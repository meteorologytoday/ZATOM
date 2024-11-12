using NCDatasets
using Formatting
using ArgParse
using JSON
using Statistics
using LinearAlgebra
using JLD2

include("ZAOM/julia/constants.jl")
include("ZAOM/julia/ZAOM.jl")
include("ZAOM/julia/ContMethod.jl")

using .ContMethod1

s = ArgParseSettings()
@add_arg_table s begin

    "--scan-dir"
        help = "Scan direction. Either positive (>0) or negative (<0)"
        arg_type = Float64
        required = true
 
    "--scan-range"
        help = "Scan range of Q in unit of Sv. Scan stops when Q is out of this range."
        nargs = 2
        arg_type = Float64
        required = true
 
    "--scales"
        help = "Scan range of Q in unit of Sv. Scan stops when Q is out of this range."
        nargs = 3
        arg_type = Float64
        default = [ 1e-2, 1e-3, 1e4 ] 
 
    "--input"
        arg_type = String
        required = true
 
    "--output"
        arg_type = String
        required = true
 
    "--ds-exp-lb"
        help = "ds_exp of continuition lower bound"
        arg_type = Float64
        default = -5.0

    "--record-limit"
        help = "Limit of count"
        arg_type = Int64
        default = 2000

      
end

parsed = parse_args(ARGS, s)

JSON.print(parsed, 4)

# original 1e-13
res_target   = 1e-13
nwt_min_iter = 3
nwt_max_iter = 10

Δs_fixed = 1.0

#scales = [ 1e-2, 1e-3, 0.01 * 1e6 ] 
scales = parsed["scales"]

Ns = Inf
Q_min = parsed["scan-range"][1] * 1e6
Q_max = parsed["scan-range"][2] * 1e6

if Q_min > Q_max
    throw(ErrorException("Values `--scan-range` should be increasing."))
end


if parsed["scan-dir"] == 0.0
    throw(ErrorException("Option `--scan-dir` must not be zero"))
elseif parsed["scan-dir"] > 0
    Q0 = Q_min
elseif parsed["scan-dir"] < 0
    Q0 = Q_max
end

init_ΔQ = sign(parsed["scan-dir"]) * 0.01e6


model = ZAOM.loadModel(parsed["input"])

# Initiate model (update buoyancy and such)
ZAOM.updateB!(model)
ZAOM.updateSaltTransForcing!(model)

# The ghost state is used because the ContMethod module
# lives in a linearly transformed space. So an intermediate
# state memory will be convenient.
ghost = Dict(
    :_X => copy(model.state._X),
    :Q  => zeros(Float64, 1),
    :dXds_ => 0.0 * model.state.X_,
    :dQds  => zeros(Float64, 1), 
)
ghost[:_dXds] = reshape(ghost[:dXds_], :)

function _F_dFdx_dFdp(x, p)
   
    F, dFdx = ZAOM.cal_F_dFdX!(model; do_dFdX=true)

    dSfdp = reshape(model.core.op_mtx[:dSfdQ], :, 1)
    dFdp = vcat(dSfdp*0, dSfdp)  # first part is temperature

    return F, dFdx, dFdp

end

# This callback is expected to be called whenever the model's
# state (T, S, Q) is updated.
function _callback(cmi)

    global ghost
    setXP!(cmi; s=cmi.s, x=ghost[:_X], p=ghost[:Q])
    model.state._X .= ghost[:_X]
    ZAOM.updateB!(model)
    ZAOM.updateSaltTransForcing!(model, ghost[:Q][1])
    ZAOM.diagnoseTransport!(model)
    
    setXP!(cmi; s=cmi.ṡ, x=ghost[:_dXds], p=ghost[:dQds])
end


Xscale_vec = ones(Float64, length(model.state._X))
T_pts = size(model.state._X_, 1)
Xscale_vec[1:T_pts]         .= scales[1]
Xscale_vec[T_pts+1:2*T_pts] .= scales[2]

cmi = ContMethod1.CMInfo(
    Nx = length(model.state._X),
    F_dFdx_dFdp     = _F_dFdx_dFdp,
    mx              = Xscale_vec,
    mp              = scales[3],
    nwt_min_iter    = nwt_min_iter,
    nwt_max_iter    = nwt_max_iter,
    res             = res_target,
    newton_callback = _callback,
    ds_exp_lb = parsed["ds-exp-lb"],
    skip_unconverge = true,
)

println("1. Initial guess the first tangent vector.")
println("(a) We do newton method to ensure that now is a steady state solution.")
println("Q0 = ", Q0/1e6, " Sv")
ZAOM.updateSaltTransForcing!(model, Q0)
ZAOM.stepModel_newton!(
    model;
    res_target  = res_target,
    min_iter    = nwt_min_iter,
    max_iter    = 20,#nwt_max_iter,
    verbose     = true,
    fastforward = true,
)

# Copy the first state
X_copy = copy(model.state._X)
Q0 = model.env.Q

println("(b) Natually perturbe Q to Q+ΔQ and solve steady state again")
println("init_ΔQ = ", init_ΔQ/1e6, " Sv")

ZAOM.updateSaltTransForcing!(model, Q0 + init_ΔQ)
ZAOM.stepModel_newton!(
    model;
    res_target  = res_target,
    min_iter    = nwt_min_iter,
    max_iter    = 20,#nwt_max_iter,
    verbose     = true,
    fastforward = true,
)


ghost[:_dXds] .= model.state._X - X_copy
ghost[:dQds][1] = init_ΔQ


# Set the tangent vector ṡ
setS!(cmi; s=cmi.ṡ, x=ghost[:_dXds], p=ghost[:dQds]) 

# Assign initial states.
model.state._X .= X_copy
ZAOM.updateB!(model)
ZAOM.updateSaltTransForcing!(model, Q0)

ghost[:_X]     .= X_copy
ghost[:Q][1]    = Q0
setS!(cmi; s=cmi.s, x=ghost[:_X], p=ghost[:Q]) 



JLD2.save("old_cont.jld2", "cmi", cmi, "model", model)

Dataset(parsed["output"], "c") do ds

    gd = model.env.gd_bib
    
    defDim(ds, "Ns", Inf)
    defDim(ds, "Nz", gd.Nz)
    defDim(ds, "Ny", gd.Ny)
    defDim(ds, "Nx",  2   )
    defDim(ds, "NX", model.env.NX)
    
    defDim(ds, "Nzp1", gd.Nz+1)
    defDim(ds, "Nyp1", gd.Ny+1)
    defDim(ds, "Nxp1", 3)
    
    # coordinate variables 
    z_T    = defVar(ds, "z_T", Float64, ("Nz",))
    z_W    = defVar(ds, "z_W", Float64, ("Nzp1",))
    y_T    = defVar(ds, "y_T", Float64, ("Ny",))
    y_V    = defVar(ds, "y_V", Float64, ("Nyp1",))
    Δx_T   = defVar(ds, "dx_T", Float64, ("Nx",))

    z_T[:] = model.env.gd_bib.z_T[:, 1, 1] 
    z_W[:] = model.env.gd_bib.z_W[:, 1, 1]
    y_T[:] = model.env.gd_bib.ϕ_T[1, :, 1]
    y_V[:] = model.env.gd_bib.ϕ_V[1, :, 1]
    Δx_T[:] = model.env.gd_bib.Δx_T[1, 1, :]

    #
    X_SS_   = defVar(ds, "X_SS_", Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))

    # simulation variables            

    Ψb      = defVar(ds, "Psib",   Float64, ("Nzp1", "Nyp1", "Ns"))
    X_      = defVar(ds, "X_",     Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
    bw      = defVar(ds, "bw",     Float64, ("Nz", "Ny", "Ns"))
    be      = defVar(ds, "be",     Float64, ("Nz", "Ny", "Ns"))
    diff_b  = defVar(ds, "diff_b", Float64, ("Nz", "Ny", "Ns"))
    ww      = defVar(ds, "ww",     Float64, ("Nzp1", "Ny", "Ns"))
    we      = defVar(ds, "we",     Float64, ("Nzp1", "Ny", "Ns"))
    vw      = defVar(ds, "vw",     Float64, ("Nz", "Nyp1", "Ns"))
    ui      = defVar(ds, "ui",     Float64, ("Nz", "Ny", "Ns"))
    Q       = defVar(ds, "Q",      Float64, ("Ns",))
    cvg     = defVar(ds, "cvg", Float64, ("Ns",))
    res     = defVar(ds, "res", Float64, ("Ns",))
    
    dXds_   = defVar(ds, "dXds_",   Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
    dQds    = defVar(ds, "dQds",    Float64, ("Ns",))
    
    # Diagnostic information
    ∫X      = defVar(ds, "int_X",    Float64, ("NX", "Ns",))
    ∫ADV    = defVar(ds, "int_ADV",  Float64, ("NX", "Ns",))
    ∫HDIF   = defVar(ds, "int_HDIF", Float64, ("NX", "Ns",))
    ∫VDIF   = defVar(ds, "int_VDIF", Float64, ("NX", "Ns",))
    ∫SS     = defVar(ds, "int_SS",   Float64, ("NX", "Ns",))
    ∫ERR    = defVar(ds, "int_ERR",  Float64, ("NX", "Ns",))
    ρc∫vT   = defVar(ds, "heat_transport", Float64, ("Nyp1", "Ns"))
    ρ∫vS    = defVar(ds, "salt_transport", Float64, ("Nyp1", "Ns"))

    ∫Δv = model.env.gd_bib.∫Δv
    function record!(_t; converge::Bool, residue::Float64)

        X_SS_[:, :, :, :, _t] = model.state.X_SS_
        Ψb[:, :, _t]      = model.state.Ψb
        X_[:, :, :, :, _t] = model.state.X_[:, :, :, :]
        bw[:, :, _t]      = model.state.b[:, :, 1]
        be[:, :, _t]      = model.state.b[:, :, 2]
        diff_b[:, :, _t]  = model.state.b[:, :, 2] - model.state.b[:, :, 1]
        ww[:, :, _t]      = model.state.w[:, :, 1]
        we[:, :, _t]      = model.state.w[:, :, 2]
        vw[:, :, _t]      = model.state.v[:, :, 1]
        ui[:, :, _t]      = model.state.u[:, :, 2]
        Q[_t]             = model.env.Q
        cvg[_t]           = ( converge ) ? 1.0 : 0.0
        res[_t]           = residue

        dXds_[:, :, :, :, _t] = ghost[:dXds_]
        dQds[_t]              = ghost[:dQds][1]

        ∫X[:, _t]         = model.state.diag[:∫X]     / ∫Δv
        ∫ADV[:, _t]       = model.state.diag[:∫ADV]   / ∫Δv
        ∫HDIF[:, _t]      = model.state.diag[:∫HDIF]  / ∫Δv
        ∫VDIF[:, _t]      = model.state.diag[:∫VDIF]  / ∫Δv
        ∫SS[:, _t]        = model.state.diag[:∫SS]    / ∫Δv
        ∫ERR[:, _t]       = model.state.diag[:∫ERR]   / ∫Δv
        
        ρc∫vT[:, _t]      = model.state.diag[:ρc∫vT]        
        ρ∫vS[:, _t]       = model.state.diag[:ρ∫vS]        
        
        println(format("Record the {:d}-th result. ", _t))
    end

    cnt = 1
#    record!(cnt); cnt+=1;

    for i = 1:parsed["record-limit"]

        global ghost
        
        println(format("The {}-th search. Current Q = {:.2e} Sv. max Ψb = {:.2e} Sv", i, model.env.Q / 1e6, maximum(model.state.Ψb) / 1e6))

        if model.env.Q > Q_max || model.env.Q < Q_min
            println("Now Q is not within range [$(Q_min), $(Q_max)]")
            break
        end

        #= Natural continuition 
        ZAOM.updateSaltTransForcing!(model, model.env.Q + 0.1*1e6)
        ZAOM.stepModel_newton!(
            model;
            res_target  = res,
            min_iter    = nwt_min_iter,
            max_iter    = nwt_max_iter,
            verbose     = true,
            fastforward = true,
        )
        =#

        _converge, _residue = ContMethod1.doContinuition!(
            cmi,
            verbose=true, 
        )

        record!(cnt; converge = _converge, residue=_residue); cnt+=1;
    end
end

println("End Program.")




