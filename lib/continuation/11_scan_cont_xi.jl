using NCDatasets
using Formatting
using ArgParse
using JSON
using Statistics
using LinearAlgebra
println("###########################################")
println("LinearAlgebra.BLAS.get_num_threads() = ", LinearAlgebra.BLAS.get_num_threads())
println("###########################################")
include("helper_functions.jl")

println("""
This program takes a snapshot file (`--input-snapshot`) to begin continuition. The initial
direction is given by `--scan-dir`. Each scan contains `--scan-count` points. 

The output file are in the folder `[CASENAME]` and file name is "dddd.nc" CASENAME
is provided as `--casename`.

On initialization, if filename  "dddd.nc" is found, then by default the last record of
the file with the largest 4 digits number is used to initialize the model and continue
the scanning; if not, then start the entire new scan from "0001.nc"

""")

s = ArgParseSettings()
@add_arg_table s begin

    "--scan-dir"
        help = "Scan direction. Either positive (>0) or negative (<0)"
        arg_type = Float64
        required = true
 
    "--scan-range"
        help = "Scan range of ξ (unitless). Scan stops when ξ is out of this range."
        nargs = 2
        arg_type = Float64
        required = true
 
    "--scales"
        help = "Scales of variables (T, S, ξ)."
        nargs = 3
        arg_type = Float64
        default = [ 1.0, 1.0, 1.0 ] 
 
    "--init-snapshot"
        help = "The snapshot file to begin with."
        arg_type = String
        required = true
 
    "--casename"
        help = "The casename of this scan. The output folder will be created using this name in the same folder as --input-snapshot"
        arg_type = String
        required = true
 
    "--ds-exp-lb"
        help = "ds_exp of continuition lower bound"
        arg_type = Float64
        default = -5.0

    "--scan-count"
        help = "Limit of count"
        arg_type = Int64
        default = 100

    "--newton-min-iter"
        help = "Newton method param."
        arg_type = Int64
        default = 3 

    "--newton-max-iter"
        help = "Newton method param."
        arg_type = Int64
        default = 20 

    "--newton-res-target"
        help = "Newton method param."
        arg_type = Float64
        default = 1e-13

end


parsed = parse_args(ARGS, s)

JSON.print(parsed, 4)

ξ_min = parsed["scan-range"][1]
ξ_max = parsed["scan-range"][2]

println("(ξ_min, ξ_max) = ($ξ_min, $ξ_max)")


if ξ_min > ξ_max
    throw(ErrorException("Values `--scan-range` should be increasing."))
end


include("../../julia/constants.jl")
include("../../julia/ZAOM.jl")
include("../../julia/ContMethod.jl")

using .ContMethod1
using .ZAOM

function getGhost(model)
    local ghost = Dict(
        :X_ => copy(model.state.X_),
        :ξ  => zeros(Float64, 1),
        :dXds_ => 0.0 * model.state.X_,
        :dξds  => zeros(Float64, 1), 
    )
    ghost[:_dXds] = reshape(ghost[:dXds_], :)
    ghost[:_X] = reshape(ghost[:X_], :)

    return ghost
end

function getXScaleVec(model, scale_T, scale_S)
    local Xscale_vec = ones(Float64, length(model.state._X))
    local T_pts = size(model.state._X_, 1)
    Xscale_vec[1:T_pts]         .= scale_T
    Xscale_vec[T_pts+1:2*T_pts] .= scale_S
    return Xscale_vec
end



function execContinuition(;
    snapshot_file  :: String,         # This file is only used to construct the model
    init_file      :: String,         # This file will contain _X, xi, _dXds, dxids
    scales         :: Array{Float64}, # Scales of T, S, ξ
    scan_steps     :: Int64,
    ds_exp_lb      :: Float64,          # This parameter controls the detail of the stepping
    output_file    :: String,
    newton_min_iter :: Int64,
    newton_max_iter :: Int64,
    newton_res_target :: Float64,
    ξ_min :: Float64,
    ξ_max :: Float64,
    record_init :: Bool = false,
)

    println("Init_file: ", init_file)
    println("Scales : ", scales)
    println("(ξ_min, ξ_max) = ($ξ_min, $ξ_max)")

    local model = ZAOM.loadModel(snapshot_file)
    #global model, ghost

    # The ghost state is used because the ContMethod module
    # lives in a linearly transformed space. So an intermediate
    # state memory will be convenient.
    local ghost = getGhost(model)

    function _F_dFdx_dFdp(x, p)
       
        F, dFdx = ZAOM.cal_F_dFdX!(model; do_dFdX=true)

        dSfdp = reshape(model.core.op_mtx[:dSfdξ], :, 1)
        dFdp = vcat(dSfdp*0, dSfdp)  # first part is temperature so that dTdp = 0

        return F, dFdx, dFdp

    end

    # This callback is expected to be called whenever the model's
    # state (T, S, ξ) is updated.
    function _callback(cmi)

        setXP!(cmi; s=cmi.s, x=ghost[:_X], p=ghost[:ξ])
        model.state._X .= ghost[:_X]
        ZAOM.updateB!(model)
        ZAOM.updateSaltTransForcing!(model; ξ=ghost[:ξ][1])
        ZAOM.diagnoseTransport!(model)
        
        setXP!(cmi; s=cmi.ṡ, x=ghost[:_dXds], p=ghost[:dξds])
    end


    local Xscale_vec =  getXScaleVec(model, scales[1], scales[2])
    local cmi = ContMethod1.CMInfo(
        Nx = length(model.state._X),
        F_dFdx_dFdp     = _F_dFdx_dFdp,
        mx              = Xscale_vec,
        mp              = scales[3],
        nwt_min_iter    = newton_min_iter,
        nwt_max_iter    = newton_max_iter,
        res             = newton_res_target,
        newton_callback = _callback,
        ds_exp_lb = ds_exp_lb,
        skip_unconverge = true,
    )

    # Loading initial state and tangent vector
    Dataset(init_file, "r") do ds

        Ns = ds.dim["Ns"]        
        println("Ns: $Ns")

        ghost[:_X]      .= reshape(ds["X_"][:, :, :, :, Ns], :)
        ghost[:ξ][1]     = ds["xi"][Ns]
        dxds = reshape(ds["dxds_"][:, :, :, :, Ns], :)
        dqds = ds["dqds"][Ns]

        model.state._X .= ghost[:_X]
        ZAOM.updateB!(model)
        ZAOM.updateSaltTransForcing!(model; ξ=ghost[:ξ][1])
       
        setS!(cmi; s=cmi.s, x=ghost[:_X],    p=ghost[:ξ])
        setS!(cmi; s=cmi.ṡ, x=dxds, p=[dqds], no_scale=true)

        cmi.ds_exp = ds["ds_exp"][Ns]
    end
        
    if model.env.ξ > ξ_max || model.env.ξ < ξ_min
        println("Now ξ = $(model.env.ξ) is not within range [$(ξ_min), $(ξ_max)]")
        println("Not doing anything. Stop program.")
        return 1
    end

    println("Output file: $(output_file)")
    Dataset(output_file, "c") do ds

        local gd = model.env.gd_bib
        
        defDim(ds, "Ns", Inf)
        defDim(ds, "Nz", gd.Nz)
        defDim(ds, "Ny", gd.Ny)
        defDim(ds, "Nx",  2   )
        defDim(ds, "NX", model.env.NX)
        
        defDim(ds, "Nzp1", gd.Nz+1)
        defDim(ds, "Nyp1", gd.Ny+1)
        defDim(ds, "Nxp1", 3)

        ds.attrib["scale_T"] = parsed["scales"][1]
        ds.attrib["scale_S"] = parsed["scales"][2]
        ds.attrib["scale_xi"] = parsed["scales"][3]
        ds.attrib["dlambda_b"] = model.env.Δλb
        ds.attrib["dlambda_i"] = model.env.Δλi
        ds.attrib["Q"]        = model.env.Q

        # coordinate variables 
        z_T    = defVar(ds, "z_T", Float64, ("Nz",))
        z_W    = defVar(ds, "z_W", Float64, ("Nzp1",))
        y_T    = defVar(ds, "y_T", Float64, ("Ny",))
        y_V    = defVar(ds, "y_V", Float64, ("Nyp1",))
        Δx_T   = defVar(ds, "dx_T", Float64, ("Ny", "Nx",))

        z_T[:] = model.env.gd_bib.z_T[:, 1, 1] 
        z_W[:] = model.env.gd_bib.z_W[:, 1, 1]
        y_T[:] = model.env.gd_bib.ϕ_T[1, :, 1]
        y_V[:] = model.env.gd_bib.ϕ_V[1, :, 1]
        Δx_T[:] = model.env.gd_bib.Δx_T[1, :, :]

        #
        X_SS_        = defVar(ds, "X_SS_",       Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
        X_VDIFU_     = defVar(ds, "X_VDIFU_",    Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
        X_HDIFU_     = defVar(ds, "X_HDIFU_",    Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))
        X_ZNLHDIFU_  = defVar(ds, "X_ZNLHDIFU_", Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))
        X_CVA_       = defVar(ds, "X_CVA_",      Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))
        X_ADV_ZOC_   = defVar(ds, "X_ADV_ZOC_",  Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))
        X_ADV_MOC_   = defVar(ds, "X_ADV_MOC_",  Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))
        X_FRC_       = defVar(ds, "X_FRC_",      Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))
        dXdt_SUM_    = defVar(ds, "dXdt_SUM_",   Float64, ("Nz", "Ny", "Nx", "NX", "Ns" ))

        u_ZOC   = defVar(ds, "u_ZOC",     Float64, ("Nz", "Ny", "Nxp1", "Ns"))
        v_ZOC   = defVar(ds, "v_ZOC",     Float64, ("Nz", "Nyp1", "Nx", "Ns"))
        w_ZOC   = defVar(ds, "w_ZOC",     Float64, ("Nzp1", "Ny", "Nx", "Ns"))

        u_MOC   = defVar(ds, "u_MOC",     Float64, ("Nz", "Ny", "Nxp1", "Ns"))
        v_MOC   = defVar(ds, "v_MOC",     Float64, ("Nz", "Nyp1", "Nx", "Ns"))
        w_MOC   = defVar(ds, "w_MOC",     Float64, ("Nzp1", "Ny", "Nx", "Ns"))



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
        χ       = defVar(ds, "chi",    Float64, ("Nzp1", "Ny", "Ns"))
        ξ       = defVar(ds, "xi",      Float64, ("Ns",))
        qe      = defVar(ds, "qe",     Float64, ("Nz", "Ny", "NX", "Ns"))
        qw      = defVar(ds, "qw",     Float64, ("Nz", "Ny", "NX", "Ns"))
        cvg     = defVar(ds, "cvg",    Float64, ("Ns",))
        res         = defVar(ds, "res",      Float64, ("Ns",))
        var_logAbsDetJ      = defVar(ds, "logAbsDetJ",  Float64, ("Ns",))
        var_logAbsDetJ_sign = defVar(ds, "logAbsDetJ_sign",  Float64, ("Ns",))
        var_stable          = defVar(ds, "stable",  Float64, ("Ns",))
        
        # Lower cases means in s space. The original data 
        # should be dxds_ * scale. This is because we need
        # to keep the precision of ds vector to achieve bit-
        # to-bit repreducability.
        dxds_   = defVar(ds, "dxds_",   Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
        dqds    = defVar(ds, "dqds",    Float64, ("Ns",))
        ds_exp  = defVar(ds, "ds_exp",  Float64, ("Ns",))
        
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
        function record!(_t; converge::Bool, residue::Float64, logAbsDetJ::Float64, logAbsDetJ_sign::Float64, stable::Float64)

            X_SS_[:, :, :, :, _t]    = model.state.X_SS_
            X_VDIFU_[:, :, :, :, _t] = model.state.X_VDIFU_
            X_HDIFU_[:, :, :, :, _t] = model.state.X_HDIFU_
            X_ZNLHDIFU_[:, :, :, :, _t] = model.state.X_ZNLHDIFU_
            X_CVA_[:, :, :, :, _t] = model.state.X_CVA_
            X_ADV_ZOC_[:, :, :, :, _t] = model.state.X_ADV_ZOC_
            X_ADV_MOC_[:, :, :, :, _t] = model.state.X_ADV_MOC_
            X_FRC_[:, :, :, :, _t] = model.state.X_FRC_
            dXdt_SUM_[:, :, :, :, _t] = model.state.dXdt_SUM_
            u_ZOC[:, :, :, _t]  = model.state.u_ZOC
            v_ZOC[:, :, :, _t]  = model.state.v_ZOC
            w_ZOC[:, :, :, _t]  = model.state.w_ZOC
            u_MOC[:, :, :, _t]  = model.state.u_MOC
            v_MOC[:, :, :, _t]  = model.state.v_MOC
            w_MOC[:, :, :, _t]  = model.state.w_MOC
     

            Ψb[:, :, _t]      = model.state.Ψb
            X_[:, :, :, :, _t] = model.state.X_[:, :, :, :]
            bw[:, :, _t]      = model.state.b[:, :, 1]
            be[:, :, _t]      = model.state.b[:, :, 2]
            diff_b[:, :, _t]  = model.state.b[:, :, 2] - model.state.b[:, :, 1]
            ww[:, :, _t]      = model.state.w[:, :, 1]
            we[:, :, _t]      = model.state.w[:, :, 2]
            vw[:, :, _t]      = model.state.v[:, :, 1]
            ui[:, :, _t]      = model.state.u[:, :, 2]
            χ[:, :, _t]       = model.state.χ[:, :]
            ξ[_t]             = model.env.ξ
            qw[:, :, :, _t]   = model.state.diag[:q_cva][:, :, 1, :]
            qe[:, :, :, _t]   = model.state.diag[:q_cva][:, :, 2, :]
 
            cvg[_t]           = ( converge ) ? 1.0 : 0.0
            res[_t]           = residue
            var_logAbsDetJ[_t]       = logAbsDetJ
            var_logAbsDetJ_sign[_t]  = logAbsDetJ_sign
            var_stable[_t]  = stable

            dxds_[:, :, :, :, _t] = reshape(cmi.ṡ[1:end-1], size(model.state.X_))
            dqds[_t]              = cmi.ṡ[end]
            ds_exp[_t]            = cmi.ds_exp

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

        local cnt = 1
        
        if record_init
            println("The flag `record_init` is on. Record the initial status.")

            _, dFdx = ZAOM.cal_F_dFdX!(model; do_dFdX=true)
            dFdx = Matrix(dFdx)
            evs = eigvals(dFdx)
            #println(evs)
            stable = ( all(real.(evs) .< 0.0) ) ? 1.0 : 0.0
            println("Stability of the new steady state: ", ( stable == 1.0 ) ? "stable" : "unstable") 

            record!(cnt; converge = true, residue=NaN, logAbsDetJ=NaN, logAbsDetJ_sign=NaN, stable=stable); cnt+=1;
        end
        
        for i = 1:scan_steps

            println(format("The {}-th search. Current ξ = {:.5f}. (max Ψb, minΨb) = ({:.3e}, {:3e}) Sv", i, model.env.ξ, maximum(model.state.Ψb) / 1e6, minimum(model.state.Ψb/1e6)))

            _converge, _residue, (_logAbsDetJ, _logAbsDetJ_sign) = ContMethod1.doContinuition!(
                cmi,
                verbose=true, 
            )

            # Determine the stability
            _, dFdx = ZAOM.cal_F_dFdX!(model; do_dFdX=true)
            dFdx = Matrix(dFdx)
            evs = eigvals(dFdx)
            #println(evs)
            stable = ( all(real.(evs) .< 0.0) ) ? 1.0 : 0.0
            println("Stability of the new steady state: ", ( stable == 1.0 ) ? "stable" : "unstable") 
        
            record!(cnt; converge = _converge, residue=_residue, logAbsDetJ=_logAbsDetJ, logAbsDetJ_sign=_logAbsDetJ_sign, stable=stable); cnt+=1;

            if model.env.ξ > ξ_max || model.env.ξ < ξ_min
                println("Now ξ = $(model.env.ξ) is not within range [$(ξ_min), $(ξ_max)]")
                break
            end

            #if i == 5
            #    JLD2.save("new_5.jld2", "cmi", cmi, "model", model)
            #end
        end
    end

    return 0
end

println("Start Program")
println("Detecting initialization...")

while true

    global status, init_file, next_file = detectInitialization(
        dirname(parsed["init-snapshot"]),
        parsed["casename"],
    )

    global next_file_tmp = joinpath(dirname(next_file), "$(basename(next_file)).tmp")

    if status == :INIT

        println("No previous scan. Create a new one: $(next_file_tmp)")

        local ξ0

        if parsed["scan-dir"] == 0.0
            throw(ErrorException("Option `--scan-dir` must not be zero"))
        elseif parsed["scan-dir"] > 0
            ξ0 = ξ_min
        elseif parsed["scan-dir"] < 0
            ξ0 = ξ_max
        end

        local init_Δξ = sign(parsed["scan-dir"]) * 1e-3


        global model = ZAOM.loadModel(parsed["init-snapshot"])
        global ghost = getGhost(model)
        
        println("1. Initial guess the first tangent vector.")
        println("(a) We do newton method to ensure that now is a steady state solution.")
        println(format("ξ0 = {:.5f}", ξ0))
        ZAOM.updateSaltTransForcing!(model; ξ=ξ0)
        ZAOM.stepModel_newton!(
            model;
            res_target  = parsed["newton-res-target"],
            min_iter    = parsed["newton-min-iter"],
            max_iter    = parsed["newton-max-iter"],
            verbose     = true,
            fastforward = true,
        )

        # Copy the first state
        local X_copy = copy(model.state.X_)

        println("(b) Natually perturbe ξ to ξ+Δξ and solve steady state again")
        println(format("init_Δξ = {:.5f}.", init_Δξ))

        ZAOM.updateSaltTransForcing!(model; ξ=ξ0 + init_Δξ)
        ZAOM.stepModel_newton!(
            model;
            res_target  = parsed["newton-res-target"],
            min_iter    = parsed["newton-min-iter"],
            max_iter    = parsed["newton-max-iter"],
            verbose     = true,
            fastforward = true,
        )


        ghost[:dXds_] .= model.state.X_ - X_copy
        ghost[:dξds][1] = init_Δξ
        
        ghost[:dXds_] .= model.state.X_ - X_copy

        # Restore model state
        println("Outputting: $(next_file_tmp)") 
        Dataset(next_file_tmp, "c") do ds

            local gd = model.env.gd_bib
            
            defDim(ds, "Ns", Inf)
            defDim(ds, "Nz", gd.Nz)
            defDim(ds, "Ny", gd.Ny)
            defDim(ds, "Nx",  2   )
            defDim(ds, "NX", model.env.NX)
 
            ds.attrib["scale_T"] = parsed["scales"][1]
            ds.attrib["scale_S"] = parsed["scales"][2]
            ds.attrib["scale_xi"] = parsed["scales"][3]
            ds.attrib["dlambda_b"] = model.env.Δλb
            ds.attrib["dlambda_i"] = model.env.Δλi
            ds.attrib["Q"]         = model.env.Q
                
            X_      = defVar(ds, "X_",      Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
            ξ       = defVar(ds, "xi",       Float64, ("Ns",))
            dxds_   = defVar(ds, "dxds_",   Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
            dqds    = defVar(ds, "dqds",    Float64, ("Ns",))
            ds_exp  = defVar(ds, "ds_exp",  Float64, ("Ns",))
           
            Xscale_vec = reshape(getXScaleVec(model, parsed["scales"][1], parsed["scales"][2]), size(X_copy))
             
            X_[:, :, :, :, 1]      = X_copy
            ξ[1]                   = ξ0
            dxds_[:, :, :, :, 1]   = ghost[:dXds_] ./ Xscale_vec
            dqds[1]                = ghost[:dξds][1] ./ parsed["scales"][3]
            ds_exp[1]              = 0.0
        end

        mv(next_file_tmp, next_file; force=true)

    elseif status == :CONT
       
        println("Initial file determined") 
        break

    end

end

println("Init the continuition with : $init_file")
println("Running continuition...")
status = execContinuition(;
    snapshot_file     = parsed["init-snapshot"],
    init_file         = init_file,
    output_file       = next_file_tmp,
    scales            = parsed["scales"],
    scan_steps        = parsed["scan-count"],
    ds_exp_lb         = parsed["ds-exp-lb"],
    newton_res_target = parsed["newton-res-target"],
    newton_min_iter   = parsed["newton-min-iter"],
    newton_max_iter   = parsed["newton-max-iter"],
    ξ_min = ξ_min,
    ξ_max = ξ_max,
    record_init = basename(init_file) == "0000.nc",
)

if status == 0
    println("Moving $(next_file_tmp) to $(next_file)")        
    mv(next_file_tmp, next_file; force=true)
else
    println("Error: execContinuition returns nonzero code $(status). Please check.")
end

println("End Program.")
