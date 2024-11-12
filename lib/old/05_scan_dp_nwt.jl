using NCDatasets
using Formatting
using ArgParse
include("ZAOM/julia/constants.jl")
include("ZAOM/julia/ZAOM.jl")

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
 
    "--scan-pts"
        help = "The points of Q_vec"
        arg_type = Int64
        default = 21

    "--scan-days"
        help = "Simulation time of each given Q in days."
        arg_type = Int64
        required = true
 
    "--input"
        arg_type = String
        required = true
 
    "--output"
        arg_type = String
        required = true
        
end

parsed = parse_args(ARGS, s)

nwt_residue  = 1e-13
nwt_min_iter = 3
nwt_max_iter = 20

Δs_fixed = 1.0
scales = [ 1e-2, 1e-3, 0.01 * 1e6 ] 

Ns = Inf
Q_min = parsed["scan-range"][1] * 1e6
Q_max = parsed["scan-range"][2] * 1e6

if Q_min > Q_max
    throw(ErrorException("Values `--scan-range` should be increasing."))
end


if parsed["scan-dir"] == 0.0
    throw(ErrorException("Option `--scan-dir` must not be zero"))
end

Q_vec = collect(range(Q_min, Q_max, length=parsed["scan-pts"]))

if sign(parsed["scan-dir"]) < 0
    reverse!(Q_vec)
end


model = ZAOM.loadModel(parsed["input"])

# Initiate model (update buoyancy and such)
ZAOM.updateB!(model)
ZAOM.updateSaltTransForcing!(model)



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
    function record!(_t)

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

    callback = function(k)
        ZAOM.diagnoseTransport!(model)
        record!(k)
        println(format("Now Q = {:.3f} Sv", model.env.Q / 1e6))
    end

    break_test = function(k)
        println(format("Now Q = {:.3f} Sv", model.env.Q / 1e6))
        return model.env.Q > Q_max || model.env.Q < Q_min
    end


    model.env.Δt = 86400.0
    for i=1:length(Q_vec)
        
        Q_now = Q_vec[i]
        print(format("Step: {:d} with Q = {:.3f} Sv \n", i, Q_now / 1e6))
        model.env.Q = Q_now
        ZAOM.updateB!(model)
        ZAOM.updateSaltTransForcing!(model)

        _X = copy(model.state._X) 
 
        for trial=1:5

            local flag, cnt, res

            for t=1:parsed["scan-days"]
                print(format("Step: {:d}/{:d}\r", t, parsed["scan-days"]))
                ZAOM.stepModel_rk2!(model)
            end


            try


                flag, cnt, res = ZAOM.stepModel_newton!(
                    model;
                    res_target = nwt_residue,
                    min_iter = nwt_min_iter,
                    max_iter = nwt_max_iter,
                    verbose = true,
                    fastforward = true,
                )

            if flag == :converge
                break
            else
                model.state._X .= _X            
                continue
            end


            catch e
                if isa(err, LinearAlgebra.SingularException)
                    model.state._X .= _X            
                    continue
                else
                    throw(e)
                end
            end

        end

        record!(i)

    end
end

println("End Program.")

