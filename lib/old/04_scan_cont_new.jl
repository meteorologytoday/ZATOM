using NCDatasets
using Formatting
using ArgParse
using JSON
using Statistics
using LinearAlgebra

include("ZAOM/julia/constants.jl")
include("ZAOM/julia/ZAOM.jl")
include("ZAOM/julia/ContMethod.jl")

using .ContMethod1
using .ZAOM

function execContinuition(;
    snapshot_file  :: String,         # This file is only used to construct the model
    init_file      :: String,         # This file will contain _X, Q, _dXds, dQds
    init_idx       :: Int64,          # Specify which scan step is used as the initial state. Negative means counting from the end
    scales         :: Array{Float64}, # Scales of T, S, Q
    scan_steps     :: Int64,
    ds_exp_lb      :: Int64,          # This parameter controls the detail of the stepping
    output_file    :: String,
    nwt_params     :: Tuple,          # (iter_min, iter_max, target_residue)
)

    local model = ZAOM.loadModel(snapshot_file)

    # Initiate model (update buoyancy and such)
    ZAOM.updateB!(model)
    ZAOM.updateSaltTransForcing!(model)

    # The ghost state is used because the ContMethod module
    # lives in a linearly transformed space. So an intermediate
    # state memory will be convenient.
    local ghost = Dict(
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

        setXP!(cmi; s=cmi.s, x=ghost[:_X], p=ghost[:Q])
        model.state._X .= ghost[:_X]
        ZAOM.updateB!(model)
        ZAOM.updateSaltTransForcing!(model, ghost[:Q][1])
        ZAOM.diagnoseTransport!(model)
        
        setXP!(cmi; s=cmi.ṡ, x=ghost[:_dXds], p=ghost[:dQds])
    end


    local Xscale_vec = ones(Float64, length(model.state._X))
    local T_pts = size(model.state._X_, 1)
    Xscale_vec[1:T_pts]         .= scales[1]
    Xscale_vec[T_pts+1:2*T_pts] .= scales[2]

    local cmi = ContMethod1.CMInfo(
        Nx = length(model.state._X),
        F_dFdx_dFdp     = _F_dFdx_dFdp,
        mx              = Xscale_vec,
        mp              = scales[3],
        nwt_min_iter    = nwt_params[1],
        nwt_max_iter    = nwt_params[2],
        res             = nwt_params[3],
        newton_callback = _callback,
        ds_exp_lb = ds_exp_lb,
        skip_unconverge = true,
    )


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

        local cnt = 1
        
        for i = 1:scan_steps

            println(format("The {}-th search. Current Q = {:.2e} Sv. (max Ψb, minΨb) = ({:.3e}, {:3e}) Sv", i, model.env.Q / 1e6, maximum(model.state.Ψb) / 1e6i, minimum(model.state.Ψb/1e6)))

            _converge, _residue = ContMethod1.doContinuition!(
                cmi,
                verbose=true, 
            )

            record!(cnt; converge = _converge, residue=_residue); cnt+=1;
        end
    end

end
println("End Program.")




