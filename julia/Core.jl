mutable struct Core

    sop_b         :: MatrixSpatialOperators
    sop_bb        :: MatrixSpatialOperators
    sop_bib       :: MatrixSpatialOperators

    es            :: EllipticSolver
    vds_iso       :: VerticalDiffusionSolver
    vds_cva       :: VerticalDiffusionSolver
    hds           :: HorizontalDiffusionSolver
    wksp          :: Workspace

    op_mtx        :: Any         
    diag_mtx      :: Any 

    function Core(env)
        
        NX = env.NX
        Nx = 2
        Ny = env.Ny_bsn
        Nz = env.Nz

        D_factor = env.Δλb / (env.Δλb + env.Δλi)
        println("Dilution factor: ", D_factor)

        println("Making matricies")
        @time sop_b   = MatrixSpatialOperators(env.gd_b)
        @time sop_bb  = MatrixSpatialOperators(env.gd_bb)

        println("Making bib")
        @time sop_bib = MatrixSpatialOperators(env.gd_bib)

        es = EllipticSolver(env.gd_b, sop_b, sop_bb, env.μ)
        vds_iso = VerticalDiffusionSolver(sop_bb; K_iso=env.Kv, K_cva=env.Kv, Δ=env.cva_Δ, dilution_factor=D_factor)
        vds_cva = VerticalDiffusionSolver(sop_bb; K_iso=0.0,    K_cva=env.Kv_cva - env.Kv, Δ=env.cva_Δ, dilution_factor=D_factor)
        hds = HorizontalDiffusionSolver(sop_bib; K=env.Kh, dilution_factor=D_factor)
        
        wksp = Workspace(;
            Ny = Ny,
            Nz = Nz,
        ) 

        hdiv     = zeros(Float64, 1, Ny)

        # construct forcing
        forcing_coes = Array{Any, 1}(undef, NX)
            
        forcing_coes[1] = genForcingCoe(sop_bib.op.T_dim, env.tracer_forcings[1]; MLT=env.MLT_T, MLT_shape=env.MLT_shape, gd=env.gd_bib)
        forcing_coes[2] = genForcingCoe(sop_bib.op.T_dim, env.tracer_forcings[2]; MLT=env.MLT_S, MLT_shape=env.MLT_shape, gd=env.gd_bib)

        if NX >= 3
            for i=3:NX
                forcing_coes[i] = genForcingCoe(sop_bib.op.T_dim, env.tracer_forcings[i]; MLT=env.MLT_T, MLT_shape=env.MLT_shape, gd=env.gd_bib)
            end
        end
        forcing = blockdiag( Tuple(forcing_coes)... )
        
        # Some additional matrices for Newton method 

        # calculate projection matrix T-grid to V- and W-grid
        # TODO: make your own U_interp_T
#                sop_bb.op.U_W_T   ;
#                sop_bb.U_interp_T ;
        P = [
                sop_bb.op.U_W_T   ;
                sop_bb.V_interp_T ;
                sop_bb.W_interp_T
        ]

        ∇ = sparse([ sop_bb.T_DIVx_U sop_bb.T_DIVy_V  sop_bb.T_DIVz_W ])
        
        # calculate the dilute matrix which only
        # operates on the eastern (i.e. second slab)
        D = ones(Float64, sop_bb.op.T_pts)
        D[(1+Nz*Ny):(Nz*Ny*2)] .= D_factor
        D = spdiagm(0 => D)
        
        D∇ = sparse(D * ∇)
        
        S, G, Q = calSGQ(es)
        EXTRACT_ZOC, EXTRACT_MOC = calEXTRACT_MOCZOC(sop_b)


        #### Old formulation #####
        # project from U,V,W grid onto T
        P_old = blockdiag(spdiagm(sop_bb.op.T_pts, sop_bb.op.U_pts), sop_bb.T_interp_V, sop_bb.T_interp_W)
        Σ_old = spdiagm(0 => ones(Float64, sop_bb.op.T_pts))
        Σ_old = [ Σ_old Σ_old Σ_old ]

        ΣP_old = Σ_old * P_old
        ∇_old = sparse([ spdiagm(sop_bb.op.U_pts, sop_bb.op.T_pts) ; sop_bb.V_∂y_T ; sop_bb.W_∂z_T ])



        op_vdiff_iso, _ = calOp_vdiff(vds_iso, zeros(Float64, sop_bb.op.T_dim...); cal_Jacobian=false)

        op_mtx = Dict{Symbol, Any}( 
#                :zonal_diffusion => calOp_zonal_diff(sopz, env.Khx),
#                :∇       => ∇,
                :P           => P,
                :D∇          => D∇,
                :S           => S, 
                :G           => G, 
                :Q           => Q, 
                :hdiff       => calOp_hdiff(hds),
                :hdiff_zonal => calOp_zonal_diff(sop_bib, env.Khx, env.gd_bib.Δx_T[1, :, 1]),
                :forcing     => forcing,
                :∇_old       => ∇_old,
                :ΣP_old      => ΣP_old,
                :vdiff_iso   => op_vdiff_iso,
                :EXTRACT_MOC => EXTRACT_MOC,
                :EXTRACT_ZOC => EXTRACT_ZOC,
        )
        
         
        _tmp = sparse(reshape(env.gd_bib.Δv_T, 1, :))
        diag_mtx = Dict(
            :∫dv => blockdiag( ntuple((a,)-> _tmp, NX)... )
        )

       

        new(
            sop_b,
            sop_bb,
            sop_bib,
#            sopz,
            es,
            vds_iso,
            vds_cva,
            hds,
            wksp,
            op_mtx,
            diag_mtx,
        )
        
    end
end


function genForcingCoe(
    dim        :: Tuple,
    tracer_forcing :: TracerForcing;
    gd  :: Union{Grid, Nothing} = nothing,
    MLT :: Float64 = -1.0,
    MLT_shape :: String,
)
    forcing_coe = - ones(Float64, dim...) / tracer_forcing.τ
    if tracer_forcing.style == :NONE
        forcing_coe .= 0.0
    elseif tracer_forcing.style == :SURFACE
        forcing_coe .*= genMixedlayerWeight(gd.z_T, MLT, MLT_shape)
        #forcing_coe[2:end, :, :] .= 0  # Only restore the surface
    elseif tracer_forcing.style == :EVERYWHERE
        # nothing to do
    else
        throw(ErrorException("Unrecognized `style` parameter: " * string(style)))
    end

    return spdiagm( 0 => reshape(forcing_coe, :) )
end

function genMixedlayerWeight(z :: Union{AbstractArray, Float64}, MLT :: Float64, MLT_shape :: String)

    if MLT <= 0.0
        throw(ErrorException("MLT should be a positive number"))
    end

    if MLT_shape == "linear"
        pos = (x,) -> max(x, 0.0)
        w = pos.(1.0 .+ z ./ MLT)
    elseif MLT_shape == "step"
        pos = (x,) -> (x > - MLT) ? 1.0 : 0.0
        w = pos.(z)
    else
        throw(ErrorException("Unknown MLT_shape: $MLT_shape (only `linear` and `step` are allowed)."))
    end
    
    return w

end
