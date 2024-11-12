function stepModel_ab3!(
   m :: Model;
)

    @fast_extract m

    reset!(co.wksp)
    
    solveVelocity!(m)

    stepSlice!(m, st.bw, st.vw, st.ww, st.G_bw, st.bw_forcing)
    stepSlice!(m, st.be, st.ve, st.we, st.G_be, st.be_forcing)

    
    # Zonal diffusion between bw and be
    # This is relevant to the boundary condition
    stepZonalDiffusion!(m, st.b)

    # Important: Adam-Bashforth
    rotateG_idx!(co.G_idx)
    
end

function stepSlice!(
    m   :: Model,
    b   :: AbstractArray{Float64, 2}, 
    v   :: AbstractArray{Float64, 2}, 
    w   :: AbstractArray{Float64, 2}, 
    G_b :: AbstractArray{Float64, 3}, 
    b_forcing :: AbstractArray{Float64, 2}, 
)
    @fast_extract m
    

    # Advection

    doAdvection!(;
        sop   = co.sop,
        wksp  = co.wksp,
        G_b   = G_b,
        G_idx = co.G_idx,
        Δt    = ev.Δt,
        b     = b,
        v     = v,
        w     = w,
    )

    # Euler Backward Method:
    #
    # *** Diffusion and restoring ***
    # (b_t+1 - b_t) / dt =  OP1 * b_t+1 + OP2 * (b_t+1 - b_target)
    # b_t+1 - b_t =  dt * OP1 * b_t+1 + dt * OP2 * (b_t+1 - b_target)
    # (I - OP1 * dt - OP2 * dt) b_t+1 = b_t - dt * OP2 * b_target
    # b_t+1 = (I - OP1 * dt - OP2 * dt) \ (b_t - dt * OP2 * b_target)
    #

    b_flat = view(b, :)
    eT_send_T = co.vds.eT_send_T
    T_send_eT = co.vds.T_send_eT
    eT_I_eT   = co.vds.eT_I_eT

    op_hdiff = co.op_mtx[:2][:hdiff]
    op_vdiff = calOp2_vdiff(co.vds, b_flat)
    op_forcing = co.op_mtx[:2][:forcing] 
    
    EBM = lu( eT_I_eT - ev.Δt * ( op_hdiff + op_vdiff + op_forcing) )

    b_flat .= T_send_eT * ( EBM \ (eT_send_T * (b_flat - ev.Δt * op_forcing * view(b_forcing, :)) ) )

end

function stepZonalDiffusion!(
    m   :: Model,
    b   :: AbstractArray{Float64, 3}, 
)

    @fast_extract m
    
    b_flat = view(b, :)

    # b_t+1 - b_t / dt =  OP * b_t+1
    # (I - OP * dt) b_t+1 = b_t
    # b_t+1 = ( I - OP * dt ) \ b_t

    b_flat .= co.lu_mtx[:EBM_zonal_diffusion] \ b_flat

end


function solveVelocity!(
    m :: Model,
)
    @fast_extract m

    gd = ev.gd
    
    # Solve we, and Ψb
    solveΨ_star!(co.es, st.Ψb_star;  bw=st.bw, be=st.be) 
    solveW!(co.es, st.we;            be=st.be, domain=:east_boundary) 
    Ψ_star_to_vel!(co.sop, w=st.ww, v=st.vw, Ψ_star=st.Ψb_star) 
   
    # Solve ui
    st.ui[:] .= gd.Δx_T[:] .* ( co.sop.T_DIVz_W * st.we[:] )

    st.Ψb  .= st.Ψb_star .* gd.Δx_VW
    st.ww .-= st.we

end

function doAdvection!(;
    sop   :: MatrixSpatialOperators,
    wksp  :: Workspace,
    G_b   :: AbstractArray{Float64, 3},
    G_idx :: Dict,
    Δt    :: Float64,
    b     :: AbstractArray{Float64, 2},
    v     :: AbstractArray{Float64, 2},
    w     :: AbstractArray{Float64, 2},
)

    Δt0 = G_idx[:now]
    Δt1 = G_idx[:one_Δt_ago]
    Δt2 = G_idx[:two_Δt_ago]

    G_b_Δt0 = view(G_b, :, :, Δt0)
    G_b_Δt1 = view(G_b, :, :, Δt1)
    G_b_Δt2 = view(G_b, :, :, Δt2)

    cal_u∇b!(;
        sop   = sop,
        wksp  = wksp,
        u∇b   = G_b_Δt0,
        b     = b,
        v     = v,
        w     = w,
    )

    @. G_b_Δt0 *= -1
    @. b += Δt * ABIII(G_b_Δt0, G_b_Δt1, G_b_Δt2)

end


function cal_u∇b!(;
    sop   :: MatrixSpatialOperators,
    wksp  :: Workspace,
    u∇b  :: AbstractArray{Float64},
    b     :: AbstractArray{Float64},
    v     :: AbstractArray{Float64},
    w     :: AbstractArray{Float64},
)
    tmp_V = getSpace!(wksp, :V)
    tmp_W = getSpace!(wksp, :W)
    tmp_T = getSpace!(wksp, :T)

    mul_autoflat!(tmp_V, sop.V_∂y_T, b)
    mul_autoflat!(tmp_W, sop.W_∂z_T, b)

    @. tmp_V *= v
    @. tmp_W *= w

    mul_autoflat!(u∇b,   sop.T_interp_V, tmp_V)
    mul_autoflat!(tmp_T, sop.T_interp_W, tmp_W)

    @. u∇b += tmp_T

end

function Ψ_star_to_vel!(
    es      :: MatrixSpatialOperators;
    w       :: AbstractArray{Float64},
    v       :: AbstractArray{Float64},
    Ψ_star  :: AbstractArray{Float64},
)

    mul_autoflat!(w, es.W_DIVy_VW, Ψ_star)
    mul_autoflat!(v, es.V_DIVz_VW, Ψ_star)

    v .*= -1
 
end


