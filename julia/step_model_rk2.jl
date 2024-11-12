function stepModel_rk2!(
    m :: Model,
)

    @fast_extract m 

    # Step 1: Advection using RK2

    ∫X_bef = co.diag_mtx[:∫dv] * st._X
    
    _X_copy = copy(st._X)
    _b_copy = copy(st._b)

    op_adv_ZOC, op_adv_MOC = calOp_adv!(m)
    op_adv = op_adv_ZOC + op_adv_MOC
    
    op_adv_ZOC = nblockdiag(op_adv_ZOC, ev.NX)
    op_adv_MOC = nblockdiag(op_adv_MOC, ev.NX)
    op_adv     = nblockdiag(op_adv    , ev.NX)

    k1 = op_adv * st._X
    @. st._X += (ev.Δt/2.0) * k1
    
    # need to update buoyancy because velocity is a function of buoyancy
    updateB!(m)
 
    op_adv_ZOC, op_adv_MOC = calOp_adv!(m)
    op_adv = op_adv_ZOC + op_adv_MOC
    
    op_adv_ZOC = nblockdiag(op_adv_ZOC, ev.NX)
    op_adv_MOC = nblockdiag(op_adv_MOC, ev.NX)
    op_adv     = nblockdiag(op_adv    , ev.NX)
    k2 = op_adv * st._X

    @. st._X = _X_copy + ev.Δt * k2
    
    # update for consistency
    updateB!(m)
    
    ∫X_aft_adv = co.diag_mtx[:∫dv] * st._X
    st.diag[:∫ADV] .= ∫X_aft_adv - ∫X_bef
    
    #@. st._X_ADV = (st._X - _X_copy) / ev.Δt

    # Step 2: diffusion and restoring
   
    op_vdiff_cva, _ = calOp_vdiff(co.vds_cva, _b_copy)
 
    op_vdiff_cva  = nblockdiag(op_vdiff_cva, ev.NX)
    op_vdiff_iso  = nblockdiag(co.op_mtx[:vdiff_iso], ev.NX)
    op_hdiff      = nblockdiag(co.op_mtx[:hdiff], ev.NX)
    op_hdiff_zonal= nblockdiag(co.op_mtx[:hdiff_zonal], ev.NX)
    op_forcing    = co.op_mtx[:forcing]

    # *** Diffusion and restoring ***
    # (b_t+1 - b_t) / dt =  OP1 * b_t+1 + OP2 * (b_t+1 - b_target) + const
    # b_t+1 - b_t =  dt * OP1 * b_t+1 + dt * OP2 * (b_t+1 - b_target) + dt * const 
    # (I - OP1 * dt - OP2 * dt) b_t+1 = b_t - dt * OP2 * b_target + dt * const
    # b_t+1 = (I - OP1 * dt - OP2 * dt) \ (b_t - dt * OP2 * b_target + dt * const)

    # Test conservation
    op = op_vdiff_iso + op_vdiff_cva + op_hdiff + op_hdiff_zonal + op_forcing
    
    F_EBM = lu( I - ev.Δt * op )
    st._X .= F_EBM \ ( st._X - ev.Δt * op_forcing * st._X_tgt + ev.Δt * st._X_SS)
    
    updateB!(m)
    
    # Compute and save the convective adjustment and diffusion
    st.diag[:q_cva][:] .= op_vdiff_cva * st._X
    st._X_VDIFU[:] .= op_vdiff_iso * st._X
    st._X_HDIFU[:] .= op_hdiff * st._X


    st.diag[:∫HDIF] .= ev.Δt .* co.diag_mtx[:∫dv] * op_hdiff * st._X
    st.diag[:∫VDIF] .= ev.Δt .* co.diag_mtx[:∫dv] * (op_vdiff_iso + op_vdiff_cva) * st._X
    st.diag[:∫SS]   .= ev.Δt .* co.diag_mtx[:∫dv] * op_forcing * (st._X - st._X_tgt)
    ∫X_aft = co.diag_mtx[:∫dv] * st._X

    st.diag[:∫X]   .= ∫X_aft
    @. st.diag[:∫ERR] = (∫X_aft - ∫X_aft_adv) - (st.diag[:∫HDIF] + st.diag[:∫VDIF] + st.diag[:∫SS])

    diagnoseTransport!(m) 

end
