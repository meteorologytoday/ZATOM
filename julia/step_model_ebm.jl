function stepModel_ebm!(
    m :: Model,
)

    @fast_extract m 

    # Function calOp_adv! updates Ψb_star, Ψb, and vel
    op_adv        = calOp_adv!(m)

    op_vdiff      = calOp_vdiff(co.vds, st.b)
    op_hdiff      = co.op_mtx[:hdiff]
    op_forcing    = co.op_mtx[:forcing]
#    op_zonal_diff = co.op_mtx[:zonal_diffusion]

    op = op_adv + op_vdiff + op_hdiff + op_forcing 
    #op = op_adv + op_vdiff + op_hdiff + op_forcing + op_zonal_diff
    #op = op_hdiff

    b_flat = view(st.b, :)

    # *** Diffusion and restoring ***
    # (b_t+1 - b_t) / dt =  OP1 * b_t+1 + OP2 * (b_t+1 - b_target)
    # b_t+1 - b_t =  dt * OP1 * b_t+1 + dt * OP2 * (b_t+1 - b_target)
    # (I - OP1 * dt - OP2 * dt) b_t+1 = b_t - dt * OP2 * b_target
    # b_t+1 = (I - OP1 * dt - OP2 * dt) \ (b_t - dt * OP2 * b_target)
    
    F_EBM = lu( I - ev.Δt * op )
    b_flat .= F_EBM \ ( b_flat - ev.Δt * op_forcing * view(st.b_forcing, :) )

end
