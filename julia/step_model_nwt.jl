function stepModel_newton!(
    m           :: Model;
    res_target  :: Float64 = 0.0,
    min_iter    :: Integer = 0,
    max_iter    :: Integer = 0,
    verbose     :: Bool = false,
    fastforward :: Bool=false,
)
    
    @fast_extract m

    _F_dFdX = (x) -> cal_F_dFdX!(m; do_dFdX = true) 
    callback = callback=()->updateB!(m)

    if fastforward
        flag, cnt, res, Δ, J, logDetJ = doNewtonMethod!(st._X, _F_dFdX, res_target, min_iter, max_iter; verbose=verbose, callback=callback)
    else
        res, _, _, logDetJ = doNewtonMethodOnce!(st._X, _F_dFdX; verbose=verbose, callback=callback)
    end

    diagnoseTransport!(m) 


    if fastforward
        return flag, cnt, res, logDetJ
    else
        return res, logDetJ
    end
end

function cal_F_dFdX!(
    m :: Model;
    do_dFdX :: Bool = false,
)
    @fast_extract m

    op_adv_ZOC, op_adv_MOC = calOp_adv!(m)
    op_adv = op_adv_ZOC + op_adv_MOC
    
    op_adv_ZOC = nblockdiag(op_adv_ZOC, ev.NX)
    op_adv_MOC = nblockdiag(op_adv_MOC, ev.NX)
    op_adv     = nblockdiag(op_adv    , ev.NX)


    op_vdiff_cva, J_vdiff_cva = calOp_vdiff(co.vds_cva, st.b; cal_Jacobian=do_dFdX)
    op_vdiff_cva   = nblockdiag(op_vdiff_cva, ev.NX)
    
    op_vdiff_iso = nblockdiag(co.op_mtx[:vdiff_iso], ev.NX)
    
    op_hdiff   = nblockdiag(co.op_mtx[:hdiff], ev.NX)
    op_hdiff_zonal = nblockdiag(co.op_mtx[:hdiff_zonal], ev.NX)
    op_forcing = co.op_mtx[:forcing]
    
    op = op_adv_ZOC + op_adv_MOC + op_vdiff_iso + op_vdiff_cva + op_hdiff + op_hdiff_zonal + op_forcing
    F = op * st._X - op_forcing * st._X_tgt + st._X_SS
   
    # Compute and save the convective adjustment and diffusion
    st._X_CVA[:]       .= op_vdiff_cva * st._X
    st._X_VDIFU[:]     .= op_vdiff_iso * st._X
    st._X_HDIFU[:]     .= op_hdiff * st._X
    st._X_ZNLHDIFU[:]  .= op_hdiff_zonal * st._X
    st._X_ADV_ZOC[:]   .= op_adv_ZOC * st._X
    st._X_ADV_MOC[:]   .= op_adv_MOC * st._X
    st._X_FRC[:]       .= op_forcing * (st._X - st._X_tgt)
    
    @. st._dXdt_SUM = (
        st._X_ADV_ZOC 
        + st._X_ADV_MOC 
        + st._X_VDIFU 
        + st._X_CVA 
        + st._X_HDIFU 
        + st._X_ZNLHDIFU 
        + st._X_SS 
        + st._X_FRC
    )
 
    if do_dFdX

        ∂b∂X = zeros(Float64, ev.NX)
        ∂b∂X[1] =   g * α
        ∂b∂X[2] = - g * β
            
        first_term  = op_adv
        second_term = Array{Any}(undef, ev.NX)

        for k = 1:ev.NX # Each row
            row = Array{Any}(undef, ev.NX)
            PX = co.op_mtx[:P] * st._X_[:, k]
            Σ∇PXS = - co.op_mtx[:D∇] * spdiagm( 0 => PX ) * co.op_mtx[:S]
            for kk = 1:ev.NX # Same row, different columns
                row[kk] = ∂b∂X[kk] * Σ∇PXS 
            end
            second_term[k] = hcat(row...)
        end
       
        second_term = vcat(second_term...)

        J_adv = first_term + second_term
        
        J_vdiff_cva   = nblockdiag(J_vdiff_cva, ev.NX) # It has already been computed above
        J_vdiff_iso   = op_vdiff_iso
        J_hdiff       = op_hdiff
        J_hdiff_zonal   = op_hdiff_zonal
        J_forcing = co.op_mtx[:forcing]

        dFdX = J_adv + J_vdiff_iso + J_vdiff_cva + J_hdiff + J_hdiff_zonal + J_forcing
    end


    if ! do_dFdX
        return F
    else
        return F, dFdX
    end
end



