function diagnoseTransport!(
    m :: Model
)

    @fast_extract m
    

    T_V = reshape(co.sop_bib.V_interp_T * st._T, co.sop_bib.op.V_dim...)
    S_V = reshape(co.sop_bib.V_interp_T * st._S, co.sop_bib.op.V_dim...)

    flux =st.v .* ev.gd_bib.Δx_V .* ev.gd_bib.Δz_V
    # calculate the heat and salt transport
    ρc∫vT = sum( reshape(T_V .* flux * ρ0, co.sop_bib.op.V_dim...), dims=(1, 3))[1, :, 1]
    ρ∫vS  = sum( reshape(S_V .* flux * ρ0,  co.sop_bib.op.V_dim...), dims=(1, 3))[1, :, 1]
 
    st.diag[:ρc∫vT] .= ρc∫vT
    st.diag[:ρ∫vS]  .= ρ∫vS

   
    return ρc∫vT, ρ∫vS


end
