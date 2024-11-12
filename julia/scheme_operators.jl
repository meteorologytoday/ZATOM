function calOp_zonal_diff(
    sop   :: MatrixSpatialOperators, 
    K     :: AbstractArray{Float64, 2}, # (Nz, Ny)
    Δx    :: AbstractArray{Float64, 1}, # (Ny,)
)

    mop = sop.op

    if mop.Nx != 2
        throw(ErrorException("Zonal diffusion is used only with Nx = 2"))
    end
  
    invΔx = repeat( reshape(Δx, 1, mop.Ny), outer=(mop.Nz, 1) )[:].^(-1)

    U_invΔx_U = spdiagm(0 => [ invΔx ; invΔx ; invΔx ])

    op_chg = sop.T_DIVx_U * sop.U_mask_U * spdiagm(0 => [K[:]*0 ; K[:] ; K[:]*0]) * U_invΔx_U * (mop.U_W_T - mop.U_E_T)

    return op_chg

end

# This function updates Ψb_star, Ψb, and vel.
# The derived advection operator is returned.
function calOp_adv!(
    m :: Model;
    update :: Bool = true,
)

    @fast_extract m

    _vel = co.op_mtx[:S] * st._b
    _vel_ZOC = co.op_mtx[:EXTRACT_ZOC] * co.op_mtx[:S] * st._b
    _vel_MOC = co.op_mtx[:EXTRACT_MOC] * co.op_mtx[:S] * st._b


    if update
        st._vel        .= _vel
        st._vel_ZOC    .= _vel_ZOC
        st._vel_MOC    .= _vel_MOC

        st.Ψb_star[:] .= (co.op_mtx[:Q] * st._b)[1:length(st.Ψb_star)]
        st.Ψb         .= st.Ψb_star .* view(ev.gd_b.Δx_VW, :, :, 1)
        st.χ          .= - view(st.w, :, :, 2) .* view(ev.gd_b.Δx_W, :, :, 1)
    end


    #println(size(co.op_mtx[:DΣ∇]))
    #println(size(VEL))
    #println(size(co.op_mtx[:P]))


    # new formulation
    return (- co.op_mtx[:D∇] * spdiagm( 0 => _vel_ZOC ) * co.op_mtx[:P] ), (- co.op_mtx[:D∇] * spdiagm( 0 => _vel_MOC ) * co.op_mtx[:P] ) 
            

    # old formulation
    # return - co.op_mtx[:ΣP_old] * spdiagm( 0 => VEL ) * co.op_mtx[:∇_old]

end



"""
    return matrix : G 

    such that G [ ψ^* ; we ] = [ ui; vw ; ww ; ve ; we ]
"""
function calG(
    msp      :: MatrixSpatialOperators;
)

    op = msp.op
    V_I_V = spdiagm( 0 => ones(op.V_pts) )
    W_I_W = spdiagm( 0 => ones(op.W_pts) )

    # [ ψ^* ; we ] |-> [ vw ; ve (=0) ; ww + we ; we ]
    A = blockdiag( 

        [
            -msp.V_DIVz_VW ;
            spzeros(op.V_pts, op.VW_pts) ;
            msp.W_DIVy_VW ;
        ],
    
        sparse(1.0I, op.W_pts, op.W_pts)
    )

#    r = zeros(op.V_pts + op.W_pts)
#    r[op.V_pts+1:end] .= -1.0

    # [ vw ; ve ; ww + we ; we ] |-> [ vw ; ve ; ww ; we ]


    W_substract_W = spdiagm(
         0        =>   ones(op.W_pts*2),
         op.W_pts => - ones(op.W_pts  ), 
    )

    B = blockdiag(V_I_V, V_I_V, W_substract_W)
    


    # [ vw ; ve ; ww ; we ] |-> [ vw ; ve ; ww ; we ; ui ]
    C = blockdiag(V_I_V, V_I_V, W_I_W)
    C = blockdiag(
        C,
        [
            W_I_W ;
            msp.T_Δx_T * msp.T_DIVz_W
        ]
    )

    # [ vw ; ve ; ww ; we ; ui ] |-> [ uw=0; ui; ue=0 ; vw ; ve ; ww ; we ]
    m = 3 * op.T_pts + 2 * op.W_pts + 2 * op.V_pts
    n =     op.T_pts + 2 * op.W_pts + 2 * op.V_pts
    skip_vw_cnt = 2 * op.W_pts + 2 * op.V_pts
    get_ui = zeros(Float64, 2 * op.T_pts)
    get_ui[op.T_pts:end] .= 1.0
    D = spdiagm(m, n, (-op.T_pts*3) => ones(Float64, skip_vw_cnt), (skip_vw_cnt-op.T_pts) => get_ui)

    return dropzeros(D * C * B * A)
end

"""
    return matrix : Q 
    such that Q [ bw ; be ] = [ ψ^* ;  we ]
"""
function calQ(es :: EllipticSolver)

    fops = es.forcing_ops 
        
    # [ bw ; be ] |-> [ ψ^* ;  we ]
    return sparse([ 
        es.VW_send_eVW * inv(es.eVW_Ψ_starop_eVW |> lu) * ( es.eVW_send_VW * fops.Ψb_star_b ) ;
        es.W_send_eW * inv(es.eW_LAPz_eW |> lu) * ( es.eW_send_W * fops.we_b ) ;
    ])
end



"""
    return matrix : S 
    such that S [ bw ; be ] = [ uw=0; ui; ue=0 ; vw ; ve ; ww ; we ]
"""
function calS(es :: EllipticSolver)
    # [ bw ; be ] |-> [ ψ^* ;  we ] |-> [ uw=0; ui; ue=0 ; vw ; ve ; ww ; we ]
    return calG(es.sop2) * calQ(es)
end


"""
    return matrix : S, G, Q 
"""
function calSGQ(es :: EllipticSolver)
    return calS(es), calG(es.sop2), calQ(es)
end


"""
    return matrices : EXTRACT_ZOC, EXTRACT_MOC
    such that EXTRACT_ZOC [ uw ; ui ; ue ; vw ; ve ; ww ; we ] = [ uw=0; ui  ; ue=0 ; vw=0 ; ve=0 ; w_ZOC=-we   ; we   ]
              EXTRACT_MOC [ uw ; ui ; ue ; vw ; ve ; ww ; we ] = [ uw=0; ui=0; ue=0 ; vw   ; ve=0 ; w_MOC=ww+we ; we=0 ]
"""
function calEXTRACT_MOCZOC(msp::MatrixSpatialOperators)
    
    op = msp.op
    m = 3 * op.T_pts + 2 * op.W_pts + 2 * op.V_pts
    
    # EXTRACT_ZOC
    # A [ uw ; ui ; ue ; vw ; ww ; ve ; we ] = [ uw=0; ui    ; ue=0 ; vw=0 ; ww=0      ; ve=0 ; we     ]
    A = blockdiag( 
        sparse(0.0I, op.T_pts, op.T_pts), # clear uw
        sparse(1.0I, op.T_pts, op.T_pts), # keep  ui
        sparse(0.0I, op.T_pts, op.T_pts), # clear ue
        sparse(0.0I, op.V_pts, op.V_pts), # clear vw
        sparse(0.0I, op.V_pts, op.V_pts), # clear ve
        sparse(0.0I, op.W_pts, op.W_pts), # clear ww
        sparse(1.0I, op.W_pts, op.W_pts), # keep  we
    )

    # B [ uw ; ui ; ue ; vw ; ww ; ve ; we ] = [ we ]
    skip_cnt = 3 * op.T_pts + 2 * op.W_pts + 2 * op.V_pts
    only_we = [ spzeros(op.W_pts, 3 * op.T_pts + 2 * op.V_pts + op.W_pts) sparse(1.0I, op.W_pts, op.W_pts) ]
    B = [
        spzeros(3 * op.T_pts, m) ;
        spzeros(op.V_pts, m)     ;
        spzeros(op.V_pts, m)     ;
        -1.0 * only_we           ; # multiply by -1.0
        spzeros(op.W_pts, m)     ;
    ]
    
    EXTRACT_ZOC = A + B
    EXTRACT_MOC = I - EXTRACT_ZOC

    return EXTRACT_ZOC, EXTRACT_MOC

end


