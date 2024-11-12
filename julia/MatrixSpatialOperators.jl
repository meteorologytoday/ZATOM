mutable struct MatrixSpatialOperators

    op       :: MatrixOperators

    T_DIVx_U    :: AbstractArray{Float64, 2}

    T_DIVy_V    :: AbstractArray{Float64, 2}
    VW_DIVz_V   :: AbstractArray{Float64, 2}

    VW_DIVy_W   :: AbstractArray{Float64, 2}
    W_DIVy_VW   :: AbstractArray{Float64, 2}
    V_DIVz_VW   :: AbstractArray{Float64, 2}

    T_DIVz_W    :: AbstractArray{Float64, 2}
 
    V_DIVy_T    :: AbstractArray{Float64, 2}
    W_DIVz_T    :: AbstractArray{Float64, 2}
    
    V_∂y_T      :: AbstractArray{Float64, 2}
    W_∂z_T      :: AbstractArray{Float64, 2}

    T_∂x_U      :: AbstractArray{Float64, 2}
    T_∂y_V      :: AbstractArray{Float64, 2}
    T_∂z_W      :: AbstractArray{Float64, 2}
    V_∂z_VW     :: AbstractArray{Float64, 2}
    
    T_LAPy_T      :: AbstractArray{Float64, 2}
    W_LAPz_W      :: AbstractArray{Float64, 2}
    VW_LAPz_VW    :: AbstractArray{Float64, 2}

    V_m∂y_T     :: AbstractArray{Float64, 2}
    T_mDIVy_V    :: AbstractArray{Float64, 2}
    T_mLAPy_T    :: AbstractArray{Float64, 2}

    U_interp_T  :: AbstractArray{Float64, 2} 
    V_interp_T  :: AbstractArray{Float64, 2} 
    W_interp_T  :: AbstractArray{Float64, 2} 

    T_interp_V  :: AbstractArray{Float64, 2} 
    W_interp_V  :: AbstractArray{Float64, 2} 
    
    T_interp_W  :: AbstractArray{Float64, 2} 
    V_interp_W  :: AbstractArray{Float64, 2} 
    


    VW_interp_T :: AbstractArray{Float64, 2} 
    VW_interp_V :: AbstractArray{Float64, 2} 
    VW_interp_W :: AbstractArray{Float64, 2} 
    
    W_interp_VW :: AbstractArray{Float64, 2} 
    
    W_m_interp_V  :: AbstractArray{Float64, 2} 
    VW_m_interp_V  :: AbstractArray{Float64, 2} 

    T_mask_T       :: AbstractArray{Float64, 2}
    U_mask_U       :: AbstractArray{Float64, 2}
    V_mask_V       :: AbstractArray{Float64, 2}
    W_mask_W       :: AbstractArray{Float64, 2}
    VW_mask_VW       :: AbstractArray{Float64, 2}
    T_bordermask_T :: AbstractArray{Float64, 2}

    T_Δx_T :: AbstractArray{Float64, 2}

    T_invΔx_T :: AbstractArray{Float64, 2}
    V_invΔx_V :: AbstractArray{Float64, 2}
    
    T_Δv_T :: AbstractArray{Float64, 2}
    T_invΔv_T :: AbstractArray{Float64, 2}
    
    T_f_T :: AbstractArray{Float64, 2}
    V_f_V :: AbstractArray{Float64, 2}
    
    T_ϵ_T :: AbstractArray{Float64, 2}
    V_ϵ_V :: AbstractArray{Float64, 2}
 
    T_invD_T :: AbstractArray{Float64, 2}
    V_invD_V :: AbstractArray{Float64, 2}
  
    V_mΔay_V
    mΔay_V
 
    function MatrixSpatialOperators(
        gd :: Grid;
    )


        Nx = gd.Nx
        Ny = gd.Ny
        Nz = gd.Nz
        
        cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))
        
        @time op = MatrixOperators(Ny=Ny, Nz=Nz, Nx=Nx)
        
        
        mask_flat = view(gd.mask, :)

        onU_if_unblocked_west_onT = op.U_E_T  * mask_flat
        onU_if_unblocked_east_onT = op.U_W_T  * mask_flat

        onV_if_unblocked_north_onT = op.V_S_T  * mask_flat
        onV_if_unblocked_south_onT = op.V_N_T  * mask_flat
        onV_m_if_unblocked_north_onT = op.V_mS_T * mask_flat
        onV_m_if_unblocked_south_onT = op.V_mN_T * mask_flat
 
        onW_if_unblocked_up_onT    = op.W_DN_T * mask_flat
        onW_if_unblocked_dn_onT    = op.W_UP_T * mask_flat

        onVW_if_unblocked_north_dn_onT  = op.VW_S_W * op.W_UP_T * mask_flat
        onVW_if_unblocked_north_up_onT  = op.VW_S_W * op.W_DN_T * mask_flat
        onVW_if_unblocked_south_dn_onT  = op.VW_N_W * op.W_UP_T * mask_flat
        onVW_if_unblocked_south_up_onT  = op.VW_N_W * op.W_DN_T * mask_flat

        U_mask = onU_if_unblocked_east_onT .* onU_if_unblocked_west_onT
        V_mask = onV_if_unblocked_north_onT .* onV_if_unblocked_south_onT
        V_m_mask = onV_m_if_unblocked_north_onT .* onV_m_if_unblocked_south_onT
        W_mask = onW_if_unblocked_up_onT    .* onW_if_unblocked_dn_onT
        VW_mask = ( 
            onVW_if_unblocked_north_dn_onT
         .* onVW_if_unblocked_north_up_onT
         .* onVW_if_unblocked_south_dn_onT
         .* onVW_if_unblocked_south_up_onT
        )

        T_mask_T   = mask_flat |> cvtDiagOp
        U_mask_U   = U_mask    |> cvtDiagOp
        V_mask_V   = V_mask    |> cvtDiagOp
        V_m_mask_V = V_m_mask  |> cvtDiagOp
        W_mask_W   = W_mask    |> cvtDiagOp
        VW_mask_VW = VW_mask   |> cvtDiagOp

        T_bordermask_T = cvtDiagOp( 
               (op.T_N_T  * mask_flat)
            .* (op.T_S_T  * mask_flat)
            .* (op.T_UP_T * mask_flat)
            .* (op.T_DN_T * mask_flat)
        )

        # ===== [ BEGIN grid length, area and volume ] =====
        
        Δx_T = gd.Δx_T
        Δx_V = gd.Δx_V
        Δx_W = gd.Δx_W
        Δx_VW = gd.Δx_VW
 
        Δz_T = gd.Δz_T
        Δz_U = gd.Δz_U
        Δz_V = gd.Δz_V
        Δz_W = gd.Δz_W
        Δz_VW = gd.Δz_VW
 
        Δy_T = gd.Δy_T
        Δy_U = gd.Δy_U
        Δy_V = gd.Δy_V
        Δy_W = gd.Δy_W
        Δy_VW = gd.Δy_VW
        
        Δax_U  = Δy_U  .* Δz_U
        
        Δay_T  = Δx_T  .* Δz_T
        Δay_V  = Δx_V  .* Δz_V
        Δay_W  = Δx_W  .* Δz_W
        Δay_VW = Δx_VW .* Δz_VW

        Δaz_T  = Δx_T  .* Δy_T
        Δaz_V  = Δx_V  .* Δy_V
        Δaz_W  = Δx_W  .* Δy_W
        Δaz_VW = Δx_VW .* Δy_VW

        Δv_T   = Δx_T  .* Δy_T  .* Δz_T
        Δv_V   = Δx_V  .* Δy_V  .* Δz_V
        Δv_W   = Δx_W  .* Δy_W  .* Δz_W
        Δv_VW  = Δx_VW .* Δy_VW .* Δz_VW
        
        U_Δax_U   = Δax_U  |> cvtDiagOp

        T_Δay_T   = Δay_T  |> cvtDiagOp
        V_Δay_V   = Δay_V  |> cvtDiagOp
        W_Δay_W   = Δay_W  |> cvtDiagOp
        VW_Δay_VW = Δay_VW |> cvtDiagOp

        T_Δaz_T   = Δaz_T  |> cvtDiagOp
        V_Δaz_V   = Δaz_V  |> cvtDiagOp
        W_Δaz_W   = Δaz_W  |> cvtDiagOp
        VW_Δaz_VW = Δaz_VW |> cvtDiagOp

        T_Δv_T     = Δv_T  |> cvtDiagOp
        T_invΔv_T  = (Δv_T.^(-1))  |> cvtDiagOp
        V_invΔv_V  = (Δv_V.^(-1))  |> cvtDiagOp
        W_invΔv_W  = (Δv_W.^(-1))  |> cvtDiagOp
        VW_invΔv_VW  = (Δv_VW.^(-1))  |> cvtDiagOp
       
        T_Δx_T = ( Δx_T ) |> cvtDiagOp

        T_invΔx_T = ( Δx_T.^(-1) ) |> cvtDiagOp
        T_invΔy_T = ( Δy_T.^(-1) ) |> cvtDiagOp
        T_invΔz_T = ( Δz_T.^(-1) ) |> cvtDiagOp
 
        V_invΔx_V = ( Δx_V.^(-1) ) |> cvtDiagOp
        V_invΔy_V = ( Δy_V.^(-1) ) |> cvtDiagOp
        V_invΔz_V = ( Δz_V.^(-1) ) |> cvtDiagOp

        W_invΔx_W = ( Δx_W.^(-1) ) |> cvtDiagOp
        W_invΔy_W = ( Δy_W.^(-1) ) |> cvtDiagOp
        W_invΔz_W = ( Δz_W.^(-1) ) |> cvtDiagOp
 
        VW_invΔx_VW = ( Δx_VW.^(-1) ) |> cvtDiagOp
        VW_invΔy_VW = ( Δy_VW.^(-1) ) |> cvtDiagOp
        VW_invΔz_VW = ( Δz_VW.^(-1) ) |> cvtDiagOp
 
        # ===== [ END grid length, area and volume ] =====

        # ===== [ BEGIN interpolation ] =====
        function selfDivision(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse
            for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg:_end] ./= wgt
                end
            end
            
            return transpose(m_t) |> sparse
        end
        
        ones_T  = ones(Float64, op.T_pts)
        ones_V  = ones(Float64, op.V_pts)
        ones_W  = ones(Float64, op.W_pts)
        ones_VW = ones(Float64, op.VW_pts)
 
        U_interp_T = (op.U_W_T + op.U_E_T) * T_mask_T
        U_interp_T = selfDivision(U_interp_T, ones_T)
        
        V_interp_T = (op.V_S_T + op.V_N_T) * T_mask_T
        V_interp_T = selfDivision(V_interp_T, ones_T)
        
        W_interp_T = (op.W_DN_T + op.W_UP_T) * T_mask_T
        W_interp_T = selfDivision(W_interp_T, ones_T)
 
        VW_interp_T = (op.VW_N_W + op.VW_S_W) * (op.W_UP_T + op.W_DN_T) * T_mask_T
        VW_interp_T = selfDivision(VW_interp_T, ones_T)
 
        T_interp_V = (op.T_N_V + op.T_S_V) * V_mask_V
        T_interp_V = selfDivision(T_interp_V, ones_V)
       
        W_interp_V = (op.W_UP_T + op.W_DN_T) * (op.T_N_V + op.T_S_V) * V_mask_V
        W_interp_V = selfDivision(W_interp_V, ones_V)

        T_interp_W = (op.T_DN_W + op.T_UP_W) * W_mask_W
        T_interp_W = selfDivision(T_interp_W, ones_W)
 
        V_interp_W = (op.V_N_T + op.V_S_T) * (op.T_DN_W + op.T_UP_W) * W_mask_W
        V_interp_W = selfDivision(V_interp_W, ones_W)
 
        VW_interp_V = (op.VW_UP_V + op.VW_DN_V) * V_mask_V
        VW_interp_V = selfDivision(VW_interp_V, ones_V)

        
        VW_interp_W = (op.VW_S_W + op.VW_N_W) * W_mask_W
        VW_interp_W = selfDivision(VW_interp_W, ones_W)

        W_interp_VW = (op.W_S_VW + op.W_N_VW) * VW_mask_VW
        W_interp_VW = selfDivision(W_interp_VW, ones_VW)


        W_m_interp_V  = (op.W_UP_T + op.W_DN_T) * (op.T_N_V + op.T_S_V) * V_m_mask_V
        W_m_interp_V  = selfDivision(W_m_interp_V, ones_V)

        VW_m_interp_V = (op.VW_UP_V + op.VW_DN_V) * V_m_mask_V
        VW_m_interp_V = selfDivision(VW_m_interp_V, ones_V)



        
        # ===== [ END interpolation ] =====

        # MAGIC!!
        T_DIVx_U  = T_mask_T   * T_invΔv_T   * ( op.T_W_U  - op.T_E_U  ) * U_Δax_U  ; dropzeros!(T_DIVx_U);

        T_DIVy_V  = T_mask_T   * T_invΔv_T   * ( op.T_S_V  - op.T_N_V  ) * V_Δay_V  ; dropzeros!(T_DIVy_V);
        VW_DIVz_V = VW_mask_VW * VW_invΔv_VW * ( op.VW_DN_V - op.VW_UP_V ) * V_Δaz_V  ; dropzeros!(VW_DIVz_V);
        T_DIVz_W  = T_mask_T   * T_invΔv_T   * ( op.T_DN_W - op.T_UP_W ) * W_Δaz_W  ; dropzeros!(T_DIVz_W);

        VW_DIVy_W = VW_invΔv_VW * ( op.VW_S_W  - op.VW_N_W ) * W_Δay_W    ; dropzeros!(VW_DIVy_W);

        W_DIVy_VW = W_mask_W * W_invΔv_W   * ( op.W_S_VW  - op.W_N_VW ) * VW_Δay_VW  ; dropzeros!(W_DIVy_VW);
        #W_DIVy_VW = W_mask_W * W_invΔy_W * ( op.W_S_VW  - op.W_N_VW )  ; dropzeros!(W_DIVy_VW);
        #W_DIVy_VW =2*W_mask_W * W_invΔv_W   * ( I  - op.W_N_VW ) * VW_Δay_VW  ; dropzeros!(W_DIVy_VW);
        #W_DIVy_VW = 2 * W_mask_W * W_invΔv_W   * ( op.W_S_VW  - I ) * VW_Δay_VW  ; dropzeros!(W_DIVy_VW);
        V_DIVz_VW = V_invΔv_V   * ( op.V_DN_VW  - op.V_UP_VW ) * VW_Δaz_VW  ; dropzeros!(V_DIVz_VW);
        
        V_DIVy_T = V_mask_V * V_invΔv_V * ( op.V_S_T  - op.V_N_T  ) * T_Δay_T  ; dropzeros!(V_DIVy_T);
        W_DIVz_T = W_mask_W * W_invΔv_W * ( op.W_DN_T - op.W_UP_T ) * T_Δaz_T  ; dropzeros!(W_DIVz_T);
        #W_DIVz_T = W_mask_W * W_invΔz_W * ( op.W_DN_T - op.W_UP_T )  ; dropzeros!(W_DIVz_T);

        V_∂y_T = V_mask_V * V_invΔy_V * (op.V_S_T  - op.V_N_T)                 ; dropzeros!(V_∂y_T);
        W_∂z_T = W_mask_W * W_invΔz_W * (op.W_DN_T - op.W_UP_T)                ; dropzeros!(W_∂z_T);

        T_∂x_U = T_mask_T * T_invΔx_T * ( op.T_W_U - op.T_E_U )               ; dropzeros!(T_∂x_U);
        T_∂y_V = T_mask_T * T_invΔy_T * ( op.T_S_V - op.T_N_V )               ; dropzeros!(T_∂y_V);
        T_∂z_W = T_mask_T * T_invΔz_T * ( op.T_DN_W - op.T_UP_W )             ; dropzeros!(T_∂z_W);
        
        V_∂z_VW = V_mask_V * V_invΔz_V * ( op.V_DN_VW - op.V_UP_VW )          ; dropzeros!(V_∂z_VW);
        

        T_LAPy_T   =  T_DIVy_V * V_∂y_T
        W_LAPz_W   =  W_DIVz_T * T_∂z_W
        VW_LAPz_VW = VW_DIVz_V * V_∂z_VW
      
        #= 
        # Check vertical divergence
        x1 = dropzeros(W_mask_W * W_invΔz_W * (op.W_DN_T - op.W_UP_T))
        #x1 = dropzeros(W_mask_W * (op.W_DN_T - op.W_UP_T))
        x2 = dropzeros( T_invΔv_T * ( op.T_DN_W - op.T_UP_W ) * W_Δaz_W * x1)
#        x2 = dropzeros(( op.T_DN_W - op.T_UP_W ) * x1)
        x3 = dropzeros(sparse(ones(1, op.T_pts)) * T_Δv_T * x2)
        if length(x3.nzval) != 0
            println("length: ", length(x1.nzval))
            println(x2)
#            println(x2.nzval)
            throw(ErrorException("Check vertical divergence: Something is wrong"))
        end
#        throw(ErrorException("Stop here"))
        =#


        # This mLAPy is used for diffusion that needs to conserve total tracer
        # The north and south boundaries on two sides communicate through different widths
        # so we need to take extra care here
        mΔay_V = copy(Δay_V)
        
        if op.Nx == 2

            mΔay_V[:, end, 2] .= mΔay_V[:, end, 1]
            mΔay_V[:, 1,   2] .= mΔay_V[:, 1  , 1]

        end

        V_mΔay_V = cvtDiagOp(mΔay_V)
        
        V_m∂y_T   = V_invΔy_V * (op.V_mS_T  - op.V_mN_T)                          ; dropzeros!(V_m∂y_T);
        T_mDIVy_V = T_mask_T  * T_invΔv_T * ( op.T_S_V  - op.T_N_V  ) * V_mΔay_V  ; dropzeros!(T_mDIVy_V);
#        T_mLAPy_T = T_mDIVy_V * V_m∂y_T                                           ; dropzeros!(T_mLAPy_T);
        T_mLAPy_T = T_mDIVy_V * V_m∂y_T                                           ; dropzeros!(T_mLAPy_T);

        #=
        if Nx == 2
            println("!!! Nx == 2 !!!")

            x1 = dropzeros(sparse(ones(1, op.T_pts) * ( op.T_S_V  - op.T_N_V  ) * (op.V_mS_T - op.V_mN_T)))
            x2 = dropzeros(sparse(ones(1, op.T_pts) * (T_Δv_T * T_invΔv_T) * ( op.T_S_V  - op.T_N_V  ) * V_mΔay_V * V_invΔy_V * (op.V_mS_T - op.V_mN_T)))
            #x2 = dropzeros(sparse(ones(1, op.T_pts) * (T_Δv_T * T_invΔv_T) * ( op.T_S_V  - op.T_N_V  ) * V_mΔay_V * (op.V_mS_T - op.V_mN_T)))
            if length(x2.nzval) != 0
                println("length: ", length(x2.nzval))
                println(x2.nzval)
                throw(ErrorException("Something is wrong"))
            end
            #throw(ErrorException("Stop here"))
        end
        =#

        f_T   = 2 * gd.Ω * sin.(gd.ϕ_T)
        f_V   = 2 * gd.Ω * sin.(gd.ϕ_V)
        
        ϵ_T   = 0 * f_T .+ gd.ϵ
        ϵ_V   = 0 * f_V .+ gd.ϵ
        
        invD_T = (f_T.^2 + ϵ_T.^2).^(-1)
        invD_V = (f_V.^2 + ϵ_V.^2).^(-1)
         
        T_f_T = 2 * gd.Ω * sin.(gd.ϕ_T) |> cvtDiagOp
        V_f_V = 2 * gd.Ω * sin.(gd.ϕ_V) |> cvtDiagOp

        T_ϵ_T = ϵ_T |> cvtDiagOp
        V_ϵ_V = ϵ_V |> cvtDiagOp
        
        T_invD_T = invD_T |> cvtDiagOp
        V_invD_V = invD_V |> cvtDiagOp
            
        return new(
            op,
            
            T_DIVx_U,

            T_DIVy_V,
            VW_DIVz_V,

            VW_DIVy_W,
            W_DIVy_VW,
            V_DIVz_VW,

            T_DIVz_W,
         
            V_DIVy_T,
            W_DIVz_T,
            
            V_∂y_T,
            W_∂z_T,

            T_∂x_U,
            T_∂y_V,
            T_∂z_W,
            V_∂z_VW,
            
            T_LAPy_T,
            W_LAPz_W,
            VW_LAPz_VW,

            V_m∂y_T,
            T_mDIVy_V,
            T_mLAPy_T,

            U_interp_T,
            V_interp_T,
            W_interp_T,

            T_interp_V,
            W_interp_V,
            
            T_interp_W,
            V_interp_W,

            VW_interp_T,
            VW_interp_V,
            VW_interp_W,
            
            W_interp_VW,

            W_m_interp_V,
            VW_m_interp_V,
            
            T_mask_T,
            U_mask_U,
            V_mask_V,
            W_mask_W,
            VW_mask_VW,
            T_bordermask_T,

            T_Δx_T,

            T_invΔx_T,
            V_invΔx_V,
            
            T_Δv_T,
            T_invΔv_T,
            
            T_f_T,
            V_f_V,
            
            T_ϵ_T,
            V_ϵ_V,

            T_invD_T,
            V_invD_V,

            V_mΔay_V, 
            mΔay_V, 
#            Δay_V, 
        )
    end
end
