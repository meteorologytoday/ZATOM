mutable struct EllipticSolver

    sop2         :: MatrixSpatialOperators
    sop3         :: MatrixSpatialOperators

    eW_length    :: Int64
    eVW_length   :: Int64

    eW_send_W    :: AbstractArray{Float64, 2}
    W_send_eW    :: AbstractArray{Float64, 2}

    eVW_send_VW     :: AbstractArray{Float64, 2}
    VW_send_eVW     :: AbstractArray{Float64, 2}
    VW_cosϕ_eVW     :: AbstractArray{Float64, 2}

    eW_LAPz_eW          :: AbstractArray{Float64, 2}
    eVW_Ψ_starop_eVW    :: AbstractArray{Float64, 2}
    
    forcing_ops
    tool

    wksp_w
    wksp_Ψ_star

    # sop2 = operator in 2D (V, W)
    # sop3 = operator in 3D (U, V, W) = 2 layers in X direction
    function EllipticSolver(
        gd2   :: Grid,
        sop2  :: MatrixSpatialOperators,
        sop3  :: MatrixSpatialOperators,
        μ     :: Float64,
    )
        cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))
        
        Ny = sop2.op.Ny
        Nz = sop2.op.Nz

        # Create eW
        # need a mask excluding bnd points
        mask_eff_W = reshape(sop2.W_mask_W * ones(Float64, sop2.op.W_pts), sop2.op.W_dim...)

        # Create coversion matrix and its inverse
        W_num = reshape(collect(1:length(mask_eff_W)), size(mask_eff_W)...)
        active_num_eff_W = W_num[ mask_eff_W .==1 ]

        eW_send_W = sop2.op.W_I_W[active_num_eff_W, :]; dropzeros!(eW_send_W)
        W_send_eW = sparse(eW_send_W')
 
        # Create eVW
        mask_eff_VW = reshape(sop2.VW_mask_VW * ones(Float64, sop2.op.VW_pts), sop2.op.VW_dim...)
        VW_num = reshape(collect(1:length(mask_eff_VW)), size(mask_eff_VW)...)
        active_num_eff_VW = VW_num[ mask_eff_VW .==1 ]

        eVW_send_VW = sop2.op.VW_I_VW[active_num_eff_VW, :]; dropzeros!(eVW_send_VW)
        VW_send_eVW = sparse(eVW_send_VW')
    
        VW_cosϕ_VW = cos.(gd2.ϕ_VW) |> cvtDiagOp
        VW_cosϕ_eVW = VW_cosϕ_VW * VW_send_eVW
         
        # identity 
        eW_I_eW   = eW_send_W   * W_send_eW
        eVW_I_eVW = eVW_send_VW * VW_send_eVW


        # Extraction matrix
        # These matrices extract either east or west boundary points.
        # It operates on 3D grids and return a slice of 2D grids.

        A = spdiagm(0 => ones(Float64, sop2.op.W_pts))
        B = spzeros(sop2.op.W_pts, sop2.op.W_pts)
        W2_extract_west_W3 = sparse([ A B ])
        W2_extract_east_W3 = sparse([ B A ])
 
        A = spdiagm(0 => ones(Float64, sop2.op.T_pts))
        B = spzeros(sop2.op.T_pts, sop2.op.T_pts)
        T2_extract_west_T3 = sparse([ A B ])
        T2_extract_east_T3 = sparse([ B A ])
 
        A = spdiagm(0 => ones(Float64, sop2.op.VW_pts))
        B = spzeros(sop2.op.VW_pts, sop2.op.VW_pts)
        VW2_extract_west_VW3 = sparse([ A B ])
        VW2_extract_east_VW3 = sparse([ B A ])
        
        #if ! (eVW_I_eVW == VW_send_eVW == eVW_send_VW)
        #    throw(ErrorException("!!"))
        #end

        # construct operators on the left-hand-side
        eW_LAPz_eW       = sparse(eW_send_W   * sop2.W_LAPz_W   * W_send_eW)
        eVW_Ψ_starop_eVW = sparse(eVW_send_VW * sop2.VW_LAPz_VW * VW_send_eVW)

        # Remove the spurious signal along the north and south boundary
        #killNS_W = reshape(ones(Float64, sop2.op.W_pts), sop2.op.W_dim...)
        #killNS_W[:, [1,sop2.op.Ny], :] .= 0.0
        #W_killNS_W = spdiagm(0 => view(killNS_W, :))

        forcing_ops = (

            we_b = sparse( 
               #- W2_extract_east_W3 * sop3.W_m_interp_V * sop3.V_invΔx_V * sop3.V_invD_V * sop3.V_f_V * sop3.V_m∂y_T
               - W2_extract_east_W3 * sop3.W_interp_V * sop3.V_invΔx_V * sop3.V_invD_V * sop3.V_f_V * sop3.V_∂y_T
            ),

            Ψb_star_b = μ * 2 * sparse(
                VW2_extract_west_VW3 * (

                    #sop3.VW_m_interp_V * sop3.V_invD_V * sop3.V_ϵ_V * sop3.V_m∂y_T
                  #  sop3.VW_interp_V * sop3.V_invD_V * sop3.V_ϵ_V * sop3.V_∂y_T
                   sop3.VW_interp_T * sop3.T_invD_T * sop3.T_f_T * sop3.T_invΔx_T

                ) - VW2_extract_east_VW3 * (

                   sop3.VW_interp_T * sop3.T_invD_T * sop3.T_f_T * sop3.T_invΔx_T

                )
            )

        )


        eW_length = size(eW_I_eW)[1]
        eVW_length = size(eVW_I_eVW)[1]
        
        tool = (
            Ψ_starop = lu(eVW_Ψ_starop_eVW),
            Lap = lu(eW_LAPz_eW),
        )
 
        wksp_w = (
            rhs_W  = zeros(Float64, sop2.op.W_pts),
            rhs_W2 = zeros(Float64, sop2.op.W_pts),
            rhs_eW = zeros(Float64, eW_length),
            lhs_eW = zeros(Float64, eW_length),
            lhs_W  = zeros(Float64, sop2.op.W_pts),
        )

        wksp_Ψ_star = (
            rhs_VW  = zeros(Float64, sop2.op.VW_pts),
            rhs_VW2 = zeros(Float64, sop2.op.VW_pts),
            rhs_eVW = zeros(Float64, eVW_length),
            lhs_eVW = zeros(Float64, eVW_length),
            lhs_VW  = zeros(Float64, sop2.op.VW_pts),
        )


       
        return new(

            sop2,
            sop3,

            eW_length,
            eVW_length,

            eW_send_W,
            W_send_eW,

            eVW_send_VW,
            VW_send_eVW,
            VW_cosϕ_eVW,

            eW_LAPz_eW,
            eVW_Ψ_starop_eVW,

            forcing_ops,
            tool,

            wksp_w,
            wksp_Ψ_star,

        ) 

    end
end

function solveΨ_star!(
    es       :: EllipticSolver,
    Ψ_star :: AbstractArray{Float64};
    bw     :: Union{Nothing, AbstractArray{Float64}} = nothing,
    be     :: Union{Nothing, AbstractArray{Float64}} = nothing,
)

    wksp = es.wksp_Ψ_star
    fops = es.forcing_ops 

    rhs_VW = fops.Ψb_star_bw * view(bw, :) + fops.Ψb_star_be * view(be, :)
    view(Ψ_star, :) .= es.VW_send_eVW * ( es.tool.Ψ_starop \ ( es.eVW_send_VW * view(rhs_VW, :) ) )

end

function solveW!(
    es     :: EllipticSolver,
    w      :: AbstractArray{Float64};
    bw     :: Union{Nothing, AbstractArray{Float64}} = nothing,
    be     :: Union{Nothing, AbstractArray{Float64}} = nothing,
    domain :: Symbol,
)

    wksp = es.wksp_w
    fops = es.forcing_ops 

    if domain == :west_boundary
        mul!(wksp.rhs_W,  fops.west_be, view(be, :)) 
        mul!(wksp.rhs_W2, fops.west_bw, view(bw, :)) 
        @. wksp.rhs_W += wksp.rhs_W2
    elseif domain == :east_boundary
        rhs_W = fops.we_be *  view(be, :)
    else
        throw(ErrorException("Unsupported domain: " * string(domain)))
    end

    view(w, :) .= es.W_send_eW * ( es.tool.Lap \ ( es.eW_send_W * view(rhs_W, :) ) )

end





