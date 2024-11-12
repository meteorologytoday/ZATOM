mutable struct VerticalDiffusionSolver
    
    sop :: MatrixSpatialOperators

    eT_length    :: Int64
   
    eT_send_T    :: AbstractArray{Float64, 2}
    T_send_eT    :: AbstractArray{Float64, 2}
    eT_I_eT      :: AbstractArray{Float64, 2}
    
    K_iso :: Float64
    K_cva :: Float64

    Δ        :: Float64
    dilution_factor :: Float64
    dilution_matrix :: AbstractArray{Float64, 2}
    K_func   :: Function
    dK_func  :: Function

    wksp

    function VerticalDiffusionSolver(
        sop   :: MatrixSpatialOperators;
        K_iso :: Float64,
        K_cva :: Float64,
        Δ     :: Float64,
        dilution_factor :: Float64 = 1.0,
    )

        if K_cva < K_iso 
            throw(ErrorException("Kva should be equal or greater than Kv"))
        elseif K_cva == K_iso
            println("No convective adjustment") 
        end

        cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))
        Nz = sop.op.Nz
        Ny = sop.op.Ny
        Nx = sop.op.Nx
         
        # Create coversion matrix and its inverse
        T_num        = reshape(sop.T_mask_T * collect(1:sop.op.T_pts), Nz, Ny, Nx)
        active_T_num = T_num[ T_num .!= 0.0 ]
        eT_send_T    = sop.op.T_I_T[active_T_num, :]
        T_send_eT    = eT_send_T' |> sparse
        eT_I_eT      = eT_send_T * T_send_eT

        eT_length = size(eT_I_eT)[1]

        K = spdiagm(0 => ones(Float64, sop.op.W_pts))
        

        if Nx != 2
            throw(ErrorException("VerticalDiffusionSolver.jl only supports Nx = 2"))
        end

        dilution_matrix = zeros(Float64, Nz, Ny, Nx)
        dilution_matrix[:, :, 1] .= 1.0
        dilution_matrix[:, :, 2] .= dilution_factor
        dilution_matrix = cvtDiagOp(dilution_matrix)

        # iEBD = inverse of Euler Backward Diffusion ( I - Δt * ∇K∇ ) 
        wksp = (
            above_W    = zeros(Float64, sop.op.W_pts),
            below_W    = zeros(Float64, sop.op.W_pts),
            Δ_W        = zeros(Float64, sop.op.W_pts),
            W_K_W      = spdiagm(0 => ones(Float64, sop.op.W_pts)),
            W_K∇_T     = spzeros(Float64, sop.op.W_pts, sop.op.T_pts),
            T_∇K∇_T    = spzeros(Float64, sop.op.T_pts, sop.op.T_pts),
            eT_∇K∇_T   = spzeros(Float64, eT_length, sop.op.T_pts),
            eT_∇K∇_eT  = spzeros(Float64, eT_length, eT_length),
            eT_iEBD_eT = spzeros(Float64, eT_length, eT_length),
            rhs_eT     = zeros(Float64, eT_length),
            lhs_eT     = zeros(Float64, eT_length),
        )


        # diffusivity function

        # Step function
        if Δ == 0
            K_func  = (x,) -> (x >= 0) ? K_iso : K_cva
            dK_func = (x,) -> 0.0
        elseif Δ > 0.0

            ΔK = K_cva - K_iso
            K_func  = (x,) -> ΔK *  σ(x/Δ + 1.0) + K_iso
            dK_func = (x,) -> ΔK * dσ(x/Δ + 1.0) / Δ
         
            #K_func  = (x,) -> 0.5 * ( 1.0 - tanh(2.0*x/Δ) ) * (K_cva - K_iso) + K_iso
            #dK_func = (x,) -> - sech(2.0*x/Δ)^2.0 / Δ * (K_cva - K_iso)
        else
            throw(ErrorException("Δ must be non-negative."))
        end

        return new(
            sop,
            eT_length,
            eT_send_T,
            T_send_eT,
            eT_I_eT,
            K_iso,
            K_cva,
            Δ,
            dilution_factor,
            dilution_matrix,
            K_func,
            dK_func,
            wksp,
        )

    end
end

function σ(x)
    return ( x > 1.0 ) ? 0.0 : ( (x < 0) ? 1.0 : ( 2.0 * x^3 - 3.0 * x^2 + 1.0  ) )
end

function dσ(x)
    return ( x > 1.0 ) ? 0.0 : ( (x < 0) ? 0.0 : ( 6.0 * x^2 - 6.0 * x ) )
end

function genDiffFuncs(
    Δ     :: Float64,
    K_iso :: Float64,
    K_cva :: Float64,
)

    ΔK = K_cva - K_iso
    K_func  = (x,) -> ΔK *  σ(x/Δ + 1.0) + K_iso
    dK_func = (x,) -> ΔK * dσ(x/Δ + 1.0) / Δ
     
    return K_func, dK_func
     
end

function calOp_vdiff(
    vds   :: VerticalDiffusionSolver,
    input :: Union{AbstractArray{Float64}, Nothing} = nothing;
    cal_Jacobian :: Bool = false,
)
#    println("size: ", size(input))
    
    J = nothing

    if input == nothing # no convective_adjustment
        
        for i = 1:vds.sop.op.W_pts
            vds.wksp.W_K_W[i, i] = vds.K_iso
        end

    else  # with convective adjustment

        input_flat = view(input, :)
        dbdz = vds.sop.W_∂z_T * input_flat
        for i = 1:vds.sop.op.W_pts
            vds.wksp.W_K_W[i, i] = vds.K_func(dbdz[i])
        end

#=
        above_W = vds.wksp.above_W
        below_W = vds.wksp.below_W

        input_flat = view(input, :)

        mul!(vds.wksp.above_W, vds.sop.op.W_DN_T, input_flat)
        mul!(vds.wksp.below_W, vds.sop.op.W_UP_T, input_flat)

        @. vds.wksp.Δ_W = vds.wksp.above_W - vds.wksp.below_W

        for i = 1:vds.sop.op.W_pts
            vds.wksp.W_K_W[i, i] = (vds.wksp.Δ_W[i] >= 0.0 ) ? vds.K_iso : vds.K_cva
        end
=#

        

        # make sure boundary does not leak
        vds.wksp.W_K_W .= vds.sop.W_mask_W * vds.wksp.W_K_W
        op_vdiff = sparse(vds.dilution_matrix * vds.sop.T_DIVz_W * vds.wksp.W_K_W * vds.sop.W_∂z_T)

        if cal_Jacobian

            # compute the second term
            
            # compute ∇b (already done above)
            # compute f'(∇b)
            dKdb = spdiagm(0 => dbdz .* vds.dK_func.(dbdz) )

            J = op_vdiff + sparse( vds.dilution_matrix * vds.sop.T_DIVz_W * dKdb * vds.sop.W_∂z_T )
            
        end
    end

    return op_vdiff, J

end

