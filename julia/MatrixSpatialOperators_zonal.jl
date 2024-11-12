mutable struct MatrixSpatialOperators_zonal

    op        :: MatrixOperators
    
    U_mask_U   :: AbstractArray{Float64, 2}

    T_Δx_T    :: AbstractArray{Float64, 2}
    T_invΔx_T :: AbstractArray{Float64, 2}
 
    function MatrixSpatialOperators_zonal(gd :: Grid)

        Ny = gd.Ny
        Nz = gd.Nz
        Nx = 2
        
 
        cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))

        @time op = MatrixOperators(Ny=Ny, Nz=Nz, Nx=Nx)
        
        Δxw_T  = gd.Δx_T[:, :, 1] 
        Δxie_T = gd.Δx_T[:, :, 1] * (gd.Δλ + gd.Δλi) / gd.Δλ

        Δx_T = [ Δxw_T[:] ; Δxie_T[:] ]

        T_Δx_T     = Δx_T        |> cvtDiagOp
        T_invΔx_T  = Δx_T.^(-1)  |> cvtDiagOp


        mask3_u = ones(Float64, Nz, Ny, Nx)
        mask3_u[:, :, 1] .= 0 
        
        U_mask_U = spdiagm(0 => view(mask3_u, :))

        return new(
            op,

            U_mask_U,

            T_Δx_T,
            T_invΔx_T, 
         
        ) 
    end
end
