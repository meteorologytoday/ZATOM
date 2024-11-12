mutable struct HorizontalDiffusionSolver
    
    sop :: MatrixSpatialOperators

    eT_length    :: Int64
   
    eT_send_T    :: AbstractArray{Float64, 2}
    T_send_eT    :: AbstractArray{Float64, 2}
    eT_I_eT      :: AbstractArray{Float64, 2}

    K            :: Float64

    dilution_factor :: Float64
    dilution_matrix :: AbstractArray{Float64, 2}

    function HorizontalDiffusionSolver(
        sop   :: MatrixSpatialOperators;
        K     :: Float64,
        dilution_factor :: Float64 = 1.0,
    )


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

        if Nx != 2
            throw(ErrorException("VerticalDiffusionSolver.jl only supports Nx = 2"))
        end
        dilution_matrix = zeros(Float64, Nz, Ny, Nx)
        dilution_matrix[:, :, 1] .= 1.0
        dilution_matrix[:, :, 2] .= dilution_factor
        dilution_matrix = cvtDiagOp(dilution_matrix)


        return new(
            sop,
            eT_length,
            eT_send_T,
            T_send_eT,
            eT_I_eT,
            K,
            dilution_factor,
            dilution_matrix,
        )

    end
end

function calOp_hdiff(
    hds :: HorizontalDiffusionSolver,
)

    #return sparse(hds.K * hds.eT_send_T * hds.sop.T_mLAPy_T * hds.T_send_eT)
    return sparse(hds.dilution_matrix * hds.K * hds.eT_send_T * hds.sop.T_LAPy_T * hds.T_send_eT)
end
