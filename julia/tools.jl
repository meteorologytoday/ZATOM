using SparseArrays


function speye( 
    dtype::DataType = Float64,
    n :: Integer
)
    return spdiagm(0 => ones(dtype, n))
end
