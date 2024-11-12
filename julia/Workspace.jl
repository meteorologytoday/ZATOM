mutable struct Workspace

    Ny :: Int64
    Nz :: Int64

    T  :: Array
    V  :: Array
    W  :: Array
    VW :: Array

    ptr :: Dict

    function Workspace(;
        Ny :: Int64,
        Nz :: Int64,
        T :: Int64=0,
        V :: Int64=0,
        W :: Int64=0,
        VW :: Int64=0,
    )

        ptr = Dict(
            :T => 1,
            :V => 1,
            :W => 1,
            :VW => 1,
        )
        return new(
            Ny, Nz,
            [], [], [], [],
            ptr,
        )

    end
    

end

function getSpace!(
    wksp :: Workspace,
    grid :: Symbol;
    flat :: Bool = false
)
    i = wksp.ptr[grid]
    list = getfield(wksp, grid)

    if i > length(list)
        println("Running out of workspace of " * string(grid) * ", create new...")
        push!(list, genEmptyGrid(wksp, Float64, grid))
    end
    
    wksp.ptr[grid] += 1

    return ( flat ) ? view(list[i], :) : list[i]
end

function reset!(
    wksp :: Workspace,
    grid :: Symbol=:ALL,
)
    if grid == :ALL
        for k in keys(wksp.ptr)
            wksp.ptr[k] = 1
        end
    else
        wksp.ptr[grid] = 1
    end

end

function genEmptyGrid(
    wksp  :: Workspace,
    dtype :: DataType,
    grid  :: Symbol,
)
    Ny, Nz = wksp.Ny, wksp.Nz
    dim = Dict(
        :T  =>  [Nz,   Ny  ],
        :V  =>  [Nz,   Ny+1],
        :W  =>  [Nz+1, Ny  ],
        :VW =>  [Nz+1, Ny+1],
    )[grid]

    return zeros(dtype, dim...)

end
