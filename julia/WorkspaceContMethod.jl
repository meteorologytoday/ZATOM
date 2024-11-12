mutable struct Workspace

    Ns   :: Int64
    ptr  :: Int64
    arrs :: Array

    function Workspace(Ns :: Int64)

        return new(
            Ns, 1, [],
        )

    end
    

end

function getSpace!(
    wksp :: Workspace,
)

    if wksp.ptr > length(wksp.arrs)
        println("Running out of workspace, create new...")
        push!(wksp.arrs, zeros(Float64, wksp.Ns))
    end
   
    result = wksp.arrs[wksp.ptr] 
    wksp.ptr += 1

    return result
end

function reset!(
    wksp :: Workspace,
)
    wksp.ptr = 1
end

function genEmptyGrid(
    wksp  :: Workspace,
)
    return zeros(dtype, dim...)
end
