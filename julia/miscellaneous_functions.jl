
function get∫b(m :: Model)
    return sum(m.core.sop_bib.T_Δv_T * m.state._b)
end

function updateB!(m :: Model)
    m.state.b .= TS2b.(m.state.T, m.state.S)
end

function nblockdiag(a, n :: Integer)
    return blockdiag(ntuple((i,)-> a, n)...) 
end

@inline function TS2b(T::Float64, S::Float64)
    ΔT = T - T_ref
    ΔS = S - S_ref
    return g * ( α*ΔT - β*ΔS )
end
