"""
"""
abstract type AbstractGF end

"""
"""
struct TimeOrderedGF <: AbstractGF
    L   # Lesser
    G   # Greater
    R   # Retarded
    ts  # time-grid
    
    TimeOrderedGF(L, G, ts) = new(L, G, LowerTriangular(G - L), ts)
end

"""
"""
struct TimeOrderedConvolution <: AbstractGF
    A::AbstractGF # Left GF
    B::AbstractGF # Right GF
    dts::AbstractMatrix # Trapezoidal rule is achieved by a multiplication with this matrix
end

"""
"""
function conv(A::TimeOrderedGF, B::TimeOrderedGF)
    dts = B.ts[2:end] - B.ts[1:end-1]
    dts = 0.5 * Diagonal([dts; 0.0] + [0.0; dts])
    return TimeOrderedConvolution(A, B, dts)
end
conv(A::TimeOrderedConvolution, B::AbstractGF) = TimeOrderedConvolution(A, B, A.dts)
conv(A::AbstractGF, B::TimeOrderedConvolution) = TimeOrderedConvolution(A, B, B.dts)
const ⋆ = conv # left-associative operator


function Base.:*(a::Number, g::AbstractGF)
    return TimeOrderedGF(a * lesser(g), a * greater(g), g.ts)
end

function Base.:+(g1::TimeOrderedGF, g2::AbstractGF)
    return TimeOrderedGF(lesser(g1) + lesser(g2), greater(g1) + greater(g2), g1.ts)
end
