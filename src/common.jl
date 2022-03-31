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

Base.:*(a::Number, g::TimeOrderedGF) = TimeOrderedGF(a * lesser(g), a * greater(g), g.ts)
Base.:+(g1::TimeOrderedGF, g2::TimeOrderedGF) = TimeOrderedGF(lesser(g1) + lesser(g2), greater(g1) + greater(g2), g1.ts)

"""
"""
struct TimeOrderedConvolution <: AbstractGF
    A::AbstractGF # Left GF
    B::AbstractGF # Right GF
    dts::AbstractMatrix # Trapezoidal rule is achieved by a multiplication with this matrix
end

"""
"""
function conv(A::TimeOrderedGF, B::TimeOrderedGF; 
    dts = reduce(hcat, ([calculate_weights(A.ts[1:i], ones(Int64, i-1)); zeros(length(A.ts)-i)] for i in eachindex(A.ts)))
    )
    c = TimeOrderedConvolution(A, B, dts)
    return TimeOrderedGF(lesser(c), greater(c), B.ts) # no need for lazy `TimeOrderedConvolution`
end
const ⋆ = conv # left-associative operator
