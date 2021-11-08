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
    dts′::AbstractMatrix # And a correction has to be applied (yay theta functions)
end

"""
"""
function conv(A::TimeOrderedGF, B::TimeOrderedGF)
    dts_ = B.ts[2:end] - B.ts[1:end-1]
    dts = 0.5 * Diagonal([dts_; 0.0] + [0.0; dts_])
    dts′ = 0.5 * Diagonal([dts_; 0.0])
    c = TimeOrderedConvolution(A, B, dts, dts′)
    return TimeOrderedGF(lesser(c), greater(c), B.ts) # no need for lazy `TimeOrderedConvolution`
end
const ⋆ = conv # left-associative operator
