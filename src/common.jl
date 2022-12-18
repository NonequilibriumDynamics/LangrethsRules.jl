"""
"""
abstract type AbstractGreenFunction end

"""
"""
struct TimeOrderedGreenFunction <: AbstractGreenFunction
    L   # Lesser
    G   # Greater
    A   # Advanced
    
    TimeOrderedGreenFunction(L, G) = new(L, G, UpperTriangular(L - G))
end

Base.:*(a::Number, g::TimeOrderedGreenFunction) = TimeOrderedGreenFunction(a * lesser(g), a * greater(g))
Base.:+(g1::TimeOrderedGreenFunction, g2::TimeOrderedGreenFunction) = TimeOrderedGreenFunction(lesser(g1) + lesser(g2), greater(g1) + greater(g2))

"""
"""
struct TimeOrderedConvolution{TA<:AbstractGreenFunction, TB<:AbstractGreenFunction} <: AbstractGreenFunction
    A::TA
    B::TB
    dts::AbstractMatrix # Integration is achieved by multiplication with this matrix
end

"""
"""
function conv(A::AbstractGreenFunction, B::AbstractGreenFunction, dts)
    c = TimeOrderedConvolution(A, B, dts)
    return TimeOrderedGreenFunction(lesser(c), greater(c))
end