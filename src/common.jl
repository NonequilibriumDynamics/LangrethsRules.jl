abstract type AbstractGreenFunction end

"""
    A time-ordered Green function.

    This framework arises because the Heaviside functions that work very nicely 
    on paper are different beasts when the Green function is discretized in a matrix form.
"""
struct TimeOrderedGreenFunction <: AbstractGreenFunction
    L   # Lesser
    G   # Greater
    R   # Retarded (this is not the true retarded function)
    A   # Advanced (this is not the true advanced function)

    hs

    function TimeOrderedGreenFunction(L, G, hs)
        A = UpperTriangular(G - L) .* hs
        new(L, G, adjoint(A), A, hs)
    end
end

struct TimeOrderedConvolution{T<:AbstractGreenFunction, V<:AbstractGreenFunction} <: AbstractGreenFunction
    A::T
    B::V

    hs
end

function Base.:+(a::AbstractGreenFunction, b::AbstractGreenFunction)
    @assert a.hs === b.hs
    TimeOrderedGreenFunction(lesser(a) + lesser(b), greater(a) + greater(b), a.hs)
end

function Base.:*(a::AbstractGreenFunction, b::AbstractGreenFunction)
    @assert a.hs === b.hs
    TimeOrderedConvolution(a, b, a.hs)
end
