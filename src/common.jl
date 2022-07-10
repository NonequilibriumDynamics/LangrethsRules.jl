abstract type AbstractGreenFunction end

struct TimeOrderedGreenFunction <: AbstractGreenFunction
    L   # Lesser
    G   # Greater
end

struct TimeOrderedConvolution{N} <: AbstractGreenFunction
    gs::NTuple{N, AbstractGreenFunction}
    hs
end

function Base.:+(a::AbstractGreenFunction, b::AbstractGreenFunction)
    TimeOrderedGreenFunction(lesser(a) + lesser(b), greater(a) + greater(b))
end
