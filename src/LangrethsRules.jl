module LangrethsRules

using LinearAlgebra
import KadanoffBaym: calculate_weights

export TimeOrderedGF, TimeOrderedConvolution, ⋆
export greater, lesser, advanced, retarded

function skew_hermitify!(x)
    for i in 1:size(x, 1)
        for j in 1:(i - 1)
            x[j,i] = -conj(x[i,j])
        end

        x[i,i] = eltype(x) <: Real ? 0.0 : im * imag(x[i,i])
    end
    x
end

include("common.jl")
include("rules.jl")

end # module
