using LangrethsRules
using Test

N = 100

# Suppose a non-equidistant time-grid
ts = sort(N*rand(N));

# And 2 time-ordered GFs defined on that grid
g1 = TimeOrderedGF(rand(N,N), rand(N,N), ts)
g2 = TimeOrderedGF(rand(N,N), rand(N,N), ts)

# Obtaining convolutions is simple (and fast)
result1 = greater(g1 ⋆ g2)
result2 = lesser(g1 ⋆ g2 ⋆ g1) # can also chain convolutions

function integrate(x::AbstractVector, y::AbstractVector)
    if isone(length(x))
        return zero(first(y))
    end
    @inbounds retval = (x[2] - x[1]) * (y[1] + y[2])
    @inbounds @fastmath @simd for i in 2:(length(y) - 1)
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 1//2 * retval
end

function conv_greater_benchmark(g1::TimeOrderedGF, g2::TimeOrderedGF, ts)
    g = zero(g1.L)
    for i in 1:size(g,1)
        for j in 1:size(g,2)
            g[i,j] += integrate(ts, retarded(g1)[i, :] .* greater(g2)[:, j])
            g[i,j] += integrate(ts, greater(g1)[i, :] .* advanced(g2)[:, j])
        end
    end
    g
end

@test conv_greater_benchmark(g1, g2, ts) ≈ greater(g1 ⋆ g2)
@test retarded((g1 ⋆ g2) ⋆ g1) ≈ retarded(g1 ⋆ (g2 ⋆ g1))
