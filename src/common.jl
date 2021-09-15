"""
"""
abstract type AbstractGF end

"""
"""
struct TimeOrderedGF <: AbstractGF
    L   # Lesser
    G   # Greater
    R   # Retarded
    A   # Advanced
    ts  # time-grid
end

function TimeOrderedGF(L::AbstractMatrix{T}, G::AbstractMatrix{T}, ts) where {T}
    R = zero(L)
    A = zero(L)
    
    for i in 1:size(L, 1)
        for j in 1:size(L, 2)
            R[i,j] = i > j ? G[i,j] - L[i,j] : i == j ? 0.5 * (G[i,j] - L[i,j]) : zero(T)
            A[i,j] = j > i ? L[i,j] - G[i,j] : i == j ? 0.5 * (L[i,j] - G[i,j]) : zero(T)
        end
    end
    
    return TimeOrderedGF(L, G, R, A, ts)
end

function Base.:*(a::Number, g::TimeOrderedGF)
    return TimeOrderedGF(a * lesser(g), a * greater(g), g.ts)
end

function Base.:+(g1::TimeOrderedGF, g2::AbstractGF)
    return TimeOrderedGF(lesser(g1) + lesser(g2), greater(g1) + greater(g2), g1.ts)
end

"""
"""
struct TimeOrderedConvolution <: AbstractGF
    A::AbstractGF
    B::AbstractGF
    dts::AbstractMatrix
end

"""
"""
function conv(A::AbstractGF, B::AbstractGF)
    dts = B.ts[2:end] - B.ts[1:end-1]
    dts = 0.5 * Diagonal([dts; 0.0] + [0.0; dts])
    return TimeOrderedConvolution(A, B, dts)
end
conv(A::AbstractGF, B::AbstractGF...) = conv(A, conv(B...))
const ⋆ = conv
