lesser(x::TimeOrderedGreenFunction) = x.L
function lesser(x::TimeOrderedConvolution{2})
    (LowerTriangular(x.hs) .* retarded(x.gs[1])) * lesser(x.gs[2]) + lesser(x.gs[1]) * (advanced(x.gs[2]) .* UpperTriangular(x.hs))
end
function lesser(x::TimeOrderedConvolution{3})
    r1 = LowerTriangular(x.hs) .* retarded(x.gs[1])
    a3 = advanced(x.gs[3]) .* UpperTriangular(x.hs)
    l = lesser(x.gs[1]) * (advanced(x.gs[2]) .* UpperTriangular(x.hs)) * a3
    l - l' + r1 * lesser(x.gs[2]) * a3
end

# function greater(x::TimeOrderedConvolution{2})
#     (LowerTriangular(x.hs) .* retarded(x.gs[1])) * greater(x.gs[2]) + greater(x.gs[1]) * (advanced(x.gs[2]) .* UpperTriangular(x.hs))
# end
# function greater(x::TimeOrderedConvolution{3})
#     r1 = LowerTriangular(x.hs) .* retarded(x.gs[1])
#     a3 = advanced(x.gs[3]) .* UpperTriangular(x.hs)
#     g = greater(x.gs[1]) * (advanced(x.gs[2]) .* UpperTriangular(x.hs)) * a3
#     g - g' + r1 * greater(x.gs[2]) * a3
# end
greater(x::TimeOrderedGreenFunction) = x.G
greater(x::TimeOrderedConvolution{1}) = greater(first(x.gs))
function greater(x::TimeOrderedConvolution)
    (LowerTriangular(x.hs) .* retarded(first(x.gs))) * greater(TimeOrderedConvolution(tail(x.gs), x.hs)) + greater(first(x.gs)) * (advanced(TimeOrderedConvolution(tail(x.gs), x.hs)) .* UpperTriangular(x.hs))
end
# function greater(x::TimeOrderedConvolution{3})
#     r1 = LowerTriangular(x.hs) .* retarded(x.gs[1])
#     a3 = advanced(x.gs[3]) .* UpperTriangular(x.hs)
#     g = greater(x.gs[1]) * (advanced(x.gs[2]) .* UpperTriangular(x.hs)) * a3
#     g - g' + r1 * greater(x.gs[2]) * a3
# end

retarded(x::TimeOrderedGreenFunction) = LowerTriangular(lesser(x) - greater(x))
retarded(x::TimeOrderedConvolution{1}) = retarded(first(x.gs))
function retarded(x::TimeOrderedConvolution)
    (LowerTriangular(x.hs) .* retarded(first(x.gs))) * retarded(TimeOrderedConvolution(tail(x.gs), x.hs))
end

advanced(x::TimeOrderedGreenFunction) = UpperTriangular(greater(x) - lesser(x))
advanced(x::TimeOrderedConvolution{1}) = advanced(first(x.gs))
function advanced(x::TimeOrderedConvolution)
    (advanced(first(x.gs)) .* UpperTriangular(x.hs)) * advanced(TimeOrderedConvolution(tail(x.gs), x.hs))
end
