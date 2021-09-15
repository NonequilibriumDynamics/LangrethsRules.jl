greater(g::TimeOrderedGF) = g.G
lesser(g::TimeOrderedGF) = g.L
retarded(g::TimeOrderedGF) = [(g.G[i,j] - g.L[i,j]) * (i > j ? 1.0 : i == j ? 0.5 : 0.0) for i in 1:size(g.G, 1), j in 1:size(g.G, 2)]
advanced(g::TimeOrderedGF) = [(g.L[i,j] - g.G[i,j]) * (i > j ? 1.0 : i == j ? 0.5 : 0.0) for i in 1:size(g.G, 1), j in 1:size(g.G, 2)]

greater(g::TimeOrderedConvolution) = greater(g.A) * g.dts * advanced(g.B) + retarded(g.A) * g.dts * greater(g.B)
lesser(g::TimeOrderedConvolution) = lesser(g.A) * g.dts * advanced(g.B) + retarded(g.A) * g.dts * lesser(g.B)
retarded(g::TimeOrderedConvolution) = retarded(g.A) * g.dts * retarded(g.B)
advanced(g::TimeOrderedConvolution) = advanced(g.A) * g.dts * advanced(g.B)
