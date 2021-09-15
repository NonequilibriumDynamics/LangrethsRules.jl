greater(g::TimeOrderedGF) = g.G
lesser(g::TimeOrderedGF) = g.L
retarded(g::TimeOrderedGF) = g.R
advanced(g::TimeOrderedGF) = g.A

greater(g::TimeOrderedConvolution) = greater(g.A) * g.dts * advanced(g.B) + retarded(g.A) * g.dts * greater(g.B)
lesser(g::TimeOrderedConvolution) = lesser(g.A) * g.dts * advanced(g.B) + retarded(g.A) * g.dts * lesser(g.B)
retarded(g::TimeOrderedConvolution) = retarded(g.A) * g.dts * retarded(g.B)
advanced(g::TimeOrderedConvolution) = advanced(g.A) * g.dts * advanced(g.B)
