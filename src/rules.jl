greater(g::TimeOrderedGF) = g.G
lesser(g::TimeOrderedGF) = g.L
retarded(g::TimeOrderedGF) = g.R
advanced(g::TimeOrderedGF) = adjoint(retarded(g))

greater(g::TimeOrderedConvolution) = (retarded(g.A) * g.dts - Diagonal(retarded(g.A)) * g.dts′) * greater(g.B) + greater(g.A) * (g.dts * advanced(g.B) - g.dts′ * Diagonal(advanced(g.B))) |> skew_hermitify!
lesser(g::TimeOrderedConvolution) = (retarded(g.A) * g.dts - Diagonal(retarded(g.A)) * g.dts′) * lesser(g.B) + lesser(g.A) * (g.dts * advanced(g.B) - g.dts′ * Diagonal(advanced(g.B))) |> skew_hermitify!
