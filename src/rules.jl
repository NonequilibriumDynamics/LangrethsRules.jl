greater(g::TimeOrderedGreenFunction) = g.G
lesser(g::TimeOrderedGreenFunction) = g.L
advanced(g::TimeOrderedGreenFunction) = g.A
retarded(g::TimeOrderedGreenFunction) = adjoint(advanced(g))

greater(g::TimeOrderedConvolution) = (retarded(g.A) .* LowerTriangular(transpose(g.dts))) * greater(g.B) + greater(g.A) * (UpperTriangular(g.dts) .* advanced(g.B)) |> skew_hermitify!
lesser(g::TimeOrderedConvolution) = (retarded(g.A) .* LowerTriangular(transpose(g.dts))) * lesser(g.B) + lesser(g.A) * (UpperTriangular(g.dts) .* advanced(g.B)) |> skew_hermitify!