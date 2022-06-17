lesser(x::TimeOrderedGreenFunction) = x.L

function lesser(x::TimeOrderedConvolution{TimeOrderedGreenFunction, TimeOrderedGreenFunction})
    a, b = x.A, x.B
    a.R * b.L + a.L * b.A
end

function lesser(x::TimeOrderedConvolution{TimeOrderedGreenFunction, TimeOrderedConvolution{TimeOrderedGreenFunction,TimeOrderedGreenFunction}})
    a, b, c = x.A, x.B.A, x.B.B
    lesser3(a, b, c)
end

function lesser(x::TimeOrderedConvolution{TimeOrderedConvolution{TimeOrderedGreenFunction,TimeOrderedGreenFunction}, TimeOrderedGreenFunction})
    a, b, c = x.A.A, x.A.B, x.B
    lesser3(a, b, c)
end

function lesser3(a, b, c)
    x = a.L * b.A * c.A
    x - x' + a.R * b.L * c.A
end

greater(x::TimeOrderedGreenFunction) = x.G

function greater(x::TimeOrderedConvolution{TimeOrderedGreenFunction, TimeOrderedGreenFunction})
    a, b = x.A, x.B
    a.R * b.G + a.G * b.A
end

function greater(x::TimeOrderedConvolution{TimeOrderedGreenFunction, TimeOrderedConvolution{TimeOrderedGreenFunction,TimeOrderedGreenFunction}})
    a, b, c = x.A, x.B.A, x.B.B
    greater3(a, b, c)
end

function greater(x::TimeOrderedConvolution{TimeOrderedConvolution{TimeOrderedGreenFunction,TimeOrderedGreenFunction}, TimeOrderedGreenFunction})
    a, b, c = x.A.A, x.A.B, x.B
    greater3(a, b, c)
end

function greater3(a, b, c)
    x = a.G * b.A * c.A
    x - x' + a.R * b.G * c.A
end
