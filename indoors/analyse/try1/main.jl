using Distributions, StaticArrays, GLMakie, AlgebraOfGraphics
using Turing
const SV = SVector{2, Float64}

d = Truncated(Normal(0.3, 1), -π, π)
n = 1000

xs = Vector{SV}(undef, n)
xs[1] = zero(SV)
for i in 2:n
    xs[i] = xs[i-1] + SV(reverse(sincos(rand(d))))
end
lines(xs, axis = (; aspect = DataAspect()))

@model function bmodel(xs)
    μ ~ VonMises(0, 10)
    κ ~ Truncated(InverseGamma(2, 3), 1e-3, Inf)
    α ~ Truncated(Normal(μ, κ), -π, π)
    n = length(xs)
    xys = Vector{SV}(undef, n)
    xys[1] = zero(SV)
    for i in 2:n
        xys[i] = xys[i-1] + SV(reverse(sincos(α)))
    end
    xs .= xys
    return xs
end

chain = sample(bmodel(xs), NUTS(), 10_000, progress=false)

chain = sample(bmodel(xs), NUTS(), MCMCThreads(), 10_000, 5, progress=true)

fig = Figure()
for (i, var_name) in enumerate(chain.name_map.parameters)
    draw!(
        fig[i, 1],
        data(chain) *
        mapping(var_name; color=:chain => nonnumeric) *
        AlgebraOfGraphics.density() *
        visual(fillalpha=0)
    )
end





n = 100
xys = accumulate(1:n; init = zero(SV)) do xy, _
    xy + SV(reverse(sincos(rand(d))))
end
lines(xys, axis = (; aspect = DataAspect()))



xy = Vector{SV}(undef, n)
xy[1] = zero(SV)
for i in 2:n
    xy[i] = xy[i-1] + SV(reverse(sincos(rand(d))))
end

lines(xy)
