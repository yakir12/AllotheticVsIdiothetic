using LinearAlgebra
using GLMakie
using Optim, StaticArrays

xy(t) = SVector(t^2 - 1, t^3 - t)

t = range(-2, 2, 100)
lines(xy.(t), axis = (; aspect = DataAspect()))

min_step = 1
min_x = -2
max_x = 2
lx = Float64[min_x, min_step]
ux = Float64[max_x - min_step, max_x - min_x]
con_c!(c, x) = (c[1] = sum(x); c)
lc = [min_x + min_step]
uc = [max_x]
dfc = TwiceDifferentiableConstraints(con_c!, lx, ux, lc, uc, :forward)

x0 = [(min_x + max_x - min_step)/2, (max_x - (min_x + max_x - min_step)/2)/2] .+ randn(2)/10

fun(x) = norm(xy(x[1]) - xy(sum(x)))
df = TwiceDifferentiable(fun, x0; autodiff=:forward)

res = optimize(df, dfc, x0, IPNewton())

t1 = range(lx[1], ux[1], 1000)
step = range(lx[2], ux[2], 1000)
z = [lc[] < t + s < uc[] ? fun([t, s]) : NaN for t in t1, s in step]
fig = heatmap(t1, step, 1 ./ z, colorscale=log, axis = (; xlabel = "t1", ylabel = "step"))
scatter!(x0..., color = :green)
scatter!(Optim.minimizer(res)..., color = :red)
display(fig)

