using LinearAlgebra
using GLMakie
using Optim, StaticArrays, Dierckx

function detect_self_intersection(spl, t1, t2)
    min_step = 1
    lx = Float64[t1, min_step]
    ux = Float64[t2 - min_step, t2 - t1]
    con_c!(c, x) = sum!(c, x)
    lc = [t1 + min_step]
    uc = [t2]
    dfc = TwiceDifferentiableConstraints(con_c!, lx, ux, lc, uc, :forward)

    x0 = (2lx .+ ux) ./ 3

    # x0 = [(t1 + t2 - min_step)/2, (t2 - (t1 + t2 - min_step)/2)/2]

    fun(x) = norm(spl(x[1]) - spl(sum(x)))
    df = TwiceDifferentiable(fun, x0; autodiff=:forward)

    res = optimize(df, dfc, x0, IPNewton())
    t, step = Optim.minimizer(res)
    return (t1 = t, t2 = t + step)
end

function unwrap!(x, period = 2π)
    y = convert(eltype(x), period)
    v = first(x)
    for k = eachindex(x)
        x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
    end
    return x
end

spl(t) = SVector(t^2 - 1, t^3 - t)
t1 = -2
t2 = 2

t = range(t1, t2, 100)
xy = spl.(t)
lines(xy, axis = (; aspect = DataAspect()))
ti1, ti2 = detect_self_intersection(spl, t1, t2)
scatter!(spl(ti1), color = :red, markersize = 20)
scatter!(spl(ti2), color = :green, marker = :+)


t = filter(x -> !(ti1 < x < ti2), t)
push!(t, ti1)
sort!(t)
spl2 = ParametricSpline(t, stack(spl.(t)), k = 2)
# t = range(t1, t2, 1000)
xy = SVector{2, Float64}.(spl2.(t))
lines(xy, axis = (; aspect = DataAspect()))

θ = [atan(reverse(derivative(spl2, ti))...) for ti in t]
unwrap!(θ)
lines(rad2deg.(θ))

using Interpolations

n = 100
t = range(t1, t2, n)
inter = spl(ti1)
A = Vector{SVector{2, Float64}}(undef, n)
for (i, t) in enumerate(t)
    if ti1 < t < ti2
        A[i] = inter
    else
        A[i] = spl(t)
    end
end
lines(A, axis = (; aspect = DataAspect()))

itp = Interpolations.scale(interpolate(stack(A)', (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)

t = range(t1, t2, 100)
xy = Point2f.(itp.(t, 1), itp.(t, 2))
lines(xy, axis = (; aspect = DataAspect()))

θ = [atan(only(Interpolations.gradient(itp, ti, 2)), only(Interpolations.gradient(itp, ti, 1))) for ti in t if !(ti1 < ti < ti2)]
unwrap!(θ)
lines(rad2deg.(θ))


n = 100
t = collect(range(t1, t2, n))
push!(t, ti1)
sort!(t)
xy = [spl(t) for t in t if !(ti1 < t < ti2)]
itp = interpolate(stack(xy)', (BSpline(Cubic(Natural(OnGrid()))), NoInterp()))

t = range(1, length(xy), 1000)
xy = Point2f.(itp.(t, 1), itp.(t, 2))
lines(xy, axis = (; aspect = DataAspect()))

θ = [atan(only(Interpolations.gradient(itp, ti, 2)), only(Interpolations.gradient(itp, ti, 1))) for ti in t if !(ti1 < ti < ti2)]
unwrap!(θ)
lines(rad2deg.(θ))
