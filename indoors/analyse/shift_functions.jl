# function combine_factors(light, induced, run)
#     induced = induced ? " induced" : ""
#     run = run > 1 ? " $run" : ""
#     string(light, induced, run)
# end
#
# function convert_dance_by_to_binary(dance_by)
#     dance_by == "disrupt" && return true
#     dance_by == "no" && return false
#     error("third dance_by option: $dance_by")
# end


function guess_logistic(x, y)
    ymin, ymax = extrema(y)
    x₀i = last(findmax(log.(1 ./ abs.(y .- (ymax - ymin)/2 .- ymin))))
    p0 = Float64[ymax - ymin, -1^(mean(y) > y[1]), x[x₀i], -ymin]
end

fit_logistic(θ) = fit_logistic(parent(θ.l), parent(θ.θ))

function fit_logistic(x, θ)
    p0 = guess_logistic(x, θ)
    lb = Float64[0, -10, x[1], -20]
    ub = Float64[20, 10, x[end], 20]
    fit = curve_fit(logistic, x, θ, p0, lower = lb, upper = ub)
    tss = sum(abs2, θ .- mean(θ))
    (; ks = LsqFit.coef(fit), logistic_rsquare = 1 - LsqFit.rss(fit)/tss)
end

logistic(x, L, k, x₀, y₀) = L / (1 + exp(-k*(x - x₀))) - y₀
logistic(x, p) = logistic.(x, p...)

function get_turn_profile(xy, poi)
    # xy = runs.xy[1]
    # poi = runs.poi[1]
    Δ = diff(xy)
    θ = DimStack((; l = cumsum(norm.(Δ)), θ = unwrap!(splat(atan).(reverse.(Δ)))))
    lpoi = θ[Ti = Near(poi)].l
    i1 = findfirst(≥(lpoi - 5), θ.l)
    i2 = @something findfirst(≥(lpoi + 15), θ.l) length(θ)
    return θ[i1:i2]
end

# function arclength(spl, t1, t2; kws...)
#     knots = get_knots(spl)
#     filter!(t -> t1 < t < t2, knots)
#     pushfirst!(knots, t1)
#     push!(knots, t2)
#     s, _ = quadgk(t -> norm(derivative(spl, t)), knots; kws...)
#     return s
# end

# function get_pathlength(t, spl)
#     n = length(t)
#     l = zeros(n)
#     for i in 2:n
#         l[i] = arclength(spl, t[i-1], t[i])
#     end
#     return cumsum(l)
# end

# function get_pathlength(xy)
#     ll = cumsum(norm.(diff(xy)))
#     pushfirst!(ll, 0.0)
#     return ll
# end

function unwrap!(x, period = 2π)
    y = convert(eltype(x), period)
    v = first(x)
    for k = eachindex(x)
        x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
    end
    return x
end


wrap2pi(x::typeof(π)) = rem2pi(float(x), RoundNearest)
wrap2pi(x) = rem2pi(x, RoundNearest)
