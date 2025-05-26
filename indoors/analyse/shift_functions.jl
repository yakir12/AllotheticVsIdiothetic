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

function fit_logistic(l, θ)
    p0 = guess_logistic(l, θ)
    # @show p0
    lb = Float64[0, -10, l[1], -20]
    ub = Float64[20, 10, l[end], 20]
    fit = curve_fit(logistic, l, θ, p0, lower = lb, upper = ub)
    tss = sum(abs2, θ .- mean(θ))
    (; ks = LsqFit.coef(fit), logistic_rsquare = 1 - LsqFit.rss(fit)/tss)
end

logistic(x, L, k, x₀, y₀) = L / (1 + exp(-k*(x - x₀))) - y₀
logistic(x, p) = logistic.(x, p...)

function get_turn_profile(xy, poi)

    xy = runs.xy[1]
    poi = runs.poi[1]

    ll = cumsum(norm.(diff(xy)))

    l1 = ll[Ti = Near(poi)] - 5
    i1 = findfirst(≥(l1), ll)

    l2 = ll[Ti = Near(poi)] + 15
    i2 = findfirst(≥(l2), ll)

    θ = [atan(reverse(xy[i])...) for i in i1:i2]
    unwrap!(θ)

    t = lookup(xy, Ti)[i1:i2]

    DimVector(θ, (; Ti = t, Le = ll[i1:i2]))

    DimStack(DimVector(θ, (; Ti = t)), DimVector(θ, (; Ti = t)))







    l2 = ll[Ti = Near(poi)] + 15

    n = 1000
    tl = range(t[1], t[end], n)
    poi_index = findfirst(≥(poi), tl)
    l = get_pathlength(xy)
    l1 = l[poi_index] - 5 # 5 cm beofe the poi
    i1 = findfirst(≥(l1), l)
    l2 = l[poi_index] + 15 # 5 cm beofe the poi
    i2 = something(findfirst(≥(l2), l), n)
    θ = [atan(reverse(derivative(spl, tl[i]))...) for i in i1:i2]
    unwrap!(θ)
    lpoi = l[poi_index]
    lθ = l[i1:i2]
    lpoi_index = findfirst(≥(lpoi), lθ)
    return (; lpoi_index, lpoi, lθ, θ)
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


