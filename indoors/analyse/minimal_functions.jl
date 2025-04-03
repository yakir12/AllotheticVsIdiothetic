const SV = SVector{2, Float64}

function trectify(fs, xs)
    n = length(xs)
    ys = Vector{Vector{SVector{2, Float64}}}.(undef, n)
    Threads.@threads for i in 1:n
        y = fs[i].(xs[i])
        y .-= Ref(y[1])
        ys[i] = y
    end
    return ys
end



function get_calibration(calibration_id)
    c = CameraCalibrations.load(joinpath(results_dir, calibration_id))
    f = rectification(c)
    # fun(ij)::Function = f.(ij)
    return f
end

tosecond(t::T) where {T <: TimePeriod} = t / convert(T, Dates.Second(1))
tosecond(t::TimeType) = tosecond(t - Time(0))
tosecond(sec::Real) = sec

# function get_txy(tij_file, rectify)
#     tij = CSV.File(joinpath(results_dir, tij_file))
#     t = range(tij.t[1], tij.t[end], length = length(tij))
#     ij = SVector{2, Int}.(tij.i, tij.j)
#     xy = rectify.(ij)
#     (; t, xy)
# end

function get_tij(file)
    tij = CSV.File(joinpath(results_dir, file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    (; t = collect(t), ij = tuple.(tij.i, tij.j))
end

function remove_stops!(t, ij)
    tokill = map(==, ij[1:end - 1], ij[2:end])
    pushfirst!(tokill, false)
    deleteat!(t, tokill)
    deleteat!(ij, tokill)
    return (; t, ij)
end

# function remove_stops!(t, xy)
#     tokill = Int[]
#     last_xy = xy[1]
#     for i in 2:length(t)
#         Δ = norm(xy[i] - last_xy)
#         if Δ > 0.1
#             last_xy = xy[i]
#         else
#             push!(tokill, i)
#         end
#     end
#     deleteat!(t, tokill)
#     deleteat!(xy, tokill)
#     return (; t, xy)
# end
#
function remove_loops!(t, xy)
    inds, _ = self_intersections(Point2f.(xy))
    ranges = splat(UnitRange).(Iterators.partition(inds, 2))
    filter!(<(51) ∘ length, ranges)
    tokill = vcat(ranges...)
    unique!(tokill)
    sort!(tokill)
    deleteat!(t, tokill)
    deleteat!(xy, tokill)
    return (; t, xy)
end

function sparseify(t, xy, poi)
    spl = ParametricSpline(t, stack(xy))
    tl = round(Int, t[1]):round(Int, t[end])
    poi_index = findfirst(==(round(Int, poi)), tl)
    (; t = tl, xy = SV.(spl.(tl)), poi_index)
end

function smooth(t, xy)
    spl = ParametricSpline(t, stack(xy), k = 3, s = 50)
    xys = SV.(spl.(t))
    trans = Translation(-xys[1])
    xys .= trans.(xys)
    # t1 = t[1]
    # t2 = t[end]
    # xys1 = xys[1]
    # xys2 = xys[end]
    # remove_loops!(t, xys)
    # if t[1] ≠ t1
    #     pushfirst!(t, t1)
    #     pushfirst!(xys, xys1)
    # end
    # if t[end] ≠ t2
    #     push!(t, t2)
    #     push!(xys, xys2)
    # end
    # xys .-= Ref(xys[1])
    return (; xys, spl, trans)
end



# function remove_cycles(t, ij)
#     t, ij = interpolate_pixels(t, ij)
#     ids = levelsmap(ij)
#     vs = [ids[k] for k in ij]
#     g = DiGraph(length(vs))
#     for (v1, v2) in zip(vs[1:end-1], vs[2:end])
#         add_edge!(g, v1, v2)
#     end
#     c = simplecycles(g)
#     if isempty(c)
#         return (t, ij)
#     end
#     tokill = sort(unique(reduce(vcat, c)))
#     deleteat!(ij, tokill)
#     deleteat!(t, tokill)
#     return (t, ij)
# end
#
# function interpolate_pixels(t, ij)
#     n = length(t)
#     ijs = Vector{Vector{Tuple{Int, Int}}}(undef, n - 1)
#     ts = Vector{Vector{Float64}}(undef, n - 1)
#     for (i, (t1, t2, (y1, x1), (y2, x2))) in enumerate(zip(t[1:end-1], t[2:end], ij[1:end-1], ij[2:end]))
#         indices = bresenham(y1, x1, y2, x2)
#         ijs[i] = indices
#         n = length(indices)
#         if n > 1
#             ts[i] = collect(range(t1, t2, n))
#         else
#             ts[i] = [t1]
#         end
#     end
#     ij = reduce(vcat, ijs)
#     t = reduce(vcat, ts)
#     tokill = findall(==(0) ∘ norm, diff(SV.(ij))) .+ 1
#     deleteat!(ij, tokill)
#     deleteat!(t, tokill)
#     return (t, ij)
# end
#
# function bresenham(y0::Int, x0::Int, y1::Int, x1::Int)
#     dx = abs(x1 - x0)
#     dy = abs(y1 - y0)
#     sx = x0 < x1 ? 1 : -1
#     sy = y0 < y1 ? 1 : -1;
#     err = (dx > dy ? dx : -dy) / 2
#     indices = Vector{Tuple{Int, Int}}(undef, 0)
#     while true
#         push!(indices, (y0, x0))
#         (x0 != x1 || y0 != y1) || break
#         e2 = err
#         if e2 > -dx
#             err -= dy
#             x0 += sx
#         end
#         if e2 < dy
#             err += dx
#             y0 += sy
#         end
#     end
#     return indices
# end
#
# function clean_coords(t, ij)
#     for i in 1:100
#         t, ij = remove_cycles(t, ij)
#     end
#     tl = range(t[1], t[end], step = 1/5)
#     tp = ParametricSpline(t, stack(ij))
#     return (; t = tl, ij = [Tuple(round.(Int, tp(t))) for t in tl])
# end

# function clean_coords!(xy)
#     for i in 2:length(xy)
#         Δ = norm(xy[i] - xy[i - 1])
#         if Δ < 1
#             xy[i] = xy[i - 1]
#         end
#     end
#     return xy
# end



# function prune_coords!(t, xy)
#     rmax = 10 # cm
#     rmin = 1 # cm
#     for i in 1:length(xy) - 1
#         j = copy(i + 1)
#         while norm(xy[j] - xy[i]) < rmax
#             j += 1
#
#
#
#     Δ = round(Int, 1/step(t)) # Δ steps for 1 second
#     xy = xy[1:Δ:end]
#     t = range(t[1], t[end], length(xy))
#     return (; t, xy)
# end

function impute_poi_time(t, xy)
    poi_index = something(findfirst(>(10) ∘ norm, xy .- Ref(xy[1])), length(xy))
    t[poi_index]
end

function glue_intervention!(xy, t, intervention)
    inter_i = something(findfirst(≥(intervention), t), length(t))
    h = 2
    diffs = diff(xy[1:inter_i + h + 1])
    Δs = norm.(diffs)
    μ = mean(Δs[1:inter_i - h])
    σ = std(Δs[1:inter_i - h], mean = μ)
    for i in inter_i - h:inter_i + h
        if Δs[i] > μ + 1.5σ
            xy[i + 1:end] .-= Ref(diffs[i] - μ*normalize(diffs[i]))
            return Δs[i]
        end
    end
    return missing
end

# function glue_poi_index!(xy, t, poi::Float64)
#     poi_index = something(findfirst(≥(poi), t), length(t))
#     Δ = gluePOI!(xy, poi_index)
#     (; poi_index, dance_jump = Δ)
# end

# function get_spline(xy, t)
#     k, s = (3, 300)
#     k, s = (3, 100)
#     tp = ParametricSpline(t, stack(xy); k, s)
# end
#
# function smooth_center(t, spl)
#     xy = SV.(spl.(t))
#     return xy .- Ref(xy[1])
# end

function get_rotation(p2)
    θ = π/2 - atan(reverse(p2)...)
    LinearMap(Angle2d(θ))
end

function get_smooth_center_poi_rotate(t, spl, poi_index)
    f = SV ∘ spl
    p1 = f(t[1])
    p2 = f(t[poi_index])
    trans = Translation(-p2)
    p1 = trans(p1)
    rot = get_rotation(p1)
    return rot ∘ trans ∘ f
end

function get_exit_angle(xyp, r)
    i = findfirst(≥(r) ∘ norm, xyp)
    if isnothing(i)
        return missing
    end
    x, y = xyp[i]
    atan(y, x)
end

function mean_resultant_vector(θ)
    norm(mean(SV ∘ reverse ∘ sincos, θ))
end


function get_turn_profile(t, spl, poi)
    n = 1000
    tl = range(t[1], t[end], n)
    poi_index = findfirst(≥(poi), tl)
    l = get_pathlength(SV.(spl.(tl)))
    l1 = l[poi_index] - 10 # 5 cm beofe the poi
    i1 = findfirst(≥(l1), l)
    l2 = l[poi_index] + 10 # 5 cm beofe the poi
    i2 = something(findfirst(≥(l2), l), n)
    θ = [atan(reverse(derivative(spl, tl[i]))...) for i in i1:i2]
    unwrap!(θ)
    return (; lpoi = l[poi_index], lθ = l[i1:i2], θ)
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

function get_pathlength(xys)
    l = cumsum(norm.(diff(xys)))
    pushfirst!(l, 0.0)
    return l
end

function unwrap!(x, period = 2π)
    y = convert(eltype(x), period)
    v = first(x)
    for k = eachindex(x)
        x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
    end
    return x
end

function cropto(t, spl, trans, l)
    i = findlast(<(l) ∘ norm ∘ trans ∘ SV ∘ spl, t)
    t1 = t[i]
    t2 = t[i+1]
    fun(t) = abs2(norm(trans(SV(spl(t)))) - l)
    res = optimize(fun, t1, t2)
    tend = Optim.minimizer(res)
    tl = range(t[1], tend, 100)
    trans.(SV.(spl.(tl)))
end

function guess_logistic(x, y)
    ymin, ymax = extrema(y)
    x₀i = last(findmax(log.(1 ./ abs.(y .- (ymax - ymin)/2 .- ymin))))
    p0 = Float64[ymax - ymin, -1^(mean(y) > y[1]), x[x₀i], -ymin]
end

function fit_logistic(l, θ)
    p0 = guess_logistic(l, θ)
    # @show p0
    lb = Float64[0, -10, l[1], 0]
    ub = Float64[20, 10, l[end], 10]
    fit = curve_fit(logistic, l, θ, p0, lower = lb, upper = ub)
    tss = sum(abs2, θ .- mean(θ))
    (; ks = LsqFit.coef(fit), logistic_rsquare = 1 - LsqFit.rss(fit)/tss)
end

logistic(x, L, k, x₀, y₀) = L / (1 + exp(-k*(x - x₀))) - y₀
logistic(x, p) = logistic.(x, p...)

inv_logistic(y, L, k, x₀, y₀) = log(L/(y + y₀) - 1) / (-k) + x₀

norm_logistic(y, L, k, y₀) = sign(k)*((y + y₀) / L - 0.5) + 0.5

# """
#     wrap2pi(angle)
#
# Limit the angle to the range -π .. π .
# """
# wrap2pi(::typeof(pi)) = π
# function wrap2pi(angle)
#     y = rem(angle, 2π)
#     abs(y) > π && (y -= 2π * sign(y))
#     return y
# end
wrap2pi(x::typeof(π)) = rem2pi(float(x), RoundNearest)
wrap2pi(x) = rem2pi(x, RoundNearest)
#######################





# function totuple(x::AbstractString)
#     if contains(x, '(')
#         m = match(r"^\((\d+),\s*(\d+)\)$", x)
#         Tuple{Int, Int}(parse.(Int, m.captures))
#     else
#         parse(Int, x)
#     end
# end
# totuple(x) = x
# function calibrate_and_smooth(c, track, s, k)
#     xy_pixels = RowCol.(track.x, track.y)
#     xy_cm = CameraCalibrationMeta.rectification(c).(xy_pixels)
#     XY = reduce(hcat, xy_cm)
#     spl = ParametricSpline(track.t, XY; s, k)
#     return SV ∘ spl ∘ tosecond
# end
#
#
# function track_function(start, POI, c, track; s = 100, k = 2)
#     spl = calibrate_and_smooth(c, track, s, k)
#     cr = center_and_rotate(start, POI, spl)
#     return cr ∘ spl
# end
#
#
# cordlength(rotated) = cordlength(rotated.xy[rotated.i:end])
# cordlength(xy::Vector{SV}) = norm(diff([xy[1], xy[end]]))
#
# curvelength(rotated) = curvelength(rotated.xy[rotated.i:end])
# function curvelength(xy::Vector{SV})
#     p0 = xy[1]
#     s = 0.0
#     for p1 in xy[2:end]
#         s += norm(p1 - p0)
#         p0 = p1
#     end
#     return s
# end
#
# # tortuosity(rotated) = tortuosity(rotated.xy[rotated.i:end])
# # tortuosity(xy::Vector{SV}) = cordlength(xy) / curvelength(xy)
#
#
#
# function angle_between(p1, p2)
#     θ = acos(dot(p1, p2) / norm(p1) / norm(p2))
#     return -sign(cross(p1, p2))*θ
# end
#
# function cumulative_angle(f, t1, t2)
#     xy = f.(range(t1, t2; step = Second(1)))
#     δ = filter(!iszero ∘ sum, diff(xy))
#     α = sum(splat(angle_between) ∘ reverse, partition(δ, 2, 1); init=0.0)
#     rad2deg(α)
# end
#
# function plotone(run_id, xy, poi_index, center_rotate, sxy, t, intervention, spontaneous_end, l, θ, condition)
#     fig  = Figure()
#     ax = Axis(fig[1,1], aspect = DataAspect(), title = string(run_id, " ", condition))
#     # for r  in (30, 50)
#     #     lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
#     # end
#     scatter!(ax, center_rotate.(xy[1:poi_index]), markersize = 2)
#     scatter!(ax, center_rotate.(xy[poi_index:end]), markersize = 2)
#     lines!(ax, sxy[1:poi_index])
#     lines!(ax, sxy[poi_index:end])
#     intervention_i = findfirst(≥(intervention), t)
#     scatter!(ax, sxy[intervention_i], label = "intervention")
#     if !ismissing(spontaneous_end)
#         spontaneous_end_i = findfirst(≥(spontaneous_end), t)
#         scatter!(ax, sxy[spontaneous_end_i], label = "spontaneous dance")
#     end
#     ax = Axis(fig[2,1], ytickformat = "{:n}°", xlabel = "Path length (cm)", ylabel = "θ")
#     lines!(ax, l, rad2deg.(θ))
#     scatter!(ax, l[poi_index], rad2deg(θ[poi_index]))
#     rowsize!(fig.layout, 1, Relative(9/10))
#     return fig
# end
#
#
# # function get_turn_profile(t, spl, poi_index, p)
# #     i = round(Int, poi_index*p)
# #     der = derivative.(Ref(spl), t[1:i])
# #     θ = [atan(reverse(d)...) for d in der]
# #     unwrap!(θ)
# #     θ .-= θ[1]
# #     return rad2deg(abs(θ[end]))
# # end
#
#
# function get_turn_profile(t, spl)
#     θ = [atan(reverse(derivative(spl, ti))...) for ti in t]
#     θ .-= θ[1]
#     unwrap!(θ)
#     if mean(θ) < 0
#         θ .*= -1
#     end
#     return θ
# end
#
# function get_turn_profile(t, spl, poi_index)
#     Δ = round(Int, 1/step(t))
#     tθ = t[poi_index - Δ:end]
#     θ = [atan(reverse(derivative(spl, ti))...) for ti in tθ]
#     θ₀ = θ[1]
#     θ .-= θ₀
#     unwrap!(θ)
#     # m, M = extrema(θ)
#     # if abs(m) > M
#     if mean(θ) < 0
#         θ .*= -1
#     end
#     return (; tθ = tθ .- tθ[1], θ = θ)
# end
#
# function getIQR(x)
#     α = 0.05
#     d = Truncated(Distributions.fit(Normal, x), 0, Inf)
#     c1, μ, c2 = quantile(d, [α/2, 0.5, 1 - α/2])
#     (; c1, μ, c2)
# end
#
# function create_track(k)
#     xy = [zero(Point2f)]
#     t = 0
#     while last(last(xy)) > -50
#         t += 0.1
#         Δ = reverse(sincos(logistic(t, 2π, k, 0) + π/2))
#         push!(xy, xy[end] + Point2f(Δ))
#     end
#     spl = ParametricSpline(range(0, t, length(xy)), stack(xy))
#     xys = Point2f.(spl.(range(0, t, 100)))
#     xys[end] = Point2f(xys[end][1], -50)
#     return xys
# end
#
#
#
# function get_mean_speed(t, spl, poi_index)
#     t1 = t[poi_index]
#     t2 = t[end]
#     s = arclength(spl, t1, t2)
#     return s / (t2 - t1)
# end
#
# # function t2length(t, spl)
# #     n = length(t)
# #     l = Vector{Float64}(undef, n)
# #     l[1] = 0.0
# #     for i in 2:n
# #         l[i] = first(quadgk(t -> norm(derivative(spl, t)), t[1], t[i]))
# #     end
# #     return l
# # end
