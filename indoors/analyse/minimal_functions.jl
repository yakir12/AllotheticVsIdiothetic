const SV = SVector{2, Float64}

function get_calibration(calibration_id)
    c = CameraCalibrations.load(joinpath(results_dir, calibration_id))
    rectification(c)
end

tosecond(t::T) where {T <: TimePeriod} = t / convert(T, Dates.Second(1))
tosecond(t::TimeType) = tosecond(t - Time(0))
tosecond(sec::Real) = sec

function get_txy(tij_file, rectify)
    tij = CSV.File(joinpath(results_dir, tij_file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    ij = SVector{2, Int}.(tij.i, tij.j)
    xy = rectify.(ij)
    (; t, xy)
end

function clean_coords!(xy)
    for i in 2:length(xy)
        Δ = norm(xy[i] - xy[i - 1])
        if Δ < 0.1
            xy[i] = xy[i - 1]
        end
    end
    return xy
end

impute_poi_time(_, _, poi::Float64) = poi
function impute_poi_time(t, xy, ::Missing)
    poi_index = something(findfirst(>(10) ∘ norm, xy .- Ref(xy[1])), length(xy))
    t[poi_index]
end

function gluePOI!(xy, poi_index)
    diffs = norm.(diff(xy[1:poi_index - 1]))
    μ = mean(diffs)
    σ = std(diffs)
    Δ = only(diff(xy[poi_index - 1:poi_index]))
    l = norm(Δ)
    if l > μ + 1.5σ
        xy[poi_index:end] .-= Ref(Δ)
        return l
    else
        return missing
    end
end

function glue_poi_index!(xy, t, poi::Float64)
    poi_index = something(findfirst(≥(poi), t), length(t))
    Δ = gluePOI!(xy, poi_index)
    (; poi_index, dance_jump = Δ)
end

function get_spline(xy, t)
    k, s = (3, 300)
    tp = ParametricSpline(t, stack(xy); k, s)
end

function smooth_center(t, spl)
    xy = SV.(spl.(t))
    return xy .- Ref(xy[1])
end

function get_rotation(p2)
    θ = -π/2 - atan(reverse(p2)...)
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

function get_exit_angle(xyp, r = 20)
    i = findfirst(≥(r) ∘ norm, xyp)
    if isnothing(i)
        return missing
    end
    atan(reverse(xyp[i])...)
end

function mean_resultant_vector(θ)
    norm(mean(SV ∘ reverse ∘ sincos, θ))
end



#######################





function totuple(x::AbstractString)
    if contains(x, '(')
        m = match(r"^\((\d+),\s*(\d+)\)$", x)
        Tuple{Int, Int}(parse.(Int, m.captures))
    else
        parse(Int, x)
    end
end
totuple(x) = x
function calibrate_and_smooth(c, track, s, k)
    xy_pixels = RowCol.(track.x, track.y)
    xy_cm = CameraCalibrationMeta.rectification(c).(xy_pixels)
    XY = reduce(hcat, xy_cm)
    spl = ParametricSpline(track.t, XY; s, k)
    return SV ∘ spl ∘ tosecond
end


function track_function(start, POI, c, track; s = 100, k = 2)
    spl = calibrate_and_smooth(c, track, s, k)
    cr = center_and_rotate(start, POI, spl)
    return cr ∘ spl
end


cordlength(rotated) = cordlength(rotated.xy[rotated.i:end])
cordlength(xy::Vector{SV}) = norm(diff([xy[1], xy[end]]))

curvelength(rotated) = curvelength(rotated.xy[rotated.i:end])
function curvelength(xy::Vector{SV})
    p0 = xy[1]
    s = 0.0
    for p1 in xy[2:end]
        s += norm(p1 - p0)
        p0 = p1
    end
    return s
end

# tortuosity(rotated) = tortuosity(rotated.xy[rotated.i:end])
# tortuosity(xy::Vector{SV}) = cordlength(xy) / curvelength(xy)



function angle_between(p1, p2)
    θ = acos(dot(p1, p2) / norm(p1) / norm(p2))
    return -sign(cross(p1, p2))*θ
end

function cumulative_angle(f, t1, t2)
    xy = f.(range(t1, t2; step = Second(1)))
    δ = filter(!iszero ∘ sum, diff(xy))
    α = sum(splat(angle_between) ∘ reverse, partition(δ, 2, 1); init=0.0)
    rad2deg(α)
end

function plotone(run_id, xy, poi_index, center_rotate, sxy, t, intervention, spontaneous_end, l, θ, condition)
    fig  = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect(), title = string(run_id, " ", condition))
    # for r  in (30, 50)
    #     lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
    # end
    scatter!(ax, center_rotate.(xy[1:poi_index]), markersize = 2)
    scatter!(ax, center_rotate.(xy[poi_index:end]), markersize = 2)
    lines!(ax, sxy[1:poi_index])
    lines!(ax, sxy[poi_index:end])
    intervention_i = findfirst(≥(intervention), t)
    scatter!(ax, sxy[intervention_i], label = "intervention")
    if !ismissing(spontaneous_end)
        spontaneous_end_i = findfirst(≥(spontaneous_end), t)
        scatter!(ax, sxy[spontaneous_end_i], label = "spontaneous dance")
    end
    ax = Axis(fig[2,1], ytickformat = "{:n}°", xlabel = "Path length (cm)", ylabel = "θ")
    lines!(ax, l, rad2deg.(θ))
    scatter!(ax, l[poi_index], rad2deg(θ[poi_index]))
    rowsize!(fig.layout, 1, Relative(9/10))
    return fig
end


function cropto(xy, l)
    i = something(findfirst(>(l) ∘ norm, xy), length(xy))
    xy[1:i-1]
end

function unwrap!(x, period = 2π)
    y = convert(eltype(x), period)
    v = first(x)
    for k = eachindex(x)
        x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
    end
    return x
end

# function get_turn_profile(t, spl, poi_index, p)
#     i = round(Int, poi_index*p)
#     der = derivative.(Ref(spl), t[1:i])
#     θ = [atan(reverse(d)...) for d in der]
#     unwrap!(θ)
#     θ .-= θ[1]
#     return rad2deg(abs(θ[end]))
# end


function get_turn_profile(t, spl)
    θ = [atan(reverse(derivative(spl, ti))...) for ti in t]
    θ .-= θ[1]
    unwrap!(θ)
    if mean(θ) < 0
        θ .*= -1
    end
    return θ
end

function get_turn_profile(t, spl, poi_index)
    Δ = round(Int, 1/step(t))
    tθ = t[poi_index - Δ:end]
    θ = [atan(reverse(derivative(spl, ti))...) for ti in tθ]
    θ₀ = θ[1]
    θ .-= θ₀
    unwrap!(θ)
    # m, M = extrema(θ)
    # if abs(m) > M
    if mean(θ) < 0
        θ .*= -1
    end
    return (; tθ = tθ .- tθ[1], θ = θ)
end

function fit_logistic(tθ, θ)
    lb = [0.0, 0]
    ub = [100.0, 5pi]
    p0 = [0.1, pi]
    model(x, p) = logistic.(x, p[2], p[1], 0)
    fit = curve_fit(model, tθ, θ, p0, lower = lb, upper = ub)
    fit.param
end

function getIQR(x)
    α = 0.05
    d = Truncated(Distributions.fit(Normal, x), 0, Inf)
    c1, μ, c2 = quantile(d, [α/2, 0.5, 1 - α/2])
    (; c1, μ, c2)
end

function create_track(k)
    xy = [zero(Point2f)]
    t = 0
    while last(last(xy)) > -50
        t += 0.1
        Δ = reverse(sincos(logistic(t, 2π, k, 0) + π/2))
        push!(xy, xy[end] + Point2f(Δ))
    end
    spl = ParametricSpline(range(0, t, length(xy)), stack(xy))
    xys = Point2f.(spl.(range(0, t, 100)))
    xys[end] = Point2f(xys[end][1], -50)
    return xys
end

logistic(x, L, k, x₀) = L / (1 + exp(-k*(x - x₀))) - L/2

function arclength(spl, t1, t2; kws...)
    knots = get_knots(spl)
    filter!(t -> t1 < t < t2, knots)
    pushfirst!(knots, t1)
    push!(knots, t2)
    s, _ = quadgk(t -> norm(derivative(spl, t)), knots; kws...)
    return s
end

function get_pathlength(t, spl)
    n = length(t)
    l = zeros(n)
    for i in 2:n
        l[i] = arclength(spl, t[i-1], t[i])
    end
    return cumsum(l)
end

function get_mean_speed(t, spl, poi_index)
    t1 = t[poi_index]
    t2 = t[end]
    s = arclength(spl, t1, t2)
    return s / (t2 - t1)
end

# function t2length(t, spl)
#     n = length(t)
#     l = Vector{Float64}(undef, n)
#     l[1] = 0.0
#     for i in 2:n
#         l[i] = first(quadgk(t -> norm(derivative(spl, t)), t[1], t[i]))
#     end
#     return l
# end
