const SV = SVector{2, Float64}

###

function remove_loops(xy)
    inds, _ = self_intersections(Point2f.(xy))
    rngs = splat(UnitRange).(Iterators.partition(inds, 2))
    filter!(<(51) ∘ length, rngs)
    tokill = vcat(rngs...)
    unique!(tokill)
    sort!(tokill)
    keep::Vector{Int} = collect(eachindex(xy))
    setdiff!(keep, tokill)
    return xy[keep]
end

function glue_intervention!(xy, intervention::Real)
    h = 1
    diffs = diff(xy[Ti = 0.0 .. intervention + h])
    Δs = norm.(diffs)
    μ = mean(Δs[Ti = 0.0 .. intervention - h])
    σ = std(Δs[Ti = 0.0 .. intervention - h], mean = μ)
    threshold = μ + 2σ
    v = Δs[Ti = intervention - h .. intervention + h]
    i = findfirst(>(threshold), v)
    if isnothing(i)
        return 0.0
    else
        t1 = lookup(v, Ti)[i]
        step = diffs[Ti = At(t1)] - μ*normalize(diffs[Ti = At(t1)])
        t2 = lookup(v, Ti)[i + 1]
        xy[Ti = t2 .. Inf] .-= Ref(step)
        return v[i]
    end
end

function sparseify(xy)
    spl = ParametricSpline(lookup(xy, Ti), stack(xy))
    # tl = round(Int, t[1]):round(Int, t[end])
    tl = range(extrema(lookup(xy, Ti))..., step = 0.5)
    # poi_index = findfirst(==(round(Int, poi)), tl)
    DimVector(SV.(spl.(tl)), Ti(tl))
end

function has_stops(ij)
    for i in 2:length(ij)
        if ij[i] == ij[i - 1]
            return true
        end
    end
    return false
end

function remove_stops(ij)
    keep = [1]
    for i in 2:length(ij)
        if ij[i] ≠ ij[i - 1]
            push!(keep, i)
        end
    end
    ij[keep]
end

function get_tij(file)
    tij = CSV.File(file)
    t = range(tij.t[1], tij.t[end], length = length(tij))
    pixels = DimVector(CameraCalibrations.RowCol.(tij.i, tij.j), Ti(t))
end
function smooth(xy)
    t = lookup(xy, Ti)
    spl = ParametricSpline(t, stack(xy), k = 3, s = 25)
    tl = range(first(t), last(t), 100length(t))
    DimVector(SV.(spl.(tl)), Ti(tl))
end

function center2start(xy)
    trans = Translation(-first(xy))
    trans.(xy)
end

function intersection(orig, dir)
    b = -orig⋅dir
    disc = b^2 - orig⋅orig + 1
    if disc ≥ 0
        d = sqrt(disc)
        t2 = b + d
        if t2 ≥ 0
            t1 = b - d
            return t1 > 0 ? (t1, t2) : (Inf, t2)
        end
    end
    return (Inf, Inf)
end

function cropto(xy, l)
    i = findfirst(>(l) ∘ norm, xy)
    isnothing(i) && return copy(xy)
    p1 = xy[i - 1]
    p2 = xy[i]
    dir = normalize(p2 - p1)
    _, d = intersection(p1/l, normalize(p2 - p1))
    p2 = p1 + d*l*dir
    cropped = xy[1:i]
    cropped[i] = p2
    return cropped
end

function rotate2poi(xy, poi)
    p2 = xy[Ti = Near(poi)]
    θ = π/2 - atan(reverse(p2)...)
    rot = LinearMap(Angle2d(θ))
    rot.(xy)
end

function center2poi_and_crop(xy, poi)
    trans = Translation(-xy[Ti = Near(poi)])
    trans.(xy[Ti = poi..Inf])
end

###

# function trectify(fs, xs)
#     n = length(xs)
#     ys = Vector{DimVector{SV}}.(undef, n)
#     Threads.@threads for i in 1:n
#         y = fs[i].(xs[i])
#         y .-= Ref(y[1])
#         ys[i] = y
#     end
#     return ys
# end

function get_calibration(file)
    c = CameraCalibrations.load(file)
    f = rectification(c, findfirst(contains("extrinsic"), c.files))
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

function impute_poi_time(xy)
    p1 = xy[1]
    for i in eachindex(xy)
        if norm(xy[i] - p1) > 10
            return lookup(xy, Ti)[i - 1]
        end
    end
    return last(lookup(xy, Ti))
end

