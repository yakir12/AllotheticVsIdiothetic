tosecond(t::T) where {T <: TimeType} = tosecond(t - Time(0))
tosecond(t::T) where {T <: TimePeriod} = t / convert(T, Dates.Second(1))

SimpTrack.track(file, start::Time, stop::Time; kwargs...) = track(file, tosecond(start), tosecond(stop); kwargs...)

function get_smooth_calibrate(clb, trk, s, k)
    t, xy = trk
    xy = clb.itform.(xy)
    # smooth 
    XY = reshape(collect(Iterators.flatten(xy)), 2, :)
    spl = ParametricSpline(t, XY; s, k)
    return SV ∘ spl ∘ tosecond
end

function get_rotation(start, POI, spl)
    # center and rotate
    p1 = spl(start)
    p2 = spl(POI)
    trans = Translation(-p1)
    x, y = normalize(p2 - p1)
    θ = π/2 - atan(y, x)
    rot = recenter(LinearMap(Angle2d(θ)), p1)
    return trans ∘ rot
end

function get_transformation(clb, trk, start, POI; s = 100, k = 2)
    spl = get_smooth_calibrate(clb, trk, s, k)
    trans = get_rotation(start, POI, spl)
    return trans ∘ spl
end

# function rotate_trim(start, stop, POI, spl; cutoff = 50)
#     # resample
#     t = range(tosecond(start), tosecond(stop), round(Int, 30tosecond(stop - start))) # 30 fps
#     xy = SV.(spl.(t))
#
#     # center and rotate
#     trans = Translation(-xy[1])
#     p1 = xy[1]
#     i = findfirst(≥(tosecond(POI)), t)
#     p2 = xy[i]
#     x, y = normalize(p2 - p1)
#     θ = π/2 - atan(y, x)
#     rot = recenter(LinearMap(Angle2d(θ)), p1)
#     cmp = trans ∘ rot
#     xy .= cmp.(xy)
#
#     # trim
#     n = length(xy)
#     j = findfirst(>(cutoff) ∘ norm, xy)
#     if !isnothing(j)
#         @assert j ≥ i "The trimming cutoff must be larger than the POI, $i"
#         xy = xy[1:j - 1]
#         t = t[1:j - 1]
#     end
#
#     return (; t, i, xy)
# end

# function plotit(name, start, stop, spline)
#     fig = Figure();
#     ax = Axis(fig[1,1], aspect=DataAspect())
#     ts = range(start, stop, 100)
#     xy = Point2f.(spline.(ts))
#     xy .-= xy[1]
#     lines!(ax, xy)
#     lines!(ax, Circle(zero(Point2f), 350), color=:black)
#     save("$name.pdf", fig)
# end

# function plotthem(names, starts, stops, splines)
#     fig = Figure();
#     ax = Axis(fig[1,1], aspect=DataAspect())
#     lins = []
#     for (name, start, stop, spline) in zip(names, starts, stops, splines)
#         ts = range(start, stop, 100)
#         xy = Point2f.(spline.(ts))
#         xy .-= xy[1]
#         push!(lins, lines!(ax, xy))
#     end
#     lines!(ax, Circle(zero(Point2f), 350), color=:black)
#     Legend(fig[1,2], lins, names; nbanks=2)
#     # axislegend(ax, position = (1.2,1))
#     save("all.pdf", fig)
# end

function save_csv(name, track)
    t, xy = track
    df = DataFrame(second = t, x = first.(xy), y = last.(xy))
    CSV.write("$name.csv", df)
end

# function save_vid(name, file, track)
#     t, xy = track
#     vid = openvideo(file)
#     sz = out_frame_size(vid)
#     h = 200
#     fig = Figure(size=(h, round(Int, h * sz[2] / sz[1] / VideoIO.aspect_ratio(vid))), figure_padding=0)
#     ax = Axis(fig[1,1])
#     img = Observable(rotr90(read(vid)))
#     image!(ax, img)
#     y, x = xy[1]
#     point = Observable(Point2f(x, sz[2] - y))
#     scatter!(ax, point, marker='+', color=:red)
#     hidespines!(ax)
#     hidedecorations!(ax)
#     t₀ = gettime(vid)
#     t .+= t₀
#     seek(vid, t[1])
#     framerate = round(Int, 2length(t)/(t[end] - t[1]))
#     record(fig, "$name.mp4", zip(xy, vid); framerate) do (xy, frame)
#         img[] = rotr90(frame)
#         y, x = xy
#         point[] = Point2f(x, sz[2] - y)
#     end
# end

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

# cumulative_angle(rotated) = cumulative_angle(rotated.xy[rotated.i:end])
# function cumulative_angle(xy::Vector{SV})
#     δ = filter(!iszero ∘ sum, diff(xy[1:60:end])) # cause 30 fps
#     α = sum(splat(angle_between) ∘ reverse, partition(δ, 2, 1); init=0.0)
#     rad2deg(α)
# end




# function angle_between(oldu, newu)
#   x1, y1 = oldu
#   x2, y2 = newu
#   sign((x2 - x1)*(y2 + y1))*angle(oldu, newu)
# end
#
# cumulative_angle(rotated) = cumulative_angle(rotated.xy[rotated.i:end])
# function cumulative_angle(xy::Vector{SV})
#     δ = filter(!iszero ∘ sum, diff(xy))
#     if isempty(δ)
#         0.0
#     else
#         # δ = filter(!iszero ∘ sum, diff(transfrm.(range(t1, t2, 1000))))
#         rad2deg(sum(splat(angle_between), partition(δ, 2, 1); init = 0.0))
#     end
# end
