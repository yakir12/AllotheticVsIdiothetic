tosecond(t::T) where {T <: Dates.TimePeriod} = t / convert(T, Dates.Second(1))
tosecond(t::T) where {T <: Dates.TimeType} = tosecond(t - Time(0))

# function get_spline(track)
#     t, xy = track
#     XY = reshape(collect(Iterators.flatten(xy)), 2, :)
#     ParametricSpline(t, XY; s = 10000, k = 2)
# end

function smooth_rotate_split(track, POI)
    t, xy = track
    XY = reshape(collect(Iterators.flatten(xy)), 2, :)
    spl = ParametricSpline(t, XY; s = 10000, k = 2)

    trans = Translation(-spl(t[1]))
    p1 = SV(spl(t[1]))
    p2 = SV(spl(POI))
    x, y = normalize(p2 - p1)
    θ = π/2 - atan(y, x)
    rot = recenter(LinearMap(Angle2d(θ)), p1)
    t2xy = trans ∘ rot ∘ SV ∘ spl

    before = t2xy.(range(t[1], POI, step = 1/25))
    after = t2xy.(range(POI, t[end], step = 1/25))

    return (before, after)
end


function plotit(name, start, stop, spline)
    fig = Figure();
    ax = Axis(fig[1,1], aspect=DataAspect())
    ts = range(start, stop, 100)
    xy = Point2f.(spline.(ts))
    xy .-= xy[1]
    lines!(ax, xy)
    lines!(ax, Circle(zero(Point2f), 350), color=:black)
    save("$name.pdf", fig)
end

function plotthem(names, starts, stops, splines)
    fig = Figure();
    ax = Axis(fig[1,1], aspect=DataAspect())
    lins = []
    for (name, start, stop, spline) in zip(names, starts, stops, splines)
        ts = range(start, stop, 100)
        xy = Point2f.(spline.(ts))
        xy .-= xy[1]
        push!(lins, lines!(ax, xy))
    end
    lines!(ax, Circle(zero(Point2f), 350), color=:black)
    Legend(fig[1,2], lins, names; nbanks=2)
    # axislegend(ax, position = (1.2,1))
    save("all.pdf", fig)
end

function save_csv(name, track)
    t, xy = track
    df = DataFrame(second = t, x = first.(xy), y = last.(xy))
    CSV.write("$name.csv", df)
end

function save_vid(name, file, track)
    t, xy = track
    vid = openvideo(file)
    sz = out_frame_size(vid)
    h = 200
    fig = Figure(size=(h, round(Int, h * sz[2] / sz[1] / VideoIO.aspect_ratio(vid))), figure_padding=0)
    ax = Axis(fig[1,1])
    img = Observable(rotr90(read(vid)))
    image!(ax, img)
    y, x = xy[1]
    point = Observable(Point2f(x, sz[2] - y))
    scatter!(ax, point, marker='+', color=:red)
    hidespines!(ax)
    hidedecorations!(ax)
    t₀ = gettime(vid)
    t .+= t₀
    seek(vid, t[1])
    framerate = round(Int, 2length(t)/(t[end] - t[1]))
    record(fig, "$name.mp4", zip(xy, vid); framerate) do (xy, frame)
        img[] = rotr90(frame)
        y, x = xy
        point[] = Point2f(x, sz[2] - y)
    end
end

function angle_between(oldu, newu)
  x1, y1 = oldu
  x2, y2 = newu
  sign((x2 - x1)*(y2 + y1))*angle(oldu, newu)
end

function cumulative_angle(xy)
    δ = filter(!iszero ∘ sum, diff(xy))
    # δ = filter(!iszero ∘ sum, diff(transfrm.(range(t1, t2, 1000))))
    rad2deg(sum(splat(angle_between), partition(δ, 2, 1)))
end

# xy = [SV(0,0), SV(1,0.1), SV(2, -0.1)]
# save("figure.pdf", scatter(xy, axis=(;aspect=DataAspect())))
# δ = diff(xy)
# rad2deg(sum(splat(angle_between), partition(δ, 2, 1)))


cordlength(xy) = norm(diff([xy[1], xy[end]]))

function curvelength(xy)
    p0 = xy[1]
    s = 0.0
    for p1 in xy[2:end]
        s += norm(p1 - p0)
        p0 = p1
    end
    return s
end

# function straightness(transform, t1, t2)
#     xy = transform.(range(t1, t2, 1000))
#     cordlength(xy) / curvelength(xy)
# end


# function rotate(spl, t2)
#     t1 = first(spl.t)
#     trans = Translation(-spl(t1))
#     p1 = SV(spl(t1))
#     p2 = SV(spl(t2))
#     x, y = normalize(p2 - p1)
#     θ = π/2 - atan(y, x)
#     rot = recenter(LinearMap(Angle2d(θ)), p1)
#     return trans ∘ rot
# end

