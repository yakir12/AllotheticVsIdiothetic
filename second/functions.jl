tosecond(t::T) where {T <: TimeType} = tosecond(t - Time(0))
tosecond(t::T) where {T <: TimePeriod} = t / convert(T, Dates.Second(1))

function calibrate_and_smooth(c, track, s, k)
    xy_pixels = RowCol.(track.x, track.y)
    xy_cm = CameraCalibrationMeta.rectification(c).(xy_pixels)
    XY = reduce(hcat, xy_cm)
    spl = ParametricSpline(track.t, XY; s, k)
    return SV ∘ spl ∘ tosecond
end

function center_and_rotate(start, POI, spl)
    p1 = spl(start)
    p2 = spl(POI)
    trans = Translation(-p1)
    x, y = normalize(p2 - p1)
    θ = π/2 - atan(y, x)
    rot = recenter(LinearMap(Angle2d(θ)), p1)
    return trans ∘ rot
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

