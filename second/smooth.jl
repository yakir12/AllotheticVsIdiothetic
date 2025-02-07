using StaticArrays, LinearAlgebra, GLMakie, CSV, Rotations, CoordinateTransformations, Statistics, Polynomials, Optim, FresnelIntegrals

function get_rotation(p2)
    θ = π/2 - atan(reverse(p2)...)
    LinearMap(Angle2d(θ))
end

function get_xy(t, _a, b, c, c0)
    a = Complex(_a)
    st, ct = sincos(b^2/4a - c)
    C, S = fresnel((b + 2a*t)/sqrt(2π*a))
    s = sqrt(π/2a)*(ct*S - st*C)
    c = sqrt(π/2a)*(ct*C + st*S)
    SVector{2, Float64}(real(c), real(s)) - c0
end

function get_track(p, t, p1)
    scale, a, b, c = p
    c0 = get_xy(0, a, b, c, zero(SVector{2, Float64}))
    track = [scale*get_xy(ti, a, b, c, c0) + p1 for ti in t]
    return track
end

function get_guess(Δt, xy)
    Δ = diff(xy)
    θ = splat(atan).(reverse.(Δ))
    t = range(0; length = length(θ), step = Δt)
    polynom = Polynomials.fit(t, θ, 2)
    μ = mean(norm, Δ)
    step_per_second = μ/Δt
    p0 = [step_per_second; reverse(coeffs(polynom))]
    return p0
end

function get_coeffs(Δt, xy)
    p1 = xy[1]
    t = range(0; length = length(xy), step = Δt)
    function fun(p)
        track = get_track(p, t, p1)
        sqrt(mean(LinearAlgebra.norm_sqr, track .- xy))
    end
    p0 = get_guess(Δt, xy)
    o = optimize(fun, p0, Optim.Options(iterations = 100000))
    p = o.minimizer
    return p
end

function get_smooth(t, xy, poi_index)
    n = length(t)
    indices = collect(1:poi_index:n)
    if n ∉ indices
        push!(indices, n)
    end

    xy = xy[indices]
    t = t[indices]
    t .-= t[1]
    Δt = mean(diff(t))

    p = get_coeffs(Δt, xy[2:end])
    p1 = xy[2]
    tend = t[end - 1]
    ts = range(0, tend, 1000)
    track = get_track(p, ts, p1)
    pushfirst!(track, zero(SVector{2, Float64}))
    return (; ts, track)
end

