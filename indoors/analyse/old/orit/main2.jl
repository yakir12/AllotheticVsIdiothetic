# using StaticArrays, LinearAlgebra, GLMakie, Polynomials, Optim, CSV, Statistics
#
# function next_step(θ₁, crwθ)
#     θ₂ = θ₁ + crwθ
#     y, x = sincos(θ₂)
#     step = SVector{2, Float64}(x, y)
#     Δ = normalize(step) # normalize the step to 1
#     return (θ₂, Δ)
# end
#
# function get_track(θ₀, crw, radius, tΔ)
#     track = SVector{2, Float64}[]
#     xy = zero(SVector{2, Float64})
#     push!(track, xy)
#     θ = θ₀
#     t = 0.0
#     while norm(xy) < radius
#         t += tΔ
#         θ, Δ = next_step(θ, crw(t))
#         xy += Δ
#         push!(track, xy)
#     end
#     return track
# end
#
#
# xy = [SVector{2, Float64}(row...) for row in CSV.File("../centered rotataed raw/160.csv")]
# # xy = xy[100:end]
# xy .-= Ref(xy[1])
# θ₀ = -atan(reverse(xy[end])...)
# radius = norm(xy[end])
#
# function get_track(p)
#     crw = Polynomial(p[2:end])
#     get_track(p[1], crw, radius, 1)
# end
#
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# lines!(ax, xy)
# p = [3.7, 0.021, 0.0002 ]
# track = get_track(p)
# lines!(ax, track)
# scatter!(xy[1:60:end])
#

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
    polynom = fit(t, θ, 2)
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
    return track
end




file = CSV.File("../calibrated/160.csv")
xy = [SVector{2, Float64}(row.x, row.y) for row in file]
t = file.t
poi_index = 60
xy .-= Ref(xy[1])
rot = get_rotation(xy[poi_index])
xy .= rot.(xy)

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
scatter!(ax, xy)
lines!(ax, track, color = :red)
fig


