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
    unwrap!(θ)
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


function fun(t, xy, poi_index)
    n = length(t)
    indices = collect(1:poi_index:n)
    if n ∉ indices
        push!(indices, n)
    end
    xy = xy[indices]
    t = t[indices]
    t .-= t[1]
    Δt = mean(diff(t))
    xy = xy[2:end]
    p1 = xy[1]
    t = range(0; length = length(xy), step = Δt)
    Δ = diff(xy)
    θ = splat(atan).(reverse.(Δ))
    unwrap!(θ)
    t = range(0; length = length(θ), step = Δt)
    return t, θ
end

GLMakie.activate!()

row = runs[7,:]
t, xy, poi_index = (row.t, row.xy, row.poi_index)
xy1 = [p1 for (p1, p2) in zip(xy[1:end-1], xy[2:end]) if p1 ≠ p2]
xy1 = xy1[1:10:end]
n = length(xy1)
itp = interpolate(stack(xy1)', (BSpline(Cubic(Natural(OnGrid()))), NoInterp()))

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
scatter!(ax, xy)
xy2 = [SVector{2, Float64}(itp(i, 1:2)) for i in range(1, n, 100n)]
lines!(ax, xy2, color = :red)

Δxy = minimum(norm, diff(xy2))/10

function radiallength(itp, i, Δxy)
    p1 = SVector{2, Float64}(itp(i, 1:2))
    fun(t) = abs2(norm(SVector{2, Float64}(itp(t, 1:2)) - p1) - Δxy)
    o = optimize(fun, i, i + 1)
    # if !Optim.converged(o)
    #     error()
    # end
    # if o.minimum > 0.001
    #     error()
    # end
    return Optim.minimizer(o)
end

ts = Float64[]
push!(ts, 1.0)
while ts[end] < n - 1
    push!(ts, radiallength(itp, ts[end], Δxy))
end

xy2 = [SVector{2, Float64}(itp(t, 1:2)) for t in ts]
maximum(abs2.(norm.(diff(xy2)) .- Δxy))


function get_θ(itp, t)
    x = only(gradient(itp, t, 1))
    y = only(gradient(itp, t, 2))
    atan(y, x)
end

θ = get_θ.(Ref(itp), ts)
unwrap!(θ)

xy2 = cumsum(Δxy .* Point2f.(reverse.(sincos.(θ))))

fc = 0.5
flt = digitalfilter(Lowpass(fc), Butterworth(5));
# flt = digitalfilter(Lowpass(fc), Elliptic(4, 3, 2))
θ2 = filtfilt(flt, θ)
# lines(θ)
# lines!(θ2)
xy3 = cumsum(Δxy .* Point2f.(reverse.(sincos.(θ2))))
fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
# scatter!(ax, xy)
lines!(ax, xy2, color = :red)
lines!(ax, xy3, color = :green)






Δt = step(t)
Δ = diff(xy)
filter!(!iszero ∘ norm, Δ)
θ = splat(atan).(reverse.(Δ))
unwrap!(θ)

l = norm.(Δ)
L = cumsum(l)
L .-= L[1]
itp = interpolate((L, ), θ, Gridded(Linear()))
μ = mean(l)
lL = range(0, L[end], step = μ)
lθ = itp.(lL)

lines(L, θ)
lines!(lL, lθ)

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
lines!(ax, xy)
xy2 = cumsum(μ .* Point2f.(reverse.(sincos.(lθ))))
scatter!(ax, xy2)

signal = θ

# N = length(signal) - 1
# Ts = step(t)
# fs = 1 / Ts
# t0 = 0
# tmax = t0 + N * Ts
# t = t0:Ts:tmax
# F = fftshift(fft(signal))
# freqs =  fftshift(fftfreq(length(t), fs))

fc = 0.04
flt = digitalfilter(Lowpass(fc), Butterworth(5));
lθ2 = filtfilt(flt, lθ)  #filtered signal
fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
lines!(ax, xy)
xy2 = cumsum(μ .* Point2f.(reverse.(sincos.(lθ2))))
lines!(ax, xy2)







fig = Figure()
ax = Axis(fig[1,1])
lines(θ)

using Interpolations



itp = cubic_spline_interpolation(L, θ)

using FFTW

# Number of points
N = 2^12 - 1
# Sample spacing
Ts = 1 / (1.1 * N)
# Sample rate
fs = 1 / Ts
# Start time
t0 = 0
tmax = t0 + N * Ts

# time coordinate
t = t0:Ts:tmax

# The underlying signal here is the sum of a sine wave at 60 cycles per second
# and its second harmonic (120 cycles per second) at half amplitude. We have
# discrete observations (samples) of this signal at each time `t`, with `fs`
# samples per second.

signal = sin.(2π * 60 * t) + .5 * sin.(2π * 120 * t)

# The `fft` function calculates the (discrete) Fourier transform of its input.
# The first half of the returned array contains the positive frequencies, while
# the second half contains the negative ones. For visualization purposes, we
# rearrange the array to have the zero-frequency at the center.

signal = θ
N = length(signal) - 1
Ts = step(t)
fs = 1 / Ts
# Start time
t0 = 0
tmax = t0 + N * Ts
# time coordinate
t = t0:Ts:tmax


F = fftshift(fft(signal))
freqs =  fftshift(fftfreq(length(t), fs))

fc = 0.1
flt = digitalfilter(Lowpass(fc), Butterworth(6));
d2 = filtfilt(flt, signal)  #filtered signal
# Plot
fig = Figure()
ax = Axis(fig[1,1], title="Signal", xlabel="time (s)")
lines!(ax, t, signal)
lines!(ax, t, d2)

ax = Axis(fig[1,2], title="Spectrum", limits = ((0, 5), nothing), xlabel="frequency (Hz)")
lines!(ax, freqs, abs.(F))










fig = Figure()
ax = Axis(fig[1,1])
t, θ = fun(t, xy, poi_index)
polynom = Polynomials.fit(t, θ, 2)
scatter!(t, θ)
ts = range(t[1], t[end], 100)
lines!(ts, polynom.(ts))
