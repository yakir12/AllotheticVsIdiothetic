using StaticArrays, LinearAlgebra, GLMakie, Polynomials, Optim, CSV, Statistics

function next_step(θ, brwθ, crwθ, w)
    brwyx = sincos(brwθ)
    crwyx = sincos(θ + crwθ)
    y, x = @. w*brwyx + (1 - w)*crwyx
    θ = atan(y, x)
    step = SVector{2, Float64}(x, y)
    Δ = 0.3normalize(step) # normalize the step to 1
    return (θ, Δ)
end

function get_track(θ₀, nsteps, brw, crw, w)
    track = Vector{SVector{2, Float64}}(undef, nsteps)
    xy = zero(SVector{2, Float64})
    track[1] = xy
    θ = θ₀
    for i in 2:nsteps
        θ, Δ = next_step(θ, brw(i), crw(i), w)
        xy += Δ
        track[i] = xy
    end
    return track
end


xy = [SVector{2, Float64}(row...) for row in CSV.File("../centered rotataed raw/160.csv")]
xy = xy[100:end]
xy .-= Ref(xy[1])
θ₀ = atan(reverse(xy[end])...)
nsteps = length(xy)

function get_track(p)
    brw_poly = Polynomial(p[1:2])
    brw(i) = brw_poly(i/nsteps)
    crw_poly = Polynomial(p[3])
    crw(i) = crw_poly(i/nsteps)
    w = 0.3
    get_track(θ₀, nsteps, brw, crw, w)
end

function fun(p)
    track = get_track(p)
    sqrt(mean(LinearAlgebra.norm_sqr, track .- xy))
end

x0 = rand(3)
o = optimize(fun, x0, Optim.Options(iterations = 100000))
p = o.minimizer

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
lines!(ax, xy)
track = get_track(p)
lines!(ax, track)
fig
