using StaticArrays, LinearAlgebra

function next_step(θ, brwθ, crwθ, w)
    brwyx = sincos(brwθ)
    crwyx = sincos(θ + crwθ)
    y, x = @. w*brwyx + (1 - w)*crwyx
    atan(y, x)
end

w = 0.3
brwθ = 0.3
crwθ = 0.1

n = 30
θs = accumulate(1:n, init = 0) do θ, _
    next_step(θ, brwθ, crwθ, w)
end

function fun(brwθ, crwθ)
    θ̂s = accumulate(1:n, init = 0) do θ, _
        next_step(θ, brwθ, crwθ, w)
    end
    abs2.(θs .- θ̂s)
end

using GLMakie

p1 = range(0, 1, 1000)
p2 = range(0, 1, 1000)
z = [fun(p1i, p2i) for p1i in p1, p2i in p2]

fig = Figure()
for i in 1:n
    ax = Axis(fig[1,i])
    h = heatmap!(ax, p1, p2, getindex.(z, i); colorscale = log10)
    scatter!(ax, brwθ, crwθ)
    Colorbar(fig[2,i], h; vertical = false)
end

fig = Figure()
ax = Axis(fig[1,1])
h = heatmap!(ax, p1, p2, sum.(z); colorscale = log10)
scatter!(ax, brwθ, crwθ)
Colorbar(fig[2,1], h; vertical = false)
