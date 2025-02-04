using LowLevelParticleFilters, CSV, LinearAlgebra, StaticArrays, Distributions, GLMakie, Random

tbl = CSV.File("../trajectories/1.csv")
y = SVector{2, Float64}.(tbl.x, tbl.y)

nx = 4 # Dimension of state: we have position (2d), speed and angle
ny = 2 # Dimension of measurements, we can measure the x and the y
@inline pos(s) = s[SVector(1,2)]
@inline vel(s) = s[3]
@inline ϕ(s) = s[4]

dgσ = 0.1 # the deviation of the measurement noise distribution
dvσ = .03 # the deviation of the dynamics noise distribution
ϕσ  = 0.5
P = [0.995 0.005; 0.0 1] # Transition probability matrix, we model the search mode as "terminal"
μ = [1.0, 0.0] # Initial mixing probabilities
R1 = Diagonal([1e-1, 1e-1, dvσ, ϕσ].^2)
R2 = dgσ^2*I(ny) # Measurement noise covariance matrix
d0 = MvNormal(SVector(y[1]..., 0.5, atan((y[2]-y[1])...)), [3.,3,2,2])

@inline function dynamics(s,_,modegain,t,w,m)
    # current state
    v = vel(s)
    a = ϕ(s)
    p = pos(s)

    y_noise, x_noise, v_noise, ϕ_noise = w

    # next state
    v⁺ = max(0.999v + v_noise, 0.0)
    a⁺ = a + (ϕ_noise*(1 + m*modegain))/(1 + v⁺) # next state velocity is used here
    p⁺ = p + SVector(y_noise, x_noise) + SVector(sincos(a))*v⁺ # current angle but next velocity
    SVector(p⁺[1], p⁺[2], v⁺, a⁺) # all next state
end
@inline measurement(s,u,p,t) = s[SVector(1,2)] # We observe the position coordinates with the measurement

u = zeros(length(y)) # no control inputs
kffalse = UnscentedKalmanFilter{false,false,true,false}((x,u,p,t,w)->dynamics(x,u,p,t,w,false), measurement, R1, R2, d0; ny, nu=0, p=10)

measure = LowLevelParticleFilters.measurement(kffalse)
yh = []
for t in eachindex(y)
    kffalse(u[t], y[t]) # Performs both prediction and correction
    xh = state(kffalse)
    yht = measure(xh, u[t], nothing, t)
    push!(yh, yht)
end

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
lines!(ax, y, label = "raw")
lines!(ax, yh, label = "fit")
axislegend(ax, position=:rb)

display(fig)


