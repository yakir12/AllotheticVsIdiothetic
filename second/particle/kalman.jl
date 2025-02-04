using LowLevelParticleFilters, GLMakie, Random
Random.seed!(1)

# Create a time series for filtering
x = [zeros(50); 0:100]
T = length(x)
Y = x + randn(T)
lines(x, label = "Measurement")
lines!(Y, label = "True state to be tracked")


y = [[y] for y in Y] # create a vector of vectors for the KF
u = fill([], T) # No inputs in this example :(

# Define the model
Ts = 1
A = [1 Ts; 0 1]
B = zeros(2, 0)
C = [1 0]
D = zeros(0, 0)
R2 = [1;;]

σws = [1e-2, 1e-5] # Dynamics noise standard deviations

for σw in σws
    R1 = σw*[Ts^3/3 Ts^2/2; Ts^2/2 Ts] # The dynamics noise covariance matrix is σw*Bw*Bw' where Bw = [Ts^2/2; Ts]
    kf = KalmanFilter(A, B, C, D, R1, R2)
    yh = []
    measure = LowLevelParticleFilters.measurement(kf)
    for t = 1:T # Main filter loop
        kf(u[t], y[t]) # Performs both prediction and correction
        xh = state(kf)
        yht = measure(xh, u[t], nothing, t)
        push!(yh, yht)
    end
    Yh = vec(stack(yh))
    lines!(Yh, label="Estimate σ_w = $σw")
end
axislegend()





















using LowLevelParticleFilters, LinearAlgebra, StaticArrays, Distributions, GLMakie, Random
using CSV

tbl = CSV.File(joinpath("../tracks/1.csv"))
xy = SVector{2, Float64}.(tbl.x, tbl.y)

#
# using DelimitedFiles
# path = "../track.csv"
# xyt = readdlm(path)
# tosvec(y) = reinterpret(SVector{length(y[1]),Float64}, reduce(hcat,y))[:] |> copy # helper function
# y = tosvec(collect(eachrow(xyt[:,1:2])))

N = 2000 # Number of particles in the particle filter
n = 4 # Dimension of state: we have position (2D), speed and angle
p = 2 # Dimension of measurements, we can measure the x and the y
@inline pos(s) = s[SVector(1,2)]
@inline vel(s) = s[3]
@inline ϕ(s) = s[4]
@inline mymode(s) = s[5]

dgσ = 1 # the deviation of the measurement noise distribution
dvσ = 0.3 # the deviation of the dynamics noise distribution
ϕσ  = 0.5
const switch_prob = 0.03 # Probability of mode switch
const dg = MvNormal(@SVector(zeros(p)), dgσ^2) # Measurement noise Distribution
const df = LowLevelParticleFilters.TupleProduct((Normal.(0,[1e-1, 1e-1, dvσ, ϕσ])...,Binomial(1,switch_prob)))
d0 = MvNormal(SVector(xy[1]..., 0.5, atan((xy[2]-xy[1])...), 0), [3.,3,2,2,0])
const noisevec = zeros(5) # cache vector

@inline function dynamics(s,u,p,t,noise=false)
    # current state
    m = mymode(s)
    v = vel(s)
    a = ϕ(s)
    p = pos(s)
    # get noise
    if noise
        y_noise, x_noise, v_noise, ϕ_noise,_ = rand!(df, noisevec)
    else
        y_noise, x_noise, v_noise, ϕ_noise = 0.,0.,0.,0.
    end
    # next state
    v⁺ = max(0.999v + v_noise, 0.0)
    m⁺ = Float64(m == 0 ? rand() < switch_prob : true)
    a⁺ = a + (ϕ_noise*(1 + m*10))/(1 + v⁺) # next state velocity is used here
    p⁺ = p + SVector(y_noise, x_noise) + SVector(sincos(a))*v⁺ # current angle but next velocity
    SVector{5,Float64}(p⁺[1], p⁺[2], v⁺, a⁺, m⁺) # all next state
end
function measurement_likelihood(s,u,xy,p,t)
    logpdf(dg, pos(s)-xy) # A simple linear measurement model with normal additive noise
end
@inline measurement(s,u,p,t,noise=false) = s[SVector(1,2)] + noise*rand(dg) # We observe the position coordinates with the measurement

u = zeros(length(xy))
pf = AuxiliaryParticleFilter(AdvancedParticleFilter(N, dynamics, measurement, measurement_likelihood, df, d0))
T = length(xy)
sol = forward_trajectory(pf,u[1:T],xy[1:T])
(; x,w,we,ll) = sol
# plot(sol, markerstrokecolor=:auto, m=(2,0.5))

xh = mean_trajectory(x,we)

"plotting helper function"
function to1series(x::AbstractVector, xy)
    r,c = size(xy)
    y2 = vec([xy; fill(Inf, 1, c)])
    x2 = repeat([x; Inf], c)
    x2,y2
end
to1series(xy) = to1series(1:size(xy,1),xy)

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
lines!(ax, xy, label = "raw")

M = 1000 # Number of smoothing trajectories, NOTE: if this is set higher, the result will be better at the expense of linear scaling of the computational cost.
sb,ll = smooth(pf, M, u, xy) # Sample smoothing particles (b for backward-trajectory)
sbm = smoothed_mean(sb)     # Calculate the mean of smoothing trajectories
sbt = smoothed_trajs(sb)    # Get smoothing trajectories
lines!(ax, sbm[1,:],sbm[2,:], label = "smoothed")
axislegend(ax, position = :rb)

fig
# plot!(fig2, identity.(sbm'[:,4]), lab="smoothed")
#
# plot(xh[:,5], lab="Filtering")
# plot!(to1series(sbt[5,:,:]')..., lab="Smoothing", title="Mode trajectories", l=(:black,0.2))
