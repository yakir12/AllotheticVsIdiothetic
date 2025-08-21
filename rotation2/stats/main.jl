using CSV, DataFrames, GLM
using GLMakie, AlgebraOfGraphics

df = CSV.read("data.csv", DataFrame)

m = glm(@formula(lap ~ elevation), df, Binomial(), LogitLink())

n = 100
newdf = DataFrame(elevation = range(0, 90, n), lap = falses(n))
plu = predict(m, newdf, interval = :confidence)
newdf.prediction = disallowmissing(plu.prediction)
newdf.lower = disallowmissing(plu.lower)
newdf.upper = disallowmissing(plu.upper)

df2 = combine(groupby(df, :elevation), :lap => count => :lap, nrow => :n)
transform!(df2, [:lap, :n] => ByRow(/) => :p)

band(newdf.elevation, newdf.lower, newdf.upper)
lines!(newdf.elevation, newdf.prediction)
scatter!(df2.elevation, df2.p)

using Turing
using StatsFuns: logistic

@model function bmodel(xs, ys, n)
    # Priors
    α ~ Normal(-4, 2)
    β ~ Normal(0, 2)
    # σ² ~ truncated(Normal(0, 1), lower=0)
    # Likelihood
    for i in 1:n
        p = logistic(α + β*xs[i])
        ys[i] ~ Bernoulli(p)
    end
end
posterior = bmodel(df.elevation, df.lap, nrow(df))
chain = sample(posterior, NUTS(), 1000)

fig = Figure()
for (i, var_name) in enumerate(chain.name_map.parameters)
    draw!(
        fig[i, 1],
        data(chain) *
        mapping(var_name; color=:chain => nonnumeric) *
        AlgebraOfGraphics.density() *
        visual(fillalpha=0)
    )
end

# x = -10:100
# mylogistic(x, x₀, k, L) = L/(1 + exp(-k*(x - x₀)))
# lines(x, mylogistic.(2x, 90, 1/30, 180)/2; axis = (; limits = (nothing, (0, 90))))
# ablines!(0, 1)

# transθ(x) = mylogistic.(x, 45, 1/20, 90)
# transϵ(x) = mylogistic.(x, 90, 1/30, 180)
# trans(x) = mylogistic.(x, 45, 1/20, 90)

function my_link(ϵ, θ, a, b) 
    # 2atand(tand(trans(ϵ))/cosd(trans(θ)))/180
    (2atand(tand(ϵ/2)/cosd(θ)) - ϵ)/(180 - ϵ)*a + b
end



scatter(df2.elevation, df2.p)
for ϵ in 10:2:15
    # p = @. (2atand(tand(ϵ/2)/cosd(elevation)) - ϵ)/(180 - ϵ)
    # lines!(elevation, p)
    lines!(elevation, my_link.(ϵ, elevation, 1/4, 0.02))
end

@model function bmodel(xs, ys, n)
    # Priors
    ϵ ~ truncated(Normal(10, 2), 0.1, 179.9)
    a ~ truncated(Normal(0.5, 1), 0, 1)
    b ~ truncated(Normal(0.0, 1), 0, 1)
    ps = my_link.(ϵ, xs, a, b)
    if any(p -> p < 0 || p > 1, ps)
        @addlogprob! -Inf
        # Exit the model evaluation early
        return nothing
    end
    # Likelihood
    for i in 1:n
        ys[i] ~ Bernoulli(ps[i])
    end
    return nothing
end
df2 = transform(df, :elevation => ByRow(x -> x == 90 ? 89.9 : x), renamecols = false)
posterior = bmodel(df2.elevation, df2.lap, nrow(df))
chain = sample(posterior, NUTS(), 1000)
# chain = sample(posterior, Prior(), 1000)
# chain = sample(posterior, SMC(), 1000)

# sampler   = RandPermGibbs(SliceSteppingOut(.2))
# chain = sample(posterior, externalsampler(sampler), 1000)

fig = Figure()
for (i, var_name) in enumerate(chain.name_map.parameters)
    draw!(
          fig[i, 1],
          data(chain) *
          mapping(var_name; color=:chain => nonnumeric) *
          AlgebraOfGraphics.density() *
          visual(fillalpha=0)
         )
end

param = DataFrame(quantile(chain, q = [0.025, 0.5, 0.975]))
rename!(param, "2.5%" => :lower, "50.0%" => :prediction, "97.5%" => :upper)
ϵ = subset(param, :parameters => ByRow(==(:ϵ)))
a = subset(param, :parameters => ByRow(==(:a)))
b = subset(param, :parameters => ByRow(==(:b)))

n = 100
elevation = range(0, 90, n)
newdf = DataFrame(elevation = elevation, 
                  lower = my_link.(only(ϵ.lower), elevation, only(a.lower), only(b.lower)), 
                  prediction = my_link.(only(ϵ.prediction), elevation, only(a.prediction), only(b.prediction)), 
                  upper = my_link.(only(ϵ.upper), elevation, only(a.upper), only(b.upper)), 
                 )

df3 = combine(groupby(df, :elevation), :lap => count => :lap, nrow => :n)
transform!(df3, [:lap, :n] => ByRow(/) => :p)

band(newdf.elevation, newdf.lower, newdf.upper, color = :gray)
lines!(newdf.elevation, newdf.prediction, color = :white)
scatter!(df3.elevation, df3.p, color = :black)

