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

x = -10:100
lines(x, logistic.(x ./ 90 .- 1))

function my_link(ϵ, θ) 
    if 0 < ϵ < 180 && 0 < θ < 90
        2atand(tand(ϵ/2)/cosd(θ))/180
    elseif ϵ ≤ 0 || θ ≤ 0
        0.0
    else
        1.0
    end
end
@model function bmodel(xs, ys, n)
    # Priors
    ϵ ~ truncated(Normal(10, 2), 0.1, 179.9)
    # Likelihood
    for i in 1:n
        p = my_link(ϵ, xs[i])
        ys[i] ~ Bernoulli(p)
    end
end
df2 = transform(df, :elevation => ByRow(x -> x == 90 ? 89.9 : x), renamecols = false)
posterior = bmodel(df2.elevation, df2.lap, nrow(df))
# chain = sample(posterior, NUTS(), 1000)
# chain = sample(posterior, Prior(), 1000)
chain = sample(posterior, SMC(), 1000)

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
