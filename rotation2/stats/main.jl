using CSV, DataFrames, GLM, Optim
using GLMakie, AlgebraOfGraphics

n = 100
elevation = range(0, 90, n)

function plotit(df1, m, df2, file)
    fig = band(df2.elevation, df2.lower, df2.upper, color = :gray, axis = (; title = file, xlabel = "Elevation (°)", ylabel = "Probability", limits = ((0, 90),(-0.01, 0.5))))
    lines!(df2.elevation, df2.prediction, color = :white)
    scatter!(df1.elevation, df1.p, color = :black, markersize = 20)
    !ismissing(m) && text!(0, 0.5, text = string("BIC = ", round(Int, GLM.bic(m))), align = (:left, :top), offset = (10, -10))
    save("$file.png", fig)
end

function prediction(chain, x)
    p = get_params(chain[200:end, :, :])
    targets = my_link.(p.ϵ', x, p.a', p.b')
    ps = quantile.(eachrow(targets), Ref([0.025, 0.5, 0.975]))
    DataFrame(elevation = x, lower = getindex.(ps, 1), prediction = getindex.(ps, 2), upper = getindex.(ps, 3))
end

df = CSV.read("data.csv", DataFrame)

m = glm(@formula(lap ~ elevation), df, Binomial())

n = 100
newdf = DataFrame(elevation = range(0, 100, n), lap = falses(n))
plu = predict(m, newdf, interval = :confidence)
newdf.prediction = disallowmissing(plu.prediction)
newdf.lower = disallowmissing(plu.lower)
newdf.upper = disallowmissing(plu.upper)

df2 = combine(groupby(df, :elevation), :lap => count => :lap, nrow => :n)
transform!(df2, [:lap, :n] => ByRow(/) => :p)

plotit(df2, m, newdf, "glm with logit")

using Turing
using StatsFuns: logistic

@model function bmodel(xs, ys)
    α ~ Normal(0, 3)
    β ~ Normal(0, 3)
    ps = logistic.(α .+ β .* xs)
    ys ~ product_distribution(Bernoulli.(ps))
    # for i in 1:n
    #     p = logistic(α + β*xs[i])
    #     ys[i] ~ Bernoulli(p)
    # end
end
posterior = bmodel(df.elevation, df.lap)
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


p = get_params(chain[200:end, :, :])
targets = logistic.(p.α' .+ p.β' .* elevation)
ps = quantile.(eachrow(targets), Ref([0.025, 0.5, 0.975]))
newdf = DataFrame(elevation = elevation, lower = getindex.(ps, 1), prediction = getindex.(ps, 2), upper = getindex.(ps, 3))

plotit(df2, missing, newdf, "bayesian with logit")


# x = -10:100
# mylogistic(x, x₀, k, L) = L/(1 + exp(-k*(x - x₀)))
# lines(x, mylogistic.(2x, 90, 1/30, 180)/2; axis = (; limits = (nothing, (0, 90))))
# ablines!(0, 1)

# transθ(x) = mylogistic.(x, 45, 1/20, 90)
# transϵ(x) = mylogistic.(x, 90, 1/30, 180)
# trans(x) = mylogistic.(x, 45, 1/20, 90)

function my_link(ϵ, θ, a, b) 
    θ = θ < 0 ? zero(θ) : θ > 90 ? 90one(θ) : θ
    ϵ = ϵ < 0 ? zero(ϵ) : ϵ > 180 ? 180one(ϵ) : ϵ
    (2atand(tand(ϵ/2)/cosd(θ)) - ϵ)/(180 - ϵ)*a + b
end


@model function bmodel(xs, ys)
    # Priors
    ϵ ~ truncated(Normal(6, 10), 0.1, 179.9)
    # ϵ ~ Uniform(1, 179)
    a ~ truncated(Normal(0.5, 1), 0, 1)
    # a ~ Uniform(0.01, 0.99)
    b ~ truncated(Normal(0, 1), 0, 1)
    # b ~ Uniform(0.01, 0.99)
    ps = my_link.(ϵ, xs, a, b)
    if any(p -> p < 0 || p > 1, ps)
        @addlogprob! -Inf
        # Exit the model evaluation early
        return nothing
    end
    # Likelihood
    ys ~ product_distribution(Bernoulli.(ps))
    return nothing
end
df3 = transform(df, :elevation => ByRow(x -> x == 90 ? 89.9 : x), renamecols = false)
posterior = bmodel(df3.elevation, df3.lap)
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

newdf = prediction(chain, elevation)
plotit(df2, missing, newdf, "bayesian with new link")





fun(ϵ) = θ -> (2atand(tand(ϵ/2)/cosd(θ)) - ϵ)/(180 - ϵ)
ifun(ϵ) = p -> asecd(cotd(ϵ/2)*tand(1/2*(-ϵ*p + 180p + ϵ)))

bic_glm(ϵ) = glm(@formula(lap ~ p2), transform(df, :elevation => ByRow(fun(ϵ)) => :p2) , Binomial(), IdentityLink())
optimal_bic = optimize(GLM.bic ∘ bic_glm, 0.01, 179.99)
ϵ = optimal_bic.minimizer
m = bic_glm(ϵ)

newdf = DataFrame(p2 = range(0, 1, n), lap = falses(n))
plu = predict(m, newdf, interval = :confidence)
newdf.prediction = disallowmissing(plu.prediction)
newdf.lower = disallowmissing(plu.lower)
newdf.upper = disallowmissing(plu.upper)

transform!(newdf, :p2 => ByRow(ifun(ϵ)) => :elevation)

plotit(df2, m, newdf, "glm with new link")


