using DataFramesMeta, CSV, Chain
using LinearAlgebra
using GLMakie, AlgebraOfGraphics
using BetaRegression
using Random
using StandardizedPredictors
using Statistics
using Distributions

df = @chain "data.csv" begin
    CSV.read(DataFrame)
    @subset :elevation .> 40
    @select :id = string.(:elevation, "_", :nr_in_cond) :elevation :exit = :rel_bearing
    dropmissing(:exit)
    @transform :exit = deg2rad.(mod.(90 .- :exit, 360))
    # groupby(:id)
    # @transform :exit = :exit .- angle(mean(exp, :exit * 1im))
end

get_rvalue(x) = norm(mean(cis, x))

df2 = combine(groupby(df, [:elevation, :id]), :exit => get_rvalue => :rvalue)

# @chain df2 begin
#     subset(:rvalue => ByRow(==(1)), view = true)
#     @transform! :rvalue = maximum(filter(<(1), df2.rvalue))
# end

fig = data(df2) * mapping(:elevation, :rvalue) * visual(BoxPlot) |> draw()


frm = @formula(rvalue ~ 1 + elevation)
m = BetaRegression.fit(BetaRegressionModel, frm, df2)#; contrasts)

tbl = coeftable(m)
row = (; Pair.(Symbol.(tbl.rownms), tbl.cols[tbl.pvalcol])...)  # Extract p-values



# MixedModelsSim.parametricbootstrap only has a method for MixedModels.MixedModel, so it
# can't be used on a BetaRegressionModel directly. We hand-roll the same idea instead:
# simulate new rvalue from Beta(μᵢ·φ, (1-μᵢ)·φ) using this model's own fitted means μ and
# precision φ (i.e. treat the fitted model as the assumed-true population), refit, and
# record how often each coefficient's p-value clears the significance threshold.
function beta_parametric_bootstrap(rng, nsim, m, df, frm)
    μ = fitted(m)
    ϕ = precision(m)
    dfsim = copy(df)
    results = NamedTuple[]
    for i in 1:nsim
        dfsim.rvalue = rand.(rng, Beta.(μ .* ϕ, (1 .- μ) .* ϕ))
        msim = BetaRegression.fit(BetaRegressionModel, frm, dfsim)
        tblsim = coeftable(msim)
        for (name, p) in zip(tblsim.rownms, tblsim.cols[tblsim.pvalcol])
            push!(results, (; iter = i, coefname = name, p = p))
        end
    end
    return results
end

sim = beta_parametric_bootstrap(MersenneTwister(12321), 1000, m, df2, frm)

# Same Wilson-score-interval logic as used for the GLMM power tables, just built on top
# of our own (coefname, p) rows instead of a MixedModels.MixedModelBootstrap.
function power_table_ci(results, alpha=0.05; z=1.96)
    dd = Dict{String,Tuple{Int,Int}}()
    for row in results
        sig, tot = get(dd, row.coefname, (0, 0))
        dd[row.coefname] = (sig + (row.p < alpha), tot + 1)
    end
    return [begin
        p̂ = sig / tot
        denom = 1 + z^2 / tot
        centre = (p̂ + z^2 / (2tot)) / denom
        halfwidth = z / denom * sqrt(p̂ * (1 - p̂) / tot + z^2 / (4tot^2))
        (; coefname = key, nsim = tot, power = p̂, lower = max(0.0, centre - halfwidth), upper = min(1.0, centre + halfwidth))
    end for (key, (sig, tot)) in dd]
end

DataFrame(power_table_ci(sim, 0.05))
