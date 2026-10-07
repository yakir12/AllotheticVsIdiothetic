using Dates
using CSV, DataFramesMeta, GLMakie, Chain, MixedModels, AlgebraOfGraphics
using GLM
using Statistics
using StandardizedPredictors
using MixedModelsSim, Random

df = @chain "data.csv" begin
    CSV.read(DataFrame, select = ["condition", "beetle_id", "roll_duration_min", "#_dances"])
    @rename :id = :beetle_id :duration = :roll_duration_min :n = $"#_dances"
    @transform :duration = Int.((:duration .- Time(0)) ./ Nanosecond(Second(1)))
end

frm = @formula(n ~ 1 + condition + duration + (1|id))
contrasts = Dict(:duration => Center())
m = glmm(frm, df, Poisson(); contrasts) 

newdf = allcombinations(DataFrame, :condition => unique(df.condition), :duration => range(extrema(df.duration)..., 100))
newdf.n .= 0
newdf.id .= "new beetle"
newdf.n .= predict(m, newdf)

fig = (data(newdf) * mapping(:duration, :n, color = :condition) * visual(Lines) + 
 data(df) * mapping(:duration, :n, color = :condition) * visual(Scatter)) |> draw()

#Beetles dance significantly more often on natural terrain than smooth terrain overall (P < 0.001), and that gap is largest for short rolls (~4x at the shortest observed durations, 2 minutes) and shrinks toward no detectable difference for the longest rolls, 10 minutes (significant terrain type × duration interaction with P = 0.0049). We used a GLMM with a Poisson distribution family (n = 43).
#
# Beetles dance significantly more often on natural terrain than smooth terrain (beetles dance roughly twice as often on natural terrain; P < 0.001). We used a GLMM with a Poisson distribution family (n = 43).
#
#

sim = parametricbootstrap(MersenneTwister(12321), 1000, m)

power_table(sim, 0.001)

# power_table only reports the point proportion (successes / nsim); to attach a
# confidence interval we need the same counts it computes internally, then treat
# "proportion of simulations reaching significance" as a binomial proportion and
# use a Wilson score interval (better-behaved than the normal approximation when
# power is close to 0 or 1, which is common in these tables).
function power_table_ci(sim::MixedModels.MixedModelBootstrap, alpha=0.05; z=1.96)
    nsim = length(sim.objective)
    dd = Dict{Symbol,Int}()
    for row in sim.coefpvalues
        val = get!(dd, row.coefname, 0)
        dd[row.coefname] = val + (row.p < alpha)
    end
    return [begin
        p̂ = x / nsim
        denom = 1 + z^2 / nsim
        centre = (p̂ + z^2 / (2nsim)) / denom
        halfwidth = z / denom * sqrt(p̂ * (1 - p̂) / nsim + z^2 / (4nsim^2))
        (; coefname = string(key), power = p̂, lower = max(0.0, centre - halfwidth), upper = min(1.0, centre + halfwidth))
    end for (key, x) in dd]
end

DataFrame(power_table_ci(sim, 0.001))
