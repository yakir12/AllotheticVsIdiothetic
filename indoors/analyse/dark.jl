using AlgebraOfGraphics, GLMakie, CairoMakie
using GLM
# using MultivariateStats
using DataFramesMeta, Chain
using HypothesisTests
using GeometryBasics
# using StatsBase, Graphs
using Dates, LinearAlgebra, Statistics, Random
using CSV, DataFrames, CameraCalibrations
using Interpolations, StaticArrays, Dierckx, CoordinateTransformations, Rotations
using OhMyThreads
using LsqFit
using Optim
using CategoricalArrays
using Distributions
using IntervalSets
using QuadGK
using BetaRegression
using DimensionalData

GLMakie.activate!()

include("general_functions.jl")
include("dark_functions.jl")

const l = 50

output = "dark"
if isdir(output)
    rm.(readdir(output; join = true))
else
    mkdir(output)
end

const results_dir = "../track_calibrate/tracks and calibrations"

runs = @chain joinpath(results_dir, "runs.csv") begin
    CSV.read(DataFrame)
    @select Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size)
    @subset :light .≠ "shift"
    @transform :tij_file = joinpath.(results_dir, :tij_file)
    @transform :location = "Lund"
end
calibs = @chain joinpath(results_dir, "calibs.csv") begin
    CSV.read(DataFrame)
    @transform :calibration_file = joinpath.(results_dir, :calibration_id)
    @transform :rectify = get_calibration.(:calibration_file)
    @select Cols(:calibration_id, :rectify)
end
leftjoin!(runs, calibs, on = :calibration_id)
@chain runs begin
    @select! Not(:calibration_id)
    @transform! :dance_spontaneous = .!ismissing.(:spontaneous_end)
    @transform! :condition = string.(:light, " ", :dance_by, " ", :at_run)
    @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
    @rename! :intervention = :poi
    @transform! :pixels = get_tij.(:tij_file)
    @transform! :pixels = remove_stops.(:pixels)
    # @aside @assert !any(has_stops, _.pixels) "some stops remain?!" # convert to a test
    @rtransform! :xy = :rectify.(:pixels)
    @aside @chain _ begin 
        @subset(:dance_by .≠ "no"; view = true)
        @transform! :jump = glue_intervention!.(:xy, :intervention)
    end
    @aside @chain _ begin 
        @subset(:light .== "remain"; view = true)
        @transform! :intervention = impute_poi_time.(:xy)
    end
    @transform! :spontaneous_end = passmissing(tosecond).(:spontaneous_end)
    @transform! :poi = coalesce.(:spontaneous_end, :intervention)
    disallowmissing!(:poi)
    @select! Not(:intervention)
    @transform! :xy = remove_loops.(:xy)
    @transform! :xy = sparseify.(:xy)
    @transform! :smooth = smooth.(:xy)
    @transform! :centered2start = center2start.(:smooth)
    @transform! :cropped = cropto.(:centered2start, l)
    @transform! :rotated2poi = rotate2poi.(:cropped, :poi)
    @transform! :centered2poi_and_cropped = center2poi_and_crop.(:rotated2poi, :poi)
    transform!([:pixels, :xy, :smooth, :centered2start, :cropped, :rotated2poi, :centered2poi_and_cropped] .=> ByRow(parent), renamecols = false)
end

############ plot the tracks to check validity

fig = (pregrouped(runs.smooth => first => "X (cm)", runs.smooth => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines; color = :red) + pregrouped(runs.xy => first => "X (cm)", runs.xy => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines)) |> draw(; axis = (; width = 400, height = 400, limits = ((-l, l), (-l, l))));
save(joinpath(output, "overview_dark.png"), fig)

################################################### 10 random tracks

n = 10
df = @chain runs begin
    @subset norm.(last.(:cropped)) .≈ l
    @groupby :dance_by
    @transform :n = 1:length(:dance_by)
    @subset :n .≤ n
end
@assert all(==(n), combine(groupby(df, :dance_by), nrow).nrow)

fig = pregrouped(df.cropped => first => "X (cm)", df.cropped => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, l)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure1.png"), fig)

# fig = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
# for ax in fig.figure.content 
#     if ax isa Axis
#         for r  in (30, 50)
#             lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
#         end
#     end
# end
# save(joinpath(output, "figure1a.png"), fig)


######################## Figure 2

fig = pregrouped(runs.rotated2poi => first => "X (cm)", runs.rotated2poi => last => "Y (cm)", col = runs.condition, color = runs.y2025) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure2.png"), fig)

################################################### Figure 3


df = @subset runs :at_run .== 1
# df = deepcopy(runs)
nr = 3
l1 = floor(Int, minimum(norm ∘ last, df.centered2poi_and_cropped))
change minimum to 5!!!!
rl = range(5, l1, nr)
transform!(df, :centered2poi_and_cropped => ByRow(xy -> get_exit_angle.(Ref(xy), rl)) => :θs)
select!(df, Cols(:light, :dance_by, :θs, :centered2poi_and_cropped))
df.light = categorical(df.light)
levels!(df.light, ["remain", "dark"])
df.dance_by = categorical(df.dance_by)
levels!(df.dance_by, ["no", "hold", "disrupt"])

# R = 60
# fig = pregrouped(df.centered2poi_and_cropped => first => "X (cm)", df.centered2poi_and_cropped => last => "Y (cm)", col = df.light => renamer("remain" => "Lights on", "dark" => "Lights off"), row = df.dance_by => renamer("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball")) * visual(Lines) |> draw(; axis=(width = 300, height = 300, limits = ((-R, R), (-R, R))))

# save(joinpath(output, "figure2a.png"), fig)

light_rename = Dict("remain" => "Lights on", "dark" => "Lights off")
dance_by_rename = Dict("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball")
fig = Figure()
for ((light, dance_by), g) in pairs(groupby(df, [:light, :dance_by]))
    i = findfirst(==(dance_by), levels(df.dance_by))
    j = findfirst(==(light), levels(df.light))
    ax = PolarAxis(fig[i, j], rlimits = (0, 44), width = 300, height = 300)
    for xy in g.centered2poi_and_cropped
        θ = splat(atan).(xy) .+ π/2
        r = norm.(xy)
        lines!(ax, θ, r, color = :black)
    end
    if i == 1
        Label(fig[0,j], light_rename[string(light)], font = :bold)
    end
    if j == 2
        Label(fig[i,3], dance_by_rename[string(dance_by)], rotation = -π/2, font = :bold)
    end
end
resize_to_layout!(fig)

save(joinpath(output, "figure2b.png"), fig)

##################################

df = deepcopy(runs)
nr = 3
l1 = floor(Int, minimum(norm ∘ last, df.centered2poi_and_cropped))
rl = range(1e-3, l1, nr)
transform!(df, :centered2poi_and_cropped => ByRow(xy -> get_exit_angle.(Ref(xy), rl)) => :θs)
select!(df, Cols(:light, :dance_by, :at_run, :θs, :centered2poi_and_cropped))
# df.light = categorical(df.light)
# levels!(df.light, ["remain", "dark"])
# df.dance_by = categorical(df.dance_by)
# levels!(df.dance_by, ["no", "hold", "disrupt"])
@transform! df :grp = string.(:light, :dance_by, :at_run)
df.grp = categorical(df.grp)
levels!(df.grp, [
 "remainno1",
 "darkno1",
 "darkno10",
 "darkhold1",
 "darkdisrupt1"
])

fig = Figure()
for (i, label) in enumerate(["Light", "Manipulation", "At run"])
     Label(fig[i,0], label)
end
for (j, (k, g)) in enumerate(pairs(groupby(df, :grp)))
    for (i, col) in enumerate([:light, :dance_by, :at_run])
        Label(fig[i, j], string(g[1,col]))
    end
    ax = PolarAxis(fig[4, j], rlimits = (0, 44), width = 300, height = 300)
    for xy in g.centered2poi_and_cropped
        θ = splat(atan).(xy) .+ π/2
        r = norm.(xy)
        lines!(ax, θ, r, color = :black)
    end
end
resize_to_layout!(fig)

save(joinpath(output, "figure2c.png"), fig)

# select!(df, Not(:centered2poi_and_cropped, :grp))
# select!(df, Cols(:condition, :θs))
df.r .= Ref(rl)



lightdf = @chain df begin
    @subset :dance_by .== "no" :at_run .== 1
    @select Not(:dance_by, :at_run)
end
fm = @formula(mean_resultant_vector ~ light + r)
newlight = combine(groupby(lightdf, :light), :r => first ∘ first => :r)
nr2 = 100
rl2 = range(extrema(rl)..., nr2)
newlight.r .= Ref(rl2)
newlight = flatten(newlight, :r)
newlight.mean_resultant_vector .= 0.0
pc, c = bootstrap(lightdf, fm, newlight, [:light])
select!(pc, Not(Symbol("(Precision)")))

y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newlight.lower .= getindex.(y, 1)
newlight.mean_resultant_vector .= getindex.(y, 2)
newlight.upper .= getindex.(y, 3)
pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)
pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
df2 = stack(pvalues, Not(:what), variable_name = :source)
pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])

show(pvalues, show_row_number=false, eltypes=false)

fig = data(newlight) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :light => renamer("remain" => "Lights on", "dark" => "Lights off")) * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3a.png"), fig)

df2 = combine(groupby(newlight, :light), nrow)
@transform! df2 :r̄ = critical_r.(:nrow)

fig = (data(newlight) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :light => renamer("remain" => "Lights on", "dark" => "Lights off")) * visual(LinesFill) + data(df2) * mapping(:r̄) * visual(HLines)) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3a1.png"), fig)

atrun = @chain df begin
    @subset :dance_by .== "no" :light .== "dark"
    @select Not(:dance_by, :light)
end
fm = @formula(mean_resultant_vector ~ at_run + r)
newrun = combine(groupby(atrun, :at_run), :r => first ∘ first => :r)
nr2 = 100
rl2 = range(extrema(rl)..., nr2)
newrun.r .= Ref(rl2)
newrun = flatten(newrun, :r)
newrun.mean_resultant_vector .= 0.0
pc, c = bootstrap(atrun, fm, newrun, [:at_run])
select!(pc, Not(Symbol("(Precision)")))

y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newrun.lower .= getindex.(y, 1)
newrun.mean_resultant_vector .= getindex.(y, 2)
newrun.upper .= getindex.(y, 3)
pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)
pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
df2 = stack(pvalues, Not(:what), variable_name = :source)
pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])

show(pvalues, show_row_number=false, eltypes=false)

fig = data(newrun) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :at_run => nonnumeric => "At run #") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3b.png"), fig)

df2 = combine(groupby(newrun, :at_run), nrow)
@transform! df2 :r̄ = critical_r.(:nrow)

fig = (data(newrun) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :at_run => nonnumeric => "At run #") * visual(LinesFill) + data(df2) * mapping(:r̄) * visual(HLines)) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3b1.png"), fig)


dance = @chain df begin
    @subset :at_run .== 1 :light .== "dark"
    @select Not(:at_run, :light)
end
fm = @formula(mean_resultant_vector ~ dance_by + r)
newdance = combine(groupby(dance, :dance_by), :r => first ∘ first => :r)
nr2 = 100
rl2 = range(extrema(rl)..., nr2)
newdance.r .= Ref(rl2)
newdance = flatten(newdance, :r)
newdance.mean_resultant_vector .= 0.0
pc, c = bootstrap(dance, fm, newdance, [:dance_by])
select!(pc, Not(Symbol("(Precision)")))

y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdance.lower .= getindex.(y, 1)
newdance.mean_resultant_vector .= getindex.(y, 2)
newdance.upper .= getindex.(y, 3)
pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)
pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
df2 = stack(pvalues, Not(:what), variable_name = :source)
pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])

show(pvalues, show_row_number=false, eltypes=false)

fig = data(newdance) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :dance_by => renamer("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball") => "Dance induced by") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3c.png"), fig)


df2 = combine(groupby(newdance, :dance_by), nrow)
@transform! df2 :r̄ = critical_r.(:nrow)

fig = (data(newdance) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :dance_by => renamer("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball") => "Dance induced by") * visual(LinesFill) + data(df2) * mapping(:r̄) * visual(HLines)) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3c1.png"), fig)

newlight.at_run .= 1
newlight.dance_by .= Ref("no")
newlight.test .= Ref("light")
newrun.light .= Ref("dark")
newrun.dance_by .= Ref("no")
newrun.test .= Ref("run")
newdance.light .= Ref("dark")
newdance.at_run .= 1
newdance.test .= Ref("manipulation")

df3 = vcat(newlight, newrun, newdance)
@transform! df3 :grp = string.(:light, " ", :at_run, " ", :dance_by)

fig = data(df3) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, col = :test, color = :grp) * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

save(joinpath(output, "figure3d.png"), fig)


###### compare "full" model: mean_resultant_vector ~ at_run + light + dance_by + r

fm = @formula(mean_resultant_vector ~ at_run + light + dance_by + r)
newdf = combine(groupby(df, [:at_run, :light, :dance_by]), :r => first ∘ first => :r)
nr2 = 100
rl2 = range(extrema(rl)..., nr2)
newdf.r .= Ref(rl2)
newdf = flatten(newdf, :r)
newdf.mean_resultant_vector .= 0.0
pc, c = bootstrap(df, fm, newdf, [:at_run, :light, :dance_by])
select!(pc, Not(Symbol("(Precision)")))

y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdf.lower .= getindex.(y, 1)
newdf.mean_resultant_vector .= getindex.(y, 2)
newdf.upper .= getindex.(y, 3)
pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)
pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
df2 = stack(pvalues, Not(:what), variable_name = :source)
pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])

show(pvalues, show_row_number=false, eltypes=false)

@transform! newdf :grp = string.(:light, " ", :at_run, " ", :dance_by)

fig = data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, col = :grp) * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))

fig = data(newdf) * mapping(:r, :lower, :upper, col = :grp) * visual(Band; alpha = 0.5) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))


tracks2plot = combine(groupby(df, [:light, :dance_by, :at_run]), :centered2poi_and_cropped => Ref => :centered2poi_and_cropped)
stats2plot = combine(groupby(newdf, [:light, :dance_by, :at_run]), [:r, :mean_resultant_vector, :lower, :upper] .=> Ref .=> [:r, :mean_resultant_vector, :lower, :upper])
toplot = leftjoin(tracks2plot, stats2plot, on = [:light, :dance_by, :at_run])
sort!(toplot, order(:upper, by = first, rev = true))

@transform! toplot :r̄ = critical_r.(length.(:centered2poi_and_cropped))

fig = Figure()
for (i, label) in enumerate(["Light", "Manipulation", "At run"])
     Label(fig[i,0], label)
end
for (j, row) in enumerate(eachrow(toplot))
# for (j, (k, g)) in enumerate(pairs(groupby(df, [:at_run, :light, :dance_by])))
    for (i, col) in enumerate([:light, :dance_by, :at_run])
        Label(fig[i, j], string(row[col]))
    end
    ax = PolarAxis(fig[4, j], rlimits = (0, 44), width = 300, height = 300)
    for xy in row.centered2poi_and_cropped
        θ = splat(atan).(xy) .+ π/2
        r = norm.(xy)
        lines!(ax, θ, r, color = :black)
    end
    ax = Axis(fig[5, j], limits = ((0, l1), (0, 1)), width = 300, height = 300, ylabel = "Mean resultant vector length")
    linesfill!(ax, row.r, row.mean_resultant_vector, lower = row.lower, upper = row.upper)
    if j ≠ 1
        hideydecorations!(ax, grid = false, minorgrid = false)
    end
    hlines!(ax, row.r̄, color = :grey)
end
Label(fig[6, 1:5], "Radius (cm)")
resize_to_layout!(fig)

save(joinpath(output, "figure3e.png"), fig)


###### straightness

df = @chain runs begin
    @transform :path_length = path_length_at.(:centered2poi_and_cropped, l1) .- l1
    @transform :grp = string.(:light, " ", :at_run, " ", :dance_by)
end

data(df) * mapping(:grp, :path_length) * visual(RainClouds) |> draw()

########### each separate

lightdf = @chain df begin
    @subset :dance_by .== "no" :at_run .== 1
    @select Not(:dance_by, :at_run)
end
fm = @formula(path_length ~ light)
m = glm(fm, lightdf, Gamma())

atrun = @chain df begin
    @subset :dance_by .== "no" :light .== "dark"
    @select Not(:dance_by, :light)
end
fm = @formula(path_length ~ at_run)
m = glm(fm, atrun, Gamma())

dance = @chain df begin
    @subset :at_run .== 1 :light .== "dark"
    @select Not(:at_run, :light)
end
fm = @formula(path_length ~ dance_by)
m = glm(fm, dance, Gamma())

#### all together

m = glm(@formula(path_length ~ light + at_run + dance_by), df, Gamma())

# newdf = combine(groupby(df, [:at_run, :light, :dance_by]), :path_length => first => :path_length)
# newdf = hcat(newdf, predict(m, newdf, interval = :confidence))
# @transform! newdf :μ = l1 ./ (:prediction .+ l1)
# @transform! newdf :upper = l1 ./ (:upper .+ l1)
# @transform! newdf :lower = l1 ./ (:lower .+ l1)


# sdjfghsdjkfhlsfhj
#
# df = DataFrame(y = rand(10), x = repeat(["c", "b"], 5))
# fm = @formula(y ~ x)
# m = BetaRegression.fit(BetaRegressionModel, fm, df)
#
#
# # df.condition = categorical(df.condition)
# # levels!(df.condition, ["remain", "no", "hold", "disrupt"])
#
#
# # fig = Figure()
# # for (i, (k, g)) in enumerate(pairs(groupby(df, :condition)))
# #     ax = PolarAxis(fig[i,1])
# #     # _df = flatten(g, [:θs, :r])
# #     # scatter!(ax, _df.θs, _df.r)
# #     for row in eachrow(g)
# #         scatter!(ax, row.θs, row.r)
# #     end
# # end
#
# newdf = combine(groupby(df, [:light, :dance_by, :at_run]), :r => first ∘ first => :r)
# nr2 = 100
# rl2 = range(extrema(rl)..., nr2)
# newdf.r .= Ref(rl2)
# newdf = flatten(newdf, :r)
# newdf.mean_resultant_vector .= 0.0
#
# fm = @formula(mean_resultant_vector ~ light + r + dance_by + at_run)
#
# # conditions = levels(df.condition)
# # nr2 = 100
# # rl2 = range(extrema(rl)..., nr2)
# # newdf = DataFrame(r = repeat(rl2, outer = length(conditions)), condition = repeat(conditions, inner = nr2), mean_resultant_vector = zeros(nr2*length(conditions)))
# # fm = @formula(mean_resultant_vector ~ r*condition)
#
# pc, c = bootstrap(df)
# select!(pc, Not(Symbol("(Precision)")))
#
# y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
# newdf.lower .= getindex.(y, 1)
# newdf.mean_resultant_vector .= getindex.(y, 2)
# newdf.upper .= getindex.(y, 3)
# pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)
# pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
# df2 = stack(pvalues, Not(:what), variable_name = :source)
# pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])
#
#
# df2 = combine(groupby(df, [:light, :dance_by]), nrow)
# @transform! df2 :r̄ = critical_r.(:nrow)
#
#
# # n = nrow(df)
# # df2 = flatten(df, [:θs, :r])
# # df3 = combine(groupby(df2, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
# # df3.condition = categorical(df3.condition)
# # levels!(df3.condition, ["remain", "no", "hold", "disrupt"])
# # m = BetaRegression.fit(BetaRegressionModel, fm, df3)
#
# # + data(df2) * mapping(:r̄, col = :light, color = :dance_by => renamer("no" => "No disruption", "hold" => "Hold", "disrupt" => "Taken off ball") => "Light") * visual(HLines) 
# #, yticks = ([0; sort(df2.r̄); 0.5; 1], ["0", "", "", "", "", "0.5", "1"]
#
# fig = data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, group = [:light, :dance_by], col = :light => renamer("remain" => "Lights on", "dark" => "Lights off"), color = :dance_by => renamer("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball") => "Dance induced by") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l1), (0, 1))))
#
#
#
# save(joinpath(output, "figure3.png"), fig)
#
# show(pvalues, show_row_number=false, eltypes=false)
#
