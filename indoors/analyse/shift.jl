using AlgebraOfGraphics, GLMakie, CairoMakie

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

GLMakie.activate!()

include("minimal_functions.jl")

output = "figures"
mkpath("figures")

const results_dir = "../track_calibrate/tracks and calibrations"

function combine_factors(light, induced, run)
    induced = induced ? " induced" : ""
    run = run > 1 ? " $run" : ""
    string(light, induced, run)
end

runs = CSV.read(joinpath(results_dir, "runs.csv"), DataFrame)
@rsubset!(runs, :light == "shift")
select!(runs, Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size))
calibs = CSV.read(joinpath(results_dir, "calibs.csv"), DataFrame)
transform!(calibs, :calibration_id => ByRow(get_calibration) => :rectify)
select!(calibs, Cols(:calibration_id, :rectify))
leftjoin!(runs, calibs, on = :calibration_id)
select!(runs, Not(:calibration_id))
@chain runs begin
    @rename! :intervention = :poi
    @rtransform! $AsTable = get_tij(:tij_file)
    @rtransform! $AsTable = remove_stops!(:t, :ij)
    @transform! :xy = trectify(:rectify, :ij)
    @aside @chain _ begin 
        @subset(:dance_induced; view = true)
        # @aside pregrouped(_.xy => first, _.xy => last)  * visual(Lines) * pregrouped(layout = _.run_id => nonnumeric) |> draw(figure = (;size = (1000, 1000)), axis = (;aspect = DataAspect())) |> save("before.png")
        @rtransform! :jump = glue_intervention!(:xy, :t, :intervention)
    end
    @aside @chain _ begin 
        @subset(:light .== "remain"; view = true)
        @rtransform! :poi = impute_poi_time(:t, :xy)
    end
    @rtransform! :spontaneous_end = passmissing(tosecond)(:spontaneous_end)
    @rtransform! :poi = coalesce(:spontaneous_end, :intervention, :poi)
    disallowmissing!(:poi)
    @rtransform! :poi_index = something(findfirst(≥(:poi), :t), length(:t))
    @rtransform! $AsTable = remove_loops!(:t, :xy)
    @rtransform! $AsTable = sparseify(:t, :xy, :poi)
    @rtransform! $AsTable = smooth(:t, :xy)
    @rtransform! :dance_spontaneous = !ismissing(:spontaneous_end)
    @rtransform! :condition = combine_factors(:light, :dance_induced, :at_run)
    @rtransform! :rot = get_rotation(:xys[:poi_index])
    @rtransform! :l = get_pathlength(:xys)
    @rtransform! $AsTable = get_turn_profile(:t, :spl, :poi_index)
    @rtransform! $AsTable = fit_logistic(:lθ, :θ)
    @rtransform! :θs = logistic.(:t, Ref(:ks))
    @transform! :k = getindex.(:ks, 2)
    @rtransform! :Δθ = abs(wrap2pi(-(reverse(extrema(:θs))...)))
    @rtransform! :Δpoi = :ks[3] - :poi
end

################################################### 10 random tracks

df = @chain runs begin
    @subset norm.(last.(:xys)) .> 50
    @rtransform :xys = cropto(:t, :spl, :trans, 50)
    @rtransform :xysr = :rot.(:xys)
    @groupby :dance_induced
    @transform :n = 1:length(:dance_induced)
    @subset :n .≤ 10
end
@assert all(==(10), combine(groupby(df, :dance_induced), nrow).nrow)

fig = pregrouped(df.xys => first, df.xys => last, col = df.dance_induced => renamer(true => "Induced", false => "Not")) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
save(joinpath(output, "figure4.png"), fig)

fig = pregrouped(df.xysr => first, df.xysr => last, col = df.dance_induced => renamer(true => "Induced", false => "Not")) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
save(joinpath(output, "figure4a.png"), fig)


################################################### raw data for turning angles

df = @chain runs begin
    @rtransform :lθ = :lθ .- first(:lθ)
    @rtransform :θ = :θ .- first(:θ)
    @rtransform :θ = mean(:θ) < 0 ? -(:θ) : :θ
end

fig = pregrouped(df.lθ => "Path length (cm)", df.θ => rad2deg => "Turning", col = df.dance_induced => renamer(true => "Induced", false => "Not")) * visual(Lines) |> draw(; axis = (; width = 400, height = 400, ytickformat = "{:n}°", yticks = 0:90:270, limits = ((nothing, 10), (-80, nothing))))

save(joinpath(output, "figure5.png"), fig)


# function loops(xy)
#     inds, _ = self_intersections(Point2f.(xy))
#     ranges = splat(UnitRange).(Iterators.partition(inds, 2))
#     length.(ranges)
# end
#
# l = map(loops, df.xy)
# l = reduce(vcat, l)

# df = runs
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect(), limits = ((-100, 100), (-100, 100)))
# sl = SliderGrid(fig[3, 1], (range = 1:nrow(df), ))
# i = only(sl.sliders).value
# poi1 = lift(i) do i
#     df.xy[i][1:df.poi_index[i]]
# end
# poi2 = lift(i) do i
#     df.xys[i][1:df.poi_index[i]]
# end
# tra1 = lift(i) do i
#     df.xy[i][df.poi_index[i]:end]
# end
# tra2 = lift(i) do i
#     df.xys[i][df.poi_index[i]:end]
# end
# sl = SliderGrid(fig[4, 1], (range = 1:minimum(length, df.t), ))
# j = only(sl.sliders).value
# xy1 = lift(i, j) do i, j
#     df.xy[i][j]
# end
# xy2 = lift(i, j) do i, j
#     df.xys[i][j]
# end
# lines!(ax, poi1, color = :black)
# lines!(ax, poi2, color = :red)
# lines!(ax, tra1, color = :blue)
# lines!(ax, tra2, color = :green)
# scatter!(ax, xy1)
# scatter!(ax, xy2)
# ax1 = Axis(fig[2,1], limits = (nothing, (-90, 360)), yticks = -90:90:360)
# turn1 = lift(i) do i
#     Point2f.(df.tθ[i][1:df.poi_index[i]], rad2deg.(df.θ[i][1:df.poi_index[i]]))
# end
# turn2 = lift(i) do i
#     Point2f.(df.tθ[i][df.poi_index[i]:end], rad2deg.(df.θ[i][df.poi_index[i]:end]))
# end
# txy = lift(i, j) do i, j
#     Point2f(df.tθ[i][j], rad2deg(df.θ[i][j]))
# end
# lines!(ax1, turn1, color = :red)
# lines!(ax1, turn2, color = :green)
# scatter!(ax1, txy)
# on(turn1) do t
#     limits!(ax1, -1, first(last(t)) + 1, -90, 360)
# end

# n = 100
# x = range(1, 100, n)
#
# L = 10 # total height, from lowest to highest
# k = -1 # slope, lower is slower
# x₀ = 10 # x value for the midpoint
# y₀ = -3 # y displacement
# y = logistic(x, [L, k, x₀, y₀]) .+ 0.1randn(n)
# lines(x, y)
#
# p0 = guess_logistic(x, y)
# lb = Float64[1, -10, -100, -100]
# ub = Float64[100, 10, 100, 100]
# fit = curve_fit(logistic, x, y, p0, lower = lb, upper = ub)
# yl = logistic(x, fit.param)
# lines!(x, yl)

df1 = @rsubset(runs, (:dance_induced))
# @rtransform! df1 :tθ =  :tθ .- first(:tθ)

fig = (pregrouped(df1.tθ, df1.θ => rad2deg)  * visual(Lines) + pregrouped(df1.t, df1.θs => rad2deg)  * visual(Lines; color = :red)) * pregrouped(layout = string.(df1.run_id, " ", round.(df1.k, digits = 4))) |> draw(axis = (; yticks = -360:180:360), facet = (; linkxaxes = :none, linkyaxes = :all))
# display(fig)

fig = pregrouped([ts .- first(ts) for ts in runs.t], runs.θs => rad2deg, col = runs.dance_induced)  * visual(Lines) * pregrouped(group = string.(runs.run_id, " ", round.(runs.k, digits = 4))) |> draw(axis = (; yticks = -360:180:360))

fig = data(runs) * mapping(:Δθ => rad2deg, :k => abs, color = :dance_induced, marker = :dance_induced) * visual(Scatter) |> draw()
display(fig)

data(runs) * mapping(:logistic_rsquare, :k => abs) * visual(Scatter; color = :red, strokecolor = :black, strokewidth = 1) |> draw()


data(runs) * mapping(:Δθ => abs ∘ rad2deg, :k => abs) * visual(Scatter; color = :red, strokecolor = :black, strokewidth = 1) |> draw()


# df2 = only(eachrow(@rsubset df1 :run_id == 165))
# lines(df2.tθ, df2.θ)
# lines!(df2.tθ, df2.θs)
#
# df2.ks
#
# df3 = only(eachrow(@rsubset df1 :run_id == 93))
# df3.ks

data(runs) * mapping(:dance_induced, :k => abs, row = :at_run => nonnumeric) * visual(RainClouds, violin_limits = (0, Inf)) |> draw()

data(runs) * mapping(:dance_induced, :Δpoi => abs, row = :at_run => nonnumeric) * visual(RainClouds, violin_limits = (0, Inf)) |> draw()


# i = 7
# df = runs[i:i, :]
# (pregrouped(df.xy => first, df.xy => last)  * visual(Lines) + pregrouped(df.xys => first, df.xys => last)  * visual(Lines; color = :red)) * pregrouped(layout = string.(df.run_id, " ", df.condition)) |> draw(figure = (;size = (1000, 1000)), axis = (;aspect = DataAspect()))
#


# @chain runs begin 
#     @subset(:dance_induced)
#     pregrouped(_.xy => first, _.xy => last)  * visual(Lines) * pregrouped(layout = _.run_id => nonnumeric) |> draw(figure = (;size = (1000, 1000)), axis = (;aspect = DataAspect())) |> save("after.png")
# end


# x figure: 10 raw tracks: original direction
# x figure: 10 raw tracks: pointing to same direction
# x switch time to path length
# x figure: raw data for the turning angles
# - run glm (Gamma) with weights from the r square -> figure
# - maybe test variance differences on Δθ

sdjkfghasdklhfsdkh



# runs.condition .= categorical(runs.condition; levels = ["remain", "dark", "dark induced", "shift", "shift induced"], ordered = true)






df = subset(runs, :light => ByRow(==("shift")))

my_renamer = uppercasefirst ∘ string


# df1 = subset(df, :run_id => ByRow(==(8)))

# fig = pregrouped(df1.xys => first, df1.xys => last) * visual(Lines) |> draw(; axis = (; aspect = DataAspect()))
# display(fig)

# fig = (pregrouped(df1.xys => first, df1.xys => last) + pregrouped(df1.xy => first, df1.xy => last)) * visual(Lines) |> draw(; axis = (; aspect = DataAspect()))

# pregrouped(df.tθ, df.θ => rad2deg, row = df.dance_induced, col = df.at_run => nonnumeric) * visual(Lines) |> draw()

df1 = subset(df, :dance_induced => ByRow(x -> x))
(pregrouped(df1.tθ, df1.θ => rad2deg)  * visual(Lines) + pregrouped(df1.tθ, df1.θs => rad2deg)  * visual(Lines; color = :red)) * pregrouped(layout = string.(df1.run_id, " ", df1.condition)) |> draw()


(pregrouped(df.xy => first, df.xy => last)  * visual(Lines) + pregrouped(df.xys => first, df.xys => last)  * visual(Lines; color = :red)) * pregrouped(layout = string.(df.run_id, " ", df.condition)) |> draw()



sdjkhfgksdljfhsfhj

# ij = [(0,0), (1,0), (2,0), (3,1), (3,2), (1,2), (1,1), (2,0), (3,0), (4,0)]
# remove_cycles!(ij)
#
# ids = levelsmap(ij)
# vs = [ids[k] for k in ij]
# g = DiGraph(length(vs))
# for (v1, v2) in zip(vs[1:end-1], vs[2:end])
#     add_edge!(g, v1, v2)
# end
# c = only(simplecycles(g))
#
# xy = Point2f.(ij)
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# scatterlines!(ax, xy)
# scatter!(ax, xy[c], color = :red)
# deleteat!(xy, c)
# lines!(ax, xy, color = :green)
#
# deleteat!(ij, c)
#



df1 = subset(df, :run_id => ByRow(==(8)))
transform!(df1, :tij_file => ByRow(get_ij) => [:t, :ij])

# fig = pregrouped(df1.ij => first, df1.ij => last) * visual(ScatterLines) |> draw(; axis = (; aspect = DataAspect()))



ij = copy(df1.ij[1])
t = copy(df1.t[1])
tokill = findall(==(0) ∘ norm, diff(SV.(ij))) .+ 1
deleteat!(ij, tokill)
deleteat!(t, tokill)
t, ij = clean_coords(t, ij)
lines(SV.(ij), axis = (; aspect = DataAspect()))
# lines!(SV.(ij))




using Interpolations

itp = interpolate(a, BSpline(Constant()))

tokill = findall(==(0) ∘ norm, diff(SV.(ij))) .+ 1
deleteat!(ij, tokill)
lines(SV.(ij), axis = (; aspect = DataAspect()))
while remove_cycles!(ij)
end
lines!(SV.(ij))



function remove_cycles!(ij)
    ids = StatsBase.levelsmap(ij)
    vs = [ids[k] for k in ij]
    n = length(vs)
    g = DiGraph(n)
    for (v1, v2) in zip(vs[1:end-1], vs[2:end])
        add_edge!(g, v1, v2)
    end
    c = simplecycles(g)
    tokill = sort(unique(reduce(vcat, c)))
    deleteat!(ij, tokill)
    return ij
end

remove_cycles!(vs)

lines(SV.(vs), axis = (; aspect = DataAspect()))
lines!(SV.(df1.ij[1]))

dsjhflskdhfkdshflskhf


using Graphs

g = path_graph(60)

g = path_digraph(60)
graphplot(g)

using GLMakie, GraphMakie
using GraphMakie.NetworkLayout
g = smallgraph(:dodecahedral)
graphplot(g; layout=Stress(; dim=3))



ujfhsdfhjsdkh

######################## Plots 


# my_renamer = renamer("remain" => "Remain",
#                      "dark" => "Dark",
#                      "dark 10" => "Dark 10",
#                      "dark induced" => "Dark induced")

transform!(groupby(df, :condition), eachindex => :n)
combine(groupby(df, :condition), nrow).nrow



######################## Figure 1



######################## Figure 2

df1 = transform(df, [:xys, :tform] => ByRow((xy, f) -> f.(xy)) => :xyp)

fig = pregrouped(df1.xyp => first => "X (cm)", df1.xyp => last => "Y (cm)", col = df1.condition => my_renamer) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))

save(joinpath(output, "figure5.png"), fig)

################################################### Figure 3

df1 = transform(df, [:t, :poi_index, :tform] => ByRow((t, i, f) -> f.(t[i:end])) => :xyp)
nr = 3
l = floor(Int, minimum(norm ∘ last, df1.xyp))
r = range(1e-3, l, nr)
transform!(df1, :xyp => ByRow(xyp -> get_exit_angle.(Ref(xyp), r)) => :θs)
select!(df1, Cols(:condition, :θs))
df1.r .= Ref(r)

rl = range(extrema(r)..., 100)
conditions = levels(df1.condition)
newdf = DataFrame(r = repeat(rl, outer = length(conditions)), condition = repeat(conditions, inner = 100), mean_resultant_vector = zeros(100length(conditions)))
fm = @formula(mean_resultant_vector ~ r*condition)
function _bootstrap(df)
    n = nrow(df)
    df1 = flatten(df[sample(1:n, n), :], [:θs, :r])
    df2 = combine(groupby(df1, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
    try
        m = BetaRegression.fit(BetaRegressionModel, fm, df2)
        predict(m, newdf)
    catch ex
        missing
    end
end
function bootstrap(df)
    n = 10_000
    ys = Vector{Vector{Float64}}(undef, n)
    i = 0
    while i < n
        y = _bootstrap(df)
        if ismissing(y)
            continue
        else
            i += 1
            ys[i] = y
        end
    end
    return stack(ys)
end
c = bootstrap(df1)
y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdf.lower .= getindex.(y, 1)
newdf.mean_resultant_vector .= getindex.(y, 2)
newdf.upper .= getindex.(y, 3)

fig = data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :condition => my_renamer => "Light") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant length", xlabel = "Radius (cm)", limits = ((0, l), (0, 1))))

save(joinpath(output, "figure6.png"), fig)
