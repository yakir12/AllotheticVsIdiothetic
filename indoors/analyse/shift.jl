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
    # @rtransform! :tform = get_rotation(:xys[:poi_index])
    @rtransform! :l = get_pathlength(:xys)
    @rtransform! $AsTable = get_turn_profile(:t, :spl, :poi_index)
    @rtransform! :ks = fit_logistic(:tθ, :θ)
    @rtransform! :θs = logistic.(:tθ, Ref([:ks; 0]))
    @transform! :k = first.(:ks)
end

function loops(xy)
    inds, _ = self_intersections(Point2f.(xy))
    ranges = splat(UnitRange).(Iterators.partition(inds, 2))
    length.(ranges)
end

l = map(loops, df.xy)
l = reduce(vcat, l)

df = copy(runs)
fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect(), limits = ((-100, 100), (-100, 100)))
sl = SliderGrid(fig[3, 1], (range = 1:nrow(df), ))
i = only(sl.sliders).value
poi1 = lift(i) do i
    df.xy[i][1:df.poi_index[i]]
end
poi2 = lift(i) do i
    df.xys[i][1:df.poi_index[i]]
end
tra1 = lift(i) do i
    df.xy[i][df.poi_index[i]:end]
end
tra2 = lift(i) do i
    df.xys[i][df.poi_index[i]:end]
end
sl = SliderGrid(fig[4, 1], (range = 1:minimum(length, df.t), ))
j = only(sl.sliders).value
xy1 = lift(i, j) do i, j
    df.xy[i][j]
end
xy2 = lift(i, j) do i, j
    df.xys[i][j]
end
lines!(ax, poi1, color = :black)
lines!(ax, poi2, color = :red)
lines!(ax, tra1, color = :blue)
lines!(ax, tra2, color = :green)
scatter!(ax, xy1)
scatter!(ax, xy2)
ax1 = Axis(fig[2,1], limits = (nothing, (-90, 360)), yticks = -90:90:360)
turn1 = lift(i) do i
    Point2f.(df.tθ[i][1:df.poi_index[i]], rad2deg.(df.θ[i][1:df.poi_index[i]]))
end
turn2 = lift(i) do i
    Point2f.(df.tθ[i][df.poi_index[i]:end], rad2deg.(df.θ[i][df.poi_index[i]:end]))
end
txy = lift(i, j) do i, j
    Point2f(df.tθ[i][j], rad2deg(df.θ[i][j]))
end
lines!(ax1, turn1, color = :red)
lines!(ax1, turn2, color = :green)
scatter!(ax1, txy)
on(turn) do t
    limits!(ax1, -1, first(last(t)) + 1, -90, 360)
end



df1 = @rsubset(runs, !(:dance_induced))

(pregrouped(df1.tθ, df1.θ => rad2deg)  * visual(Lines) + pregrouped(df1.tθ, df1.θs => rad2deg)  * visual(Lines; color = :red)) * pregrouped(layout = string.(df1.run_id, " ", df1.condition)) |> draw(axis = (;limits = (nothing, (-10, 360))))

data(runs) * mapping(:dance_induced, :k, row = :at_run => nonnumeric) * visual(RainClouds, violin_limits = (0, Inf)) |> draw()


i = 7
df = runs[i:i, :]
(pregrouped(df.xy => first, df.xy => last)  * visual(Lines) + pregrouped(df.xys => first, df.xys => last)  * visual(Lines; color = :red)) * pregrouped(layout = string.(df.run_id, " ", df.condition)) |> draw(figure = (;size = (1000, 1000)), axis = (;aspect = DataAspect()))



# @chain runs begin 
#     @subset(:dance_induced)
#     pregrouped(_.xy => first, _.xy => last)  * visual(Lines) * pregrouped(layout = _.run_id => nonnumeric) |> draw(figure = (;size = (1000, 1000)), axis = (;aspect = DataAspect())) |> save("after.png")
# end



sdjkfghasdklhfsdkh

julia> findall(>(360), rad2deg.(last.(runs.ks)) ./ 2)
3-element Vector{Int64}:
  72
 112
 123


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

df1 = subset(df, :xys => ByRow(≥(50) ∘ norm ∘ last))
transform!(df1, :xys => ByRow(xy -> cropto(xy, 50)) => :xys)
transform!(groupby(df1, :condition), eachindex => :n)
subset!(df1, :n => ByRow(≤(10)))
@assert all(==(10), combine(groupby(df1, :condition), nrow).nrow)

fig = pregrouped(df1.xys => first, df1.xys => last, col = df1.condition => my_renamer) * visual(Lines) |> draw(; figure = (; size = (1200, 400)), axis=(aspect=DataAspect(), ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure4.png"), fig)


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
