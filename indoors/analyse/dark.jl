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

GLMakie.activate!()

include("minimal_functions.jl")

output = "figures"
# if isdir(output)
#     rm.(readdir(output; join = true))
# else
#     mkdir("figures")
# end
#
const results_dir = "../track_calibrate/tracks and calibrations"

function combine_factors(light, induced, run)
    induced = induced ? " induced" : ""
    run = run > 1 ? " $run" : ""
    string(light, induced, run)
end

function convert_dance_by_to_binary(dance_by)
    dance_by ∈ ("disrupt", "hold") && return true
    dance_by == "no" && return false
    error("third dance_by option: $dance_by")
end

runs = @chain joinpath(results_dir, "runs.csv") begin
    CSV.read(DataFrame)
    @subset :light .≠ "shift"
    @transform :dance_induced = convert_dance_by_to_binary.(:dance_by)
    @select Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size)
end

calibs = @chain joinpath(results_dir, "calibs.csv") begin
    CSV.read(DataFrame)
    @transform :rectify = get_calibration.(:calibration_id)
    @select Cols(:calibration_id, :rectify)
end
leftjoin!(runs, calibs, on = :calibration_id)
@chain runs begin
    @select! Not(:calibration_id)
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
    @rtransform! $AsTable = get_turn_profile(:t, :spl, :poi)
    @rtransform! $AsTable = fit_logistic(:lθ, :θ)
    transform!(:ks => ByRow(identity) => [:L, :k, :x₀, :y₀])
    @rtransform! :θs = logistic.(:lθ, :L, :k, :x₀, :y₀)
    @rtransform! :y2025 = Year(:start_datetime) == 2025 ? "2025" : "earlier"
end

############ plot the tyracks to check validity
CairoMakie.activate!()
path = "tmp"
mkpath(path)
rm.(readdir(path, join = true))
Threads.@threads for row in eachrow(runs)
    run_id, xy, xys, poi_index, lθ, θ, θs, lpoi, lpoi_index, dance_induced = (row.run_id, row.xy, row.xys, row.poi_index, row.lθ, row.θ, row.θs, row.lpoi, row.lpoi_index, row.dance_induced)
    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect(), title = string(dance_induced, " ", round(Int, rad2deg(row.L)), "°"))
    lines!(ax, xy[1:poi_index], color = :black)
    lines!(ax, xy[poi_index:end], color = :gray)
    lines!(ax, xys[1:poi_index], color = :red)
    lines!(ax, xys[poi_index:end], color = :orange)
    ax = Axis(fig[1,2])
    lines!(ax, lθ[1:lpoi_index], rad2deg.(θ[1:lpoi_index]), color = :black)
    lines!(ax, lθ[lpoi_index:end], rad2deg.(θ[lpoi_index:end]), color = :gray)
    lines!(ax, lθ[1:lpoi_index], rad2deg.(θs[1:lpoi_index]), color = :red)
    lines!(ax, lθ[lpoi_index:end], rad2deg.(θs[lpoi_index:end]), color = :orange)
    save(joinpath(path, "$(run_id).png"), fig)
end
GLMakie.activate!()

################################################### 10 random tracks

l = 50
n = 10
df = @chain runs begin
    @subset norm.(last.(:xys)) .> l
    @rtransform :xys = cropto(:t, :spl, :trans, l)
    @rtransform :xysr = :rot.(:xys)
    @groupby :dance_induced
    @transform :n = 1:length(:dance_induced)
    @subset :n .≤ n
end
@assert all(==(n), combine(groupby(df, :dance_induced), nrow).nrow)

fig = pregrouped(df.xys => first => "X (cm)", df.xys => last => "Y (cm)", col = df.dance_induced => renamer(true => "Induced", false => "Not"), color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
save(joinpath(output, "figure1.png"), fig)

fig = pregrouped(df.xysr => first => "X (cm)", df.xysr => last => "Y (cm)", col = df.dance_induced => renamer(true => "Induced", false => "Not"), color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
save(joinpath(output, "figure1a.png"), fig)


################################################### raw data for turning angles


######################## Figure 2

df1 = transform(df, [:t, :poi_index, :tform] => ByRow((t, i, f) -> f.(t[i:end])) => :xyp)

fig = pregrouped(df1.xyp => first => "X (cm)", df1.xyp => last => "Y (cm)", col = df1.condition => my_renamer) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))

save(joinpath(output, "figure2.png"), fig)

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

save(joinpath(output, "figure3.png"), fig)
