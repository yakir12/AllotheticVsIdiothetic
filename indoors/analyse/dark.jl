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

include("minimal_functions.jl")

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
    @rtransform! :pixels = get_tij(:tij_file)
    @rtransform! :pixels = remove_stops(:pixels)
    # @aside @assert !any(has_stops, _.pixels) "some stops remain?!" # convert to a test
    @transform! :xy = trectify(:rectify, :pixels)
    @aside @chain _ begin 
        @subset(:dance_by .≠ "no"; view = true)
        @rtransform! :jump = glue_intervention!(:xy, :intervention)
    end
    @aside @chain _ begin 
        @subset(:light .== "remain"; view = true)
        @rtransform! :intervention = impute_poi_time(:xy)
    end
    @rtransform! :spontaneous_end = passmissing(tosecond)(:spontaneous_end)
    @rtransform! :poi = coalesce(:spontaneous_end, :intervention)
    disallowmissing!(:poi)
    @select! Not(:intervention)
    @transform! :xy = remove_loops.(:xy)
    @transform! :xy = sparseify.(:xy)
    @transform! :smooth = smooth.(:xy)
    @transform! :centered2start = center2start.(:smooth)
    @transform! :cropped = cropto.(:centered2start, l)
    @transform! :rotated2poi = rotate2poi.(:cropped, :poi)
    @transform! :centered2poi_and_cropped = center2poi_and_crop.(:rotated2poi, :poi)
    @transform! :dance_spontaneous = .!ismissing.(:spontaneous_end)
    @transform! :condition = string.(:light, " ", :dance_by, " ", :at_run)
    @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
end


askjdhfksahflksahls

#     # @rtransform! :poi_index = something(findfirst(≥(:poi), :t), length(:t))
#     @rtransform! :xysr = :rot.(:xys)
#     @rtransform! :xysrc = :xysr[:poi_index:end] .- Ref(:xysr[:poi_index])
#     # @rtransform! $AsTable = get_turn_profile(:t, :spl, :poi)
#     # @rtransform! $AsTable = fit_logistic(:lθ, :θ)
#     # transform!(:ks => ByRow(identity) => [:L, :k, :x₀, :y₀])
#     # @rtransform! :θs = logistic.(:lθ, :L, :k, :x₀, :y₀)
#     @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
#     @aside pregrouped(_.xys => first, _.xys => last)  * visual(Lines) * pregrouped(layout = _.run_id => nonnumeric) |> draw(figure = (;size = (1000, 1000)), axis = (;aspect = DataAspect())) |> save("summary.png")
# end

############ plot the tyracks to check validity
CairoMakie.activate!()
path = "tmp"
mkpath(path)
rm.(readdir(path, join = true))
Threads.@threads for row in eachrow(runs)
    run_id, xy, xys, poi_index, dance_by = (row.run_id, row.xy, row.xys, row.poi_index, row.dance_by)
    # run_id, xy, xys, poi_index, lθ, θ, θs, lpoi, lpoi_index, dance_by = (row.run_id, row.xy, row.xys, row.poi_index, row.lθ, row.θ, row.θs, row.lpoi, row.lpoi_index, row.dance_by)
    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect(), title = dance_by)
    lines!(ax, xy[1:poi_index], color = :black)
    lines!(ax, xy[poi_index:end], color = :gray)
    lines!(ax, xys[1:poi_index], color = :red)
    lines!(ax, xys[poi_index:end], color = :orange)
    # ax = Axis(fig[1,2])
    # lines!(ax, lθ[1:lpoi_index], rad2deg.(θ[1:lpoi_index]), color = :black)
    # lines!(ax, lθ[lpoi_index:end], rad2deg.(θ[lpoi_index:end]), color = :gray)
    # lines!(ax, lθ[1:lpoi_index], rad2deg.(θs[1:lpoi_index]), color = :red)
    # lines!(ax, lθ[lpoi_index:end], rad2deg.(θs[lpoi_index:end]), color = :orange)
    save(joinpath(path, "$(run_id).png"), fig)
end
GLMakie.activate!()

################################################### 10 random tracks

n = 10
df = @chain runs begin
    @subset norm.(last.(:cropped)) .≈ l
    @groupby :dance_by
    @transform :n = 1:length(:dance_by)
    @subset :n .≤ n
    @transform :cropped = parent.(:cropped)
    @transform :rotated2poi = parent.(:rotated2poi)
end
@assert all(==(n), combine(groupby(df, :dance_by), nrow).nrow)


# df = DataFrame(tracks = [DimVector(rand(SVector{2, Float64}, 10), Ti(1:10)) for _ in 1:3], name = string.('a':'c')) 
# pregrouped(df.tracks => first, df.tracks => last, color = df.name) * visual(Lines) |> draw()
#
# df = flatten(df, :tracks)
# transform!(df, :tracks => [:x, :y])
# data(df) * mapping(:x, :y, color = :name) * visual(Lines) |> draw()
#
# df = DataFrame(tracks = [rand(SVector{2, Float64}, 10) for _ in 1:3], name = string.('a':'c')) 
# pregrouped(df.tracks => first, df.tracks => last, color = df.name) * visual(Lines) |> draw()

fig = pregrouped(df.cropped => first => "X (cm)", df.cropped => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, l)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure1.png"), fig)

fig = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
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

df = @transform runs :rotated2poi = parent.(:rotated2poi)
fig = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)", col = df.condition, color = df.y2025) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))
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
# @rtransform! df :condition = :light == "remain" ? "remain" : :dance_by
nr = 3
l = floor(Int, minimum(norm ∘ last, df.xysrc))
rl = range(1e-3, l, nr)
transform!(df, :xysrc => ByRow(xysrc -> get_exit_angle.(Ref(xysrc), rl)) => :θs)
select!(df, Cols(:light, :dance_by, :θs, :xysrc))
df.light = categorical(df.light)
levels!(df.light, ["remain", "dark"])
df.dance_by = categorical(df.dance_by)
levels!(df.dance_by, ["no", "hold", "disrupt"])

R = 60
fig = pregrouped(df.xysrc => first => "X (cm)", df.xysrc => last => "Y (cm)", col = df.light => renamer("remain" => "Lights on", "dark" => "Lights off"), row = df.dance_by => renamer("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball")) * visual(Lines) |> draw(; axis=(width = 300, height = 300, limits = ((-R, R), (-R, R))))
# , width = 300, height = 300

save(joinpath(output, "figure2a.png"), fig)

light_rename = Dict("remain" => "Lights on", "dark" => "Lights off")
dance_by_rename = Dict("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball")
fig = Figure()
for ((light, dance_by), g) in pairs(groupby(df, [:light, :dance_by]))
    i = findfirst(==(dance_by), levels(df.dance_by))
    j = findfirst(==(light), levels(df.light))
    ax = PolarAxis(fig[i, j], rlimits = (0, 44), width = 300, height = 300)
    for xy in g.xysrc
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


select!(df, Not(:xysrc))
# select!(df, Cols(:condition, :θs))
df.r .= Ref(rl)






# df.condition = categorical(df.condition)
# levels!(df.condition, ["remain", "no", "hold", "disrupt"])


# fig = Figure()
# for (i, (k, g)) in enumerate(pairs(groupby(df, :condition)))
#     ax = PolarAxis(fig[i,1])
#     # _df = flatten(g, [:θs, :r])
#     # scatter!(ax, _df.θs, _df.r)
#     for row in eachrow(g)
#         scatter!(ax, row.θs, row.r)
#     end
# end

newdf = combine(groupby(df, [:light, :dance_by]), :r => first ∘ first => :r)
nr2 = 100
rl2 = range(extrema(rl)..., nr2)
newdf.r .= Ref(rl2)
newdf = flatten(newdf, :r)
newdf.mean_resultant_vector .= 0.0

fm = @formula(mean_resultant_vector ~ light + r + dance_by)

# conditions = levels(df.condition)
# nr2 = 100
# rl2 = range(extrema(rl)..., nr2)
# newdf = DataFrame(r = repeat(rl2, outer = length(conditions)), condition = repeat(conditions, inner = nr2), mean_resultant_vector = zeros(nr2*length(conditions)))
# fm = @formula(mean_resultant_vector ~ r*condition)

pc, c = bootstrap(df)
select!(pc, Not(Symbol("(Precision)")))

y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdf.lower .= getindex.(y, 1)
newdf.mean_resultant_vector .= getindex.(y, 2)
newdf.upper .= getindex.(y, 3)
pvalues = combine(pc, All() .=> stats, renamecols = false)
pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
df2 = stack(pvalues, Not(:what), variable_name = :source)
pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])


df2 = combine(groupby(df, [:light, :dance_by]), nrow)
@transform! df2 :r̄ = critical_r.(:nrow)


# n = nrow(df)
# df2 = flatten(df, [:θs, :r])
# df3 = combine(groupby(df2, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
# df3.condition = categorical(df3.condition)
# levels!(df3.condition, ["remain", "no", "hold", "disrupt"])
# m = BetaRegression.fit(BetaRegressionModel, fm, df3)

# + data(df2) * mapping(:r̄, col = :light, color = :dance_by => renamer("no" => "No disruption", "hold" => "Hold", "disrupt" => "Taken off ball") => "Light") * visual(HLines) 
#, yticks = ([0; sort(df2.r̄); 0.5; 1], ["0", "", "", "", "", "0.5", "1"]

fig = data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, group = [:light, :dance_by], col = :light => renamer("remain" => "Lights on", "dark" => "Lights off"), color = :dance_by => renamer("no" => "Nothing", "hold" => "Holding ball", "disrupt" => "Removing from ball") => "Dance induced by") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant vector length", xlabel = "Radius (cm)", width = 300, height = 300, limits = ((0, l), (0, 1))))



save(joinpath(output, "figure3.png"), fig)

show(pvalues, show_row_number=false, eltypes=false)

