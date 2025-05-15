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

function combine_factors(light, dance_by, run)
    "$light $dance_by $run"
    # dance_by = dance_by == " no" ? "" : dance_by
    # run = run > 1 ? " $run" : ""
    # string(light, dance_by, run)
end

# function convert_dance_by_to_binary(dance_by)
#     dance_by ∈ ("disrupt", "hold") && return true
#     dance_by == "no" && return false
#     error("third dance_by option: $dance_by")
# end

runs = @chain joinpath(results_dir, "runs.csv") begin
    CSV.read(DataFrame)
    @select Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size)
    @subset :light .≠ "shift"
    # @transform :dance_induced = convert_dance_by_to_binary.(:dance_by)
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
    @transform! :_distrupt = :dance_by .≠ "no"
    @aside @chain _ begin 
        @subset(:_distrupt; view = true)
        # @subset(:dance_by => x -> x .≠ "no"; view = true)
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
    @rtransform! :condition = combine_factors(:light, :dance_by, :at_run)
    @rtransform! :rot = get_rotation(:xys[:poi_index])
    @rtransform! :xysr = :rot.(:xys)
    @rtransform! :xysrc = :xysr[:poi_index:end] .- Ref(:xysr[:poi_index])
    @rtransform! $AsTable = get_turn_profile(:t, :spl, :poi)
    @rtransform! $AsTable = fit_logistic(:lθ, :θ)
    transform!(:ks => ByRow(identity) => [:L, :k, :x₀, :y₀])
    @rtransform! :θs = logistic.(:lθ, :L, :k, :x₀, :y₀)
    @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
end

############ plot the tyracks to check validity
CairoMakie.activate!()
path = "tmp"
mkpath(path)
rm.(readdir(path, join = true))
Threads.@threads for row in eachrow(runs)
    run_id, xy, xys, poi_index, lθ, θ, θs, lpoi, lpoi_index, dance_by = (row.run_id, row.xy, row.xys, row.poi_index, row.lθ, row.θ, row.θs, row.lpoi, row.lpoi_index, row.dance_by)
    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect(), title = string(dance_by, " ", round(Int, rad2deg(row.L)), "°"))
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
    @groupby :dance_by
    @transform :n = 1:length(:dance_by)
    @subset :n .≤ n
end
@assert all(==(n), combine(groupby(df, :dance_by), nrow).nrow)

fig = pregrouped(df.xys => first => "X (cm)", df.xys => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
save(joinpath(output, "figure1.png"), fig)

fig = pregrouped(df.xysr => first => "X (cm)", df.xysr => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
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

fig = pregrouped(runs.xysr => first => "X (cm)", runs.xysr => last => "Y (cm)", col = runs.condition, color = runs.y2025) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))

save(joinpath(output, "figure2.png"), fig)

fig = pregrouped(runs.xysrc => first => "X (cm)", runs.xysrc => last => "Y (cm)", col = runs.light => sorter(["remain", "dark"]), row = runs.dance_by => sorter(["no", "hold", "disrupt"]),  color = runs.y2025 => "Recorded at", linestyle = runs.at_run => nonnumeric => "At run") * visual(Lines) |> draw(; figure = (; size = (1200, 1000)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))));
save(joinpath(output, "figure2a.png"), fig)

################################################### Figure 3


df = @subset runs :at_run .== 1
@rtransform! df :condition = :light == "remain" ? "remain" : :dance_by
nr = 3
l = floor(Int, minimum(norm ∘ last, df.xysrc))
rl = range(1e-3, l, nr)
transform!(df, :xysrc => ByRow(xysrc -> get_exit_angle.(Ref(xysrc), rl)) => :θs)

select!(df, Cols(:condition, :θs))
df.r .= Ref(rl)

df.condition = categorical(df.condition)
levels!(df.condition, ["remain", "no", "hold", "disrupt"])


fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df, :condition)))
    ax = PolarAxis(fig[i,1])
    # _df = flatten(g, [:θs, :r])
    # scatter!(ax, _df.θs, _df.r)
    for row in eachrow(g)
        scatter!(ax, row.θs, row.r)
    end
end


# df2 = flatten(df, [:θs, :r])
# data(df2) * mapping(:r, :θs, row = :condition => renamer("remain" => "Light on", "no" => "No disruption", "hold" => "Hold", "disrupt" => "Taken off ball")) * visual(Violin) |> draw()

conditions = levels(df.condition)
nr2 = 100
rl2 = range(extrema(rl)..., nr2)
newdf = DataFrame(r = repeat(rl2, outer = length(conditions)), condition = repeat(conditions, inner = nr2), mean_resultant_vector = zeros(nr2*length(conditions)))
fm = @formula(mean_resultant_vector ~ r*condition)

function _bootstrap(df)
    n = nrow(df)
    df2 = flatten(df[sample(1:n, n), :], [:θs, :r])
    df3 = combine(groupby(df2, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
    df3.condition = categorical(df3.condition)
    levels!(df3.condition, ["remain", "no", "hold", "disrupt"])
    try
        m = BetaRegression.fit(BetaRegressionModel, fm, df3)
        tbl = coeftable(m)
        row = (; Pair.(Symbol.(tbl.rownms), tbl.cols[tbl.pvalcol])...)
        # row = tbl.cols[tbl.pvalcol]
        (row,  predict(m, newdf))
    catch ex
        missing
    end
end
function bootstrap(df)
    n = 10_000
    ys = Vector{Vector{Float64}}(undef, n)
    ys = Matrix{Float64}(undef, nrow(newdf), n)
    rows = DataFrame(var"(Intercept)" = Float64[], r = Float64[], var"condition: hold" = Float64[], var"condition: disrupt" = Float64[], var"condition: no" = Float64[], var"r & condition: hold" = Float64[], var"r & condition: disrupt" = Float64[], var"r & condition: no" = Float64[], var"(Precision)" = Float64[])
    i = 0
    while i < n
        rowy = _bootstrap(df)
        if ismissing(rowy)
            continue
        else
            i += 1
            row, ys[:, i] = rowy
            push!(rows, row)
        end
    end
    return rows, stack(ys)
end
pc, c = bootstrap(df)
select!(pc, Not(Symbol("(Precision)")))

y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdf.lower .= getindex.(y, 1)
newdf.mean_resultant_vector .= getindex.(y, 2)
newdf.upper .= getindex.(y, 3)
function stats(x)
    q1, med, q2 = quantile(x, [0.025, 0.5, 0.975])
    mod = mode(x)
    μ = mean(x)
    [q1, med, q2, mod, μ]
end
pvalues = combine(pc, All() .=> stats, renamecols = false)
pvalues.what = ["q1", "median", "q2", "mode", "mean"]
df2 = stack(pvalues, variable_name = :source)
pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["q1", "median", "q2", "mode", "mean"])

critical R values!!!!

# n = nrow(df)
# df2 = flatten(df, [:θs, :r])
# df3 = combine(groupby(df2, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
# df3.condition = categorical(df3.condition)
# levels!(df3.condition, ["remain", "no", "hold", "disrupt"])
# m = BetaRegression.fit(BetaRegressionModel, fm, df3)



fig = data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :condition => renamer("remain" => "Light on", "no" => "No disruption", "hold" => "Hold", "disrupt" => "Taken off ball") => "Light") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant length", xlabel = "Radius (cm)", limits = ((0, l), (0, 1))));

save(joinpath(output, "figure3.png"), fig)


CSV.write("$nr.csv", pvalues)
