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
include("shift_functions.jl")

const l = 50

output = "shift"
if isdir(output)
    rm.(readdir(output; join = true))
else
    mkdir(output)
end

const results_dir = "../track_calibrate/tracks and calibrations"

runs = @chain joinpath(results_dir, "runs.csv") begin
    CSV.read(DataFrame)
    @select Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size)
    @subset :light .== "shift"
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
    @transform! :dance_spontaneous = .!ismissing.(:spontaneous_end)
    @transform! :condition = string.(:light, " ", :dance_by, " ", :at_run)
    @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
    @rename! :intervention = :poi
    @transform! :pixels = get_tij.(:tij_file)
    @transform! :pixels = remove_stops.(:pixels)
    # @aside @assert !any(has_stops, _.pixels) "some stops remain?!" # convert to a test
    @rtransform! :xy = :rectify.(:pixels)
    @aside @chain _ begin 
        @subset(:dance_by .≠ "disrupt"; view = true)
        @transform! :jump = glue_intervention!.(:xy, :intervention)
    end
    # @aside @chain _ begin 
    #     @subset(:light .== "remain"; view = true)
    #     @rtransform! :poi = impute_poi_time(:t, :xy)
    # end
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
    # transform!([:pixels, :xy, :smooth, :centered2start, :cropped, :rotated2poi, :centered2poi_and_cropped] .=> ByRow(parent), renamecols = false)

    # @rtransform! $AsTable = fit_logistic(:lθ, :θ)
    # transform!(:ks => ByRow(identity) => [:L, :k, :x₀, :y₀])
    # @rtransform! :θs = logistic.(:lθ, :L, :k, :x₀, :y₀)
end

sdjfhksdjhflasfhj

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

# (pregrouped(map(x -> fill(x, 2), runs.lpoi), fill([-360, 360], nrow(runs)))  * visual(Lines; color = :green) + pregrouped(runs.lθ, runs.θ => rad2deg)  * visual(Lines) + pregrouped(runs.lθ, runs.θs => rad2deg)  * visual(Lines; color = :red)) * pregrouped(layout = runs.run_id => nonnumeric) |> draw(facet = (; linkxaxes = :none, linkyaxes = :all)) |> display



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
save(joinpath(output, "figure4.png"), fig)

fig = pregrouped(df.xysr => first => "X (cm)", df.xysr => last => "Y (cm)", col = df.dance_induced => renamer(true => "Induced", false => "Not"), color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
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
    @subset :logistic_rsquare .> 0.95
    @transform :Δθ = wrap2pi.(:L)
    transform(:L => ByRow(l -> sincos(l .+ π/2)) => [:v, :u])
    @transform :Δl = (:x₀ .- :lpoi)
    @transform :absk = abs.(:k)
end

fig = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v, col = :dance_induced => renamer(false => "Not", true => "Induced"), row = :at_run => nonnumeric, color = :Δθ => rad2deg => "Total turn") * visual(Arrows) |> draw(scales(Color = (; colormap = :cyclic_wrwbw_40_90_c42_n256_s25, colorrange = (-180, 180))); axis = (; width = 400, height = 400), colorbar = (; ticks = -180:90:180, tickformat = "{:n}°"))

save(joinpath(output, "figure5.png"), fig)

# # PCA
@chain df begin
    @transform! :Δθnormalized = (π .- abs.(:Δθ)) ./ π 
    @transform! :kabs = abs.(:k)
    @transform! :Δlabs = abs.(:Δl)
end
#
# Xtr = Matrix(select(df, :Δθnormalized, :kabs, :Δlabs))
# M = MultivariateStats.fit(PCA, Xtr; maxoutdim=3)
# data((; at_run = df.at_run, dance_induced = df.dance_induced, pc1 = first(eachcol(M.proj)), pc2 = last(eachcol(M.proj)))) * mapping(:pc1, :pc2, col = :dance_induced, row = :at_run => nonnumeric) * visual(Scatter) |> draw()



# correlation tests
cor(df.Δlabs, df.kabs)
cor(df.Δlabs, df.Δθnormalized)
cor(df.kabs, df.Δθnormalized)

# m = glm(@formula(kabs ~ dance_induced*at_run), df, Gamma())
# m = glm(@formula(kabs ~ dance_induced + at_run), df, Gamma())

m = glm(@formula(kabs ~ dance_induced), df, Gamma())

# predictions = DataFrame(dance_induced = [true, false])
# predictions = hcat(predictions, rename(predict(m, predictions; interval = :confidence), :prediction => :k, :lower => :lower_k, :upper => :upper_k))


# m = glm(@formula(Δlabs ~ dance_induced*at_run), df, Gamma())
# m = glm(@formula(Δlabs ~ dance_induced + at_run), df, Gamma())
m = glm(@formula(Δlabs ~ dance_induced), df, Gamma())

# predictions = hcat(predictions, rename(predict(m, predictions; interval = :confidence), :prediction => :l, :lower => :lower_l, :upper => :upper_l))

# m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_induced*at_run), df)
# m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_induced + at_run), df)
m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_induced), df)

# predict(m, predictions)


@chain df begin
    @rtransform! :lθshifted = :lθ .- :lpoi # .- :x₀
    @rtransform! :θnormalized = :θ .- :θ[:lpoi_index]#:k < 0 ? :θ .+ :y₀ .- π : :θ .+ :y₀
end
fig = pregrouped(df.lθshifted => "Distance from POI (path length cm)", df.θnormalized => rad2deg) * visual(Lines) * mapping(col = df.dance_induced => renamer(false => "Not", true => "Induced"), row = df.at_run => nonnumeric, color = df.y2025) |> draw()#; axis = (; width = 400, height = 400, ytickformat = "{:n}°", yticks = -180:90:270, limits = ((-5, 5), nothing)))

save(joinpath(output, "figure6.png"), fig)

