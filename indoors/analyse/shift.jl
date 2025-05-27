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
    @transform! :θ = get_turn_profile.(:smooth, :poi)
    @transform! $AsTable = fit_logistic.(:θ)
    transform!(:ks => ByRow(identity) => [:L, :k, :x₀, :y₀])
    @rtransform! :θs = logistic.(:θ.l, :L, :k, :x₀, :y₀)
    transform!([:pixels, :xy, :smooth, :centered2start, :cropped, :rotated2poi, :centered2poi_and_cropped] .=> ByRow(parent), renamecols = false)
end


############ plot the tracks to check validity

fig = (pregrouped(runs.smooth => first => "X (cm)", runs.smooth => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines; color = :red) + pregrouped(runs.xy => first => "X (cm)", runs.xy => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines)) |> draw(; axis = (; width = 400, height = 400, limits = ((-l, l), (-l, l))));
save(joinpath(output, "overview_shift.png"), fig)

################################################### 10 random tracks

n = 10
df = @chain runs begin
    @subset norm.(last.(:cropped)) .≈ l
    @groupby :dance_by
    @transform :n = 1:length(:dance_by)
    @subset :n .≤ n
end
@assert all(==(n), combine(groupby(df, :dance_by), nrow).nrow)

fig = pregrouped(df.cropped => first => "X (cm)", df.cropped => last => "Y (cm)", col = df.dance_by) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, l)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure1.png"), fig)

fig = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)", col = df.dance_by) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure1a.png"), fig)

################################################### raw data for turning angles


df = @chain runs begin
    @subset :logistic_rsquare .> 0.95
    @transform :Δθ = wrap2pi.(:L)
    transform(:L => ByRow(l -> sincos(l .+ π/2)) => [:v, :u])
    @rtransform :Δl = (:x₀ .- :θ.l[Ti = Near(:poi)])
    @transform :absk = abs.(:k)
end

fig = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v, col = :dance_by, row = :at_run => nonnumeric, color = :Δθ => rad2deg => "Total turn") * visual(Arrows) |> draw(scales(Color = (; colormap = :cyclic_wrwbw_40_90_c42_n256_s25, colorrange = (-180, 180))); axis = (; width = 400, height = 400), colorbar = (; ticks = -180:90:180, tickformat = "{:n}°"))

save(joinpath(output, "figure5.png"), fig)

# # PCA
@chain df begin
    @transform! :Δθnormalized = (π .- abs.(:Δθ)) ./ π 
    @transform! :kabs = abs.(:k)
    @transform! :Δlabs = abs.(:Δl)
end

# Xtr = Matrix(select(df, :Δθnormalized, :kabs, :Δlabs))
# M = MultivariateStats.fit(PCA, Xtr; maxoutdim=3)
# data((; at_run = df.at_run, dance_by = df.dance_by, pc1 = first(eachcol(M.proj)), pc2 = last(eachcol(M.proj)))) * mapping(:pc1, :pc2, col = :dance_by, row = :at_run => nonnumeric) * visual(Scatter) |> draw()



# correlation tests
cor(df.Δlabs, df.kabs)
cor(df.Δlabs, df.Δθnormalized)
cor(df.kabs, df.Δθnormalized)

m = glm(@formula(kabs ~ dance_by*at_run), df, Gamma())
m = glm(@formula(kabs ~ dance_by + at_run), df, Gamma())

m = glm(@formula(kabs ~ dance_by), df, Gamma())

# predictions = DataFrame(dance_by = [true, false])
# predictions = hcat(predictions, rename(predict(m, predictions; interval = :confidence), :prediction => :k, :lower => :lower_k, :upper => :upper_k))


# m = glm(@formula(Δlabs ~ dance_by*at_run), df, Gamma())
# m = glm(@formula(Δlabs ~ dance_by + at_run), df, Gamma())
m = glm(@formula(Δlabs ~ dance_by), df, Gamma())

# predictions = hcat(predictions, rename(predict(m, predictions; interval = :confidence), :prediction => :l, :lower => :lower_l, :upper => :upper_l))

# m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_by*at_run), df)
# m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_by + at_run), df)
m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_by), df)

# predict(m, predictions)


@chain df begin
    @rtransform! :lθshifted = parent(:θ.l .- :θ.l[Ti = Near(:poi)]) # .- :x₀
    @rtransform! :θnormalized = parent(:θ.θ .- :θ.θ[Ti = Near(:poi)])#:k < 0 ? :θ .+ :y₀ .- π : :θ .+ :y₀
end

fig = pregrouped(df.lθshifted => "Distance from POI (path length cm)", df.θnormalized => rad2deg) * visual(Lines) * mapping(col = df.dance_by, row = df.at_run => nonnumeric, color = df.y2025) |> draw()#; axis = (; width = 400, height = 400, ytickformat = "{:n}°", yticks = -180:90:270, limits = ((-5, 5), nothing)))

save(joinpath(output, "figure6.png"), fig)

