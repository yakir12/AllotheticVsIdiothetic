using SimpTrack
using Dates, LinearAlgebra
using CSV, DataFrames, Chain, DataFramesMeta, IterTools, AngleBetweenVectors
using Dierckx, VideoIO, CameraCalibrations, OhMyThreads, ImageTransformations, ImageDraw, Colors
using StaticArrays, CoordinateTransformations, Rotations
# using CairoMakie
using GLMakie
using AlgebraOfGraphics

const SV = SVector{2, Float64}

include("functions.jl")
include("calibrations.jl")

data_path = "data"

# Calibrations

calibrations = CSV.read(joinpath(data_path, "calib.csv"), DataFrame)

# data quality checks
@assert allunique(calibrations.calibration_id) "Calibration IDs not identical"
for row in eachrow(calibrations)
    @assert row.start < row.stop "start must occur before stop in calibration $(row.calibration_id)"
    file = joinpath(data_path, row.file)
    @assert isfile(file) "calibration file $file does not exist"
    # @assert that the time stamps are within the period of the video
end

@rselect!(calibrations, :calibration_id, :calib = calib(joinpath(data_path, :file), :start, :stop, :extrinsic))

# # save calibrations images for quality assesment
# suspect_calibration = "20220304_calibration.mov 2nd"
# check_calibration(calibrations, suspect_calibration)

# Tracking
df = CSV.read(joinpath(data_path, "runs.csv"), DataFrame)

# data quality checks
for row in eachrow(df)
    @assert row.start ≤ row.POI ≤ row.stop "POI is not within the track in row $row"
    @assert row.calibration_id ∈ calibrations.calibration_id "Calibration $calibration_id is missing from the calibrations file."
    file = joinpath(data_path, row.path, row.file)
    @assert isfile(file) "run file $file does not exist"
    # @assert that the time stamps are within the period of the video
end

df.track = tmap(track, joinpath.(data_path, df.path, df.file), df.start, df.stop)

# Combine the two
leftjoin!(df, calibrations, on = :calibration_id)
select!(df, Not(:calibration_id))

# calibrate and smooth the track
transform!(df, [:calib, :track] => ByRow(calibrate_smooth) => :spline)

# # check the calibrated-smoothed tracks and videos
# foreach(enumerate(eachrow(df))) do (i, row)
#     t, xy = row.track
#     save_vid(string(i), joinpath(data_path, row.path, row.file), row.calib.tform, t, xy, row.spline)
# end

# rotate and trim
transform!(df, [:start, :stop, :POI, :spline] => ByRow(rotate_trim) => :rotated)


colors = Dict(zip(unique(df.condition), Makie.wong_colors()))
@rtransform! df :color = colors[:condition]


fig = Figure()
ax = Axis(fig[1,1], aspect=DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)")
for r  in (30, 50)
    lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
end
for row in eachrow(df)
    _, i, xy = row.rotated
    lines!(ax, xy[1:i], color = :gray)
    lines!(ax, xy[i:end], color = row.color, label = row.condition)
end
axislegend(ax, merge = true)
save("fig.png", fig)


tbl = @chain df begin
    @transform :cordlength = cordlength.(:rotated)
    @transform :curvelength = curvelength.(:rotated)
    @transform :straightness = 1 .- :cordlength ./ :curvelength
    @transform :cumulative_angle = cumulative_angle.(:rotated)
    @select :condition :cordlength :curvelength :straightness :cumulative_angle
end

CSV.write("stats.csv", tbl)
