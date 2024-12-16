using CSV, DataFrames
results_dir = "tracks and calibrations"

runs = CSV.read(joinpath("..", "first", results_dir, "runs.csv"), DataFrame)
calib = CSV.read(joinpath("..", "first", results_dir, "calib.csv"), DataFrame)



using SimpTrack
using Dates, LinearAlgebra
using Dierckx, VideoIO, CameraCalibrations, OhMyThreads, ImageTransformations, ImageDraw, Colors
using StaticArrays, CoordinateTransformations, Rotations
# using CairoMakie
using GLMakie
using AlgebraOfGraphics

const SV = SVector{2, Float64}

include("functions.jl")
include("calibrations.jl")

calib_file = "/home/yakir/tmp/Elin_tracks/calib.csv"
runs_file = "/home/yakir/tmp/Elin_tracks/runs.csv"

# Calibrations

calibrations = CSV.read(calib_file, DataFrame)

# data quality checks
@assert allunique(calibrations.calibration_id) "Calibration IDs not identical"
for row in eachrow(calibrations)
    @assert row.start < row.stop "start must occur before stop in calibration $(row.calibration_id)"
    file = joinpath(row.path, row.file)
    @assert isfile(file) "calibration file $file does not exist"
    # @assert that the time stamps are within the period of the video
end

@rselect!(calibrations, :calibration_id, :calib = calib(:calibration_id, joinpath(:path, :file), :start, :stop, :extrinsic))

# # save calibrations images for quality assesment
# suspect_calibration = "20220304_calibration.mov 2nd"
# check_calibration(calibrations, suspect_calibration)

# Tracking
df = CSV.read(runs_file, DataFrame)

# data quality checks
for row in eachrow(df)
    @assert row.start ≤ row.POI ≤ row.stop "POI is not within the track in row $row"
    @assert row.calibration_id ∈ calibrations.calibration_id "Calibration $calibration_id is missing from the calibrations file."
    file = joinpath(row.path, row.file)
    @assert isfile(file) "run file $file does not exist"
    # @assert that the time stamps are within the period of the video
end

df.track = tmap(track, joinpath.(df.path, df.file), df.start, df.stop)

# Combine the two
leftjoin!(df, calibrations, on = :calibration_id)
select!(df, Not(:calibration_id))



# s = 100
# row = first(eachrow(df))
# runi, clib, trck, start, stop, POI = (1, row.calib, row.track, row.start, row.stop, row.POI)
# tfm = get_transformation(clib, trck, start, POI; s)
# # xy = tfm.(range(start, stop, step = Millisecond(100)))
# # lines(xy)
#
# t = range(start, stop, step = Second(1))
#
# trajectory_length = cumsum(norm.(diff(tfm.(t))))
# ca = [cumulative_angle(tfm, start, tstop) for tstop in t[2:end]]
# total_cumulative_angle = round(Int, cumulative_angle(tfm, start, stop))


ns = 5
ss = round.(Int, range(50, 250, ns))
function explore(runi, calib, track, start, stop, POI)
    w = 300
    fig = Figure(size=(2w, ns*w))
    axs1 = []
    axs2 = []
    for (j, s) in enumerate(ss)
        tfm = get_transformation(calib, track, start, POI; s)
        t = range(start, stop, step = Second(2))
        trajectory_length = cumsum(norm.(diff(tfm.(t))))
        ca = [cumulative_angle(tfm, start, tstop) for tstop in t[2:end]]
        total_cumulative_angle = round(Int, cumulative_angle(tfm, start, stop))
        ax = Axis(fig[j, 1], aspect=DataAspect(), ylabel = "s = $s")
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
        xy = tfm.(range(start, POI, step = Millisecond(100)))
        lines!(ax, xy)
        xy = tfm.(range(POI, stop, step = Millisecond(100)))
        lines!(ax, xy)
        j ≠ ns && hidexdecorations!(ax, grid = false, minorgrid = false)#, ticks = false, minorticks = false)
        push!(axs1, ax)
        ax = Axis(fig[j, 2], xlabel = "Trajectory length (cm)", yaxisposition = :right, aspect = AxisAspect(1), ylabel = "Cumulative angle (°)", title = "Final cumulative angle = $total_cumulative_angle °")
        lines!(ax, trajectory_length, ca)
        j ≠ ns && hidexdecorations!(ax, grid = false, minorgrid = false)#, ticks = false, minorticks = false)
        push!(axs2, ax)
    end
    # linkaxes!(filter(x -> x isa Axis, fig.content)...)
    linkaxes!(axs1...)
    linkaxes!(axs2...)
    save("run $runi.png", fig)
end
foreach(enumerate(eachrow(df))) do (runi, row)
    explore(runi, row.calib, row.track, row.start, row.stop, row.POI)
end

# n = 10
# α = rand(-10:10, n)
#

# TODO: ask Vishaal about the noise effect on cumulative turn angle
function angle_between(p1, p2)
    θ = acos(dot(p1, p2) / norm(p1) / norm(p2))
    return -sign(cross(p1, p2))*θ
end
function fun(xy)
    δ = diff(xy)
    α = sum(splat(angle_between) ∘ reverse, partition(δ, 2, 1))
    rad2deg(α)
    # pushfirst!(α, atan(reverse(δ[1])...))
    # rad2deg(sum(splat(angle_between), partition(δ, 2, 1)))
end
n = 1000
xy = [SV(reverse(sincos(θ))) for θ in range(0, π, n)]
fig, ax, h = lines(xy, label = "clean",  axis = (;aspect = DataAspect()))
xy1 = xy .+ 0.01randn(SV, n)
lines!(ax, xy1, label = "with noise")
xy2 = xy1[round.(Int, range(1, n, 15))]
# scatterlines!(ax, xy2, label = "subsampled")
δ = diff(xy2)
α = map(splat(angle_between) ∘ reverse, partition(δ, 2, 1))
pushfirst!(α, atan(reverse(δ[1])...))
l = [SV(reverse(sincos(i))) for i in cumsum(α)]
arrows!(ax, Point2f.(xy2[1:end-1]), l .* norm.(δ), label = "subsampled")
axislegend(ax)
fun(xy)
fun(xy1)
fun(xy2)
#
# using Statistics
#
# steps = range(1, 100, 101)
# m = 1000
# n = 1000
# Δ = zeros(length(steps), m)
# for mi in 1:m
#     xy = [SV(reverse(sincos(θ))) for θ in range(0, π, n)]
#     xy1 = xy .+ 0.01randn(SV, n)
#     Δ[:,mi] = map(steps) do step
#         xy2 = xy1[round.(Int, range(1, n; step))]
#         abs(180 - fun(xy2))
#     end
# end
# lines(steps, vec(mean(Δ, dims=2)))
#
#









#
# δα = [0, 0, 45, 45, 45 , 45, 0]
# α = cumsum(δα)
# l = [SV(reverse(sincosd(i))) for i in α]
# xy = cumsum(l)
# arrows!(ax, Point2f.(xy[1:end-1]), l[2:end])
#
# Δ = diff(α)
# sum(Δ[2:end])
#
#
# cumulative_angle(xy)
#
