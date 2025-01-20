# using GLMakie
using CairoMakie

using Dates, LinearAlgebra
using CSV, DataFrames, DataFramesMeta, CameraCalibrations
using Interpolations, StaticArrays, Dierckx, CoordinateTransformations, Rotations
using OhMyThreads

tosecond(t::T) where {T <: TimePeriod} = t / convert(T, Dates.Second(1))
tosecond(t::TimeType) = tosecond(t - Time(0))
tosecond(sec::Real) = sec
function totuple(x::AbstractString)
    if contains(x, '(')
        m = match(r"^\((\d+),\s*(\d+)\)$", x)
        Tuple{Int, Int}(parse.(Int, m.captures))
    else
        parse(Int, x)
    end
end
totuple(x) = x
function get_calibration(calibration_id)
    c = CameraCalibrations.load(joinpath(results_dir, calibration_id))
    rectification(c)
end

const results_dir = "../first/tracks and calibrations"

# TODO: make sure you're not doing things double (like for the same row_number)
# maybe appply the register earlier
# this might look a lot better with dataframesmeta or something
# maybe save in first the name of the csv files, complete, and the ones for the calibration too, just for ease of loading

runs = CSV.read(joinpath(results_dir, "runs.csv"), DataFrame)
transform!(runs, :POI => ByRow(passmissing(tosecond)); renamecols = false)
calibs = CSV.read(joinpath(results_dir, "calibs.csv"), DataFrame)
# minimal work requred
transform!(calibs, :calibration_id => ByRow(get_calibration) => :rectify)
transform!(calibs, [:rectify, :center_ij] => ByRow((f, c) -> passmissing(f)(totuple(c))) => :center)
transform!(calibs, [:rectify, :north_ij] => ByRow((f, c) -> passmissing(f)(totuple(c))) => :north)
select!(calibs, Cols(:calibration_id, :rectify, :center, :north))
leftjoin!(runs, calibs, on = :calibration_id)
select!(runs, Not(:runs_path, :start_location, :calibration_id, :fps, :target_width, :runs_file, :window_size))
rename!(runs, :runs_start => :start, :runs_stop => :stop)

# trying to remove that weird run
deleteat!(runs, 7:8)

function get_txy(tij_file, start, ::Missing, stop, rectify)
    tij = CSV.File(joinpath(results_dir, tij_file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    xy = rectify.(SVector{2, Int}.(tij.i, tij.j))
    return [(; t, xy)]
end
function get_txy(tij_file, start, POI, stop, rectify)
    tij = CSV.File(joinpath(results_dir, tij_file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    xy = rectify.(SVector{2, Int}.(tij.i, tij.j))
    i1 = something(findfirst(≥(start), t), 1)
    i2 = something(findfirst(≥(POI), t), length(t))
    i3 = something(findfirst(≥(stop), t), length(t))
    t1 = t[i1:i2]
    xy1 = xy[i1:i2]
    t2 = t[i2:i3]
    xy2 = xy[i2:i3]
    return [(; t = t1, xy = xy1), (; t = t2, xy = xy2)]
end

function smooth(t, xy)
    tp = ParametricSpline(t, stack(xy); k = 1, s = 5000)
    ts = range(t[1], t[end]; step = 1/1)
    sxy = SVector{2, Float64}.(tp.(ts))
    return (; t, xy = sxy)
end

function center_and_rotate(p1, p2)
    trans = Translation(-p1)
    x, y = normalize(p2 - p1)
    θ = π/2 - atan(y, x)
    rot = recenter(LinearMap(Angle2d(θ)), p1)
    return  trans ∘ rot
end
function register!(norths, centers, xys)
    p1 = first(centers)
    p2 = last(first(xys))
    trans = center_and_rotate(p1, p2)
    norths .= trans.(norths)
    centers .= trans.(centers)
    for xy in xys
        xy .= trans.(xy)
    end
end

function get_txy(tij_file, start, stop, rectify)
    tij = CSV.File(joinpath(results_dir, tij_file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    xy = rectify.(SVector{2, Int}.(tij.i, tij.j))
    return (; t, xy)
end
function concatenate(tij_files, starts, POIs, stops, rectifies)
    if length(POIs) > 1
        a = get_txy.(tij_files, starts, stops, rectifies)
        xy = vcat(getfield.(a, :xy)...)
        t = range(0; length = length(xy), step = step(a[1].t))
        l = length(a[1].t)
        (; t, xy, POI = t[l])
    else
        a = only(get_txy.(tij_files, starts, stops, rectifies))
        t = range(0; step = step(a.t), length = length(a.t))
        (; t, xy = a.xy, POI = POIs[])
    end
end


df = combine(groupby(runs, :run_id), [:north, :center] .=> Ref ∘ first, [:tij_file, :start, :POI, :stop, :rectify] => Ref ∘ concatenate => [:t, :xy, :poi], keepkeys = true, renamecols = false)
dropmissing!(df, :poi)
transform!(df, [:t, :xy] => ByRow(smooth); renamecols = false)
transform!(df, [:north, :center, :xy] => passmissing(register!); renamecols = false)


if isdir("tracks")
    rm("tracks", recursive=true)
end
mkpath("tracks")

@tasks for row in eachrow(df)
    fig  = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect())
    for r  in (300, 500)
        lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
    end
    lines!(ax, row.xy)
    save(joinpath("tracks", string(row.run_id, ".png")), fig)
end


#
#
#
# transform(groupby(runs, :run_id), [:tij_file, :start, :POI, :stop, :rectify] => ((f, s, p, ss, r) -> (@show length(p))))# => [:t, :xy, :POI])
#
#
# df = flatten(select(transform(runs, [:tij_file, :start, :POI, :stop, :rectify] => ByRow(get_txy) => :txy), Not(:start, :stop, :POI)), :txy)
# transform!(df, :txy => [:t, :xy])
# select!(df, Not(:txy))
# #
# # just ignore the runs that are split differently, there's just one
# subset!(transform!(groupby(df, :run_id), nrow), :nrow => ByRow(==(2)))
# @assert all(==(2), combine(groupby(df, :run_id), nrow).nrow)
# #
# transform!(df, [:t, :xy] => ByRow(smooth) => [:t, :xy])
# transform!(groupby(df, :run_id), [:north, :center, :xy] => register!; renamecols = false)
# #
#
# mkpath("tracks")
# gs = collect(groupby(df, :run_id))
# @tasks for i in 1:length(gs)
#     g = gs[i]
#     k = (; run_id = g.run_id[1])
#     fig  = Figure()
#     ax = Axis(fig[1,1], aspect = DataAspect())
#     for r  in (300, 500)
#         lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
#     end
#     for g1 in eachrow(g)
#         lines!(ax, g1.xy)
#     end
#     save(joinpath("tracks", string(k.run_id, ".png")), fig)
# end
#
#

















#
# using SimpTrack
# using Dates, LinearAlgebra
# using Dierckx, VideoIO, CameraCalibrations, OhMyThreads, ImageTransformations, ImageDraw, Colors
# using StaticArrays, CoordinateTransformations, Rotations
# # using CairoMakie
# using GLMakie
# using AlgebraOfGraphics
#
# const SV = SVector{2, Float64}
#
# include("functions.jl")
# include("calibrations.jl")
#
# calib_file = "/home/yakir/tmp/Elin_tracks/calib.csv"
# runs_file = "/home/yakir/tmp/Elin_tracks/runs.csv"
#
# # Calibrations
#
# calibrations = CSV.read(calib_file, DataFrame)
#
# # data quality checks
# @assert allunique(calibrations.calibration_id) "Calibration IDs not identical"
# for row in eachrow(calibrations)
#     @assert row.start < row.stop "start must occur before stop in calibration $(row.calibration_id)"
#     file = joinpath(row.path, row.file)
#     @assert isfile(file) "calibration file $file does not exist"
#     # @assert that the time stamps are within the period of the video
# end
#
# @rselect!(calibrations, :calibration_id, :calib = calib(:calibration_id, joinpath(:path, :file), :start, :stop, :extrinsic))
#
# # # save calibrations images for quality assesment
# # suspect_calibration = "20220304_calibration.mov 2nd"
# # check_calibration(calibrations, suspect_calibration)
#
# # Tracking
# df = CSV.read(runs_file, DataFrame)
#
# # data quality checks
# for row in eachrow(df)
#     @assert row.start ≤ row.POI ≤ row.stop "POI is not within the track in row $row"
#     @assert row.calibration_id ∈ calibrations.calibration_id "Calibration $calibration_id is missing from the calibrations file."
#     file = joinpath(row.path, row.file)
#     @assert isfile(file) "run file $file does not exist"
#     # @assert that the time stamps are within the period of the video
# end
#
# df.track = tmap(track, joinpath.(df.path, df.file), df.start, df.stop)
#
# # Combine the two
# leftjoin!(df, calibrations, on = :calibration_id)
# select!(df, Not(:calibration_id))
#
#
#
#
#
# # s = 100
# # row = first(eachrow(df))
# # runi, clib, trck, start, stop, POI = (1, row.calib, row.track, row.start, row.stop, row.POI)
# # tfm = get_transformation(clib, trck, start, POI; s)
# # # xy = tfm.(range(start, stop, step = Millisecond(100)))
# # # lines(xy)
# #
# # t = range(start, stop, step = Second(1))
# #
# # trajectory_length = cumsum(norm.(diff(tfm.(t))))
# # ca = [cumulative_angle(tfm, start, tstop) for tstop in t[2:end]]
# # total_cumulative_angle = round(Int, cumulative_angle(tfm, start, stop))
#
#
# ns = 5
# ss = round.(Int, range(50, 250, ns))
# function explore(runi, calib, track, start, stop, POI)
#     w = 300
#     fig = Figure(size=(2w, ns*w))
#     axs1 = []
#     axs2 = []
#     for (j, s) in enumerate(ss)
#         tfm = get_transformation(calib, track, start, POI; s)
#         t = range(start, stop, step = Second(2))
#         trajectory_length = cumsum(norm.(diff(tfm.(t))))
#         ca = [cumulative_angle(tfm, start, tstop) for tstop in t[2:end]]
#         total_cumulative_angle = round(Int, cumulative_angle(tfm, start, stop))
#         ax = Axis(fig[j, 1], aspect=DataAspect(), ylabel = "s = $s")
#         for r  in (30, 50)
#             lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
#         end
#         xy = tfm.(range(start, POI, step = Millisecond(100)))
#         lines!(ax, xy)
#         xy = tfm.(range(POI, stop, step = Millisecond(100)))
#         lines!(ax, xy)
#         j ≠ ns && hidexdecorations!(ax, grid = false, minorgrid = false)#, ticks = false, minorticks = false)
#         push!(axs1, ax)
#         ax = Axis(fig[j, 2], xlabel = "Trajectory length (cm)", yaxisposition = :right, aspect = AxisAspect(1), ylabel = "Cumulative angle (°)", title = "Final cumulative angle = $total_cumulative_angle °")
#         lines!(ax, trajectory_length, ca)
#         j ≠ ns && hidexdecorations!(ax, grid = false, minorgrid = false)#, ticks = false, minorticks = false)
#         push!(axs2, ax)
#     end
#     # linkaxes!(filter(x -> x isa Axis, fig.content)...)
#     linkaxes!(axs1...)
#     linkaxes!(axs2...)
#     save("run $runi.png", fig)
# end
# foreach(enumerate(eachrow(df))) do (runi, row)
#     explore(runi, row.calib, row.track, row.start, row.stop, row.POI)
# end
#
# # n = 10
# # α = rand(-10:10, n)
# #
#
# # TODO: ask Vishaal about the noise effect on cumulative turn angle
# function angle_between(p1, p2)
#     θ = acos(dot(p1, p2) / norm(p1) / norm(p2))
#     return -sign(cross(p1, p2))*θ
# end
# function fun(xy)
#     δ = diff(xy)
#     α = sum(splat(angle_between) ∘ reverse, partition(δ, 2, 1))
#     rad2deg(α)
#     # pushfirst!(α, atan(reverse(δ[1])...))
#     # rad2deg(sum(splat(angle_between), partition(δ, 2, 1)))
# end
# n = 1000
# xy = [SV(reverse(sincos(θ))) for θ in range(0, π, n)]
# fig, ax, h = lines(xy, label = "clean",  axis = (;aspect = DataAspect()))
# xy1 = xy .+ 0.01randn(SV, n)
# lines!(ax, xy1, label = "with noise")
# xy2 = xy1[round.(Int, range(1, n, 15))]
# # scatterlines!(ax, xy2, label = "subsampled")
# δ = diff(xy2)
# α = map(splat(angle_between) ∘ reverse, partition(δ, 2, 1))
# pushfirst!(α, atan(reverse(δ[1])...))
# l = [SV(reverse(sincos(i))) for i in cumsum(α)]
# arrows!(ax, Point2f.(xy2[1:end-1]), l .* norm.(δ), label = "subsampled")
# axislegend(ax)
# fun(xy)
# fun(xy1)
# fun(xy2)
# #
# # using Statistics
# #
# # steps = range(1, 100, 101)
# # m = 1000
# # n = 1000
# # Δ = zeros(length(steps), m)
# # for mi in 1:m
# #     xy = [SV(reverse(sincos(θ))) for θ in range(0, π, n)]
# #     xy1 = xy .+ 0.01randn(SV, n)
# #     Δ[:,mi] = map(steps) do step
# #         xy2 = xy1[round.(Int, range(1, n; step))]
# #         abs(180 - fun(xy2))
# #     end
# # end
# # lines(steps, vec(mean(Δ, dims=2)))
# #
# #
#
#
#
#
#
#
#
#
#
# #
# # δα = [0, 0, 45, 45, 45 , 45, 0]
# # α = cumsum(δα)
# # l = [SV(reverse(sincosd(i))) for i in α]
# # xy = cumsum(l)
# # arrows!(ax, Point2f.(xy[1:end-1]), l[2:end])
# #
# # Δ = diff(α)
# # sum(Δ[2:end])
# #
# #
# # cumulative_angle(xy)
# #
