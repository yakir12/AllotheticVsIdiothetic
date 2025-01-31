using AlgebraOfGraphics, GLMakie
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

const results_dir = "../indoors/tracks and calibrations"

# TODO: make sure you're not doing things double (like for the same row_number)
# maybe appply the register earlier
# this might look a lot better with dataframesmeta or something
# maybe save in first the name of the csv files, complete, and the ones for the calibration too, just for ease of loading

runs = CSV.read(joinpath(results_dir, "runs.csv"), DataFrame)
# transform!(runs, :POI => ByRow(passmissing(tosecond)); renamecols = false)
calibs = CSV.read(joinpath(results_dir, "calibs.csv"), DataFrame)
# minimal work requred
transform!(calibs, :calibration_id => ByRow(get_calibration) => :rectify)
transform!(calibs, [:rectify, :center_ij] => ByRow((f, c) -> passmissing(f)(totuple(c))) => :center)
transform!(calibs, [:rectify, :north_ij] => ByRow((f, c) -> passmissing(f)(totuple(c))) => :north)
select!(calibs, Cols(:calibration_id, :rectify, :center, :north))
leftjoin!(runs, calibs, on = :calibration_id)
select!(runs, Not(:runs_path, :start_location, :calibration_id, :fps, :target_width, :runs_file, :window_size))

function get_rotation(xy)
    # θ = -atan(reverse(sum(normalize, xy))...)
    # θ = -atan(reverse(sum(normalize, xy))...)
    p1 = first(xy)
    p2 = last(xy)
    x, y = normalize(p2 - p1)
    θ = π/2 - atan(y, x)
    # p = first.(xy) \ last.(xy)
    # θ = atan(p)
    # if sum(last, xy) > 0
    #     θ += π
    # end
    # θ = π/2 - θ
    LinearMap(Angle2d(θ))
end
function get_txy(tij_file, rectify, poi)
    tij = CSV.File(joinpath(results_dir, tij_file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    poi_index = something(findfirst(>(poi), t), length(t))
    ij = SVector{2, Int}.(tij.i, tij.j)
    xy = rectify.(ij)
    tp = ParametricSpline(t, stack(xy); k = 2, s = 300)
    xy .= SVector{2, Float64}.(tp.(t))
    trans = Translation(-xy[1])
    xy .= trans.(xy)
    rot = get_rotation(xy[2:poi_index])
    xy .= rot.(xy)
    (; poi_index, t, xy)
end
transform!(runs, [:tij_file, :rectify, :poi] => ByRow(get_txy) => [:poi_index, :t, :xy])
function plotone(run_id, xy, poi_index)
    fig  = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect(), autolimitaspect = 1, title = string(run_id), limits = ((-60, 60), (-60, 60)))
    for r  in (30, 50)
        lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
    end
    lines!(ax, xy[1:poi_index])
    lines!(ax, xy[poi_index:end])
    # θ = atan(reverse(sum(normalize, xy[2:poi_index]))...)
    # arrows!(ax, [0], [0], [cos(θ)], [sin(θ)], lengthscale = 10)
    # a, b = [first.(xy[1:poi_index]) ones(poi_index)] \ last.(xy[1:poi_index])
    # ablines!(ax, b, a, color = :gray)
    save(joinpath("tracks", string(run_id, ".png")), fig)
end
if isdir("tracks")
    rm("tracks", recursive=true)
end
mkpath("tracks")
CairoMakie.activate!()
@tasks for row in eachrow(runs)
    plotone(row.run_id, row.xy, row.poi_index)
end






GLMakie.activate!()

words = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"]

transform!(runs, :dance => ByRow(x -> x ? "dance" : "no dance"), :at_run => ByRow(x -> words[x]), renamecols = false)

function cropto(xy, l)
    i = something(findfirst(>(l) ∘ norm, xy), length(xy))
    xy[1:i-1]
end

transform!(runs, :xy => ByRow(xy -> cropto(xy, 50)); renamecols = false)

df1 = flatten(runs, :xy)
transform!(df1, :xy => [:x, :y])

plt = data(df1) * mapping(:x => "X (cm)", :y => "Y (cm)", group=:run_id => nonnumeric, col = :dance => nonnumeric, row = :at_run => nonnumeric => "at run", color = :light) * visual(Lines)
fig = draw(plt; axis=(aspect=1, ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save("figure1.png", fig)

# display(fig)



plt = data(df1) * mapping(:x => "X (cm)", :y => "Y (cm)", group=:run_id => nonnumeric, col = :dance => nonnumeric, row = :light => nonnumeric => "at run", color = :at_run) * visual(Lines)
fig = draw(plt; axis=(aspect=1, ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save("figure2.png", fig)

function get_spline(tij_file, rectify, poi)
    tij = CSV.File(joinpath(results_dir, tij_file))
    t = range(tij.t[1], tij.t[end], length = length(tij))
    poi_index = something(findfirst(>(poi), t), length(t))
    ij = SVector{2, Int}.(tij.i, tij.j)
    xy = rectify.(ij)
    (;t, spl = ParametricSpline(t, stack(xy); k = 2, s = 300), poi_index)
end
df = select(runs, [:tij_file, :rectify, :poi] => ByRow(get_spline) => [:t, :spl, :poi_index], Cols(:run_id, :dance, :light, :at_run))

using Statistics
mean_angle(θ) = angle(mean(exp, θ*im))

function unwrap!(x, period = 2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end

function plot_direction(poi_index, spl, t, run_id)
    # row = df[2,:]
    # poi_index = row.poi_index
    # spl = row.spl
    # t = row.t
    # run_id = row.run_id
    # GLMakie.activate!()
    der = derivative.(Ref(spl), t)
    fig = Figure()
    ax = Axis(fig[1,1], autolimitaspect = 1, aspect = DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)")
    lines!(ax, Point2f.(spl.(t[1:poi_index])))
    lines!(ax, Point2f.(spl.(t[poi_index:end])))
    ax = Axis(fig[2,1], limits = (extrema(t), nothing), xlabel = "Time (sec)", ylabel = "Direction (°)")
    θ = [atan(reverse(d)...) for d in der]
    unwrap!(θ)
    lines!(ax, t[1:poi_index], rad2deg.(θ[1:poi_index]))
    lines!(ax, t[poi_index:end], rad2deg.(θ[poi_index:end]))
    θs = [mean_angle(θ[i]) for i in (1:poi_index, poi_index:length(t))]
    lines!(ax, t[[1, poi_index, poi_index, end]], [fill(rad2deg(θs[1]), 2); fill(rad2deg(θs[2]), 2)], color=:gray)
    ax = Axis(fig[3,1], limits = (extrema(t), nothing), xlabel = "Time (sec)", ylabel = "Difference in direction (°)")
    lines!(ax, t[2:end], diff(θ))
    display(fig)
    # save(joinpath("directions", string(run_id, ".png")), fig)
end

row = df[15,:]
plot_direction(row.poi_index, row.spl, row.t, row.run_id)


if isdir("directions")
    rm("directions", recursive=true)
end
mkpath("directions")
CairoMakie.activate!()
@tasks for row in eachrow(df)
    plot_direction(row.poi_index, row.spl, row.t, row.run_id)
end

    GLMakie.activate!()
function turning_event(t, spl, poi_index)
    der = derivative.(Ref(spl), t)
    θ = [atan(reverse(d)...) for d in der]
    # d, i = findmax(diff(θ[poi_index:end]))
    d, i = findmax(abs.(diff(θ)))
    # (; d, dt = t[i + poi_index - 1] - t[poi_index])
    (; d = rad2deg(d), dt = t[i] - t[poi_index])
end
df = select(runs, [:tij_file, :rectify, :poi] => ByRow(get_spline) => [:t, :spl, :poi_index], Cols(:run_id, :dance, :light, :at_run));
transform!(df, [:t, :spl] => ByRow(turning_event) => [:d, :dt]);
plt = data(df) * mapping(:dt => "Time from POI (sec)", :d => "Turn (°)", col = :dance => nonnumeric, row = :light => nonnumeric => "at run", color = :at_run) * visual(Scatter)
fig = draw(plt)





# display(fig)


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
