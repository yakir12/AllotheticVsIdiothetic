using Dates, LinearAlgebra
using DataFrames, CSV, CameraCalibrationMeta, OhMyThreads, Dierckx, StaticArrays, CoordinateTransformations, Rotations
using GLMakie
const SV = SVector{2, Float64}

include("functions.jl")

data_path = "../first"

runs_file = "../../data/runs.csv"
df = CSV.read(runs_file, DataFrame)
select!(df, Not(:path, :file))

foreach(df.calibration_id) do calibration_id
    @assert isfile(joinpath(data_path, "$calibration_id.calib")) "Calibration $calibration_id is missing from the calibrations file."
end

transform!(df, :calibration_id => ByRow(name -> CameraCalibrationMeta.load(joinpath(data_path, "$name.calib"))) => :calib)
select!(df, Not(:calibration_id))

df.track = tmap(eachrow(df)) do row
    name = rownumber(row)
    CSV.read(joinpath(data_path, "$name.csv"), DataFrame)
end

s, k = (100, 2)
step = Millisecond(33)
radii = (30, 50)

transform!(df, [:start, :POI, :calib, :track] => ByRow(track_function) => :track)
select!(df, Not(:calib))


fig = Figure()
for (i, (k, grp)) in enumerate(pairs(groupby(df, :condition)))
    ax = Axis(fig[i,1], aspect=DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)", title = k.condition)
    for r  in radii
        lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
    end
    for row in eachrow(grp)
        ts = range(row.start, row.POI; step)
        before = row.track.(ts)
        lines!(ax, before, color = :gray)
        ts = range(row.POI, row.stop; step)
        after = row.track.(ts)
        filter!(≤(last(radii)) ∘ norm, after)
        lines!(ax, after, color = :black)
    end
end

save("fig.png", fig)



tbl = @chain df begin
    @transform :cordlength = cordlength.(:rotated)
    @transform :curvelength = curvelength.(:rotated)
    @transform :straightness = 1 .- :cordlength ./ :curvelength
    @transform :cumulative_angle = cumulative_angle.(:rotated)
    @select :condition :cordlength :curvelength :straightness :cumulative_angle
end

CSV.write("stats.csv", tbl)




