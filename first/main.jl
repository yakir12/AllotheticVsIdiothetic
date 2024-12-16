# First
using Revise
using Fromage

# data_path = "/media/yakir/205D-185A/Bastien's"
# data_path = "/home/yakir/mnt/dacke_lab_data/Data/Elin/Project_AllotheticVsIdiothetic_indoors"
data_path = "/media/yakir/205D-185A/Bastien's NEW"
main(data_path)


# Second
using DataFrames, CSV, CameraCalibrations, StaticArrays, GLMakie
using AstroLib, TimeZone, TOMLs

path = "/home/yakir/new_projects/Elin/code/first/fromage/tracks and calibrations"
runs = CSV.read(joinpath(path, "runs.csv"), DataFrame)

# subset!(runs, :run_number => ByRow(==(4)))

# calibs = CSV.read(joinpath(path, "calib.csv"), DataFrame)
# leftjoin!(runs, calibs, on = :calibration_id, makeunique = true)


fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
for r in (15, 30, 150, 300)
    lines!(Circle(zero(Point2f), r), color = :gray)
end
for row in eachrow(runs)
    c = CameraCalibrations.load(joinpath(path, string(row.calibration_id, ".calib")))
    xyt = CSV.read(joinpath(path, string(row.row_number, ".csv")), DataFrame)
    sort!(xyt, :t)
    xy = Point2f.(pop.(c.(RowCol.(xyt.x, xyt.y))))
    xy .-= first(xy)
    lines!(ax, xy)
end
display(fig)
