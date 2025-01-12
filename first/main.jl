# First
using Revise
using Fromage

# data_path = "/media/yakir/205D-185A/Bastien's"
data_path = "/home/yakir/mnt/dacke_lab_data/Data/Elin/Project_AllotheticVsIdiothetic_indoors" # took about 17 minutes
# data_path = "/home/yakir/mnt/dacke_lab_data/Data/Elin/Project_AllotheticVsIdiothetic_outdoors_50cm"
# data_path = "/home/yakir/mnt/dacke_lab_data/Data/2024/South Africa 2024 Nov/Bastien"
# data_path = "/media/yakir/3A82C7D782C79633/Bastien's backup"
# data_path = "/media/yakir/205D-185A/Bastien's NEW"
# data_path = "/media/yakir/E/data_lamarcki/Group_B/Round_1/12_11"

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end

main(data_path)


#
# files = readdir("jl_5Tly1A", join = true)
# n_corners = (6, 9)
# checker_size = 4
# aspect = 1
#
# c, ϵ = fit(files, n_corners, checker_size; aspect)
# @assert !isnothing(findfirst(contains(r"extrinsic"), c.files)) "fail"
#
# # Second
# using DataFrames, CSV, CameraCalibrations, StaticArrays, GLMakie
# using AstroLib, TimeZones, TOML
#
# path = "/home/yakir/new_projects/Elin/code/first/fromage/tracks and calibrations"
# runs = CSV.read(joinpath(path, "runs.csv"), DataFrame)
#
# # subset!(runs, :run_number => ByRow(==(4)))
#
# # calibs = CSV.read(joinpath(path, "calib.csv"), DataFrame)
# # leftjoin!(runs, calibs, on = :calibration_id, makeunique = true)
#
#
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# for r in (15, 30, 150, 300)
#     lines!(Circle(zero(Point2f), r), color = :gray)
# end
# for row in eachrow(runs)
#     c = CameraCalibrations.load(joinpath(path, string(row.calibration_id, ".calib")))
#     xyt = CSV.read(joinpath(path, string(row.row_number, ".csv")), DataFrame)
#     sort!(xyt, :t)
#     xy = Point2f.(pop.(c.(RowCol.(xyt.x, xyt.y))))
#     xy .-= first(xy)
#     lines!(ax, xy)
# end
# display(fig)
#
# x,y = (439, 44)
# w, h = (988, 932)
# 988:932:439:44
# scale=320:-1
#
# crop=745:920:328:50
