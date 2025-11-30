using Revise
using Fromage

data_path = "/home/yakir/mnt/Data/Elin/Project_AllotheticVsIdiothetic_outdoors"

rm.(filter(startswith("jl_"), readdir(".")), force = true, recursive = true)

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end
if isdir("tracks")
    rm("tracks", recursive=true)
end
if isdir("calibrations")
    rm("calibrations", recursive=true)
end

Fromage.only_calibrate(data_path, todo = ["20221124_mirror_dance10_spontaneous_Bela-Bela_23.MTS"])
