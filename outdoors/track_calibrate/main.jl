using Revise
using Fromage

data_path = "/home/yakir/mnt/Data/Elin/Project_AllotheticVsIdiothetic_outdoors"

rm.(filter(startswith("jl_"), readdir(".")), force = true, recursive = true)

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end

# main(data_path)
main(data_path, todo = ["20230111_mirror_dance01_spontaneous_ID03.MTS"])
