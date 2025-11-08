using Revise
using Fromage

data_path = "/home/yakir/mnt/Data/Elin/Project_AllotheticVsIdiothetic_indoors"

rm.(filter(startswith("jl_"), readdir(".")), force = true, recursive = true)

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end

main(data_path)
