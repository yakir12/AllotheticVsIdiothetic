using Revise
using Fromage

data_path = "/home/yakir/mnt/Data/Elin/AllotheticVsIdiothetic"

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end

main(data_path)


