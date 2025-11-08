using Revise
using Fromage

data_path = "/home/yakir/mnt/Data/Elin/Project_AllotheticVsIdiothetic_indoors"

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end

main(data_path; delim = ';')
