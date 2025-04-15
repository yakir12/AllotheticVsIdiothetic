using Revise
using Fromage

data_path = "/home/yakir/mnt/dacke_lab_data/Data/Elin/Project_AllotheticVsIdiothetic_outdoors_50cm"

if isdir("tracks and calibrations")
    rm("tracks and calibrations", recursive=true)
end

main(data_path; delim = ';')
