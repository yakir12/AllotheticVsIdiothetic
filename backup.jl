#! julia

file = only(ARGS)

root = "/home/yakir/mnt/dacke_lab_data/Data/Elin"
rp = relpath(file, root)
org = joinpath(root, "backup", rp)
mv(file, org)
run(`ffmpeg -i $org -c:v libx264 -crf 18 -preset slow -c:a copy -movflags use_metadata_tags -map_metadata 0 $file`)
