using Revise
using ProgressMeter
using CameraCalibrations
using DataFrames, Combinatorics
using CSV

video = "/home/yakir/mnt/dacke_lab_data/Data/Elin/Project_AllotheticVsIdiothetic_outdoors_50cm/20221125_outdoors_1.MTS"
ss = 17
stop = 50
t = stop - ss
r = 1
files = "intrinsic%03d.png"
run(`ffmpeg -loglevel 8 -ss $ss -i $video -t $t -r $r -vf format=gray,yadif=1,gblur=sigma=1 -pix_fmt gray $files`)

run(`ffmpeg -loglevel 8 -ss $ss -i $video -t $t -r $r -vf format=gray,yadif=1,scale=sar"*"iw:ih -pix_fmt gray $files`)

function detect_corners(gray, n_corners, flags)
    ret, _ = CameraCalibrations.cv2.findChessboardCorners(gray, n_corners, flags)
    return Bool(ret)
end

flags = [CameraCalibrations.cv2.CALIB_CB_ADAPTIVE_THRESH, CameraCalibrations.cv2.CALIB_CB_NORMALIZE_IMAGE, CameraCalibrations.cv2.CALIB_CB_FILTER_QUADS, CameraCalibrations.cv2.CALIB_CB_FAST_CHECK]
df = allcombinations(DataFrame, n_corners = [(i, j) for i in 3:8 for j in 3:8 if i ≥ j], flags = sum.(combinations(flags)))
df.n_files .= 0

files = filter(x -> last(splitext(x)) == ".png", readdir("."))

grays = map(files) do file
    img = CameraCalibrations.cv2.imread(file)
    CameraCalibrations.cv2.cvtColor(img, CameraCalibrations.cv2.COLOR_BGR2GRAY)
end

@showprogress for row in eachrow(df)
    row.n_files = count(grays) do gray
        detect_corners(gray, row.n_corners, row.flags)
    end
end

CSV.write("results.csv", df)

n_corners = (3, 3)
i = findfirst(gray -> detect_corners(gray, n_corners, flags[1]), grays)


c, e = fit(files, (4,6), 1)
