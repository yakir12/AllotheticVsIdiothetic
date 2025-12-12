using Test

get_total_rotationd(start, stop, cw, fullturns) = rad2deg(get_total_rotation(deg2rad(start), deg2rad(stop), cw, fullturns))

start = 0
stop = 90
cw = false
fullturns = 0
@test get_total_rotationd(start, stop, cw, fullturns) == 90

start = 0
stop = 90
cw = true
fullturns = 0
@test get_total_rotationd(start, stop, cw, fullturns) == -270

start = 90
stop = 0
cw = true
fullturns = 0
@test get_total_rotationd(start, stop, cw, fullturns) == -90

start = 90
stop = 0
cw = false
fullturns = 0
@test get_total_rotationd(start, stop, cw, fullturns) == 270

#####

start = convert2cartesiand(235)
stop = convert2cartesiand(0)
cw = true
fullturns = 0
@test get_total_rotationd(start, stop, cw, fullturns) == -125

convert2cartesian(rad) = 2π - rad
convert2cartesiand(deg) = 360 - deg

convert2cartesiand(90)


