function seek_snap(path, vid, t, aspect, name = t)
    seek(vid, t)
    img = read(vid)
    save(joinpath(path, "$name.png"), warp(img, aspect))
end

calib(calibration_id, file, start::Time, stop::Time, extrinsic::Time; checker_size::Real = 3.9, temporal_step::Real = 2, n_corners::NTuple{2, Int} = (5,8)) = calib(calibration_id, file, tosecond(start), tosecond(stop), tosecond(extrinsic), checker_size, temporal_step, n_corners)

function calib(calibration_id, file, start::Real, stop::Real, extrinsic::Real, checker_size, temporal_step, n_corners)
    vid = openvideo(file, target_format=VideoIO.AV_PIX_FMT_GRAY8)
    aspect = LinearMap(SDiagonal{2, Float64}([1, 1/VideoIO.aspect_ratio(vid)]))
    read(vid)
    t₀ = gettime(vid)
    start += t₀
    stop += t₀
    extrinsic += t₀
    ts = range(start, stop, step = temporal_step)
    # mktempdir(@__DIR__; cleanup=false) do path
    path = mktempdir(@__DIR__; cleanup=false)
    foreach(t -> seek_snap(path, vid, t, aspect), ts)
    seek_snap(path, vid, extrinsic, aspect, "extrinsic")

    files = readdir(path; join=true)

    c = Calibration(files, n_corners, checker_size)
    extrinsic_index = findfirst(contains(r"extrinsic"), c.files)

    if isnothing(extrinsic_index)
        @error "The extrinsic image in calibration $calibration_id, is unusable! Please choose a different time point for the extrinsic image."
    end

    tform =  aspect ∘ c.real2image[extrinsic_index] ∘ Base.Fix2(push, 0)
    itform = pop ∘ c.image2real[extrinsic_index] ∘ inv(aspect)

    return (; tform, itform, c)

end

check_calibration(suspect_calibration) = CameraCalibrations.plot(only(subset(calibrations, :calibration_id => ByRow(==(suspect_calibration))).calib).c)

function save_vid(name, file, tform, t, xy, spl; axs = (-100:100, -100:100))

    smooth_xy = SV.(spl.(t))

    vid = openvideo(file)
    fig = Figure()#size=(h, round(Int, h * sz[2] / sz[1] / VideoIO.aspect_ratio(vid))), figure_padding=0)
    ax = Axis(fig[1,1], aspect = DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)")
    imgw = warp(read(vid), tform, axs)
    img = Observable(collect(imgw))
    image!(ax, GLMakie.ClosedInterval.(axs)..., img)
    point = Observable(smooth_xy[1])
    scatter!(ax, point, marker='+', color=:red)
    t₀ = gettime(vid)
    t .+= t₀
    seek(vid, t[1])
    framerate = round(Int, 2length(t)/(t[end] - t[1]))
    record(fig, "$name.mp4", zip(vid, xy, smooth_xy); framerate) do (frame, (x, y), smooth_xyi)
        ImageDraw.draw!(frame, CirclePointRadius(y, x, 20, thickness = 10, fill = false), colorant"blue")
        img[] = collect(warp(frame, tform, axs))
        point[] = smooth_xyi
    end

end


