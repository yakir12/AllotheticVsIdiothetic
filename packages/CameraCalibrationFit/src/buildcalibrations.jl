function obj2img(k, Rs, ts, frow, fcol, crow, ccol, checker_size)
    intrinsic = AffineMap(SDiagonal(frow, fcol), SVector(crow, ccol))
    distort(rc) = lens_distortion(rc, k)
    extrinsics = AffineMap.(Base.splat(RotationVec).(Rs), SVector{3, Float64}.(ts))
    scale = LinearMap(SDiagonal{3}(I/checker_size))
    return intrinsic, distort, extrinsics, scale
end

struct CalibrationFit
    files # image files
    n_files # how many
    n_corners # the number of corners in each of the two dimensions of the checkerboard
    checker_size # physical size of the checker (e.g. in cm)
    objpoints # ideal real world coordinates of the checkerboard corners
    imgpointss # the detetcted corners in each of the images, in pixel coordinates
    k # radial lens distortion coefficient
    real2image # the transformation that converts real world to image coordinates
    image2real # the transformation that converts image coordinates to real world
    sz # the dimensions of the images
    intrinsic::CameraCalibrationMeta.AM
    extrinsics::Vector{CameraCalibrationMeta.AMext}
    scale::LinearMap{SDiagonal{3, Float64}}
end

function CalibrationFit(files, n_corners, checker_size, with_distortion)
    files, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol = detect_fit(unique(files), n_corners, with_distortion)
    objpoints .*= checker_size
    intrinsic, distort, extrinsics, scale = obj2img(k, Rs, ts, frow, fcol, crow, ccol, checker_size)
    real2image = .∘(Ref(intrinsic), distort, Ref(PerspectiveMap()), extrinsics, Ref(scale))
    inv_scale, inv_extrinsics, inv_perspective_maps, inv_distort, inv_intrinsic = img2obj(intrinsic, extrinsics, scale, k)
    image2real = .∘(Ref(inv_scale), inv_extrinsics, inv_perspective_maps, inv_distort, Ref(inv_intrinsic))
    return CalibrationFit(files, length(files), n_corners, checker_size, objpoints, imgpointss, k, real2image, image2real, sz, intrinsic, extrinsics, scale)
end

function CameraCalibrationMeta.Calibration(cf::CalibrationFit)
    CameraCalibrationMeta.Calibration(cf.intrinsic, cf.extrinsics, cf.scale, cf.k, cf.files, cf.real2image, cf.image2real)
end

"""
    Calibration(files, n_corners, checker_size, extrinsic_index)
Build a calibration object. `files` are the image files of the checkerboard. `n_corners` is a tuple of the number of corners in each of the sides of the checkerboard. `checker_size` is the physical size of the checker (e.g. in cm). `with_distortion` controls if radial lens distortion is included in the model or not.
"""
function fit(files, n_corners, checker_size; with_distortion = true, improve_n = nothing, improve_threshold = nothing, inverse_samples = 100, with_plot = false)
    cf = CalibrationFit(files, n_corners, checker_size, with_distortion)
    cf = improve(cf, improve_n, improve_threshold)
    # ϵ = calculate_errors(cf, inverse_samples)
    if with_plot
        plot(cf)
    end
    return CameraCalibrationMeta.Calibration(cf)
end
#
# """
#     c(i::RowCol, extrinsic)
# Convert the row column StaticArray `i` to its real-world equivalent, `XYZ`, for the extrinsic parameters from the `extrinsic` image (given as an index or file name).
# """
# (c::Calibration)(i::RowCol, extrinsic_index::Int) = c.image2real[extrinsic_index](i)
#
# """
#     c(xyz::XYZ, extrinsic)
# Convert the x, y, z StaticArray `xyz` to its pixel-coordinate equivalent, `RowCol` for the extrinsic parameters from the `extrinsic` image (given as an index or file name).
# """
# (c::Calibration)(xyz::XYZ, extrinsic_index::Int) = c.real2image[extrinsic_index](xyz)
#
# function (c::Calibration)(coordinate, extrinsic_file::AbstractString)
#     extrinsic_index = findfirst(==(extrinsic_file), c.files)
#     c(coordinate, extrinsic_index)
# end
#
# """
#     rectification(c, extrinsic_index)
# Return a function that accepts an instance of `::RowCol` and converts it to its real-world equivalent for the extrinsic parameters from the `extrinsic_index` image, without its third dimension (which would be ≈ 0): an `xy` coordinate.
# """
# rectification(c, extrinsic_index) = pop ∘ c.image2real[extrinsic_index]

function _reprojection(c, i)
    imgpoints = c.imgpointss[i]
    reprojected = c.(c.objpoints, i)
    sum(LinearAlgebra.norm_sqr, reprojected .- imgpoints)
end

"""
    calculate_errors(c)
Calculate reprojection, projection, distance, and inverse errors for the calibration `c`. `distance` measures the mean error of the distance between all adjacent checkerboard corners from the expected `checker_size`. `inverse` measures the mean error of applying the calibration's transformation and its inverse `inverse_samples` times.
"""
function calculate_errors(c, inverse_samples=100)
    reprojection = 0.0
    projection = 0.0
    distance = 0.0
    inverse = 0.0
    for (i, imgpoints) in pairs(c.imgpointss)
        reprojection += _reprojection(c, i)

        projected = c.(imgpoints, i)
        projection += sum(LinearAlgebra.norm_sqr, projected .- c.objpoints)

        distance += sum(1:2) do dims
            sum(abs2, norm.(diff(projected; dims)) .- c.checker_size)
        end

        inverse += sum(1:inverse_samples) do _
            rc = rand(RowCol{Float64}) .* (c.sz .- 1) .+ 1
            projected = c(rc, i)
            reprojected = c(projected, i)
            LinearAlgebra.norm_sqr(rc .- reprojected)
        end
    end
    n = prod(c.n_corners)*c.n_files
    reprojection = sqrt(reprojection/n)
    projection = sqrt(projection/n)
    distance = sqrt(distance/prod(c.n_corners .- 1)/c.n_files)
    inverse = sqrt(inverse/inverse_samples/c.n_files)
    return (; n = c.n_files, reprojection, projection, distance, inverse)
end

"""
    improve
Identify all the images that had relatively high reprojection errors, and rerun the calibration without them. Include a maximum of `n` images with the lowest reprojection error, or all the images with an error lower than `threshold`.
"""
function improve(c, n, threshold)
    c.n_files ≤ n && return c
    reprojection = sqrt.(_reprojection.(Ref(c), 1:c.n_files) ./ prod(c.n_corners))
    cutoff = max(threshold, sort(reprojection)[n])
    files = [file for (file, ϵ) in zip(c.files, reprojection) if ϵ ≤ cutoff]
    CalibrationFit(files, c.n_corners, c.checker_size; with_distortion=c.k ≠ 0)
end
improve(c, n::Int, threshold::Nothing) = improve(c, n, 2)
improve(c, n::Nothing, threshold::Int) = improve(c, 15, threshold)
improve(c, n::Nothing, threshold::Nothing) = c
