module CameraCalibrationIO

using CameraCalibrationMeta
using Serde, Rotations, CoordinateTransformations, StaticArrays, LinearAlgebra

function Serde.SerJson.ser_type(::Type{CameraCalibrationMeta.AM}, v::SDiagonal{2, Float64})
    return v.diag
end

function Serde.SerJson.ser_type(::Type{CameraCalibrationMeta.LM}, v::SDiagonal{3, Float64})
    return v.diag
end

function Serde.SerJson.ser_type(::Type{CameraCalibrationMeta.AMext}, v::RotationVec{Float64})
    return [v.sx, v.sy, v.sz]
end

function Serde.deser(::Type{CameraCalibrationMeta.AM}, ::Type{SVector{2, Float64}}, v::Vector{Any})
    return SVector{2, Float64}(v)
end

function Serde.deser(::Type{CameraCalibrationMeta.AMext}, ::Type{SVector{2, Float64}}, v::Vector{Any})
    return SVector{2, Float64}(v)
end

function Serde.deser(::Type{CameraCalibrationMeta.AM}, ::Type{SDiagonal{2, Float64}}, v::Vector{Any})
    return Diagonal(SVector{2, Float64}(v))
end

function Serde.deser(::Type{CameraCalibrationMeta.LM}, ::Type{SDiagonal{3, Float64}}, v)
    return Diagonal(SVector{3, Float64}(v))
end

function Serde.deser(::Type{CameraCalibrationMeta.AMext}, ::Type{RotationVec{Float64}}, v::Vector{Any})
    return RotationVec(v...)
end

struct CalibrationIO
    intrinsic::CameraCalibrationMeta.AM
    extrinsics::Vector{CameraCalibrationMeta.AMext}
    scale::CameraCalibrationMeta.LM
    k::Float64
    files::Vector{String}
end

function CameraCalibrationMeta.Calibration(cio::CalibrationIO)
    distort(rc) = lens_distortion(rc, cio.k)
    real2image = .∘(Ref(cio.intrinsic), distort, Ref(PerspectiveMap()), cio.extrinsics, Ref(cio.scale))
    inv_scale, inv_extrinsics, inv_perspective_maps, inv_distort, inv_intrinsic = img2obj(cio.intrinsic, cio.extrinsics, cio.scale, cio.k)
    image2real = .∘(Ref(inv_scale), inv_extrinsics, inv_perspective_maps, inv_distort, Ref(inv_intrinsic))
    CameraCalibrationMeta.Calibration(cio.intrinsic, cio.extrinsics, cio.scale, cio.k, cio.files, real2image, image2real)
end

function CalibrationIO(c::CameraCalibrationMeta.Calibration)
    CalibrationIO(c.intrinsic, c.extrinsics, c.scale, c.k, c.files)
end

function load(file)
    cio = deser_json(CalibrationIO, read(file, String))
    return CameraCalibrationMeta.Calibration(cio)
end

save(file, cio::CalibrationIO) = open(file, "w") do io
    print(io, to_json(cio))
end
save(file, c::CameraCalibrationMeta.Calibration) = save(file, CalibrationIO(c))

end # module CameraCalibrationIO
