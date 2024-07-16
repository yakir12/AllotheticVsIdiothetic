using SimpTrack
using Dates, LinearAlgebra
using CSV, DataFrames, Chain, DataFramesMeta
using Dierckx, VideoIO, OhMyThreads
using StaticArrays, CoordinateTransformations, Rotations
using CairoMakie

const SV = SVector{2, Float64}

include("functions.jl")

data_path = "data"

# Tracking

df = @chain joinpath(data_path, "runs.csv") begin
    CSV.read(DataFrame)
    @transform :calibration_id = :calibration
    @rtransform :name = first(splitext(:file))
    @rtransform :file = joinpath(data_path, :path, :file)
    @rtransform :start = tosecond(:start - Time(0))
    @rtransform :stop = tosecond(:stop - Time(0))
end
df.track = tmap(track, df.file, df.start, df.stop)

# save mini video for debugging purposes
tforeach(save_vid, df.name, df.file, df.track)


# rotate them

function rotate(xy, i)
    trans = Translation(-xy[1])
    p1 = xy[1]
    p2 = xy[i]
    x, y = normalize(p2 - p1)
    θ = π/2 - atan(y, x)
    rot = recenter(LinearMap(Angle2d(θ)), p1)
    cmp = trans ∘ rot
    cmp.(xy)
end

colors = Dict(zip(unique(df.condition), Makie.wong_colors()))

@chain df begin
    @rtransform! :POI_i = findfirst(≥(tosecond(:POI - Time(0))), first(:track))
    @rtransform! :t = first(:track)
    @rtransform! :xy = SV.(last(:track))
    @rtransform! :rotated = rotate(:xy, :POI_i)
    @rtransform! :before = :rotated[1: :POI_i]
    @rtransform! :after = :rotated[:POI_i:end]
    @rtransform! :color = colors[:condition]
end

fig = Figure()
ax = Axis(fig[1,1], aspect=DataAspect())
for row in eachrow(df)
    lines!(ax, row.before, color = :gray)
    lines!(ax, row.after, color = row.color, label = row.condition)
end
axislegend(ax, merge = true)
save("figure.pdf", fig)

