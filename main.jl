using SimpTrack
using Dates, LinearAlgebra
using CSV, DataFrames, Chain, DataFramesMeta
using Dierckx, VideoIO, OhMyThreads
using StaticArrays, CoordinateTransformations, Rotations, AngleBetweenVectors, IterTools
using CairoMakie

const SV = SVector{2, Float64}

include("functions.jl")

data_path = "data"
fps = 25 # frames per second

# Tracking

df = @chain joinpath(data_path, "runs.csv") begin
    CSV.read(DataFrame)
    @rtransform :name = first(splitext(:file))
    @rtransform :file = realpath(joinpath("..", data_path, :path, :file))
    transform!([:start, :POI, :stop] .=> ByRow(tosecond); renamecols = false)
end

df.track = tmap(track, df.file, df.start, df.stop)

# # save mini video for debugging purposes
tforeach(save_vid, df.name, df.file, df.track)


# rotate them

@chain df begin
    transform!([:track, :POI] => ByRow(smooth_rotate_split) => [:before, :after])
    transform!([:before, :after] .=> ByRow(cordlength))
    transform!([:before, :after] .=> ByRow(curvelength))
    transform!([:before_cordlength, :before_curvelength] => ByRow(/) => :before_straightness)
    transform!([:after_cordlength, :after_curvelength] => ByRow(/) => :after_straightness)
    transform!([:before, :after] .=> ByRow(cumulative_angle))
end

fig = Figure();
axs = []
for (i, (k, grp)) in enumerate(pairs(DataFrames.groupby(df, :condition)))
    ax = Axis(fig[1,i], aspect=DataAspect(), title = k.condition)
    push!(axs, ax)
    for row in eachrow(grp)
        lines!(ax, row.before, color = :gray)
        lines!(ax, row.after, color = :black)
    end
end
linkaxes!(axs...)
save("figure.pdf", fig)

CSV.write("summary.csv", select(df, Cols("name", "condition", r"cord", r"curve", r"straight", r"angle")))
