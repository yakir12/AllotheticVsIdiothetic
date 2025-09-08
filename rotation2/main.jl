using Dates, LinearAlgebra, Statistics, Random
using CSV
using DataFramesMeta
using AlgebraOfGraphics
using GLMakie
using CairoMakie
# using StatsBase
using GLM

cutoff = 10

const pt = 4/3
const inch = 96
const cm = inch / 2.54

set_theme!(Theme(Label = (;fontsize = 8pt), Figure = (size = (5cm, 10cm), fontsize = 8pt, fonts = (; regular = "Helvetica")), Axis = (xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt), Legend = (labelsize = 8pt, labelfont = "Helvetica")))

function angular_distance(θ₁, θ₂)
    Δ = abs(θ₂ - θ₁)
    Δ = mod(Δ, 360)
    return min(Δ, 360 - Δ)
end


function normalize_angle(degrees)
    normalized = degrees % 360
    if normalized < 0
        normalized += 360
    end
    return normalized
end

function get_total_rotation(start, stop, clockwise, turns)
    Δ = stop - start
    if Δ > 180
        Δ -= 360
    elseif Δ < -180
        Δ += 360
    end
    if clockwise
        if Δ < 0
            Δ += 360
        end
        Δ += 360turns
    else
        if Δ > 0
            Δ -= 360
        end
        Δ -= 360turns
    end
    return Δ
end


function get_residual(placed, dance, clockwise)
    down = placed + dance
    ϵ = rem(down, 360)
    if clockwise
        if ϵ < 180
            ϵ
        else
            ϵ - 360
        end
    else
        if ϵ < 180
            ϵ - 180
        else
            360 - ϵ
        end
    end
end

function plotit!(ax, df, column)
    w = 3#90/100*h
    for (k, g) in pairs(groupby(sort(df, column, rev = true), :elevation))
        n = nrow(g)
        elevation = k.elevation
        c = g[!, column]
        y = range(0, 1, n + 1)[1:end-1]
        h = step(y)
        for (y, c) in zip(y, c)
            poly!(ax, Rect(elevation - w*0.5, y, w, h), color = c ? :black : :transparent, strokecolor = :lightgray, strokewidth = 1)
        end
    end
end


# condition               elevation               exit_angle_BUTT
# full_lap_1              full_lap_2              full_lap_3
# go_down_angle_FACE      index                   individual
# nr_direction_changes    nr_in_cond              person_extracted
# placed_from_angle_FACE  rotation_1_direction    rotation_2_direction
# rotation_3_direction    rotation_category       stop_1_angle_FACE
# stop_2_angle_FACE       stop_3_angle_FACE       total_rotation
# year_extracted

df = @chain "rotation_elevation_compiled.csv" begin
    CSV.read(DataFrame; missingstring = ["", "NA"], select = ["elevation",
                                                              "condition",
                                                              "total_rotation",
                                                              "person_extracted",
                                                              "go_down_angle_FACE",
                                                              "individual",
                                                              "exit_angle_BUTT",
                                                              "full_lap_1",
                                                              "full_lap_2",
                                                              "full_lap_3",
                                                              "nr_direction_changes",
                                                              "placed_from_angle_FACE",
                                                              "rotation_1_direction",
                                                              "rotation_2_direction",
                                                              "rotation_3_direction",
                                                              "rotation_category",
                                                              "stop_1_angle_FACE",
                                                              "stop_2_angle_FACE",
                                                              "stop_3_angle_FACE"])
    @rename begin
        :id = $"individual"
        :placed = $"placed_from_angle_FACE"
        # :direction = $"rotation_1_direction"
        :category = $"rotation_category"
        :exit = $"exit_angle_BUTT"
        :down = $"go_down_angle_FACE"
        :abs_total = $"total_rotation"
    end
    # @transform :total = 360 .- :total
    @aside begin @chain _ begin
            @subset :condition .== "dark"
            @select :abs_total
            dropmissing(:abs_total)
            dark = _
        end
    end
    @subset :person_extracted .== "Elin" :condition .== "LED"
    @select Not(:condition, :person_extracted)
    # dropmissing(:id)
    dropmissing(:down)
    disallowmissing!(Cols(:elevation, :exit, :full_lap_1, :down, :id, :nr_direction_changes, :placed, :rotation_1_direction, :category, :stop_1_angle_FACE, :abs_total))
    @aside begin
        @assert all(∈(("cw", "ccw", "both")), _.category)
        for col in (:rotation_1_direction, :rotation_2_direction, :rotation_3_direction)
            @assert all(∈(("cw", "ccw")), skipmissing(_[!, col]))
        end
        @assert all(==("both"), @subset(_, :nr_direction_changes .> 0).category)
        @assert all(>(0), @subset(_, :category .== "both").nr_direction_changes)
        @assert all(ismissing, @subset(_, ismissing.(:stop_2_angle_FACE)).full_lap_2)
        @assert all(ismissing, @subset(_, ismissing.(:full_lap_2)).stop_2_angle_FACE)
        @assert all(ismissing, @subset(_, ismissing.(:stop_3_angle_FACE)).full_lap_3)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_3)).stop_3_angle_FACE)
        @assert all(!ismissing, @subset(_, .!ismissing.(:stop_2_angle_FACE)).full_lap_2)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_2)).stop_2_angle_FACE)
        @assert all(!ismissing, @subset(_, .!ismissing.(:stop_3_angle_FACE)).full_lap_3)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_3)).stop_3_angle_FACE)
    end
    @aside tbl = copy(_)

    @transform begin
        :cw1 = passmissing(==("cw")).(:rotation_1_direction)
        :cw2 = passmissing(==("cw")).(:rotation_2_direction)
        :cw3 = passmissing(==("cw")).(:rotation_3_direction)
    end

    @transform :dance1 = get_total_rotation.(:placed, :stop_1_angle_FACE, :cw1, :full_lap_1)
    @transform :dance2 = passmissing(get_total_rotation).(:stop_1_angle_FACE, :stop_2_angle_FACE, :cw2, :full_lap_2)
    @transform :dance3 = passmissing(get_total_rotation).(:stop_2_angle_FACE, :stop_3_angle_FACE, :cw3, :full_lap_3)

    @aside @assert all(select(select(_, r"dance" => ByRow((xs...) -> sum(abs, skipmissing(xs))) => :abs_total2, :abs_total), All() => ByRow(==) => :same).same)

    @aside @assert all(select(select(_, Cols(r"dance", :placed) => ByRow((xs...) -> normalize_angle(sum(skipmissing(xs)))) => :down2, :down), All() => ByRow(==) => :same).same)
    # @select Not(:down2, 

    @select begin
        :placed = normalize_angle.(:placed .- :down)  # fix this shit
        :exit = normalize_angle.(:exit .- :down .- 180) 
        Not(:down)
    end

    @transform :placed_from_left = :placed .> 180
    @transform :shorter_direction1 = :placed_from_left .== :cw1
    transform(r"full_lap" => ByRow((xs...) -> any(>(0), skipmissing(xs))) => :lap)
    @transform :changed_direction = :nr_direction_changes .> 0
    @transform :residual1 = get_residual.(:placed, :dance1, :cw1)
    @transform :too_close = angular_distance.(:placed, 0) .< cutoff
    # # @transform :residual_magnitude = ifelse.(:direction .== "longer", -:residual, :residual)
end


# scatter(df.exit)






# GLMakie.activate!()

fig = Figure(size = (12cm, 11cm))
row1 = fig[1,1] = GridLayout()
xlimits = (-5, 95)
width = 5
relative = 0.1
ax11 = Axis(row1[1,1], ylabel = "Absolut total rotation (°)", xticks = 0:30:90, yticks = 0:180:2000, limits = (xlimits, nothing))
boxplot!(ax11, df.elevation, df.abs_total; width, color = :gray, show_outliers = false)
ax12 = Axis(row1[1,2], yticks = 0:180:2000, limits = (xlimits, nothing), xticks = ([45], ["dark"]))
boxplot!(ax12, 45ones(nrow(dark)), dark.abs_total; width = width/relative, color = :gray, show_outliers = false)
linkyaxes!(ax11, ax12)
hidexdecorations!(ax12, ticklabels = false, ticks = false)
hideydecorations!(ax12, grid = false, minorgrid = false)
colsize!(row1, 2, Relative(relative))


row2 = fig[2,1] = GridLayout()


limits = (nothing, (nothing, 0.74))
ax2 = Axis(row2[1,1]; xticks = 0:30:90, yticks = 0:0.25:1, ylabel = "Proportion of population", title = "Shorter direction", limits)
plotit!(ax2, df, :shorter_direction1)



ax3 = Axis(row2[1,2]; xticks = 0:30:90, title = "Direction changes", limits)
plotit!(ax3, df, :changed_direction)

hideydecorations!(ax3, label = false, grid = false, minorgrid = false)

linkaxes!(ax3, ax2)

m = glm(@formula(lap ~ elevation), df, Binomial())

n = 100
newdf = DataFrame(elevation = range(0, 90, n), lap = falses(n))
plu = predict(m, newdf, interval = :confidence)
newdf.prediction = disallowmissing(plu.prediction)
newdf.lower = disallowmissing(plu.lower)
newdf.upper = disallowmissing(plu.upper)

ax4 = Axis(row2[1,3]; xticks = 0:30:90, title = "Extra lap", limits)

plotit!(ax4, df, :lap)
band!(ax4, newdf.elevation, newdf.lower, newdf.upper, color = (Makie.wong_colors()[1], 0.2))
lines!(newdf.elevation, newdf.prediction, color = Makie.wong_colors()[1])


hideydecorations!(ax4, label = false, grid = false, minorgrid = false)

linkaxes!(ax4, ax3)


Label(row2[2, :], "Elevation (°)")

Label(row1[1, 1, TopLeft()], "A", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
Label(row2[1, 1, TopLeft()], "B", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
Label(row2[1, 2, TopLeft()], "C", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
Label(row2[1, 3, TopLeft()], "D", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)


GLMakie.activate!()
save("figure.png", fig)

CairoMakie.activate!()
save("figure.pdf", fig)





                                                                                                                                                                   # color = :cw1 => renamer(true => "clockwise", false => "counterclockwise")

# GLMakie.activate!()


_df = @chain df begin
    @subset :lap #.!:changed_direction
    @select :elevation :placed :too_close :changed_direction
end
fig = Figure()
ax = PolarAxis(fig[1,1], theta_0 = -π/2, rticks = sort(unique(_df.elevation)), direction = -1)
h = []
l = []
for (k, g) in pairs(groupby(_df, :changed_direction))
    _h = scatter!(ax, deg2rad.(g.placed), g.elevation, label = k.changed_direction)
    push!(h, _h)
    push!(l, k.changed_direction ? "changed direction" : "didn't")
end
@subset! _df :too_close
c = scatter!(ax, deg2rad.(_df.placed), _df.elevation, markersize = 20, color = :transparent, strokecolor = :red, strokewidth = 1)
Legend(fig[1,2], [h; c], [l; "too close?"])

save("figure2.pdf", fig)

m = glm(@formula(lap ~ elevation), df, Binomial())

m = glm(@formula(lap ~ elevation), @subset(df, .!:too_close), Binomial())

m = glm(@formula(lap ~ elevation), @subset(df, .!:changed_direction), Binomial())

m = glm(@formula(lap ~ elevation), @subset(df, .!:changed_direction, .!:too_close), Binomial())


m = glm(@formula(changed_direction ~ elevation), df, Binomial())

m = glm(@formula(changed_direction ~ elevation), @subset(df, .!:too_close), Binomial())

# using CategoricalArrays
#
# df2 = @transform df :absplaced = abs.(:placed .- 180)
# fmt(from, to, i; leftclosed, rightclosed) = (from + to)/2
# @transform! df2 :placed_bin = cut(:absplaced, range(0, 180, 10), labels = fmt)
# @transform! groupby(df2, [:placed_bin, :elevation]) :nlapped = count(:lap)
#
# data(df2) * mapping(:placed_bin, :nlapped, row = :elevation) * visual(Scatter) |> draw()
#
# m = glm(@formula(lap ~ absplaced*elevation), df2, Binomial())




#
# sdfjhsdjfhlksjhflsdjahf
#
#
#
#
# colors = Makie.wong_colors()
# fig = Figure()
# ax = Axis(fig[1,1], title = "First rotation", ylabel = "Number of individuals", xlabel = "Elevation (°)", xticks = 0:30:90)
# barplot!(ax, df2.elevation, df2.value, stack = df2.variable, color = colors[df2.variable])
# df3 = subset(df2, :variable => ByRow(==(1)))
# text!(ax, df3.elevation, df3.value, text = string.(df3.proportion_shorter, " %"), rotation = π/2, align = (:right, :center), color = :white)
# labels = ["Shorter", "Longer"]
# elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
# title = "Direction"
# Legend(fig[1,2], elements, labels, title)
# save("shorter direction.png", fig)
#
#
# fig = Figure()
# ax = Axis(fig[1,1], title = "First rotation", ylabel = "Rotation (°)", xlabel = "Elevation (°)", xticks = 0:30:90, yticks = 0:180:720)
# boxplot!(ax, df.elevation, abs.(rad2deg.(df.placed2down)); width = 5)
# save("first turning.png", fig)
#
#
# df2 = combine(groupby(df, :elevation), :full_lap_1 => sum => :full1, :full_lap_2 => sum ∘ skipmissing => :full2, :full_lap_3 => sum ∘ skipmissing => :full3, :full_lap_1 => length => :n)
# @transform! df2 :proportion_full = round.(Int, 100 .* (:full1 .+ :full2 .+ :full3) ./ :n)
# df2 = stack(df2, [:full1, :full2, :full3])
# @transform! df2 :variable = parse.(Int, last.(:variable))
#
#
# colors = Makie.wong_colors()
# fig = Figure()
# ax = Axis(fig[1,1], title = "Full lap", ylabel = "Number of laps", xlabel = "Elevation (°)", xticks = 0:30:90)
# barplot!(ax, df2.elevation, df2.value, stack = df2.variable, color = colors[df2.variable])
# df3 = subset(df2, :variable => ByRow(==(1)))
# text!(ax, df3.elevation, df3.value, text = string.(df3.proportion_full, " %"), rotation = π/2, align = (:right, :center), color = :white)
# labels = ["First", "Second"]
# elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
# title = "Lap"
# Legend(fig[1,2], elements, labels, title)
# save("full laps.png", fig)
#
#
#
#
# # df2 = combine(groupby(df, :elevation), :full_lap_1 => sum => :full, :full_lap_1 => length => :n)
# # @transform! df2 :proportion_full_lap = :full ./ :n
# # lines(df2.elevation, df2.proportion_full_lap)
# #
# #
# #
# # df2 = combine(groupby(df, :elevation), :full_lap_1 => sum => :full1, :full_lap_2 => sum ∘ skipmissing => :full2, :full_lap_3 => sum ∘ skipmissing => :full3, :full_lap_1 => length => :n)
# # transform!(df2, r"full" => ByRow((xs...) -> sum(xs)) => :full)
# # @transform! df2 :proportion_full_lap = :full ./ :n
# #
# # lines(df2.elevation, df2.proportion_full_lap)
#
# df2 = combine(groupby(df, :elevation), :nr_direction_changes => sum => :nr, :nr_direction_changes => length => :n)
# @transform! df2 :proportion_nr = :nr ./ :n
#
# fig = Figure()
# ax = Axis(fig[1,1], title = "Direction changes", ylabel = "Number of changes", xlabel = "Elevation (°)", xticks = 0:30:90)
# lines!(ax, df2.elevation, df2.nr)
# save("direction changes.png", fig)
#
#
#
# df2 = combine(groupby(df, [:elevation, :id]), r"full_lap" => ((xs...) -> Ref(only.(xs))) => [:one, :two, :three])
# df3 = @rsubset df2 !iszero(coalesce(:one, 0)) || !iszero(coalesce(:two, 0)) || !iszero(coalesce(:three, 0))
#
# @show df3
#
# @rtransform! df2 :lap = (coalesce(:one, 0) + coalesce(:two, 0) + coalesce(:three, 0)) > 0
#
# @select! df2 Not(:id, :one, :two, :three)
#
# CSV.write("binomial data.csv", df2)
#
# sort!(df2, [:elevation, order(:lap, rev = true)])
#
# df3 = combine(groupby(df2, :elevation), nrow => :n, :lap => count => :lap)
# @transform! df3 :proportion_lap = 100 .* :lap ./ :n
#
#
#
# fig = Figure()
# # w = 3
# gap = 0.0
# ax = Axis(fig[1,1], xticks = 0:30:90, xlabel = "Elevation (°)", ylabel = "Proportion (%)", aspect = 1)
# for (k, g) in pairs(groupby(df2, :elevation))
#     n = nrow(g)
#     elevation = k.elevation
#     c = g.lap
#     y = range(0, 100, n + 1)[1:end-1]
#     h = step(y)
#     w = 3#90/100*h
#     for (y, c) in zip(y, c)
#         poly!(ax, Rect(elevation - w*0.5, y, w, h), color = c ? :black : :white, strokecolor = :lightgray, strokewidth = 2)
#     end
# end
# # lines!(df3.elevation, df3.proportion_lap)
# scatter!(ax, 0, 1, color = :transparent)
# save("full laps2.png", fig)
#
#
#
#
#
#
#
#
#
# # df2 = combine(groupby(df, [:elevation, :id]), Cols(r"full_lap", r"stop_") => ((f1, f2, f3, s1, s2, s3) -> Ref((; one = only(xs[1]), two = only(xs[2])))) => [:one, :two])
#
# df2 = transform(df, Cols(r"full_lap") .=> ByRow(x -> coalesce(x, 0)), renamecols = false)
# df2 = combine(groupby(df2, [:elevation, :id]), r"full_lap" => ((xs...) -> Ref((; one = only(xs[1]), two = only(xs[2])))) => [:one, :two])
# @rsubset! df2 !iszero(:one) || !iszero(:two) 
# select!(df2, Not(:id))
# sort!(df2, [:elevation, :one, :two])
# transform!(groupby(df2, :elevation), :elevation => (xs -> 1:length(xs)) => :y)
#
#
# colors = Makie.wong_colors()
# pushfirst!(colors, colorant"lightgray")
# fig = Figure()
# w = 1
# gap = 0.1
# ax = Axis(fig[1,1], xticks = 0:30:90, yticks = 1:8, xlabel = "Elevation (°)", ylabel = "Individuals")
# for (elevation, c1, c2, y) in eachrow(df2)
#     poly!(ax, Rect(elevation - w, y - 0.5 + gap, w, 1 - gap), color = colors[c1 + 1], strokecolor = :black, strokewidth = 1)
#     poly!(ax, Rect(elevation, y - 0.5 + gap, w, 1 - gap), color = colors[c2 + 1], strokecolor = :black, strokewidth = 1)
# end
# scatter!(ax, 0, 1, color = :transparent)
# # text!(ax, df3.elevation, df3.value, text = string.(df3.proportion_full, " %"), rotation = π/2, align = (:right, :center), color = :white)
# labels = string.(0:2)
# elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
# title = "Full laps"
# Legend(fig[1,2], elements, labels, title)
#
# save("shorter direction2.png", fig)
#
#
#
#
#
#
#
#
# # df2 = transform(df, Cols(r"full_lap") .=> ByRow(x -> coalesce(x, 0)), renamecols = false)
# # df2 = combine(groupby(df2, [:elevation, :id]), r"full_lap" => ((xs...) -> any(>(0), only.(xs))) => :one)
# # # sort!(df2, [:elevation, :one], rev = [false, true])
# # df2 = combine(groupby(df2, :elevation), :one => count => :lap, :one => length => :n)
# # @transform! df2 :proportion_lap = round.(Int, 100 .* :lap ./ :n)
# #
# # colors = Makie.wong_colors()
# # fig = Figure()
# # ax = Axis(fig[1,1], title = "Full lap", ylabel = "Number of laps", xlabel = "Elevation (°)", xticks = 0:30:90)
# # barplot!(ax, df2.elevation, df2.value, stack = df2.variable, color = colors[df2.variable])
# # df3 = subset(df2, :variable => ByRow(==(1)))
# # text!(ax, df3.elevation, df3.value, text = string.(df3.proportion_full, " %"), rotation = π/2, align = (:right, :center), color = :white)
# # labels = ["First", "Second"]
# # elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
# # title = "Lap"
# # Legend(fig[1,2], elements, labels, title)
# # save("full laps.png", fig)
# #
# #
# #
# #
# #
# # transform!(groupby(df2, :elevation), :id => (xs -> 1:length(xs)) => :y)
# #
# # fig = Figure()
# # ax = Axis(fig[1,1])
# # for row in eachrow(df2)
# #     poly!(ax, Rect(row.elevation, row.y, 1, 1), color = :transparent, strokecolor = :black, strokewidth = 1)
# #     laps = filter(!iszero, [row.full_lap_1, row.full_lap_2, row.full_lap_3])
# #     if !isempty(laps)
# #         text!(ax, row.elevation + 0.5, row.y + 0.5, text = join(laps, ", "), align = (:center, :center))
# #     end
# # end
# #
# # @transform! df2 :proportion_full = round.(Int, 100 .* (:full1 .+ :full2 .+ :full3) ./ :n)
# # df2 = stack(df2, [:full1, :full2, :full3])
# # @transform! df2 :variable = parse.(Int, last.(:variable))
