using Dates, LinearAlgebra, Statistics, Random
using CSV
using DataFramesMeta
using AlgebraOfGraphics
using GLMakie
using CairoMakie
# using StatsBase
using GLM


const pt = 4/3
const inch = 96
const cm = inch / 2.54

set_theme!(Theme(Label = (;fontsize = 8pt), Figure = (size = (5cm, 10cm), fontsize = 8pt, fonts = (; regular = "Helvetica")), Axis = (xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt), Legend = (labelsize = 8pt, labelfont = "Helvetica")))


function get_total_rotation(start, stop, cw, fullturns)
    if start < stop
        if cw 
            stop -= (fullturns + 1)*2π #
        else
            stop += fullturns*2π #
        end
    elseif start > stop
        if cw 
            stop += fullturns*2π #
        else
            stop += (fullturns + 1)*2π #
        end
    else
        throw(ArgumentError("start and stop can't be equal"))
    end
    return stop - start
end

get_total_rotation(start, stop, lastcw) = get_total_rotation(start, stop, start > stop, 0)

function fix_stop_equals_start(start, stop, lastcw)
    start ≠ stop && return stop
    Δ = lastcw ? -0.01 : 0.01
    return stop + Δ
end

function plotit!(ax, df, column)
    w = 3#90/100*h
    for (k, g) in pairs(groupby(sort(df, column, rev = true), :elevation))
        n = nrow(g)
        elevation = k.elevation
        c = g[!, column]
        y = range(0, 100, n + 1)[1:end-1]
        h = step(y)
        for (y, c) in zip(y, c)
            poly!(ax, Rect(elevation - w*0.5, y, w, h), color = c ? :black : :transparent, strokecolor = :lightgray, strokewidth = 1)
        end
    end
end




df = @chain "rotation_elevation_compiled.csv" begin
    CSV.read(DataFrame; missingstring = ["", "NA"])
    @aside begin @chain _ begin
            @subset :condition .== "dark"
            dropmissing(:total_rotation)
            dark = _
        end
    end
    @subset :person_extracted .== "Elin" :condition .== "LED"
    @select Not(:condition, :index, :nr_in_cond, :person_extracted, :year_extracted)
    dropmissing(:go_down_angle_FACE)
    dropmissing(:individual)
    disallowmissing!(Cols(:elevation, :exit_angle_BUTT, :full_lap_1, :go_down_angle_FACE, :individual, :nr_direction_changes, :placed_from_angle_FACE, :rotation_1_direction, :rotation_category, :stop_1_angle_FACE, :total_rotation))
    @aside begin
        @assert all(∈(("cw", "ccw", "both")), _.rotation_category)
        for col in (:rotation_1_direction, :rotation_2_direction, :rotation_3_direction)
            @assert all(∈(("cw", "ccw")), skipmissing(_[!, col]))
        end
        @assert all(==("both"), @subset(_, :nr_direction_changes .> 0).rotation_category)
        @assert all(>(0), @subset(_, :rotation_category .== "both").nr_direction_changes)
        @assert all(ismissing, @subset(_, ismissing.(:stop_2_angle_FACE)).full_lap_2)
        @assert all(ismissing, @subset(_, ismissing.(:full_lap_2)).stop_2_angle_FACE)
        @assert all(ismissing, @subset(_, ismissing.(:stop_3_angle_FACE)).full_lap_3)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_3)).stop_3_angle_FACE)
        @assert all(!ismissing, @subset(_, .!ismissing.(:stop_2_angle_FACE)).full_lap_2)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_2)).stop_2_angle_FACE)
        @assert all(!ismissing, @subset(_, .!ismissing.(:stop_3_angle_FACE)).full_lap_3)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_3)).stop_3_angle_FACE)
    end
    @transform :cw = :rotation_1_direction .== "cw"
    transform([:exit_angle_BUTT, :go_down_angle_FACE, :placed_from_angle_FACE, :stop_1_angle_FACE, :stop_2_angle_FACE, :stop_3_angle_FACE, :total_rotation] .=> ByRow(passmissing(deg2rad)); renamecols = false)
    @transform :go_down_angle_FACE = fix_stop_equals_start.(:placed_from_angle_FACE, :go_down_angle_FACE, :cw)
    @transform :placed2down = get_total_rotation.(:placed_from_angle_FACE, :go_down_angle_FACE, :cw, :full_lap_1)
    @transform :placed_from_left = sign.(:placed_from_angle_FACE .- :go_down_angle_FACE) .> 0
    @transform :shorter_direction = :placed_from_left .== :cw
    @rtransform! :lap = :full_lap_1 > 0 || coalesce(:full_lap_2, 0) > 0 || coalesce(:full_lap_3, 0) > 0
    @transform :changed_direction = :nr_direction_changes .> 0
end


fig = Figure(size = (12cm, 11cm))
row1 = fig[1,1] = GridLayout()
xlimits = (-5, 95)
width = 5
relative = 0.1
ax11 = Axis(row1[1,1], ylabel = "Absolut total rotation (°)", xticks = 0:30:90, yticks = 0:180:2000, limits = (xlimits, nothing))
boxplot!(ax11, df.elevation, rad2deg.(df.total_rotation); width, color = :gray, show_outliers = false)
ax12 = Axis(row1[1,2], yticks = 0:180:2000, limits = (xlimits, nothing), xticks = ([45], ["dark"]))
boxplot!(ax12, 45ones(nrow(dark)), dark.total_rotation; width = width/relative, color = :gray, show_outliers = false)
linkyaxes!(ax11, ax12)
hidexdecorations!(ax12, ticklabels = false, ticks = false)
hideydecorations!(ax12, grid = false, minorgrid = false)
colsize!(row1, 2, Relative(relative))


row2 = fig[2,1] = GridLayout()


limits = (nothing, (nothing, 75))
ax2 = Axis(row2[1,1]; xticks = 0:30:90, ylabel = "Shorter direction (%)", limits)
plotit!(ax2, df, :shorter_direction)



ax3 = Axis(row2[1,2]; xticks = 0:30:90, ylabel = "Direction changes (%)", limits)
plotit!(ax3, df, :changed_direction)

hideydecorations!(ax3, label = false, grid = false, minorgrid = false)

linkaxes!(ax3, ax2)


m = glm(@formula(lap ~ elevation), df, Binomial())

n = 100
newdf = DataFrame(elevation = range(0, 90, n), lap = falses(n))
plu = predict(m, newdf, interval = :confidence)
newdf.prediction = 100disallowmissing(plu.prediction)
newdf.lower = 100disallowmissing(plu.lower)
newdf.upper = 100disallowmissing(plu.upper)

ax4 = Axis(row2[1,3]; xticks = 0:30:90, ylabel = "Extra lap (%)", limits)

band!(ax4, newdf.elevation, newdf.lower, newdf.upper, color = :lightgray)
lines!(newdf.elevation, newdf.prediction, color = :white)

plotit!(ax4, df, :lap)

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
# df2 = combine(groupby(df, [:elevation, :individual]), r"full_lap" => ((xs...) -> Ref(only.(xs))) => [:one, :two, :three])
# df3 = @rsubset df2 !iszero(coalesce(:one, 0)) || !iszero(coalesce(:two, 0)) || !iszero(coalesce(:three, 0))
#
# @show df3
#
# @rtransform! df2 :lap = (coalesce(:one, 0) + coalesce(:two, 0) + coalesce(:three, 0)) > 0
#
# @select! df2 Not(:individual, :one, :two, :three)
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
# # df2 = combine(groupby(df, [:elevation, :individual]), Cols(r"full_lap", r"stop_") => ((f1, f2, f3, s1, s2, s3) -> Ref((; one = only(xs[1]), two = only(xs[2])))) => [:one, :two])
#
# df2 = transform(df, Cols(r"full_lap") .=> ByRow(x -> coalesce(x, 0)), renamecols = false)
# df2 = combine(groupby(df2, [:elevation, :individual]), r"full_lap" => ((xs...) -> Ref((; one = only(xs[1]), two = only(xs[2])))) => [:one, :two])
# @rsubset! df2 !iszero(:one) || !iszero(:two) 
# select!(df2, Not(:individual))
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
# # df2 = combine(groupby(df2, [:elevation, :individual]), r"full_lap" => ((xs...) -> any(>(0), only.(xs))) => :one)
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
# # transform!(groupby(df2, :elevation), :individual => (xs -> 1:length(xs)) => :y)
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
