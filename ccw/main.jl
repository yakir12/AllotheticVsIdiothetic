using Statistics#, Random
using CSV
using DataFramesMeta
# using PrettyTables
using AlgebraOfGraphics
using GLMakie
using CairoMakie
using GLM
# using Optim
using HypothesisTests
using MixedModels

const pt = 4/3
const inch = 96
const cm = inch / 2.54

cutoff = 10

set_theme!(Theme(Figure = (size = (5cm, 10cm), fontsize = 8pt, fonts = (; regular = "Helvetica")), Axis = (xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt), Legend = (labelsize = 8pt, labelfont = "Helvetica")))

# rm.(filter(==(".png") ∘ last ∘ splitext, readdir()))

function report(df, file)
    header = """Out of a total of $(nrow(df)) runs, $(count(==("both"), df.category)) turned both directions and $(count(>(0), df.lap)) turned more than 360°. Out of the ones that didn't turn multiple times nor turned more than 360°, $(count(==("shorter"), df.direction1)) turned towards the shorter direction and $(count(==("longer"), df.direction1)) turned towards the longer direction."""

    gs = groupby(df, :direction1)
    df2 = combine(gs, :residual1 .=> [mean, std] .=> [:μ, :σ]) 
    g1, g2 = [DataFrame(g) for g in gs]
    if nrow(g1) < 2 || nrow(g2) < 2
        open("stats $file.md", "w") do io
            print(io, """# Summary
                  $header
                  """
                 )
        end
        return nothing
    end
    vartest = VarianceFTest(g1.residual1, g2.residual1)
    meantest = UnequalVarianceTTest(g1.residual1, g2.residual1)
    meanmagtest = UnequalVarianceTTest(g1.residual_magnitude1, g2.residual_magnitude1)
    df3 = DataFrame(μ = mean(df.residual_magnitude1), σ = std(df.residual_magnitude1))
    t = OneSampleTTest(df.residual_magnitude1)
    μ = mean(df.residual_magnitude1)
    σ = std(df.residual_magnitude1)
    open("stats $file.md", "w") do io
        print(io, """# Summary
              $header

              The residuals—defined as the difference between the expected and actual dance directions—show the following mean and standard deviation for beetles turning in the shorter versus longer direction:

              $df2

              An F-test assessing the null hypothesis that the two groups of residuals have equal variances revealed a significant difference in variance:

              $vartest

              Given this, we applied Welch’s t-test to determine whether the group means differ significantly. The test confirmed a significant difference:

              $meantest

              This indicates that beetles turning in the shorter direction (i.e., placed on the left half-circle and turning clockwise, or placed on the right and turning counterclockwise) performed significantly longer dances than those turning in the longer direction—and vice versa.

              To assess whether the degree of over- or undershooting differed between the dance direction groups, we used Welch’s t-test on the magnitude of the turn (ignoring direction). The result showed no significant difference:

              $meanmagtest

              This suggests that beetles overshot or undershot by approximately the same amount (mean: $(μ)°, standard deviation: $(σ)°). Finally, we tested whether the magnitude of the residuals differed significantly from zero using a one-sample t-test:

              $t

              The result confirmed that the residuals are not centered around zero—indicating that beetles significantly over- or undershot their intended goal direction.""")
    end
    return nothing
end


function angular_distance(θ₁, θ₂)
    Δ = θ₂ - θ₁
    mod(Δ, 180)
end

function wrap_degrees(deg)
    if deg > 180
        deg -= 360
    end
    return deg
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


function get_residual(ϵ, clockwise)
    clockwise || return -ϵ
    return ϵ
end

convert2lap(x) = abs(x) ÷ 360

angular_mean(deg) = atand(mean(sind, deg), mean(cosd, deg))

df = @chain "repeated_dances.csv" begin
    CSV.read(DataFrame; select = ["individual_number",
                                  "placed from angle (degrees)",
                                  "rotation 1 direction",
                                  "rotation 2 direction",
                                  "rotation 3 direction",
                                  "angle before direction change 1 (degrees)",
                                  "angle before direction change 2 (degrees)",
                                  "rotation 1 full lap? (degrees)",
                                  "rotation 2 full lap? (degrees)",
                                  "rotation 3 full lap? (degrees)",
                                  "rotation category measured",
                                  "exit angle (degrees)",
                                  "go down angle (degrees)",
                                  "total absolute degrees of rotation",
                                  "index"])
    @rename begin
        :id = $"individual_number"
        :placed = $"placed from angle (degrees)"
        :category = $"rotation category measured"
        :exit = $"exit angle (degrees)"
        :down = $"go down angle (degrees)"
        :abs_total = $"total absolute degrees of rotation"
        :rotation_1_direction = $"rotation 1 direction"
        :rotation_2_direction = $"rotation 2 direction"
        :rotation_3_direction = $"rotation 3 direction"
        :stop_1_angle = $"angle before direction change 1 (degrees)"
        :stop_2_angle = $"angle before direction change 2 (degrees)"
    end
    @aside @assert all(∈(-720:360:720), skipmissing(unique(Matrix(select(_, r"full lap")))))
    @select begin
        :full_lap_1 = convert2lap.(coalesce.($"rotation 1 full lap? (degrees)", 0))
        :full_lap_2 = ifelse.(ismissing.("rotation_2_direction"), missing, convert2lap.(coalesce.($"rotation 2 full lap? (degrees)", 0)))
        :full_lap_3 = ifelse.(ismissing.("rotation_3_direction"), missing, convert2lap.(coalesce.($"rotation 3 full lap? (degrees)", 0)))
        :stop_1_angle = coalesce.(:stop_1_angle, :down)
        :stop_2_angle = ifelse.(ismissing.("rotation_2_direction"), missing, coalesce.(:stop_2_angle, :down))
        :stop_3_angle = ifelse.(ismissing.("rotation_3_direction"), missing, :down)
        Not("rotation 1 full lap? (degrees)", "rotation 2 full lap? (degrees)", "rotation 3 full lap? (degrees)")
    end
    @transform :full_lap_2 = ifelse.(ismissing.(:stop_2_angle), :full_lap_2, coalesce.(:full_lap_2, 0))
    @transform :full_lap_3 = ifelse.(ismissing.(:stop_3_angle), :full_lap_3, coalesce.(:full_lap_3, 0))
    transform(r"rotation_\d_direction" => ByRow((xs...) -> count(!ismissing, xs) - 1) => :nr_direction_changes)

    dropmissing(:down)
    disallowmissing!(Cols(:exit, :full_lap_1, :down, :id, :nr_direction_changes, :placed, :rotation_1_direction, :category, :stop_1_angle, :abs_total))
    @aside begin
        @assert all(∈(("cw", "ccw", "both")), _.category)
        for col in (:rotation_1_direction, :rotation_2_direction, :rotation_3_direction)
            @assert all(∈(("cw", "ccw")), skipmissing(_[!, col]))
        end
        @assert all(==("both"), @subset(_, :nr_direction_changes .> 0).category)
        @assert all(>(0), @subset(_, :category .== "both").nr_direction_changes)
        @assert all(ismissing, @subset(_, ismissing.(:stop_2_angle)).full_lap_2)
        @assert all(ismissing, @subset(_, ismissing.(:full_lap_2)).stop_2_angle)
        @assert all(ismissing, @subset(_, ismissing.(:stop_3_angle)).full_lap_3)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_3)).stop_3_angle)
        @assert all(!ismissing, @subset(_, .!ismissing.(:stop_2_angle)).full_lap_2)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_2)).stop_2_angle)
        @assert all(!ismissing, @subset(_, .!ismissing.(:stop_3_angle)).full_lap_3)
        @assert all(!ismissing, @subset(_, .!ismissing.(:full_lap_3)).stop_3_angle)
    end

    @transform begin
        :cw1 = passmissing(==("cw")).(:rotation_1_direction)
        :cw2 = passmissing(==("cw")).(:rotation_2_direction)
        :cw3 = passmissing(==("cw")).(:rotation_3_direction)
    end

    @transform :dance1 = get_total_rotation.(:placed, :stop_1_angle, :cw1, :full_lap_1)
    @transform :dance2 = passmissing(get_total_rotation).(:stop_1_angle, :stop_2_angle, :cw2, :full_lap_2)
    @transform :dance3 = passmissing(get_total_rotation).(:stop_2_angle, :stop_3_angle, :cw3, :full_lap_3)

    @aside @assert all(select(select(_, r"dance" => ByRow((xs...) -> sum(abs, skipmissing(xs))) => :abs_total2, :abs_total), All() => ByRow(==) => :same).same)

    @aside @assert all(select(select(_, Cols(r"dance", :placed) => ByRow((xs...) -> normalize_angle(sum(skipmissing(xs)))) => :down2, :down), All() => ByRow(==) => :same).same)

    transform(r"dance" => ByRow((xs...) -> sum(skipmissing(xs))) => :total_dance) 
    transform(r"full_lap_" => ByRow((xs...) -> sum(skipmissing(xs))) => :lap) 

    @groupby :id
    @transform :mean_down = angular_mean(:down) :n = 1:length(:down)
    @transform :placed = wrap_degrees.(:placed .- :mean_down)
    @transform :down = wrap_degrees.(:down .- :mean_down) 
    @rtransform :direction1 = (:placed .< 0) == :cw1 ? "shorter" : "longer" 
    @transform :residual1 = get_residual.(:down, :cw1)
    @transform :residual_magnitude1 = ifelse.(:direction1 .== "longer", -:residual1, :residual1)
    @transform :overshoot = :residual1 .> 0
    @transform :too_close = angular_distance.(:placed, 0) .< cutoff
    @transform :changed_direction = :nr_direction_changes .> 0
    @transform :minimized_rotation = :category .≠ "both" .&& abs.(:dance1) .< 180

end


m = glmm(@formula(minimized_rotation ~ 1 + (1|id)), df, Binomial())

@show m
n = 100
newdf = DataFrame(elevation = range(0, 90, n), minimized_rotation = falses(n), id = 9999)
plu = predict(m, newdf, new_re_levels=:population, type=:response)

newdf.prediction = disallowmissing(plu.prediction)
newdf.lower = disallowmissing(plu.lower)
newdf.upper = disallowmissing(plu.upper)

# using MixedModels
#
# m = lmm(@formula(exit ~ down + (1|id)), df)
#
# corr = @chain df begin
#     @groupby :id
#     transform(:down => angular_mean => :mean_down, 
#             :exit => angular_mean => :mean_exit)
#     @transform :down = rem2pi.(:down .- :mean_down, RoundNearest) :exit = rem2pi.(:exit .- :mean_exit, RoundNearest)
# end

# AlgebraOfGraphics.set_aog_theme!()


# (pregrouped([0], [1]) * visual(ABLines; color = :gray) + data(corr) * mapping(:exit => rad2deg => "Mean exit (°)", :down => rad2deg => "Mean down (°)", color = :id => nonnumeric) * visual(Scatter)) |> draw(; axis = (; xticks = -180:90:180, yticks = -180:90:180, aspect = DataAspect()))

report(df, "nothing excluded")

notboth = @subset df :category .≠ "both"
report(notboth, "both excluded")

reversing = (; xreversed = false, yreversed = false)

example_individual = 24

colors = reverse(Makie.wong_colors())
gap = 40
max_y = maximum(abs, notboth.dance1)
fig = Figure(size = (12cm, 12cm))
ax2 = Axis(fig[1:3,1]; reversing..., limits = (-180 - gap, 180 + gap, -max_y - gap, max_y + gap), aspect = AxisAspect((180 + gap)/(max_y + gap)), xaxisposition = :top, yaxisposition = :right, xticks = ([-90, 90], [rich("Left", color = colors[5]), rich("Right", color = colors[4])]), yticks = (max_y ./ [-2, 2], [rich("Counterclockwise", color = colors[7]), rich("Clockwise", color = colors[6])]), yticklabelrotation = -π/2)
hidespines!(ax2)
hidedecorations!(ax2, ticklabels = false)
ax = Axis(fig[1:3,1]; reversing..., limits = (-180 - gap, 180 + gap, -max_y - gap, max_y + gap), aspect = AxisAspect((180 + gap)/(max_y + gap)), yticks = -720:180:720, xticks = -180:180:180, ylabel = "Total rotation (°)", xlabel = "Initial orientation relative to intended bearing (°)")
for (i, label) in zip([0, 1, -1], ["shorter rotation direction", "longer rotation direction", "longer rotation direction"])
    ablines!(ax, 360i, -1; label, color = abs(i), colorrange = (0, 2), colormap = colors)#, linestyle = :dash)
end
_df = @subset notboth :id .≠ example_individual
scatter!(ax, _df.placed, _df.dance1, color = (:black, 0.5))
_df = @subset notboth :id .== example_individual
scatter!(ax, _df.placed, _df.dance1, color = (:red, 0.5))
# _df = @subset notboth :lap .> 0
# scatter!(ax, _df.placed, _df.dance1, color = :transparent, markersize = 20, strokecolor = :red, strokewidth = 1)
poly!(ax, Rect(-180 - gap/2, -max_y - gap/2, 175 + gap/2, 2*max_y + gap), color = (colors[5], 0.2))
poly!(ax, Rect(5, -max_y - gap/2, 180 + gap/2, 2*max_y + gap), color = (colors[4], 0.2))
poly!(ax, Rect(-180 - 0.75gap, -max_y - 0.75gap, 360 + 1.5gap, max_y - 5 + 0.75gap), color = :transparent, strokecolor = colors[7], strokewidth = 2)
poly!(ax, Rect(-180 - 0.75gap, 5, 360 + 1.5gap, max_y - 5 + 0.75gap), color = :transparent, strokecolor = colors[6], strokewidth = 2)
Legend(fig[1,2], ax, merge = true)
ax = Axis(fig[2,2], xlabel = "Residuals (°)", ylabel = "Counts", xticks = -180:90:180)
hist!(ax, notboth.residual_magnitude1, color = :black)
ax = Axis(fig[3,2], xlabel = "Sequential rotation events", ylabel = "Absolute total rotation (°)", yticks = 0:180:1000)
for g in groupby(df, :id)
    lines!(ax, g.n, g.abs_total, color = (:black, 0.1))
end
boxplot!(ax, df.n, df.abs_total, color = :gray)

Label(fig[1:3, 1,  TopLeft()], "A", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
Label(fig[2, 2, TopLeft()], "B", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
Label(fig[3, 2, TopLeft()], "C", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)

# display(fig)

# GLMakie.activate!()
CairoMakie.activate!()

save("scatter.pdf", fig)



# GLMakie.activate!()
fig = Figure(size = (10cm, 20cm))
xy = data(subset(notboth, :id => ByRow(≠(example_individual)))) * mapping(:placed, :dance1, layout = :id => nonnumeric) * visual(Scatter; color = (:black, 0.5))
xy42 = data(subset(notboth, :id => ByRow(==(example_individual)))) * mapping(:placed, :dance1, layout = :id => nonnumeric) * visual(Scatter; color = (:red, 0.5))
abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = [2,1,3,1,2])) * mapping(:a, :b, color = :color) * visual(ABLines; color = colors)
toplot = abline + xy + xy42
g = draw!(fig[1,1], toplot, scales(; Layout = (; legend = false) ); axis = (; reversing..., xlabel = "Initial orientation relative to intended bearing (°)", ylabel = "Total rotation (°)", xticks = [-180, 180], yticks = -720:360:720, aspect = DataAspect()))
resize_to_layout!(fig)

save("inidividual relationship.pdf", fig)
save("inidividual relationship.png", fig)






# data(df2) * mapping(:residual, color = :direction) * visual(Hist) |> draw()


@subset! notboth :lap .== 0

report(notboth, "both excluded and extra laps excluded")

both = @subset df :category .== "both"

fig = Figure(size = (12cm, 12cm))
ax2 = Axis(fig[1:2,1]; reversing..., limits = (-180 - gap, 180 + gap, -max_y - gap, max_y + gap), aspect = AxisAspect((180 + gap)/(max_y + gap)), xaxisposition = :top, yaxisposition = :right, xticks = ([-90, 90], [rich("Right", color = colors[5]), rich("Left", color = colors[4])]), yticks = (max_y ./ [-2, 2], [rich("Clockwise", color = colors[7]), rich("Counterclockwise", color = colors[6])]), yticklabelrotation = -π/2)
hidespines!(ax2)
hidedecorations!(ax2, ticklabels = false)
ax = Axis(fig[1:2,1]; reversing..., limits = (-180 - gap, 180 + gap, -max_y - gap, max_y + gap), aspect = AxisAspect((180 + gap)/(max_y + gap)), yticks = -720:180:720, xticks = -180:180:180, ylabel = "Total rotation (°)", xlabel = "Initial orientation relative to intended bearing (°)")
for (i, label) in zip([0, 1, -1], ["shorter rotation direction", "longer rotation direction", "longer rotation direction"])
    ablines!(ax, 360i, -1; label, color = abs(i), colorrange = (0, 2), colormap = colors)#, linestyle = :dash)
end
scatter!(ax, both.placed, both.total_dance, color = (:black, 0.5))
poly!(ax, Rect(-180 - gap/2, -max_y - gap/2, 175 + gap/2, 2*max_y + gap), color = (colors[5], 0.2))
poly!(ax, Rect(5, -max_y - gap/2, 180 + gap/2, 2*max_y + gap), color = (colors[4], 0.2))
poly!(ax, Rect(-180 - 0.75gap, -max_y - 0.75gap, 360 + 1.5gap, max_y - 5 + 0.75gap), color = :transparent, strokecolor = colors[7], strokewidth = 2)
poly!(ax, Rect(-180 - 0.75gap, 5, 360 + 1.5gap, max_y - 5 + 0.75gap), color = :transparent, strokecolor = colors[6], strokewidth = 2)
Legend(fig[1,2], ax, merge = true)
ax = Axis(fig[2,2], xlabel = "Residuals (°)", ylabel = "Counts", xticks = -180:90:180)
hist!(ax, both.residual_magnitude1, color = :black)

save("scatter both.pdf", fig)


report(both, "single turning directions excluded")


morelaps = @subset df :lap .> 0 :category .≠ "both"

report(morelaps, "both excluded and laps shorter than 360 excluded")




_df = @chain df begin
    @subset :lap .> 0 #.!:changed_direction
    @select :placed :too_close :changed_direction
end

fig = Figure()
ax = PolarAxis(fig[1,1], theta_0 = -π/2, direction = -1)
h = []
l = []
for (k, g) in pairs(groupby(_df, :changed_direction))
    _h = scatter!(ax, deg2rad.(g.placed), 9ones(nrow(g)), label = k.changed_direction)
    push!(h, _h)
    push!(l, k.changed_direction ? "changed direction" : "didn't")
end
@subset! _df :too_close
c = scatter!(ax, deg2rad.(_df.placed), 9ones(nrow(_df)), markersize = 20, color = :transparent, strokecolor = :red, strokewidth = 1)
Legend(fig[1,2], [h; c], [l; "too close?"])

save("figure2.pdf", fig)

m = glm(@formula(overshoot ~ direction1), @subset(df, .!:changed_direction), Binomial())
newdf = DataFrame(direction1 = ["shorter", "longer"], overshoot = falses(2))
plu = predict(m, newdf, interval = :confidence)
newdf.prediction = 100disallowmissing(plu.prediction)
newdf.lower = 100disallowmissing(plu.lower)
newdf.upper = 100disallowmissing(plu.upper)



# using GLM
#
# using CategoricalArrays
#
# fmt(from, to, i; leftclosed, rightclosed) = (from + to)/2
# @transform! df :absplaced = abs.(:placed)
# @transform! df begin
#     :absplaced_bin = cut(:absplaced, range(0, 180, 10), labels = fmt)
#     :lap_bool = :lap .> 0
#     :both = :category .== "both"
# end
# @transform! groupby(df, [:absplaced_bin, :category]) :lapn = count(:lap_bool)
#
# data(df) * mapping(:absplaced_bin, :lapn, row = :both) * visual(Scatter) |> draw()
#
# data(subset(df, :both)) * mapping(:absplaced_bin, :lapn) * visual(Scatter) |> draw()
#
#
#
# df2 = @transform df :placed = abs.(rad2deg.(:placed))
# fmt(from, to, i; leftclosed, rightclosed) = (from + to)/2
# @select! df2 begin
#     :placed
#     :placed_bin = cut(:placed, range(0, 180, 10), labels = fmt)
#     :lap = :lap .> 0
# end
# df3 = combine(groupby(df2, :placed_bin), nrow,  :lap => count => :lap)
# @transform! df3 :percent = 100*(:lap ./ :nrow)
# m = glm(@formula(lap ~ placed), df2, Binomial())
#
# n = 100
# newdf = DataFrame(placed = range(0, 180, n), lap = falses(n))
# plu = predict(m, newdf, interval = :confidence)
# newdf.prediction = 100disallowmissing(plu.prediction)
# newdf.lower = 100disallowmissing(plu.lower)
# newdf.upper = 100disallowmissing(plu.upper)
# GLMakie.activate!()
# band(newdf.placed, newdf.lower, newdf.upper)
# lines!(newdf.placed, newdf.prediction, color = :white)
# scatter!(unwrap.(df3.placed_bin), df3.percent)
#
#
#
# df2 = @transform df :both = :category .== "both"
#
# m = glm(@formula(both ~ placed), df2, Binomial())
#
#
# m = glm(@formula(abs(residual_magnitude) ~ placed), df, Gamma())
#
#
#
# df2 = transform(df, [:placed, :residual_magnitude] .=> ByRow(rad2deg ∘ abs), renamecols = false)
#
# m = glm(@formula(residual_magnitude ~ placed), df2, Gamma())
# m = lm(@formula(log(residual_magnitude) ~ placed), df2)
#
# n = 100
# newdf = DataFrame(placed = range(0, 180, n), lap = falses(n))
# plu = predict(m, newdf, interval = :confidence)
# newdf.prediction = disallowmissing(plu.prediction)
# newdf.lower = disallowmissing(plu.lower)
# newdf.upper = disallowmissing(plu.upper)
# GLMakie.activate!()
# band(newdf.placed, exp.(newdf.lower), exp.(newdf.upper))
# lines!(newdf.placed, exp.(newdf.prediction), color = :white)
# scatter!(df2.placed, df2.residual_magnitude)
#

#
# dfgjhsdflkhjgdlfghjsdflghjdlghj
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# both = @chain both begin
#     @rsubset :category ∉ ("cw", "ccw")
#     @aside @assert all(_.direction .≠ _.category)
#     @select Not(:category)
#     @aside @assert all(c -> all(0 .≤ _[!,c] .≤ 360), [:placed, :exit, :down])
#     transform([:placed, :exit, :down] .=> ByRow(deg2rad), renamecols = false)
#     @select :cw = :direction .== "cw" Not(:direction)
#     @aside @assert all(_.placed .≠ _.down)
#     @select :dance = get_rotation.(:placed, :down, :cw, :lap) Not(:lap)
#     @aside @assert all((sign.(_.dance) .< 0) .== _.cw)
#     @aside @assert all(isapprox.(rem2pi.(_.placed .+ _.dance .- _.down, RoundNearest), 0, atol = 1e-10))
#     @transform :exit = fix_stop_equals_start.(:down, :exit, :cw)
#     @aside @assert all(_.exit .≠ _.down)
#     @select Not(:exit)
#     @groupby :id
#     @select :mean_down = angular_mean(:down) :n = 1:length(:down) Not(:down)
#     @select :placed = rem2pi.(:placed .- :mean_down, RoundNearest) Not(:mean_down)
#     @rtransform :direction = (:placed .> 0) == :cw ? "shorter" : "longer" 
#     @transform :residual = get_residual.(:placed, :dance)
# end
#
# colors = reverse(Makie.wong_colors())
# gap = 40
# fig = Figure(size = (12cm, 10cm))
# ax2 = Axis(fig[1:3,1], limits = (-180 - gap, 180 + gap, -720 - gap, 720 + gap), aspect = AxisAspect((180 + gap)/(720 + gap)), xaxisposition = :top, yaxisposition = :right, xticks = ([-90, 90], [rich("Left", color = colors[5]), rich("Right", color = colors[4])]), yticks = ([-360, 360], [rich("Counterclockwise", color = colors[7]), rich("Clockwise", color = colors[6])]), yticklabelrotation = -π/2)
# hidespines!(ax2)
# hidedecorations!(ax2, ticklabels = false)
# ax = Axis(fig[1:3,1], xreversed = true, yreversed = true, limits = (-180 - gap, 180 + gap, -720 - gap, 720 + gap), aspect = AxisAspect((180 + gap)/(720 + gap)), yticks = -720:180:720, xticks = -180:180:180, ylabel = "Accumulated yaw (°)", xlabel = "Initial orientation relative to intended bearing (°)")
# for (i, label) in zip([0, 1, 2, -1, -2], ["shorter rotation direction", "longer rotation direction", "additional lap", "longer rotation direction", "additional lap"])
#     ablines!(ax, 360i, -1; label, color = abs(i), colorrange = (0, 2), colormap = colors)#, linestyle = :dash)
# end
# scatter!(ax, rad2deg.(both.placed), rad2deg.(both.dance), color = (:black, 0.5))
# poly!(ax, Rect(-180 - gap/2, -720 - gap/2, 175 + gap/2, 2*720 + gap), color = (colors[4], 0.2))
# poly!(ax, Rect(5, -720 - gap/2, 180 + gap/2, 2*720 + gap), color = (colors[5], 0.2))
# poly!(ax, Rect(-180 - 0.75gap, -720 - 0.75gap, 360 + 1.5gap, 715 + 0.75gap), color = :transparent, strokecolor = colors[6], strokewidth = 2)
# poly!(ax, Rect(-180 - 0.75gap, 5, 360 + 1.5gap, 720 + 0.75gap), color = :transparent, strokecolor = colors[7], strokewidth = 2)
# Legend(fig[1,2], ax, merge = true)
# ax = Axis(fig[2,2], xlabel = "Residuals (°)", ylabel = "Counts", xticks = -180:90:180)
# hist!(ax, rad2deg.(both.residual), color = :black)
#
# Label(fig[1:3, 1,  TopLeft()], "A", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
# Label(fig[2, 2, TopLeft()], "B", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
#
#
#
#
# sdjkfhsdjfhasdklfhjl
#
#
# OneSampleTTest(df.residual)
#
#
# df2 = @chain df begin
#     @subset deg2rad(20) .< abs.(:placed) .< deg2rad(160) abs.(:dance) .< deg2rad(400)
#     @transform :residual = get_residual.(:placed, :dance) :shortdirection = :cw .== (:placed .> 0)
# end
#
# data(df2) * mapping(:residual => rad2deg, color = :shortdirection, stack=:shortdirection) * histogram(bins = 20) |> draw()
#
# combine(groupby(df2, :shortdirection), :residual .=> [rad2deg ∘ mean, rad2deg ∘ std])
#
#
#
# sdjkfhsdjfhasdklfhjl
#
# function to_minimize(placed, dance, μ)
#     s = 0.0
#     for (_x, y) in zip(placed, dance)
#         x = rem2pi(_x .- μ, RoundNearest)
#         r, _ = findmin(i -> abs(y - (i*2π - x)), -2:2)
#         s += r
#     end
#     return s
# end
#
# function get_min_residual(placed, dance, cw)
#     intercepts = 360(-2:2)
#     dance2 = -placed .+ deg2rad.(intercepts)
#     Δ = dance .- dance2
#     _,i = findmin(abs.(Δ))
#     down = placed + dance
#     @assert isapprox(Δ[i], rem2pi(down, RoundNearest), atol = 1e-10) Δ[i] , rem2pi(placed + dance, RoundNearest)
#     shoot = if cw == Δ[i] > 0
#         down > deg2rad(intercepts[i]) ? "undershoot" : "overshoot"
#     else
#         down < deg2rad(intercepts[i]) ? "undershoot" : "overshoot"
#     end
#     (; residual = Δ[i], line = intercepts[i], shoot)
# end
#
# df = @chain "data.csv" begin
#     CSV.read(DataFrame; select = ["individual_number",
#                                   "placed from angle (degrees)",
#                                   "rotation 1 direction",
#                                   "rotation category measured",
#                                   "exit angle (degrees)",
#                                   "go down angle (degrees)",
#                                   "full lap"])
#     @rename begin
#         :id = $"individual_number"
#         :placed = $"placed from angle (degrees)"
#         :direction = $"rotation 1 direction"
#         :category = $"rotation category measured"
#         :exit = $"exit angle (degrees)"
#         :down = $"go down angle (degrees)"
#         :lap = $"full lap"
#     end
#     @rsubset :category ∈ ("cw", "ccw")
#     @aside @assert all(_.direction .== _.category)
#     @select Not(:category)
#     @aside @assert all(c -> all(0 .≤ _[!,c] .≤ 360), [:placed, :exit, :down])
#     # transform([:placed, :exit, :down] .=> ByRow(x -> x - 180), renamecols = false)
#     transform([:placed, :exit, :down] .=> ByRow(deg2rad), renamecols = false)
#     @select :cw = :direction .== "cw" Not(:direction)
#     @aside @assert all(_.placed .≠ _.down)
#     @select :dance = get_rotation.(:placed, :down, :cw, :lap) Not(:lap)
#     @aside @assert all((sign.(_.dance) .< 0) .== _.cw)
#     @aside @assert all(isapprox.(rem2pi.(_.placed .+ _.dance .- _.down, RoundNearest), 0, atol = 1e-10))
#     @transform :exit = fix_stop_equals_start.(:down, :exit, :cw)
#     @aside @assert all(_.exit .≠ _.down)
#     @transform :down2exit = get_rotation.(:down, :exit, :cw)
#     @groupby :id
#     select(:down => angular_mean => :mean_down, 
#            :exit => angular_mean => :mean_exit, 
#            [:placed, :dance] => ((placed, dance) -> optimize(μ -> to_minimize(placed, dance, μ), -π, π).minimizer) => :mean_optimal,
#            Not(:down, :exit))
#     # @select :means = [angular_mean(:down), angular_mean(:exit)] Not(:down)
#     @select begin
#         :placed_down = rem2pi.(:placed .- :mean_down, RoundNearest) 
#         :placed_exit = rem2pi.(:placed .- :mean_exit, RoundNearest) 
#         :placed_optimal = rem2pi.(:placed .- :mean_optimal, RoundNearest) 
#         Not(:mean_down, :mean_exit, :mean_optimal)
#     end
#     @rtransform begin 
#         :direction_down = (:placed_down .> 0) == :cw ? "shorter" : "longer" 
#         :direction_exit = (:placed_exit .> 0) == :cw ? "shorter" : "longer"
#         :direction_optimal = (:placed_optimal .> 0) == :cw ? "shorter" : "longer"
#     end
#     transform([:placed_down, :dance, :cw] => ByRow(get_min_residual) => [:residual_down, :line_down, :shoot_down], 
#               [:placed_exit, :dance, :cw] => ByRow(get_min_residual) => [:residual_exit, :line_exit, :shoot_exit],
#               [:placed_optimal, :dance, :cw] => ByRow(get_min_residual) => [:residual_optimal, :line_optimal, :shoot_optimal])
#     # @transform :norm_resid = abs.(:residual_down ./ :dance)
# end
#
#
# for μ in (:_down, :_exit, :_optimal)
#     scatter_layer = data(df) * mapping(Symbol(:placed, μ) => rad2deg, :dance => rad2deg) * visual(Scatter)
#     abline_layer = data(DataFrame(a = 360(-2:2), b = -1, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
#     fig = draw(scatter_layer + abline_layer; axis = (; aspect = DataAspect(), yticks = -720:360:720, xticks = -180:180:180))
#     save("scatter$μ.png", fig)
# end
#
# μ = :_down
# df2 = @transform df :norm_dance = :dance .- (-:placed_down .+ deg2rad.(:line_down))
# @rtransform! df2 :norm_dance = :placed_down > 0 ? -:norm_dance : :norm_dance
# @rtransform! df2 :shorter = :cw == (:placed_down > 0)
# scatter_layer = data(df2) * mapping(:shorter, :norm_dance => rad2deg, color = :shorter) * visual(BoxPlot)
# # abline_layer = data(DataFrame(a = 360(-2:2), b = 0, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
# fig = draw(scatter_layer)#; axis = (; yticks = -720:360:720, xticks = -180:180:180))
#
#
# df2 = @rsubset df (:direction_down == "longer" && :line_down == 360) || (:direction_down == "shorter" && :line_down == 0)
# violin_layer = data(df2) * mapping(:direction_down => sorter("shorter", "longer") => "Turning direction", :residual_down => rad2deg => "Residuals (°)") * visual(BoxPlot)
# # violin_layer = data(df2) * mapping(:direction_down => sorter("shorter", "longer") => "Turning direction", :residual_down => rad2deg => "Residuals (°)", color = :shoot_down => sorter("undershoot", "overshoot") => "", dodge = :shoot_down) * visual(BoxPlot)
# fig = draw(violin_layer)#; axis = (; yticks = 0:30:180, width = 200, height = 200))
# resize_to_layout!(fig)
#
#
# for μ in (:_down, :_exit, :_optimal)
#
#     μ = :_down
#     scatter_layer = data(df) * mapping(Symbol(:placed, μ) => rad2deg, Symbol(:residual, μ) => rad2deg, col = Symbol(:direction, μ), row = Symbol(:line, μ) => nonnumeric, color = Symbol(:shoot, μ)) * visual(Scatter)
#     # band_layer = data(rename(intervals, Dict(x => Symbol(x, μ) for x in [:direction, :line, :shoot]))) * mapping(:placed, :lower, :upper, col = Symbol(:direction, μ), row = Symbol(:line, μ) => nonnumeric, color = Symbol(:shoot, μ)) * visual(Band; alpha = 0.2)
#     fig = draw(scatter_layer)#; axis = (; xticks = 0:90:180, yticks = 0:90:180, limits = (-10, 190, -10, 190), width = 200, height = 200))
#     resize_to_layout!(fig)
#     save("residual$μ.png", fig)
#
# end
#
# sdhgfksdhjgksdjhfgksdfjhg
#
# function get_intervals(direction, line, shoot, placed)
#     if direction == "longer"
#         if line == 0
#             if shoot == "undershoot"
#                 (missing, missing)
#             else
#                 (placed, 180)
#             end
#         else
#             (0, 180)
#         end
#     else
#         if line == 0
#             if shoot == "undershoot"
#                 (0, placed)
#             else
#                 (0, 180)
#             end
#         else
#             (0, 180)
#         end
#     end
# end
#
# intervals = @chain DataFrame(direction = ["longer", "shorter"]) begin
#     @rtransform :line = 0:360:720
#     flatten(:line)
#     @rtransform :shoot = ["undershoot", "overshoot"]
#     flatten(:shoot)
#     @rtransform :placed = 0:180
#     flatten(:placed)
#     transform([:direction, :line, :shoot, :placed] => ByRow(get_intervals) => [:lower, :upper])
#     dropmissing(:lower)
# end
#
#
#
# for mean in (:_down, :_exit, :_optimal)
#     scatter_layer = data(df) * mapping(Symbol(:placed, mean) => rad2deg, :dance => rad2deg) * visual(Scatter)
#     abline_layer = data(DataFrame(a = 360(-2:2), b = -1, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
#     fig = draw(scatter_layer + abline_layer; axis = (; aspect = DataAspect(), yticks = -720:360:720, xticks = -180:180:180))
#     save("scatter$mean.png", fig)
#     scatter_layer = data(df) * mapping(Symbol(:placed, mean) => rad2deg ∘ abs, Symbol(:residual, mean) => rad2deg, col = Symbol(:direction, mean), row = Symbol(:line, mean) => nonnumeric, color = Symbol(:shoot, mean)) * visual(Scatter)
#     band_layer = data(rename(intervals, Dict(x => Symbol(x, mean) for x in [:direction, :line, :shoot]))) * mapping(:placed, :lower, :upper, col = Symbol(:direction, mean), row = Symbol(:line, mean) => nonnumeric, color = Symbol(:shoot, mean)) * visual(Band; alpha = 0.2)
#     fig = draw(band_layer + scatter_layer; axis = (; xticks = 0:90:180, yticks = 0:90:180, limits = (-10, 190, -10, 190), width = 200, height = 200))
#     resize_to_layout!(fig)
#     save("residual$mean.png", fig)
# end
#
#
# # μ = :_down
# # scatter_layer = data(df) * mapping(:dance => rad2deg ∘ abs, :norm_resid, col = Symbol(:direction, μ), row = Symbol(:line, μ) => nonnumeric, color = Symbol(:shoot, μ)) * visual(Scatter)
# # # abline_layer = pregrouped([0], [1]) * visual(ABLines; color = :gray)
# # fig = draw(scatter_layer; axis = (; width = 200, height = 200))
# # resize_to_layout!(fig)
#
# # fig = scatter(rad2deg.(abs.(df.placed_down)), rad2deg.(abs.(df.dance)), axis = (; xlabel = "placed", ylabel = "danced", aspect = DataAspect()))
# # ablines!(0, 1, color = :black)
#
#
#
# skdjfhsljkhflsdhjdfsdhjk
#
#
# df = CSV.read("data.csv", DataFrame, select = ["individual_number",
#                                                "placed from angle (degrees)",
#                                                "rotation 1 direction",
#                                                "rotation category measured",
#                                                "exit angle (degrees)",
#                                                "go down angle (degrees)",
#                                                "full lap",
#                                                "total absolute degrees of rotation"])
#
# transform!(groupby(df, :individual_number), :individual_number => (g -> 1:length(g)) => :n)
#
# dforg = copy(df)
#
# # transform!(groupby(df, :individual_number), :individual_number => (x -> 1:length(x)) => :n)
#
# subset!(df, "rotation category measured" => ByRow(∈(("cw", "ccw"))))
# @assert all(df[!, "rotation 1 direction"] .== df[!, "rotation category measured"])
# select!(df, Not("rotation category measured"))
#
#
# transform!(df, ["placed from angle (degrees)",
#                 "exit angle (degrees)",
#                 "go down angle (degrees)"] .=> ByRow(deg2rad) .=> 
#            ["placed from angle",
#             "exit angle",
#             "go down angle"])
# select!(df, Not("placed from angle (degrees)",
#                 "exit angle (degrees)",
#                 "go down angle (degrees)"))
#
# transform!(df, "rotation 1 direction" => ByRow(==("cw")) => :cw)
# select!(df, Not("rotation 1 direction"))
#
# transform!(df, ["placed from angle", "go down angle", "cw", "full lap"] => ByRow(get_total_rotation) => :total_dance)
# rename!(df, "placed from angle" => :start)
#
# transform!(df, ["go down angle", "exit angle", "cw"] => ByRow(fix_stop_equals_start) => "exit angle")
# transform!(df, ["go down angle", "exit angle", "cw"] => ByRow(get_total_rotation) => :down2exit)
# # rename!(df, "go down angle" => :down)
# rename!(df, "exit angle" => :exit)
#
# select!(df, Not("go down angle", "full lap"))
#
# function to_minimize(start, total_dance, μ)
#     s = 0.0
#     for (_x, y) in zip(start, total_dance)
#         x = rem2pi(_x .- μ, RoundNearest)
#         r, _ = findmin(i -> abs(y - (i*2π - x)), -2:2)
#         s += r
#     end
#     return s
# end
#
# angular_mean(θs) = angle(sum(exp, 1im*θs))
#
# transform!(groupby(df, :individual_number), :exit => angular_mean => :mean_exit, [:start, :total_dance] => ((s, p) -> angular_mean(s .+ p)) => :mean_down, [:start, :total_dance] => ((start, total_dance) -> optimize(μ -> to_minimize(start, total_dance, μ), -pi, pi).minimizer) => :optimal)
#
# fig1 = data(df) * (mapping(:mean_exit => rad2deg => "Mean exit", :mean_down => rad2deg) * visual(Scatter; label = "Mean down") + mapping(:mean_exit => rad2deg => "Mean exit", :optimal => rad2deg) * visual(Scatter; color = :red, label = "Mean optimal") + pregrouped([0], [1]) * visual(ABLines; color = :gray, label = "y = x")) |> draw(; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -180:180:180, ytickformat = "{:n}°", aspect = DataAspect(), width = 200, limits = ((-190, 190),(-190, 190))))
# resize_to_layout!(fig1)
#
# save("means.pdf", fig1)
#
#
# function get_abs_residuals(shortdirection, placed_from_left, start, total_dance, cw)
#     # Δ = asin(sin(total_dance + start))
#     # placed_from_left ? Δ : -Δ
#     intercepts = -720:360:720
#     Δs = deg2rad.(intercepts) .- start .- total_dance
#     _, i = findmin(abs, Δs)
#     Δ = Δs[i]
#     return (intercept = intercepts[i], residual = cw ? Δ : -Δ)
#     # n = abs(total_dance ÷ 2π)
#     # intercept = shortdirection ? 0 : placed_from_left ? n*2π : -n*2π
#     # ŷ = intercept - start
#     # Δ = ŷ - total_dance 
#     # placed_from_left ? Δ : -Δ
# end
#
# example_individual = 24
#
# for μ in (:mean_exit, :mean_down, :optimal)
#     df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)
#     transform!(df2, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
#     df3 = df2#subset(df2, :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
#     # df3 = subset(df2, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
#     transform!(df3, [:placed_from_left, :cw] => ByRow(==) => :shortdirection)
#     transform!(df3, [:shortdirection, :placed_from_left, :start, :total_dance, :cw] => ByRow(get_abs_residuals) => [:intercept, :residual], :start => (s -> abs.(s)) => :x)
#
#     (data(combine(groupby(df3, :shortdirection), :residual => mean => :μ)) * mapping(:μ => rad2deg => "Residuals (°)", col = :shortdirection => renamer(true => "short direction", false => "long direction")) * visual(HLines)) + (data(df3) * mapping(:start => rad2deg ∘ abs => "Initial body orientation (°)", :residual => rad2deg => "Residuals (°)", col = :shortdirection => renamer(true => "short direction", false => "long direction"), color = :total_dance => nonnumeric ∘ (x -> round(Int, abs(x)÷2π)) => "Additional lap") * visual(Scatter)) |> draw(axis = (aspect = DataAspect(), xticks = 0:30:180, yticks = -180:30:180, limits = ((0, 180), (-180, 180)))) |> save("scatter $μ.png")
#     # data(df3) * mapping(:start => abs ∘ rad2deg, :residual => rad2deg, row = :shortdirection => renamer(true => "short direction", false => "long direction")) * (linear() + visual(Scatter)) |> draw() |> save("scatter $μ.png")
#     data(df3) * mapping(:residual => rad2deg, color = :shortdirection => renamer(true => "short direction", false => "long direction")) * histogram(Stairs; bins = 10) |> draw() |> save("hist $μ.png")
#     @show combine(groupby(df3, :shortdirection), :residual => rad2deg ∘ mean)
#
#     # data(df3) * mapping(:shortdirection => renamer(true => "short direction", false => "long direction"), :residual => rad2deg) * visual(Violin; datalimits = extrema) |> draw() |> save("violin $μ.png")
#
#     fig = Figure(size = (12cm, 10cm), fontsize = 8pt, fonts = (; regular = "Helvetica"))
#     subgl_left = GridLayout()
#     subgl_right = GridLayout()
#     fig.layout[1, 1] = subgl_left
#     fig.layout[1, 2] = subgl_right
#     xy = data(subset(df2, :individual_number => ByRow(≠(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation") * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5)
#     xy42 = data(subset(df2, :individual_number => ByRow(==(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation") * visual(Scatter; color = :red, strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5)
#     abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
#     vline = data(DataFrame(m = [-160, 20], M = [-20, 160])) * mapping(:m, :M) * visual(VSpan; color = (:gray, 0.2))
#     # vline = data(DataFrame(geometry = [Rect(-160, -360, 140, 720), Rect(20, -360, 140, 720)])) * mapping(:geometry) * visual(Poly; color = (:gray, 0.2), label = "included")
#     toplot = vline + abline + xy + xy42
#     g = draw!(subgl_left[1,1], toplot; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xlabel = "Initial body orientation", xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect()))
#
#     legend!(subgl_right[1,1], g; framevisible = false)
#     draw!(subgl_right[2,1], data(df3) * mapping(:residual => rad2deg => "Residuals") * visual(Hist; bins = 15); axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, ylabel = "#", xtickformat = "{:n}°"))
#     dance_number = data(dforg) * (mapping(:n => "Dance number", "total absolute degrees of rotation" => "Total absolute rotation") * visual(BoxPlot) + mapping(:n => "Dance number", "total absolute degrees of rotation" => "Total absolute rotation", group = :individual_number => nonnumeric) * visual(Lines; alpha = 0.1, linewidth = 1))
#     draw!(subgl_right[3,1], dance_number; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, yticks = 0:360:1000, ytickformat = "{:n}°"))
#
#     Label(subgl_left[1, 1,  TopLeft()], "A", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
#     Label(subgl_right[2, 1, TopLeft()], "B", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
#     Label(subgl_right[3, 1, TopLeft()], "C", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
#
#     resize_to_layout!(fig)
#     save("relationship $μ.pdf", fig)
#     save("relationship $μ.png", fig)
#     resid = rad2deg.(df3.residual)
#     @show μ, round(Int, mean(resid)), round(Int, std(resid))
# end
#
#
#
#
# μ = :mean_down
# df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)
# transform!(df2, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
# df3 = subset(df2, :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
# # df3 = subset(df2, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
# transform!(df3, [:placed_from_left, :cw] => ByRow(==) => :overshoot)
# transform!(df3, [:overshoot, :placed_from_left, :start, :total_dance] => ByRow(get_abs_residuals) => :residual, :start => (s -> abs.(s)) => :x)
#
# fig = Figure(size = (8cm, 18cm), fontsize = 12pt, fonts = (; regular = "Helvetica"))
# xy = data(subset(df2, :individual_number => ByRow(≠(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation", layout = :individual_number => nonnumeric) * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
# xy42 = data(subset(df2, :individual_number => ByRow(==(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation", layout = :individual_number => nonnumeric) * visual(Scatter; color = :red, strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
# abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
# toplot = abline + xy + xy42
# g = draw!(fig[1,1], toplot; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xlabel = "Initial body orientation", xticks = [-180, 180], xtickformat = "{:n}°", yticks = -720:360:720, ytickformat = "{:n}°", aspect = DataAspect()))
# resize_to_layout!(fig)
# save("inidividual relationship $μ.pdf", fig)
# save("inidividual relationship $μ.png", fig)
#
# # mean_exit, 25, 54
# # mean_down, 26, 53
# # optimal, 27, 57
#
# shuffle!(df2)
# sort!(df2, [:placed_from_left, :cw, :start])
#
# fig = Figure()
# for (i, (k, g)) in enumerate(pairs(groupby(df2, :individual_number)))
#     ij = CartesianIndices((5, 6))[i]
#     ax = Axis(fig[Tuple(ij)...], aspect = DataAspect(), height = 200, width = 200)
#     for (j, row) in enumerate(eachrow(g))
#         radius = j + 1
#         poly!(ax, Makie.GeometryBasics.Polygon(Circle(zero(Point2f), radius+0.5), [Circle(zero(Point2f), radius-0.5)]), color = row.placed_from_left ? :blue : :green, alpha = 0.25)
#         color = (; color = row.cw ? :blue : :green)
#         ts = row.start .+ range(0, row.total_dance, length = 100)
#         lines!(ax, radius*Point2f.(reverse.(sincos.(ts .+ π/2))); color...)
#         α = row.start + row.total_dance
#         α += row.cw ? -π/2 : π/2
#         scatter!(ax, radius*Point2f(reverse(sincos(row.start + π/2))); color..., marker = '|', markersize=10, rotation = row.start + pi/2 + π/2)
#         scatter!(ax, radius*Point2f(reverse(sincos(row.start + row.total_dance + π/2))); color..., marker = :utriangle, markersize=10, rotation=α)
#     end
#     text!(ax, 0, 0; text = string(k.individual_number), align = (:center, :center))
#     hidedecorations!(ax)
#     hidespines!(ax)
# end
# resize_to_layout!(fig)
# save("circles.pdf", fig)
