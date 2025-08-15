using Dates, LinearAlgebra, Statistics, Random
using CSV, DataFrames
using DataFramesMeta, Chain
using AlgebraOfGraphics, CairoMakie
# using StatsBase
using LsqFit

pt = 4/3
inch = 96
cm = inch / 2.54

# function get_iqr(data)
#     iqr_value = iqr(data)
#     q1 = quantile(data, 0.25)  # 25th percentile
#     q3 = quantile(data, 0.75)  # 75th percentile
#     lower_bound = q1 - 1.5 * iqr_value
#     upper_bound = q3 + 1.5 * iqr_value
#     return lower_bound, upper_bound
# end


df = @chain "rotation_all_data.csv" begin
    CSV.read(DataFrame; missingstring="NA")
    @subset :person_extracted .== "Elin" :year_extracted .== 2025 :condition .== "LED"
    dropmissing(:total_rotation)
    # @subset :total_rotation .< 1000
end


# data(df) * mapping(:elevation, :total_rotation) * visual(BoxPlot; width = 5) |> draw(; axis = (; xticks = 0:10:90))

# @rsubset! df rng[1] .< :total_rotation .< rng[2]

get_α(ϵ, C, θ) = C + 2atand(tand(ϵ/2)/cosd(θ))

# fig = Figure()
# ax = Axis(fig[1,1])
# scatter!(ax, df.elevation, df.total_rotation)
# sg = SliderGrid(
#     fig[1, 2],
#     (label = "ϵ", range = range(0, 90, 1000), startvalue = 1),
#     (label = "Factor", range = range(0, 300, 1000), startvalue = 1),
#     tellheight = false)
# θ = 0:90
# l = lift(sg.sliders[1].value, sg.sliders[2].value) do ϵ, C
#     get_α.(ϵ, C, θ)
# end
# lines!(ax, θ, l)


df2 = @chain df begin
    @groupby :elevation
    @combine :total_rotation = mean(:total_rotation)
end

μ = mean(df2.total_rotation)
resid = df2.total_rotation .- μ
SSE1 = round(Int, sum(abs2, resid))
RMSE1 = round(Int, sqrt(mean(abs2, resid)))

@. model(x, p) = get_α.(p[1], p[2], x)
xdata = Float64.(df2.elevation)
ydata = Float64.(df2.total_rotation)
lower = Float64[0, 0]
upper = Float64[90, 400]
p0 = Float64[4, 150]
fit = curve_fit(model, xdata, ydata, p0; lower, upper)
# for i in 1:2
#     set_close_to!(sg.sliders[i], fit.param[i])
# end


SSE2 = round(Int, sum(abs2, fit.resid))
RMSE2 = round(Int, sqrt(mean(abs2, fit.resid)))

ϵ, C = round.(Int, fit.param)

@show SSE1, SSE2, RMSE1, RMSE2, ϵ, C

# bp = data(df) * mapping(:elevation, :total_rotation) * visual(RainClouds; )
bp = data(df) * mapping(:elevation, :total_rotation) * visual(BoxPlot; show_outliers = false, width = 5)
elevation = 0:90
avg = data((; elevation, total_rotation = fill(μ, length(elevation)))) * mapping(:elevation, :total_rotation) * visual(Lines, color = :blue, label = "μ")
mdl = data((; elevation, total_rotation = get_α.(fit.param[1], fit.param[2], elevation))) * mapping(:elevation, :total_rotation) * visual(Lines, color = :red, label = "model")

fig = Figure(size = (8cm, 5cm), fontsize = 8pt, fonts = (; regular = "Helvetica"));
ax = Axis(fig[1, 1], xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xlabel = "Elevation (°)", ylabel = "Yaw rotation (°)", xticks = 0:30:90, yticks = 0:180:700)
grid = draw!(ax, bp + avg + mdl)
legend!(fig[1, 1], grid; tellheight=false, tellwidth=false, halign=:left, valign=:top, margin=(10,10,10,10), framevisible = false, orientation = :horizontal)
save("figure1.pdf", fig)

