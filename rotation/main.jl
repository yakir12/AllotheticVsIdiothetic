using Dates, LinearAlgebra, Statistics, Random
using CSV, DataFrames
using DataFramesMeta, Chain
using AlgebraOfGraphics, GLMakie

df = @chain "rotation_all_data.csv" begin
    CSV.read(DataFrame; missingstring="NA")
    @subset :person_extracted .== "Elin" :year_extracted .== 2025 :condition .== "LED"
    dropmissing(:total_rotation)
    # @rtransform :elevation = :elevation == 90 ? 89.9 : :elevation
end

# data(df) * mapping(:elevation, :total_rotation) * visual(BoxPlot; width = 5) |> draw(; axis = (; xticks = 0:10:90))


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

using LsqFit


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





bp = data(df) * mapping(:elevation, :total_rotation) * visual(BoxPlot; show_outliers = false, width = 5, label = "data")
elevation = 0:90
avg = data((; elevation, total_rotation = fill(μ, length(elevation)))) * mapping(:elevation, :total_rotation) * visual(Lines, color = :blue, label = "μ")
mdl = data((; elevation, total_rotation = get_α.(fit.param[1], fit.param[2], elevation))) * mapping(:elevation, :total_rotation) * visual(Lines, color = :red, label = "model")
fig = (bp + avg + mdl) |> draw(; axis = (; xticks = 0:30:90))
save("figure1.png", fig)

SSE1, SSE2, RMSE1, RMSE2
