using AlgebraOfGraphics, GLMakie, CairoMakie

using Dates, LinearAlgebra, Statistics
using CSV, DataFrames, CameraCalibrations
using Interpolations, StaticArrays, Dierckx, CoordinateTransformations, Rotations
using OhMyThreads
using LsqFit
using CategoricalArrays
using Distributions
using IntervalSets
using QuadGK

include("functions.jl")
# include("smooth.jl")

output = "dark"

const results_dir = "../track_calibrate/tracks and calibrations"

runs = CSV.read(joinpath(results_dir, "runs.csv"), DataFrame)
calibs = CSV.read(joinpath(results_dir, "calibs.csv"), DataFrame)
transform!(calibs, :calibration_id => ByRow(get_calibration) => :rectify)
transform!(calibs, [:rectify, :center_ij] => ByRow((f, c) -> passmissing(f)(totuple(c))) => :center)
transform!(calibs, [:rectify, :north_ij] => ByRow((f, c) -> passmissing(f)(totuple(c))) => :north)
select!(calibs, Cols(:calibration_id, :rectify, :center, :north))
leftjoin!(runs, calibs, on = :calibration_id)
select!(runs, Not(:runs_path, :start_location, :calibration_id, :fps, :target_width, :runs_file, :window_size))
transform!(runs, :dance => ByRow(x -> x ? "dance" : "no dance"), renamecols = false)
rename!(runs, :poi => :intervention)
transform!(runs, :spontaneous_end => ByRow(passmissing(tosecond)), renamecols = false)
transform!(runs, [:spontaneous_end, :intervention] => ByRow(coalesce) => :poi)
transform!(runs, [:tij_file, :rectify, :poi] => ByRow(get_txy) => [:t, :xy, :poi_index, :dance_jump])
transform!(runs, [:xy, :t] => ByRow(get_spline) => :spl)
transform!(runs, [:t, :spl, :poi_index] => ByRow(get_center_rotate) => :center_rotate)
transform!(runs, [:t, :spl, :center_rotate] => ByRow((t, spl, tform) -> tform.(SVector{2, Float64}.(spl.(t)))) => :sxy)
transform!(runs, :sxy => ByRow(xy -> cropto(xy, 50)); renamecols = false)
transform!(runs, [:t, :spl, :poi_index] => ByRow(get_turn_profile) => [:tθ, :θ])
transform!(runs, [:tθ, :θ] => ByRow(fit_logistic) => :ks)
transform!(runs, [:tθ, :ks] => ByRow((x, p) -> logistic.(x, p[2], p[1], 0)) => :θs, :ks => ByRow(first) => :k)
transform!(runs, :spontaneous_end => ByRow(!ismissing) => :spontaneous, :runs_dance => ByRow(!ismissing) => :induced)
transform!(runs, [:induced, :spontaneous] => ByRow((e, s) -> e ? "induced" : s ? "spontaneous" : "no") => :danced)
transform!(runs, [:t, :spl, :poi_index] => ByRow(get_mean_speed) => :speed)

runs.danced .= categorical(runs.danced; levels = ["no", "spontaneous", "induced"], ordered = true)



###########################################
##### Only shift ##########################

df = subset(runs, :light => ByRow(==("dark")))

###########################################
# spontaneous

df1 = subset(runs, :spontaneous)

df1 = flatten(df1, :sxy)
transform!(df1, :sxy => [:x, :y])

words = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"]

plt = data(df1) * mapping(:x => "X (cm)", :y => "Y (cm)", group=:run_id => nonnumeric, row = :at_run => renamer(words)) * visual(Lines)
fig = draw(plt; axis=(aspect=1, ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "spontaneous.png"), fig)

###########################################

df1 = flatten(df, :sxy)
transform!(df1, :sxy => [:x, :y])

plt = data(df1) * mapping(:x => "X (cm)", :y => "Y (cm)", group=:run_id => nonnumeric, col = :danced => sorter("no", "spontaneous", "induced"), row = :at_run => renamer(words)) * visual(Lines)
fig = draw(plt; axis=(aspect=1, ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure1.png"), fig)

###########################################

df1 = flatten(df, :sxy)
transform!(df1, :sxy => [:x, :y])

plt = data(df1) * mapping(:x => "X (cm)", :y => "Y (cm)", group=:run_id => nonnumeric, col = :danced => sorter("no", "spontaneous", "induced"), color = :at_run => renamer(words) => "at run") * visual(Lines)
fig = draw(plt; axis=(aspect=1, ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure2.png"), fig)

     #)##########################################

df1 = transform(df, [:xy, :poi_index] => ByRow((x, i) -> x[1:i]) => :xy, [:sxy, :poi_index] => ByRow((x, i) -> x[1:i]) => :sxy)
df1 = flatten(df1, [:sxy, :xy])
transform!(df1, [:xy, :center_rotate] => ByRow((p, f) -> f(p)) => :xy)
transform!(df1, :sxy => [:sx, :sy], :xy => [:x, :y])
d = data(df1)
m = mapping(group=:run_id => nonnumeric, col = :danced => sorter("no", "spontaneous", "induced"), row = :at_run => renamer(words))

plt = d * m * (mapping(:x, :y) * visual(Lines) + mapping(:sx, :sy) * visual(Lines, color = :red))
fig = draw(plt; axis=(aspect=DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)"), figure = (;size = (1000, 4000)))

save(joinpath(output, "start to POI.png"), fig)

###########################################

# df1 = deepcopy(df)
# df1.p .= Ref((0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1))
# df1 = flatten(df1, :p)
# transform!(df1, [:t, :spl, :poi_index, :p] => ByRow(get_turn_profile) => :turn)
#
# plt = data(df1) * mapping(:danced => sorter("no", "spontaneous", "induced"), :turn => "Turn (°)", row = :p => (x -> string(round(Int, 100x), " %")), col = :at_run => renamer(words)) * visual(Violin, datalimits=(0, Inf))
#
# fig = draw(plt; figure = (;size = (800, 1600)))
#
# save("figure3.png", fig)

###########################################

df2 = flatten(df, [:θ, :tθ])
plt = data(df2) * mapping(:tθ => "Time from POI (sec)", :θ => rad2deg => "Turn (°)", group=:run_id => nonnumeric, col = :danced => sorter("no", "spontaneous", "induced"), row = :at_run => renamer(words)) * visual(Lines)
fig = draw(plt; axis=(; limits = ((0,50), (-200, 200)), yticks = [-180, 0, 180]))
save(joinpath(output, "from 1 sec before POI.png"), fig)

###########################################

# df1 = vcat(df, transform(df, :at_run => ByRow(_ -> "all"), renamecols = false))
df1 = deepcopy(df)

df1.row_number .= 1:nrow(df1)
n = ceil(Int, sqrt(nrow(df1)))
transform!(df1, :row_number => ByRow(i -> Tuple(CartesianIndices((n, n))[i])) => [:row, :col])
df2 = flatten(df1, [:tθ, :θ, :θs])

plt = data(df2) * (mapping(:tθ => "Time from POI (sec)", :θ => rad2deg => "Turn (°)", col = :col => nonnumeric, row = :row => nonnumeric) * visual(Lines) + mapping(:tθ => "Time from POI (sec)", :θs => rad2deg => "Turn (°)", col = :col => nonnumeric, row = :row => nonnumeric) * visual(Lines, color = :red))
fig = draw(plt; figure = (;size = (2000, 2000)), axis = (; limits = ((0, 60), (nothing, 400))))

save(joinpath(output, "fits.png"), fig)

###########################################

plt = data(df1) * mapping(:danced => sorter("no", "spontaneous", "induced"), :k, row = :at_run => renamer(words)) * visual(RainClouds, violin_limits = extrema)
fig = draw(plt; axis = (; limits=(nothing, (nothing, 4))), figure =(; size = (400, 1000)))

save(joinpath(output, "k.png"), fig)

###########################################

plt = data(df) * mapping(:danced => sorter("no", "spontaneous", "induced"), :speed => "Speed (cm/s)", row = :at_run => renamer(words)) * visual(RainClouds, violin_limits = (0, Inf))
fig = draw(plt; figure =(; size = (400, 1000)))

save(joinpath(output, "speed1.png"), fig)

hist(df.speed; axis = ( ;xlabel = "Speed (cm/s)")) |> save(joinpath(output, "speed2.png"))

###########################################

df1 = combine(groupby(df, [:at_run, :danced]), :k => getIQR => [:c1, :μ, :c2])
tθ = range(0, 20, 1000)
transform!(df1, [:c1, :μ, :c2] .=> ByRow(k -> rad2deg.(logistic.(tθ, 2π, k, 0))) .=> [:θc1, :θμ, :θc2])
df1.tθ .= Ref(tθ)
df2 = flatten(df1, [:tθ, :θc1, :θμ, :θc2])
plt = data(df2) * mapping(:tθ => "Time from POI (sec)", :θμ, lower = :θc1, upper = :θc2, row = :at_run => renamer(words), color = :danced => sorter("no", "spontaneous", "induced") => "Dance") * visual(LinesFill) 
fig = draw(plt; axis = (; xscale = Makie.pseudolog10,  ylabel = "Turn", yticks = 0:45:180, ytickformat = "{:n}°"))#, limits = (extrema(tθ), (0, 180))))
save(joinpath(output, "mean k.png"), fig)

# plt = data(df2) * mapping(:tθ => "Time from POI (sec)", :θμ, lower = :θc1, upper = :θc2, color = :at_run => renamer("1" => "first", "4" => "forth", "10" => "tenth"), row = :danced => sorter("no", "spontaneous", "induced") => "Dance") * visual(LinesFill) 
# fig = draw(plt; axis = (; xscale = Makie.pseudolog10,  ylabel = "Turn", yticks = 0:45:180, ytickformat = "{:n}°"))#, limits = (extrema(tθ), (0, 180))))
# save("mean k.png", fig)

###########################################

df1 = combine(groupby(df, [:danced, :at_run]), :k => getIQR => [:c1, :μ, :c2])
transform!(df1, [:c1, :μ, :c2] .=> ByRow(create_track) .=> [:xyc1, :xyμ, :xyc2])
transform!(df1, [:xyc1, :xyc2] => ByRow((l, u) -> Makie.Polygon([l; reverse(u)])) => :poly)
df2 = flatten(df1, [:xyc1, :xyμ, :xyc2])
transform!(df2, :xyμ => [:xμ, :yμ])

# allwords = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"]
# words = [string(i) => allwords[i] for i in [1,4,10]]
# df2.at_run .= categorical(df2.at_run; levels = first.(words), ordered = true)

plt = data(df2) * (
                   mapping(:poly, row = :at_run => renamer(words), color = :danced => sorter("no", "spontaneous", "induced") => "Dance") * visual(Poly, alpha = 0.15)
                   + 
                   mapping(:xμ, :yμ, row = :at_run => renamer(words), color = :danced => sorter("no", "spontaneous", "induced") => "Dance") * visual(Lines)
                  )
fig = draw(plt; axis = (; aspect = DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)"))

save(joinpath(output, "mean tracks.png"), fig)



###########################################

df1 = subset(df, :danced => ByRow(==("no")))

df2 = sort!(df1, :k, rev = true)[1:3,:]
df2 = flatten(df2, [:θ, :tθ])
fig = data(df2) * mapping(:tθ => "Time from POI (sec)", :θ => rad2deg => "Turn (°)", group = :run_id => nonnumeric) * visual(Lines) |> draw

save(joinpath(output, "fastest no dance.png"), fig)


#############################


using GLM

df1 = select(subset(df, :danced => ByRow(≠("spontaneous"))), Cols(:danced, :k))
data(df1) * mapping(:danced, :k) * visual(RainClouds, violin_limits = (0, Inf)) |> draw

m1 = glm(@formula(k ~ danced), df1, Gamma())

m2 = glm(@formula(k ~ danced + at_run), df, Gamma())
m3 = glm(@formula(k ~ danced), df, Gamma())
lrtest(m1.model, m2.model)
lrtest(m2.model, m3.model)

newx = select(combine(groupby(df, :danced), :k => mean => :k), Not(:k))
k = predict(m3, newx; interval = :confidence, level = 0.95)
df1 = hcat(newx, k)

tθ = range(0, 20, 1000)
transform!(df1, [:lower, :prediction, :upper] .=> ByRow(k -> rad2deg.(logistic.(tθ, 2π, k, 0))) .=> [:θc1, :θμ, :θc2])
df1.tθ .= Ref(tθ)
df2 = flatten(df1, [:tθ, :θc1, :θμ, :θc2])

plt = data(df2) * mapping(:tθ => "Time from POI (sec)", :θμ, lower = :θc1, upper = :θc2, color = :danced => sorter("no", "spontaneous", "induced") => "Dance") * visual(LinesFill) 
fig = draw(plt; axis = (; xscale = Makie.pseudolog10,  ylabel = "Turn", yticks = 0:45:180, ytickformat = "{:n}°"))#, limits = (extrema(tθ), (0, 180))))

save(joinpath(output, "mean k.png"), fig)








##################################### DARK



# df1 = transform(runs, [:t, :spl, :poi_index] => ByRow(get_turn_profile) => [:tθ, :θ])
# df2 = flatten(df1, [:θ, :tθ])
# plt = data(df2) * mapping(:tθ => "Time from POI (sec)", :θ => rad2deg => "Turn (°)", group=:run_id => nonnumeric, col = :dance => nonnumeric, row = :at_run => nonnumeric) * visual(Lines)
# fig = draw(plt; figure = (;size = (700, 1000)), axis = (; limits = ((0, 20), (nothing, 400))));

# save("dark from 0.83 sec before POI.png", fig)







# f(x, L, k, x₀) = L / (1 + exp(-k*(x - x₀)))
#
# n = 100
# t = range(0, 10, n)
# L = 180
# k = 0.5
# x₀ = 0
# y = f.(t, L, k, x₀) .+ 0.02randn(n)
# lines(t, y)
#
# lb = Float64[1, 0.5, -2]
# ub = Float64[3, 2, 2]
# p0 = Float64[3, 1, 1]
# fit = curve_fit(model, t, y, p0, lower = lb, upper = ub)
#
# confidence_inter = confint(fit; level=0.95)
#
#
# row = runs[1,:]
# xy = row.xy
# t = row.t
# poi_index = row.poi_index
#
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# lines!(ax, xy)
#
# accumulate(t[poi_index - 24 : end]) do ti
#
#
#     kasjdhgfkhgfksdhak
#
#     function get_turn_profile(t, spl, poi_index)
#         # l = -3
#         # row = df[1,:]
#         # t = row.t
#         # spl = row.spl
#         # poi_index = row.poi_index
#         o = optimize(t1 -> abs2(first(quadgk(t -> norm(derivative(spl, t)), t1, t[poi_index])) - 3), t[1], t[poi_index])
#         t1 = o.minimizer
#         o = optimize(t2 -> abs2(first(quadgk(t -> norm(derivative(spl, t)), t[poi_index], t2)) - l), t[poi_index], t[end])
#         t2 = o.minimizer
#         ts = range(t1, t[end], step = 1/30)
#         der = derivative.(Ref(spl), ts)
#         θ = [atan(reverse(d)...) for d in der]
#         unwrap!(θ)
#         θ .-= θ[1]
#         # lines(ts, θ, axis = (; limits=((t1, t[end]), nothing)))
#         spl1 = Spline1D(ts, θ)
#         θ1 = spl1(t2)
#         θtotal = spl1(ts[end])
#         (; θ1, θtotal)
#     end
#
#
#
#     abtrace = data((; intercept = [0], slope = [1]))  * mapping(:intercept, :slope) * visual(ABLines, color = :gray)
#     df1 = subset(df, :light => ByRow(==("shift")))
#     # df1 = vcat(df1, transform(df1, :at_run => ByRow(_ -> "pooled"); renamecols = false))
#     df1 = map(0:3:12) do h
#         df2 = transform(df1, [:t, :spl, :poi_index] => ByRow((args...) -> get_turn_profile(args..., h)) => [:θ1, :θtotal])
#         df2.h .= h
#         df2
#     end |> splat(vcat)
#
#
#     plt = abtrace + data(df1) * mapping(:θ1 => rad2deg => "Turn between 3 cm before and h cm after the POI (°)", :θtotal => rad2deg => "Total turn (°)", col = :dance => nonnumeric, row = :h => nonnumeric, color = :at_run) * visual(Scatter)
#     fig = draw(plt; axis = (; yticks = [-180, 0, 180], xticks = [-180, 0, 180]))
#
#     save("figure1.png", fig)
#
#
#
#     function get_turn_profile(t, spl, poi_index)
#         # l = -3
#         # row = df[1,:]
#         # t = row.t
#         # spl = row.spl
#         # poi_index = row.poi_index
#         o = optimize(t1 -> abs2(first(quadgk(t -> norm(derivative(spl, t)), t1, t[poi_index])) - 3), t[1], t[poi_index])
#         t1 = o.minimizer
#         ts = range(t1, t[end], step = 1/30)
#         der = derivative.(Ref(spl), ts)
#         θ = [atan(reverse(d)...) for d in der]
#         unwrap!(θ)
#         θ .-= θ[1]
#         # lines(ts, θ, axis = (; limits=((t1, t[end]), nothing)))
#         (; tθ = ts .- ts[1], θ = abs.(θ))
#     end
#     df1 = subset(df, :light => ByRow(==("shift")))
#     transform!(df1, [:t, :spl, :poi_index] => ByRow(get_turn_profile) => [:tθ, :θ])
#
#     df1 = flatten(df1, [:θ, :tθ])
#
#     plt = data(df1) * mapping(:tθ => "Time from 3 cm before the POI (s)", :θ => rad2deg => "Turn (°)", group=:run_id => nonnumeric, col = :dance => nonnumeric, row = :at_run => nonnumeric) * visual(Lines)
#     fig = draw(plt; figure = (;size = (700, 1000)), axis=(; limits = ((0, 30),(0, 300)) ))
#
#     save("figure2.png", fig)
#




# function method1(spl, t1, t2)
#     tl = LinRange(t1, t2, 100_000)
#     xy = spl(tl)
#     Δxy = xy[:,2:end] - xy[:,1:end-1]
#     sum(norm.(eachcol(Δxy)))
# end
# function method2(spl, t1, t2)
#     res, _ = quadgk(t -> norm(derivative(spl, t)), t1, t2, rtol = 1e-8)
#     res
# end
# function method2a(spl)
#     knots = Dierckx.get_knots(spl)
#     s = 0.0
#     for (k1, k2) in zip(knots[1:end-1], knots[2:end])
#         res, _ = quadgk(t -> norm(derivative(spl, t)), k1, k2)
#         s += res
#     end
#     return s
# end
# function method3(spl, x_gauss, w_gauss)
#     knots = Dierckx.get_knots(spl)
#     integral = sum(eachindex(knots)[2:end]) do i
#         t1, t2 = knots[i-1], knots[i]
#         scale = (t2 - t1) / 2
#         sum(zip(x_gauss, w_gauss)) do ((xg, wg))
#             t = (xg + 1) * scale + t1 # map from (-1,1) to (t1,t2)
#             norm(derivative(spl, t)) * wg
#         end * scale
#     end
# end
# function method4(spl; kws...)
#     knots = Dierckx.get_knots(spl)
#     s, _ = quadgk(t -> norm(derivative(spl, t)), knots; kws...)
#     return s
# end
#
#
#
# row = collect(eachrow(df))[5]
# spl = row.spl
# t = row.t
# s, _ = quadgk(t -> norm(derivative(spl, t)), extrema(t)..., rtol = 1e-9)
# rel(x) = (x - s)/s
#
# rel(method1(spl, extrema(t)...))
# rel(method2(spl, extrema(t)...))
# rel(method2a(spl))
# rel(method3(spl, gauss(15)...))
# rel(method4(spl, order=4, rtol=1e-5))
#
#
# @btime method1($spl, extrema(t)...)
# @btime method2($spl, extrema(t)...)
# @btime method2a($spl)
# @btime method3($spl, $gauss(15)...)
# @btime method4($spl, order=4, rtol=1e-5)
#


# row = collect(eachrow(df))[5]
# spl = row.spl
# t = row.t
# poi_index = row.poi_index
# lines(Point2f.(spl.(t)), axis = (;aspect = DataAspect()))
# scatter!(Point2f.(spl.(t[poi_index])))
#
# t2 = t[poi_index:end]
# der = derivative.(Ref(spl), t2)
# θ = [atan(reverse(d)...) for d in der]
# unwrap!(θ)
# spl1 = Spline1D(t2, θ)
# S = Chebyshev(ApproxFun.ClosedInterval(extrema(t2)...));
# p = points(S, length(t2));
# v = spl1.(p);  
# f = Fun(S,ApproxFun.transform(S,v));
# lines(t2, θ)
# lines!(t2, f.(t2))
#
#
# ff = cumsum(f')
# lines(t2, ff.(t2))
#
# ff(t2[1] + 2) - ff(t2[1])
# ff(t2[end])
#
# f = Fun...
# F = cumsum(f)
# a, b = (1, 2)
# s = F(b) - F(a)
#
# S = 0..π
# n = 100
# pts = ApproxFun.points(S, n)
# t = range(0, pi, 50)
# xy = reverse.(sincos.(t)) 
# spl = ParametricSpline(t, stack(xy); k = 2, s = 0)
# der = derivative.(Ref(spl), pts)
# θ = [atan(reverse(d)...) for d in der]
# unwrap!(θ)
# lines(rad2deg.(θ))
#
# S = Chebyshev()
#
# nfun = Fun(S, ApproxFun.transform(S, θ))
#
# t = range(0, π, 100)
# lines(nfun)
#
#
#
# using Statistics
# mean_angle(θ) = angle(mean(exp, θ*im))
#
#
# function plot_direction(poi_index, spl, t, run_id)
#     # row = df[2,:]
#     # poi_index = row.poi_index
#     # spl = row.spl
#     # t = row.t
#     # run_id = row.run_id
#     # GLMakie.activate!()
#     der = derivative.(Ref(spl), t)
#     fig = Figure()
#     ax = Axis(fig[1,1], autolimitaspect = 1, aspect = DataAspect(), xlabel = "X (cm)", ylabel = "Y (cm)")
#     lines!(ax, Point2f.(spl.(t[1:poi_index])))
#     lines!(ax, Point2f.(spl.(t[poi_index:end])))
#     ax = Axis(fig[2,1], limits = (extrema(t), nothing), xlabel = "Time (sec)", ylabel = "Direction (°)")
#     θ = [atan(reverse(d)...) for d in der]
#     unwrap!(θ)
#     lines!(ax, t[1:poi_index], rad2deg.(θ[1:poi_index]))
#     lines!(ax, t[poi_index:end], rad2deg.(θ[poi_index:end]))
#     θs = [mean_angle(θ[i]) for i in (1:poi_index, poi_index:length(t))]
#     lines!(ax, t[[1, poi_index, poi_index, end]], [fill(rad2deg(θs[1]), 2); fill(rad2deg(θs[2]), 2)], color=:gray)
#     ax = Axis(fig[3,1], limits = (extrema(t), nothing), xlabel = "Time (sec)", ylabel = "Difference in direction (°)")
#     lines!(ax, t[2:end], diff(θ))
#     display(fig)
#     # save(joinpath("directions", string(run_id, ".png")), fig)
# end
#
# row = df[15,:]
# plot_direction(row.poi_index, row.spl, row.t, row.run_id)
#
#
# if isdir("directions")
#     rm("directions", recursive=true)
# end
# mkpath("directions")
# CairoMakie.activate!()
# @tasks for row in eachrow(df)
#     plot_direction(row.poi_index, row.spl, row.t, row.run_id)
# end
#
# GLMakie.activate!()
# function turning_event(t, spl, poi_index)
#     der = derivative.(Ref(spl), t)
#     θ = [atan(reverse(d)...) for d in der]
#     unwrap!(θ)
#     # d, i = findmax(diff(θ[poi_index:end]))
#     d, i = findmax(abs.(diff(θ)))
#     # (; d, dt = t[i + poi_index - 1] - t[poi_index])
#     (; d = rad2deg(d), dt = t[i] - t[poi_index])
# end
# df = select(runs, [:tij_file, :rectify, :poi] => ByRow(get_spline) => [:t, :spl, :poi_index], Cols(:run_id, :dance, :light, :at_run));
# transform!(df, [:t, :spl, :poi_index] => ByRow(turning_event) => [:d, :dt]);
#
# plt = data(df) * mapping(:dt => "Time from POI (sec)", :d => "Turn (°)", col = :dance => nonnumeric, row = :light => nonnumeric) * AlgebraOfGraphics.density(npoints=30) * visual(Contour)
# fig = draw(plt)
#
# plt = data(df) * mapping(:dance, :d => "Turn (°)", color = :light, dodge = :light, row = :at_run) * visual(Violin, datalimits = (0, 180))
# fig = draw(plt)
#



# display(fig)


#
#
#
# transform(groupby(runs, :run_id), [:tij_file, :start, :POI, :stop, :rectify] => ((f, s, p, ss, r) -> (@show length(p))))# => [:t, :xy, :POI])
#
#
# df = flatten(select(transform(runs, [:tij_file, :start, :POI, :stop, :rectify] => ByRow(get_txy) => :txy), Not(:start, :stop, :POI)), :txy)
# transform!(df, :txy => [:t, :xy])
# select!(df, Not(:txy))
# #
# # just ignore the runs that are split differently, there's just one
# subset!(transform!(groupby(df, :run_id), nrow), :nrow => ByRow(==(2)))
# @assert all(==(2), combine(groupby(df, :run_id), nrow).nrow)
# #
# transform!(df, [:t, :xy] => ByRow(smooth) => [:t, :xy])
# transform!(groupby(df, :run_id), [:north, :center, :xy] => register!; renamecols = false)
# #
#
# mkpath("tracks")
# gs = collect(groupby(df, :run_id))
# @tasks for i in 1:length(gs)
#     g = gs[i]
#     k = (; run_id = g.run_id[1])
#     fig  = Figure()
#     ax = Axis(fig[1,1], aspect = DataAspect())
#     for r  in (300, 500)
#         lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
#     end
#     for g1 in eachrow(g)
#         lines!(ax, g1.xy)
#     end
#     save(joinpath("tracks", string(k.run_id, ".png")), fig)
# end
#
#

















#
# using SimpTrack
# using Dates, LinearAlgebra
# using Dierckx, VideoIO, CameraCalibrations, OhMyThreads, ImageTransformations, ImageDraw, Colors
# using StaticArrays, CoordinateTransformations, Rotations
# # using CairoMakie
# using GLMakie
# using AlgebraOfGraphics
#
# const SV = SVector{2, Float64}
#
# include("calibrations.jl")
#
# calib_file = "/home/yakir/tmp/Elin_tracks/calib.csv"
# runs_file = "/home/yakir/tmp/Elin_tracks/runs.csv"
#
# # Calibrations
#
# calibrations = CSV.read(calib_file, DataFrame)
#
# # data quality checks
# @assert allunique(calibrations.calibration_id) "Calibration IDs not identical"
# for row in eachrow(calibrations)
#     @assert row.start < row.stop "start must occur before stop in calibration $(row.calibration_id)"
#     file = joinpath(row.path, row.file)
#     @assert isfile(file) "calibration file $file does not exist"
#     # @assert that the time stamps are within the period of the video
# end
#
# @rselect!(calibrations, :calibration_id, :calib = calib(:calibration_id, joinpath(:path, :file), :start, :stop, :extrinsic))
#
# # # save calibrations images for quality assesment
# # suspect_calibration = "20220304_calibration.mov 2nd"
# # check_calibration(calibrations, suspect_calibration)
#
# # Tracking
# df = CSV.read(runs_file, DataFrame)
#
# # data quality checks
# for row in eachrow(df)
#     @assert row.start ≤ row.POI ≤ row.stop "POI is not within the track in row $row"
#     @assert row.calibration_id ∈ calibrations.calibration_id "Calibration $calibration_id is missing from the calibrations file."
#     file = joinpath(row.path, row.file)
#     @assert isfile(file) "run file $file does not exist"
#     # @assert that the time stamps are within the period of the video
# end
#
# df.track = tmap(track, joinpath.(df.path, df.file), df.start, df.stop)
#
# # Combine the two
# leftjoin!(df, calibrations, on = :calibration_id)
# select!(df, Not(:calibration_id))
#
#
#
#
#
# # s = 100
# # row = first(eachrow(df))
# # runi, clib, trck, start, stop, POI = (1, row.calib, row.track, row.start, row.stop, row.POI)
# # tfm = get_transformation(clib, trck, start, POI; s)
# # # xy = tfm.(range(start, stop, step = Millisecond(100)))
# # # lines(xy)
# #
# # t = range(start, stop, step = Second(1))
# #
# # trajectory_length = cumsum(norm.(diff(tfm.(t))))
# # ca = [cumulative_angle(tfm, start, tstop) for tstop in t[2:end]]
# # total_cumulative_angle = round(Int, cumulative_angle(tfm, start, stop))
#
#
# ns = 5
# ss = round.(Int, range(50, 250, ns))
# function explore(runi, calib, track, start, stop, POI)
#     w = 300
#     fig = Figure(size=(2w, ns*w))
#     axs1 = []
#     axs2 = []
#     for (j, s) in enumerate(ss)
#         tfm = get_transformation(calib, track, start, POI; s)
#         t = range(start, stop, step = Second(2))
#         trajectory_length = cumsum(norm.(diff(tfm.(t))))
#         ca = [cumulative_angle(tfm, start, tstop) for tstop in t[2:end]]
#         total_cumulative_angle = round(Int, cumulative_angle(tfm, start, stop))
#         ax = Axis(fig[j, 1], aspect=DataAspect(), ylabel = "s = $s")
#         for r  in (30, 50)
#             lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
#         end
#         xy = tfm.(range(start, POI, step = Millisecond(100)))
#         lines!(ax, xy)
#         xy = tfm.(range(POI, stop, step = Millisecond(100)))
#         lines!(ax, xy)
#         j ≠ ns && hidexdecorations!(ax, grid = false, minorgrid = false)#, ticks = false, minorticks = false)
#         push!(axs1, ax)
#         ax = Axis(fig[j, 2], xlabel = "Trajectory length (cm)", yaxisposition = :right, aspect = AxisAspect(1), ylabel = "Cumulative angle (°)", title = "Final cumulative angle = $total_cumulative_angle °")
#         lines!(ax, trajectory_length, ca)
#         j ≠ ns && hidexdecorations!(ax, grid = false, minorgrid = false)#, ticks = false, minorticks = false)
#         push!(axs2, ax)
#     end
#     # linkaxes!(filter(x -> x isa Axis, fig.content)...)
#     linkaxes!(axs1...)
#     linkaxes!(axs2...)
#     save("run $runi.png", fig)
# end
# foreach(enumerate(eachrow(df))) do (runi, row)
#     explore(runi, row.calib, row.track, row.start, row.stop, row.POI)
# end
#
# # n = 10
# # α = rand(-10:10, n)
# #
#
# # TODO: ask Vishaal about the noise effect on cumulative turn angle
# function angle_between(p1, p2)
#     θ = acos(dot(p1, p2) / norm(p1) / norm(p2))
#     return -sign(cross(p1, p2))*θ
# end
# function fun(xy)
#     δ = diff(xy)
#     α = sum(splat(angle_between) ∘ reverse, partition(δ, 2, 1))
#     rad2deg(α)
#     # pushfirst!(α, atan(reverse(δ[1])...))
#     # rad2deg(sum(splat(angle_between), partition(δ, 2, 1)))
# end
# n = 1000
# xy = [SV(reverse(sincos(θ))) for θ in range(0, π, n)]
# fig, ax, h = lines(xy, label = "clean",  axis = (;aspect = DataAspect()))
# xy1 = xy .+ 0.01randn(SV, n)
# lines!(ax, xy1, label = "with noise")
# xy2 = xy1[round.(Int, range(1, n, 15))]
# # scatterlines!(ax, xy2, label = "subsampled")
# δ = diff(xy2)
# α = map(splat(angle_between) ∘ reverse, partition(δ, 2, 1))
# pushfirst!(α, atan(reverse(δ[1])...))
# l = [SV(reverse(sincos(i))) for i in cumsum(α)]
# arrows!(ax, Point2f.(xy2[1:end-1]), l .* norm.(δ), label = "subsampled")
# axislegend(ax)
# fun(xy)
# fun(xy1)
# fun(xy2)
# #
# # using Statistics
# #
# # steps = range(1, 100, 101)
# # m = 1000
# # n = 1000
# # Δ = zeros(length(steps), m)
# # for mi in 1:m
# #     xy = [SV(reverse(sincos(θ))) for θ in range(0, π, n)]
# #     xy1 = xy .+ 0.01randn(SV, n)
# #     Δ[:,mi] = map(steps) do step
# #         xy2 = xy1[round.(Int, range(1, n; step))]
# #         abs(180 - fun(xy2))
# #     end
# # end
# # lines(steps, vec(mean(Δ, dims=2)))
# #
# #
#
#
#
#
#
#
#
#
#
# #
# # δα = [0, 0, 45, 45, 45 , 45, 0]
# # α = cumsum(δα)
# # l = [SV(reverse(sincosd(i))) for i in α]
# # xy = cumsum(l)
# # arrows!(ax, Point2f.(xy[1:end-1]), l[2:end])
# #
# # Δ = diff(α)
# # sum(Δ[2:end])
# #
# #
# # cumulative_angle(xy)
# #
