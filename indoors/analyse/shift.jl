using AlgebraOfGraphics, GLMakie, CairoMakie

using HypothesisTests

# using StatsBase, Graphs
using Dates, LinearAlgebra, Statistics, Random
using CSV, DataFrames, CameraCalibrations
using Interpolations, StaticArrays, Dierckx, CoordinateTransformations, Rotations
using OhMyThreads
using LsqFit
using CategoricalArrays
using Distributions
using IntervalSets
using QuadGK
using BetaRegression

GLMakie.activate!()

include("minimal_functions.jl")

output = "figures"
mkpath("figures")

const results_dir = "../track_calibrate/tracks and calibrations"

function combine_factors(light, induced, run)
    induced = induced ? " induced" : ""
    run = run > 1 ? " $run" : ""
    string(light, induced, run)
end

runs = CSV.read(joinpath(results_dir, "runs.csv"), DataFrame)
transform!(runs, :spontaneous_end => ByRow(!ismissing) => :dance_spontaneous)
transform!(runs, [:light, :dance_induced, :at_run] => ByRow(combine_factors) => :condition)
# runs.condition .= categorical(runs.condition; levels = ["remain", "dark", "dark induced", "shift", "shift induced"], ordered = true)

calibs = CSV.read(joinpath(results_dir, "calibs.csv"), DataFrame)
transform!(calibs, :calibration_id => ByRow(get_calibration) => :rectify)
select!(calibs, Cols(:calibration_id, :rectify))
leftjoin!(runs, calibs, on = :calibration_id)
select!(runs, Not(:runs_path, :start_location, :calibration_id, :fps, :target_width, :runs_file, :window_size))

rename!(runs, :poi => :intervention)
transform!(runs, :spontaneous_end => ByRow(passmissing(tosecond)), renamecols = false)
transform!(runs, [:spontaneous_end, :intervention] => ByRow(coalesce) => :poi)

transform!(runs, :tij_file => ByRow(get_tij) => [:t, :ij])
# transform!(runs, [:t, :ij] => ByRow(clean_coords) => [:t, :ij])
transform!(runs, [:ij, :rectify] => ByRow((ij, fun) -> fun.(ij)) => :xy)



df = subset(runs, :light => ByRow(==("shift")))
df1 = subset(df, :run_id => ByRow(==(8)))

function remove_stops(t_range::AbstractRange,, xy)
    tokill = Int[]
    last_xy = xy[1]
    for i in 2:length(t)
        Δ = norm(xy[i] - last_xy)
        if Δ > 0.1
            last_xy = xy[i]
        else
            push!(tokill, i)
        end
    end
    t = collect(t_range)
    deleteat!(t, tokill)
    deleteat!(xy, tokill)
    return (t, xy)
end

function clean_track(t, xy)
    t, xy = remove_stops(t, xy)
    n = length(xy)
    t = 1:n
    dierckx_spline = ParametricSpline(t, stack(xy), k = 3, s = 100)
    t = range(1, n, 100n)
    n = length(t)
    xys = Point2.(dierckx_spline.(t))


end

xy = copy(df1.xy[])
δ = norm.(diff(xy))
tokill = findall(δ .< 0.1)
deleteat!(xy, tokill)



using GeometryBasics

n = length(xy)
t = 1:n
dierckx_spline = ParametricSpline(t, stack(xy), k = 3, s = 100)
t = range(1, n, 100n)
n = length(t)
xys = Point2.(dierckx_spline.(t))

lines(xys, axis = (; aspect = DataAspect(), limits = ((-50, 30), (0, 30))))

inds, _ = self_intersections(xys)
pushfirst!(inds, 1)
push!(inds, n)

t = t[vcat(splat(UnitRange).(Iterators.partition(inds, 2))...)]

n = length(t)
xys = dierckx_spline.(t)

scatter!(Point2f.(xys))

interpolations_spline = interpolate(stack(xys)', (BSpline(Cubic(Natural(OnGrid()))), NoInterp()))

t = range(1, n, 100000)

xys = [SV(interpolations_spline(i, 1:2)) for i in t]

lines!(xys)

kdsahgsahglsh



interpolations_spline = interpolate(stack(xys)', (BSpline(Cubic(Natural(OnGrid()))), NoInterp()))
interpolations_spl(t) = interpolations_spline(t, 1:2)

xys2 = interpolations_spl.(t)

lines(xy[1:10:end], axis = (; aspect = DataAspect(), limits = ((-50, 30), (0, 30))))

lines!(SV.(xys))
lines!(SV.(xys2))

# scatter!(xys[todo])

# θ = [atan(reverse(derivative(spl, ti))...) for ti in t]
# unwrap!(θ)

# todo = findall(>(π) ∘ abs, θ)

# lines(t, rad2deg.(θ))
# scatter!(t[todo], rad2deg.(θ[todo]))



using Optim, StaticArrays, Dierckx

function detect_self_intersection(spl, t1, t2)
    min_step = 10
    lx = Float64[t1, min_step]
    ux = Float64[t2 - min_step, t2 - t1]
    con_c!(c, x) = sum!(c, x)
    lc = [t1 + min_step]
    uc = [t2]
    dfc = TwiceDifferentiableConstraints(con_c!, lx, ux, lc, uc, :forward)

    x0 = (2lx .+ ux) ./ 3

    # x0 = [(t1 + t2 - min_step)/2, (t2 - (t1 + t2 - min_step)/2)/2]

    fun(x) = norm(spl(x[1]) - spl(sum(x)))
    df = TwiceDifferentiable(fun, x0; autodiff=:forward)

    res = optimize(df, dfc, x0, IPNewton())
    t, step = Optim.minimizer(res)
    return (t1 = t, t2 = t + step)
end

t1, t2 = detect_self_intersection(interpolations_spl, 1, 100)





# transform!(runs, [:t, :xy] => ByRow(prune_coords!) => [:t, :xy])
transform!(runs, [:t, :xy, :poi] => ByRow(impute_poi_time) => :poi)
disallowmissing!(runs, :poi)
transform!(runs, [:xy, :t, :poi] => ByRow(glue_poi_index!) => [:poi_index, :dance_jump])
transform!(runs, [:xy, :t] => ByRow(get_spline) => :spl)

transform!(runs, :t => ByRow(x -> similar(x, Float64)) => :l)
Threads.@threads for row in eachrow(runs)
    row.l = get_pathlength(row.t, row.spl)
end

transform!(runs, [:t, :spl] => ByRow(smooth_center) => :xyc)
transform!(runs, [:t, :spl, :poi_index] => ByRow(get_smooth_center_poi_rotate) => :tform)

transform!(runs, [:t, :spl, :poi_index] => ByRow(get_turn_profile) => [:tθ, :θ])
# transform!(runs, [:tθ, :θ] => ByRow(fit_logistic) => :ks)
# transform!(runs, [:tθ, :ks] => ByRow((x, p) -> logistic.(x, p[2], p[1], 0)) => :θs, :ks => ByRow(first) => :k)

df = subset(runs, :light => ByRow(==("shift")))

my_renamer = uppercasefirst ∘ string


df1 = subset(df, :run_id => ByRow(==(8)))

fig = pregrouped(df1.xyc => first, df1.xyc => last) * visual(Lines) |> draw(; axis = (; aspect = DataAspect()))
display(fig)

sdjkhfgksdljfhsfhj
# fig = (pregrouped(df1.xyc => first, df1.xyc => last) + pregrouped(df1.xy => first, df1.xy => last)) * visual(Lines) |> draw(; axis = (; aspect = DataAspect()))

pregrouped(df.tθ, df.θ => rad2deg, row = df.dance_induced, col = df.at_run => nonnumeric) * visual(Lines) |> draw()

pregrouped(df.tθ, df.θ => rad2deg, layout = string.(df.run_id, " ", df.condition)) * visual(Lines) |> draw()

fig



t = df1.t[1]
spl = df1.spl[1]
t1 = 30:0.04:40
xy = SV.(spl.(t1))
lines(xy)
scatter!(SV(spl(31.31)))
scatter!(SV(spl(39.52)))

f = Figure()
ax1 = Axis(f[1, 1], yticklabelcolor = :blue, ylabel = "x")
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, ylabel = "y")
hidespines!(ax2)
hidexdecorations!(ax2)
lines!(ax1, t1, first.(xy), color = :blue)
lines!(ax2, t1, last.(xy), color = :red)
vlines!(ax2, [31.31, 39.52])


function arclength(spl, t1, t2; kws...)
    knots = get_knots(spl)
    filter!(t -> t1 < t < t2, knots)
    pushfirst!(knots, t1)
    push!(knots, t2)
    s, _ = quadgk(t -> norm(derivative(spl, t)), knots; kws...)
    return s
end
xy = SV.(spl.(t))

sdjkhfgksdljfhsfhj

# ij = [(0,0), (1,0), (2,0), (3,1), (3,2), (1,2), (1,1), (2,0), (3,0), (4,0)]
# remove_cycles!(ij)
#
# ids = levelsmap(ij)
# vs = [ids[k] for k in ij]
# g = DiGraph(length(vs))
# for (v1, v2) in zip(vs[1:end-1], vs[2:end])
#     add_edge!(g, v1, v2)
# end
# c = only(simplecycles(g))
#
# xy = Point2f.(ij)
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# scatterlines!(ax, xy)
# scatter!(ax, xy[c], color = :red)
# deleteat!(xy, c)
# lines!(ax, xy, color = :green)
#
# deleteat!(ij, c)
#



df1 = subset(df, :run_id => ByRow(==(8)))
transform!(df1, :tij_file => ByRow(get_ij) => [:t, :ij])

# fig = pregrouped(df1.ij => first, df1.ij => last) * visual(ScatterLines) |> draw(; axis = (; aspect = DataAspect()))



ij = copy(df1.ij[1])
t = copy(df1.t[1])
tokill = findall(==(0) ∘ norm, diff(SV.(ij))) .+ 1
deleteat!(ij, tokill)
deleteat!(t, tokill)
t, ij = clean_coords(t, ij)
lines(SV.(ij), axis = (; aspect = DataAspect()))
# lines!(SV.(ij))




using Interpolations

itp = interpolate(a, BSpline(Constant()))

tokill = findall(==(0) ∘ norm, diff(SV.(ij))) .+ 1
deleteat!(ij, tokill)
lines(SV.(ij), axis = (; aspect = DataAspect()))
while remove_cycles!(ij)
end
lines!(SV.(ij))



function remove_cycles!(ij)
    ids = StatsBase.levelsmap(ij)
    vs = [ids[k] for k in ij]
    n = length(vs)
    g = DiGraph(n)
    for (v1, v2) in zip(vs[1:end-1], vs[2:end])
        add_edge!(g, v1, v2)
    end
    c = simplecycles(g)
    tokill = sort(unique(reduce(vcat, c)))
    deleteat!(ij, tokill)
    return ij
end

remove_cycles!(vs)

lines(SV.(vs), axis = (; aspect = DataAspect()))
lines!(SV.(df1.ij[1]))

dsjhflskdhfkdshflskhf


using Graphs

g = path_graph(60)

g = path_digraph(60)
graphplot(g)

using GLMakie, GraphMakie
using GraphMakie.NetworkLayout
g = smallgraph(:dodecahedral)
graphplot(g; layout=Stress(; dim=3))



ujfhsdfhjsdkh

######################## Plots 


# my_renamer = renamer("remain" => "Remain",
#                      "dark" => "Dark",
#                      "dark 10" => "Dark 10",
#                      "dark induced" => "Dark induced")

transform!(groupby(df, :condition), eachindex => :n)
combine(groupby(df, :condition), nrow).nrow



######################## Figure 1

df1 = subset(df, :xyc => ByRow(≥(50) ∘ norm ∘ last))
transform!(df1, :xyc => ByRow(xy -> cropto(xy, 50)) => :xyc)
transform!(groupby(df1, :condition), eachindex => :n)
subset!(df1, :n => ByRow(≤(10)))
@assert all(==(10), combine(groupby(df1, :condition), nrow).nrow)

fig = pregrouped(df1.xyc => first, df1.xyc => last, col = df1.condition => my_renamer) * visual(Lines) |> draw(; figure = (; size = (1200, 400)), axis=(aspect=DataAspect(), ))
for ax in fig.figure.content 
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save(joinpath(output, "figure4.png"), fig)


######################## Figure 2

df1 = transform(df, [:t, :poi_index, :tform] => ByRow((t, i, f) -> f.(t[i:end])) => :xyp)

fig = pregrouped(df1.xyp => first => "X (cm)", df1.xyp => last => "Y (cm)", col = df1.condition => my_renamer) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))

save(joinpath(output, "figure5.png"), fig)

################################################### Figure 3

df1 = transform(df, [:t, :poi_index, :tform] => ByRow((t, i, f) -> f.(t[i:end])) => :xyp)
nr = 3
l = floor(Int, minimum(norm ∘ last, df1.xyp))
r = range(1e-3, l, nr)
transform!(df1, :xyp => ByRow(xyp -> get_exit_angle.(Ref(xyp), r)) => :θs)
select!(df1, Cols(:condition, :θs))
df1.r .= Ref(r)

rl = range(extrema(r)..., 100)
conditions = levels(df1.condition)
newdf = DataFrame(r = repeat(rl, outer = length(conditions)), condition = repeat(conditions, inner = 100), mean_resultant_vector = zeros(100length(conditions)))
fm = @formula(mean_resultant_vector ~ r*condition)
function _bootstrap(df)
    n = nrow(df)
    df1 = flatten(df[sample(1:n, n), :], [:θs, :r])
    df2 = combine(groupby(df1, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
    try
        m = BetaRegression.fit(BetaRegressionModel, fm, df2)
        predict(m, newdf)
    catch ex
        missing
    end
end
function bootstrap(df)
    n = 10_000
    ys = Vector{Vector{Float64}}(undef, n)
    i = 0
    while i < n
        y = _bootstrap(df)
        if ismissing(y)
            continue
        else
            i += 1
            ys[i] = y
        end
    end
    return stack(ys)
end
c = bootstrap(df1)
y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdf.lower .= getindex.(y, 1)
newdf.mean_resultant_vector .= getindex.(y, 2)
newdf.upper .= getindex.(y, 3)

fig = data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :condition => my_renamer => "Light") * visual(LinesFill) |> draw(; axis = (; ylabel = "Mean resultant length", xlabel = "Radius (cm)", limits = ((0, l), (0, 1))))

save(joinpath(output, "figure6.png"), fig)
