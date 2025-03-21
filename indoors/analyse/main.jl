using AlgebraOfGraphics, GLMakie, CairoMakie

using HypothesisTests

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

transform!(runs, [:tij_file, :rectify] => ByRow(get_txy) => [:t, :xy])
transform!(runs, :xy => ByRow(clean_coords!) => :xy)
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

######################## Plots 

df = subset(runs, :light => ByRow(≠("shift")))

# my_renamer = uppercasefirst ∘ string
my_renamer = renamer("remain" => "Remain",
                     "dark" => "Dark",
                     "dark 10" => "Dark 10",
                     "dark induced" => "Dark induced")

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

save(joinpath(output, "figure1.png"), fig)


######################## Figure 2

df1 = transform(df, [:t, :poi_index, :tform] => ByRow((t, i, f) -> f.(t[i:end])) => :xyp)

fig = pregrouped(df1.xyp => first => "X (cm)", df1.xyp => last => "Y (cm)", col = df1.condition => my_renamer) * visual(Lines) |> draw(; figure = (; size = (1200, 300)), axis=(aspect=DataAspect(), limits = (nothing, (-5, nothing))))

save(joinpath(output, "figure2.png"), fig)

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

save(joinpath(output, "figure3.png"), fig)


################################################## Figure 4

df = subset(runs, :light => ByRow(==("shift")))

my_renamer = uppercasefirst ∘ string









sdjfdshflsjhfshdflhdlshflkdsa












df2 = flatten(df1, [:θs, :r])
df3 = combine(groupby(df2, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
df3.condition .= categorical(df3.condition; levels = ["remain", "dark", "dark induced"], ordered = true)

data(df3) * mapping(:r, :mean_resultant_vector, color = :condition => my_renamer) * visual(Lines) |> draw


CSV.write("data.csv", df3)

df4 = CSV.read("fitted mean_resultant_vector.csv", DataFrame)

data(df4) * mapping(:r, Symbol("q_0.5"), lower = Symbol("q_0.025"), upper = Symbol("q_0.975"), color = :condition => my_renamer) * visual(LinesFill) |> draw

fig = Figure()
ax = Axis(fig[1,1])
band!(ax, df4, lower, upper, alpha = 0.3)


rl = range(extrema(df3.r)..., 100)
newdf = DataFrame(r = repeat(rl, outer = 3), condition = repeat(levels(df3.condition), inner = 100), mean_resultant_vector = zeros(300))
fm = @formula(mean_resultant_vector ~ r*condition)
function bootstrap(df)
    # for i in 1:3
    try
        m = BetaRegression.fit(BetaRegressionModel, fm, df[sample(1:nrow(df), nrow(df)), :])
        predict(m, newdf)
    catch
        missings(nrow(newdf))
    end
    # end
end
c = stack(bootstrap(df3) for _ in 1:10_000)
y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
newdf.lower .= getindex.(y, 1)
newdf.mean_resultant_vector .= getindex.(y, 2)
newdf.upper .= getindex.(y, 3)

data(newdf) * mapping(:r, :mean_resultant_vector, lower = :lower, upper = :upper, color = :condition => my_renamer) * visual(LinesFill) |> draw





fig = Figure()
ax = Axis(fig[1,1])
for l in c
    lines!(ax, rl, l[1:100])
end




fm = @formula(mean_resultant_vector ~ r*condition)
m = BetaRegression.fit(BetaRegressionModel, fm, df3)

rl = range(extrema(df3.r)..., 100)
newdf = DataFrame(r = repeat(rl, outer = 3), condition = repeat(levels(df3.condition), inner = 100), mean_resultant_vector = zeros(300))
newdf.mean_resultant_vector = predict(m, newdf)

data(newdf) * mapping(:r, :mean_resultant_vector, color = :condition => my_renamer) * visual(Lines) |> draw







b, a... = coef(m)

yl = BetaRegression.linkinv.(LogitLink(), b .+ a[1] .* xl)
ci = confint(m)
cib1, cib2 = ci[1, :]
cia1, cia2 = ci[2, :]
lower = BetaRegression.linkinv.(LogitLink(), cia1 .* xl .+ cib1)
upper = BetaRegression.linkinv.(LogitLink(), cia2 .* xl .+ cib2)




sort(combine(groupby(df2, [:condition, :r]), :θs => Ref => :θs), [:condition, :r])


α = 0.05
n = 10_000
nsamples = 50
fm = @formula(mean_resultant_vector ~ r)


function sample_mrvl(θs)
    y = mean_resultant_vector.(eachrow(stack(sample(θs, nsamples))))
    df = DataFrame(r = r, mean_resultant_vector = y)
    m = BetaRegression.fit(BetaRegressionModel, fm, df)
    return coef(m)
end

function bootstrap(θs)
    ba = [zeros(2) for _ in 1:n]
    Threads.@threads for i in 1:n
        for j in 1:3
            try
                ba[i] = sample_mrvl(θs)
                break
            catch ex
                @warn "failed to converge"
                if j == 3
                    @error "failed 3 times"
                end
            end
        end
    end
    μb, μa = mean(ba)
    σb, σa = std(ba)
    cia1, cia2 = μa .+ 1.96σa .* [-1, 1]
    cib1, cib2 = μb .+ 1.96σb .* [-1, 1]

    r = range(0, l, 100)
    ci1 = BetaRegression.linkinv.(LogitLink(), cia1 .* r .+ cib1)
    μ = BetaRegression.linkinv.(LogitLink(), μa .* r .+ μb)
    ci2 = BetaRegression.linkinv.(LogitLink(), cia2 .* r .+ cib2)

    return (; r, ci1, μ, ci2)
end

df2 = combine(groupby(df1, :condition), :θs => bootstrap => [:r, :ci1, :μ, :ci2])


data(df2) * mapping(:r, :μ, lower = :ci1, upper = :ci2, color = :condition) * visual(LinesFill) |> draw()





newdf = DataFrame(x = xl, y = zeros(length(xl)))
newdf.y = predict(m, newdf)
disallowmissing!(newdf)
scatter(x, y)
band!(xl, lower, upper, alpha = 0.3)
lines!(xl, y2)
lines!(xl, newdf.y)



















lsdjkfhlsdjfhlsdfjhlasdfhslf


d = Normal(1,3)
h = rand(d, 1_000_000)
sample_sizes = [10, 100, 1000, 10_000]
σ = map(sample_sizes) do sample_size
    mean((std(sample(h, sample_size)) for _ in 1:10_000))
end

lines(sample_sizes, σ, axis = (;yscale = log, xscale = log))

m = [mean(sample(h, 1000)) for _ in 1:10_000]
var(m)






transform!(runs, [:t, :spl, :poi_index] => ByRow(get_center_rotate) => :center_rotate)













# here
transform!(runs, [:t, :spl] => ByRow(smooth_center) => :center_smooth_xy)
transform!(runs, [:l, :center_smooth_xy] => ByRow((l, xy) -> norm.(xy) ./ l) => :straightness)
transform!(runs, [:t, :spl, :center_rotate] => ByRow((t, spl, tform) -> tform.(SVector{2, Float64}.(spl.(t)))) => :sxy)
# transform!(runs, :sxy => ByRow(xy -> cropto(xy, 50)) => :sxy)
transform!(runs, [:t, :spl, :poi_index] => ByRow(get_turn_profile) => [:tθ, :θ])
transform!(runs, [:t, :spl] => ByRow(get_turn_profile) => :Θ)
transform!(runs, [:tθ, :θ] => ByRow(fit_logistic) => :ks)
transform!(runs, [:tθ, :ks] => ByRow((x, p) -> logistic.(x, p[2], p[1], 0)) => :θs, :ks => ByRow(first) => :k)
# transform!(runs, :spontaneous_end => ByRow(!ismissing) => :spontaneous, :runs_dance => ByRow(!ismissing) => :induced)
# transform!(runs, [:induced, :spontaneous] => ByRow((e, s) -> e ? "induced" : s ? "spontaneous" : "no") => :danced)
transform!(runs, [:t, :spl, :poi_index] => ByRow(get_mean_speed) => :speed)

# runs.danced .= categorical(runs.danced; levels = ["no", "spontaneous", "induced"], ordered = true)


######################### Figure 2


df1 = deepcopy(df)
df2 = flatten(df1, [:l, :Θ])

plt = data(df2) * mapping(:l => "Path length (cm)", :Θ => rad2deg => "θ", group = :run_id => nonnumeric, col = :condition => my_renamer) * visual(Lines)
fig = draw(plt, axis = (;  yticks = -90:90:2*360, ytickformat = "{:n}°", limits = ((0, 60), (-90, 2*360))))

save(joinpath(output, "figure2.png"), fig)


function getCI(x; α = 0.05)
    d = Truncated(Distributions.fit(Normal, x), 0, 1)
    c1, μ, c2 = quantile(d, [α/2, 0.5, 1 - α/2])
    (; c1, μ, c2)
end

df1 = subset(df, :light => ByRow(==("remain")))
df2 = flatten(df1, [:l, :straightness])
subset!(df2, :straightness => ByRow(!isnan), :l => ByRow(<(60)))
transform!(df2, :l => ByRow(x -> round(Int, x)) => :L)
df3 = combine(groupby(df2, :L), :straightness => getCI => [:c1, :μ, :c2])
sort!(df3, :L)
# subset!(df2, :L => ByRow(x -> 10 ≤ x ≤ 50))

L = copy(df3.L)
c1 = copy(df3.c1)
μ = copy(df3.μ)
c2 = copy(df3.c2)

const K = 11

function findi2(straightness)
    k = copy(K)
    for i in k:51
        if i < length(straightness) && c1[i] < straightness[i] < c2[i]
            k += 1
        else
            return k
        end
    end
    return k
end

crop(straightness) = straightness[1:findi2(straightness)]

df1 = subset(df, :light => ByRow(==("dark")))
df2 = flatten(df1, [:l, :straightness])
# subset!(df2, :straightness => ByRow(!isnan), :l => ByRow(<(60)))
transform!(df2, :l => ByRow(x -> round(Int, x)) => :L)
df3 = combine(groupby(df2, [:run_id, :L]), :straightness => mean => :straightness)
sort!(df3, [:run_id, :L])
df4 = combine(groupby(df3, :run_id), :L => Ref => :L, :straightness => Ref => :straightness)
transform!(df4, :straightness => ByRow(findi2) => :k)
transform!(df4, [:L, :k] => ByRow((x, i) -> x[K:i]) => :L)
transform!(df4, [:straightness, :k] => ByRow((x, i) -> x[K:i]) => :straightness)
subset!(df4, :k => ByRow(>(K + 2)))

fig = Figure()
ax = Axis(fig[1,1], xlabel = "Path langth (cm)", ylabel = "Straightness", limits = ((10, 50), (0.7, 1)))
band!(ax, L, c1, c2, color = :black,  alpha = 0.5)
lines!(ax, L, μ, color = :black)
for row in eachrow(df4)
    lines!(ax, row.L, row.straightness, color = :red)
end

save(joinpath(output, "figure2b.png"), fig)


df1 = subset(df, :light => ByRow(∈(("dark", "remain"))))
transform!(groupby(df1, :light), eachindex => :group_count)
subset!(df1, :group_count => ByRow(<(41)))
df2 = flatten(df1, [:l, :straightness])
subset!(df2, :straightness => ByRow(!isnan), :l => ByRow(<(60)))
transform!(df2, :l => ByRow(x -> round(Int, x)) => :L)
df3 = combine(groupby(df2, [:light, :L]), :straightness => getCI => [:c1, :μ, :c2])

plt = data(df3) * mapping(:L => "Path length (cm)", :μ, lower = :c1, upper = :c2, color = :light) * visual(LinesFill)
fig = draw(plt; axis = (; limits = ((0, 50), (0, 1))))
save(joinpath(output, "marie.png"), fig)


##############################

path = "tracks"
if isdir(path)
    rm(path, recursive=true)
end
mkpath(path)
CairoMakie.activate!()
@tasks for row in eachrow(runs)
    fig = plotone(row.run_id, row.xy, row.poi_index, row.center_rotate, row.sxy, row.t, row.poi, row.spontaneous_end, row.l, row.Θ, row.condition)
    save(joinpath(path, string(row.run_id, ".png")), fig)
end
GLMakie.activate!()

