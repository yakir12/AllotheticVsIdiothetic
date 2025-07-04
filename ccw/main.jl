using Statistics
using CSV, DataFrames
using GLMakie, AlgebraOfGraphics
# using MixedModels, GLM
using Random
using Optim

rm.(filter(==(".png") ∘ last ∘ splitext, readdir()))

function dance(start, stop, cw, fullturns)
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

function dance(start, stop, lastcw)
    if start < stop
        dance(start, stop, false, 0)
    elseif start > stop
        dance(start, stop, true, 0)
    else
        if lastcw
            dance(start, stop - 0.01, lastcw, 0)
        else
            dance(start, stop + 0.01, lastcw, 0)
        end
    end
end

# df = CSV.read("data.csv", DataFrame)

# function calc_turn(α, β)
#     a = [sincosd(α)...]
#     b = [sincosd(β)...]
#     if sign(a' * b) < 0
#         "ccw"
#     else
#         "cw"
#     end
# end


# function calc_turn(init, turn)
#     Δ = mod(turn - init, 360)
#     if Δ > 180
#         (Δ - 180, "ccw")
#     else
#         (Δ, "cw")
#     end
# end

# norm2init(init, turn) = rem2pi(turn - init, RoundNearest)

# function total_rotation(placed, direction, full, exited)
#     Δ = direction == "ccw" ? placed - exited : exited - placed
#     Δ += Δ < 0 ? 2pi : 0
#     Δ += full*2pi
#     Δ = direction == "cw" ? Δ : -Δ
#     return Δ
# end

df = CSV.read("data.csv", DataFrame, select = ["individual_number",
                                               "placed from angle (degrees)",
                                               "rotation 1 direction",
                                               "rotation category measured",
                                               "exit angle (degrees)",
                                               "go down angle (degrees)",
                                               "full lap", 
                                               "total absolute degrees of rotation"])

# transform!(groupby(df, :individual_number), :individual_number => (x -> 1:length(x)) => :n)

subset!(df, "rotation category measured" => ByRow(∈(("cw", "ccw"))))
@assert all(df[!, "rotation 1 direction"] .== df[!, "rotation category measured"])
select!(df, Not("rotation category measured"))

transform!(df, ["placed from angle (degrees)",
                "exit angle (degrees)",
                "go down angle (degrees)",
                "total absolute degrees of rotation"] .=> ByRow(deg2rad) .=> 
           ["placed from angle",
            "exit angle",
            "go down angle",
            "total absolute of rotation"])
select!(df, Not("placed from angle (degrees)",
                "exit angle (degrees)",
                "go down angle (degrees)",
                "total absolute degrees of rotation"))

transform!(df, "rotation 1 direction" => ByRow(==("cw")) => :cw)

transform!(df, ["placed from angle", "go down angle", "cw", "full lap"] => ByRow(dance) => :placed2down)
rename!(df, "placed from angle" => :start)

transform!(df, ["go down angle", "exit angle", "cw"] => ByRow(dance) => :down2exit)
# rename!(df, "go down angle" => :down)
rename!(df, "exit angle" => :exit)

select!(df, Not("go down angle", "full lap"))

# @assert all(

# for row in eachrow(df)
#     θ1 = row.start + row.placed2down
#     θ2 = row.down
#     if !isapprox(rem2pi(θ1, RoundNearest), rem2pi(θ2, RoundNearest), atol = 1e-10)
#         @show rem2pi(θ1, RoundNearest), rem2pi(θ2, RoundNearest)
#     end
# end


# transform!(df, [:start, :μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)
# select!(df, Not(:μ, :exit))


# transform!(df, :start => ByRow(x -> sign(sin(x)) > 0) => :placed_from_left)

# transform!(df, :placed2down => ByRow(x -> sum(abs, diff(x))) => :dance)
# transform!(df, :down2exit => ByRow(x -> sum(abs, diff(x))) => :work)
# transform!(df, :down2exit => ByRow(x -> abs(x[end])) => :deviation)


# display(fig)

# save("centered to mean exit.png", fig)

function to_minimize(start, placed2down, μ)
    start2 = rem2pi.(start .- μ, RoundNearest)
    r = [first(findmin(i -> abs(y - (i*2π - x)), -2:2)) for (x, y) in zip(start2, placed2down)]
    sum(r)
end

angular_mean(θs) = angle(sum(exp, 1im*θs))

transform!(groupby(df, :individual_number), :exit => angular_mean => :mean_exit, [:start, :placed2down] => ((s, p) -> angular_mean(s .+ p)) => :mean_down, [:start, :placed2down] => ((start, placed2down) -> optimize(μ -> to_minimize(start, placed2down, μ), -pi, pi).minimizer) => :optimal)

fig = data(df) * (mapping(:mean_exit => rad2deg => "Mean exit", :mean_down => rad2deg) * visual(Scatter; label = "Mean down") + mapping(:mean_exit => rad2deg => "Mean exit", :optimal => rad2deg) * visual(Scatter; color = :red, label = "Mean optimal") + pregrouped([0], [1]) * visual(ABLines; color = :gray, label = "y = x")) |> draw(; axis = (; xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -180:180:180, ytickformat = "{:n}°", aspect = DataAspect(), width = 200, limits = ((-190, 190),(-190, 190))))
resize_to_layout!(fig)

save("means.png", fig)


for μ in (:mean_exit, :mean_down, :optimal)
    df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)

    fig = (data(df2) * mapping(:start => rad2deg =>  "Placed", :placed2down => rad2deg => "Danced") * visual(Scatter; alpha = 0.5, label = "data") + data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["longest", "longer", "shortest", "longer", "longest"])) * mapping(:a, :b, color = :color => sorter("shortest", "longer", "longest") => "") * visual(ABLines; label = "y = -x") + data(DataFrame(x = [-160, 160])) * mapping(:x) * visual(VLines; color = :gray, linestyle = :dash, label = "Error margins")) |> draw(; axis = (; xticks = [-160, 0, 160], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect(), width = 200))
    resize_to_layout!(fig)

    save("relationship $μ.png", fig)
end


transform!(df, [:start, :mean_down] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)

transform!(df, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
# transform!(df, [:placed_from_left, :cw] => ByRow((left, cw) -> string(left ? "left" : "right", " ", cw ? "cw" : "ccw")) => :quadrant)

# df2 = copy(df)
df2 = subset(df, :placed2down => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))

function fun(left, cw, start, placed2down)
    # @assert length(unique(left)) == 1
    # @assert length(unique(cw)) == 1
    if left[1]
        if cw[1]
            p = 0
        else
            p = 2pi
        end
        tf = @. placed2down < -start + p
    else
        if cw[1]
            p = -2pi
        else
            p = 0
        end
        tf = @. placed2down > -start + p
    end
    count(tf)
end

tbl = combine(groupby(df2, [:placed_from_left, :cw]), [:placed_from_left, :cw, :start, :placed2down] => fun => :shoot, nrow)

function fun(overshoot, placed_from_left, start, placed2down)
    intercept = overshoot ? 0 : placed_from_left ? 2pi : -2pi
    ŷ = intercept - start
    Δ = ŷ - placed2down 
    placed_from_left ? Δ : -Δ
end

transform!(df2, [:placed_from_left, :cw] => ByRow(==) => :overshoot)

transform!(df2, [:overshoot, :placed_from_left, :start, :placed2down] => ByRow(fun) => :residuals, :start => (s -> abs.(s)) => :x)

data(df2) * mapping(:overshoot, :residuals, color = :overshoot) * visual(BoxPlot) |> draw()

rad2deg(mean(df2.residuals))

# lmmod = lmm(@formula(residuals ~ 1 + overshoot + (1|individual_number)), df2)
#
# lmmod = lm(@formula(residuals ~ 1 + overshoot), df2)
#
# data(df3) * mapping(:x, :residuals) * visual(Scatter) + pregrouped([0], [1]) * visual(ABLines) + pregrouped([[c] for c in coef(lmod)]...) * visual(ABLines) |> draw()#; axis = (; aspect = DataAspect()))

# using Turing
#
#
# @model function bmodel(g, x, y)
#     σ ~ InverseGamma(2, 2)
#     slope ~ Normal(1, 3)
#     intercept ~ filldist(Normal(0, 3), length(unique(g)))
#     b = intercept[g]
#     for i in eachindex(x)
#         μ = b[i] + slope*x[i]
#         y[i] ~ Normal(μ, σ)
#     end
# end
#
# model = bmodel(df3.individual_number, df3.x, df3.residuals)
# chain = sample(model, NUTS(), MCMCThreads(), 100000, 4)
#
# # chain = sample(model, Gibbs(), 1000)
#
# describe(chain)
#
# fig = Figure()
# for (i, var_name) in enumerate(chain.name_map.parameters)
#     draw!(
#           fig[i, 1],
#           data(chain) *
#           mapping(var_name; color=:chain => nonnumeric) *
#           AlgebraOfGraphics.density() *
#           visual(fillalpha=0); axis = (; height = 200))
# end
# resize_to_layout!(fig)
#
# @model function bmodel(g, x, y)
#     σ ~ InverseGamma(2, 2)
#     intercept ~ filldist(Normal(0, 3), length(unique(g)))
#     b = intercept[g]
#     for i in eachindex(x)
#         μ = b[i] + x[i]
#         y[i] ~ Normal(μ, σ)
#     end
# end
#
# model = bmodel(df3.individual_number, df3.x, df3.residuals)
# chain = sample(model, NUTS(), MCMCThreads(), 100000, 4)
#
# # chain = sample(model, Gibbs(), 1000)
#
# describe(chain)
#
# fig = Figure()
# for (i, var_name) in enumerate(chain.name_map.parameters)
#     draw!(
#           fig[i, 1],
#           data(chain) *
#           mapping(var_name; color=:chain => nonnumeric) *
#           AlgebraOfGraphics.density() *
#           visual(fillalpha=0); axis = (; height = 200))
# end
# resize_to_layout!(fig)





# AlgebraOfGraphics.density(bandwidth=0.5) |> draw()


# df2 = transform(df, [:start, :cw] => ByRow((s, c) -> string(s < 0 ? "right" : "left", " ", c ? "cw" : "ccw")) => :quadrant)
# # subset!(df2, :placed2down => ByRow(<(2pi) ∘ abs), :start => ByRow(>(0.5) ∘ abs))
# @. model(x, p) = -x + p[2]*exp(x*p[1] + p[3])
# @. model(x, p) = -p[1]* x + p[2]
# df3 = combine(groupby(df2, :quadrant), [:start, :placed2down] => ((x, y) -> Ref(coef(curve_fit(model, x, y, ones(2))))) => :coef)
# transform!(df3, :quadrant => ByRow(q -> contains(q, "right") ? range(-pi, 0, 100) : range(0, pi, 100)) => :xl)
# transform!(df3, [:xl, :coef] => ByRow(model) => :yl)
# select!(df3, Not(:coef))
# df4 = flatten(df3, [:xl, :yl])
# xy = data(df2) * mapping(:start, :placed2down) * visual(Scatter)
# line = data(df4) * mapping(:xl, :yl, group = :quadrant) * visual(Lines; color = :red)
# draw(xy + line)


fig = Figure()
subgl_left = GridLayout()
subgl_right = GridLayout()
fig.layout[1, 1] = subgl_left
fig.layout[1, 2] = subgl_right
xy = data(df) * mapping(:start => rad2deg => "Placed", :placed2down => rad2deg => "Danced") * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["extra loop", "longer", "shorter", "longer", "extra loop"])) * mapping(:a, :b, color = :color => sorter("shorter", "longer", "extra loop") => "") * visual(ABLines)
vline = data(DataFrame(geometry = [Rect(-160, -360, 140, 720), Rect(20, -360, 140, 720)])) * mapping(:geometry) * visual(Poly; color = (:gray, 0.2), label = "included")
toplot = vline + abline + xy
g = draw!(subgl_left[1,1], toplot; axis = (; xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect(), width = 200))
legend!(subgl_right[1,1], g)
draw!(subgl_right[2,1], data(df2) * mapping(:residuals => rad2deg) * visual(Hist); axis = (; ylabel = "#", xticks = [25, 180], xtickformat = "{:n}°", width = 200))
resize_to_layout!(fig)


save("relationship.png", fig)

# using LsqFit


# lines!(xdata, model(xdata, coef(fit)))

# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# hexbin!(ax, df.start, df.placed2down, cellsize = (1, 1), threshold = 0, colormap = [Makie.to_color(:transparent); Makie.to_colormap(:viridis)],
#     strokewidth = 0.5,
#     strokecolor = :gray50,
#     colorscale = Makie.pseudolog10)
#
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect(), xticks = -180:90:180, yticks = -720:180:720)
# datashader!(ax, Point2f.(rad2deg.(df.start), rad2deg.(df.placed2down)), binsize = 10, colormap=[:transparent, :black], async = false)
# ablines!(360 .* (-2:2), -1, label = ["longest", "longer", "shortest", "longer", "longest"])

# df2 = combine(groupby(df, [:start, :placed2down]), nrow)
# tricontourf(df2.start, df2.placed2down, df2.nrow)

# df2 = combine(groupby(df, [:start, :placed2down]), nrow)
# heatmap(df2.start, df2.placed2down, df2.nrow)


# transform!(df, [:start, :placed2down] => ByRow((x, y) -> findmin(i -> abs(y - (i*2π - x)), -2:2)) => [:r, :i])
#
# fig = (data(df) * mapping(:start => rad2deg =>  "Placed down (°)", :placed2down => rad2deg => "Danced (°)", layout = :individual_number => nonnumeric) * visual(Scatter; label = "data") + data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["longest", "longer", "shortest", "longer", "longest"])) * mapping(:a, :b, color = :color => sorter("shortest", "longer", "longest")) * visual(ABLines; label = "y = -x")) |> draw(; axis = (; aspect = DataAspect(), xticks = -180:180:180, yticks = -720:180:720, height = 200))
# resize_to_layout!(fig)
# save("relationship2.png", fig)


# fig = (data(df) * mapping(:start => rad2deg =>  "Placed down (°)", :placed2down => rad2deg => "Danced (°)", layout = :individual_number => nonnumeric) * visual(Scatter; label = "data") + data(DataFrame(a = 360 .* (-2:2), b = -1)) * mapping(:a, :b) * visual(ABLines; color = :red, label = "y = -x")) |> draw(; axis = (; xticks = -180:180:180, yticks = -720:180:720, aspect = DataAspect(), width = 200))


shuffle!(df)
sort!(df, [:placed_from_left, :cw, :start])

fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df, :individual_number)))
    ij = CartesianIndices((5, 6))[i]
    ax = Axis(fig[Tuple(ij)...], aspect = DataAspect(), height = 200, width = 200)
    for (j, row) in enumerate(eachrow(g))
        radius = j + 1
        poly!(ax, Makie.GeometryBasics.Polygon(Circle(zero(Point2f), radius+0.5), [Circle(zero(Point2f), radius-0.5)]), color = row.placed_from_left ? :blue : :green, alpha = 0.25)
        color = (; color = row.cw ? :blue : :green)
        ts = row.start .+ range(0, row.placed2down, length = 100)
        lines!(ax, radius*Point2f.(reverse.(sincos.(ts .+ π/2))); color...)
        α = row.start + row.placed2down
        α += row.cw ? -π/2 : π/2
        scatter!(ax, radius*Point2f(reverse(sincos(row.start + π/2))); color..., marker = '|', markersize=10, rotation = row.start + pi/2 + π/2)
        scatter!(ax, radius*Point2f(reverse(sincos(row.start + row.placed2down + π/2))); color..., marker = :utriangle, markersize=10, rotation=α)
    end
    text!(ax, 0, 0; text = string(k.individual_number), align = (:center, :center))
    hidedecorations!(ax)
    hidespines!(ax)
end
resize_to_layout!(fig)

save("centered to mean go down.png", fig)



wejrhlwehrwkehrlwhrwlwjeh


below is the bayesian stuff

# df2 = combine(groupby(df, :individual_number), [:placed_from_left, :cw] => ((x1, x2) -> round(Int, 100count(x1 .≠ x2)/length(x1))) => :longest, [:placed_from_left, :cw] => ((x1, x2) -> round(Int, 100count(x1 .== x2)/length(x1))) => :shortest, :cw => (x -> round(Int, 100count(!, x)/length(x))) => :ccw, :cw => (x -> round(Int, 100count(x)/length(x))) => :cw, :placed2down => (x -> round(Int, 100count(x -> abs(x) < π, x)/length(x))) => :shortest_dance)
#
# mean.(eachcol(df2))

place2weight(x) = 1 - abs((abs(x) - π/2)/π/2)
df2 = select(df, :individual_number, :start, :cw)
# transform!(groupby(df2, :individual_number), groupindices => :id)

# @assert all(≠(0), df2.side)
# @assert all(x -> 0 ≤ x ≤ 1, df2.w)


n = 30
handedness = rand(n)
shortest = 0.99
nreps = 10
df3 = DataFrame(start = rand(-π..π, nreps*n), individual_number = repeat(1:n, inner = nreps))
function generate(start, individual_number)
    w = place2weight(start)
    s = sign(start)
    p2 = (1 + s*(2shortest - 1))/2
    z = handedness[individual_number] + w*p2
    p = 1 / (1 + exp(-z))
    d = Bernoulli(p)
    rand(d)
end
transform!(df3, [:start, :individual_number] => ByRow(generate) => :cw)
hist(df3.cw)

using Turing

@model function bmodel(starts, cws)
    # handedness ~ Beta(1, 1)
    # shortest ~ Beta(1, 1)
    handedness ~ Uniform(0, 1)
    shortest ~ Uniform(0, 1)
    ratio ~ Beta(1, 1)
    for i in eachindex(starts)
        w = place2weight(starts[i])
        s = sign(starts[i])
        p2 = (1 + s*(2shortest - 1))/2
        z = ratio*handedness + (1 - ratio)*w*p2
        cws[i] ~ BernoulliLogit(z)
    end
end

model = bmodel(df2.start, df2.cw)
chain = sample(model, NUTS(), MCMCThreads(), 100000, 4)

# chain = sample(model, Gibbs(), 1000)

describe(chain)

fig = Figure()
for (i, var_name) in enumerate(chain.name_map.parameters)
    draw!(
          fig[i, 1],
          data(chain) *
          mapping(var_name; color=:chain => nonnumeric) *
          AlgebraOfGraphics.density() *
          visual(fillalpha=0)
          ; axis = (; limits = ((0, 1), nothing)))
end
# resize_to_layout!(fig)



@model function bmodel(individual_number, starts, cws)
    handedness ~ filldist(Uniform(0, 1), maximum(individual_number))
    p1 = handedness[individual_number]
    shortest ~ Uniform(0, 1)
    ratio ~ filldist(Beta(1, 1), maximum(individual_number))
    r = ratio[individual_number]
    # ratio ~ Beta(1, 1)
    for i in eachindex(starts)
        w = place2weight(starts[i])
        s = sign(starts[i])
        p2 = (1 + s*(2shortest - 1))/2
        # z = ratio*p1[i] + (1 - ratio)*w*p2
        z = r[i]*p1[i] + (1 - r[i])*w*p2
        cws[i] ~ BernoulliLogit(z)
    end
end

model = bmodel(df2.individual_number, df2.start, df2.cw)
chain = sample(model, NUTS(), MCMCThreads(), 10000, 4)

# chain = sample(model, Gibbs(), 1000)

describe(chain)

fig = Figure()
for (i, var_name) in enumerate(chain.name_map.parameters)
    draw!(
          fig[i, 1],
          data(chain) *
          mapping(var_name; color=:chain => nonnumeric) *
          AlgebraOfGraphics.density() *
          visual(fillalpha=0)
          ; axis = (; width = 100, height = 100))
end
resize_to_layout!(fig)


fjhsdkfhlsakhfs

transform!(groupby(df, :individual_number), :placed2down => (θs -> angle(sum(exp.(1im*last(θs))))) => :μ)

transform!(df, [:μ, :placed2down] => ByRow((x, xs) -> xs .- x) => :placed2down)
transform!(df, [:μ, :down2exit] => ByRow((x, xs) -> xs .- x) => :down2exit)
select!(df, Not(:μ))

transform!(df, :placed2down => ByRow(x -> sign(sin(x[1])) > 0) => :placed_from_left)

shuffle!(df)
sort!(df, [:placed_from_left, :cw, order(:placed2down, by = x -> first(x) + 10000)])

fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df, :individual_number)))
    ij = CartesianIndices((5, 6))[i]
    ax = Axis(fig[Tuple(ij)...], aspect = DataAspect(), height = 200, width = 200)
    for (j, row) in enumerate(eachrow(g))
        radius = j + 1
        poly!(ax, Makie.GeometryBasics.Polygon(Circle(zero(Point2f), radius+0.5), [Circle(zero(Point2f), radius-0.5)]), color = row.placed_from_left ? :blue : :green, alpha = 0.25)
        color = (; color = row.cw ? :blue : :green)
        lines!(ax, radius*Point2f.(reverse.(sincos.(row.placed2down .+ π/2))); color...)
        α = row.placed2down[1]
        α += row.cw ? π : 0
        scatter!(ax, radius*Point2f(reverse(sincos(row.placed2down[1] + π/2))); color..., marker = :utriangle, markersize=10, rotation = α + π/2)
        color = (; color = :red)
        scatter!(ax, radius*Point2f(reverse(sincos(row.down2exit[end] + π/2))); color..., marker = '|', markersize=10, rotation=row.down2exit[end] + pi/2 + π/2)
        # text!(ax, radius, 0; text = string(row.n), align = (:center, :center))
    end
    text!(ax, 0, 0; text = string(k.individual_number), align = (:center, :center))
    hidedecorations!(ax)
    hidespines!(ax)
end
resize_to_layout!(fig)

# display(fig)

save("centered to mean down from ball.png", fig)



transform!(df, :placed2down => ByRow(xs -> xs .- last(xs)) => :placed2down)
transform!(df, :down2exit => ByRow(xs -> xs .- first(xs)) => :down2exit)
transform!(df, :placed2down => ByRow(x -> sign(sin(x[1])) > 0) => :placed_from_left)

shuffle!(df)
sort!(df, [:placed_from_left, :cw, order(:placed2down, by = x -> first(x) + 10000)])

fig = Figure(size = (2000, 2000))
ax = Axis(fig[1,1], aspect = DataAspect())
for (j, row) in enumerate(eachrow(df))
    radius = j + 1
    poly!(ax, Makie.GeometryBasics.Polygon(Circle(zero(Point2f), radius+0.5), [Circle(zero(Point2f), radius-0.5)]), color = row.placed_from_left ? :blue : :green, alpha = 0.25)
    color = (; color = row.cw ? :blue : :green)
    lines!(ax, radius*Point2f.(reverse.(sincos.(row.placed2down .+ π/2))); color...)
    α = row.placed2down[1]
    α += row.cw ? π : 0
    scatter!(ax, radius*Point2f(reverse(sincos(row.placed2down[1] + π/2))); color..., marker = :utriangle, markersize=10, rotation = α + π/2)
    color = (; color = :red)
    scatter!(ax, radius*Point2f(reverse(sincos(row.down2exit[end] + π/2))); color..., marker = '|', markersize=10, rotation=row.down2exit[end] + pi/2 + π/2)
end
hidedecorations!(ax)
hidespines!(ax)

save("aligned to down from ball.png", fig)









sdjjfhgksdhfgskdhfgskdfg




fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df, :individual_number)))
    # i, (k, g) = first(enumerate(pairs(groupby(df, :individual_number))))
    ij = CartesianIndices((5, 6))[i]
    ax = PolarAxis(fig[Tuple(ij)...], rgridvisible = false, spinevisible = false, theta_0 = pi/2, height = 200, width = 200)#, thetaticks = (range(0, length = 4, step = pi/2), string.(range(0, length = 4, step = 90)) .* "°"))
    for (j, row) in enumerate(eachrow(g))
        r = fill(j, length(row.placed2down))
        lines!(ax, range(0, 2pi, 100), r, color = row.placed_from_left ? :blue : :green, linewidth = 5, alpha = 0.25)
        color = (; color = row.cw ? :blue : :green)
        lines!(ax, row.placed2down, r; color...)
        scatter!(ax, row.placed2down[1], j; color..., marker = '|', markersize=10, rotation=row.placed2down[1])
        # θ = min(row.placed2down[1] + 0.2, row.placed2down[end])
        ls = range(j*row.placed2down[1], j*row.placed2down[end], step = (pi)*sign(step(row.placed2down)))[2:end-1]
        if isempty(ls)
            ls = range(j*row.placed2down[50], j*row.placed2down[50])
        end
        for l in ls[1:1]
            θ = l / j
            α = θ + π/2
            α += row.cw ? pi : 0
            scatter!(ax, θ, j; color..., marker = :utriangle, markersize=10, rotation = α)
        end
        color = (; color = :red)
        # lines!(ax, row.down2exit, r; color...)
        scatter!(ax, row.down2exit[end], j; color..., marker = '|', markersize=10, rotation=row.down2exit[end])
    end
    hidedecorations!(ax)
end
resize_to_layout!(fig)

display(fig)


i, (k, g) = first(enumerate(pairs(groupby(df, :individual_number))))
(j, row) = first(enumerate(eachrow(g)))
fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
lines!(Circle(zero(Point2f), j), color = row.placed_from_left ? :blue : :green, linewidth = 5, alpha = 0.25)
color = (; color = row.cw ? :blue : :green)
lines!(ax, j*Point2f.(reverse.(sincos.(row.placed2down))); color...)
scatter!(ax, j*Point2f(reverse(sincos(row.placed2down[1]))); color..., marker = '|', markersize=10, rotation=row.placed2down[1] + pi/2)
ls = range(j*row.placed2down[1], j*row.placed2down[end], step = (pi)*sign(step(row.placed2down)))
if length(ls) < 2
    l = j*row.placed2down[50]
else
    l = ls[2]
end
θ = l / j
α = θ
α += row.cw ? pi : 0
scatter!(ax, j*Point2f(reverse(sincos(θ))); color..., marker = :utriangle, markersize=10, rotation = α)
color = (; color = :red)
scatter!(ax, j*Point2f(reverse(sincos(row.down2exit[end]))); color..., marker = '|', markersize=10, rotation=row.down2exit[end] + pi/2)


lkhfdshflshfsh




transform!(df, :placed => ByRow(x -> x < 0 ? "ccw" : "cw") => "placed from")

transform!(df, ["placed", "rotation 1 direction", "full lap", "descended"] => ByRow(total_rotation) => :rotated2descent)
transform!(df, ["placed", "rotation 1 direction", "full lap", "exited"] => ByRow(total_rotation) => :rotated2exit)

# @assert all(df."total absolute degrees of rotation" .== round.(Int, rad2deg.(abs.(df.rotated))))
select!(df, Not("total absolute degrees of rotation"))

# for (i, row) in enumerate(eachrow(df))
#     if !isapprox(rem2pi(row.placed - row.rotated, RoundNearest),row.exited, atol = 1e-14)
#         @show i
#     end
# end

# @assert all(isapprox.(rem2pi.(df.placed .- df.rotated, RoundNearest), df.exited, atol = 1e-14))

# transform!(df, ["initial exit", "placed from angle (degrees)"] => ByRow(calc_turn) => ["absolute placed from angle (degrees)", "calculated placed from"])

# transform!(df, :w => x -> x/sum(x), renamecols = false)
# This fails because the data is wrong!
# @assert all(df[!, "placed from"] .== df[!, "calculated placed from"])
# df.rownumber .= 1:nrow(df)
# for row in eachrow(df)
#     if row."placed from" ≠ row."calculated placed from"
#         @show row.rownumber
#     end
# end

# transform!(groupby(df, ["individual_number", "rotation 1 direction"]), x -> 1:nrow(x))
# rename!(df, :x1 => :r)



function get_trs(θ₁, rotated, r)
    a = 0.2
    θ₂ = θ₁ + rotated
    θ = range(θ₁, θ₂, 100)
    Point2f.(θ, r .+ a .* θ)
end
fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df, :individual_number)))
    # i, (k, g) = first(enumerate(pairs(groupby(df, :individual_number))))
    ij = CartesianIndices((5, 6))[i]
    ax = PolarAxis(fig[Tuple(ij)...], rgridvisible = false, theta_0 = pi/2, height = 200, width = 200)#, thetaticks = (range(0, length = 4, step = pi/2), string.(range(0, length = 4, step = 90)) .* "°"))
    for (j, row) in enumerate(eachrow(g))
        trs = get_trs(row.placed, row.rotated2exit, 2j)
        # trs = get_trs(row.placed, row.rotated, 2j)
        lines!(ax, trs, color = :black)
        scatter!(ax, trs[end], color = :black)
        trs = get_trs(row.placed, row.rotated2descent, 2j)
        scatter!(ax, trs[end], color = :gray)
    end
    hidedecorations!(ax)
end
resize_to_layout!(fig)

display(fig)


sjdhfsdhjflsdhlfsdkjh

# select!(df, Not("calculated placed from", "initial exit", "placed from angle (degrees)", "absolute placed from angle (degrees)"))

rename!(df, Dict("rotation 1 direction" => "y", "placed from" => "x", "individual_number" => "id"))
# transform!(df, :id => ByRow(string), renamecols = false)
transform!(df, [:x, :y] .=> ByRow(==("cw")), renamecols = false) # equal to cw






transform!(df, :placed => ByRow(x -> (pi - abs(x))/pi) => "w")




# using MixedModels

df = CSV.read("repeated_direction.csv", DataFrame, select = ["individual_number", # (ID)
                                                             # "run_column", # (trial number)
                                                             "rotation category measured", # (rotationen som utfördes: "cw", "ccw", "both" or "no")
                                                             "placed from", # (hur dyngbaggen placerades i relation till dess exit angle)
                                                             # "smallest possible direction", # (motsatt riktning till "placed from" - denna riktning som du vill jämföra mot "rotation category measured")
                                                             # "Minimizing yaw rotation?" # (svarar på om dyngbaggen tog minst möjliga rotation under denna trial: 1/0) (edited)
                                                            ])
rename!(df, Dict("rotation category measured" => "y", "placed from" => "x", "individual_number" => "id"))
subset!(df, :y => ByRow(∈(("cw", "ccw"))))
transform!(df, :id => ByRow(string), renamecols = false)
transform!(df, [:x, :y] .=> ByRow(==("cw")), renamecols = false) # equal to cw

# transform!(df, ["placed from", "rotation category measured"] .=> ByRow(txt -> Float64(txt == "cw") + rand()/10) .=> [:x, :y])

data(df) * mapping(:x, :y, layout = :id) * visual(Scatter) |> draw()



function fun(g)
    fm = @formula(y ~ 1 + x)
    gm = glm(fm, g, Binomial(), wts = g.w/sum(g.w))
    newdf = DataFrame(x = range(0, 1, 100), y =  false)
    newdf = hcat(select(newdf, Not(:y)), predict(gm, newdf; interval = :confidence))
    tbl = coeftable(gm)
    p = tbl.cols[tbl.pvalcol]
    newdf.intercept .= p[1]
    newdf.slope .= p[2]
    newdf
end


df2 = combine(groupby(df, :id), fun)

df3 = combine(groupby(df2, :id), :intercept => first => :intercept, :slope => first => :slope)

data(df2) * mapping(:x => "Placed from", :prediction, lower = :lower, upper = :upper, layout = :id) * visual(LinesFill) |> draw(; axis = (height = 200, width = 200, ylabel = "Rotated", xticks = ([0, 1], ["ccw", "cw"]), yticks = ([0, 1], ["ccw", "cw"])))


fm = @formula(y ~ 1 + x + (1|id))
gm = fit(MixedModel, fm, df, Bernoulli())

newdf = DataFrame(x = ["ccw", "cw"], y =  ["ccw", "cw"], id = "NEW")

newdf.y .= predict(gm, newdf; new_re_levels = :population, type = :response)

n = 1000
newdf = DataFrame(x = repeat(["ccw", "cw"], outer = n), y = repeat(["ccw", "cw"], outer = n), id = "NEW")
simulate(MersenneTwister(42), gm, df)
