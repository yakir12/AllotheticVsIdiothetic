kakakakakakakaka

transform!(df, [:start, :mean_down] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)

transform!(df, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
# transform!(df, [:placed_from_left, :cw] => ByRow((left, cw) -> string(left ? "left" : "right", " ", cw ? "cw" : "ccw")) => :quadrant)

# df2 = copy(df)
df2 = subset(df, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))

function fun(left, cw, start, total_dance)
    # @assert length(unique(left)) == 1
    # @assert length(unique(cw)) == 1
    if left[1]
        if cw[1]
            p = 0
        else
            p = 2pi
        end
        tf = @. total_dance < -start + p
    else
        if cw[1]
            p = -2pi
        else
            p = 0
        end
        tf = @. total_dance > -start + p
    end
    count(tf)
end

tbl = combine(groupby(df2, [:placed_from_left, :cw]), [:placed_from_left, :cw, :start, :total_dance] => fun => :shoot, nrow)

function fun(overshoot, placed_from_left, start, total_dance)
    intercept = overshoot ? 0 : placed_from_left ? 2pi : -2pi
    ŷ = intercept - start
    Δ = ŷ - total_dance 
    placed_from_left ? Δ : -Δ
end

transform!(df2, [:placed_from_left, :cw] => ByRow(==) => :overshoot)

transform!(df2, [:overshoot, :placed_from_left, :start, :total_dance] => ByRow(fun) => :residuals, :start => (s -> abs.(s)) => :x)

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
# # subset!(df2, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(>(0.5) ∘ abs))
# @. model(x, p) = -x + p[2]*exp(x*p[1] + p[3])
# @. model(x, p) = -p[1]* x + p[2]
# df3 = combine(groupby(df2, :quadrant), [:start, :total_dance] => ((x, y) -> Ref(coef(curve_fit(model, x, y, ones(2))))) => :coef)
# transform!(df3, :quadrant => ByRow(q -> contains(q, "right") ? range(-pi, 0, 100) : range(0, pi, 100)) => :xl)
# transform!(df3, [:xl, :coef] => ByRow(model) => :yl)
# select!(df3, Not(:coef))
# df4 = flatten(df3, [:xl, :yl])
# xy = data(df2) * mapping(:start, :total_dance) * visual(Scatter)
# line = data(df4) * mapping(:xl, :yl, group = :quadrant) * visual(Lines; color = :red)
# draw(xy + line)


fig = Figure()
subgl_left = GridLayout()
subgl_right = GridLayout()
fig.layout[1, 1] = subgl_left
fig.layout[1, 2] = subgl_right
xy = data(df) * mapping(:start => rad2deg => "Placed", :total_dance => rad2deg => "Danced") * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["extra loop", "longer", "shorter", "longer", "extra loop"])) * mapping(:a, :b, color = :color => sorter("shorter", "longer", "extra loop") => "") * visual(ABLines)
vline = data(DataFrame(geometry = [Rect(-160, -360, 140, 720), Rect(20, -360, 140, 720)])) * mapping(:geometry) * visual(Poly; color = (:gray, 0.2), label = "included")
toplot = vline + abline + xy
g = draw!(subgl_left[1,1], toplot; axis = (; xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect(), width = 200))
legend!(subgl_right[1,1], g)
draw!(subgl_right[2,1], data(df2) * mapping(:residuals => rad2deg) * visual(Hist); axis = (; ylabel = "#", xticks = [25, 180], xtickformat = "{:n}°", width = 200))
resize_to_layout!(fig)


save("relationship.png", fig)


for μ in (:mean_exit, :mean_down, :optimal)
    df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)

    fig = Figure()
    subgl_left = GridLayout()
    subgl_right = GridLayout()
    fig.layout[1, 1] = subgl_left
    fig.layout[1, 2] = subgl_right
    xy = data(df2) * mapping(:start => rad2deg => "Placed", :total_dance => rad2deg => "Danced") * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
    abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["extra loop", "longer", "shorter", "longer", "extra loop"])) * mapping(:a, :b, color = :color => sorter("shorter", "longer", "extra loop") => "") * visual(ABLines)
    vline = data(DataFrame(geometry = [Rect(-160, -360, 140, 720), Rect(20, -360, 140, 720)])) * mapping(:geometry) * visual(Poly; color = (:gray, 0.2), label = "included")
    toplot = vline + abline + xy
    g = draw!(subgl_left[1,1], toplot; axis = (; xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect(), width = 200))
    legend!(subgl_right[1,1], g)


    df3 = subset(df, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))

    transform!(df3, [:placed_from_left, :cw] => ByRow(==) => :overshoot)

    transform!(df3, [:overshoot, :placed_from_left, :start, :total_dance] => ByRow(fun) => :residuals, :start => (s -> abs.(s)) => :x)

    draw!(subgl_right[2,1], data(df3) * mapping(:residuals => rad2deg) * visual(Hist); axis = (; ylabel = "#", xticks = [25, 180], xtickformat = "{:n}°", width = 200))
    resize_to_layout!(fig)

    save("relationship $μ.png", fig)
end




# using LsqFit


# lines!(xdata, model(xdata, coef(fit)))

# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect())
# hexbin!(ax, df.start, df.total_dance, cellsize = (1, 1), threshold = 0, colormap = [Makie.to_color(:transparent); Makie.to_colormap(:viridis)],
#     strokewidth = 0.5,
#     strokecolor = :gray50,
#     colorscale = Makie.pseudolog10)
#
# fig = Figure()
# ax = Axis(fig[1,1], aspect = DataAspect(), xticks = -180:90:180, yticks = -720:180:720)
# datashader!(ax, Point2f.(rad2deg.(df.start), rad2deg.(df.total_dance)), binsize = 10, colormap=[:transparent, :black], async = false)
# ablines!(360 .* (-2:2), -1, label = ["longest", "longer", "shortest", "longer", "longest"])

# df2 = combine(groupby(df, [:start, :total_dance]), nrow)
# tricontourf(df2.start, df2.total_dance, df2.nrow)

# df2 = combine(groupby(df, [:start, :total_dance]), nrow)
# heatmap(df2.start, df2.total_dance, df2.nrow)


# transform!(df, [:start, :total_dance] => ByRow((x, y) -> findmin(i -> abs(y - (i*2π - x)), -2:2)) => [:r, :i])
#
# fig = (data(df) * mapping(:start => rad2deg =>  "Placed down (°)", :total_dance => rad2deg => "Danced (°)", layout = :individual_number => nonnumeric) * visual(Scatter; label = "data") + data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["longest", "longer", "shortest", "longer", "longest"])) * mapping(:a, :b, color = :color => sorter("shortest", "longer", "longest")) * visual(ABLines; label = "y = -x")) |> draw(; axis = (; aspect = DataAspect(), xticks = -180:180:180, yticks = -720:180:720, height = 200))
# resize_to_layout!(fig)
# save("relationship2.png", fig)


# fig = (data(df) * mapping(:start => rad2deg =>  "Placed down (°)", :total_dance => rad2deg => "Danced (°)", layout = :individual_number => nonnumeric) * visual(Scatter; label = "data") + data(DataFrame(a = 360 .* (-2:2), b = -1)) * mapping(:a, :b) * visual(ABLines; color = :red, label = "y = -x")) |> draw(; axis = (; xticks = -180:180:180, yticks = -720:180:720, aspect = DataAspect(), width = 200))


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
        ts = row.start .+ range(0, row.total_dance, length = 100)
        lines!(ax, radius*Point2f.(reverse.(sincos.(ts .+ π/2))); color...)
        α = row.start + row.total_dance
        α += row.cw ? -π/2 : π/2
        scatter!(ax, radius*Point2f(reverse(sincos(row.start + π/2))); color..., marker = '|', markersize=10, rotation = row.start + pi/2 + π/2)
        scatter!(ax, radius*Point2f(reverse(sincos(row.start + row.total_dance + π/2))); color..., marker = :utriangle, markersize=10, rotation=α)
    end
    text!(ax, 0, 0; text = string(k.individual_number), align = (:center, :center))
    hidedecorations!(ax)
    hidespines!(ax)
end
resize_to_layout!(fig)

save("centered to mean go down.png", fig)



wejrhlwehrwkehrlwhrwlwjeh


below is the bayesian stuff

# df2 = combine(groupby(df, :individual_number), [:placed_from_left, :cw] => ((x1, x2) -> round(Int, 100count(x1 .≠ x2)/length(x1))) => :longest, [:placed_from_left, :cw] => ((x1, x2) -> round(Int, 100count(x1 .== x2)/length(x1))) => :shortest, :cw => (x -> round(Int, 100count(!, x)/length(x))) => :ccw, :cw => (x -> round(Int, 100count(x)/length(x))) => :cw, :total_dance => (x -> round(Int, 100count(x -> abs(x) < π, x)/length(x))) => :shortest_dance)
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

transform!(groupby(df, :individual_number), :total_dance => (θs -> angle(sum(exp.(1im*last(θs))))) => :μ)

transform!(df, [:μ, :total_dance] => ByRow((x, xs) -> xs .- x) => :total_dance)
transform!(df, [:μ, :down2exit] => ByRow((x, xs) -> xs .- x) => :down2exit)
select!(df, Not(:μ))

transform!(df, :total_dance => ByRow(x -> sign(sin(x[1])) > 0) => :placed_from_left)

shuffle!(df)
sort!(df, [:placed_from_left, :cw, order(:total_dance, by = x -> first(x) + 10000)])

fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df, :individual_number)))
    ij = CartesianIndices((5, 6))[i]
    ax = Axis(fig[Tuple(ij)...], aspect = DataAspect(), height = 200, width = 200)
    for (j, row) in enumerate(eachrow(g))
        radius = j + 1
        poly!(ax, Makie.GeometryBasics.Polygon(Circle(zero(Point2f), radius+0.5), [Circle(zero(Point2f), radius-0.5)]), color = row.placed_from_left ? :blue : :green, alpha = 0.25)
        color = (; color = row.cw ? :blue : :green)
        lines!(ax, radius*Point2f.(reverse.(sincos.(row.total_dance .+ π/2))); color...)
        α = row.total_dance[1]
        α += row.cw ? π : 0
        scatter!(ax, radius*Point2f(reverse(sincos(row.total_dance[1] + π/2))); color..., marker = :utriangle, markersize=10, rotation = α + π/2)
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



transform!(df, :total_dance => ByRow(xs -> xs .- last(xs)) => :total_dance)
transform!(df, :down2exit => ByRow(xs -> xs .- first(xs)) => :down2exit)
transform!(df, :total_dance => ByRow(x -> sign(sin(x[1])) > 0) => :placed_from_left)

shuffle!(df)
sort!(df, [:placed_from_left, :cw, order(:total_dance, by = x -> first(x) + 10000)])

fig = Figure(size = (2000, 2000))
ax = Axis(fig[1,1], aspect = DataAspect())
for (j, row) in enumerate(eachrow(df))
    radius = j + 1
    poly!(ax, Makie.GeometryBasics.Polygon(Circle(zero(Point2f), radius+0.5), [Circle(zero(Point2f), radius-0.5)]), color = row.placed_from_left ? :blue : :green, alpha = 0.25)
    color = (; color = row.cw ? :blue : :green)
    lines!(ax, radius*Point2f.(reverse.(sincos.(row.total_dance .+ π/2))); color...)
    α = row.total_dance[1]
    α += row.cw ? π : 0
    scatter!(ax, radius*Point2f(reverse(sincos(row.total_dance[1] + π/2))); color..., marker = :utriangle, markersize=10, rotation = α + π/2)
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
        r = fill(j, length(row.total_dance))
        lines!(ax, range(0, 2pi, 100), r, color = row.placed_from_left ? :blue : :green, linewidth = 5, alpha = 0.25)
        color = (; color = row.cw ? :blue : :green)
        lines!(ax, row.total_dance, r; color...)
        scatter!(ax, row.total_dance[1], j; color..., marker = '|', markersize=10, rotation=row.total_dance[1])
        # θ = min(row.total_dance[1] + 0.2, row.total_dance[end])
        ls = range(j*row.total_dance[1], j*row.total_dance[end], step = (pi)*sign(step(row.total_dance)))[2:end-1]
        if isempty(ls)
            ls = range(j*row.total_dance[50], j*row.total_dance[50])
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
lines!(ax, j*Point2f.(reverse.(sincos.(row.total_dance))); color...)
scatter!(ax, j*Point2f(reverse(sincos(row.total_dance[1]))); color..., marker = '|', markersize=10, rotation=row.total_dance[1] + pi/2)
ls = range(j*row.total_dance[1], j*row.total_dance[end], step = (pi)*sign(step(row.total_dance)))
if length(ls) < 2
    l = j*row.total_danctotal_danctotal_dance[50]
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

