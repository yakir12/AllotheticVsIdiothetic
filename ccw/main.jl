using Statistics
using CSV, DataFrames
using AlgebraOfGraphics
# using GLMakie
using CairoMakie
# using MixedModels, GLM
using Random
using Optim

pt = 4/3
inch = 96
cm = inch / 2.54

# rm.(filter(==(".png") ∘ last ∘ splitext, readdir()))

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

function fix_stop(start, stop, lastcw)
    start ≠ stop && return stop
    Δ = lastcw ? -0.01 : 0.01
    return stop + Δ
end

df = CSV.read("data.csv", DataFrame, select = ["individual_number",
                                               "placed from angle (degrees)",
                                               "rotation 1 direction",
                                               "rotation category measured",
                                               "exit angle (degrees)",
                                               "go down angle (degrees)",
                                               "full lap",
                                              "total absolute degrees of rotation"])

transform!(groupby(df, :individual_number), :individual_number => (g -> 1:length(g)) => :n)

dforg = copy(df)

# transform!(groupby(df, :individual_number), :individual_number => (x -> 1:length(x)) => :n)

subset!(df, "rotation category measured" => ByRow(∈(("cw", "ccw"))))
@assert all(df[!, "rotation 1 direction"] .== df[!, "rotation category measured"])
select!(df, Not("rotation category measured"))

transform!(df, ["placed from angle (degrees)",
                "exit angle (degrees)",
                "go down angle (degrees)"] .=> ByRow(deg2rad) .=> 
           ["placed from angle",
            "exit angle",
            "go down angle"])
select!(df, Not("placed from angle (degrees)",
                "exit angle (degrees)",
                "go down angle (degrees)"))

transform!(df, "rotation 1 direction" => ByRow(==("cw")) => :cw)
select!(df, Not("rotation 1 direction"))

transform!(df, ["placed from angle", "go down angle", "cw", "full lap"] => ByRow(get_total_rotation) => :total_dance)
rename!(df, "placed from angle" => :start)

transform!(df, ["go down angle", "exit angle", "cw"] => ByRow(fix_stop) => "exit angle")
transform!(df, ["go down angle", "exit angle", "cw"] => ByRow(get_total_rotation) => :down2exit)
# rename!(df, "go down angle" => :down)
rename!(df, "exit angle" => :exit)

select!(df, Not("go down angle", "full lap"))

function to_minimize(start, total_dance, μ)
    s = 0.0
    for (_x, y) in zip(start, total_dance)
        x = rem2pi(_x .- μ, RoundNearest)
        r, _ = findmin(i -> abs(y - (i*2π - x)), -2:2)
        s += r
    end
    return s
end

angular_mean(θs) = angle(sum(exp, 1im*θs))

transform!(groupby(df, :individual_number), :exit => angular_mean => :mean_exit, [:start, :total_dance] => ((s, p) -> angular_mean(s .+ p)) => :mean_down, [:start, :total_dance] => ((start, total_dance) -> optimize(μ -> to_minimize(start, total_dance, μ), -pi, pi).minimizer) => :optimal)

fig1 = data(df) * (mapping(:mean_exit => rad2deg => "Mean exit", :mean_down => rad2deg) * visual(Scatter; label = "Mean down") + mapping(:mean_exit => rad2deg => "Mean exit", :optimal => rad2deg) * visual(Scatter; color = :red, label = "Mean optimal") + pregrouped([0], [1]) * visual(ABLines; color = :gray, label = "y = x")) |> draw(; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -180:180:180, ytickformat = "{:n}°", aspect = DataAspect(), width = 200, limits = ((-190, 190),(-190, 190))))
resize_to_layout!(fig1)

save("means.pdf", fig1)


function get_abs_residuals(overshoot, placed_from_left, start, total_dance)
    n = abs(total_dance) < 2π ? 1 : 2
    intercept = overshoot ? 0 : placed_from_left ? n*2π : -n*2π
    ŷ = intercept - start
    Δ = ŷ - total_dance 
    placed_from_left ? Δ : -Δ
end

example_individual = 24

for μ in (:mean_exit, :mean_down, :optimal)
    df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)
    transform!(df2, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
    df3 = subset(df2, :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
    # df3 = subset(df2, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
    transform!(df3, [:placed_from_left, :cw] => ByRow(==) => :overshoot)
    transform!(df3, [:overshoot, :placed_from_left, :start, :total_dance] => ByRow(get_abs_residuals) => :residuals, :start => (s -> abs.(s)) => :x)

    data(df3) * mapping(:start => abs ∘ rad2deg, :residuals => rad2deg, row = :overshoot => renamer(true => "overshut", false => "undershut"), col = :placed_from_left => renamer(true => "placed from left", false => "placed from right")) * (linear() + visual(Scatter)) |> draw() |> save("scatter $μ.png")

    fig = Figure(size = (12cm, 10cm), fontsize = 8pt, fonts = (; regular = "Helvetica"))
    subgl_left = GridLayout()
    subgl_right = GridLayout()
    fig.layout[1, 1] = subgl_left
    fig.layout[1, 2] = subgl_right
    xy = data(subset(df2, :individual_number => ByRow(≠(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation") * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5)
    xy42 = data(subset(df2, :individual_number => ByRow(==(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation") * visual(Scatter; color = :red, strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5)
    abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
    vline = data(DataFrame(m = [-160, 20], M = [-20, 160])) * mapping(:m, :M) * visual(VSpan; color = (:gray, 0.2))
    # vline = data(DataFrame(geometry = [Rect(-160, -360, 140, 720), Rect(20, -360, 140, 720)])) * mapping(:geometry) * visual(Poly; color = (:gray, 0.2), label = "included")
    toplot = vline + abline + xy + xy42
    g = draw!(subgl_left[1,1], toplot; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xlabel = "Initial body orientation", xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect()))

    legend!(subgl_right[1,1], g; framevisible = false)
    draw!(subgl_right[2,1], data(df3) * mapping(:residuals => rad2deg => "Residuals") * visual(Hist; bins = 25); axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, ylabel = "#", xticks = 0:180:1000, xtickformat = "{:n}°"))
    dance_number = data(dforg) * (mapping(:n => "Dance number", "total absolute degrees of rotation" => "Total absolute rotation") * visual(BoxPlot) + mapping(:n => "Dance number", "total absolute degrees of rotation" => "Total absolute rotation", group = :individual_number => nonnumeric) * visual(Lines; alpha = 0.1, linewidth = 1))
    draw!(subgl_right[3,1], dance_number; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, yticks = 0:360:1000, ytickformat = "{:n}°"))

    Label(subgl_left[1, 1,  TopLeft()], "A", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
    Label(subgl_right[2, 1, TopLeft()], "B", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)
    Label(subgl_right[3, 1, TopLeft()], "C", fontsize = 12pt, padding = (0, 5, 5, 0), halign = :right)

    resize_to_layout!(fig)
    save("relationship $μ.pdf", fig)
    save("relationship $μ.png", fig)
    resid = rad2deg.(df3.residuals)
    @show μ, round(Int, mean(resid)), round(Int, std(resid))
end




μ = :mean_down
df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)
transform!(df2, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
df3 = subset(df2, :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
# df3 = subset(df2, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
transform!(df3, [:placed_from_left, :cw] => ByRow(==) => :overshoot)
transform!(df3, [:overshoot, :placed_from_left, :start, :total_dance] => ByRow(get_abs_residuals) => :residuals, :start => (s -> abs.(s)) => :x)

fig = Figure(size = (8cm, 18cm), fontsize = 12pt, fonts = (; regular = "Helvetica"))
xy = data(subset(df2, :individual_number => ByRow(≠(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation", layout = :individual_number => nonnumeric) * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
xy42 = data(subset(df2, :individual_number => ByRow(==(example_individual)))) * mapping(:start => rad2deg => "Initial body orientation", :total_dance => rad2deg => "Total rotation", layout = :individual_number => nonnumeric) * visual(Scatter; color = :red, strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["additional lap", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "additional lap"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "additional lap") => "") * visual(ABLines)
toplot = abline + xy + xy42
g = draw!(fig[1,1], toplot; axis = (; xticksize = 3, yticksize = 3, titlesize = 8pt, titlefont = :regular, xlabelsize = 8pt, ylabelsize = 8pt, xticklabelsize = 6pt, yticklabelsize = 6pt, xlabel = "Initial body orientation", xticks = [-180, 180], xtickformat = "{:n}°", yticks = -720:360:720, ytickformat = "{:n}°", aspect = DataAspect()))
resize_to_layout!(fig)
save("inidividual relationship $μ.pdf", fig)
save("inidividual relationship $μ.png", fig)

# mean_exit, 25, 54
# mean_down, 26, 53
# optimal, 27, 57

shuffle!(df2)
sort!(df2, [:placed_from_left, :cw, :start])

fig = Figure()
for (i, (k, g)) in enumerate(pairs(groupby(df2, :individual_number)))
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
save("circles.pdf", fig)
