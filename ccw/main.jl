using Statistics
using CSV, DataFrames
using GLMakie, AlgebraOfGraphics
# using MixedModels, GLM
using Random
using Optim

rm.(filter(==(".png") ∘ last ∘ splitext, readdir()))

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
                                               "full lap"])

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

fig1 = data(df) * (mapping(:mean_exit => rad2deg => "Mean exit", :mean_down => rad2deg) * visual(Scatter; label = "Mean down") + mapping(:mean_exit => rad2deg => "Mean exit", :optimal => rad2deg) * visual(Scatter; color = :red, label = "Mean optimal") + pregrouped([0], [1]) * visual(ABLines; color = :gray, label = "y = x")) |> draw(; axis = (; xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -180:180:180, ytickformat = "{:n}°", aspect = DataAspect(), width = 200, limits = ((-190, 190),(-190, 190))))
resize_to_layout!(fig1)

save("means.png", fig1)


function get_abs_residuals(overshoot, placed_from_left, start, total_dance)
    intercept = overshoot ? 0 : placed_from_left ? 2pi : -2pi
    ŷ = intercept - start
    Δ = ŷ - total_dance 
    placed_from_left ? Δ : -Δ
end

for μ in (:mean_exit, :mean_down, :optimal)
    df2 = transform(df, [:start, μ] => ByRow((start, μ) -> rem2pi(start - μ, RoundNearest)) => :start)
    transform!(df2, :start => ByRow(x -> sign(x) > 0) => :placed_from_left)
    df3 = subset(df2, :total_dance => ByRow(<(2pi) ∘ abs), :start => ByRow(x -> deg2rad(20) < abs(x) < deg2rad(160)))
    transform!(df3, [:placed_from_left, :cw] => ByRow(==) => :overshoot)
    transform!(df3, [:overshoot, :placed_from_left, :start, :total_dance] => ByRow(get_abs_residuals) => :residuals, :start => (s -> abs.(s)) => :x)
    fig = Figure()
    subgl_left = GridLayout()
    subgl_right = GridLayout()
    fig.layout[1, 1] = subgl_left
    fig.layout[1, 2] = subgl_right
    xy = data(df2) * mapping(:start => rad2deg => "Placement relative to intended direction", :total_dance => rad2deg => "Degree of rotation") * visual(Scatter; strokewidth = 1, strokecolor = :white, markersize = 6, alpha = 0.5, label = "data", legend = (; alpha = 1))
    abline = data(DataFrame(a = 360 .* (-2:2), b = -1, color = ["extra loop", "longer rotation direction", "shorter rotation direction", "longer rotation direction", "extra loop"])) * mapping(:a, :b, color = :color => sorter("shorter rotation direction", "longer rotation direction", "extra loop") => "") * visual(ABLines)
    vline = data(DataFrame(geometry = [Rect(-160, -360, 140, 720), Rect(20, -360, 140, 720)])) * mapping(:geometry) * visual(Poly; color = (:gray, 0.2), label = "included")
    toplot = vline + abline + xy
    g = draw!(subgl_left[1,1], toplot; axis = (; xticks = [-180, 0, 180], xtickformat = "{:n}°", yticks = -720:180:720, ytickformat = "{:n}°", aspect = DataAspect(), width = 200))
    legend!(subgl_right[1,1], g)
    draw!(subgl_right[2,1], data(df3) * mapping(:residuals => rad2deg) * visual(Hist); axis = (; ylabel = "#", xticks = [25, 180], xtickformat = "{:n}°", width = 200))
    resize_to_layout!(fig)
    save("relationship $μ.png", fig)
end




