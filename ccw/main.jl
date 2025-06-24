using Statistics
using CSV, DataFrames
using GLMakie, AlgebraOfGraphics
using GLM


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

norm2init(init, turn) = rem2pi(turn - init, RoundNearest)

function total_rotation(placed, direction, full, exited)
    Δ = direction == "ccw" ? placed - exited : exited - placed
    Δ += Δ < 0 ? 2pi : 0
    Δ += full*2pi
    Δ = direction == "cw" ? Δ : -Δ
    return Δ
end

df = CSV.read("data.csv", DataFrame, select = ["individual_number",
                                               "placed from",
                                               "placed from angle (degrees)",
                                               "rotation 1 direction",
                                               "rotation category measured",
                                               "exit angle (degrees)",
                                               "go down angle (degrees)",
                                               "full lap", 
                                               "total absolute degrees of rotation"])

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

transform!(groupby(df, :individual_number), "exit angle" => (θs -> angle(sum(exp.(1im*θs)))) => :μ)

transform!(df, ["μ", "placed from angle"] => ByRow(norm2init) => :placed)
transform!(df, ["μ", "go down angle"] => ByRow(norm2init) => :descended)
transform!(df, ["μ", "exit angle"] => ByRow(norm2init) => :exited)
select!(df, Not("μ", "placed from angle", "go down angle", "exit angle"))

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


sort!(df, ["individual_number", "placed from", "rotated2descent"])

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
