using CSV, DataFrames
using GLMakie, AlgebraOfGraphics
using GLM

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
    gm = glm(fm, g, Binomial())
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
