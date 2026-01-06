using BetaRegression, DataFrames, AlgebraOfGraphics, GLMakie

d1 = Beta(10, 1.5)
lines(d1)
d2 = Beta(10, 2.5)
lines!(d2)

n = 200
df = DataFrame(grp = rand(string.('a':'c'), n))
transform!(df, :grp => ByRow(x -> x == "a" ? rand(d1) : rand(d2)) => :y) 
data(df) * mapping(:grp, :y) * visual(BoxPlot) |> draw()
# transform!(df, :grp => ByRow(x -> x == "b" ? "0" : x), renamecols = false) 
m = BetaRegression.fit(BetaRegressionModel, @formula(y ~ grp), df)


