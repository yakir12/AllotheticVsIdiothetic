critical_r(n, p = 0.05) = sqrt(-4n * log(p) + 4n - log(p)^2 + 1)/(2n)

function mean_resultant_vector(θ)
    norm(mean(SV ∘ sincos, θ))
end

function _bootstrap(df)
    n = nrow(df)
    df2 = flatten(df[sample(1:n, n), :], [:θs, :r])
    df3 = combine(groupby(df2, [:light, :dance_by, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
    df3.light = categorical(df3.light)
    levels!(df3.light, ["remain", "dark"])
    df3.dance_by = categorical(df3.dance_by)
    levels!(df3.dance_by, ["no", "hold", "disrupt"])
    try
        m = BetaRegression.fit(BetaRegressionModel, fm, df3)
        tbl = coeftable(m)
        row = (; Pair.(Symbol.(tbl.rownms), tbl.cols[tbl.pvalcol])...)
        # row = tbl.cols[tbl.pvalcol]
        (row,  predict(m, newdf))
    catch ex
        missing
    end
end

# function _bootstrap(df)
#     n = nrow(df)
#     df2 = flatten(df[sample(1:n, n), :], [:θs, :r])
#     df3 = combine(groupby(df2, [:condition, :r]), :θs => mean_resultant_vector => :mean_resultant_vector)
#     df3.condition = categorical(df3.condition)
#     levels!(df3.condition, ["remain", "no", "hold", "disrupt"])
#     try
#         m = BetaRegression.fit(BetaRegressionModel, fm, df3)
#         tbl = coeftable(m)
#         row = (; Pair.(Symbol.(tbl.rownms), tbl.cols[tbl.pvalcol])...)
#         # row = tbl.cols[tbl.pvalcol]
#         (row,  predict(m, newdf))
#     catch ex
#         missing
#     end
# end
function bootstrap(df)
    n = 10_000
    row, y = _bootstrap(df)
    v1 = (; pairs(row)...)
    rows = Vector{typeof(v1)}(undef, n)
    # rows = DataFrame(Dict(pairs(row)))
    # empty!(rows)
    ys = Matrix{Float64}(undef, nrow(newdf), n)
    i = 0
    Threads.@threads for i in 1:n
    # while i < n
        rowy = _bootstrap(df)
        if ismissing(rowy)
            continue
        else
            # i += 1
            rows[i], ys[:, i] = rowy
            # row, ys[:, i] = rowy
            # push!(rows, row)
        end
    end
    keep = findall(i -> isassigned(rows, i), eachindex(rows))
    return DataFrame(rows[keep]), ys[:,keep]
    # return rows, stack(ys)
end

# function bootstrap(df)
#     n = 10_000
#     row, y = _bootstrap(df)
#     rows = DataFrame(Dict(pairs(row)))
#     empty!(rows)
#     ys = Matrix{Float64}(undef, nrow(newdf), n)
#     i = 0
#     while i < n
#         rowy = _bootstrap(df)
#         if ismissing(rowy)
#             continue
#         else
#             i += 1
#             row, ys[:, i] = rowy
#             push!(rows, row)
#         end
#     end
#     return rows, stack(ys)
# end



function formatp(p)
    if p > 0.05
        ">0.05"
    elseif p > 0.001
        string(round(p, sigdigits = 2))
    else
        "<0.001"
    end
end

function stats(x)
    p = something(findfirst(>(0.05), sort!(x)), length(x))/length(x)
    q1, med, q2 = quantile(x, [0.025, 0.5, 0.975])
    mod = mode(x)
    μ = mean(x)
    [string(round(Int, 100p), "%"), formatp(q1), formatp(med), formatp(q2), formatp(mod), formatp(μ)]
end


function get_exit_angle(xyp, r)
    i = findfirst(≥(r) ∘ norm, xyp)
    if isnothing(i)
        return missing
    end
    x, y = xyp[i]
    atan(y, x)
end

