function path_length_at(xy, l1)
    p1 = xy[1]
    s = 0.0
    for p2 in xy[2:end]
        s += norm(p2 - p1)
        norm(p1) ≥ l1 && return s
        p1 = p2
    end
    return s
end

critical_r(n, p = 0.95) = critical_r̄(n, p)

function sample_r̄(n)
    p = 0.0im
    for _ in 1:n
        α = 2π*rand()
        p += exp(α*im)
    end
    return norm(p)/n
end

function critical_r̄(n, p, nn = 10^7)
    r̄s = Vector{Float64}(undef, nn)
    Threads.@threads for i in 1:nn
        r̄s[i] = sample_r̄(n)
    end
    return quantile(r̄s, p)
end

function mean_resultant_vector(θ)
    norm(mean(SV ∘ sincos, θ))
end

function _bootstrap(df, fm, factors, newdf)

    # df = lightdf
    # newdf = newlight

    n = nrow(df)
    df2 = flatten(df[sample(1:n, n), :], [:θs, :r])
    df3 = combine(groupby(df2, [factors..., :r]), :θs => mean_resultant_vector => :mean_resultant_vector)

    # df3[:, factor] = categorical(df3[:, factor])
    # levels!(df3[:, factor], levels)
    sort!(df3, factors)
    # df3.dance_by = categorical(df3.dance_by)
    # levels!(df3.dance_by, ["no", "hold", "disrupt"])
    try
        m = BetaRegression.fit(BetaRegressionModel, fm, df3)
        # @show df3
        tbl = coeftable(m)
        row = (; Pair.(Symbol.(tbl.rownms), tbl.cols[tbl.pvalcol])...)
        # row = tbl.cols[tbl.pvalcol]
        (row,  predict(m, newdf))
    catch ex
        # @show df3, fm, newdf, ex
        missing
    end
end

function __bootstrap(df, fm, factors, newdf)
    n = 100
    for i in 1:n
        rowy = _bootstrap(df, fm, factors, newdf)
        if !ismissing(rowy)
            return rowy
        end
    end
    error("tried $n times!")
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
function bootstrap(df, fm, newdf, factors)
    n = 10_000
    row, y = __bootstrap(df, fm, factors, newdf)
    v1 = (; pairs(row)...)
    rows = Vector{typeof(v1)}(undef, n)
    # rows = DataFrame(Dict(pairs(row)))
    # empty!(rows)
    ys = Matrix{Float64}(undef, nrow(newdf), n)
    i = 0
    # Threads.@threads 
    for i in 1:n
        # while i < n
        rowy = _bootstrap(df, fm, factors, newdf)
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

