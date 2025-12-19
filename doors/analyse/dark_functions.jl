# ============================================================================
# CIRCULAR STATISTICS AND BOOTSTRAP FUNCTIONS
# ============================================================================
# Functions for analyzing directional data using circular statistics.
# Primary metric: Mean Resultant Vector Length (r̄)
# Statistical inference: Parametric bootstrap with beta regression
# ============================================================================

# ============================================================================
# GEOMETRIC CALCULATIONS
# ============================================================================

"""
    path_length_at(xy, l1)

Calculate total path length traveled until reaching distance l1 from origin.
Accumulates Euclidean distances between consecutive positions, stopping when
the beetle reaches l1 cm from the start point.

Returns cumulative path length (cm).
"""
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

"""
    get_exit_angle(xyp, r)

Extract heading angle when trajectory first reaches radius r from origin.

This captures the direction in which the beetle exits a circular region
centered at the POI. Used to assess directional bias at multiple scales.

Parameters:
- xyp: Trajectory (post-POI, centered on POI)
- r: Radius threshold (cm)

Returns angle in radians, or missing if trajectory never reaches radius r.
"""
function get_exit_angle(xyp, r)
    i = findfirst(≥(r) ∘ norm, xyp)
    if isnothing(i)
        return missing
    end
    x, y = xyp[i]
    atan(y, x)  # Angle in radians
end

# ============================================================================
# CIRCULAR STATISTICS
# ============================================================================

"""
    mean_resultant_vector(θ)

Compute mean resultant vector length from a set of angles.

This is the primary measure of directional consistency:
- r̄ = 0: angles uniformly distributed (no preferred direction)
- r̄ = 1: all angles identical (perfect consensus)
- 0 < r̄ < 1: partial directional bias

Method: Represent each angle as unit vector, compute vector mean, take length.
"""
function mean_resultant_vector(θ)
    norm(mean(SV ∘ sincos, θ))
end

"""
    critical_r(n, p = 0.95)

Compute critical value for mean resultant vector length under null hypothesis
of uniform circular distribution. Used for significance testing.

Parameters:
- n: Sample size
- p: Confidence level (default 95%)

Returns critical r̄ value. If observed r̄ > critical_r̄, reject uniformity.
"""
critical_r(n, p = 0.95) = critical_r̄(n, p)

"""
    sample_r̄(n)

Simulate mean resultant vector length from n uniformly distributed angles.
This generates one random sample under the null hypothesis of circular uniformity.

Returns r̄ ∈ [0, 1] for this random sample.
"""
function sample_r̄(n)
    p = 0.0im
    for _ in 1:n
        α = 2π*rand()  # Random angle uniformly distributed
        p += exp(α*im)  # Complex representation: sum unit vectors
    end
    return norm(p)/n  # Mean resultant vector length
end

"""
    critical_r̄(n, p, nn = 10^7)

Monte Carlo estimation of critical r̄ value.
Generates nn random samples from uniform circular distribution and computes
the p-th quantile (typically 95th percentile).

This provides the threshold for statistical significance testing.
"""
function critical_r̄(n, p, nn = 10^7)
    r̄s = Vector{Float64}(undef, nn)
    Threads.@threads for i in 1:nn
        r̄s[i] = sample_r̄(n)
    end
    return quantile(r̄s, p)
end

# ============================================================================
# BOOTSTRAP INFERENCE
# ============================================================================

"""
    _bootstrap(df, fm, factors, newdf)

Single bootstrap resample iteration.

Procedure:
1. Resample rows with replacement from original data
2. Compute mean resultant vector for each (factor, radius) combination
3. Fit beta regression model to resampled data
4. Generate predictions for newdf grid

Returns (p-values, predictions) or missing if model fitting fails.
"""
function _bootstrap(df, fm, factors, newdf)

    # df = lightdf
    # newdf = newlight

    n = nrow(df)
    df2 = flatten(df[sample(1:n, n), :], [:θs, :r])  # Resample with replacement
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
        row = (; Pair.(Symbol.(tbl.rownms), tbl.cols[tbl.pvalcol])...)  # Extract p-values
        # row = tbl.cols[tbl.pvalcol]
        (row,  predict(m, newdf))
    catch ex
        # @show df3, fm, newdf, ex
        missing  # Return missing if model fails to converge
    end
end

"""
    __bootstrap(df, fm, factors, newdf)

Bootstrap wrapper with retry logic.
Attempts _bootstrap up to 100 times, retrying on failures.
This handles occasional convergence issues in beta regression.
"""
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

"""
    bootstrap(df, fm, newdf, factors)

Perform parametric bootstrap with n=10,000 resamples.

Procedure:
1. Repeatedly resample data and fit beta regression model
2. Collect p-values from all successful model fits
3. Collect predictions across bootstrap resamples
4. Return distributions for confidence interval construction

Parameters:
- df: Original data with angles at multiple radii
- fm: Formula for beta regression (typically: mean_resultant_vector ~ factor + r)
- newdf: Prediction grid (fine-resolution radius values)
- factors: Experimental factors (e.g., [:light], [:dance_by], [:at_run])

Returns:
- DataFrame of p-value distributions (one row per bootstrap sample)
- Matrix of predictions (nrow(newdf) × n_bootstrap)
"""
function bootstrap(df, fm, newdf, factors)
    n = 10_000  # Number of bootstrap resamples
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
            continue  # Skip failed model fits
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

# ============================================================================
# STATISTICAL UTILITIES
# ============================================================================

"""
    formatp(p)

Format p-value for display.
- p > 0.05: Display as ">0.05"
- 0.001 < p ≤ 0.05: Display rounded value
- p ≤ 0.001: Display as "<0.001"
"""
function formatp(p)
    if p > 0.05
        ">0.05"
    elseif p > 0.001
        string(round(p, sigdigits = 2))
    else
        "<0.001"
    end
end

"""
    stats(x)

Compute summary statistics from bootstrap p-value distribution.

Returns vector of formatted strings:
1. Proportion of resamples with p > 0.05 (non-significant)
2. 2.5th percentile (lower bound of 95% CI)
3. Median
4. 97.5th percentile (upper bound of 95% CI)
5. Mode
6. Mean

These quantify uncertainty in p-values across bootstrap resamples.
"""
function stats(x)
    p = count(>(0), x)/length(x)
    if p > 0.5
        1 - p
    else
        p
    end
    #
    # q1, q2 = quantile(x, [0.025, 0.975])
    # if q1 < 0 && q2 > 0
    #     ">0.05"
    # else
    #     "<0.05"
    # end
    #
    # p = something(findfirst(>(0.05), sort!(x)), length(x))/length(x)
    # q1, med, q2 = quantile(x, [0.025, 0.5, 0.975])
    # mod = mode(x)
    # μ = mean(x)
    # [string(round(Int, 100p), "%"), formatp(q1), formatp(med), formatp(q2), formatp(mod), formatp(μ)]
end

# ============================================================================
# HIGH-LEVEL ANALYSIS FUNCTIONS
# ============================================================================

"""
    analyze_factor_bootstrap(df, factor; n_radii=100, n_bootstrap=10_000)

Perform bootstrap analysis for a single experimental factor.

This function encapsulates the complete workflow for testing differences in
directional consistency (mean resultant vector length) across experimental
conditions using parametric bootstrap with beta regression.

Workflow:
1. Create high-resolution prediction grid across radii
2. Run parametric bootstrap (default: 10,000 resamples)
3. Fit beta regression: mean_resultant_vector ~ factor + r
4. Extract 95% confidence intervals from bootstrap distribution
5. Compute p-value statistics across bootstrap resamples

Parameters:
- df: DataFrame with columns: mean_resultant_vector, r (radii), factor, color
- factor: Symbol for experimental factor (:light, :dance_by, :at_run)
- n_radii: Number of radius points for smooth prediction curves (default: 100)
- n_bootstrap: Number of bootstrap resamples (default: 10,000)

Returns:
- newdata: DataFrame with interpolated predictions and confidence bands
  Columns: factor levels, r, color, mean_resultant_vector, lower, upper
- pvalues: DataFrame with bootstrap p-value distributions
  Columns: source (parameter names), proportion, Q2.5, median, Q97.5, mode, mean
"""
function analyze_factor_bootstrap(df, factor; n_radii=100, n_bootstrap=10_000)
    # Construct beta regression formula
    fm = FormulaTerm(Term(:mean_resultant_vector), (Term(factor), Term(:r)))

    # Create prediction grid with high resolution for smooth curves
    rl_range = extrema(df.r)
    rl2 = range(rl_range..., n_radii)
    newdata = combine(groupby(df, factor),
                     :r => first => :r,
                     :color => first => :color)
    newdata.r .= Ref(rl2)
    newdata = flatten(newdata, :r)
    newdata.mean_resultant_vector .= 0.0

    # Run bootstrap and extract confidence intervals
    pc, c = bootstrap(df, fm, newdata, [factor])
    select!(pc, Not(Symbol("(Precision)")))

    # Extract 2.5%, 50%, 97.5% quantiles for confidence bands
    y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
    newdata.lower .= getindex.(y, 1)
    newdata.mean_resultant_vector .= getindex.(y, 2)
    newdata.upper .= getindex.(y, 3)

    # Format p-value statistics from bootstrap distribution
    pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)

    return newdata, pvalues
end
