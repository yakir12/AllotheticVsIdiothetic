# ============================================================================
# DARKNESS EXPERIMENT ANALYSIS
# ============================================================================
# This script analyzes beetle movement trajectories in darkness experiments
# using circular statistics to assess directional orientation.
#
# Scientific Goal:
#   Determine whether beetles maintain directional consistency when:
#   - Lights turn off vs remain on
#   - Dance is induced (disrupt/hold) vs spontaneous behavior
#   - At first run (novel) vs 10th run (familiar)
#
# Key Metric: Mean Resultant Vector Length (r̄)
#   - Circular statistics measure of directional consistency
#   - r̄ = 0: uniform/random directions
#   - r̄ = 1: all trajectories point same direction
#   - Computed at multiple radii from Point Of Interest (POI)
#
# Analysis Pipeline:
#   1. Load and filter non-shift experiments (exclude spontaneous dances)
#   2. Process trajectories: rectify, smooth, center, rotate, crop
#   3. Extract exit angles at multiple radii from POI
#   4. Compute mean resultant vector with bootstrap confidence intervals
#   5. Generate polar plots and statistical comparisons
# ============================================================================

include("common_imports.jl")
using Interpolations
using OhMyThreads
using LsqFit
using Optim
using CategoricalArrays
using IntervalSets
using QuadGK
using BetaRegression

include("general_functions.jl")
include("dark_functions.jl")

# Prepare output directory for figures and results
output = setup_output_dir("dark")

const results_dir = "../track_calibrate/tracks and calibrations"

# ============================================================================
# SECTION 1: DATA LOADING AND FILTERING
# ============================================================================
# Load tracking data, excluding:
#   - "shift" experiments (analyzed separately in shift.jl)
#   - Spontaneous dances (focus on controlled interventions)
# ============================================================================

runs = @chain joinpath(results_dir, "runs.csv") begin
    CSV.read(DataFrame)
    @select Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size)
    @subset :light .≠ "shift"  # Exclude shift experiments
    @subset! ismissing.(:spontaneous_end)  # Only induced dances, no spontaneous
    @transform :tij_file = joinpath.(results_dir, :tij_file)
    @transform :location = "Lund"
end
calibs = @chain joinpath(results_dir, "calibs.csv") begin
    CSV.read(DataFrame)
    @transform :calibration_file = joinpath.(results_dir, :calibration_id)
    @transform :rectify = get_calibration.(:calibration_file)  # Create rectification function
    @select Cols(:calibration_id, :rectify)
end
leftjoin!(runs, calibs, on = :calibration_id)

# ============================================================================
# SECTION 2: TRAJECTORY PROCESSING PIPELINE
# ============================================================================
# Transform raw pixel coordinates into standardized trajectory representations.
# Key difference from shift.jl: for "remain" condition (lights stay on),
# POI time must be imputed since there's no actual intervention.
# ============================================================================

@chain runs begin
    @select! Not(:calibration_id)
    @transform! :condition = string.(:light, " ", :dance_by, " ", :at_run)
    @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
    @rename! :intervention = :poi

    # Load tracking data: time (t), pixel row (i), pixel column (j)
    @transform! :pixels = get_tij.(:tij_file)
    # Remove consecutive duplicate positions (stops)
    @transform! :pixels = remove_stops.(:pixels)

    # Rectify: convert pixel coordinates to real-world coordinates (cm)
    @rtransform! :xy = :rectify.(:pixels)

    # Fix discontinuities for induced dances (skip "no" dance condition)
    @aside @chain _ begin
        @subset(:dance_by .≠ "no"; view = true)
        @transform! :jump = glue_intervention!.(:xy, :intervention)
    end

    # For "remain" condition: no actual intervention, so impute POI time
    # (estimated as when beetle has moved 10 cm from start)
    @aside @chain _ begin
        @subset(:light .== "remain"; view = true)
        @transform! :intervention = impute_poi_time.(:xy)
    end

    @transform! :poi = :intervention
    disallowmissing!(:poi)
    @select! Not(:intervention)

    # Clean and standardize trajectories
    @transform! :xy = remove_loops.(:xy)  # Remove self-intersecting loops
    @transform! :xy = sparseify.(:xy)  # Resample at 0.5s intervals
    @transform! :smooth = smooth.(:xy)  # Apply smoothing spline

    # Standardize trajectory orientation and position:
    @transform! :centered2start = center2start.(:smooth)  # Origin at start
    @transform! :cropped = cropto.(:centered2start, MAX_TRAJECTORY_LENGTH)  # Crop to 50 cm radius
    @transform! :rotated2poi = rotate2poi.(:cropped, :poi)  # POI at 90°
    @transform! :centered2poi_and_cropped = center2poi_and_crop.(:rotated2poi, :poi)  # Origin at POI

    # Convert DimensionalData structures to plain arrays
    transform!([:pixels, :xy, :smooth, :centered2start, :cropped, :rotated2poi, :centered2poi_and_cropped] .=> ByRow(parent), renamecols = false)
end

# ============================================================================
# SECTION 3: OVERVIEW FIGURES
# ============================================================================

# Overview: Plot all tracks to visually check data quality
fig = (pregrouped(runs.smooth => first => "X (cm)", runs.smooth => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines; color = :red) + pregrouped(runs.xy => first => "X (cm)", runs.xy => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines)) |> draw(; axis = (; width = 400, height = 400, limits = ((-MAX_TRAJECTORY_LENGTH, MAX_TRAJECTORY_LENGTH), (-MAX_TRAJECTORY_LENGTH, MAX_TRAJECTORY_LENGTH))));
save_figure(fig, output, "overview_dark")

# Figure 1: Sample trajectories (10 per condition)
n = 10
df = @chain runs begin
    @subset norm.(last.(:cropped)) .≈ MAX_TRAJECTORY_LENGTH
    @groupby :dance_by
    @transform :n = 1:length(:dance_by)
    @subset :n .≤ n
end
@assert all(==(n), combine(groupby(df, :dance_by), nrow).nrow)

fig = pregrouped(df.cropped => first => "X (cm)", df.cropped => last => "Y (cm)", col = df.dance_by, color = df.y2025) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content
    if ax isa Axis
        for r  in (30, MAX_TRAJECTORY_LENGTH)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
save_figure(fig, output, "figure1")

# ============================================================================
# SECTION 4: MAIN ANALYSIS - THREE PAIRWISE COMPARISONS
# ============================================================================
# Test three hypotheses about directional consistency:
#   1. Light effect: Do beetles orient differently when lights turn off?
#   2. Dance induction effect: Does induced dance affect orientation?
#   3. Familiarity effect: Does repeated exposure (run 10) change orientation?
#
# Statistical Approach:
#   - Extract exit angles at multiple radii (5-29 cm from POI)
#   - Compute mean resultant vector length (circular statistics)
#   - Bootstrap 10,000 resamples for confidence intervals
#   - Beta regression to model r̄ as function of radius + experimental factors
# ============================================================================

using Colors

# Define distinguishable colors for each experimental group
colors = distinguishable_colors(4, [colorant"white", colorant"gray",  colorant"black"], dropseed = true, transform = deuteranopic)

df = copy(runs)

# Create human-readable labels for plotting
@chain df begin
    @rtransform! :light_styled = :light == "dark" ? "lights turn off" : "lights remain on"
    @rtransform! :dance_by_styled = :dance_by == "no" ? "dance not induced" : "dance induced"
    @rtransform! :at_run_styled = :at_run == 1 ? rich("at the 1", superscript("st"), " run") : rich("at the 10", superscript("th"), " run")
end

# Extract exit angles at multiple radii from POI
nr = 3  # Number of radii to sample
l1 = floor(Int, minimum(norm ∘ last, df.centered2poi_and_cropped))
rl = range(5, l1, nr)  # Radii from 5 cm to maximum common distance
transform!(df, :centered2poi_and_cropped => ByRow(xy -> get_exit_angle.(Ref(xy), rl)) => :θs)
df.r .= Ref(rl)

# Assign colors to experimental groups for visualization
@transform! df :color = colorant"gray"  # Default grey for reference
light = @subset df :dance_by .== "no" :at_run .== 1  # Light comparison subset
dance_by = @subset df :light .== "dark" :at_run .== 1  # Dance induction comparison
at_run = @subset df :dance_by .== "no" :light .== "dark"  # Familiarity comparison

# Color coding for each comparison:
@transform!(subset(light, :light => ByRow(==("remain")), view = true), :color = colors[1])  # Lights remain
@rtransform!(groupby(subset(dance_by, :dance_by => ByRow(≠("no")), view = true), :dance_by), :color = :dance_by == "disrupt" ? colors[2] : colors[3])  # Induced dances
@transform!(subset(at_run, :at_run => ByRow(==(10)), view = true), :color = colors[4])  # 10th run

"""
    fun(lightdf, factor)

Perform bootstrap analysis for a single experimental factor.
Fits beta regression model with formula: mean_resultant_vector ~ factor + r
Returns interpolated predictions with confidence bands and p-value statistics.
"""
function fun(lightdf, factor)
    fm = FormulaTerm(Term(:mean_resultant_vector), (Term(factor), Term(:r)))
    newlight = combine(groupby(lightdf, factor), :r => first ∘ first => :r, :color => first => :color)
    nr2 = 100
    rl2 = range(extrema(rl)..., nr2)  # High-resolution grid for smooth curves
    newlight.r .= Ref(rl2)
    newlight = flatten(newlight, :r)
    newlight.mean_resultant_vector .= 0.0
    pc, c = bootstrap(lightdf, fm, newlight, [factor])  # 10,000 bootstrap resamples
    select!(pc, Not(Symbol("(Precision)")))
    # Extract 95% confidence intervals
    y = quantile.(skipmissing.(eachrow(c)), Ref([0.025, 0.5, 0.975]))
    newlight.lower .= getindex.(y, 1)
    newlight.mean_resultant_vector .= getindex.(y, 2)
    newlight.upper .= getindex.(y, 3)
    pvalues = combine(pc, DataFrames.All() .=> stats, renamecols = false)
    pvalues.what = ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"]
    df2 = stack(pvalues, Not(:what), variable_name = :source)
    pvalues = combine(groupby(df2, :source), [:what, :value] => ((what, value) -> (; Pair.(Symbol.(what), value)...)) => ["proportion", "Q2.5", "median", "Q97.5", "mode", "mean"])
    return newlight, pvalues
end

# Run bootstrap analysis for each of the three comparisons
newlight, pvalueslight = fun(light, :light)
newdance, pvaluesdance_by = fun(dance_by, :dance_by)
newat_run, pvaluesat_run = fun(at_run, :at_run)

# Save p-value tables to CSV
CSV.write(joinpath(output, "light.csv"), pvalueslight)
CSV.write(joinpath(output, "dance_by.csv"), pvaluesdance_by)
CSV.write(joinpath(output, "at_run.csv"), pvaluesat_run)

# Save fitted data (for replotting without rerunning bootstrap)
CSV.write(joinpath(output, "light_data.csv"), newlight)
CSV.write(joinpath(output, "dance_by_data.csv"), newdance)
CSV.write(joinpath(output, "at_run_data.csv"), newat_run)

# ============================================================================
# SECTION 5: FIGURE GENERATION - THREE COMPARISONS
# ============================================================================
# Multi-panel figure showing:
#   - Top 3 rows: Experimental condition labels for each comparison
#   - Row 4: Polar plots of trajectories (black lines on polar axes)
#   - Row 5: Mean resultant vector vs radius with confidence bands
# ============================================================================

g1 = @subset light :light .== "dark"
g2 = @subset light :light .== "remain"
g3 = @subset dance_by :dance_by .≠ "no"
g4 = @subset at_run :at_run .== 10

width = 250
height = 250
fig = Figure()

# Top panels: Labels for experimental conditions
for (i, g) in enumerate((g1, g2, g3, g4))
    for (j, col) in enumerate([:light_styled, :dance_by_styled, :at_run_styled])
        Label(fig[j, i], g[1, col])
    end
    # Polar plot: trajectories emanating from center (POI)
    ax = PolarAxis(fig[4, i], rlimits = (0, 44); width, height, thetaticklabelsvisible = false)
    for row in eachrow(g)
        θ = splat(atan).(row.centered2poi_and_cropped) .+ π/2  # Convert to polar angle
        r = norm.(row.centered2poi_and_cropped)  # Radial distance
        lines!(ax, θ, r, color = row.color)
    end
end

# Bottom panels: Mean resultant vector plots with confidence bands
leg = []
for (i, g) in enumerate((newlight, newdance, newat_run))
    ax = Axis(fig[5, i + 1]; limits = ((5, 29),(0,1)), ylabel = "Resultant mean vector length",  height)
    for (k, g) in pairs(groupby(g, Not(:color, :lower, :mean_resultant_vector, :r, :upper)))
        b = band!(ax, g.r, g.lower, g.upper, color = RGBA(g.color[1], 0.25))  # Confidence band
        ll = lines!(ax, g.r, g.mean_resultant_vector, color = g.color[1])  # Median line
        push!(leg, k => [b, ll])
    end
    if i > 1
        hideydecorations!(ax, grid = false, minorgrid = false)
    end
end
Label(fig[6,2:end], "Radial distance from POI (cm)")
Legend(fig[5,1], last.(leg[[1,2,3,5,7]]), ["Lights turn off", "Lights remain on", "Dance induced by disruption", "Dance induced by holding", rich("at the 10", superscript("th"), " run")])
for i in 1:3
    rowgap!(fig.layout, i, 0)
end
resize_to_layout!(fig)

save_figure(fig, output, "darkness 3 comparisons")

# ============================================================================
# SECTION 6: STRAIGHTNESS METRIC
# ============================================================================
# Compute trajectory straightness as ratio: straight-line distance / path length
# Lower values = more tortuous path
# Higher values = straighter path
# ============================================================================

l1 = floor(Int, minimum(norm ∘ last, df.centered2poi_and_cropped))
@transform! df :straightness = l1 ./ path_length_at.(:centered2poi_and_cropped, l1)
CSV.write(joinpath(output, "straightness.csv"), select(df, [:run_id, :straightness]))
