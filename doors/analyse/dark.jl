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

runs = load_runs_and_calibs(results_dir; exclude_spontaneous=true)
@subset! runs :light .≠ "shift"  # Exclude shift experiments

# Identify spontaneous dances (vs induced dances) before processing
@transform! runs :dance_spontaneous = .!ismissing.(:spontaneous_end)

# filter out all the runs that are dance=no and have a spontaneous dance in them
@rsubset! runs :dance_by ≠ "no" || !:dance_spontaneous

# ============================================================================
# SECTION 2: TRAJECTORY PROCESSING PIPELINE
# ============================================================================
# Transform raw pixel coordinates into standardized trajectory representations.
# Key difference from shift.jl: for "remain" condition (lights stay on),
# POI time must be imputed since there's no actual intervention.
# ============================================================================

process_trajectories!(runs;
                      fix_intervention_jumps=true,
                      impute_poi_for_remain=true)

# Convert DimensionalData structures to plain arrays for easier manipulation
transform!(runs, [:pixels, :xy, :smooth, :centered2start, :cropped,
                  :rotated2poi, :centered2poi_and_cropped] .=> ByRow(parent),
           renamecols = false)

# ============================================================================
# SECTION 3: OVERVIEW FIGURES
# ============================================================================

@chain runs begin
    @rtransform! :path_length = sum(norm, diff(:centered2poi_and_cropped))
    @rtransform! :cord_length = norm(last(:centered2poi_and_cropped))
    @transform! :tortuosity = :path_length ./ :cord_length
    @transform! :straightness = 1 ./ :tortuosity
    @rtransform! :exit_angle = splat(atan)(reverse(last(:centered2poi_and_cropped)))
end

# angular_mean(θs) = atand(mean(sind, θs), mean(cosd, θs))
data(runs) * mapping(:straightness, :exit_angle, color = :condition) * visual(Scatter, strokewidth = 1, strokecolor = :white) 

function summerize_stats(df)
    n = nrow(df)
    df2 = df[sample(1:n, n), :]
    combine(groupby(df2, :condition), :exit_angle => mean_resultant_vector => :mean_resultant_vector, :straightness => mean => :straightness)
end

df = select(runs, :condition, :exit_angle, :straightness)
n_boot = 100000
d = Dict(c => Vector{Point2f}(undef, n_boot) for c in unique(runs.condition))
for i in 1:n_boot
    _df = summerize_stats(df)
    for row in eachrow(_df)
        d[row.condition][i] = Point2f(row.straightness, row.mean_resultant_vector)
    end
end
fig, ax, pl = datashader(d; async = false)

_df = select(runs, :condition, :exit_angle, :straightness)
n_boot = 100000
df = vcat((summerize_stats(_df) for _ in 1:n_boot)...)

contour_layer = AlgebraOfGraphics.density() * visual(Contour; levels = 4, linewidth = 2)
fig = data(df) * mapping(:straightness, :mean_resultant_vector, color = :condition) * contour_layer |> draw()
save_figure(fig, output, "straightness with mean resultant vector")

fig = data(df) * mapping(:straightness, :mean_resultant_vector, color = :condition) * visual(Scatter; markersize = 1, alpha = 0.1) |> draw()

# data(combine(groupby(runs, :condition), :straightness => mean => :straightness, :exit_angle => angular_mean => :exit_angle)) * mapping(:straightness, :exit_angle, color = :condition) * visual(Scatter; markersize = 25, strokewidth = 2, strokecolor = :gray) |> draw()

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
for _ax in fig.figure.content
    if _ax isa Axis
        for r  in (30, MAX_TRAJECTORY_LENGTH)
            lines!(_ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
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

p = bootstrap_dance_method(df)
CSV.write(joinpath(output, "dance_by_method.csv"), p)
p = bootstrap_dance(df)
CSV.write(joinpath(output, "dance_by.csv"), p)
p = bootstrap_light(df)
CSV.write(joinpath(output, "light.csv"), p)
p = bootstrap_at_run(df)
CSV.write(joinpath(output, "at_run.csv"), p)


# Assign colors to experimental groups for visualization
@transform! df :color = colorant"gray"  # Default grey for reference
light = @subset df :dance_by .== "no" :at_run .== 1  # Light comparison subset
dance_by = @subset df :light .== "dark" :at_run .== 1 :dance_by .≠ "hold" # Dance induction comparison
@transform! dance_by :danced = ifelse.(:dance_by .≠ "no", "yes", "no")
at_run = @subset df :dance_by .== "no" :light .== "dark"  # Familiarity comparison

# Color coding for each comparison:
@transform!(subset(light, :light => ByRow(==("remain")), view = true), :color = colors[1])  # Lights remain
@transform!(subset(dance_by, :danced => ByRow(==("yes")), view = true), :color = colors[2])  # Lights remain
@transform!(subset(at_run, :at_run => ByRow(==(10)), view = true), :color = colors[4])  # 10th run

# Compute mean resultant vector for each (factor, radius) combination
light_summary = flatten(light, [:θs, :r])
dance_by_summary = flatten(dance_by, [:θs, :r])
at_run_summary = flatten(at_run, [:θs, :r])

light_summary1 = combine(groupby(light_summary, [:light, :r]),
                       :θs => mean_resultant_vector => :mean_resultant_vector,
                       :color => first => :color)
dance_by_summary1 = combine(groupby(dance_by_summary, [:danced, :r]),
                          :θs => mean_resultant_vector => :mean_resultant_vector,
                          :color => first => :color)
at_run_summary1 = combine(groupby(at_run_summary, [:at_run, :r]),
                        :θs => mean_resultant_vector => :mean_resultant_vector,
                        :color => first => :color)

# Run bootstrap analysis for each of the three comparisons
newlight, pvalueslight = analyze_factor_bootstrap(light_summary, :light)
newdance, pvaluesdance_by = analyze_factor_bootstrap(dance_by_summary, :danced)
newat_run, pvaluesat_run = analyze_factor_bootstrap(at_run_summary, :at_run)

# Save p-value tables to CSV
# CSV.write(joinpath(output, "light.csv"), pvalueslight)
# CSV.write(joinpath(output, "dance_by.csv"), pvaluesdance_by)
# CSV.write(joinpath(output, "at_run.csv"), pvaluesat_run)

# Save fitted data (for replotting without rerunning bootstrap)
# CSV.write(joinpath(output, "light_data.csv"), newlight)
# CSV.write(joinpath(output, "dance_by_data.csv"), newdance)
# CSV.write(joinpath(output, "at_run_data.csv"), newat_run)

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
g3 = @subset dance_by :danced .≠ "no"
g4 = @subset at_run :at_run .== 10

width = 250
height = 250

fig = Figure();
# Top panels: Labels for experimental conditions
for (i, g) in enumerate((g1, g2, g4, g3))
    for (j, col) in enumerate([:light_styled, :dance_by_styled, :at_run_styled])
        Label(fig[j, i], g[1, col])
    end
    # Polar plot: trajectories emanating from center (POI)
    _ax = PolarAxis(fig[4, i], rlimits = (0, 44); width, height, thetaticklabelsvisible = false)
    for row in eachrow(g)
        θ = splat(atan).(row.centered2poi_and_cropped) .+ π/2  # Convert to polar angle
        r = norm.(row.centered2poi_and_cropped)  # Radial distance
        lines!(_ax, θ, r, color = row.color)
    end
end
# Bottom panels: Mean resultant vector plots with confidence bands
leg = []
for (i, g) in enumerate((newlight, newat_run, newdance))
    _ax = Axis(fig[5, i + 1]; limits = ((5, 29),(0,1)), ylabel = "Resultant mean vector length",  height)
    for (k, g) in pairs(groupby(g, Not(:color, :lower, :mean_resultant_vector, :r, :upper)))
        b = band!(_ax, g.r, g.lower, g.upper, color = RGBA(g.color[1], 0.25))  # Confidence band
        ll = lines!(_ax, g.r, g.mean_resultant_vector, color = g.color[1])  # Median line
        push!(leg, k => [b, ll])
    end
    if i > 1
        hideydecorations!(_ax, grid = false, minorgrid = false)
    end
end
Label(fig[6,2:end], "Radial distance from POI (cm)")
# Legend(fig[5,1], last.(leg[[2,1,4,6]]), ["Lights turn off", "Dance induced", "Lights remain on", rich("at the 10", superscript("th"), " run")])
# Legend(fig[5,1], last.(leg[[2,4,6,1]]), ["Lights turn off", "Lights remain on", rich("at the 10", superscript("th"), " run"), "Dance induced"])
Legend(fig[5,1], last.(leg[[1,2,4,5]]), ["Lights turn off", "Lights remain on", rich("at the 10", superscript("th"), " run"), "Dance induced"])
for i in 1:3
    rowgap!(fig.layout, i, 0)
end
resize_to_layout!(fig)
save_figure(fig, output, "darkness 3 comparisons");

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
