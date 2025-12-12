# ============================================================================
# SHIFT EXPERIMENT ANALYSIS
# ============================================================================
# This script analyzes beetle movement trajectories in light shift experiments.
#
# Scientific Goal:
#   Compare induced vs spontaneous "dance" behavior in response to light shifts.
#   Test whether turning patterns differ between:
#   - Dance induction method (induced by stimulus vs spontaneous)
#   - Light source type (LED in Lund vs Sun in Bela-Bela)
#
# Analysis Pipeline:
#   1. Load and merge tracking data from indoor and outdoor50 experiments
#   2. Process trajectories: rectify, smooth, center, rotate, crop
#   3. Extract turning angle profiles and fit logistic curves
#   4. Statistical modeling: GLM and beta regression on turn parameters
#   5. Generate publication-quality figures
# ============================================================================

include("common_imports.jl")
using LsqFit
using BetaRegression

include("general_functions.jl")
include("shift_functions.jl")

# Prepare output directory for figures
output = setup_output_dir("shift")

# ============================================================================
# SECTION 1: DATA LOADING AND MERGING
# ============================================================================
# Load tracking data from both indoor (Lund, LED stimulus) and outdoor50
# (Bela-Bela, Sun stimulus) experiments. Merge with camera calibration data
# to enable pixel-to-world coordinate transformation.
# ============================================================================

# Load indoor experiment data (Lund location, LED stimulus)
runs = load_runs_and_calibs("../track_calibrate/tracks and calibrations"; filter_light="shift")

# ============================================================================
# SECTION 2: TRAJECTORY PROCESSING PIPELINE
# ============================================================================
# Transform raw pixel coordinates into standardized trajectory representations
# suitable for statistical analysis. Key steps:
#   1. Rectification: pixel → world coordinates (cm)
#   2. Cleaning: remove stops, loops, discontinuities
#   3. Standardization: smooth, center, crop, rotate trajectories
#   4. Feature extraction: fit logistic curves to turning behavior
# ============================================================================

# Identify spontaneous dances (vs induced dances) before processing
@transform! runs :dance_spontaneous = .!ismissing.(:spontaneous_end)

# add another category of dance_by that codes the spontaneous dances
@rtransform! runs :dance_by = :dance_spontaneous ? "spontaneous" : :dance_by

# Apply standard trajectory processing pipeline with spontaneous dance handling
process_trajectories!(runs;
                      fix_intervention_jumps=true,
                      handle_spontaneous=true)

# Extract turning angle profile around POI (±5 to +15 cm path length)
@transform! runs :θ = get_turn_profile.(:smooth, :poi)

# Fit logistic curve to turning behavior: θ(l) = L/(1 + exp(-k(l - x₀))) - y₀
# Parameters: L = total turn, k = steepness, x₀ = inflection point, y₀ = offset
@transform! runs $AsTable = fit_logistic.(:θ)
transform!(runs, :ks => ByRow(identity) => [:L, :k, :x₀, :y₀])
@rtransform! runs :θs = logistic.(:θ.l, :L, :k, :x₀, :y₀)

# Convert to plain arrays
transform!(runs, [:pixels, :xy, :smooth, :centered2start, :cropped, :rotated2poi, :centered2poi_and_cropped] .=> ByRow(parent), renamecols = false)
# ============================================================================
# SECTION 3: PRELIMINARY TESTS
# ============================================================================
# Quick exploratory tests on exit angle and turn magnitude
# ============================================================================

#### Additional tests

# @chain runs begin
#     @rsubset :dance_spontaneous
#     @select :comment
# end
# @chain runs begin
#     @rsubset !ismissing(:comment) && contains(:comment, "spontaneous")
#     @select :dance_spontaneous
# end

@show combine(groupby(runs, [:location, :light, :dance_by, :at_run]), nrow => :n)

# @show combine(groupby(runs, [:location, :light, :dance_by2, :at_run]), nrow => :n)

# ============================================================================
# SECTION 4: FIGURE GENERATION - TRAJECTORY VISUALIZATIONS
# ============================================================================

############ Overview: Plot all tracks to visually check data quality

fig = (pregrouped(runs.smooth => first => "X (cm)", runs.smooth => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines; color = :red) + pregrouped(runs.xy => first => "X (cm)", runs.xy => last => "Y (cm)", layout = runs.run_id => nonnumeric) * visual(Lines)) |> draw(; axis = (; width = 400, height = 400, limits = ((-MAX_TRAJECTORY_LENGTH, MAX_TRAJECTORY_LENGTH), (-MAX_TRAJECTORY_LENGTH, MAX_TRAJECTORY_LENGTH))));
save_figure(fig, output, "overview_shift")

############ Figure 1a-c: Sample trajectories (10 per condition)

n = 10
df = @chain runs begin
    @subset norm.(last.(:cropped)) .≈ MAX_TRAJECTORY_LENGTH
    @groupby :dance_by
    @transform :n = 1:length(:dance_by)
    @subset :n .≤ n
end
# @assert all(==(n), combine(groupby(df, :dance_by), nrow).nrow)

fig = pregrouped(df.cropped => first => "X (cm)", df.cropped => last => "Y (cm)", col = df.dance_by) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content
    if ax isa Axis
        for r  in (30, MAX_TRAJECTORY_LENGTH)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save_figure(fig, output, "figure1a")

fig = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)", col = df.dance_by) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save_figure(fig, output, "figure1b")

df2 = stack(select(df, [:cropped, :dance_by, :rotated2poi]), [:cropped, :rotated2poi])
fig = pregrouped(df2.value => first => "X (cm)", df2.value => last => "Y (cm)", col = df2.dance_by, row = df2.variable => renamer("cropped" => "unrotated", "rotated2poi" => "rotated")) * visual(Lines) |> draw(; axis = (; width = 400, height = 400))
for ax in fig.figure.content
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end

save_figure(fig, output, "figure1c")

# ============================================================================
# SECTION 5: STATISTICAL ANALYSIS - TURNING BEHAVIOR
# ============================================================================
# Analyze the logistic curve parameters that characterize turning behavior:
#   - L: Total turn angle (radians)
#   - k: Steepness of turn (rate parameter)
#   - x₀: Distance from POI where turn begins (cm)
#   - Δl: Distance metric (x₀ relative to POI)
#   - Δθ: Total turn wrapped to [-π, π]
# ============================================================================

############ Prepare data: Filter for high-quality logistic fits (R² > 0.95)


df = @chain runs begin
    @subset :logistic_rsquare .> 0.95  # Keep only high-quality logistic fits
    @transform :Δθ = wrap2pi.(:L)  # Total turn wrapped to [-π, π]
    transform(:L => ByRow(l -> sincos(l .+ π/2)) => [:v, :u])  # Convert to vector components for arrows
    @rtransform :Δl = (:x₀ .- :θ.l[Ti = Near(:poi)])  # Distance from POI to turn inflection point
    @transform :absk = abs.(:k)  # Absolute steepness (rate of turn)
end

# fig = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v, col = :dance_by, row = :at_run => nonnumeric, color = :Δθ => rad2deg => "Total turn") * visual(Arrows) |> draw(scales(Color = (; colormap = :cyclic_wrwbw_40_90_c42_n256_s25, colorrange = (-180, 180))); axis = (; width = 200, height = 200), colorbar = (; ticks = -180:90:180, tickformat = "{:n}°"))

############ Figure 5: Arrow plot summarizing turn parameters


plt = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v, row = :location, col = :dance_by, color = :Δθ => rad2deg => "Total turn") * visual(Arrows2D)
fig = draw(plt, scales(Color = (; colormap = :cyclic_wrwbw_40_90_c42_n256_s25, colorrange = (-180, 180))); axis = (; width = 200, height = 200))

save_figure(fig, output, "figure5")


fig = Figure()
ax = PolarAxis(fig[1,1])
ts = []
ls = []
for (k, g) in pairs(groupby(df, :dance_by))
    _t = scatter!(g.L, g.absk, label = k.dance_by)
    push!(ts, _t)
    push!(ls, k.dance_by)
end
Legend(fig[1, 2], ts, ls)
save_figure(fig, output, "figure5a")


df2 = select(df, [:at_run, :dance_by, :L])
@transform! df2 :L = round.(Int, rad2deg.(:L .+ π/2))
CSV.write("L.csv", df2)

############ Compute normalized metrics for statistical modeling
# # PCA
@chain df begin
    @transform! :Δθnormalized = (π .- abs.(:Δθ)) ./ π  # Normalize turn angle to [0,1]
    @transform! :kabs = abs.(:k)  # Absolute steepness
    @transform! :Δlabs = abs.(:Δl)  # Absolute distance from POI
end

# Xtr = Matrix(select(df, :Δθnormalized, :kabs, :Δlabs))
# M = MultivariateStats.fit(PCA, Xtr; maxoutdim=3)
# data((; at_run = df.at_run, dance_by = df.dance_by, pc1 = first(eachcol(M.proj)), pc2 = last(eachcol(M.proj)))) * mapping(:pc1, :pc2, col = :dance_by, row = :at_run => nonnumeric) * visual(Scatter) |> draw()

############ Test location effect (LED vs Sun stimulus)

df2 = @chain df begin
    @subset :at_run .== 1  # First run only (to avoid repeated measures)
    # @rtransform :dance_by = :dance_by == "no" ? "no" : "yes"
end

############ Figure 7: Rainclouds plot comparing turn steepness across locations

# TODO: fix this
plt = data(df2) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v, col = :dance_by, row = :at_run => nonnumeric, color = :Δθ => rad2deg => "Total turn") * visual(Arrows2D)
fig = draw(plt, scales(Color = (; colormap = :cyclic_wrwbw_40_90_c42_n256_s25, colorrange = (-180, 180))); axis = (; width = 200, height = 200))


# fig = data(df2) * mapping(:location, :absk => "k", col = :dance_by) * visual(RainClouds, violin_limits = extrema) |> draw(; axis = (; width = 200, height = 200))
# save_figure(fig, output, "figure7")

# ============================================================================
# STATISTICAL MODELS: Test hypotheses about turning behavior
# ============================================================================

#### Model 1: Does turn steepness differ by location × dance type?
m = glm(@formula(kabs ~ location*dance_by), df2, Gamma())

#### Correlation tests: Are turn parameters independent?
cor(df.Δlabs, df.kabs)  # Distance vs steepness
cor(df.Δlabs, df.Δθnormalized)  # Distance vs total turn
cor(df.kabs, df.Δθnormalized)  # Steepness vs total turn

#### Model 2: Does turn steepness vary with dance type and run number?
m = glm(@formula(kabs ~ dance_by*at_run), df, Gamma())
m = glm(@formula(kabs ~ dance_by + at_run), df, Gamma())

#### Model 3: Simplify dance_by to binary (induced vs not)
@rtransform! df :dance_by = :dance_by == "no" ? "no" : "yes"
m = glm(@formula(kabs ~ location*dance_by), df, Gamma())

#### Model 4: Location + dance type (first run only)
m = glm(@formula(kabs ~ location+dance_by), @subset(df, :at_run .== 1), Gamma())

# predictions = DataFrame(dance_by = [true, false])
# predictions = hcat(predictions, rename(predict(m, predictions; interval = :confidence), :prediction => :k, :lower => :lower_k, :upper => :upper_k))

#### Model 5: Does distance from POI vary with dance type?
# m = glm(@formula(Δlabs ~ dance_by*at_run), df, Gamma())
# m = glm(@formula(Δlabs ~ dance_by + at_run), df, Gamma())
m = glm(@formula(Δlabs ~ dance_by), df, Gamma())

# predictions = hcat(predictions, rename(predict(m, predictions; interval = :confidence), :prediction => :l, :lower => :lower_l, :upper => :upper_l))

#### Model 6: Beta regression on normalized turn angle (bounded 0-1)
# m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_by*at_run), df)
# m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_by + at_run), df)
m = BetaRegression.fit(BetaRegressionModel, @formula(Δθnormalized ~ dance_by), df)

# predict(m, predictions)

# ============================================================================
# SECTION 6: FIGURE GENERATION - TURN ANGLE PROFILES
# ============================================================================

#### Prepare turn angle profiles: normalize so POI is at distance=0, angle=0
@chain df begin
    @rtransform! :lθshifted = parent(:θ.l .- :θ.l[Ti = Near(:poi)]) # Distance relative to POI
    @rtransform! :θnormalized = parent(:θ.θ .- :θ.θ[Ti = Near(:poi)])  # Angle relative to POI
end

############ Figure 6: Turn angle profiles over path length
fig = pregrouped(df.lθshifted => "Distance from POI (path length cm)", df.θnormalized => rad2deg) * visual(Lines) * mapping(color = df.location, col = df.dance_by, row = df.at_run => nonnumeric) |> draw()#; axis = (; width = 400, height = 400, ytickformat = "{:n}°", yticks = -180:90:270, limits = ((-5, 5), nothing)))

save_figure(fig, output, "figure6")

# ============================================================================
# SECTION 7: MAIN PUBLICATION FIGURES
# ============================================================================
# Generate multi-panel figures combining trajectory plots, turn profiles,
# and parameter summaries for first-run data only (at_run == 1)
# ============================================================================

#### Prepare dataset: First run only, high-quality fits
df = @chain runs begin
    @subset :at_run .== 1  # First run to avoid repeated measures
    @subset :logistic_rsquare .> 0.95  # High-quality fits only
    @transform :induced = :dance_by .≠ "no"  # Binary: induced vs spontaneous
    @transform :Δθ = wrap2pi.(:L)
    transform(:L => ByRow(l -> sincos(l .+ π/2)) => [:v, :u])
    @rtransform :Δl = (:x₀ .- :θ.l[Ti = Near(:poi)])
    @transform :absk = abs.(:k)
    @rtransform :lθshifted = parent(:θ.l .- :θ.l[Ti = Near(:poi)])
    @rtransform :θnormalized = parent(:θ.θ .- :θ.θ[Ti = Near(:poi)])
end

# fig = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v, col = :dance_by, row = :at_run => nonnumeric, color = :Δθ => rad2deg => "Total turn") * visual(Arrows) |> draw(scales(Color = (; colormap = :cyclic_wrwbw_40_90_c42_n256_s25, colorrange = (-180, 180))); axis = (; width = 200, height = 200), colorbar = (; ticks = -180:90:180, tickformat = "{:n}°"))

############ Figure 1: Multi-panel summary (trajectories + angles + arrows)

# Panel A: Spatial trajectories rotated so POI is north
m = mapping(col = df.induced => renamer(false => "Dance not induced", true => "Dance induced"), color = df.location => renamer("Lund" => "LED", "Bela-Bela" => "Sun") => "Stimulus")
tracks = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)") * visual(Lines)
# Panel B: Accumulated turn angle vs distance from POI
angles = pregrouped(df.lθshifted => "Distance from POI (path length cm)", df.θnormalized => rad2deg => "Accumulated turn (°)") * visual(Lines)
# Panel C: Arrow summary of turn parameters
m3 = mapping(col = :induced => renamer(false => "not", true => "induced"), color = :location => renamer("Lund" => "LED", "Bela-Bela" => "Sun") => "Stimulus")
plt = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v) * visual(Arrows2D)

fig = Figure()
h = draw!(fig[2,1], m * tracks; axis = (; width = 250, height = 250))
# Add reference circles at 30 cm and 50 cm
for ax in fig.content
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
draw!(fig[3,1], m * angles; axis = (; xticklabelsvisible = false, xticksvisible = false, xlabelvisible = false, xlabel = "", width = 250, height = 250, limits = ((-7, 17), (-200, 200)), yticks = -180:90:180))
draw!(fig[4,1], m3 * plt; axis = (; width = 250, height = 250, limits = ((-7, 17), nothing)))
legend!(fig[1, 1], h, orientation = :horizontal, titleposition = :left)
resize_to_layout!(fig)

save_figure(fig, output, "figure1")

############ Figure 2: Alternative layout grouping by location instead of induction
m = mapping(color = df.induced => renamer(false => "Dance not induced", true => "Dance induced"), col = df.location => renamer("Lund" => "LED", "Bela-Bela" => "Sun") => "Stimulus")
tracks = pregrouped(df.rotated2poi => first => "X (cm)", df.rotated2poi => last => "Y (cm)") * visual(Lines)
angles = pregrouped(df.lθshifted => "Distance from POI (path length cm)", df.θnormalized => rad2deg => "Accumulated turn (°)") * visual(Lines)
m3 = mapping(color = :induced, col = :location)
plt = data(df) * mapping(:Δl => "Distance from POI (path length cm)", :absk => "k", :u, :v) * visual(Arrows2D)
fig = Figure()
h = draw!(fig[1,1], m * tracks; axis = (; width = 250, height = 250))
for ax in fig.content
    if ax isa Axis
        for r  in (30, 50)
            lines!(ax, Circle(zero(Point2f), r), color=:gray, linewidth = 0.5)
        end
    end
end
draw!(fig[2,1], m * angles * mapping(col = df.location => "Stimulus"); axis = (; xticklabelsvisible = false, xticksvisible = false, xlabelvisible = false, xlabel = "", width = 250, height = 250, limits = ((-7, 17), (-200, 200)), yticks = -180:90:180))
draw!(fig[3,1], m3 * plt; axis = (; width = 250, height = 250, limits = ((-7, 17), nothing)))
legend!(fig[0, 1], h, orientation = :horizontal, titleposition = :left)
resize_to_layout!(fig)

save_figure(fig, output, "figure2")
