# ============================================================================
# GENERAL TRAJECTORY PROCESSING FUNCTIONS
# ============================================================================
# Helper functions for cleaning, transforming, and standardizing movement
# trajectories in animal tracking experiments.
# ============================================================================

const SV = SVector{2, Float64}

# ============================================================================
# TRAJECTORY CLEANING AND PROCESSING FUNCTIONS
# ============================================================================

###

"""
    remove_loops(xy)

Remove self-intersecting loops from trajectory by detecting intersection points
and removing all points within the loop. This cleans trajectories where the
animal retraces its path or makes tight circles.

Returns trajectory with loop segments removed (keeps loops > 50 points).
"""
function remove_loops(xy)
    inds, _ = self_intersections(Point2f.(xy))
    rngs = splat(UnitRange).(Iterators.partition(inds, 2))
    filter!(<(51) ∘ length, rngs)
    tokill = vcat(rngs...)
    unique!(tokill)
    sort!(tokill)
    keep::Vector{Int} = collect(eachindex(xy))
    setdiff!(keep, tokill)
    return xy[keep]
end

"""
    glue_intervention!(xy, intervention::Real)

Detect and correct spatial discontinuities caused by experimental interventions
(e.g., sudden light shifts that cause brief disorientation). Identifies sudden
jumps in position by comparing movement speeds before and after the intervention
time (using mean + 2σ threshold), then shifts the post-intervention trajectory
to reconnect it smoothly.

Returns the magnitude of the detected jump (0.0 if no jump detected).
Modifies xy in place.
"""
function glue_intervention!(xy, intervention::Real)
    h = 1
    diffs = diff(xy[Ti = 0.0 .. intervention + h])
    Δs = norm.(diffs)
    μ = mean(Δs[Ti = 0.0 .. intervention - h])
    σ = std(Δs[Ti = 0.0 .. intervention - h], mean = μ)
    threshold = μ + 2σ
    v = Δs[Ti = intervention - h .. intervention + h]
    i = findfirst(>(threshold), v)
    if isnothing(i)
        return 0.0
    else
        t1 = lookup(v, Ti)[i]
        step = diffs[Ti = At(t1)] - μ*normalize(diffs[Ti = At(t1)])
        t2 = lookup(v, Ti)[i + 1]
        xy[Ti = t2 .. Inf] .-= Ref(step)
        return v[i]
    end
end

"""
    sparseify(xy)

Resample trajectory at uniform 0.5 second intervals using parametric splines.
This converts irregularly sampled tracking data into evenly spaced points,
making temporal analyses more straightforward.
"""
function sparseify(xy)
    spl = ParametricSpline(lookup(xy, Ti), stack(xy))
    # tl = round(Int, t[1]):round(Int, t[end])
    tl = range(extrema(lookup(xy, Ti))..., step = 0.5)
    # poi_index = findfirst(==(round(Int, poi)), tl)
    DimVector(SV.(spl.(tl)), Ti(tl))
end

"""
    has_stops(ij)

Check if trajectory contains consecutive duplicate positions (stops).
Returns true if any two consecutive positions are identical.
"""
function has_stops(ij)
    for i in 2:length(ij)
        if ij[i] == ij[i - 1]
            return true
        end
    end
    return false
end

"""
    remove_stops(ij)

Remove consecutive duplicate positions from trajectory. This filters out
frames where the animal was stationary, keeping only frames with movement.
"""
function remove_stops(ij)
    keep = [1]
    for i in 2:length(ij)
        if ij[i] ≠ ij[i - 1]
            push!(keep, i)
        end
    end
    ij[keep]
end

"""
    get_tij(file)

Load tracking data from CSV file. Reads time (t), pixel row (i), and pixel
column (j) coordinates. Returns DimVector with time dimension.
"""
function get_tij(file)
    tij = CSV.File(file)
    t = range(tij.t[1], tij.t[end], length = length(tij))
    pixels = DimVector(CameraCalibrations.RowCol.(tij.i, tij.j), Ti(t))
end

"""
    smooth(xy)

Apply smoothing spline to trajectory to reduce noise while preserving shape.
Uses cubic spline (k=3) with smoothing parameter s=25. Resamples at higher
temporal resolution (100x denser) for smoother curves.
"""
function smooth(xy)
    t = lookup(xy, Ti)
    spl = ParametricSpline(t, stack(xy), k = 3, s = 25)
    tl = range(first(t), last(t), 100length(t))
    DimVector(SV.(spl.(tl)), Ti(tl))
end

# ============================================================================
# TRAJECTORY STANDARDIZATION FUNCTIONS
# ============================================================================
# Functions for normalizing trajectory position and orientation
# ============================================================================

"""
    center2start(xy)

Translate trajectory so starting position is at origin [0, 0].
This standardizes trajectories for comparison regardless of where they began.
"""
function center2start(xy)
    trans = Translation(-first(xy))
    trans.(xy)
end

"""
    intersection(orig, dir)

Find intersection of a ray with unit circle. Used internally by cropto()
to determine where trajectory crosses the radius threshold.

Returns tuple (t1, t2) of ray parameters for near and far intersections.
"""
function intersection(orig, dir)
    b = -orig⋅dir
    disc = b^2 - orig⋅orig + 1
    if disc ≥ 0
        d = sqrt(disc)
        t2 = b + d
        if t2 ≥ 0
            t1 = b - d
            return t1 > 0 ? (t1, t2) : (Inf, t2)
        end
    end
    return (Inf, Inf)
end

"""
    cropto(xy, l)

Crop trajectory at maximum distance l from origin. If trajectory extends beyond
radius l, it is truncated at exactly distance l along the path direction.
This ensures all trajectories have the same spatial scale for comparison.
"""
function cropto(xy, l)
    i = findfirst(>(l) ∘ norm, xy)
    isnothing(i) && return copy(xy)
    p1 = xy[i - 1]
    p2 = xy[i]
    dir = normalize(p2 - p1)
    _, d = intersection(p1/l, normalize(p2 - p1))
    p2 = p1 + d*l*dir
    cropped = xy[1:i]
    cropped[i] = p2
    return cropped
end

"""
    rotate2poi(xy, poi)

Rotate trajectory so the Point Of Interest (POI) is at 90° (north/up).
This aligns all trajectories to a common reference frame where the experimental
intervention or behavioral event occurs at the same angular position.
"""
function rotate2poi(xy, poi)
    p2 = xy[Ti = Near(poi)]
    θ = π/2 - atan(reverse(p2)...)
    rot = LinearMap(Angle2d(θ))
    rot.(xy)
end

"""
    center2poi_and_crop(xy, poi)

Center trajectory on POI position and return only the post-POI segment.
This isolates the behavioral response after the experimental intervention.
"""
function center2poi_and_crop(xy, poi)
    trans = Translation(-xy[Ti = Near(poi)])
    trans.(xy[Ti = poi..Inf])
end

# ============================================================================
# CAMERA CALIBRATION AND COORDINATE TRANSFORMATION
# ============================================================================

###

# function trectify(fs, xs)
#     n = length(xs)
#     ys = Vector{DimVector{SV}}.(undef, n)
#     Threads.@threads for i in 1:n
#         y = fs[i].(xs[i])
#         y .-= Ref(y[1])
#         ys[i] = y
#     end
#     return ys
# end

"""
    get_calibration(file)

Load camera calibration from file and return rectification function.
The returned function maps pixel coordinates (i,j) to real-world coordinates (x,y)
in centimeters, correcting for lens distortion and camera perspective.
"""
function get_calibration(file)
    c = CameraCalibrations.load(file)
    f = rectification(c, findfirst(contains("extrinsic"), c.files))
    # fun(ij)::Function = f.(ij)
    return f
end

# ============================================================================
# TIME CONVERSION UTILITIES
# ============================================================================

"""
    tosecond(t)

Convert various time representations to seconds (Float64).
Handles Time objects, TimePeriod, and numeric seconds.
"""
tosecond(t::T) where {T <: TimePeriod} = t / convert(T, Dates.Second(1))
tosecond(t::TimeType) = tosecond(t - Time(0))
tosecond(sec::Real) = sec

# function get_txy(tij_file, rectify)
#     tij = CSV.File(joinpath(results_dir, tij_file))
#     t = range(tij.t[1], tij.t[end], length = length(tij))
#     ij = SVector{2, Int}.(tij.i, tij.j)
#     xy = rectify.(ij)
#     (; t, xy)
# end

# function remove_stops!(t, xy)
#     tokill = Int[]
#     last_xy = xy[1]
#     for i in 2:length(t)
#         Δ = norm(xy[i] - last_xy)
#         if Δ > 0.1
#             last_xy = xy[i]
#         else
#             push!(tokill, i)
#         end
#     end
#     deleteat!(t, tokill)
#     deleteat!(xy, tokill)
#     return (; t, xy)
# end
#

"""
    impute_poi_time(xy)

Estimate Point Of Interest time by detecting when animal has moved more than
10 cm from starting position. Used when POI time is not explicitly recorded.
"""
function impute_poi_time(xy)
    p1 = xy[1]
    for i in eachindex(xy)
        if norm(xy[i] - p1) > 10
            return lookup(xy, Ti)[i - 1]
        end
    end
    return last(lookup(xy, Ti))
end

# ============================================================================
# ANALYSIS PIPELINE HELPER FUNCTIONS
# ============================================================================
# High-level functions for common analysis workflows
# ============================================================================

# Trajectory processing parameters
const SMOOTHING_K = 3        # Spline degree
const SMOOTHING_S = 25       # Smoothing parameter
const RESAMPLE_INTERVAL = 0.5 # Seconds
const MAX_TRAJECTORY_LENGTH = 50  # Centimeters

"""
    setup_output_dir(name)

Create or clear output directory for analysis results.
If directory exists, removes all contents. If not, creates it.

Returns the directory name.
"""
function setup_output_dir(name)
    if isdir(name)
        rm.(readdir(name; join = true))
    else
        mkdir(name)
    end
    return name
end

"""
    save_figure(fig, output_dir, filename)

Save figure in both PNG (GLMakie) and PDF (CairoMakie) formats.
Filename should not include extension - both .png and .pdf will be created.

Returns the figure object.
"""
function save_figure(fig, output_dir, filename)
    GLMakie.activate!()
    save(joinpath(output_dir, "$filename.png"), fig)
    CairoMakie.activate!()
    save(joinpath(output_dir, "$filename.pdf"), fig)
    return fig
end

"""
    load_runs_and_calibs(results_dir; filter_light=nothing, exclude_spontaneous=false)

Load tracking data and calibrations, joining them together.

Parameters:
- results_dir: Path to "tracks and calibrations" directory
- filter_light: If specified, filter for this light condition (e.g., "shift", "dark", "remain")
- exclude_spontaneous: If true, exclude spontaneous dances (keep only induced)

Returns DataFrame with runs and calibration data joined.
"""
function load_runs_and_calibs(results_dir; filter_light=nothing, exclude_spontaneous=false)
    runs = @chain joinpath(results_dir, "runs.csv") begin
        CSV.read(DataFrame)
        @select Not(:runs_path, :start_location, :fps, :target_width, :runs_file, :window_size)
        @transform :tij_file = joinpath.(results_dir, :tij_file)
    end

    if !isnothing(filter_light)
        @subset!(runs, :light .== filter_light)
    end

    if exclude_spontaneous
        @subset!(runs, ismissing.(:spontaneous_end))
    end

    calibs = @chain joinpath(results_dir, "calibs.csv") begin
        CSV.read(DataFrame)
        @transform :calibration_file = joinpath.(results_dir, :calibration_id)
        @transform :rectify = get_calibration.(:calibration_file)
        @select Cols(:calibration_id, :rectify)
    end

    leftjoin!(runs, calibs, on = :calibration_id)
    @select! runs Not(:calibration_id)
    return runs
end

"""
    process_trajectories!(df;
                          fix_intervention_jumps=true,
                          impute_poi_for_remain=false,
                          handle_spontaneous=false,
                          l=MAX_TRAJECTORY_LENGTH)

Apply standard trajectory processing pipeline in-place.

Processing steps:
1. Load pixel coordinates and convert to real-world (cm)
2. Remove stops (consecutive duplicate positions)
3. Fix discontinuities at intervention time (optional)
4. Impute POI time for "remain" condition (optional)
5. Handle spontaneous dances (optional)
6. Remove self-intersecting loops
7. Resample at uniform intervals
8. Smooth with spline
9. Center, crop, and rotate trajectories

Parameters:
- df: DataFrame with tracking data (modified in place)
- fix_intervention_jumps: Detect and correct spatial jumps at intervention
- impute_poi_for_remain: For "remain" light condition, estimate POI time
- handle_spontaneous: Process spontaneous dance timing
- l: Maximum trajectory length in cm

Returns the modified DataFrame.
"""
function process_trajectories!(df;
        fix_intervention_jumps=true,
        impute_poi_for_remain=false,
        handle_spontaneous=false,
        l=MAX_TRAJECTORY_LENGTH)
    @chain df begin
        @transform! :condition = string.(:light, " ", :dance_by, " ", :at_run)
        @rtransform! :y2025 = Year(:start_datetime) == Year(2025) ? "2025" : "earlier"
        @rename! :intervention = :poi

        # Load and clean pixel data
        @transform! :pixels = get_tij.(:tij_file)
        @transform! :pixels = remove_stops.(:pixels)
        @rtransform! :xy = :rectify.(:pixels)
    end

    # Optional: fix intervention discontinuities
    if fix_intervention_jumps
        @chain df begin
            @subset(:dance_by .≠ "no"; view = true)
            @transform! :jump = glue_intervention!.(:xy, :intervention)
        end
    end

    # Optional: impute POI for "remain" condition
    if impute_poi_for_remain
        @chain df begin
            @subset(:light .== "remain"; view = true)
            @transform! :intervention = impute_poi_time.(:xy)
        end
    end

    # Optional: handle spontaneous dances
    if handle_spontaneous
        @transform! df :spontaneous_end = passmissing(tosecond).(:spontaneous_end)
        @transform! df :poi = coalesce.(:spontaneous_end, :intervention)
    else
        @transform! df :poi = :intervention
    end

    @chain df begin
        disallowmissing!(:poi)
        @select! Not(:intervention)

        # Clean and standardize trajectories
        @transform! :xy = remove_loops.(:xy)
        @transform! :xy = sparseify.(:xy)
        @transform! :smooth = smooth.(:xy)
        @transform! :centered2start = center2start.(:smooth)
        @transform! :cropped = cropto.(:centered2start, l)
        @transform! :rotated2poi = rotate2poi.(:cropped, :poi)
        @transform! :centered2poi_and_cropped = center2poi_and_crop.(:rotated2poi, :poi)

    end
end

