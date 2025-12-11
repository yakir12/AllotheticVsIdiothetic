# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Julia-based research codebase for analyzing animal tracking and movement data. The project focuses on tracking animals (beetles) in various experimental conditions, analyzing their movement trajectories, rotational behavior, and calibrating camera systems.

The codebase is organized into multiple experimental contexts (indoors, outdoors, ccw, rotation, etc.), each with its own tracking, calibration, and analysis pipelines.

## Key Dependencies

- **Fromage**: Custom package for video tracking and calibration (provides `main()`, `only_track()`, `only_calibrate()` functions)
- **PawsomeTracker**: Video-based animal tracking
- **CameraCalibrations**: Camera calibration and image rectification
- **DataFrames/DataFramesMeta**: Data manipulation with `@chain`, `@transform`, `@subset`, etc.
- **AlgebraOfGraphics/Makie**: Plotting (GLMakie for interactive, CairoMakie for publication-quality PDFs)
- **Dierckx**: Spline interpolation (ParametricSpline)
- **DimensionalData**: Arrays with named dimensions (DimVector with Ti time dimension)
- **Statistics packages**: GLM, MixedModels, HypothesisTests, BetaRegression

## Working with Julia Scripts

### Running Analysis Scripts

Each analysis directory (e.g., `outdoors/analyse/`, `indoors/analyse/`, `ccw/`) contains a `main.jl` file that performs the full analysis pipeline:

```bash
# Navigate to the analysis directory
cd outdoors/analyse

# Run the analysis in Julia REPL
julia --project=. main.jl
```

The scripts use `Revise` for hot-reloading during development and typically:
1. Load tracking data from CSV files (`runs.csv`, `calibs.csv`)
2. Apply transformations (smoothing, centering, rotation)
3. Fit statistical models
4. Generate figures (both PNG and PDF)

### Running Tracking/Calibration

The tracking and calibration scripts are in directories like `outdoors/track_calibrate/`, `outdoors/track/`, `outdoors/calibrate/`:

```bash
cd outdoors/track_calibrate
julia --project=. main.jl
```

These scripts:
- Call functions from the `Fromage` package (`main()`, `only_track()`, `only_calibrate()`)
- Process video files from mounted data directories (e.g., `/home/yakir/mnt/Data/Elin/...`)
- Generate output in `tracks and calibrations/` directories
- Clean up temporary `jl_*` directories before running

### Project-Specific Configuration

Each subdirectory has its own Julia environment:
- `Project.toml`: Lists package dependencies
- `Manifest.toml`: Pins exact package versions
- `LocalPreferences.toml`: Configuration for packages like Fromage (e.g., `target_width`, `temporal_step`, checker pattern parameters)

To activate a specific environment:
```bash
cd outdoors/track_calibrate
julia --project=.
```

## Codebase Architecture

### Directory Structure

- **`outdoors/`, `outdoors50/`, `indoors/`**: Main experimental datasets
  - `track_calibrate/`: Combined tracking and calibration pipeline
  - `track/`: Tracking only
  - `calibrate/`: Calibration only
  - `analyse/`: Data analysis and figure generation
- **`ccw/`**: Counter-clockwise rotation experiments
- **`rotation/`, `rotation2/`**: Rotation-specific analyses
- **`second/`**: Secondary analyses
- **Root directory**: Contains `backup.jl` and `Project.toml`

### Data Pipeline

1. **Tracking**: Extract pixel coordinates (i, j) and time (t) from videos → `tij` files
2. **Calibration**: Generate camera calibration using checkerboard patterns → `.calib` files
3. **Rectification**: Convert pixel coordinates to real-world coordinates (x, y in cm)
4. **Transformation**: Apply smoothing, centering, rotation transformations
5. **Analysis**: Fit statistical models, compute metrics, generate visualizations

### Common Data Structures

- **`DimVector`**: Time-indexed vectors from DimensionalData package
  - Example: `DimVector(points, Ti(times))`
  - Access with `xy[Ti = Near(poi)]` or `xy[Ti = 0.0 .. 10.0]`
- **`SV`**: Type alias for `SVector{2, Float64}` (2D static vectors)
- **DataFrame transformations**: Heavy use of `@chain`, `@transform`, `@rtransform`, `@subset`

### Key Functions

#### In `general_functions.jl`:
- `get_tij(file)`: Load tracking data from CSV
- `get_calibration(file)`: Load camera calibration and return rectification function
- `smooth(xy)`: Apply spline smoothing with parameters `k=3, s=25`
- `center2start(xy)`: Translate trajectory to start at origin
- `cropto(xy, l)`: Crop trajectory at distance `l` (typically 50 cm)
- `rotate2poi(xy, poi)`: Rotate trajectory so POI (point of interest) is at 90°
- `remove_loops(xy)`: Remove self-intersecting trajectory segments
- `sparseify(xy)`: Resample trajectory at 0.5s intervals
- `glue_intervention!(xy, intervention)`: Detect and correct discontinuities at intervention time

#### In `shift_functions.jl`:
- `fit_logistic(x, θ)`: Fit logistic curve to turning angle data
- `get_turn_profile(xy, poi)`: Extract angle profile around POI
- `unwrap!(x, period=2π)`: Unwrap angles to avoid 2π discontinuities

## Common Development Tasks

### Modifying Analysis Parameters

Key constants are often defined at the top of `main.jl` files:
```julia
const l = 50  # Maximum trajectory length in cm
```

Fromage parameters are in `Project.toml` under `[preferences.Fromage]`:
```toml
target_width = 60
temporal_step = 2.0
checker_size = 4
n_corners = "(5, 8)"
```

### Running Tests

Tests are inline in some analysis scripts (e.g., `outdoors/analyse/runtests.jl`):
```julia
# Run with @test macros from Test.jl
julia> include("runtests.jl")
```

### Data Paths

- Production data is typically mounted at `/home/yakir/mnt/Data/Elin/Project_*`
- Output directories are named `tracks and calibrations/`, `shift/`, etc.
- The scripts automatically create output directories if they don't exist

### Figure Generation

Figures are saved in both formats for different purposes:
```julia
GLMakie.activate!()      # For interactive viewing and PNG
save(joinpath(output, "figure1.png"), fig)

CairoMakie.activate!()   # For publication-quality PDFs
save(joinpath(output, "figure1.pdf"), fig)
```

## Code Style Notes

- Use `@chain` for data transformation pipelines
- Prefer `@rtransform` for row-wise operations, `@transform` for vectorized operations
- Functions ending in `!` modify their arguments (e.g., `glue_intervention!`, `unwrap!`)
- Angle convention: Typically work in radians, convert to degrees for plotting
- Time dimension: Use `Ti` dimension from DimensionalData
- Use `passmissing()` to propagate missing values through transformations
- Statistical models use GLM.jl (`glm()`) and MixedModels.jl (`glmm()`, `lmm()`)

## Common Workflows

### Adding a New Analysis

1. Navigate to the appropriate `analyse/` directory
2. Ensure `runs.csv` and `calibs.csv` exist in `../track_calibrate/tracks and calibrations/`
3. Load data with `@chain` pipeline
4. Apply transformations using helper functions
5. Fit statistical models
6. Generate figures with AlgebraOfGraphics

### Debugging Tracking Issues

1. Check temporary error logs in `jl_*/error.log` directories (usually cleaned up)
2. Verify data paths point to correct mounted directories
3. Ensure video files have correct extensions (`.MTS`)
4. Check calibration parameters in `LocalPreferences.toml`

### Processing New Videos

1. Place videos in data directory (e.g., `/home/yakir/mnt/Data/Elin/...`)
2. Update `todo` parameter in `main.jl` to specify which videos to process
3. Run `main(data_path, todo = ["video_file.MTS"])`
4. Output will be in `tracks and calibrations/` subdirectory
