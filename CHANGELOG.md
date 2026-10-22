# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to (well, attempts to :-)) [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Closed-form Freedman–Lane permutation engine in `fit_mass_lmm`, parallelized via `Distributed`.
- Convergence tracking and full parametric statistics (`z_values`, `p_parametric`, `fit_status`) in `LmmStatsResult`.
- Support for factor contrasts (e.g., `EffectsCoding`) in `fit_mass_lmm`.
- Permutation memory guard and sample-size validation in `fit_mass_lmm`.
- Memory leak fixes and performance optimizations in `plot_databrowser` (cached ICA activations, coordinate buffers).
- Complete `Base` array and `Tables.jl` interfaces for `EegData` types.
- Canonical accessors (`time_vector`, `frequency_vector`, etc.) and `show` methods across stats types.
- Expanded unit test suites for statistics, decoding, and signal processing.

### Fixed

- LLM-guided refactor and fixes across the ICA pipeline (`infomax.jl`, `picard.jl`, `amica.jl`, `artifact_id.jl`, `components.jl`):
  - **Algorithm Corrections**: Aligned Extended Infomax natural gradient with Lee et al. (1999) / EEGLAB `runica.m` for correct sub-Gaussian separation; fixed Extended Picard line-search loss evaluation and stabilized the block Hessian solve; added parameter snapshot rollback, non-divergent scale updates, and kurtosis-based sub-Gaussian detection to AMICA.
  - **Component Subtraction & Restoration**: Fixed rank-deficient component restoration in `add_ica_components!` ($X + s A_c^T M_c^T$) preserving variance orthogonal to the ICA subspace; guarded against duplicate indices and re-subtraction overwrite; added `EpochData` and `Vector{DataFrame}` support with arbitrary numeric column types (`Float32`, `Float64`); optimized memory layout in `samples × channels`, eliminating `permutedims` copies and ~75% peak memory.
  - **Artifact Identification**: Implemented top-level keyword argument routing; robust MAD-scaled z-scoring and directional thresholding; ECG polarity invariance and continuous sample coordinate mapping for IBIs; Welch's periodogram for line noise estimation.
  - **Architecture & Performance**: Unified input preprocessing/postprocessing (`_ica_preprocess`, `_ica_postprocess`) with rank-deficiency detection, input copy protection, and activation power variance sorting; accelerated CPU/GPU execution via BLAS `mul!`, register accumulation, and tensor preallocations.
- State consistency, observer lifecycle, and pipeline ordering in `plot_databrowser`.
- Isolated flat channels, bad pixels, and non-converged fits as `NaN` in `fit_mass_lmm`, decoupled fit status from permutation inference eligibility, initialized shared outputs to `NaN`, and added valid-field tracking.
- Eliminated Windows `mmap: Access is denied` failures by restricting `SharedArray` allocations strictly to multi-process distributed execution (`nprocs() > 1`), using standard Julia `Array`s and zero-copy inputs in single-process mode.
- Ensured worker processes inherit package directory in `LOAD_PATH` with `--startup-file=no` during distributed mass-univariate LMM inference.
- Channel line coloring under ICA overlay in `plot_databrowser`, defaulting to showing only original data in standard channel color, with dynamic switching to green upon enabling "Data - ICA", and fixing callback `MethodError` on toggle.
- `BoundsError` and trial power averaging across time-frequency pipelines.
- `MethodError` in out-of-place `channel_repairable`.
- CI documentation build failure by pinning readabilities plugin.
- Respected `display_plot = false` in `plot_tf`.
- Guarded error logging across GUI plotting routines.

### Changed

- Bumped minimum supported Julia version to 1.13 across `Project.toml`, GitHub Actions CI, documentation, and badges.
- Replaced iterative LMM permutation refits with closed-form Woodbury GLS for instant permutation inference.
- Modularized repository source structure (`ica/`, `utils/data/`, `artifact_detection/`).
- Exported and documented missing functions in `public_api.jl`, standardizing on `@minimal_error`.
- Updated to `AnovaFun` 0.5 with explicit Makie backend loading.
- Optimized allocations across time-frequency, CSD, permutations, decoding, and GPU ICA.
- Accelerated topography rendering with batched interpolation and observable binding.

### Removed

- Deprecated mutating batch string aliases (`realign!`, `calculate_trigger_interval!`).
- Legacy compatibility overloads (old `Layout` JLD2 conversion, obsolete decoding layouts).

## [0.8.0] - 2026-08-13

### Added

- Implemented `Base` Array interface, conversions, iteration, and `vcat` methods for EEG data types.
- Implemented `StatsBase.mean` methods for EEG data types.
- Implemented `Tables.jl` interface for EEG data types.

### Changed

- Refactored internal code to replace direct data field access with `n_epochs` and `n_samples` helper functions.
- Consolidated layout, legend, and colorbar parameters into kwargs named tuples.
- Improved theme kwarg integration and removed the default jet colormap.
- Updated documentation and README to include Makie backend requirements.
- Refactor of the plotting API to standardize title arguments across all plots: `window_title` controls the OS window, `figure_title` (and `figure_title_fontsize`) controls the main super-title, and `plot_title` (and `plot_title_fontsize`) controls individual subplot titles.
- Added `plot_title_position` and `plot_title_align` to allow rendering subplot titles *inside* the axis (e.g., `plot_title_position=(0.5, 0.95)`), which saves vertical space for dense EEG channel grids.
- Increased the default distance criterion for spatial layout neighbours (`neighbour_criterion`) from `0.25` to `0.35` to better support standard 10-20/10-10 electrode layouts and prevent empty connectivity matrices during cluster permutation tests.
- Formatted the entire codebase using `JuliaFormatter.jl` to adhere strictly to `.JuliaFormatter.toml` styling rules.
- 100% of the public API is now fully documented with comprehensive docstrings.

### Fixed

- Fixed an issue where `EegFun.conditions(...)` array subsets (like `conditions([3, 1])`) would discard the requested array order and always extract conditions numerically/alphabetically. This ensures paired statistics (`prepare_stats`) and plotting functions accurately reflect the requested data order.
- Fixed messy legends in `plot_erp_stats` by hiding redundant labels for channels.
- Added a warning in `plot_erp_stats` when `average_channels=true` is passed, as spatial averaging is statistically invalid for cluster permutation statistics.
- Fixed `DimensionMismatch` bug in `rereference` when using array-based reference selections (like `[:M1, :M2]`).
- Fixed default arguments and added a missing dispatch for `plot_erp_errorbar(dat::ErpData)`.
- Removed invalid keyword arguments in `apply_analysis_settings` and `apply_analysis_settings!`.
- Cleaned up unused variables and arguments across `erp_measurements.jl` and `jackknife_average.jl`.

## [0.7.0] - 2026-08-07

### Changed

- Made `EegFun.jl` plotting backend-agnostic by removing direct `GLMakie` and `CairoMakie` dependencies. Users must now explicitly load their preferred Makie backend (e.g., `using GLMakie`) alongside `EegFun.jl`.
- Updated all documentation, tutorials, and quickstart examples to reflect the new `using GLMakie` requirement for interactive plotting.
- Transitioned filter configurations and keyword arguments to use Julia `Symbol`s instead of strings for better performance and consistency.
- Minor code tidy and formatting improvements.

## [0.6.0] - 2026-08-03

### Added

- Integrated Native AMICA algorithm with GPU support.
- Added GPU support for time-frequency (TF) analysis and statistics.
- Added `check_raw_data` utility function.
- Added validation to ensure the data file contains all channels required by the layout.

### Changed

- Changed saving default flags
- Improved and optimized PICARD ICA implementation, including GPU support.
- Updated tutorial code examples and documentation formatting.

## [0.5.0] - 2026-07-22

### Added

- Added `skip_existing` flag to the `preprocess` pipeline, allowing users to safely recover and resume interrupted batch jobs without reprocessing files.
- Added `dry_run` mode to the `preprocess` pipeline to validate configurations and file paths before executing time-consuming processing.
- Added elapsed time and dynamic ETA tracking to the batch processing loop to improve the quality of life for large dataset analyses.
- Added an `ARCHITECTURE.md` guide outlining the directory structure and pipeline data flow to aid new contributor onboarding.
- Ported the Chronux/EEGLAB adaptive multi-taper line noise removal algorithm (`cleanline!`).
- Integrated `cleanline!` into the automated processing pipeline (`pipeline.jl`).
- Added multi-threading support and DC spectral leakage protection (window mean-centering) for optimal `cleanline!` performance on un-filtered raw data.
- Added `[preprocess.resample]` block to the automated `pipeline.toml` configuration to allow downsampling data early in the pipeline, reducing memory usage and speeding up ICA.
- Added `[preprocess.eeg.artifact_interval_start]` and `[preprocess.eeg.artifact_interval_end]` to the `pipeline.toml` to allow users to define a smaller time-window for artifact rejection within an oversized epoch.
- Ported the Residue Iteration Decomposition (RIDE) algorithm for latency-variable ERP component separation (`compute_ride`).
- Implemented Woody filtering (`_woody_filter`) and component iterative updating for RIDE using a zero-allocation windowing approach.
- Implemented Current Source Density (CSD) / Surface Laplacian (`compute_csd!`) using Perrin et al. (1989) spherical splines, verified for parity against MNE-Python.
- zero-allocation virtual padding in time-frequency pipelines (`tf_morlet`, `tf_stft`, `tf_multitaper`) to reduce memory usage and improve speed.
- Added interactive "Topoplot 3D" plotting option to the Data Browser's selected region popup menu.
- Added keyboard interactivity (Up/Down arrow keys) to dynamically scale colormap limits in 3D topographic plots (`LScene`/`Mesh` support).

### Changed

- Consolidated the EOG detection API by replacing the redundant `detect_eog_signals!` wrapper with `detect_eog_onsets!` utilizing Julia's multiple dispatch for a cleaner interface.
- Renamed ICA component modification functions for clarity and mathematical intuition: `remove_ica_components` is now `subtract_ica_components`, and `restore_ica_components` is now `add_ica_components`.
- Refactored `interval_selection` predicate logic to use type-stable `AbstractSelection` types (`AllSelection`, `TimeSelection`) instead of returning implicit `nothing` or `Tuple` values via `times()`.
- Decoupled `baseline_interval` parameters from the `times()` default across all plotting functions (`plot_erp`, `plot_topography`, `plot_tf`, etc.), standardizing `nothing` as the explicit default for "no baseline correction."
- `resample` API has transitioned to using an explicit `target_rate` argument (e.g., target sample rate in Hz) instead of a downsampling `factor`.

### Deprecated

### Removed

- Removed old `factor`-based resampling logic and corresponding test validations.

### Fixed

- Eliminated a silent participant ID collision risk in the batch pipeline's ID fallback generator by increasing the hash modulo from 10,000 to one billion.
- Fixed TOML serialization logic for nested array configurations (e.g., `vEOG_channels`), which previously generated invalid config templates.
- Resolved `MethodError` exceptions during artifact rejection by relaxing type constraints from `Vector` to `AbstractVector` across core analytical functions (e.g., `_calculate_epoch_metrics`, `channel_repairable!`), allowing zero-allocation subsets and ranges (like `UnitRange`) to be passed without crashes.
- Fixed TOML serialization crashes for `nothing` configuration values during `pipeline.toml` template generation, correctly generating commented-out defaults to ensure valid `.toml` output parsing.
- Fixed metadata preservation (`:time` and `condition` columns) edge cases during Time-Frequency data extraction.
- Fixed `plot_topography` contour scaling issues that previously resulted in white rendering artifacts for peak values outside of standard ranges.
