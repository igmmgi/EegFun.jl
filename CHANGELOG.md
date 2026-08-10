# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to (well, attempts to :-)) [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.1] - 2026-08-10

### Fixed

- Fixed `DimensionMismatch` bug in `rereference` when using array-based reference selections (like `[:M1, :M2]`).
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
