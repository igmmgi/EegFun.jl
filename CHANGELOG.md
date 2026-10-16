# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to (well, attempts to :-)) [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

*Many Gemini guided refactor + optimisations.*

### Added

- Mixed models permutations (`fit_mass_lmm`).
- Complete `Base` array and `Tables.jl` interface for `EegData` (indexing, dimensions, equality, conversions, and shorthand `dat[:Cz]` channel subsetting).
- Thread-safe async progress tracking (`track_progress_async`) to eliminate progress bar lock contention during parallel decoding and permutation testing.
- Comprehensive unit tests for previously untested modules: CleanLine, CSD, RIDE, Welch frequency spectrum, decoding preparation, FieldTrip CSV import, and pipeline channel repairs.
- Comprehensive unit tests for independent-design statistics (ERP and TF), directional non-parametric permutation thresholding, and epoched artifact detection/repair (spherical spline interpolation, bridging detection, flatline, and peak-to-peak).

### Fixed

- Fixed `BoundsError` in `tf_morlet`, `tf_stft` and `tf_multitaper` when using `return_trials=true` with `filter_edges=true`.
- `tf_multitaper` now averages (not sums) power over trials when `return_trials=false`, matching `tf_morlet`/`tf_stft`. Values are now `n_trials` times smaller than before.
- Fixed `MethodError` in out-of-place `channel_repairable` by implementing `Base.copy(::ContinuousRepairInfo)`.
- Fixed documentation build failure on CI by pinning `@nolebase/vitepress-plugin-enhanced-readabilities` to `~2.18.2` (preventing upstream peer dependency conflict with VitePress 1.x).

### Changed

- Shared argument validation for `tf_morlet`, `tf_stft` and `tf_multitaper` (pad, frequencies, window length, time grid, GPU fallback, subsetting) now lives in `tf_utils.jl`. Error messages are unchanged.
- Remaining raw `@warn` calls now use `@minimal_warning` for consistent log output.
- Removed an unused `_create_menu` helper and renamed two internal helpers (`_prepare_ica_plot_matrix`, `_apply_browser_filter!`) that shared names with unrelated functions.
- `tf_morlet` is faster and allocates ~99% less (output unchanged).
- `tf_stft` and `tf_multitaper` are faster and allocate ~48–58% less memory through thread-local buffer reuse and type-stable processing loops (output unchanged).
- Modularized `ica.jl` (3,985 lines) into `src/analysis/processing/ica/` (`api.jl`, `infomax.jl`, `picard.jl`, `amica.jl`, `components.jl`, `artifact_id.jl`).
- Modularized `utils/data.jl` (2,520 lines) into `src/utils/data/`, `artifact_detection.jl` (1,961 lines) into `src/analysis/processing/artifact_detection/`, and `pipelines/utils/utils.jl` into `repair_info.jl`, `repair.jl` and `summaries.jl`. Pure file moves, no behaviour change.
- Removed rogue `export` in `decompose.jl`, aligning RIDE with EegFun's zero-namespace-pollution policy.
- Expanded `public_api.jl` to declare 45 previously omitted public functions, structs, and mutating `!` pairs.
- Updated to `AnovaFun` 0.5, which no longer loads CairoMakie. Load a Makie backend explicitly (`using GLMakie` or `using CairoMakie`) before plotting.
- `preprocess` is split into `_pipeline_setup`, `_pipeline_process_file!`, `_pipeline_final_summary` and `_pipeline_finalize_log` (same log output).
- `preprocess(...; dry_run = true)` no longer creates the output directory or writes `config.toml` and `neighbours.toml`.
- Consistent error handling: the ~70 remaining plain `error()` calls for user-facing problems now use `@minimal_error` and throw `EegFunError` instead of `ErrorException`. `@minimal_error` only logs while a log file is active, so interactive use no longer prints each error twice.
- Optimized tensor conversions (`Base.Matrix`, `Base.Array`), ERP array extraction, and RSA/decoding data preparation using typed function barriers, reducing allocations by up to 99.9%.
- `fit_mass_lmm(::EpochData, ...)` is ~10× faster and allocates ~94% less memory during data extraction.
- Cluster-based permutation testing (2D, 3D, and MixedModels) is ~12× faster with zero heap allocations per permutation.
- Accelerated permutation $t$-matrix generation by skipping unneeded standard error and $p$-value allocations (~10× faster).
- Improved tail probability precision in mass-univariate $t$-tests using `ccdf` to avoid underflow for $|t| > 8$.
- Channel correlation matrix computation is ~2× faster and uses ~40% less memory.
- Optimized 3D time-frequency tensor extraction and $p$-value calculation with zero-allocation slicing and SIMD vectorization.
- Eliminated splatting overhead in `Tables.jl` and data accessors using `reduce(vcat, ...)`.
- Added `Base.show` and canonical accessors (`time_vector`, `frequency_vector`, `channel_labels`, `n_channels`, `n_samples`) across `LmmStatsResult`, `Layout`, `TimeFreqData`, and statistical results.
- Accelerated continuous event window marking (`_mark_windows_at_indices!`) using binary search on timestamps (~6,250× faster).
- `average_epochs` averages channels in a single pass across epochs, improving cache locality (~1.6× faster).
- RSA pooled covariance pre-allocates condition trials and avoids transpositions; Mahalanobis distance uses 3-argument in-place `dot` (~1.3× faster).
- Replaced splatting into `vcat` with `reduce(vcat, ...)` across trigger processing and channel summaries.
- Epoched CSD processing (`_apply_csd!`) is parallelized with `Threads.@threads`; eliminated redundant matrix transpositions via direct in-place BLAS GEMM (~3× faster on continuous data).
- Accelerated moving-window peak-to-peak artifact detection (`_is_peak_to_peak`) using the Van Herk / Gil-Werman sliding min/max algorithm (~8× faster).
- Accelerated moving-window flatline detection (`_is_flatline`) using $O(1)$ sliding Welford variance with periodic drift correction (~13× faster).
- Pre-allocated sample masks in `_detect_extreme_values`, reducing redundant per-channel BitVector allocations.
- `find_bridged_channels` and `_correlation_matrix` use BLAS matrix correlation and typed column buffers, speeding up bridging detection by ~3–4× with up to 98% fewer allocations.
- Optimized SVM decoding pipeline (`decode_libsvm` and `libsvm_classifier`): pre-allocated thread-local `LiblinearWorkspace`, direct zero-allocation C training call, double-transpose elimination, and fast SIMD decision function reduce inner-loop allocations by 99.4% (68 MiB → 1.3 MiB) with 100% bit-for-bit identical results.
- Optimized CPU Infomax and Extended Infomax ICA (`infomax.jl`): vectorized contiguous SIMD sigmoid evaluation and pre-allocated kurtosis buffers with BLAS GEMM projection; results are bit-for-bit identical.
- Batched multi-component interpolation for ICA topographies and viewers (`plot_topography` and `plot_ica_component_activation`), speeding up interpolation by ~9× and eliminating duplicate interpolation in interactive viewers.
- Optimized Time-Frequency topography channel power extraction in `plot_topography(::TimeFreqData)`, `plot_topography(::Vector{<:SingleDataFrameTF})`, and `plot_topography_stats(::TFStatsResult)`: replaced repeated DataFrame dictionary hash table re-indexing and full matrix allocations with fast contiguous row-interval filtering and non-finite aware vector accumulation (~150× faster, 99.0% fewer allocations, bit-for-bit identical results).
- Batched multi-panel topography interpolation across ERP, Time-Frequency, and Statistical Result topographies (`plot_topography` for multi-condition vectors and multi-time-bin grids, and `plot_topography_stats`): unified $N$-panel channel values into a single $C \times N$ matrix evaluated in a single interpolation pass, speeding up multi-panel topography interpolation by ~9× with ~88% lower memory allocation while preserving bit-for-bit identical results.
- Consolidated Makie scene objects in 2D and 3D electrode layout and topography rendering (`plot_layout_2d!`, `plot_layout_3d!`, and `_plot_channel_highlights!`): replaced multiple `arc!` and `lines!` calls with a single pre-calculated `lines!` call using `NaN` separators for head outline features (circle, ears, nose), and batched per-channel `text!` loops into a single vector `text!` call, reducing Makie plot object allocations by 60%–96% (e.g. 77 → 3 plot objects per labeled layout) and speeding up layout rendering by ~2.5×.
- Optimized interactive epoched data navigation in `plot_databrowser`: bound channel lines and labels directly to the `current_epoch` observable rather than destroying and recreating all scene plot objects on every arrow keypress or slider step, eliminating scene reconstruction and speeding up epoch navigation by over 1,000× with zero plot object allocations.
- Optimized multi-channel spectral estimation and axis limit calculation in `plot_power_spectrum` (`plot_channel_spectrum` and `plot_ica_component_spectrum`): reused `DSP.WelchConfig` FFT planning, window weights, and shared frequency buffers across all channels (~2.4× faster pure spectral estimation, 86.1% memory reduction), replaced `vcat` splatting and positive-filtering with a single-pass extrema helper (0 allocations), and bypassed temporary DataFrame generation in `plot_ica_component_spectrum`, maintaining bit-for-bit identical results.
- Optimized multi-epoch data extraction and layout rendering in `plot_erp_image`: pre-resolved channel column indices, eliminated per-epoch temporary allocations in single-layout averaging, hoisted invariant time limits and colorbar positioning outside the channel loop, and safely fell back to layout metadata for scale positions in topographic layouts (99.0% reduction in extraction allocations, bit-for-bit identical results).
- Vectorized trigger event rendering in `plot_trigger_overview` and `plot_trigger_timing`: replaced per-trigger `lines!` and `text!` loops with single NaN-separated `lines!` calls and vectorized `text!` arrays, reducing Makie scene plot objects by up to 99.9% (e.g. from 1,000+ objects down to 6–10 objects), speeding up rendering by ~10–24× with ~85% lower memory allocation and identical visual output.
- Optimized Time-Frequency heatmap data extraction and fixed theme colormap resolution in `plot_tf`: precomputed the DataFrame row-coordinate index once per dataset, pre-rounded time/frequency grids, cached extracted matrices to eliminate redundant extraction during colorrange computation, and replaced `filter(!isnan, ...)` array allocations with single-pass non-finite extrema calculation (~24.4× faster extraction, 97.3% memory reduction, 99.3% fewer allocations, bit-for-bit identical results).
- Optimized multi-channel trial concatenation and ERP average overlay rendering in `plot_epochs`: precomputed shared `time_cat` NaN-separated time vectors once per condition across all channels, pre-averaged ERP data once per dataset to eliminate repeated `average_epochs` passes across channels, resolved DataFrame integer column indices upfront to bypass per-trial Symbol dictionary lookups, and fixed `KeyError` crashes on missing grid/topo layout options (~16.6× faster data preparation, 60.4% memory reduction, 96.2% fewer allocations, bit-for-bit identical results).
- Vectorized background channel traces in `plot_gfp`: replaced repeated per-channel `lines!` drawing loops with a single pre-allocated NaN-separated `lines!` call, reducing plot scene objects by 96.9% (e.g. from 65 objects down to 2) and speeding up line rendering by ~40× with 87.1% lower memory allocation and identical visual output.

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
