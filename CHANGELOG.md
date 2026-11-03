# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to (well, attempts to :-)) [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Deep `Base.copy` methods across all core, rejection, statistical, decoding, and RSA data containers (`ContinuousData`, `EpochData`, `ErpData`, `TimeFreqData`, `TimeFreqEpochData`, `SpectrumData`, `EpochRejectionInfo`, `PermutationResult`, `DecodedData`, `RsaData`, etc.).
- Complete `Base` array and `Tables.jl` interfaces, canonical accessors (`time_vector`, `frequency_vector`, etc.), and `show` methods across data and statistics types.
- `format_version = EEGFUN_FORMAT_VERSION` (version `1`) metadata in `.jld2` outputs with load-time compatibility checks in `read_data`, plus `JLD2.rconvert` shims for legacy `Neighbours` and 5-field `Layout` archives and `"neighbour_criterion"` `.toml` configs.
- `list_bids_runs` and BIDS `electrodes.tsv` / `coordsystem.json` layout import in `read_bids`.
- `DimensionalData.DimArray` support for `SpectrumData` and non-contiguous/subsetted `TimeFreqData` and `TimeFreqEpochData` (`EegFunDimensionalDataExt`).
- `tf_baseline` / `baseline!` support for `TimeFreqEpochData`.
- Smulders (2010) jackknife pseudovalue support via `pseudovalues::Bool` in `jackknife_average` and `jackknife_average(::ErpMeasurementsResult)`.
- `area_boundary` (`:rectified`, `:positive`, `:negative`, `:signed`) and `peak_polarity` (`:abs`, `:positive`, `:negative`) options in `erp_measurements`.
- `RideComponent` `time_window` support and reaction-time median centering (`RideResult.latency_offsets`) in `compute_ride`.
- Separate `alpha` (cluster-level FWER) and `cluster_alpha` (point-wise cluster-forming threshold) keywords in `permutation_test`, and explicit `rng` keywords across `permutation_test`, `decode_libsvm`, `test_against_chance_cluster`, `rsa_crossvalidated`, and `compare_models!`.
- Non-mutating `add_noise_ceiling`, `calculate_eog_channels`, and `detect_eog_onsets`, mutating `compare_models!`, and `create_rdm`.

### Changed

- Standardized British `neighbour*` identifiers, types/fields (`Neighbors`, `Layout.neighbors`), keyword arguments (`neighbors`, `neighbor_criterion`), and `.toml` keys to US English `neighbor*` repo-wide.
- Enforced a strict zero-export API (`0` `export` declarations; 400 curated `public` declarations in `src/public_api.jl`), updating `generate_pipeline_template()`, docstrings, and tutorials to qualify calls with `EegFun.`.
- Standardized `plot_*` return values to `NamedTuple`s starting with `(; fig, axes, ...)`, and renamed the TF significance mask keyword in `plot_tf_stats` from `significance` to `mask_style`.
- Aligned `tf_morlet`, `tf_stft`, and `tf_multitaper` with MNE discrete L2 wavelet normalization (`||W||_2 = sqrt(2)`, `power = |c|^2 / fs`) and zero-padded linear FFT convolution.
- Defaulted `analytic_test` multiple-comparison `correction_method` to `:none` (accepting `:no` as an alias).
- Bumped minimum supported Julia version to `1.13` across `Project.toml`, CI workflows, and documentation, updated `AnovaFun` to `0.5`, and tightened `Project.toml` `[compat]` bounds (`EuropeanDataFormat` `0.1.1`, `CUDA` `5, 6`, `Tables` `1`, `FFTA` `0.3.0 - 0.3.1, 0.3.3 - 0.3` excluding broken upstream `0.3.2`).
- Modularized repository source structure (`ica/`, `utils/data/`, `artifact_detection/`) and unified N-dimensional ERP/TF statistical routines and position-based selector negation (`_negate_selector`).
- `using EegFun` no longer overrides the global Makie theme; `theme_eegfun()` is applied per figure so user `set_theme!` / `with_theme` settings take precedence.
- Seeded synthetic-data helpers, interactive demos, and `src/precompile.jl` use local RNGs instead of mutating `Random.default_rng()`.
- Adopted EEGLAB `pop_eegfiltnew` default FIR transition-width heuristics (`min(max(0.25*fc, 2 Hz), fc)` for highpass; `min(max(0.25*fc, 2 Hz), fs/2 - fc)` for lowpass).

### Fixed

- **MVPA Decoding (`decode_libsvm`)**: Added strictly within-fold training-set feature z-scoring (eliminating test-set leakage); enabled linear SVM bias/intercept (`bias = 1.0`); scored imbalanced folds (`equalize_trials = false`) with macro-averaged balanced accuracy so true chance equals `1 / n_classes`; moved `equalize_trials = true` subsampling inside each Monte Carlo iteration; removed erroneous `1 / sqrt(n_iterations)` scaling of iteration SD; and freed C-heap `LIBLINEAR` models via `free_and_destroy_model` under `GC.@preserve`.
- **RSA (`rsa_crossvalidated`, `compute_noise_ceiling`, `rsa`)**: Replaced biased within-half RDM averaging with unbiased cross-fold inner products (Crossnobis / cross-validated squared Euclidean; Walther et al., 2016) and demeaned Ledoit–Wolf covariance shrinkage; added Nili et al. (2014) rank/z-normalization to group noise ceilings; implemented tie-aware Spearman `tiedrank`; and switched model comparison permutation tests to Mantel condition-label permutations with `(extreme_count + 1) / (n_permutations + 1)` p-values.
- **Time-Frequency & Spectrum**: Eliminated circular FFT wrap-around bleed across epoch boundaries in `tf_morlet`, `tf_stft`, and `tf_multitaper`; extracted primary-taper (`k = 1`) phase in `tf_multitaper` (avoiding cancellation across antisymmetric odd DPSS tapers); replaced arithmetic phase averaging with circular means (`atan(mean(sin), mean(cos))`) in `condition_average` and `channel_average`; and changed `freq_spectrum(::EpochData)` from concatenated-epoch Welch to per-epoch Welch PSD averaging.
- **Statistics (`permutation_test`, `analytic_test`, `test_against_chance_cluster`)**: Replaced single-pass `sum(x^2) - sum(x)^2 / n` variance with centered two-pass variance on CPU and GPU; doubled stored `cluster.p_value` (`min(1.0, 2 * p_one_sided)`) when `tail == :both` to match two-sided `alpha`; stopped flagging `NaN` p-values with `t > 0` as significant in `_apply_bonferroni_correction`; evaluated one-tailed non-parametric thresholds across the full permutation `t`-distribution; made permutation RNG pre-drawing and chunk-owned scratch buffers thread-count-independent and safe under `:dynamic` task migration; and avoided mutating input `Layout.neighbors`.
- **ICA (`infomax.jl`, `picard.jl`, `amica.jl`, `artifact_id.jl`, `components.jl`)**: Aligned Extended Infomax natural gradient with Lee et al. (1999) / EEGLAB `runica.m` for sub-Gaussian sources; fixed Extended Picard line-search loss and block Hessian conditioning; added parameter rollback and kurtosis sub-Gaussian updates to AMICA; preserved variance orthogonal to rank-deficient ICA subspaces in `add_ica_components!`; added `AbstractVector{<:EpochData}` support to `subtract_ica_components!` / `subtract_ica_components`; and added MAD-scaled z-scoring, polarity-invariant ECG IBI detection, and Welch line-noise estimation in artifact identification.
- **Signal Processing, Artifact Repair, ERPs & RIDE**: Fixed Thomson F-test harmonic phase reconstruction in `cleanline!`; stabilized Perrin spherical-spline matrix conditioning and metadata-column isolation in `repair_channels!`, `repair_artifacts!`, and `compute_csd!`, with warnings on operator failure in `autoreject!`; fixed `jackknife_average` stamping `first_erp.file` onto all excluded participants and misaligning batch outputs when conditions were missing; added `MultiDataFrameEeg` overloads for `subset_bad_data` and out-of-place `channel_repairable`; and updated `AnalysisInfo.sample_rate` to `Float64` with synchronization during `resample!`.
- **Data Containers & I/O**: Eliminated shallow-copy aliasing of `DataFrame`, `Layout`, and `AnalysisInfo` fields in non-mutating functions; fixed `isequal`/`hash` (`NaN`-safe equality and value-based `Neighbors` hashing) and sample-rate checks in `isapprox`; and fixed non-integer `srate` rounding, missing epoch triggers, `:sample`/`:epoch` columns, and empty layout coordinate initialization in `read_eeglab`, `read_fieldtrip`, and `read_bids`.
- **Plotting (`plot_databrowser`, `plot_tf`, GUIs)**: Fixed `plot_databrowser` ICA activation caching, multi-channel `Ctrl`-click selection, legend channel linking, and "Data - ICA" trace toggling/coloring; respected `display_plot = false` across plotting routines; and guarded GUI observer lifecycles and error logging.

### Removed

- Deprecated mutating batch string aliases (`realign!`, `calculate_trigger_interval!`).
- Internal helper symbols and `StatsAPI` re-declarations from `src/public_api.jl` (`Component`, `coef`, `confint`, `nobs`, `stderror`, `signal_to_data`, `from_dict`, `allcomb`, `consecutive`, `repair_channels_neighbor!`, `repair_artifacts_neighbor!`, `repair_artifacts_spherical_spline!`).
- Inconsistent `Base` methods on EEG data types: `IndexStyle`, scalar `eltype` (`Float64`), linear `firstindex(dat)`/`lastindex(dat)`, and implicit `convert(Matrix/Array/DataFrame, dat)` (use explicit `Matrix(dat)`, `Array(dat)`, or `DataFrame(dat)`; epoch `eltype` is now `DataFrame`).
- Unused internal helpers (`_empty_spatial_kurtosis_metrics`, `_rereference!`, `_get_time_bounds`, `_notify_data_update`, `_reset_to_original!`, `_add_to_collection!`).

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
