# EegFun.jl Architecture

## 1. Core Data Structures (`src/types/`)

All core EEG containers subtype `EegData <: EegFunData` and bundle tabular signal data (`DataFrame`s), electrode `Layout`s, sampling rates, and `AnalysisInfo` metadata:

- **`ContinuousData <: SingleDataFrameEeg`**: Continuous recordings (`data::DataFrame`).
- **`ErpData <: SingleDataFrameEeg`**: Condition-averaged Event-Related Potentials (`data::DataFrame`, `n_epochs::Int`).
- **`EpochData <: MultiDataFrameEeg`**: Single-trial segmented data for one condition (`data::Vector{DataFrame}`).
- **`TimeFreqData <: SingleDataFrameEeg` / `TimeFreqEpochData <: MultiDataFrameEeg`**: Condition-averaged and single-trial time-frequency power (`data_power`, aliased by `.data`) and phase (`data_phase`).
- **`SpectrumData <: SingleDataFrameEeg`**: Frequency-domain power spectral density (`data::DataFrame`).
- **`Layout`**: Polar/2D/3D sensor coordinates and spatial `Neighbors` graphs.

## 2. Repository & Module Hierarchy

`src/EegFun.jl` defines the public API via `src/public_api.jl` and includes subsystem hubs (`src/<subsystem>.jl`), mirrored 1-to-1 by `test/tests/`:

- `src/public_api.jl` — Alphabetized, categorized `public` symbol declarations (Julia 1.11+ `public` API; nothing is `export`ed by default except when explicitly imported).
- `src/types/` (`src/types.jl`) — Core EEG containers (`core.jl`), pipeline configs (`pipeline.jl`), statistics (`statistics.jl`), decoding (`decoding.jl`), and RSA (`rsa.jl`) types.
- `src/utils/` (`src/utils.jl`) — Error handling (`error.jl`), data access/subsetting (`data/`), JLD2 batch processing (`batch.jl`), logging, printing, GPU kernels, and synthetic test-data generators.
- `src/layouts/` — Electrode layout loading, coordinate projections, and spatial neighbor graphs.
- `src/io/` (`src/io.jl`) — Raw EEG readers (`read_data.jl`: BioSemi `.bdf`, `.edf`, BrainVision, `.fif`), EEGLAB, FieldTrip, XDF, and BIDS import/export.
- `src/config/` — Typed `PARAMETERS` registry, TOML validation, and template generation.
- `src/analysis/` — Analytical subsystems:
  - `processing/` (`src/analysis/processing.jl`):
    - `artifact_detection/` — EOG onset, extreme/step/flatline/peak-to-peak detectors, bridged channels, epoch rejection, channel repair, and cross-validated Autoreject.
    - `signal/` — Baseline correction, detrending, IIR/FIR filtering, CleanLine, mirroring, and polyphase resampling.
    - `channels/` — Channel averaging, deletion, bipolar differences, metrics/z-scores, spherical-spline & neighbor repair, summaries, rereferencing, and surface Laplacian (CSD).
    - `epochs/` — Trigger inspection, sequence-matched epoch extraction, and per-epoch transforms.
    - `conditions/` — Condition combination, difference waves, condition averaging, and grand averaging.
    - `ica/` — Infomax, Extended Infomax, Picard, and AMICA decomposition plus automated EOG/ECG/line-noise/channel-noise IC identification.
  - `erps/` (`src/analysis/erps.jl`) — Amplitude/area/fractional-latency measurements, GFP, leave-one-out jackknifing & pseudovalues, LRP, and trial realignment.
  - `ride/` (`src/analysis/ride.jl`) — Residue Iteration Decomposition (RIDE) and Woody filtering.
  - `time_frequency/` (`src/analysis/time_frequency.jl`) — Morlet wavelets, multitaper (DPSS), STFT, ITPC, Welch/FFT spectra, and TF baseline normalization.
  - `statistics/` (`src/analysis/statistics.jl`) — Analytic tests, spatiotemporal cluster-mass & TFCE permutation tests, and mass-univariate linear mixed models (`fit_mass_lmm`).
  - `decoding/` (`src/analysis/decoding.jl`) — Time-resolved and temporal-generalization MVPA (LIBLINEAR/LIBSVM) and cluster statistics.
  - `rsa/` (`src/analysis/rsa.jl`) — Representational dissimilarity matrices, cross-nobis/Mahalanobis distances, noise ceilings, and model comparisons.
- `src/pipelines/` (`src/pipelines.jl`) — Configurable end-to-end preprocessing pipelines, repair tracking, summaries, and script template generation.
- `src/plots/` (`src/plots.jl`) — Makie visualization and interactive GUIs (`utils/`, `processing/`, `erps/`, `time_frequency/`, `decoding/`, `rsa/`, `plot_gui.jl`).
- `src/demos/` (`src/demos.jl`) — Interactive Makie teaching simulations (signal processing, ICA, ERP averaging, MVPA, EOG, EMG, ECG, and EDA).
- `ext/` — Package extensions: `EegFunCUDAExt` (GPU acceleration), `EegFunDimensionalDataExt` (`DimArray`/`DimStack`), and `EegFunMixedModelsExt` (mass-univariate LMMs with Bell–McCaffrey CR2 + Freedman–Lane block sign-flip inference).
- `test/` — Self-contained unit/integration suites under `test/tests/` mirroring `src/` and `ext/`, orchestrated by `test/runtests.jl` and guarded by structural invariant checks in `test/tests/code_organization.jl`.

## 3. Code Organization & API Conventions

All files in `src/`, `ext/`, and `test/` follow strict structural invariants verified automatically in `test/tests/code_organization.jl`:

1. **File Header & Section Order**: Every `.jl` file begins with a 1-line `# ...` purpose comment. Implementation files in `src/` and `ext/` order top-level definitions from top to bottom under 38-character section banners (`Types & Constants` → `User Functions` → `Internal Functions`). Pure `include(...)` hub files contain no banners.
2. **Public vs. Internal Naming**: Every non-underscored function defined in `EegFun` is declared `public` in `src/public_api.jl`; every internal helper function starts with a leading underscore (`_foo`).
3. **Dispatch Ladder & Primary Docstrings**: Within `User Functions`, methods of each function family are ordered `(1) non-mutating in-memory` → `(2) in-place mutating !` → `(3) ::String file/batch wrappers`, followed by `Base.`/`Makie.`/`Tables.` interface extensions. Full docstrings live on the primary non-mutating in-memory method (`foo`), while `foo!` cross-references `foo`.
4. **Canonical Docstrings & Formatting**: Every docstring starts with a 4-space indented signature/binding header and uses canonical `#` Markdown headings in fixed order (`# Fields` / `# Arguments` → `# Keyword Arguments` → `# Returns` → `# Details` / `# Notes` → `# Examples` → `# References` → `# See Also`), with per-field `# Fields` bullets on all concrete `public` structs. All `.jl` files conform to `.JuliaFormatter.toml`.
5. **Predicate Selection & Fail-Fast Errors**:
   - Subsetting uses composable selector predicates (`channels(:Fz)`, `times(-0.2, 0.8)`, `epochs(1:10)`, `conditions(1:2)`, `participants(1:20)`).
   - User-fixable input/config/file/selection errors raise `EegFunError` via `@minimal_error` (logging to file only when batch/pipeline file logging is active). Plain `error()` is reserved for internal invariant bugs.

## 4. Pipeline Data Flow (`src/pipelines/`)

The `preprocess(config_path)` batch pipeline executes the following sequence:

1. Parse and validate `pipeline.toml` against the `PARAMETERS` registry into a `PipelineConfig`.
2. Load raw continuous EEG data (`ContinuousData`) and electrode `Layout` / spatial neighbors.
3. Apply high-pass/low-pass filters, optional CleanLine line-noise removal, and resampling.
4. Derive bipolar EOG channels, flag continuous artifacts, and interpolate bad channels.
5. Train ICA on clean data segments, identify ocular/cardiac/line/channel-noise ICs, and subtract artifactual components.
6. Extract condition-specific `EpochData` from trigger sequences and apply baseline correction.
7. Flag/repair/reject bad epochs (via threshold/z-score criteria, Autoreject, or interactive review).
8. Save cleaned `EpochData`, condition-averaged `ErpData`, `InfoIca`, and `ArtifactInfo` archives (`.jld2`) plus structured summary logs.
