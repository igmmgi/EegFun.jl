# EegFun.jl Architecture

## 1. Core Data Structures (`src/types/`)

All core structures inherit from `EegData` and wrap signal data (via a `DataFrame` in the `.data` field), spatial layouts, and metadata.
- **`ContinuousData`**: Continuous time-series data.
- **`EpochData`**: Trial-segmented data grouped by condition. 
- **`ErpData`**: Averaged trials (Event-Related Potentials).
- **`TimeFreqData` / `SpectrumData`**: Frequency-domain data.
- **`Layout`**: Sensor array coordinates and meshes (2D/3D).

## 2. Directory Structure (`src/`)A

- `types/`: Struct definitions.
- `io/`: Parsers (EEGLAB, FieldTrip, BDF, BIDS).
- `config/`: The `PARAMETERS` registry and TOML generation.
- `layouts/`: Coordinate projection and neighbor mesh calculations.
- `utils/`: Logging, error handling, and indexing tools.
- `analysis/`: Analytical operations.
  - `processing/`: Artifact rejection, ICA, CSD, filtering, channel repair.
  - `erps/`: Grand averaging, peak detection, jackknifing.
  - `time_frequency/`: Wavelet and multitaper decompositions.
  - `statistics/`: Permutation testing and thresholding.
  - `decoding/` & `rsa/`: MVPA implementations.
- `pipelines/`: Batch processing implementation (`preprocess()`).
- `plots/`: Makie visualization tools and GUIs.

## 3. Core Design Patterns

- **Mutating Defaults (`!`):** Core algorithms are implemented as in-place mutating functions (e.g., `filter!(dat, ...)`). Non-mutating wrappers (`filter(...)`) are automatically generated via the internal `@add_nonmutating` macro.
- **Predicate Selection:** Data subsetting uses predicate functions instead of array indices (e.g., `channels(:Fz)`, `epochs(1:5)`).
- **Fail-Fast Error Handling:** Data integrity and types are validated immediately at function entry, throwing domain-specific errors before propagating to lower-level libraries.

## 4. Pipeline Data Flow (`src/pipelines/`)

The `preprocess(config_path)` batch pipeline executes the following sequence:

1. Parse and validate `pipeline.toml` against the `PARAMETERS` registry.
2. Load raw data (`ContinuousData`).
3. Apply filters, CleanLine, and downsampling.
4. Detect EOG events, threshold continuous artifacts, and interpolate bad channels.
5. Compute and subtract Independent Components (ICA).
6. Segment continuous data into epochs via trigger events.
7. Reject epochs exceeding amplitude/variance thresholds.
8. Output cleaned `EpochData` and `ErpData` as `.jld2` archives.
