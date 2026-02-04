# EEG Data Structures

Understanding EegFun.jl's core data structures.

## Type Hierarchy

```
EegFunData (abstract)
├── EegData (abstract) - All EEG data types
│   ├── SingleDataFrameEeg (abstract) - Data in a single DataFrame
│   │   ├── ContinuousData - Raw continuous recordings
│   │   ├── ErpData - Averaged event-related potentials
│   │   ├── TimeFreqData - Time-frequency decomposition
│   │   └── SpectrumData - Power spectrum
│   └── MultiDataFrameEeg (abstract) - Data across multiple DataFrames
│       ├── EpochData - Trial-segmented data
│       └── TimeFreqEpochData - TF with individual trials
└── StatsResult (abstract) - Statistical analysis results
```

## ContinuousData

Raw, unsegmented EEG recordings.

| Field | Type | Description |
|-------|------|-------------|
| `file` | `String` | Source filename |
| `data` | `DataFrame` | Continuous data (time × channels) |
| `layout` | `Layout` | Electrode positions and neighbors |
| `sample_rate` | `Int64` | Sampling frequency in Hz |
| `analysis_info` | `AnalysisInfo` | Preprocessing metadata |

## EpochData

Trial-segmented data for event-related analysis.

| Field | Type | Description |
|-------|------|-------------|
| `file` | `String` | Source filename |
| `condition` | `Int64` | Condition number |
| `condition_name` | `String` | Human-readable condition name |
| `data` | `Vector{DataFrame}` | One DataFrame per epoch |
| `layout` | `Layout` | Electrode positions |
| `sample_rate` | `Int64` | Sampling frequency in Hz |
| `analysis_info` | `AnalysisInfo` | Preprocessing metadata |

## ErpData

Averaged event-related potentials.

| Field | Type | Description |
|-------|------|-------------|
| `file` | `String` | Source filename |
| `condition` | `Int64` | Condition number |
| `condition_name` | `String` | Human-readable condition name |
| `data` | `DataFrame` | Averaged ERP (time × channels) |
| `layout` | `Layout` | Electrode positions |
| `sample_rate` | `Int64` | Sampling frequency in Hz |
| `analysis_info` | `AnalysisInfo` | Preprocessing metadata |
| `n_epochs` | `Int64` | Number of epochs averaged |

## TimeFreqData

Time-frequency decomposition results.

| Field | Type | Description |
|-------|------|-------------|
| `file` | `String` | Source filename |
| `condition` | `Int64` | Condition number |
| `condition_name` | `String` | Human-readable condition name |
| `data_power` | `DataFrame` | Power values (time × freq × channels) |
| `data_phase` | `DataFrame` | Phase values in radians |
| `layout` | `Layout` | Electrode positions |
| `sample_rate` | `Int64` | Original sampling frequency |
| `method` | `Symbol` | Analysis method (`:wavelet`, `:multitaper`, etc.) |
| `baseline` | `BaselineInfo` | Baseline correction info (if applied) |
| `analysis_info` | `AnalysisInfo` | Preprocessing metadata |

## SpectrumData

Power spectral density (frequency domain only).

| Field | Type | Description |
|-------|------|-------------|
| `file` | `String` | Source filename |
| `condition` | `Int64` | Condition number |
| `condition_name` | `String` | Human-readable condition name |
| `data` | `DataFrame` | PSD values (freq × channels) |
| `layout` | `Layout` | Electrode positions |
| `sample_rate` | `Int64` | Original sampling frequency |
| `method` | `Symbol` | Analysis method (`:welch`, etc.) |
| `analysis_info` | `AnalysisInfo` | Preprocessing metadata |

## Layout

Electrode spatial information and neighbor relationships.

| Field | Type | Description |
|-------|------|-------------|
| `data` | `DataFrame` | Positions with columns: `label`, `x`, `y`, `z` |
| `neighbours` | `OrderedDict{Symbol,Neighbours}` | Neighbor info per electrode |
| `criterion` | `Float64` | Distance criterion for neighbors (mm) |

## Data Flow

```
Raw BDF/BrainVision File
    ↓ read_bdf() / read_brainvision()
ContinuousData
    ↓ epochs()
EpochData
    ↓ average()
ErpData
    ↓ tf_morlet() / tf_multitaper()
TimeFreqData
```

## See Also

- [Types Reference](../reference/types.md)
- [Getting started tutorial](../tutorials/getting-started.md)
