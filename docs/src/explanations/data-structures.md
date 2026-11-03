# EEG Data Structures

Understanding EegFun.jl's core data structures.

## Julia Types

If you're coming from Python or MATLAB, Julia's type system works a little differently. The full details are in the [Julia documentation](https://docs.julialang.org/en/v1/manual/types/), but the key idea for everyday EegFun.jl use is **multiple dispatch**.

This means a function can have the same name but do different things depending on the *type* of data passed to it.

In Julia function signatures, the `::` notation means *"this argument must be of this type"*:

```julia
EegFun.plot_erp(erp::EegFun.ErpData)         # plots a single ERP
EegFun.plot_erp(erp::Vector{EegFun.ErpData}) # plots multiple ERPs overlaid
```

Here `erp::EegFun.ErpData` means "an argument called `erp` of type `EegFun.ErpData`". Julia automatically selects the right function version based on what you pass in — you never need to call a different function name.

More examples from EegFun.jl:

```julia
# Same function, different input types → different behaviour
EegFun.plot_topography(erp::EegFun.ErpData, ...)          # plots ERP topography
EegFun.plot_topography(tf::EegFun.TimeFreqData, ...)      # plots time-frequency topography

EegFun.average_epochs(dat::EegFun.EpochData)              # averages a single condition
EegFun.average_epochs(dat::Vector{EegFun.EpochData})      # averages multiple conditions

EegFun.extract_epochs(dat::EegFun.ContinuousData, ...)    # extracts from continuous recording
```

Both are called `plot_erp` — Julia picks the right function version automatically based on what you pass in. You never need to call a different function name.

If you're ever unsure what type a variable is, use `typeof()`:

```julia
typeof(my_data)   # e.g. EegFun.ContinuousData
```

This is useful when reading error messages — if Julia says *"no method matching `plot_erp(::ContinuousData)`"* it means you passed the wrong type, and `typeof()` can help you check what you actually have.

### Abstract vs Concrete Types

Julia also distinguishes between **concrete** and **abstract** types:

- **Concrete types** are the actual data structures you create and work with — for example `ContinuousData`, `ErpData`, or `EpochData`. These hold real data.
- **Abstract types** are categories or groupings — you can never create one directly, but they let functions accept a *family* of related types. For example, `EegData` is abstract; it simply says "any EEG data type".

This matters in practice because a function written for `EegData` will work on *any* concrete EEG type — `ContinuousData`, `ErpData`, etc. — without needing a separate version for each one.

---

## Type Hierarchy


```
EegFunData (abstract)
└── EegData (abstract) - All EEG data types
    ├── SingleDataFrameEeg (abstract) - Data in a single DataFrame
    │   ├── ContinuousData - Raw continuous recordings
    │   ├── ErpData - Averaged event-related potentials
    │   ├── TimeFreqData - Time-frequency decomposition
    │   └── SpectrumData - Power spectrum
    └── MultiDataFrameEeg (abstract) - Data across multiple DataFrames
        ├── EpochData - Trial-segmented data
        └── TimeFreqEpochData - TF with individual trials
```

See the [Types Reference](../reference/types.md) for full field-by-field documentation of each type.

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
| `neighbors` | `Union{Nothing, OrderedDict{Symbol,Neighbors}}` | Neighbor info per electrode |
| `criterion` | `Union{Nothing, Float64}` | Distance criterion used to compute neighbors (mm) |
| `criterion_type` | `Union{Nothing, Symbol}` | Type of criterion (e.g. `:xy`, `:xyz`) |

## Data Flow

```
Raw EEG File (.bdf, .vhdr, .set, .mat)
    ↓ read_raw_data()
ContinuousData
    ↓ extract_epochs()
EpochData
    ↓ average_epochs()
ErpData
```

## See Also

- [Types Reference](../reference/types.md)
- [Getting started tutorial](../tutorials/getting-started.md)
