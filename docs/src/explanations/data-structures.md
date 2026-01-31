# EEG Data Structures

Understanding EegFun.jl's core data structures.

## Overview

EegFun.jl uses several key types to organize EEG data:

- `EegDataFrame` - Continuous EEG data
- `EpochData` - Time-locked trial data
- `ErpData` - Averaged event-related potentials
- `Layout` - Electrode spatial information

## EegDataFrame

The primary structure for continuous EEG data:

```julia
struct EegDataFrame
    data::DataFrame           # Time × channels matrix with metadata
    layout::Layout           # Electrode positions
    sampling_rate::Float64   # Sampling frequency in Hz
    # ... other fields
end
```

### Key Features

- **DataFrame backbone**: Time points as rows, channels as columns
- **Metadata columns**: `:time`, `:sample`, trigger information
- **Artifact tracking**: Boolean columns for quality control

> **TODO**: Add detailed field descriptions

## EpochData vs ErpData

### EpochData

Individual trials extracted from continuous data:

```julia
# TODO: Add EpochData structure definition
```

**Use when**: You need trial-level data for statistical analysis or quality control

### ErpData

Averaged responses across trials:

```julia
# TODO: Add ErpData structure definition
```

**Use when**: Computing grand averages or group-level ERPs

## Layout

Electrode spatial information:

```julia
struct Layout
    labels::Vector{String}    # Electrode names
    x::Vector{Float64}        # X coordinates
    y::Vector{Float64}        # Y coordinates
    # ... other fields
end
```

> **TODO**: Add coordinate system explanation

## Metadata Handling

EegFun.jl uses a structured approach to metadata:

### Time-Only Preservation

> **TODO**: Explain the "time-only metadata preservation" rule used in averaging

### Result Traceability

> **TODO**: Explain the `:file`-first schema for analysis results

## Data Flow

```
Raw BDF File
    ↓ read_raw_data()
EegDataFrame (continuous)
    ↓ extract_epochs()
EpochData (trials)
    ↓ average_epochs()
ErpData (averaged)
```

## See Also

- [API design patterns](../reference/patterns.md)
- [Getting started tutorial](../tutorials/getting-started.md)
