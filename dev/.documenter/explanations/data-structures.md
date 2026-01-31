
# EEG Data Structures {#EEG-Data-Structures}

Understanding EegFun.jl&#39;s core data structures.

## Overview {#Overview}

EegFun.jl uses several key types to organize EEG data:
- `EegDataFrame` - Continuous EEG data
  
- `EpochData` - Time-locked trial data
  
- `ErpData` - Averaged event-related potentials
  
- `Layout` - Electrode spatial information
  

## EegDataFrame {#EegDataFrame}

The primary structure for continuous EEG data:

```julia
struct EegDataFrame
    data::DataFrame           # Time × channels matrix with metadata
    layout::Layout           # Electrode positions
    sampling_rate::Float64   # Sampling frequency in Hz
    # ... other fields
end
```


### Key Features {#Key-Features}
- **DataFrame backbone**: Time points as rows, channels as columns
  
- **Metadata columns**: `:time`, `:sample`, trigger information
  
- **Artifact tracking**: Boolean columns for quality control
  
> 
> **TODO**: Add detailed field descriptions
> 


## EpochData vs ErpData {#EpochData-vs-ErpData}

### EpochData {#EpochData}

Individual trials extracted from continuous data:

```julia
# TODO: Add EpochData structure definition
```


**Use when**: You need trial-level data for statistical analysis or quality control

### ErpData {#ErpData}

Averaged responses across trials:

```julia
# TODO: Add ErpData structure definition
```


**Use when**: Computing grand averages or group-level ERPs

## Layout {#Layout}

Electrode spatial information:

```julia
struct Layout
    labels::Vector{String}    # Electrode names
    x::Vector{Float64}        # X coordinates
    y::Vector{Float64}        # Y coordinates
    # ... other fields
end
```

> 
> **TODO**: Add coordinate system explanation
> 


## Metadata Handling {#Metadata-Handling}

EegFun.jl uses a structured approach to metadata:

### Time-Only Preservation {#Time-Only-Preservation}
> 
> **TODO**: Explain the &quot;time-only metadata preservation&quot; rule used in averaging
> 


### Result Traceability {#Result-Traceability}
> 
> **TODO**: Explain the `:file`-first schema for analysis results
> 


## Data Flow {#Data-Flow}

```
Raw BDF File
    ↓ read_raw_data()
EegDataFrame (continuous)
    ↓ extract_epochs()
EpochData (trials)
    ↓ average_epochs()
ErpData (averaged)
```


## See Also {#See-Also}
- [API design patterns](../reference/patterns.md)
  
- [Getting started tutorial](../tutorials/getting-started.md)
  
