# Data Types

This section documents the main data structures used in eegfun.

## Core Data Types

### ContinuousData
Raw continuous EEG data structure containing:
- `data::DataFrame`: Raw EEG data with time and channel columns
- `layout::DataFrame`: Electrode layout information
- `sample_rate::Real`: Sampling rate in Hz
- `analysis_info::AnalysisInfo`: Analysis metadata

### EpochData
Segmented EEG data around events:
- `data::Vector{DataFrame}`: Vector of epoch DataFrames
- `layout::DataFrame`: Electrode layout information
- `sample_rate::Real`: Sampling rate in Hz
- `analysis_info::AnalysisInfo`: Analysis metadata

### ErpData
Event-related potential data:
- `data::DataFrame`: Averaged ERP data
- `layout::DataFrame`: Electrode layout information
- `sample_rate::Real`: Sampling rate in Hz
- `analysis_info::AnalysisInfo`: Analysis metadata
- `n_epochs::Int`: Number of epochs averaged

### AnalysisInfo
Analysis metadata structure:
- `reference::Symbol`: Reference type used
- `hp_filter::Float64`: High-pass filter cutoff
- `lp_filter::Float64`: Low-pass filter cutoff

## ICA Types

### InfoIca
Independent Component Analysis results structure containing ICA weights and components.

### IcaPrms
ICA parameters structure for configuring ICA analysis.

## Interval Types

### IntervalIdx
Time interval specified by indices:
- `interval_start::Int`: Start index
- `interval_end::Int`: End index

### IntervalTime
Time interval specified by time values:
- `interval_start::Real`: Start time
- `interval_end::Real`: End time 