# Utilities

This section documents utility functions for data manipulation, validation, and helper operations.

## Channel Operations

### `channel_number_to_channel_label(channel_labels, channel_numbers)`
Convert channel numbers to their corresponding labels.

**Arguments:**
- `channel_labels::Vector{Symbol}`: List of all channel labels
- `channel_numbers`: Channel number(s) to convert (can be Int, Vector{Int}, or UnitRange)

**Returns:** Vector{Symbol} - Channel labels corresponding to the input numbers

### `get_channel_indices(dat::DataFrame, channel_labels)`
Get column indices for specified channel labels.

**Arguments:**
- `dat::DataFrame`: Data containing channels
- `channel_labels::AbstractVector{<:AbstractString}`: Channel labels to find

**Returns:** Vector{Int} - Column indices for requested channels

## Data Validation

### `validate_baseline_interval(time, baseline_interval)`
Validate and convert baseline interval to index format.

**Arguments:**
- `time::AbstractVector`: Time points vector
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Interval to validate

**Returns:** IntervalIdx - Validated interval in index format

## Time and Index Operations

### `find_idx_range(time, start_time, end_time)`
Find index range corresponding to time interval.

**Returns:** UnitRange{Int} - Range of indices

### `find_idx_start_end(time, start_time, end_time)`
Find start and end indices corresponding to time interval.

**Returns:** Tuple{Int,Int} - Start and end indices

## Array Utilities

### `consecutive(f::Function, A::AbstractVector; step::Int=1)`
Apply function f to consecutive pairs of elements in vector A.

### `splitgroups(v::AbstractVector)`
Split vector into groups based on consecutive numbers.

### `search_sequence(array::AbstractVector, sequence::Int)`
Find starting indices of a sequence in an array.

## Geometry and Visualization

### `create_convex_hull(xpos, ypos, border_size)`
Create a convex hull around a set of 2D points with a specified border size using Graham's Scan algorithm.

### `best_rect(n)`
Calculate optimal rectangular dimensions for plotting n items.

## Data Processing

### `detrend(x, y)`
Remove linear trend from data.

### `extract_int(s::String)`
Extract integer from string, returning nothing if no digits found.

## Copy Methods

Custom copy functions are provided for all main data types to avoid deepcopy overhead:
- `Base.copy(::ContinuousData)`
- `Base.copy(::EpochData)` 
- `Base.copy(::ErpData)`
- `Base.copy(::AnalysisInfo)` 