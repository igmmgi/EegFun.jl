# =============================================================================
# ARTIFACT DETECTION
# =============================================================================

"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)

Detect EOG (eye movement) onsets in the data and mark them as extreme values.

# Arguments
- `dat::ContinuousData`: The continuous EEG data object
- `criterion::Real`: Threshold criterion for EOG detection
- `channel_in::Symbol`: Input channel symbol for EOG detection
- `channel_out::Symbol`: Output channel symbol to mark detected EOG events

# Modifies
- `dat`: Adds extreme value markers to the specified output channel

# Examples
```julia
# Detect EOG onsets using HEOG channel
detect_eog_onsets!(dat, 50.0, :HEOG, :is_eog_onset)
```
"""
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)
    if !(channel_in in propertynames(dat.data))
        @minimal_error("Channel $channel_in not found in data")
    end
    
    # Get the input channel data
    input_data = dat.data[!, channel_in]
    
    # Detect onsets (positive-going threshold crossings)
    onsets = _is_extreme_value(input_data, criterion)
    
    # Add the output channel to the data
    dat.data[!, channel_out] = onsets
    
    return nothing
end

"""
    _is_extreme_value(signal::AbstractVector{Float64}, threshold::Float64)

Detect extreme values in a signal using threshold crossing.

# Arguments
- `signal::AbstractVector{Float64}`: Input signal vector
- `threshold::Float64`: Threshold for extreme value detection

# Returns
- `Vector{Bool}`: Boolean vector indicating extreme values

# Examples
```julia
# Detect values above 50 μV
extreme_mask = _is_extreme_value(signal, 50.0)
```
"""
function _is_extreme_value(signal::AbstractVector{Float64}, threshold::Float64)
    return abs.(signal) .> threshold
end

"""
    _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Float64)

Detect extreme values in a signal and store results in a pre-allocated mask.

# Arguments
- `mask::Vector{Bool}`: Pre-allocated boolean vector to store results
- `signal::AbstractVector{Float64}`: Input signal vector
- `threshold::Float64`: Threshold for extreme value detection

# Modifies
- `mask`: Filled with extreme value detection results

# Examples
```julia
# Pre-allocate mask and detect extreme values
mask = Vector{Bool}(undef, length(signal))
_is_extreme_value!(mask, signal, 50.0)
```
"""
function _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Float64)
    @assert length(mask) == length(signal) "Mask and signal must have the same length"
    
    @inbounds for i in eachindex(signal)
        mask[i] = abs(signal[i]) > threshold
    end
    
    return nothing
end

"""
    is_extreme_value(dat::SingleDataFrameEeg, threshold::Real; channel_selection::Function = channels())

Detect extreme values across selected channels in EEG data.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)

# Returns
- `DataFrame`: DataFrame with extreme value detection results for each channel

# Examples
```julia
# Detect extreme values in all channels
extreme_df = is_extreme_value(dat, 100.0)

# Detect extreme values in specific channels
extreme_df = is_extreme_value(dat, 100.0, channel_selection = channels([:Fp1, :Fp2]))
```
"""
function is_extreme_value(dat::SingleDataFrameEeg, threshold::Real; channel_selection::Function = channels())
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error("No channels selected for extreme value detection")
    end
    
    # Create result DataFrame
    result_df = DataFrame()
    result_df[!, :sample] = 1:size(dat.data, 1)
    
    # Detect extreme values for each channel
    for ch in selected_channels
        channel_data = dat.data[!, ch]
        extreme_mask = _is_extreme_value(channel_data, Float64(threshold))
        result_df[!, ch] = extreme_mask
    end
    
    return result_df
end

"""
    is_extreme_value!(dat::SingleDataFrameEeg, threshold::Real; channel_selection::Function = channels())

Detect extreme values across selected channels and add results as new columns to the data.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)

# Modifies
- `dat`: Adds extreme value detection columns to the data

# Examples
```julia
# Detect extreme values and add to data
is_extreme_value!(dat, 100.0)

# Detect extreme values for specific channels
is_extreme_value!(dat, 100.0, channel_selection = channels([:Fp1, :Fp2]))
```
"""
function is_extreme_value!(dat::SingleDataFrameEeg, threshold::Real; channel_selection::Function = channels())
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error("No channels selected for extreme value detection")
    end
    
    # Detect extreme values for each channel and add to data
    for ch in selected_channels
        channel_data = dat.data[!, ch]
        extreme_mask = _is_extreme_value(channel_data, Float64(threshold))
        
        # Add column with extreme value detection results for this specific channel
        column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
        dat.data[!, column_name] = extreme_mask
    end
    
    return nothing
end

"""
    n_extreme_value(dat::SingleDataFrameEeg, threshold::Real; channel_selection::Function = channels())

Count the number of extreme values across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)

# Returns
- `DataFrame`: DataFrame with extreme value counts for each channel

# Examples
```julia
# Count extreme values in all channels
count_df = n_extreme_value(dat, 100.0)

# Count extreme values in specific channels
count_df = n_extreme_value(dat, 100.0, channel_selection = channels([:Fp1, :Fp2]))
```
"""
function n_extreme_value(dat::SingleDataFrameEeg, threshold::Real; channel_selection::Function = channels())
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error("No channels selected for extreme value counting")
    end
    
    # Count extreme values for each channel
    counts = _n_extreme_value(dat.data, selected_channels, Float64(threshold))
    
    # Create result DataFrame
    result_df = DataFrame(
        channel = selected_channels,
        n_extreme = counts,
        threshold = fill(threshold, length(selected_channels))
    )
    
    return result_df
end

"""
    _n_extreme_value(df::DataFrame, channels::Vector{Symbol}, threshold::Float64)

Internal function to count extreme values for specified channels.

# Arguments
- `df::DataFrame`: DataFrame containing the data
- `channels::Vector{Symbol}`: Vector of channel symbols to analyze
- `threshold::Float64`: Threshold for extreme value detection

# Returns
- `Vector{Int}`: Count of extreme values for each channel
"""
function _n_extreme_value(df::DataFrame, channels::Vector{Symbol}, threshold::Float64)
    counts = Int[]
    
    for ch in channels
        channel_data = df[!, ch]
        extreme_mask = _is_extreme_value(channel_data, threshold)
        push!(counts, sum(extreme_mask))
    end
    
    return counts
end
