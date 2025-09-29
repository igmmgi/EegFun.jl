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
        @minimal_error_throw("Channel $channel_in not found in data")
    end
    
    # Get the input channel data
    input_data = dat.data[!, channel_in]
    
    # Detect onsets (positive-going threshold crossings)
    onsets = _is_extreme_value(input_data, Float64(criterion))
    
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
    is_extreme_value!(dat::SingleDataFrameEeg, threshold::Int; 
                     channel_selection::Function = channels(), 
                     sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                     mode::Symbol = :combined,
                     channel_out::Union{Symbol, Nothing} = nothing)

Detect extreme values across selected channels and add results to the data.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Union{Vector{Bool}, Nothing}`: Boolean vector for sample selection (default: nothing - all samples)
- `mode::Symbol`: Mode of operation - `:combined` (single combined column, default) or `:separate` (separate columns per channel)
- `channel_out::Union{Symbol, Nothing}`: Output channel name for combined mode (default: auto-generated as `is_extreme_value_<threshold>`)

# Modifies
- `dat`: Adds extreme value detection columns to the data

# Examples
```julia
# Detect extreme values and combine into single output channel (default, auto-named)
is_extreme_value!(dat, 100)  # Creates column :is_extreme_value_100

# Detect extreme values with custom output channel name
is_extreme_value!(dat, 100, channel_out = :is_extreme_value)

# Detect extreme values and add separate columns for each channel
is_extreme_value!(dat, 100, mode = :separate)

# Detect extreme values for specific channels (auto-named)
is_extreme_value!(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Detect extreme values only for selected samples (auto-named)
is_extreme_value!(dat, 100, sample_selection = sample_mask)
```
"""
function is_extreme_value!(dat::SingleDataFrameEeg, threshold::Int; 
                          channel_selection::Function = channels(), 
                          sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                          mode::Symbol = :combined,
                          channel_out::Union{Symbol, Nothing} = nothing)
    
    # Validate mode
    if mode ∉ [:separate, :combined]
        @minimal_error_throw("mode must be :separate or :combined, got :$mode")
    end
    
    if mode == :combined
        # Combined mode - use all channels (same as original behavior)
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)
        
        # Use provided channel_out or generate default name
        output_channel = channel_out === nothing ? Symbol("is_extreme_value_$(threshold)") : channel_out
        
        # Initialize the output channel with false values
        dat.data[!, output_channel] = falses(nrow(dat.data))
        
        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            dat.data[!, output_channel] .|= extreme_mask
        end
        
    else  # mode == :separate
        # Separate mode - use specified channel selection
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)
        
        # Add results to data object
        for (ch, extreme_mask) in results
            column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
            dat.data[!, column_name] = extreme_mask
        end
    end
    
    return nothing
end

# Helper function to detect extreme values for selected channels
function _detect_extreme_values(dat::SingleDataFrameEeg, threshold::Real; 
                               channel_selection::Function = channels(), 
                               sample_selection::Union{Vector{Bool}, Nothing} = nothing)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for extreme value detection")
    end
    
    results = Dict{Symbol, Vector{Bool}}()
    
    for ch in selected_channels
        channel_data = dat.data[!, ch]
        extreme_mask = _is_extreme_value(channel_data, Float64(threshold))
        
        # Apply sample selection if provided
        if sample_selection !== nothing
            extreme_mask = extreme_mask .& sample_selection
        end
        
        results[ch] = extreme_mask
    end
    
    return results
end

"""
    is_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                    channel_selection::Function = channels(), 
                    sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                    mode::Symbol = :combined)

Detect extreme values across selected channels and return results.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Union{Vector{Bool}, Nothing}`: Boolean vector for sample selection (default: nothing - all samples)
- `mode::Symbol`: Mode of operation - `:combined` (boolean vector, default) or `:separate` (DataFrame with separate columns per channel)

# Returns
- `DataFrame` or `Vector{Bool}`: Results depending on mode

# Examples
```julia
# Detect extreme values and return combined boolean vector (default)
extreme_mask = is_extreme_value(dat, 100)

# Detect extreme values and return DataFrame with separate columns
results = is_extreme_value(dat, 100, mode = :separate)

# Detect extreme values for specific channels
extreme_mask = is_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Detect extreme values only for selected samples
extreme_mask = is_extreme_value(dat, 100, sample_selection = sample_mask)
```
"""
function is_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                         channel_selection::Function = channels(), 
                         sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                         mode::Symbol = :combined)
    
    if mode == :combined
        # Combined mode - return boolean vector directly
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)
        
        # Initialize the result with false values
        combined_mask = Vector{Bool}(falses(nrow(dat.data)))
        
        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            combined_mask .|= extreme_mask
        end
        
        return combined_mask
        
    else  # mode == :separate
        # Separate mode - create temporary data object and use mutating version
        temp_dat = deepcopy(dat)
        is_extreme_value!(temp_dat, threshold; channel_selection, sample_selection, mode = :separate)
        
        # Extract the extreme value columns in the same order as the original channels
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        extreme_cols = [Symbol("is_extreme_value_$(ch)_$(threshold)") for ch in selected_channels]
        return temp_dat.data[:, extreme_cols]
    end
end

"""
    n_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                   channel_selection::Function = channels(), 
                   sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                   mode::Symbol = :combined)

Count the number of extreme values across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Union{Vector{Bool}, Nothing}`: Boolean vector for sample selection (default: nothing - all samples)
- `mode::Symbol`: Mode for extreme value detection (:separate or :combined, default: :combined)

# Returns
- `DataFrame`: DataFrame with extreme value counts for each channel (separate mode) or total count (combined mode)

# Examples
```julia
# Count total extreme values across all channels (combined mode, default)
total_count = n_extreme_value(dat, 100)

# Count extreme values in specific channels (combined mode)
total_count = n_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Count extreme values only for selected samples (combined mode)
total_count = n_extreme_value(dat, 100, sample_selection = sample_mask)

# Count extreme values in all channels (separate mode)
count_df = n_extreme_value(dat, 100, mode = :separate)
```
"""
function n_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                        channel_selection::Function = channels(), 
                        sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                        mode::Symbol = :combined)
    # Validate mode
    if mode ∉ [:separate, :combined]
        @minimal_error_throw("mode must be :separate or :combined, got :$mode")
    end

    if mode == :combined
        # Combined mode - use all channels (same as original behavior)
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)
        
        # Initialize the result with false values
        combined_mask = Vector{Bool}(falses(nrow(dat.data)))
        
        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            combined_mask .|= extreme_mask
        end
        
        # Return total count
        return sum(combined_mask)
        
    else  # mode == :separate
        # Separate mode - count for each channel
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)
        
        # Get channels in original order
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        
        # Count extreme values for each channel in original order
        counts = [sum(results[ch]) for ch in selected_channels]
        
        # Create result DataFrame
        result_df = DataFrame(
            channel = selected_channels,
            n_extreme = counts,
            threshold = fill(threshold, length(selected_channels))
        )
        
        return result_df
    end
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
