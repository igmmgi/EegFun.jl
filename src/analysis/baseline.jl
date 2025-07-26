
using Logging

"""
    _apply_baseline!(dat::DataFrame, channels, baseline_interval)

Internal function that applies baseline correction to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to baseline correct
- `channels::Vector{Symbol}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Subtracts the mean of the baseline interval from each specified channel
- Modifies the input DataFrame in-place
"""
function _apply_baseline!(
    dat::DataFrame,
    channels::Vector{Symbol},
    baseline_interval::Union{IntervalIdx,IntervalTime}
)
    # Compute the mean once for baseline interval
    baseline_means = mean.(eachcol(dat[baseline_interval.interval_start:baseline_interval.interval_end, channels]))

    # Apply baseline correction
    @inbounds for (channel, mean_val) in zip(channels, baseline_means)
        @views dat[!, channel] .-= mean_val
    end
end

"""
    _apply_baseline!(dat::Vector{DataFrame}, channels, baseline_interval)

Internal function that applies baseline correction to each DataFrame in a vector using broadcasting.
"""
function _apply_baseline!(
    dat::Vector{DataFrame},
    channels::Vector{Symbol},
    baseline_interval::Union{IntervalIdx,IntervalTime}
)
    _apply_baseline!.(dat, Ref(channels), Ref(baseline_interval))
end

"""
    baseline!(dat::EegData, baseline_interval; channel_selection=channels())

Apply baseline correction in-place to EEG data.

# Arguments
- `dat::EegData`: The data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Effects
- Modifies the input data in-place by subtracting the baseline mean
- Uses the specified time interval for baseline calculation
"""
function baseline!(
    dat::EegData,
    baseline_interval::Union{IntervalIdx,IntervalTime};
    channel_selection::Function = channels()
)
    # Validate baseline interval
    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection)
    
    if isempty(selected_channels)
        @minimal_warning "No channels selected for baseline correction"
        return
    end
    
    @info "Applying baseline correction to $(length(selected_channels)) channels over interval: $(baseline_interval.interval_start) to $(baseline_interval.interval_end)"
    
    # Apply baseline correction (dispatch handles DataFrame vs Vector{DataFrame})
    _apply_baseline!(dat.data, selected_channels, baseline_interval)
end

"""
    baseline!(dat::EegData; channel_selection=channels())

Apply baseline correction in-place to EEG data using the entire time range.

# Arguments
- `dat::EegData`: The data to baseline correct
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Effects
- Modifies the input data in-place by subtracting the baseline mean
- Uses the entire time range for baseline calculation
"""
function baseline!(
    dat::EegData;
    channel_selection::Function = channels()
)
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, baseline_interval; channel_selection)
end

# generates all non-mutating versions
@add_nonmutating baseline!
