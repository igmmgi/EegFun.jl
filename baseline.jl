using Logging

"""
    _apply_baseline!(dat::DataFrame, channel_labels, baseline_interval)

Internal function that applies baseline correction to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Subtracts the mean of the baseline interval from each specified channel
- Modifies the input DataFrame in-place

# Throws
- `ArgumentError`: If channel_labels is empty or channels not found in data
- `ArgumentError`: If baseline_interval is invalid
"""
function _apply_baseline!(
    dat::DataFrame,
    channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx,IntervalTime},
)

    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    channel_indices = get_channel_indices(dat, channel_labels)

    # Compute the mean once for baseline interval
    baseline_means =
        mean.(eachcol(dat[baseline_interval.interval_start:baseline_interval.interval_end, channel_indices]))

    # Compute baseline means
    baseline_data = dat[baseline_interval.interval_start:baseline_interval.interval_end, channel_labels]
    baseline_means = mean.(eachcol(baseline_data))

    @info "Applying baseline correction to channels: $(join(channel_labels[1:min(3, end)], ", ")) over IntervalIdx: $(baseline_interval.interval_start) to $(baseline_interval.interval_end)"

    # Apply baseline correction
    for (col, mean_val) in zip(channel_labels, baseline_means)
        @inbounds dat[!, col] .-= mean_val
    end

end

"""
    baseline!(dat::EEGData, channel_labels, baseline_interval)

Apply baseline correction in-place to EEG data. Works with different EEG data types
(EpochData, DataFrame, etc.) through multiple dispatch.

# Arguments
- `dat::EEGData`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Modifies the input data in-place by subtracting the baseline mean
"""
function baseline!(
    dat::EegData,
    channel_labels::Union{Vector{Symbol}, Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx, IntervalTime}
)
    if dat isa EpochData
        apply_baseline_to_epochs!(dat, channel_labels, baseline_interval)
    elseif dat isa DataFrame
        _apply_baseline!(dat, channel_labels, baseline_interval)
    else
        _apply_baseline!(dat.data, channel_labels, baseline_interval)
    end
end

"""
    apply_baseline_to_epochs!(dat::EpochData, channel_labels, baseline_interval)

Apply baseline correction to each epoch in epoched EEG data.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Modifies each epoch in-place by subtracting its baseline mean
"""
function apply_baseline_to_epochs!(
    dat::EpochData,
    channel_labels::Union{Vector{Symbol}, Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx, IntervalTime}
)
    for epoch in eachindex(dat.data)
        _apply_baseline!(dat.data[epoch], channel_labels, baseline_interval)
    end
end

"""
    baseline(dat::EEGData, channel_labels, baseline_interval)

Create a baseline-corrected copy of the EEG data.

# Arguments
- `dat::EEGData`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Returns
- `EEGData`: A new copy of the input data with baseline correction applied
"""
function baseline(
    dat::EegData,
    channel_labels::Union{Vector{Symbol}, Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx, IntervalTime}
)
    dat_out = deepcopy(dat)
    baseline!(dat_out, channel_labels, baseline_interval)
    return dat_out
end

"""
    baseline(dat::EEGData, baseline_interval)

Create a baseline-corrected copy using all channels from the data's layout.

# Arguments
- `dat::EEGData`: The data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Returns
- `EEGData`: A new copy of the input data with baseline correction applied to all channels
"""
function baseline(
    dat::EegData,
    baseline_interval::Union{IntervalIdx, IntervalTime}
)
    dat_out = deepcopy(dat)
    baseline(dat, dat.layout.label, baseline_interval)
    return dat_out
end

"""
    baseline!(dat::EEGData, baseline_interval)

Apply baseline correction in-place using all channels from the data's layout.

# Arguments
- `dat::EEGData`: The data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Modifies the input data in-place by subtracting the baseline mean from all channels
"""
function baseline!(
    dat::EegData,
    baseline_interval::Union{IntervalIdx, IntervalTime}
)
    baseline!(dat, dat.layout.label, baseline_interval)
end

