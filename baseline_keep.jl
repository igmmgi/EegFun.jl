using Logging

"""
    _apply_baseline!(dat::DataFrame, channel_indices::Vector{Int}, baseline_interval)

Internal function that applies baseline correction to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to baseline correct
- `channel_indices::Vector{Int}`: Indices of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Subtracts the mean of the baseline interval from each specified channel
- Modifies the input DataFrame in-place
"""
function _apply_baseline!(
    dat::DataFrame,
    channel_indices::Vector{Int},
    baseline_interval
)

    # Compute the mean once for baseline interval
    baseline_means =
        mean.(eachcol(dat[baseline_interval.interval_start:baseline_interval.interval_end, channel_indices]))

    # Apply baseline correction
    for (col, mean_val) in zip(channel_labels, baseline_means)
        @inbounds dat[!, col] .-= mean_val
    end

end

"""
    baseline!(dat::Union{ContinuousData,ErpData}, channel_labels, baseline_interval)

Apply baseline correction in-place to continuous or ERP data.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Modifies the input data in-place by subtracting the baseline mean
"""
function baseline!(
    dat::Union{ContinuousData,ErpData},
    channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx,IntervalTime},
)
    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    channel_indices = get_channel_indices(dat, channel_labels)
    @info "Applying baseline correction to channel(s) $(print_vector_(channel_labels)) over IntervalIdx: $(baseline_interval.interval_start) to $(baseline_interval.interval_end)"
    _apply_baseline!(dat.data, channel_indices, baseline_interval)
end

"""
    baseline!(dat::Union{ContinuousData,ErpData}, channel_labels)

Apply baseline correction using entire time range as baseline interval.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct

# Effects
- Uses full time range for baseline calculation
"""
function baseline!(dat::Union{ContinuousData,ErpData}, channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}})
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, channel_labels, baseline_interval)
end

"""
    baseline!(dat::Union{ContinuousData,ErpData}, baseline_interval)

Apply baseline correction to all channels over specified interval.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Applies baseline correction to all channels
"""
function baseline!(dat::Union{ContinuousData,ErpData}, baseline_interval::Union{IntervalIdx,IntervalTime})
    baseline!(dat, dat.layout.label, baseline_interval)
end

"""
    baseline!(dat::Union{ContinuousData,ErpData})

Apply baseline correction to all channels using entire time range.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct

# Effects
- Uses all channels and full time range for baseline
"""
function baseline!(dat::Union{ContinuousData,ErpData})
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, dat.layout.label, baseline_interval)
end

"""
    baseline(dat::Union{ContinuousData,ErpData}, channel_labels, baseline_interval)

Create baseline-corrected copy of continuous or ERP data.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Returns
- New baseline-corrected copy of input data
"""
function baseline(
    dat::Union{ContinuousData,ErpData},
    channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx,IntervalTime},
)
    dat_out = deepcopy(dat)
    baseline!(dat_out, channel_labels, baseline_interval)
    return dat_out
end

"""
    baseline(dat::Union{ContinuousData,ErpData}, channel_labels)

Create baseline-corrected copy using entire time range as baseline interval.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct

# Returns
- New baseline-corrected copy of input data
"""
function baseline(dat::Union{ContinuousData,ErpData}, channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}})
    dat_out = deepcopy(dat)
    baseline!(dat_out, channel_labels)
    return dat_out
end

"""
    baseline(dat::Union{ContinuousData,ErpData}, baseline_interval)

Create baseline-corrected copy using all channels over specified interval.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Returns
- New baseline-corrected copy of input data
"""
function baseline(dat::Union{ContinuousData,ErpData}, baseline_interval::Union{IntervalIdx,IntervalTime})
    dat_out = deepcopy(dat)
    baseline!(dat_out, baseline_interval)
    return dat_out
end

"""
    baseline(dat::Union{ContinuousData,ErpData})

Create baseline-corrected copy using all channels and entire time range.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct

# Returns
- New baseline-corrected copy of input data
"""
function baseline(dat::Union{ContinuousData,ErpData})
    dat_out = deepcopy(dat)
    baseline!(dat_out)
    return dat_out
end

"""
    baseline!(dat::EpochData, channel_labels, baseline_interval)

Apply baseline correction to specified channels in epoched data.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Modifies each epoch in-place
"""
function baseline!(
    dat::EpochData,
    channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx,IntervalTime},
)
    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    channel_indices = get_channel_indices(dat, channel_labels)
    @info "Applying baseline correction to channel(s) $(print_vector_(channel_labels)) over over IntervalIdx: $(baseline_interval.interval_start) to $(baseline_interval.interval_end)"

    for epoch in eachindex(dat.data)
        _apply_baseline!(dat.data[epoch], channel_labels, baseline_interval)
    end
end

"""
    baseline!(dat::EpochData, channel_labels)

Apply baseline correction to specified channels in epoched data using entire time range.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct

# Effects
- Modifies each epoch in-place using full time range for baseline
"""
function baseline!(dat::EpochData, channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}})
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, channel_labels, baseline_interval)
end

"""
    baseline!(dat::EpochData, baseline_interval)

Apply baseline correction to all channels in epoched data over specified interval.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Effects
- Modifies all channels in each epoch in-place
"""
function baseline!(dat::EpochData, baseline_interval::Union{IntervalIdx,IntervalTime})
    baseline!(dat, dat.layout.label, baseline_interval)
end

"""
    baseline!(dat::EpochData)

Apply baseline correction to all channels in epoched data using entire time range.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct

# Effects
- Modifies all channels in each epoch in-place using full time range
"""
function baseline!(dat::EpochData)
    baseline!(dat, dat.layout.label, IntervalIdx(1, nrow(dat.data)))
end

"""
    baseline(dat::EpochData, channel_labels, baseline_interval)

Create baseline-corrected copy of epoched data.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Returns
- New baseline-corrected copy of input data
"""
function baseline(
    dat::EpochData,
    channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}},
    baseline_interval::Union{IntervalIdx,IntervalTime}
)
    dat_out = deepcopy(dat)
    baseline!(dat_out, channel_labels, baseline_interval)
    return dat_out
end

"""
    baseline(dat::EpochData, channel_labels)

Create baseline-corrected copy of epoched data using entire time range.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct

# Returns
- New baseline-corrected copy of input data
"""
function baseline(dat::EpochData, channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}})
    dat_out = deepcopy(dat)
    baseline!(dat_out, channel_labels, IntervalIdx(1, nrow(dat_out.data)))
    return dat_out
end

"""
    baseline(dat::EpochData, baseline_interval)

Create baseline-corrected copy of epoched data using all channels.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation

# Returns
- New baseline-corrected copy of input data
"""
function baseline(dat::EpochData, baseline_interval::Union{IntervalIdx,IntervalTime})
    dat_out = deepcopy(dat)
    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    baseline!(dat_out, dat_out.layout.label, baseline_interval)
    return dat_out
end

"""
    baseline(dat::EpochData)

Create baseline-corrected copy of epoched data using all channels and entire time range.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct

# Returns
- New baseline-corrected copy of input data
"""
function baseline(dat::EpochData)
    dat_out = deepcopy(dat)
    baseline!(dat_out, dat_out.layout.label, IntervalIdx(1, nrow(dat_out.data)))
    return dat_out
end

"""
    remove_mean!(dat::ContinuousData)

Remove the mean (DC offset) from each channel in continuous EEG data.

# Arguments
- `dat::ContinuousData`: The continuous EEG data to process

# Effects
- Modifies the input data in-place by subtracting the mean from each channel
"""
function remove_mean!(dat::ContinuousData)
    for channel in dat.layout.label
        dat.data[:, channel] .-= mean(dat.data[:, channel])
    end
    return dat
end