
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
    baseline!(dat::Union{ContinuousData,ErpData}, channel_labels)
    baseline!(dat::Union{ContinuousData,ErpData}, baseline_interval)
    baseline!(dat::Union{ContinuousData,ErpData})

Apply baseline correction in-place to continuous or ERP data.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: The data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct (optional)
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation (optional)

# Effects
- Modifies the input data in-place by subtracting the baseline mean
- If channel_labels omitted, uses all channels
- If baseline_interval omitted, uses entire time range
"""
function baseline!(
    dat::Union{ContinuousData,ErpData},
    channel_labels::Vector{Symbol},
    baseline_interval::Union{IntervalIdx,IntervalTime},
)
    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    channel_indices = get_channel_indices(dat, channel_labels)
    @info "Applying baseline correction to channel(s) $(print_vector_(channel_labels)) over IntervalIdx: $(baseline_interval.interval_start) to $(baseline_interval.interval_end)"
    _apply_baseline!(dat.data, channel_indices, baseline_interval)
end

function baseline!(dat::Union{ContinuousData,ErpData}, channel_labels::Vector{Symbol})
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, channel_labels, baseline_interval)
end

function baseline!(dat::Union{ContinuousData,ErpData}, baseline_interval::Union{IntervalIdx,IntervalTime})
    baseline!(dat, dat.layout.label, baseline_interval)
end

function baseline!(dat::Union{ContinuousData,ErpData})
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, dat.layout.label, baseline_interval)
end


"""
    baseline!(dat::EpochData, channel_labels, baseline_interval)
    baseline!(dat::EpochData, channel_labels)
    baseline!(dat::EpochData, baseline_interval)
    baseline!(dat::EpochData)

Apply baseline correction to epoched data.

# Arguments
- `dat::EpochData`: The epoched data to baseline correct
- `channel_labels::Union{Vector{Symbol},Vector{<:AbstractString}}`: Names of channels to correct (optional)
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Time interval for baseline calculation (optional)

# Effects
- Modifies each epoch in-place
- If channel_labels omitted, uses all channels
- If baseline_interval omitted, uses entire time range
"""
function baseline!(
    dat::EpochData,
    channel_labels::Vector{Symbol},
    baseline_interval::Union{IntervalIdx,IntervalTime},
)
    baseline_interval = validate_baseline_interval(dat.time, baseline_interval)
    channel_indices = get_channel_indices(dat, channel_labels)
    @info "Applying baseline correction to channel(s) $(print_vector_(channel_labels)) over over IntervalIdx: $(baseline_interval.interval_start) to $(baseline_interval.interval_end)"

    for epoch in eachindex(dat.data)
        _apply_baseline!(dat.data[epoch], channel_labels, baseline_interval)
    end
end

function baseline!(dat::EpochData, channel_labels::Vector{Symbol})
    baseline_interval = IntervalIdx(1, nrow(dat.data))
    baseline!(dat, channel_labels, baseline_interval)
end

function baseline!(dat::EpochData, baseline_interval::Union{IntervalIdx,IntervalTime})
    baseline!(dat, dat.layout.label, baseline_interval)
end

function baseline!(dat::EpochData)
    baseline!(dat, dat.layout.label, IntervalIdx(1, nrow(dat.data)))
end

# generates all non-mutating versions
@add_nonmutating baseline!
