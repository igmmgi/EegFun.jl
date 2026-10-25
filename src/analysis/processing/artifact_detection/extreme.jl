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
extreme_mask = _is_extreme_value(signal, 50.0) # Detect values ±50 μV
```
"""
_is_extreme_value(signal::AbstractVector{<:Real}, threshold::Real) = abs.(signal) .> threshold

"""
    _is_extreme_value!(mask::AbstractVector{Bool}, signal::AbstractVector{<:Real}, threshold::Real)

Detect extreme values in a signal and store results in a pre-allocated mask.

# Arguments
- `mask::AbstractVector{Bool}`: Pre-allocated boolean vector to store results
- `signal::AbstractVector{<:Real}`: Input signal vector
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
function _is_extreme_value!(mask::AbstractVector{Bool}, signal::AbstractVector{<:Real}, threshold::Real)
    length(mask) == length(signal) ||
        @minimal_error "Mask and signal must have the same length (mask: $(length(mask)), signal: $(length(signal)))"
    @inbounds for i in eachindex(signal)
        mask[i] = abs(signal[i]) > threshold
    end
    return nothing
end

"""
    is_extreme_value!(dat::SingleDataFrameEeg, threshold::Real; channel_selection = channels(),
    sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined, channel_out = nothing)
    is_extreme_value!(dat::MultiDataFrameEeg, threshold::Real; channel_selection = channels(),
    sample_selection = samples(), interval_selection = times(), epoch_selection = epochs(), mode::Symbol = :combined, channel_out = nothing)
    is_extreme_value!(epochs_list::Vector{EpochData}, threshold::Real; kwargs...)

Detect extreme values and flag them in-place. Adds a boolean column (`:is_extreme_value_<threshold>` by default)
to the data. The `Vector{EpochData}` form broadcasts across conditions.

# Keyword Arguments
- `channel_selection`: Channels to check (default: all layout channels)
- `mode::Symbol`: `:combined` (single OR'd column) or `:separate` (one column per channel)
- `channel_out`: Output column name (default: auto-generated)

# Examples
```julia
is_extreme_value!(dat, 100)             # → :is_extreme_value_100
is_extreme_value!(dat, 100, mode = :separate)
is_extreme_value!(epoch_data, 100, epoch_selection = epochs([1, 3]))
is_extreme_value!(epochs_unrejected, 100)  # batch over conditions
```
"""
function is_extreme_value!(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)

    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")

    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection = combined_sel)

    if mode == :combined  # any channel
        channel_out = something(channel_out, Symbol("is_extreme_value_$(threshold)"))
        dat.data[!, channel_out] = falses(n_samples(dat))
        for (_, extreme_mask) in results
            dat.data[!, channel_out] .|= extreme_mask
        end
    elseif mode == :separate # separate columns for each channel
        for (ch, extreme_mask) in results
            column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
            dat.data[!, column_name] = extreme_mask
        end
    end

    return dat
end


function is_extreme_value!(
    dat::MultiDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)

    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for extreme value detection")

    # Get selected epochs
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected for extreme value detection")

    # Use provided channel_out or generate default name
    channel_out = something(channel_out, Symbol("is_extreme_value_$(threshold)"))

    # Combine interval and sample selection once (constant across all epochs)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    # Process each selected epoch
    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]
        selected_samples = get_selected_samples(epoch_df, combined_sel)

        # Initialize artifact flag column for this epoch
        epoch_df[!, channel_out] = falses(nrow(epoch_df))

        # Create sample mask once (same for all channels)
        sample_mask = falses(nrow(epoch_df))
        sample_mask[selected_samples] .= true

        if mode == :combined
            for ch in selected_channels
                extreme_mask = _is_extreme_value(epoch_df[!, ch], threshold) .& sample_mask
                epoch_df[!, channel_out] .|= extreme_mask
            end
        else
            for ch in selected_channels
                extreme_mask = _is_extreme_value(epoch_df[!, ch], threshold) .& sample_mask
                column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
                epoch_df[!, column_name] = extreme_mask
            end
        end
    end

    return dat
end


function is_extreme_value!(dat::Vector{EpochData}, threshold::Real; kwargs...)
    for d in dat
        is_extreme_value!(d, threshold; kwargs...)
    end
    return dat
end


"""Detect extreme values for selected channels and return a `Dict{Symbol, Vector{Bool}}`."""
function _detect_extreme_values(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for extreme value detection")

    selected_samples = get_selected_samples(dat, sample_selection)
    isempty(selected_samples) && @minimal_error("No samples selected for extreme value detection")

    n_samp = n_samples(dat)
    sample_mask = falses(n_samp)
    sample_mask[selected_samples] .= true

    results = Dict{Symbol,Vector{Bool}}()
    thresh_f = Float64(threshold)

    for ch in selected_channels
        extreme_mask = _is_extreme_value(dat.data[!, ch], thresh_f)
        results[ch] = extreme_mask .& sample_mask
    end

    return results
end

"""
    is_extreme_value(dat::SingleDataFrameEeg, threshold::Real; channel_selection = channels(),
    sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined)

Detect extreme values without modifying data. Returns `Vector{Bool}` (`:combined` mode)
or a `DataFrame` (`:separate` mode).

# Examples
```julia
extreme_mask = is_extreme_value(dat, 100)
results = is_extreme_value(dat, 100, mode = :separate)
```
"""
function is_extreme_value(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)

    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection = combined_sel)

    if mode == :combined

        combined_mask = Vector{Bool}(falses(n_samples(dat)))
        # Combine results from all channels (OR operation)
        for (_, extreme_mask) in results
            combined_mask .|= extreme_mask
        end

        return combined_mask

    elseif mode == :separate
        # Separate mode - create temporary data object and use mutating version
        temp_dat = copy(dat)
        is_extreme_value!(temp_dat, threshold; channel_selection, sample_selection, interval_selection, mode = :separate)

        # Extract the extreme value columns in the same order as the original channels
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        extreme_cols = [Symbol("is_extreme_value_$(ch)_$(threshold)") for ch in selected_channels]
        return temp_dat.data[!, extreme_cols]
    end
end

"""
    n_extreme_value(dat::SingleDataFrameEeg, threshold::Real; 
                   channel_selection::Function = channels(), 
                   sample_selection::Function = samples(),
    interval_selection::Interval = times(),
                   mode::Symbol = :combined)

Count the number of extreme values across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for extreme value detection
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
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
function n_extreme_value(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)

    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection = combined_sel)

    if mode == :combined
        combined_mask = Vector{Bool}(falses(n_samples(dat)))
        for (_, extreme_mask) in results
            combined_mask .|= extreme_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        counts = [sum(results[ch]) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_extreme = counts)
    end
end

"""
    n_values(df::DataFrame, column::Symbol)
    n_values(dat::SingleDataFrameEeg, column::Symbol)
    n_values(dat::MultiDataFrameEeg, column::Symbol)
    n_values(dat::Vector{<:EegData}, column::Symbol)

Count the number of `true` values in a boolean column. Returns 0 if the column is not found.
For `MultiDataFrameEeg` and `Vector` inputs, counts are summed across all epochs/items.
"""
function n_values(df::DataFrame, column::Symbol)
    column ∉ propertynames(df) && return 0
    return sum(df[!, column])
end

n_values(dat::SingleDataFrameEeg, column::Symbol) = n_values(dat.data, column)

function n_values(dat::MultiDataFrameEeg, column::Symbol)
    return sum(n_values(epoch, column) for epoch in dat)
end

n_values(dat::Vector{<:EegData}, column::Symbol) = sum(n_values(d, column) for d in dat)

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
function _n_extreme_value(df::DataFrame, channels::AbstractVector{Symbol}, threshold::Float64)
    counts = Int[]
    for ch in channels
        channel_data = df[!, ch]
        extreme_mask = _is_extreme_value(channel_data, threshold)
        push!(counts, sum(extreme_mask))
    end
    return counts
end
