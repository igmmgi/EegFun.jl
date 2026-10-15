# === STEP/JUMP ARTIFACT DETECTION ===

"""
    _is_step_value(signal::AbstractVector{Float64}, threshold::Real)

Detect step artifacts (sudden jumps) in a signal by computing differences between consecutive samples.

A step artifact is detected when the absolute difference between consecutive samples exceeds the threshold.
This is useful for detecting cable disconnections, amplifier saturation, or sudden movement artifacts.

# Arguments
- `signal::AbstractVector{Float64}`: Input signal vector
- `threshold::Float64`: Threshold for step detection (in same units as signal, typically μV)

# Returns
- `Vector{Bool}`: Boolean vector indicating samples where a step occurred (second sample of the pair)

# Examples
```julia
step_mask = _is_step_value(signal, 50.0) # Detect jumps > 50 μV between samples
```

# Notes
- The first sample is always false (no previous sample to compare)
- A step at index i means the jump from i-1 to i exceeded threshold
"""
function _is_step_value(signal::AbstractVector{<:Real}, threshold::Real)
    n = length(signal)
    mask = Vector{Bool}(undef, n)
    mask[1] = false  # First sample has no previous sample

    @inbounds for i = 2:n
        mask[i] = abs(signal[i] - signal[i-1]) > threshold
    end

    return mask
end

"""
    is_step_value!(dat::SingleDataFrameEeg, threshold::Real; channel_selection = channels(),
    sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined, channel_out = nothing)
    is_step_value!(dat::MultiDataFrameEeg, threshold::Real; channel_selection = channels(),
    sample_selection = samples(), interval_selection = times(), epoch_selection = epochs(), mode::Symbol = :combined, channel_out = nothing)
    is_step_value!(epochs_list::Vector{EpochData}, threshold::Real; kwargs...)

Detect step artifacts (sudden voltage jumps) and flag them in-place. Adds a boolean column
(`:is_step_value_<threshold>` by default). The `Vector{EpochData}` form broadcasts across conditions.

# Examples
```julia
is_step_value!(dat, 50.0)             # → :is_step_value_50.0
is_step_value!(dat, 50.0, mode = :separate)
is_step_value!(epoch_data, 50.0, epoch_selection = epochs([1, 3]))
is_step_value!(all_epochs, 50.0)      # batch over conditions
```
"""
function is_step_value!(
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

    # Get selected channels (exclude metadata and extra columns like triggers, EOG flags)
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    # Get selected samples (combining interval and sample selection)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)

    # Use provided channel_out or generate default name
    channel_out = something(channel_out, Symbol("is_step_value_$(threshold)"))

    if mode == :combined
        # Initialize combined output column
        dat.data[!, channel_out] = falses(n_samples(dat))

        # Build boolean sample mask from indices
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true

        # Check each channel and combine results
        for ch in selected_channels
            step_mask = _is_step_value(dat.data[!, ch], threshold) .& sample_mask
            dat.data[!, channel_out] .|= step_mask
        end
    elseif mode == :separate
        # Build boolean sample mask from indices
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true

        # Create separate column for each channel
        for ch in selected_channels
            step_mask = _is_step_value(dat.data[!, ch], threshold) .& sample_mask
            column_name = Symbol("is_step_value_$(ch)_$(threshold)")
            dat.data[!, column_name] = step_mask
        end
    end

    return nothing
end


function is_step_value!(
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

    # Get selected channels (exclude metadata and extra columns) and epochs
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    # Use provided channel_out or generate default name
    channel_out = something(channel_out, Symbol("is_step_value_$(threshold)"))

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
                step_mask = _is_step_value(epoch_df[!, ch], threshold) .& sample_mask
                epoch_df[!, channel_out] .|= step_mask
            end
        else
            for ch in selected_channels
                step_mask = _is_step_value(epoch_df[!, ch], threshold) .& sample_mask
                column_name = Symbol("is_step_value_$(ch)_$(threshold)")
                epoch_df[!, column_name] = step_mask
            end
        end
    end

    return nothing
end


is_step_value!(dat::Vector{EpochData}, threshold::Real; kwargs...) = is_step_value!.(dat, threshold; kwargs...)


"""
    is_step_value(dat::SingleDataFrameEeg, threshold::Real; 
                  channel_selection::Function = channels(), 
                  sample_selection::Function = samples(),
                  interval_selection::Interval = times(),
                  mode::Symbol = :combined)

Detect step values (sudden voltage jumps) across selected channels and return results
without modifying the input data.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for step detection (in μV)
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
- `interval_selection::Interval`: Time interval for selection (default: all times)
- `mode::Symbol`: Mode of operation - `:combined` (boolean vector, default) or `:separate` (DataFrame with separate columns per channel)

# Returns
- `Vector{Bool}` (combined mode) or `DataFrame` (separate mode)

# Examples
```julia
# Detect step values and return combined boolean vector (default)
step_mask = is_step_value(dat, 50.0)

# Detect step values and return DataFrame with separate columns
results = is_step_value(dat, 50.0, mode = :separate)

# Detect step values for specific channels
step_mask = is_step_value(dat, 50.0, channel_selection = channels([:Fp1, :Fp2]))
```
"""
function is_step_value(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)

    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")

    # Get selected channels (exclude metadata and extra columns) and samples
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)

    if mode == :combined
        combined_mask = falses(n_samples(dat))
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true
        for ch in selected_channels
            step_mask = _is_step_value(dat.data[!, ch], threshold) .& sample_mask
            combined_mask .|= step_mask
        end
        return Vector{Bool}(combined_mask)

    elseif mode == :separate
        temp_dat = copy(dat)
        is_step_value!(temp_dat, threshold; channel_selection, sample_selection, interval_selection, mode = :separate)

        extreme_cols = [Symbol("is_step_value_$(ch)_$(threshold)") for ch in selected_channels]
        return temp_dat.data[!, extreme_cols]
    end
end


"""
    n_step_value(dat::SingleDataFrameEeg, threshold::Real; 
                 channel_selection::Function = channels(), 
                 sample_selection::Function = samples(),
                 interval_selection::Interval = times(),
                 mode::Symbol = :combined)

Count the number of step values (sudden voltage jumps) across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for step value detection
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
- `mode::Symbol`: Mode for step value detection (`:separate` or `:combined`, default: `:combined`)

# Returns
- `Int` (combined mode): Total number of samples with step values across all selected channels
- `DataFrame` (separate mode): DataFrame with step value counts for each channel

# Examples
```julia
# Count total step values across all channels (combined mode, default)
total_count = n_step_value(dat, 50)

# Count step values per channel (separate mode)
count_df = n_step_value(dat, 50, mode = :separate)

# Count step values in specific channels
total_count = n_step_value(dat, 50, channel_selection = channels([:Fp1, :Fp2]))
```
"""
function n_step_value(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)

    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)
    sample_mask = falses(n_samples(dat))
    sample_mask[selected_samples] .= true

    if mode == :combined
        combined_mask = falses(n_samples(dat))
        for ch in selected_channels
            step_mask = _is_step_value(dat.data[!, ch], threshold) .& sample_mask
            combined_mask .|= step_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        counts = [sum(_is_step_value(dat.data[!, ch], threshold) .& sample_mask) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_step = counts)
    end
end

