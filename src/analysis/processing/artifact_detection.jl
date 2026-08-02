"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Float64, channel_in::Symbol, channel_out::Symbol)

Detects EOG (electrooculogram) onsets in the EEG data based on a specified criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Real`: The threshold for detecting EOG onsets.
- `channel_in::Symbol`: The channel from which to detect EOG onsets.
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event.

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG) # Detect vertical EOG onsets
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG) # Detect horizontal EOG onsets
```
"""

function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol; step_size::Int = 20)
    @info "Detecting EOG onsets in channel $(channel_in) with stepsize criterion $(criterion) μV"
    channel_in ∉ propertynames(dat.data) && @minimal_error("channel $(channel_in) not found in data")

    step_samples = div(dat.sample_rate, step_size)
    eog_diff = diff(dat.data[1:step_samples:end, channel_in])
    eog_idx = findall(x -> abs(x) >= criterion, eog_diff)
    eog_idx = [idx for (i, idx) in enumerate(eog_idx) if i == 1 || (idx - eog_idx[i-1] > 2)] .* step_samples
    dat.data[!, channel_out] = falses(nrow(dat.data))
    dat.data[eog_idx, channel_out] .= true
    return nothing
end


"""
    detect_eog_onsets!(dat::EegData, eog_cfg::Dict)

Detect EOG onsets for both vertical and horizontal EOG channels based on configuration.

# Arguments
- `dat::EegData`: The EEG data object
- `eog_cfg::Dict`: EOG configuration dictionary containing vEOG and hEOG settings

# Example
```julia
eog_cfg = Dict(
    "vEOG_criterion" => 50.0,
    "hEOG_criterion" => 50.0,
    "vEOG_channels" => [["Fp1"], ["IO1"], ["vEOG"]],
    "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]]
)
detect_eog_onsets!(dat, eog_cfg)
```
"""
function detect_eog_onsets!(dat::EegData, eog_cfg::Dict)
    vEOG_cfg = eog_cfg["vEOG_channels"]
    detect_eog_onsets!(dat, eog_cfg["vEOG_criterion"], Symbol(vEOG_cfg[3][1]), Symbol("is_" * vEOG_cfg[3][1]))
    hEOG_cfg = eog_cfg["hEOG_channels"]
    detect_eog_onsets!(dat, eog_cfg["hEOG_criterion"], Symbol(hEOG_cfg[3][1]), Symbol("is_" * hEOG_cfg[3][1]))
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
extreme_mask = _is_extreme_value(signal, 50.0) # Detect values ±50 μV
```
"""
_is_extreme_value(signal::AbstractVector{Float64}, threshold::Real) = abs.(signal) .> threshold

"""
    _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Real)

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
function _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Real)
    @assert length(mask) == length(signal) "Mask and signal must have the same length"
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
        dat.data[!, channel_out] = falses(nrow(dat.data))
        for (_, extreme_mask) in results
            dat.data[!, channel_out] .|= extreme_mask
        end
    elseif mode == :separate # separate columns for each channel
        for (ch, extreme_mask) in results
            column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
            dat.data[!, column_name] = extreme_mask
        end
    end

    return nothing
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

    # Process each selected epoch
    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]

        # Get selected samples for this epoch (combining interval and sample selection)
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
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

    return nothing
end


is_extreme_value!(dat::Vector{EpochData}, threshold::Real; kwargs...) = is_extreme_value!.(dat, threshold; kwargs...)

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
        dat.data[!, channel_out] = falses(nrow(dat.data))

        # Build boolean sample mask from indices
        sample_mask = falses(nrow(dat.data))
        sample_mask[selected_samples] .= true

        # Check each channel and combine results
        for ch in selected_channels
            step_mask = _is_step_value(dat.data[!, ch], threshold) .& sample_mask
            dat.data[!, channel_out] .|= step_mask
        end
    elseif mode == :separate
        # Build boolean sample mask from indices
        sample_mask = falses(nrow(dat.data))
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

    # Process each selected epoch
    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]

        # Get selected samples for this epoch (combining interval and sample selection)
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
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
        combined_mask = falses(nrow(dat.data))
        sample_mask = falses(nrow(dat.data))
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
    sample_mask = falses(nrow(dat.data))
    sample_mask[selected_samples] .= true

    if mode == :combined
        combined_mask = falses(nrow(dat.data))
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

    results = Dict{Symbol,Vector{Bool}}()

    for ch in selected_channels
        extreme_mask = _is_extreme_value(dat.data[!, ch], Float64(threshold))

        # Apply sample selection - only keep extreme values for selected samples
        sample_mask = falses(length(extreme_mask))
        sample_mask[selected_samples] .= true
        extreme_mask = extreme_mask .& sample_mask

        results[ch] = extreme_mask
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

        combined_mask = Vector{Bool}(falses(nrow(dat.data)))
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
        combined_mask = Vector{Bool}(falses(nrow(dat.data)))
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
    return sum(n_values(epoch, column) for epoch in dat.data)
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


# === FLATLINE (LOW VARIANCE) DETECTION ===

"""
    _is_flatline(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)

Detect flatline artifacts (signal stuck with extremely low variance) using a sliding window.
Returns a boolean mask where all samples within any window whose std dev is below the threshold are true.
"""
function _is_flatline(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)
    n = length(signal)
    mask = falses(n)
    w = max(1, window_samples)

    @inbounds for i = 1:(n-w+1)
        w_view = view(signal, i:(i+w-1))
        if std(w_view) < threshold
            mask[i:(i+w-1)] .= true
        end
    end
    return mask
end

"""
    is_flatline!(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                 channel_selection::Function = channels(),
                 sample_selection::Function = samples(),
                 interval_selection::Interval = times(),
                 mode::Symbol = :combined, channel_out = nothing)
    is_flatline!(dat::MultiDataFrameEeg, ...)
    is_flatline!(dat::Vector{EpochData}, ...)

Detect flatline artifacts using a sliding window of `window_size` (in seconds).
Flags all samples in windows where standard deviation < `threshold`.
"""
function is_flatline!(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")
    window_size <= 0 && @minimal_error("window_size must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)

    window_samples = round(Int, window_size * dat.sample_rate)
    channel_out = something(channel_out, :is_flatline)

    if mode == :combined
        dat.data[!, channel_out] = falses(nrow(dat.data))
        sample_mask = falses(nrow(dat.data))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            flat_mask = _is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask
            dat.data[!, channel_out] .|= flat_mask
        end
    elseif mode == :separate
        sample_mask = falses(nrow(dat.data))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            flat_mask = _is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask
            column_name = Symbol("$(channel_out)_$(ch)")
            dat.data[!, column_name] = flat_mask
        end
    end

    return nothing
end

function is_flatline!(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")
    window_size <= 0 && @minimal_error("window_size must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    window_samples = round(Int, window_size * dat.sample_rate)
    channel_out = something(channel_out, :is_flatline)

    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        selected_samples = get_selected_samples(epoch_df, combined_sel)

        epoch_df[!, channel_out] = falses(nrow(epoch_df))
        sample_mask = falses(nrow(epoch_df))
        sample_mask[selected_samples] .= true

        if mode == :combined
            for ch in selected_channels
                flat_mask = _is_flatline(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                epoch_df[!, channel_out] .|= flat_mask
            end
        else
            for ch in selected_channels
                flat_mask = _is_flatline(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                column_name = Symbol("$(channel_out)_$(ch)")
                epoch_df[!, column_name] = flat_mask
            end
        end
    end
    return nothing
end

is_flatline!(dat::Vector{EpochData}, threshold::Real, window_size::Real; kwargs...) = is_flatline!.(dat, threshold, window_size; kwargs...)

"""
    is_flatline(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                channel_selection = channels(), sample_selection = samples(),
                interval_selection = times(), mode::Symbol = :combined)

Detect flatline artifacts without modifying data. Returns `Vector{Bool}` (`:combined` mode)
or a `DataFrame` (`:separate` mode).
"""
function is_flatline(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    temp_dat = copy(dat)
    is_flatline!(temp_dat, threshold, window_size; channel_selection, sample_selection, interval_selection, mode)

    if mode == :combined
        return temp_dat.data[!, :is_flatline]
    elseif mode == :separate
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        cols = [Symbol("is_flatline_$(ch)") for ch in selected_channels]
        return temp_dat.data[!, cols]
    end
end

"""
    n_flatline(dat::SingleDataFrameEeg, threshold::Real, window_size::Real; ...)

Count the number of flatline values across selected channels.
"""
function n_flatline(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)
    sample_mask = falses(nrow(dat.data))
    sample_mask[selected_samples] .= true
    window_samples = round(Int, window_size * dat.sample_rate)

    if mode == :combined
        combined_mask = falses(nrow(dat.data))
        for ch in selected_channels
            flat_mask = _is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask
            combined_mask .|= flat_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        counts = [sum(_is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_flatline = counts)
    end
end

# === PEAK TO PEAK ARTIFACT DETECTION ===

"""
    _is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)

Detect moving-window peak-to-peak artifacts (max - min > threshold).
Returns a boolean mask where all samples within any offending window are true.
"""
function _is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)
    n = length(signal)
    mask = falses(n)
    w = max(1, window_samples)

    @inbounds for i = 1:(n-w+1)
        w_view = view(signal, i:(i+w-1))
        if (maximum(w_view) - minimum(w_view)) > threshold
            mask[i:(i+w-1)] .= true
        end
    end
    return mask
end

"""
    is_peak_to_peak!(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                     channel_selection::Function = channels(),
                     sample_selection::Function = samples(),
                     interval_selection::Interval = times(),
                     mode::Symbol = :combined, channel_out = nothing)
    is_peak_to_peak!(dat::MultiDataFrameEeg, ...)
    is_peak_to_peak!(dat::Vector{EpochData}, ...)

Detect moving window peak-to-peak artifacts (where max-min > `threshold` inside `window_size` in seconds).
"""
function is_peak_to_peak!(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")
    window_size <= 0 && @minimal_error("window_size must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)

    window_samples = round(Int, window_size * dat.sample_rate)
    channel_out = something(channel_out, :is_peak_to_peak)

    if mode == :combined
        dat.data[!, channel_out] = falses(nrow(dat.data))
        sample_mask = falses(nrow(dat.data))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            dat.data[!, channel_out] .|= p2p_mask
        end
    elseif mode == :separate
        sample_mask = falses(nrow(dat.data))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            column_name = Symbol("$(channel_out)_$(ch)")
            dat.data[!, column_name] = p2p_mask
        end
    end

    return nothing
end

function is_peak_to_peak!(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")
    window_size <= 0 && @minimal_error("window_size must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    window_samples = round(Int, window_size * dat.sample_rate)
    channel_out = something(channel_out, :is_peak_to_peak)

    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        selected_samples = get_selected_samples(epoch_df, combined_sel)

        epoch_df[!, channel_out] = falses(nrow(epoch_df))
        sample_mask = falses(nrow(epoch_df))
        sample_mask[selected_samples] .= true

        if mode == :combined
            for ch in selected_channels
                p2p_mask = _is_peak_to_peak(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                epoch_df[!, channel_out] .|= p2p_mask
            end
        else
            for ch in selected_channels
                p2p_mask = _is_peak_to_peak(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                column_name = Symbol("$(channel_out)_$(ch)")
                epoch_df[!, column_name] = p2p_mask
            end
        end
    end
    return nothing
end

is_peak_to_peak!(dat::Vector{EpochData}, threshold::Real, window_size::Real; kwargs...) =
    is_peak_to_peak!.(dat, threshold, window_size; kwargs...)

"""
    is_peak_to_peak(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                    channel_selection = channels(), sample_selection = samples(),
                    interval_selection = times(), mode::Symbol = :combined)

Detect moving window peak-to-peak artifacts without modifying data.
Returns `Vector{Bool}` (`:combined` mode) or a `DataFrame` (`:separate` mode).
"""
function is_peak_to_peak(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    temp_dat = copy(dat)
    is_peak_to_peak!(temp_dat, threshold, window_size; channel_selection, sample_selection, interval_selection, mode)

    if mode == :combined
        return temp_dat.data[!, :is_peak_to_peak]
    elseif mode == :separate
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        cols = [Symbol("is_peak_to_peak_$(ch)") for ch in selected_channels]
        return temp_dat.data[!, cols]
    end
end

"""
    n_peak_to_peak(dat::SingleDataFrameEeg, threshold::Real, window_size::Real; ...)

Count the number of peak-to-peak artifact values across selected channels.
"""
function n_peak_to_peak(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)
    sample_mask = falses(nrow(dat.data))
    sample_mask[selected_samples] .= true
    window_samples = round(Int, window_size * dat.sample_rate)

    if mode == :combined
        combined_mask = falses(nrow(dat.data))
        for ch in selected_channels
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            combined_mask .|= p2p_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        counts = [sum(_is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_peak_to_peak = counts)
    end
end

# === BRIDGED CHANNEL DETECTION ===

"""
    find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())
    find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())

Find channels that are highly correlated with each other, indicative of an electrolyte gel bridge.
Returns a `DataFrame` of the bridged pairs and their correlation coefficients.
"""
function find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection::Function = channels())
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)

    pairs = Tuple{Symbol,Symbol}[]
    correlations = Float64[]

    if n_ch < 2
        return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])
    end

    # Calculate pairwise correlation
    for i = 1:(n_ch-1)
        for j = (i+1):n_ch
            ch1 = selected_channels[i]
            ch2 = selected_channels[j]
            corr = cor(dat.data[!, ch1], dat.data[!, ch2])
            if corr > correlation_threshold
                push!(pairs, (ch1, ch2))
                push!(correlations, corr)
            end
        end
    end

    return DataFrame(channel_1 = first.(pairs), channel_2 = last.(pairs), correlation = correlations)
end

function find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection::Function = channels())
    # For epoched data, we flatten the data vertically and compute correlation across all epochs
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)

    if n_ch < 2
        return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])
    end

    # Extract channel data across all epochs
    flattened_data = Dict{Symbol,Vector{Float64}}()
    for ch in selected_channels
        flattened_data[ch] = Float64[]
    end

    for epoch_df in dat.data
        for ch in selected_channels
            append!(flattened_data[ch], epoch_df[!, ch])
        end
    end

    pairs = Tuple{Symbol,Symbol}[]
    correlations = Float64[]

    for i = 1:(n_ch-1)
        for j = (i+1):n_ch
            ch1 = selected_channels[i]
            ch2 = selected_channels[j]
            corr = cor(flattened_data[ch1], flattened_data[ch2])
            if corr > correlation_threshold
                push!(pairs, (ch1, ch2))
                push!(correlations, corr)
            end
        end
    end

    return DataFrame(channel_1 = first.(pairs), channel_2 = last.(pairs), correlation = correlations)
end

# Automatic epoch rejection based on statistical criteria.
#
# Automatically rejects epochs based on statistical measures 
# (variance, max, min, absolute max, range, kurtosis) using z-score thresholds. 
# This is useful for removing epochs with artifacts without manual inspection.
struct Rejection
    channel::Symbol
    epoch::Int
end

"""Pretty-print a `Rejection`, showing channel and epoch."""
Base.show(io::IO, r::Rejection) = print(io, "Rejection(:$(r.channel), $(r.epoch))")

"""Check whether two `Rejection` objects refer to the same channel and epoch."""
is_equal_rejection(a::Rejection, b::Rejection) = a.channel == b.channel && a.epoch == b.epoch

"""De-duplicate a vector of `Rejection`s by `(channel, epoch)`."""
function unique_rejections(rejections::Vector{Rejection})
    seen = Set{Tuple{Symbol,Int}}()
    out = Rejection[]
    for rejection in rejections
        key = (rejection.channel, rejection.epoch)
        if key ∉ seen
            push!(seen, key)
            push!(out, rejection)
        end
    end
    return out
end

"""Return the distinct channel symbols across a vector of rejections."""
unique_channels(rejections::Vector{Rejection}) = unique(map(x -> x.channel, rejections))
"""Return the distinct epoch indices across a vector of rejections."""
unique_epochs(rejections::Vector{Rejection}) = unique(map(x -> x.epoch, rejections))


"""
    EpochInfo

Stores metadata about the epoch condition.

# Fields
- `number::Int`: Condition number from the epoch data
- `name::String`: Condition name from the epoch data
- `n::Int`: Number of epochs before rejection
"""
struct EpochInfo
    number::Int
    name::String
    n::Int
end

"""
    ZScoreRejectionInfo

Stores z-score based rejection information.

# Fields
- `z_measures::Vector{Symbol}`: Which z-score measures were evaluated
- `z_variance::Vector{Rejection}`: Rejections due to high variance
- `z_max::Vector{Rejection}`: Rejections due to high maximum values
- `z_min::Vector{Rejection}`: Rejections due to low minimum values
- `z_abs::Vector{Rejection}`: Rejections due to high absolute values
- `z_range::Vector{Rejection}`: Rejections due to large range
- `z_kurtosis::Vector{Rejection}`: Rejections due to high kurtosis
"""
struct ZScoreRejectionInfo
    z_measures::Vector{Symbol}
    z_variance::Vector{Rejection}
    z_max::Vector{Rejection}
    z_min::Vector{Rejection}
    z_abs::Vector{Rejection}
    z_range::Vector{Rejection}
    z_kurtosis::Vector{Rejection}
end

"""
    EpochRejectionInfo

Stores information about which epochs were rejected and why, and optionally tracks channel repairs.

# Fields
- `name::String`: Name/identifier for this rejection info (e.g., "rejection_step1", "rejection_step2")
- `info::EpochInfo`: Condition metadata (number, name, n_epochs)
- `n_artifacts::Int`: Total number of artifact detections (channel-epoch pairs)
- `abs_criterion::Real`: Absolute voltage threshold (μV) used for rejection
- `abs_rejections::Union{Vector{Rejection}, Nothing}`: Rejections due to absolute voltage threshold (Nothing if abs_criterion = 0)
- `z_criterion::Real`: Z-score criterion used for rejection
- `z_rejections::Union{ZScoreRejectionInfo, Nothing}`: Z-score based rejection info (Nothing if z_criterion = 0)
- `rejected::Vector{Rejection}`: All rejected epochs
- `repaired::Union{OrderedDict{Int, Vector{Symbol}}, Nothing}`: Channels repaired per epoch, ordered by epoch number (populated during repair, Nothing if no repairs)
- `skipped::Union{OrderedDict{Int, Vector{Symbol}}, Nothing}`: Channels skipped per epoch, ordered by epoch number (populated during repair, Nothing if no repairs)
"""
mutable struct EpochRejectionInfo
    name::String
    info::EpochInfo
    n_artifacts::Int
    abs_criterion::Real
    abs_rejections::Union{Vector{Rejection},Nothing}
    z_criterion::Real
    z_rejections::Union{ZScoreRejectionInfo,Nothing}
    rejected::Vector{Rejection}
    repaired::Union{OrderedDict{Int,Vector{Symbol}},Nothing}
    skipped::Union{OrderedDict{Int,Vector{Symbol}},Nothing}
end

"""Return unique rejections from an `EpochRejectionInfo`."""
unique_rejections(info::EpochRejectionInfo) = unique_rejections(info.rejected)
"""Return unique rejections for each `EpochRejectionInfo` in a vector."""
function unique_rejections(infos::Vector{EpochRejectionInfo})
    results = Vector{Rejection}[]
    for info in infos
        push!(results, unique_rejections(info))
    end
    return results
end

"""Return unique channels from an `EpochRejectionInfo` or vector thereof."""
unique_channels(rejections::EpochRejectionInfo) = unique_channels(rejections.rejected)
unique_channels(info::Vector{EpochRejectionInfo}) = unique_channels.(info)

"""Return unique epochs from an `EpochRejectionInfo` or vector thereof."""
unique_epochs(rejections::EpochRejectionInfo) = unique_epochs(rejections.rejected)
unique_epochs(info::Vector{EpochRejectionInfo}) = unique_epochs.(info)

"""Flatten all rejection data (abs + z-score) into a single DataFrame with reason tags."""
function all_data(info::EpochRejectionInfo)::DataFrame
    rows = NamedTuple{(:condition_number, :condition_name, :epoch, :channel, :reason),Tuple{Int,String,Int,Symbol,Symbol}}[]
    cnum = info.info.number
    cname = info.info.name
    _push!(rs, tag, list) =
        foreach(r -> push!(rs, (condition_number = cnum, condition_name = cname, epoch = r.epoch, channel = r.channel, reason = tag)), list)
    !isnothing(info.abs_rejections) && _push!(rows, :abs_threshold, info.abs_rejections)
    if !isnothing(info.z_rejections)
        z = info.z_rejections
        _push!(rows, :z_variance, z.z_variance)
        _push!(rows, :z_max, z.z_max)
        _push!(rows, :z_min, z.z_min)
        _push!(rows, :z_abs, z.z_abs)
        _push!(rows, :z_range, z.z_range)
        _push!(rows, :z_kurtosis, z.z_kurtosis)
    end
    return DataFrame(rows)
end
all_data(infos::Vector{EpochRejectionInfo})::DataFrame = vcat(all_data.(infos)...)


"""
    detect_bad_epochs_automatic(dat::EpochData; z_criterion::Real = 3,
                     abs_criterion::Real = 100,
                     channel_selection::Function = channels(),
                     z_measures::Vector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis])::EpochRejectionInfo

Detect bad epochs using statistical criteria and/or absolute voltage thresholds.

Supports two types of criteria for epoch rejection:
1. **Z-score criteria**: Calculates six statistical measures for each epoch (variance, max, min,
   absolute max, range, kurtosis) across selected channels. For each measure, the maximum
   across channels is taken for each epoch, then z-scored across epochs. Epochs exceeding
   the z-criterion for any selected measure are rejected. You can choose which measures to
   apply using `measures` (default: `[:variance, :max, :min, :abs, :range, :kurtosis]`).
2. **Absolute voltage criteria**: Epochs with any channel exceeding the absolute voltage
   threshold (in μV) are rejected.

# Arguments
- `dat::EpochData`: Epoched EEG data to process
- `z_criterion::Real`: Z-score threshold for rejection (default: 3.0). Set to 0 to disable z-score based rejection.
- `abs_criterion::Real`: Absolute voltage threshold in μV for rejection (default: 100.0). Set to 0 to disable absolute threshold rejection.
- `channel_selection::Function`: Channel predicate for selecting channels to analyze (default: all channels)
- `z_measures::Vector{Symbol}`: Which z-score measures to apply (default: all: `[:variance, :max, :min, :abs, :range, :kurtosis]`)

# Returns
- `EpochRejectionInfo`: Information about which epochs were rejected and why

# Requirements
- At least one of `z_criterion` or `abs_criterion` must be greater than 0

# Examples
```julia
using EegFun, JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2")

# Use default criteria (z_criterion=3.0, abs_criterion=100.0)
rejection_info = detect_bad_epochs_automatic(epochs)

# Customize z-score criterion only
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 2.0)

# Customize absolute voltage threshold only
rejection_info = detect_bad_epochs_automatic(epochs, abs_criterion = 80.0)  # 80 μV

# Use both criteria with custom values
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 2.5, abs_criterion = 80.0)

# Disable z-score rejection (use only absolute threshold)
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100.0)

# Disable absolute threshold rejection (use only z-score)
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 3.0, abs_criterion = 0)

# Check results
println("Original epochs: \$(rejection_info.n_epochs)")
println("Rejected epochs: \$(rejection_info.n_artifacts)")
println("Rejected epochs: \$(rejection_info.rejected)")
```

# Notes
- Default z-criterion: 3.0 
- Default abs-criterion: 100 μV 
- Common z-criteria: 2.0 (more aggressive), 2.5, 3.0 (more conservative)
- Common absolute criteria: 50-100 μV 
- Set either criterion to 0 to disable that type of rejection
- Rejection is based on ANY metric exceeding the criterion
- All metrics are calculated independently and combined with OR logic
"""
function detect_bad_epochs_automatic(
    dat::EpochData;
    z_criterion::Real = 3,
    abs_criterion::Real = 100,
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    z_measures::AbstractVector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis],
    name::String = "rejection_info",
)::EpochRejectionInfo

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Detecting bad epochs"

    # Validate inputs
    z_criterion < 0 && @minimal_error("Z-criterion must be non-negative")
    abs_criterion < 0 && @minimal_error("Absolute criterion must be non-negative")
    z_criterion == 0 && abs_criterion == 0 && @minimal_error("One of z_criterion or abs_criterion must be > 0")

    # Validate measures
    allowed_measures = Set([:variance, :max, :min, :abs, :range, :kurtosis])
    selected_measures = Set(z_measures)
    if !issubset(selected_measures, allowed_measures)
        invalid = collect(setdiff(selected_measures, allowed_measures))
        @minimal_error("Invalid measures: $(invalid).")
    end

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    @info "Selected channels: $(_print_vector(selected_channels))"
    isempty(selected_channels) && @minimal_error("No channels selected for epoch rejection")

    # Get sample indices for the selected interval
    selected_samples = get_selected_samples(dat, interval_selection)
    @info "Selected interval: $(length(selected_samples)) samples"
    isempty(selected_samples) && @minimal_error("No samples selected for epoch rejection")

    # Calculate metrics and identify rejected epochs
    metrics = _calculate_epoch_metrics(dat, selected_channels, selected_samples, Float64(z_criterion), Float64(abs_criterion))

    # Initialize rejection lists (all needed for EpochRejectionInfo struct)
    rejected_info = Rejection[]
    z_variance = Rejection[]
    z_max = Rejection[]
    z_min = Rejection[]
    z_abs = Rejection[]
    z_range = Rejection[]
    z_kurtosis = Rejection[]

    # Map measure directly to (rejection_list, metrics_key)
    measure_to_list_and_key = Dict(
        :variance => (z_variance, :z_variance),
        :max => (z_max, :z_max),
        :min => (z_min, :z_min),
        :abs => (z_abs, :z_abs),
        :range => (z_range, :z_range),
        :kurtosis => (z_kurtosis, :z_kurtosis),
    )

    # Build Rejection objects directly from metric results
    # Handle z-score metrics
    if z_criterion > 0
        for measure in z_measures
            rejection_list, metrics_key = measure_to_list_and_key[measure]
            for channel in selected_channels
                for epoch_idx in metrics[metrics_key][channel]
                    rejection = Rejection(channel, epoch_idx)
                    push!(rejected_info, rejection)
                    push!(rejection_list, rejection)
                end
            end
        end
    end

    # Handle absolute threshold
    abs_rejections = abs_criterion > 0 ? Rejection[] : nothing
    if abs_criterion > 0
        for channel in selected_channels
            for epoch_idx in metrics[:absolute_threshold][channel]
                rejection = Rejection(channel, epoch_idx)
                push!(rejected_info, rejection)
                push!(abs_rejections, rejection)
            end
        end
    end

    # Create z-score rejection info if z_criterion > 0
    z_rejections = z_criterion > 0 ? ZScoreRejectionInfo(z_measures, z_variance, z_max, z_min, z_abs, z_range, z_kurtosis) : nothing

    # Create rejection info
    info = EpochInfo(dat.condition, dat.condition_name, length(dat.data))

    rejection_info = EpochRejectionInfo(
        name,
        info,
        length(rejected_info),
        abs_criterion,
        abs_rejections,
        z_criterion,
        z_rejections,
        unique_rejections(rejected_info),
        nothing,  # repaired - populated during repair
        nothing,  # skipped - populated during repair
    )

    return rejection_info
end

detect_bad_epochs_automatic(dat::Vector{EpochData}; kwargs...) = detect_bad_epochs_automatic.(dat; kwargs...)


"""
    get_rejected(info::EpochRejectionInfo)::Vector{Rejection}

Get the list of rejected channel-epoch pairs from the rejection info.

# Examples
```julia
info = detect_bad_epochs_automatic(epochs)
rejected = get_rejected(info)
```
"""
get_rejected(info::EpochRejectionInfo)::Vector{Rejection} = info.rejected
get_rejected(info::Vector{EpochRejectionInfo})::Vector{Vector{Rejection}} = get_rejected.(info)



"""
Calculate statistical metrics for all epochs.

Returns a Dict with keys: `:z_variance`, `:z_max`, `:z_min`, `:z_abs`, `:z_range`, 
`:z_kurtosis`, `:absolute_threshold`. Each value is a Dict mapping channel symbols 
to vectors of epoch indices that exceeded the criteria.
"""
function _calculate_epoch_metrics(
    dat::EpochData,
    selected_channels::AbstractVector{Symbol},
    selected_samples::AbstractVector{<:Integer},
    z_criterion::Real,
    abs_criterion::Real,
)::Dict{Symbol,Dict{Symbol,Vector{Int}}}

    metric_keys = [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis, :absolute_threshold]
    metrics = Dict(k => Dict(ch => Int[] for ch in selected_channels) for k in metric_keys)

    Threads.@threads for ch in selected_channels
        # subset the data by the selected sample interval for the metric calculation
        channel_data_all = [epoch[selected_samples, ch]::Vector{Float64} for epoch in dat.data]

        if abs_criterion > 0
            abs_threshold_violations = findall(epoch_data -> maximum(abs, epoch_data) > abs_criterion, channel_data_all)
            append!(metrics[:absolute_threshold][ch], abs_threshold_violations)
        end

        if z_criterion > 0

            variances = var.(channel_data_all)
            max_values = maximum.(channel_data_all)
            min_values = minimum.(channel_data_all)
            abs_values = [maximum(abs, epoch_data) for epoch_data in channel_data_all]
            ranges = max_values .- min_values
            kurtoses = kurtosis.(channel_data_all)

            z_scores = zscore.([variances, max_values, min_values, abs_values, ranges, kurtoses])
            z_metric_keys = [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis]

            for (z_score, metric_key) in zip(z_scores, z_metric_keys)
                bad_epochs = findall(abs.(z_score) .> z_criterion)
                append!(metrics[metric_key][ch], bad_epochs)
            end

        end

    end

    return metrics
end


# === REPORTING FUNCTIONS ===

"""Pretty-print an `EpochRejectionInfo` summary with criterion, counts, and breakdowns."""
function Base.show(io::IO, info::EpochRejectionInfo)
    println(io, "EpochRejectionInfo: $(info.name)")
    println(io, "Condition: $(info.info.number): $(info.info.name)")
    println(io, "  Abs criterion: $(info.abs_criterion > 0 ? string(info.abs_criterion,  " μV") : "disabled")")
    println(io, "  Z-criterion: $(info.z_criterion > 0 ? string(info.z_criterion) : "disabled")")
    println(io, "  Epochs total: $(info.info.n), Epochs rejected: $(length(unique_epochs(info.rejected)))")
    println(io, "  Artifacts total: $(info.n_artifacts)")
    println(io, "  Rejected epochs: $(_print_vector(unique_epochs(info.rejected)))")

    if !isnothing(info.abs_rejections)
        println(io, "  Rejection breakdown (absolute):")
        println(
            io,
            "    Abs threshold: $(length(unique_epochs(info.abs_rejections))) unique epochs, $(length(unique_channels(info.abs_rejections))) unique channels",
        )
    end

    if !isnothing(info.z_rejections)
        z_info = info.z_rejections
        # Map measures to fields and labels
        field_map = Dict(
            :variance => (z_info.z_variance, "Z-Variance"),
            :max => (z_info.z_max, "Z-Maximum"),
            :min => (z_info.z_min, "Z-Minimum"),
            :abs => (z_info.z_abs, "Z-Absolute"),
            :range => (z_info.z_range, "Z-Range"),
            :kurtosis => (z_info.z_kurtosis, "Z-Kurtosis"),
        )

        # Determine which selected measures actually have entries
        nonempty_selected = Tuple{Vector{Rejection},String}[]
        for m in z_info.z_measures
            vec, label = field_map[m]
            if !isempty(vec)
                push!(nonempty_selected, (vec, label))
            end
        end
        if !isempty(nonempty_selected)
            println(io, "  Rejection breakdown (z-score):")
            for (vec, label) in nonempty_selected
                println(io, "    $(label):  $(length(unique_epochs(vec))) unique epochs, $(length(unique_channels(vec))) unique channels")
            end
        end
    end
    println(io, "")

end

"""Show each `EpochRejectionInfo` in a vector."""
Base.show(io::IO, infos::Vector{EpochRejectionInfo}) = Base.show.(Ref(io), infos)




"""
    repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Repair detected artifacts using the specified method.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use

# Available Methods
- `:neighbor_interpolation` - Weighted neighbor interpolation (default). Uses `dat.layout.neighbours` for neighbor information.
- `:spherical_spline` - Spherical spline interpolation

# Keyword Arguments (for :spherical_spline method)
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)

# Returns
Nothing (mutates dat in-place)
"""
function repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)

    if method ∉ [:neighbor_interpolation, :spherical_spline]
        throw(ArgumentError("Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"))
    end

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Repairing artifacts using method: $method"
    if method == :neighbor_interpolation
        # Determine which channels can be repaired if not already done
        if isnothing(artifacts.repaired)
            channel_repairable!(artifacts, dat.layout)
        end
        repair_artifacts_neighbor!(dat, artifacts; kwargs...)
    elseif method == :spherical_spline
        repair_artifacts_spherical_spline!(dat, artifacts; kwargs...)
    end
    return nothing
end


"""
    repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Non-mutating version of repair_artifacts!. Creates a copy of the data and repairs artifacts without modifying the original.

# Arguments
- `dat::EpochData`: The epoch data to repair (NOT modified)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use (default: :neighbor_interpolation)

# Keyword Arguments
Same as repair_artifacts! for the respective methods.

# Returns
- `EpochData`: A new EpochData object with repaired artifacts

# Examples
```julia
# Basic artifact repair (creates new object)
repaired_epochs = repair_artifacts(epochs, artifacts)

# Using spherical spline method
repaired_epochs = repair_artifacts(epochs, artifacts, :spherical_spline)

# Rejecting bad epochs entirely (use reject_epochs instead)
clean_epochs = reject_epochs(epochs, artifacts)
```
"""
function repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)
    dat_copy = copy(dat)
    repair_artifacts!(dat_copy, artifacts; method, kwargs...)
    return dat_copy
end

function repair_artifacts(dat::Vector{EpochData}, artifacts::Vector{EpochRejectionInfo}; kwargs...)
    return repair_artifacts.(dat, artifacts; kwargs...)
end

function repair_artifacts!(dat::Vector{EpochData}, artifacts::Vector{EpochRejectionInfo}; kwargs...)
    repair_artifacts!.(dat, artifacts; kwargs...)
    return nothing
end



"""
    channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)

Analyze which channels can be repaired and which cannot, based on neighbor availability.
Populates `artifacts.repaired` and `artifacts.skipped` with the analysis.

# Arguments
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts (mutated to add repair analysis)
- `layout::Layout`: Layout object containing neighbor information

# Returns
- `EpochRejectionInfo`: The same artifacts object (modified in-place)

# Notes
Only analyzes repairability - does not perform any repairs.
Use `repair_artifacts_neighbor!` to actually perform the repairs after this analysis.
"""
function channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)
    rejected = unique([r.epoch for r in artifacts.rejected])

    # Initialize tracking dictionaries in artifacts struct (OrderedDict to maintain sorted order)
    artifacts.repaired = OrderedDict{Int,Vector{Symbol}}()
    artifacts.skipped = OrderedDict{Int,Vector{Symbol}}()

    # Process epochs in sorted order to maintain ordering in OrderedDict
    for epoch_idx in sort(rejected)
        bad_channels = [artifact.channel for artifact in artifacts.rejected if artifact.epoch == epoch_idx]
        isempty(bad_channels) && continue

        repairable_channels = check_channel_neighbors(bad_channels, layout)

        if isempty(repairable_channels)
            if length(bad_channels) == 1
                @info "Epoch $epoch_idx: Cannot repair channel $(bad_channels[1]) (fewer than 2 neighbors)"
            else
                @info "Epoch $epoch_idx: Cannot repair channels $(bad_channels) (bad neighbors and/or fewer than 2 neighbors)"
            end
            # Track that all were skipped
            artifacts.skipped[epoch_idx] = bad_channels
        else
            skipped_channels = setdiff(bad_channels, repairable_channels)
            artifacts.repaired[epoch_idx] = repairable_channels
            if !isempty(skipped_channels)
                @info "Epoch $epoch_idx: Skipping repair of $(length(skipped_channels)) channel(s) with bad neighbors: $skipped_channels"
                artifacts.skipped[epoch_idx] = skipped_channels
            end
        end
    end

    return artifacts
end

channel_repairable!(artifacts::Vector{EpochRejectionInfo}, layout::Layout) = channel_repairable!.(artifacts, Ref(layout))


"""
    repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)

Repair artifacts using weighted neighbor interpolation.
Uses `artifacts.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Uses `dat.layout.neighbours` for neighbor information.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information with `repaired` already populated

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)
    # Check if repaired has been populated
    if isnothing(artifacts.repaired)
        throw(ArgumentError("repaired not populated. Call channel_repairable! first."))
    end

    if isempty(artifacts.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return dat
    end

    # Process epochs in sorted order (already sorted in OrderedDict)
    for (epoch_idx, repairable_channels) in artifacts.repaired
        @info "Epoch $epoch_idx: Repairing channels $(repairable_channels) using neighbor interpolation"

        # Use unified channel repair function with epoch selection
        repair_channels!(
            dat,
            repairable_channels;
            method = :neighbor_interpolation,
            epoch_selection = epochs([epoch_idx]),
            neighbours_dict = dat.layout.neighbours,
        )
        @info "" # formatting
    end

    return nothing
end

"""
    repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int=4, lambda::Real=1e-5)

Repair artifacts using spherical spline interpolation.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Real`: Regularization parameter (default: 1e-5)

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)
"""
function repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int = 4, lambda::Real = 1e-5)

    _ensure_coordinates_3d!(dat.layout)

    # Get all rejected epochs with their bad channels
    rejected = unique([r.epoch for r in artifacts.rejected])

    for epoch_idx in rejected
        # Get bad channels for this epoch
        bad_channels = [artifact.channel for artifact in artifacts.rejected if artifact.epoch == epoch_idx]
        isempty(bad_channels) && continue

        @info "Repairing epoch $epoch_idx channels $(bad_channels) using spherical spline interpolation"

        # Use unified channel repair function with epoch selection
        repair_channels!(dat, bad_channels; method = :spherical_spline, epoch_selection = epochs([epoch_idx]), m = m, lambda = lambda)
    end

    return dat
end



"""
    subset_bad_data(data_path::String, threshold::Float64; 
                    subset_directory::String = "excluded")

Identify and move files from participants with low data retention to a separate directory.

Participants are considered "bad" if they have less than the threshold percentage
of data remaining in ANY condition. All files associated with bad participants
are moved to a subdirectory (default: "excluded").

The function searches for "epoch_summary.jld2" in the specified directory.

# Arguments
- `data_path::String`: Path to the directory containing epoch_summary.jld2 and preprocessed files
- `threshold::Float64`: Minimum percentage threshold (e.g., 75.0 means 75% retention required)

# Keyword Arguments
- `subset_directory::String`: Name of subdirectory for excluded participants (default: "excluded")

# Examples
```julia
# Move participants with < 75% data in any condition to "excluded" subdirectory
subset_bad_data("preprocessed_files", 75.0)

# Specify custom subset directory name
subset_bad_data("/path/to/preprocessed", 80.0, subset_directory="excluded")
```
"""
function subset_bad_data(data_path::String, threshold::Real; subset_directory::String = "excluded")

    # Validate inputs/outputs
    (threshold < 0.0 || threshold > 100.0) && @minimal_error("threshold must be 0 < threshold < 100, got $threshold")
    !isdir(data_path) && @minimal_error("data_path must be a directory: $data_path")

    epoch_summary_path = joinpath(data_path, "epoch_summary.jld2")
    !isfile(epoch_summary_path) && @minimal_error("epoch_summary.jld2 not found: $data_path")

    # Load epoch summary (plain DataFrame, not EegFunData — can't use read_data)
    epoch_summary = load(epoch_summary_path, "data")

    # Use data_path as output directory
    output_directory = abspath(data_path)
    subset_dir_path = joinpath(output_directory, subset_directory)
    !isdir(subset_dir_path) && mkpath(subset_dir_path)

    # Find participants with any condition below threshold
    bad_participants = unique(epoch_summary.file[epoch_summary.percentage .< threshold])
    println("Subsetting data: $(length(bad_participants))")
    println("   N remaining: $(length(unique(epoch_summary.file)) - length(bad_participants))")
    println("   N removed: $(length(bad_participants))")

    # Create subset summary files with only non-excluded participants
    # Filter epoch_summary to exclude bad participants
    epoch_summary_subset = epoch_summary[.!in.(epoch_summary.file, Ref(bad_participants)), :]

    # Load and filter file_summary_subset
    file_summary_path = joinpath(data_path, "file_summary.jld2")
    !isfile(file_summary_path) && @minimal_error("file_summary.jld2 not found: $data_path")

    file_summary = load(file_summary_path, "data")
    file_summary_subset = file_summary[.!in.(file_summary.file, Ref(bad_participants)), :]

    # Save subset summary files
    epoch_summary_subset_path = joinpath(output_directory, "epoch_summary_subset.jld2")
    file_summary_subset_path = joinpath(output_directory, "file_summary_subset.jld2")

    jldsave(epoch_summary_subset_path; data = epoch_summary_subset)
    jldsave(file_summary_subset_path; data = file_summary_subset)

    # Find and move all files for bad participants
    all_files = readdir(output_directory)
    for participant_id in bad_participants
        # Find all files that contain participant_id
        matching_files = filter(all_files) do filename
            occursin(participant_id, filename)
        end
        # Move each matching file
        for filename in matching_files
            src_path = joinpath(output_directory, filename)
            dst_path = joinpath(subset_dir_path, filename)
            if isfile(src_path)
                mv(src_path, dst_path, force = true)
            end
        end
    end
end
