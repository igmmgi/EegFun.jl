"""
Batch ERP measurements (amplitude, latency) for EEG/ERP data.
"""

#=============================================================================
    MEASUREMENT TYPES (single source of truth)
=============================================================================#

"""Valid ERP measurement type strings."""
const VALID_MEASUREMENT_TYPES = [
    "mean_amplitude",
    "max_peak_amplitude",
    "min_peak_amplitude",
    "max_peak_latency",
    "min_peak_latency",
    "peak_to_peak_amplitude",
    "peak_to_peak_latency",
    "rectified_area",
    "integral",
    "positive_area",
    "negative_area",
    "fractional_area_latency",
    "fractional_peak_latency",
]

"""Human-readable labels paired with measurement type strings, for use in GUI menus."""
const MEASUREMENT_TYPE_LABELS = [
    "Mean Amplitude" => "mean_amplitude",
    "Max Peak Amplitude" => "max_peak_amplitude",
    "Min Peak Amplitude" => "min_peak_amplitude",
    "Max Peak Latency" => "max_peak_latency",
    "Min Peak Latency" => "min_peak_latency",
    "Peak-to-Peak Amplitude" => "peak_to_peak_amplitude",
    "Peak-to-Peak Latency" => "peak_to_peak_latency",
    "Rectified Area" => "rectified_area",
    "Integral" => "integral",
    "Positive Area" => "positive_area",
    "Negative Area" => "negative_area",
    "Fractional Area Latency" => "fractional_area_latency",
    "Fractional Peak Latency" => "fractional_peak_latency",
]

"""Set of measurement types that return latency values (seconds) rather than amplitudes (μV)."""
const _LATENCY_MEASUREMENT_TYPES =
    Set(["max_peak_latency", "min_peak_latency", "peak_to_peak_latency", "fractional_area_latency", "fractional_peak_latency"])

"""Check whether a measurement type returns a latency value."""
_is_latency_measurement(t::String) = t in _LATENCY_MEASUREMENT_TYPES

#=============================================================================
    DEFAULT KEYWORD ARGUMENTS
=============================================================================#
const ERP_MEASUREMENTS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Robust peak detection
    :local_interval => (
        3,
        "Number of samples on each side of peak (total interval = 2*local_interval + 1). Peak must be larger than neighbors and local averages within this window.",
    ),

    # Fractional area latency
    :fractional_area_fraction =>
        (0.5, "Fraction for fractional area latency (0.0-1.0). Finds latency where this fraction of area is to the left."),

    # Fractional peak latency
    :fractional_peak_fraction =>
        (0.5, "Fraction for fractional peak latency (0.0-1.0). Finds latency where amplitude is this fraction of peak."),
    :fractional_peak_direction => (:onset, "Direction for fractional peak latency: :onset (before peak) or :offset (after peak)"),
)

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC VALIDATION
=============================================================================#

"""Validate analysis type, returning error message or nothing."""
function _validate_analysis_type(analysis_type::String)
    if analysis_type ∉ VALID_MEASUREMENT_TYPES
        return "Analysis type must be one of: $(join(VALID_MEASUREMENT_TYPES, ", ")). Got: $analysis_type"
    end
    return nothing
end


"""Generate default output directory name for ERP measurements."""
function _default_measurements_output_dir(input_dir::String, analysis_type::String)
    joinpath(input_dir, "measurements_$(analysis_type)")
end

"""Normalize interval to tuple format (handles AbstractRange and Tuple)."""
function _normalize_interval(interval::Interval)
    isnothing(interval) && return nothing
    interval isa AbstractRange ? (first(interval), last(interval)) : interval
end

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC HELPERS
=============================================================================#

"""
Find robust peak (peak that is larger than neighbors and local averages).
Returns (peak_value, peak_index) or (nothing, nothing) if no robust peak found.
"""
function _find_robust_peak(
    data::AbstractVector,
    peak_type::Symbol;  # :max or :min
    local_interval::Int = 3,  # Number of samples on each side (ERPLAB "Neighborhood" parameter)
)
    n = length(data)
    total_window = 2 * local_interval + 1

    if total_window > n
        @minimal_warning "local_interval ($local_interval) requires $(total_window) samples but data has only $n, cannot detect robust peak"
        return (nothing, nothing)
    elseif local_interval > n ÷ 2
        @minimal_warning "local_interval ($local_interval) is > 50% of data length ($n), may be too large for robust peak detection"
    end

    valid_peaks = Tuple{Float64,Int}[]  # Store (value, index) for all valid peaks

    for i = (local_interval+1):(n-local_interval)
        val = data[i]

        # Use views to avoid allocations
        left_view = @view data[(i-local_interval):(i-1)]
        right_view = @view data[(i+1):(i+local_interval)]

        # Check if this is a local peak (unified logic for max/min)
        # ERPLAB checks: (i) peak > both adjacent points, (ii) peak > average of neighbors
        is_peak = if peak_type == :max
            val > data[i-1] && val > data[i+1] && val > mean(left_view) && val > mean(right_view)
        else  # :min
            val < data[i-1] && val < data[i+1] && val < mean(left_view) && val < mean(right_view)
        end

        if is_peak
            push!(valid_peaks, (val, i))
        end
    end

    # Return the most extreme peak (largest for max, smallest for min)
    if !isempty(valid_peaks)
        if peak_type == :max
            return argmax(p -> p[1], valid_peaks)  # Returns tuple with maximum value
        else
            return argmin(p -> p[1], valid_peaks)  # Returns tuple with minimum value
        end
    end

    return (nothing, nothing)
end

"""
Compute fractional area latency - finds the point that divides area into specified fraction.
Returns latency at which fraction of area is to the left.
"""
function _fractional_area_latency(data::AbstractVector, time_col::AbstractVector, fraction::Float64)

    # Edge cases
    (isempty(data) || isempty(time_col)) && return NaN
    length(data) == 1 && return time_col[1]

    # Compute total area using cumulative sum for O(N) efficiency
    # Rectangular integration: area = sum(data) * dt
    dt = mean(diff(time_col))
    cum_area = cumsum(data) .* dt
    total_area = cum_area[end]

    # Handle zero or very small area
    if abs(total_area) < 1e-12
        return time_col[1]  # Return start if no area
    end

    target_area = total_area * fraction

    # Handle fraction boundaries
    if fraction <= 0.0
        return time_col[1]
    elseif fraction >= 1.0
        return time_col[end]
    end

    # Find crossing point (handle both positive and negative total area)
    if total_area > 0
        idx = findfirst(a -> a >= target_area, cum_area)
    else
        idx = findfirst(a -> a <= target_area, cum_area)
    end

    if isnothing(idx)
        return time_col[end]
    elseif idx == 1
        return time_col[1]
    end

    # Linear interpolation for sub-sample precision
    # Latency (x) such that area(x) = target_area
    y0 = cum_area[idx-1]
    y1 = cum_area[idx]
    x0 = time_col[idx-1]
    x1 = time_col[idx]

    if abs(y1 - y0) < 1e-14
        return x0
    end

    return x0 + (target_area - y0) * (x1 - x0) / (y1 - y0)
end

"""
Compute fractional peak latency - finds point where amplitude is fraction of peak.
Returns latency before (onset) or after (offset) peak.
"""
function _fractional_peak_latency(
    data::AbstractVector,
    time_col::AbstractVector,
    peak_idx::Int,
    fraction::Float64,
    direction::Symbol,  # :onset (before peak) or :offset (after peak)
)
    # Edge cases
    if isempty(data) || isempty(time_col) || peak_idx < 1 || peak_idx > length(data)
        return NaN
    end

    peak_val = data[peak_idx]

    # Handle zero peak
    if abs(peak_val) < 1e-12
        return time_col[peak_idx]
    end

    target_val = peak_val * fraction

    if direction == :onset
        # Work backward from peak
        for i = (peak_idx-1):-1:1
            # Check for crossing (handle both positive and negative peaks)
            if (peak_val > 0 && data[i] <= target_val) || (peak_val < 0 && data[i] >= target_val)
                # Linear interpolation between i and i+1
                y0, y1 = data[i], data[i+1]
                x0, x1 = time_col[i], time_col[i+1]

                if abs(y1 - y0) < 1e-14
                    return x0
                end
                return x0 + (target_val - y0) * (x1 - x0) / (y1 - y0)
            end
        end
        return time_col[1]
    else  # :offset
        # Work forward from peak
        for i = (peak_idx+1):length(data)
            # Check for crossing
            if (peak_val > 0 && data[i] <= target_val) || (peak_val < 0 && data[i] >= target_val)
                # Linear interpolation between i-1 and i
                y0, y1 = data[i-1], data[i]
                x0, x1 = time_col[i-1], time_col[i]

                if abs(y1 - y0) < 1e-14
                    return x1
                end
                return x0 + (target_val - y0) * (x1 - x0) / (y1 - y0)
            end
        end
        return time_col[end]
    end
end

"""
Helper function for peak measurements with robust detection and fallback.
Returns (peak_value, peak_index) or (nothing, nothing) if no robust peak found.
"""
function _compute_peak_measurement(chan_data::AbstractVector, peak_type::Symbol, local_interval::Int, channel_name::Symbol)
    total_window = 2 * local_interval + 1
    if total_window > length(chan_data)
        peak_name = peak_type == :max ? "maximum" : "minimum"
        @minimal_warning "Channel $channel_name: local_interval ($local_interval, total interval $total_window) > analysis interval length ($(length(chan_data)) samples). Cannot detect robust peak, using simple $peak_name."
    end

    peak_val, peak_idx = _find_robust_peak(chan_data, peak_type; local_interval = local_interval)
    if isnothing(peak_val)
        peak_name = peak_type == :max ? "maximum" : "minimum"
        @minimal_warning "Channel $channel_name: No robust $peak_name peak found, using simple $peak_name"
        peak_idx = peak_type == :max ? argmax(chan_data) : argmin(chan_data)
        peak_val = peak_type == :max ? maximum(chan_data) : minimum(chan_data)
    end

    return (peak_val, peak_idx)
end

"""
Compute measurement for a single channel in a time interval.
"""
function _compute_measurement(
    chan_data::AbstractVector,
    selected_times::AbstractVector,
    analysis_type::String,
    measurement_kwargs::Dict{Symbol,Any},
    channel_name::Symbol,
)
    if analysis_type == "mean_amplitude"
        return mean(chan_data)

        # Peak measurements (use robust detection with fallback to simple)
    elseif analysis_type in ["max_peak_amplitude", "min_peak_amplitude", "max_peak_latency", "min_peak_latency"]
        peak_type = startswith(analysis_type, "max") ? :max : :min
        local_interval = measurement_kwargs[:local_interval]

        peak_val, peak_idx = _compute_peak_measurement(chan_data, peak_type, local_interval, channel_name)

        if analysis_type in ["max_peak_amplitude", "min_peak_amplitude"]
            return peak_val
        else  # latency measurements
            return selected_times[peak_idx]
        end

        # Peak-to-peak measurements
    elseif analysis_type in ["peak_to_peak_amplitude", "peak_to_peak_latency"]
        local_interval = measurement_kwargs[:local_interval]

        # Find both max and min peaks
        max_val, max_idx = _compute_peak_measurement(chan_data, :max, local_interval, channel_name)
        min_val, min_idx = _compute_peak_measurement(chan_data, :min, local_interval, channel_name)

        # Handle cases where peaks weren't found
        if isnothing(max_val) || isnothing(min_val)
            @minimal_warning "Channel $channel_name: Could not find both maximum and minimum peaks for peak-to-peak measurement"
            return NaN
        end

        if analysis_type == "peak_to_peak_amplitude"
            return max_val - min_val
        else  # peak_to_peak_latency
            max_time = selected_times[max_idx]
            min_time = selected_times[min_idx]
            return abs(max_time - min_time)
        end

        # Area/Integral measurements (in µVs) - compute dt only when needed
    elseif analysis_type in ["rectified_area", "integral", "positive_area", "negative_area"]
        # Edge case: single sample or empty
        if length(selected_times) <= 1 || length(chan_data) <= 1
            return 0.0
        end

        dt = mean(diff(selected_times))

        # Handle case where dt is zero or very small (shouldn't happen with valid time data)
        if dt <= 0.0
            @minimal_warning "Channel $channel_name: time step is <= 0, returning 0.0 for area measurement"
            return 0.0
        end

        # Use rectangular integration (sum * dt) - perfectly fine for uniformly sampled data
        # For uniform sampling (which EEG/ERP data always is), this is equivalent to trapezoidal
        if analysis_type == "rectified_area"
            return sum(abs, chan_data) * dt
        elseif analysis_type == "integral"
            return sum(chan_data) * dt
        elseif analysis_type == "positive_area"
            return sum(x -> max(x, 0.0), chan_data) * dt
        else  # negative_area
            return sum(x -> abs(min(x, 0.0)), chan_data) * dt
        end

        # Fractional latency measurements
    elseif analysis_type == "fractional_area_latency"
        fraction = measurement_kwargs[:fractional_area_fraction]
        return _fractional_area_latency(chan_data, selected_times, fraction)
    elseif analysis_type == "fractional_peak_latency"
        # Find the peak with maximum absolute value using robust detection
        local_interval = measurement_kwargs[:local_interval]
        max_val, max_idx = _compute_peak_measurement(chan_data, :max, local_interval, channel_name)
        min_val, min_idx = _compute_peak_measurement(chan_data, :min, local_interval, channel_name)

        # Use the peak with larger absolute value
        peak_idx = abs(max_val) >= abs(min_val) ? max_idx : min_idx
        fraction = measurement_kwargs[:fractional_peak_fraction]
        direction = measurement_kwargs[:fractional_peak_direction]
        return _fractional_peak_latency(chan_data, selected_times, peak_idx, fraction, direction)
    end

    return nothing
end

"""
Process a single DataFrame (epoch or ERP) to extract measurements.
Returns DataFrame row with metadata and channel measurements.
"""
function _process_dataframe_measurements(
    df::DataFrame,
    selected_channels::Vector{Symbol},
    analysis_interval::Interval,
    baseline_interval::Interval,
    analysis_type::String,
    file::String,
    participant::Int,
    measurement_kwargs::Dict{Symbol,Any},
)
    # Basic validation
    if isempty(df)
        @minimal_warning "DataFrame is empty"
        return nothing
    end

    if !hasproperty(df, :time)
        @minimal_error "DataFrame must have a :time column"
    end

    if isempty(selected_channels)
        @minimal_warning "No channels selected"
        return nothing
    end

    # Find time indices within analysis interval
    time_col = df[!, :time]

    if isnothing(analysis_interval)
        # Use all samples
        time_idx = 1:nrow(df)
    else
        # Normalize to tuple
        interval = _normalize_interval(analysis_interval)

        # Validate interval
        if interval[1] > interval[2]
            @minimal_error "Analysis interval start must be <= stop. Got: ($(interval[1]), $(interval[2]))"
        end

        # Find time range
        start_idx, stop_idx = _find_idx_start_end(time_col, interval[1], interval[2])
        time_idx = start_idx:stop_idx
    end

    if isempty(time_idx)
        @minimal_warning "No time points found matching analysis interval"
        return nothing
    end

    # Pre-compute selected times once (used by latency, area, and fractional measurements)
    # This allows us to avoid passing time_idx to _compute_measurement
    selected_times = @view time_col[time_idx]

    # Build metadata pairs in fixed order to ensure column order in DataFrame
    metadata_vals = Any[file, participant]
    metadata_keys = Symbol[:file, :participant]

    if hasproperty(df, :condition)
        push!(metadata_vals, df[1, :condition])
        push!(metadata_keys, :condition)
    end
    if hasproperty(df, :condition_name)
        push!(metadata_vals, df[1, :condition_name])
        push!(metadata_keys, :condition_name)
    end
    if hasproperty(df, :epoch)
        push!(metadata_vals, df[1, :epoch])
        push!(metadata_keys, :epoch)
    end

    # Add analysis type
    push!(metadata_vals, analysis_type)
    push!(metadata_keys, :analysis_type)

    # Add baseline interval columns
    if !isnothing(baseline_interval)
        baseline = _normalize_interval(baseline_interval)
        push!(metadata_vals, baseline[1])
        push!(metadata_keys, :baseline_interval_start)
        push!(metadata_vals, baseline[2])
        push!(metadata_keys, :baseline_interval_end)
    else
        push!(metadata_vals, missing)
        push!(metadata_keys, :baseline_interval_start)
        push!(metadata_vals, missing)
        push!(metadata_keys, :baseline_interval_end)
    end

    # Add analysis interval columns
    if !isnothing(analysis_interval)
        ai = _normalize_interval(analysis_interval)
        push!(metadata_vals, ai[1])
        push!(metadata_keys, :analysis_interval_start)
        push!(metadata_vals, ai[2])
        push!(metadata_keys, :analysis_interval_end)
    else
        # Use actual time range from data
        push!(metadata_vals, time_col[time_idx[1]])
        push!(metadata_keys, :analysis_interval_start)
        push!(metadata_vals, time_col[time_idx[end]])
        push!(metadata_keys, :analysis_interval_end)
    end


    metadata_final = NamedTuple{Tuple(metadata_keys)}(Tuple(metadata_vals))

    # Compute measurements for all channels
    channel_pairs = Vector{Pair{Symbol,Float64}}(undef, length(selected_channels))
    for (i, chan_symbol) in enumerate(selected_channels)
        chan_data = @view df[time_idx, chan_symbol]
        value = _compute_measurement(chan_data, selected_times, analysis_type, measurement_kwargs, chan_symbol)
        channel_pairs[i] = chan_symbol => value
    end

    # Combine into single NamedTuple
    return merge(metadata_final, NamedTuple(channel_pairs))
end

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through ERP measurements pipeline.
Returns Vector of NamedTuples (one per condition/epoch).
"""
function _process_measurements_file(
    filepath::String,
    analysis_interval::Interval,
    analysis_type::String,
    baseline_interval::Interval,
    average_channels::Bool,
    condition_selection::Function,
    channel_selection::Function,
    measurement_kwargs::Dict{Symbol,Any},
)
    filename = basename(filepath)

    # Extract participant ID
    participant = _extract_participant_id(filename)

    # Read data (using read_data which finds by type)
    data_var = read_data(filepath)
    if isnothing(data_var)
        @minimal_warning "No data variables found in $filename"
        return nothing
    end

    # Validate that data is Vector of ErpData or EpochData
    if !(data_var isa Vector{<:Union{ErpData,EpochData}})
        @minimal_warning "Invalid data type in $filename: expected Vector{ErpData} or Vector{EpochData}, got $(typeof(data_var))"
        return nothing
    end

    # Select conditions using get_selected_conditions (more general than _condition_select)
    selected_indices = get_selected_conditions(data_var, condition_selection)
    isempty(selected_indices) && @minimal_error "No conditions left following condition selection in $filename"
    data_var = data_var[selected_indices]

    # Dispatch to the Vector method which will handle each data object appropriately
    return erp_measurements!(
        data_var,
        analysis_type;
        analysis_interval = analysis_interval,
        baseline_interval = baseline_interval,
        channel_selection = channel_selection,
        average_channels = average_channels,
        participant = participant,
        measurement_kwargs...,
    )
end


#=============================================================================
    TYPE-SPECIFIC ERP MEASUREMENT METHODS
=============================================================================#

"""
Helper function to apply baseline correction to data.
Returns the DataFrames vector (modified in place if baseline applied).
"""
function _apply_baseline_correction!(dfs::Vector{DataFrame}, baseline_interval::Interval, all_channels::Vector{Symbol})
    # Skip if no baseline specified
    if isnothing(baseline_interval)
        return dfs
    end

    # Get EEG channels for baseline
    eeg_channels = all_channels[channels()(all_channels)]

    if !isempty(eeg_channels)
        # Normalize to tuple
        baseline_int = _normalize_interval(baseline_interval)

        # Convert time tuple to index tuple using first dataframe's time column
        # (assume all dataframes have the same time points for epochs/ERPs)
        time_col = dfs isa Vector ? dfs[1].time : dfs.time
        start_idx, stop_idx = _find_idx_start_end(time_col, baseline_int[1], baseline_int[2])
        baseline_idx = (start_idx, stop_idx)

        # Try to apply baseline - don't fail measurements if baseline fails
        try
            @info "Applying baseline correction to $(length(eeg_channels)) channels using interval: $baseline_int"
            _apply_baseline!(dfs, eeg_channels, baseline_idx)
        catch e
            @minimal_warning "Baseline correction failed: $(sprint(showerror, e)). Continuing without baseline."
        end
    else
        @minimal_warning "No EEG channels found for baseline correction"
    end

    return dfs
end

"""
Average channel measurement columns into a single `:mean_channels` column.
Preserves all metadata columns and replaces individual channel columns with their mean.
"""
function _average_channel_columns(row::NamedTuple, selected_channels::Vector{Symbol})
    # Separate metadata from channel data
    channel_values = [Float64(row[ch]) for ch in selected_channels]
    avg_value = mean(channel_values)

    # Build new NamedTuple with metadata + single averaged column
    metadata_keys = filter(k -> k ∉ selected_channels, keys(row))
    metadata = NamedTuple{Tuple(metadata_keys)}(Tuple(row[k] for k in metadata_keys))
    return merge(metadata, (; dv_channels = avg_value))
end

"""
    erp_measurements!(dat::ErpData, analysis_type::String; analysis_interval = times(), baseline_interval = times(), channel_selection = channels(), participant = 0, kwargs...)
    erp_measurements!(dat::EpochData, analysis_type::String; kwargs...)
    erp_measurements!(data::Vector{<:Union{ErpData,EpochData}}, analysis_type::String; kwargs...)
    erp_measurements(file_pattern::String, analysis_type::String; analysis_interval = times(), baseline_interval = times(),
    participant_selection = participants(), condition_selection = conditions(), channel_selection = channels(),
    input_dir = pwd(), output_dir = nothing, output_file = "erp_measurements", kwargs...)

Extract ERP measurements from EEG data. The mutating (`!`) forms accumulate results in Julia;
the `file_pattern` form batch-processes JLD2 files and saves a CSV.

# `analysis_type` Options
**Essential:**
- `"mean_amplitude"` — average amplitude in the analysis interval
- `"max_peak_amplitude"` / `"min_peak_amplitude"` — robust peak detection (fallback to simple max/min)
- `"max_peak_latency"` / `"min_peak_latency"` — timing of peak (in seconds)

**Specialized:**
- `"peak_to_peak_amplitude"` / `"peak_to_peak_latency"`
- `"integral"` / `"rectified_area"` / `"positive_area"` / `"negative_area"` — area in µV·s
- `"fractional_area_latency"` / `"fractional_peak_latency"` — onset/offset-style latencies

$(_generate_kwargs_doc(ERP_MEASUREMENTS_KWARGS))

# Examples
```julia
# In-memory, single ERP
erp_measurements!(dat, "mean_amplitude", analysis_interval = (0.1, 0.2), baseline_interval = (-0.2, 0.0))

# Batch across participants
erp_measurements("erps", "max_peak_amplitude",
    analysis_interval = (0.3, 0.5),
    baseline_interval = (-0.2, 0.0),
    participant_selection = participants(1:20),
    condition_selection = conditions([1, 2]))
```
"""
function erp_measurements!(
    dat::ErpData,
    analysis_type::String;
    analysis_interval::Interval = times(),
    baseline_interval::Interval = times(),
    channel_selection::Function = channels(),
    participant::Int = 0,
    average_channels::Bool = false,
    kwargs...,
)
    # Merge measurement kwargs with defaults
    measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)

    # Get selected channels
    metadata_cols = meta_labels(dat)
    all_channels = setdiff(propertynames(dat.data), metadata_cols)
    channel_mask = channel_selection(all_channels)
    selected_channels = all_channels[channel_mask]

    if isempty(selected_channels)
        @minimal_warning "No channels selected"
        return nothing
    end

    # Apply baseline correction
    dfs = [dat.data]  # Wrap in vector for baseline function
    _apply_baseline_correction!(dfs, baseline_interval, all_channels)

    # Add condition metadata if not present
    df = dfs[1]
    if !hasproperty(df, :condition)
        insertcols!(df, 1, :condition => dat.condition)
    end
    if !hasproperty(df, :condition_name)
        insertcols!(df, 2, :condition_name => dat.condition_name)
    end

    # Process the single ERP DataFrame
    result = _process_dataframe_measurements(
        df,
        selected_channels,
        analysis_interval,
        baseline_interval,
        analysis_type,
        dat.file,
        participant,
        measurement_kwargs,
    )

    # Average channels if requested
    if average_channels && !isnothing(result)
        result = _average_channel_columns(result, selected_channels)
    end

    return result
end


function erp_measurements!(
    dat::EpochData,
    analysis_type::String;
    analysis_interval::Interval = times(),
    baseline_interval::Interval = times(),
    channel_selection::Function = channels(),
    participant::Int = 0,
    average_channels::Bool = false,
    kwargs...,
)
    # Merge measurement kwargs with defaults
    measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)

    # Get selected channels
    metadata_cols = meta_labels(dat)
    first_df = dat.data[1]
    all_channels = setdiff(propertynames(first_df), metadata_cols)
    channel_mask = channel_selection(all_channels)
    selected_channels = all_channels[channel_mask]

    if isempty(selected_channels)
        @minimal_warning "No channels selected"
        return nothing
    end

    # Apply baseline correction to all epochs
    _apply_baseline_correction!(dat.data, baseline_interval, all_channels)

    # Add condition metadata if not present
    for df in dat.data
        if !hasproperty(df, :condition)
            insertcols!(df, 1, :condition => dat.condition)
        end
        if !hasproperty(df, :condition_name)
            insertcols!(df, 2, :condition_name => dat.condition_name)
        end
    end

    # Process each epoch DataFrame
    results = Vector{NamedTuple}()
    for df in dat.data
        row_data = _process_dataframe_measurements(
            df,
            selected_channels,
            analysis_interval,
            baseline_interval,
            analysis_type,
            dat.file,
            participant,
            measurement_kwargs,
        )

        if !isnothing(row_data)
            if average_channels
                row_data = _average_channel_columns(row_data, selected_channels)
            end
            push!(results, row_data)
        end
    end

    return isempty(results) ? nothing : results
end


function erp_measurements!(data::Vector{<:Union{ErpData,EpochData}}, analysis_type::String; average_channels::Bool = false, kwargs...)
    all_results = Vector{Any}()

    for dat in data
        result = erp_measurements!(dat, analysis_type; average_channels = average_channels, kwargs...)

        if !isnothing(result)
            # ErpData returns single NamedTuple, EpochData returns Vector{NamedTuple}
            if result isa Vector
                append!(all_results, result)
            else
                push!(all_results, result)
            end
        end
    end

    return isempty(all_results) ? nothing : all_results
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#


function erp_measurements(
    file_pattern::String,
    analysis_type::String;
    analysis_interval::Interval = times(),
    baseline_interval::Interval = times(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    average_channels::Bool = false,
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    output_file::String = "erp_measurements",
    kwargs...,
)

    # Setup logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)

    try

        @info "ERP measurement analysis started at $(now())"
        @log_call "erp_measurements"

        # Validation 
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error error_msg
        end
        error_msg = _validate_analysis_type(analysis_type)
        if !isnothing(error_msg)
            @minimal_error error_msg
        end

        # Validate analysis_interval
        if !isnothing(analysis_interval)
            interval = _normalize_interval(analysis_interval)
            if interval[1] > interval[2]
                @minimal_error "Analysis interval start must be <= stop. Got: ($(interval[1]), $(interval[2]))"
            end
        end

        # Merge measurement kwargs with defaults
        measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)

        # Validate measurement kwargs
        local_interval = measurement_kwargs[:local_interval]
        if local_interval < 1
            @minimal_error "local_interval must be >= 1, got: $local_interval"
        end

        fractional_area_fraction = measurement_kwargs[:fractional_area_fraction]
        if fractional_area_fraction < 0.0 || fractional_area_fraction > 1.0
            @minimal_error "fractional_area_fraction must be in [0.0, 1.0], got: $fractional_area_fraction"
        end

        fractional_peak_fraction = measurement_kwargs[:fractional_peak_fraction]
        if fractional_peak_fraction < 0.0 || fractional_peak_fraction > 1.0
            @minimal_error "fractional_peak_fraction must be in [0.0, 1.0], got: $fractional_peak_fraction"
        end

        fractional_peak_direction = measurement_kwargs[:fractional_peak_direction]
        if fractional_peak_direction !== :onset && fractional_peak_direction !== :offset
            @minimal_error "fractional_peak_direction must be :onset or :offset, got: $fractional_peak_direction"
        end

        # Setup directories
        output_dir = something(output_dir, _default_measurements_output_dir(input_dir, analysis_type))
        mkpath(output_dir)

        # Log analysis parameters
        @info "Analysis type: $analysis_type"

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files"

        # Process all files and collect results as NamedTuples
        all_results = Vector{Any}()
        processed_count = 0
        error_count = 0

        p = Progress(length(files); desc = "Measuring ERPs...")
        for (file_idx, file) in enumerate(files)
            input_path = joinpath(input_dir, file)
            @info "Processing: $file ($file_idx/$(length(files)))"

            try
                file_results = _process_measurements_file(
                    input_path,
                    analysis_interval,
                    analysis_type,
                    baseline_interval,
                    average_channels,
                    condition_selection,
                    channel_selection,
                    measurement_kwargs,
                )

                if !isnothing(file_results) && !isempty(file_results)
                    # file_results is a vector of NamedTuples
                    append!(all_results, file_results)
                    processed_count += 1
                    @info "  ✓ Extracted $(length(file_results)) measurement(s)"
                else
                    error_count += 1
                    @minimal_warning "  ✗ No measurements extracted from $file"
                end
            catch e
                # Re-throw invalid window errors (user errors that should propagate)
                # @minimal_error throws ErrorException
                if e isa ErrorException
                    error_msg = e.msg
                    if occursin("invalid", lowercase(error_msg)) && occursin("window", lowercase(error_msg))
                        rethrow(e)
                    end
                end
                @minimal_warning "Error processing $file: $(sprint(showerror, e))"
                error_count += 1
            end
            next!(p)
        end

        # Check if we have results
        if isempty(all_results)
            @minimal_warning "No results to save"
            @info "Analysis complete! Processed $processed_count files successfully, $error_count errors"
            return nothing
        end

        # Convert to DataFrame once at the end
        results_df = DataFrame(all_results)

        # Sort by participant, then condition (and epoch if present)
        sort_cols = [:participant, :condition]
        if hasproperty(results_df, :epoch)
            push!(sort_cols, :epoch)
        end
        sort!(results_df, sort_cols)

        # Save results
        output_csv = joinpath(output_dir, "$(output_file).csv")
        CSV.write(output_csv, results_df)

        @info "Saved $(nrow(results_df)) measurement(s) to: $output_csv"
        @info "Analysis complete! Processed $processed_count files successfully, $error_count errors"

        # Return as ErpMeasurementsResult with metadata
        return ErpMeasurementsResult(results_df, analysis_type, analysis_interval, baseline_interval)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
