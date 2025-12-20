"""
Batch ERP measurements (amplitude, latency) for EEG/ERP data.
"""

#=============================================================================
    DEFAULT KEYWORD ARGUMENTS
=============================================================================#
const ERP_MEASUREMENTS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Robust peak detection
    :local_window => (
        3,
        "Number of samples on each side of peak for robust peak detection (total window = 2*local_window + 1). Peak must be larger than neighbors and local averages within this window.",
    ),

    # Fractional area latency
    :fractional_area_fraction => (
        0.5,
        "Fraction for fractional area latency (0.0-1.0). Finds latency where this fraction of area is to the left.",
    ),

    # Fractional peak latency
    :fractional_peak_fraction => (
        0.5,
        "Fraction for fractional peak latency (0.0-1.0). Finds latency where amplitude is this fraction of peak.",
    ),
    :fractional_peak_direction =>
        (:onset, "Direction for fractional peak latency: :onset (before peak) or :offset (after peak)"),
)

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC VALIDATION
=============================================================================#

"""Validate analysis type, returning error message or nothing."""
function _validate_analysis_type(analysis_type::String)
    valid_types = [
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
    analysis_type ∉ valid_types &&
        return "Analysis type must be one of: $(join(valid_types, ", ")). Got: $analysis_type"
    return nothing
end


"""Generate default output directory name for ERP measurements."""
function _default_measurements_output_dir(input_dir::String, analysis_type::String)
    joinpath(input_dir, "measurements_$(analysis_type)")
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
    local_window::Int = 3,  # Number of samples on each side (ERPLAB "Neighborhood" parameter)
)
    n = length(data)
    total_window = 2 * local_window + 1

    if total_window > n
        @minimal_warning "local_window ($local_window) requires $(total_window) samples but data has only $n, cannot detect robust peak"
        return (nothing, nothing)
    elseif local_window > n ÷ 2
        @minimal_warning "local_window ($local_window) is > 50% of data length ($n), may be too large for robust peak detection"
    end

    valid_peaks = Tuple{Float64,Int}[]  # Store (value, index) for all valid peaks

    for i = (local_window+1):(n-local_window)
        val = data[i]

        # Use views to avoid allocations
        left_view = @view data[(i-local_window):(i-1)]
        right_view = @view data[(i+1):(i+local_window)]

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
    if isempty(data) || isempty(time_col)
        return NaN
    end

    if length(data) == 1
        return time_col[1]
    end

    # Compute total area using rectangular integration (fine for uniform sampling)
    dt = mean(diff(time_col))
    total_area = sum(data) * dt

    # Handle zero or very small area
    if abs(total_area) < 1e-10
        return time_col[1]  # Return start if no area
    end

    target_area = total_area * fraction

    # Handle fraction = 0.0 or 1.0
    if fraction <= 0.0
        return time_col[1]
    elseif fraction >= 1.0
        return time_col[end]
    end

    # Use binary search (like ERPLAB) for more accurate fractional area latency
    n = length(time_col)
    plow = 1
    phigh = n

    while plow <= phigh
        pmid = round(Int, (plow + phigh) / 2)
        # Compute cumulative area up to pmid
        cumulative_area = sum(data[1:pmid]) * dt

        if cumulative_area > target_area
            phigh = pmid - 1
        elseif cumulative_area < target_area
            plow = pmid + 1
        else
            break  # Exact match found
        end
    end

    # Return the latency at the found index
    return time_col[min(plow, n)]
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
    if abs(peak_val) < 1e-10
        return time_col[peak_idx]
    end

    target_val = peak_val * fraction

    if direction == :onset
        # Work backward from peak
        for i = (peak_idx-1):-1:1
            if data[i] <= target_val
                return time_col[i]
            end
        end
        return time_col[1]
    else  # :offset
        # Work forward from peak
        for i = (peak_idx+1):length(data)
            if data[i] <= target_val
                return time_col[i]
            end
        end
        return time_col[end]
    end
end

"""
Helper function for peak measurements with robust detection and fallback.
Returns (peak_value, peak_index) or (nothing, nothing) if no robust peak found.
"""
function _compute_peak_measurement(
    chan_data::AbstractVector,
    peak_type::Symbol,
    local_window::Int,
    channel_name::Symbol,
    measurement_name::String,
)
    n_samples = length(chan_data)
    total_window = 2 * local_window + 1
    if total_window > n_samples
        @minimal_warning "Channel $channel_name: local_window ($local_window, total window $total_window) > analysis window length ($n_samples samples). Cannot detect robust peak, using simple $measurement_name."
    end

    peak_val, peak_idx = _find_robust_peak(chan_data, peak_type; local_window = local_window)
    if isnothing(peak_val)
        @minimal_warning "Channel $channel_name: No robust $(peak_type == :max ? "maximum" : "minimum") peak found, using simple $measurement_name"
        peak_idx = peak_type == :max ? argmax(chan_data) : argmin(chan_data)
        peak_val = peak_type == :max ? maximum(chan_data) : minimum(chan_data)
    end

    return (peak_val, peak_idx)
end

"""
Compute measurement for a single channel in a time window.
"""
function _compute_measurement(
    chan_data::AbstractVector,
    time_col::AbstractVector,
    time_idx::AbstractVector,
    selected_times::Union{AbstractVector,Nothing},
    analysis_type::String,
    measurement_kwargs::Dict{Symbol,Any},
    channel_name::Symbol,
)
    if analysis_type == "mean_amplitude"
        return mean(chan_data)

        # Peak measurements (use robust detection with fallback to simple)
    elseif analysis_type in ["max_peak_amplitude", "min_peak_amplitude", "max_peak_latency", "min_peak_latency"]
        peak_type = startswith(analysis_type, "max") ? :max : :min
        local_window = measurement_kwargs[:local_window]
        measurement_name =
            analysis_type == "max_peak_amplitude" ? "maximum" :
            analysis_type == "min_peak_amplitude" ? "minimum" :
            analysis_type == "max_peak_latency" ? "maximum latency" : "minimum latency"

        peak_val, peak_idx =
            _compute_peak_measurement(chan_data, peak_type, local_window, channel_name, measurement_name)

        if analysis_type in ["max_peak_amplitude", "min_peak_amplitude"]
            return peak_val
        else  # latency measurements
            return time_col[time_idx[peak_idx]]
        end

        # Peak-to-peak measurements
    elseif analysis_type in ["peak_to_peak_amplitude", "peak_to_peak_latency"]
        local_window = measurement_kwargs[:local_window]

        # Find both max and min peaks
        max_val, max_idx = _compute_peak_measurement(chan_data, :max, local_window, channel_name, "maximum")
        min_val, min_idx = _compute_peak_measurement(chan_data, :min, local_window, channel_name, "minimum")

        # Handle cases where peaks weren't found
        if isnothing(max_val) || isnothing(min_val)
            @minimal_warning "Channel $channel_name: Could not find both maximum and minimum peaks for peak-to-peak measurement"
            return NaN
        end

        if analysis_type == "peak_to_peak_amplitude"
            return max_val - min_val
        else  # peak_to_peak_latency
            max_time = time_col[time_idx[max_idx]]
            min_time = time_col[time_idx[min_idx]]
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
            return sum(abs.(chan_data)) * dt
        elseif analysis_type == "integral"
            return sum(chan_data) * dt
        elseif analysis_type == "positive_area"
            return sum(max.(chan_data, 0.0)) * dt
        else  # negative_area
            return sum(abs.(min.(chan_data, 0.0))) * dt
        end

        # Fractional latency measurements
    elseif analysis_type == "fractional_area_latency"
        isnothing(selected_times) && error("selected_times should not be nothing for fractional_area_latency")
        fraction = measurement_kwargs[:fractional_area_fraction]
        return _fractional_area_latency(chan_data, selected_times, fraction)
    elseif analysis_type == "fractional_peak_latency"
        isnothing(selected_times) && error("selected_times should not be nothing for fractional_peak_latency")
        # Find the peak with maximum absolute value using robust detection
        local_window = measurement_kwargs[:local_window]
        max_val, max_idx = _compute_peak_measurement(chan_data, :max, local_window, channel_name, "maximum")
        min_val, min_idx = _compute_peak_measurement(chan_data, :min, local_window, channel_name, "minimum")

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
    analysis_window::Function,
    analysis_type::String,
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

    # Find time column and apply analysis window predicate
    time_col = df[!, :time]
    sample_mask = analysis_window(df)
    time_idx = findall(sample_mask)

    if isempty(time_idx)
        @minimal_warning "No time points found matching analysis window"
        return nothing
    end

    # Pre-compute selected times once (used by area/fractional measurements)
    selected_times = nothing
    if analysis_type in [
        "rectified_area",
        "integral",
        "positive_area",
        "negative_area",
        "fractional_area_latency",
        "fractional_peak_latency",
    ]
        selected_times = @view time_col[time_idx]
    end

    metadata_pairs = Pair{Symbol,Any}[:participant=>participant]
    if hasproperty(df, :condition)
        push!(metadata_pairs, :condition => df[1, :condition])
    end
    if hasproperty(df, :condition_name)
        push!(metadata_pairs, :condition_name => df[1, :condition_name])
    end
    if hasproperty(df, :epoch)
        push!(metadata_pairs, :epoch => df[1, :epoch])
    end

    # Compute measurements for all channels
    channel_pairs = Vector{Pair{Symbol,Float64}}(undef, length(selected_channels))
    for (i, chan_symbol) in enumerate(selected_channels)
        chan_data = @view df[time_idx, chan_symbol]
        value = _compute_measurement(
            chan_data,
            time_col,
            time_idx,
            selected_times,
            analysis_type,
            measurement_kwargs,
            chan_symbol,
        )
        channel_pairs[i] = chan_symbol => value
    end

    # Combine into single NamedTuple
    return merge(NamedTuple(metadata_pairs), NamedTuple(channel_pairs))
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
    analysis_window::Function,
    analysis_type::String,
    baseline_window::Function,
    condition_selection::Function,
    channel_selection::Function,
    measurement_kwargs::Dict{Symbol,Any},
)
    filename = basename(filepath)

    # Extract participant ID
    participant = _extract_participant_id(filename)

    # Load data (using load_data which finds by type)
    data_var = load_data(filepath)
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


    # Get selected channels from first data
    first_data = data_var[1]
    metadata_cols = meta_labels(first_data)
    first_df = if first_data isa DataFrame
        first_data
    elseif hasproperty(first_data, :data)
        first_data.data isa Vector{DataFrame} ? first_data.data[1] : first_data.data
    else
        first_data.data
    end
    all_channels = setdiff(propertynames(first_df), metadata_cols)
    channel_mask = channel_selection(all_channels)
    selected_channels = all_channels[channel_mask]

    if isempty(selected_channels)
        @minimal_warning "No channels selected in $filename"
        return nothing
    end

    # Collect results as NamedTuples
    results = Vector{Any}()

    for data in data_var
        # Get metadata columns from the original data object
        metadata_cols = meta_labels(data)

        # Get DataFrame(s) - ErpData has single df, EpochData has vector
        dfs_to_process = if data isa DataFrame
            [data]
        elseif hasproperty(data, :data)
            data.data isa Vector{DataFrame} ? data.data : [data.data]
        else
            [data.data]
        end

        # Apply baseline correction if baseline window is specified (not default all samples)
        first_df = dfs_to_process[1]
        baseline_mask = baseline_window(first_df)
        baseline_indices = findall(baseline_mask)

        # Skip baseline if window matches all samples (default) or no samples
        if length(baseline_indices) < nrow(first_df) && !isempty(baseline_indices)
            # Get channel columns (exclude metadata)
            all_channels = setdiff(propertynames(first_df), metadata_cols)
            # channels() returns a predicate that operates on a vector, not individual elements
            eeg_channels = all_channels[channels()(all_channels)]

            if !isempty(eeg_channels)
                # Use existing baseline infrastructure from baseline.jl
                interval = IntervalIndex(start = first(baseline_indices), stop = last(baseline_indices))
                @info "Applying baseline correction to $(length(eeg_channels)) channels over interval: $(first(baseline_indices)) to $(last(baseline_indices))"
                _apply_baseline!(dfs_to_process, eeg_channels, interval)
            else
                @minimal_warning "No EEG channels found for baseline correction"
            end
        elseif isempty(baseline_indices)
            @minimal_warning "Baseline window matched no samples. Skipping baseline correction."
        end

        # Use pre-determined selected_channels (already computed for type stability)

        # Add condition and condition_name to DataFrames if not present and available from data object
        if (data isa ErpData || data isa EpochData)
            for df in dfs_to_process
                if !hasproperty(df, :condition)
                    insertcols!(df, 1, :condition => data.condition)
                end
                if !hasproperty(df, :condition_name)
                    insertcols!(df, hasproperty(df, :condition) ? 2 : 1, :condition_name => data.condition_name)
                end
            end
        end

        # Process each dataframe (epoch or single ERP)
        for df in dfs_to_process
            row_data = _process_dataframe_measurements(
                df,
                selected_channels,
                analysis_window,
                analysis_type,
                participant,
                measurement_kwargs,
            )

            if !isnothing(row_data)
                push!(results, row_data)
            end
        end
    end

    # Return vector of NamedTuples (will be converted to DataFrame at top level)
    return isempty(results) ? nothing : results
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    erp_measurements(file_pattern::String, analysis_type::String;
                    analysis_window::Function = samples(),
                    baseline_window::Function = samples(),
                    participant_selection::Function = participants(),
                    condition_selection::Function = conditions(),
                    channel_selection::Function = channels(),
                    input_dir::String = pwd(),
                    output_dir::Union{String, Nothing} = nothing,
                    output_file::String = "erp_measurements")

Perform standard ERP measurements on averaged or epoched EEG data.

This function computes basic ERP measurements (mean amplitude, peak amplitude, peak latency) 
across specified time windows and saves results to CSV files.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps", "epochs_cleaned")
- `analysis_type::String`: Type of measurement:
  - **Amplitude**: "mean_amplitude", "max_peak_amplitude", "min_peak_amplitude", "peak_to_peak_amplitude"
  - **Latency**: "max_peak_latency", "min_peak_latency", "peak_to_peak_latency", "fractional_area_latency", "fractional_peak_latency"
  - **Area/Integral** (in µVs): 
    - "rectified_area": Sum of absolute values (always positive)
    - "integral": Signed integral (net area, can be positive or negative)
    - "positive_area": Area of positive values only
    - "negative_area": Area of negative values only (as absolute value)
  
  Note: Peak measurements use robust detection (local peak with neighbor/average checks) and fall back to simple peak if no robust peak is found.
- `analysis_window::Function`: Analysis window sample selection predicate (default: samples() - all samples)
- `baseline_window::Function`: Baseline window sample selection predicate (default: samples() - all samples, baseline skipped)
- `participant_selection::Function`: Participant selection predicate (default: participants() - all)
- `condition_selection::Function`: Condition selection predicate (default: conditions() - all)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `output_file::String`: Base name for output files (default: "erp_measurements")
- `kwargs...`: Additional keyword arguments for measurement parameters:
  - `local_window::Int`: Window size (in samples) for robust peak detection (default: 3)
  - `fractional_area_fraction::Float64`: Fraction for fractional area latency (default: 0.5)
  - `fractional_peak_fraction::Float64`: Fraction for fractional peak latency (default: 0.5)
  - `fractional_peak_direction::Symbol`: Direction for fractional peak latency - :onset (before peak) or :offset (after peak) (default: :onset)

# Examples
```julia
# Mean amplitude between 100-200 ms with baseline
erp_measurements("erps", "mean_amplitude", analysis_window=samples((0.1, 0.2)), baseline_window=samples((-0.2, 0.0)))

# Maximum peak between 300-500 ms for specific participants
erp_measurements("erps", "max_peak_amplitude", analysis_window=samples((0.3, 0.5)), participant_selection=participants([1, 2, 3]), condition_selection=conditions([1, 2]))

# Exclude specific participants
erp_measurements("erps", "max_peak_amplitude", analysis_window=samples((0.3, 0.5)), participant_selection=participants_not([10, 11]))

# Minimum peak for specific channels
erp_measurements("erps", "min_peak_amplitude", analysis_window=samples((0.0, 0.6)), channel_selection=channels([:Fz, :Cz, :Pz]))
```
"""
function erp_measurements(
    file_pattern::String,
    analysis_type::String;
    analysis_window::Function = samples(),
    baseline_window::Function = samples(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
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
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw error_msg
        end
        if (error_msg = _validate_analysis_type(analysis_type)) !== nothing
            @minimal_error_throw error_msg
        end

        # Merge measurement kwargs with defaults
        measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)

        # Validate measurement kwargs
        local_window = measurement_kwargs[:local_window]
        if local_window < 1
            @minimal_error "local_window must be >= 1, got: $local_window"
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

        for (file_idx, file) in enumerate(files)
            input_path = joinpath(input_dir, file)
            @info "Processing: $file ($file_idx/$(length(files)))"

            try
                file_results = _process_measurements_file(
                    input_path,
                    analysis_window,
                    analysis_type,
                    baseline_window,
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
                # @minimal_error_throw throws ErrorException
                if e isa ErrorException
                    error_msg = e.msg
                    if occursin("invalid", lowercase(error_msg)) && occursin("window", lowercase(error_msg))
                        rethrow(e)
                    end
                end
                @error "Error processing $file" exception=(e, catch_backtrace())
                error_count += 1
            end
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

        # Evaluate predicates on actual data to get readable descriptions
        # Load the first file to get actual time range
        analysis_window_desc = "custom"
        baseline_window_desc = baseline_window === nothing ? "none" : "custom"

        if !isempty(files)
            try
                println("Loading data from $(files[1])")
                first_file = joinpath(input_dir, files[1])
                data = load_data(first_file)

                mask = analysis_window(data[1].data)



                selected_times = time(data)[mask]
                time_min = round(minimum(selected_times), digits = 3)
                time_max = round(maximum(selected_times), digits = 3)
                analysis_window_desc = "$time_min:$time_max S"
                mask = baseline_window(data[1].data)
                selected_times = time(data)[mask]
                time_min = round(minimum(selected_times), digits = 3)
                time_max = round(maximum(selected_times), digits = 3)
                baseline_window_desc = "$time_min:$time_max S"

            catch
                println("Error loading data: ")
                # If loading fails, use "custom" descriptions
            end
        else
            @minimal_warning "No files found"
        end

        # Return as ErpMeasurementsResult with metadata
        return ErpMeasurementsResult(
            results_df,
            analysis_type,
            analysis_window,
            analysis_window_desc,
            baseline_window,
            baseline_window_desc,
        )

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
