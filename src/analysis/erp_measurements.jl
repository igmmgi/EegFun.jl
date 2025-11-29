"""
Batch ERP measurements (amplitude, latency) for EEG/ERP data.
"""

#=============================================================================
    DEFAULT KEYWORD ARGUMENTS
=============================================================================#
const ERP_MEASUREMENTS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Robust peak detection
    :local_window => (3, "Window size (in samples) for robust peak detection. Peak must be larger than neighbors and local averages within this window."),
    
    # Fractional area latency
    :fractional_area_fraction => (0.5, "Fraction for fractional area latency (0.0-1.0). Finds latency where this fraction of area is to the left."),
    
    # Fractional peak latency
    :fractional_peak_fraction => (0.5, "Fraction for fractional peak latency (0.0-1.0). Finds latency where amplitude is this fraction of peak."),
    :fractional_peak_direction => (:onset, "Direction for fractional peak latency: :onset (before peak) or :offset (after peak)"),
)

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC VALIDATION
=============================================================================#

"""Validate analysis type, returning error message or nothing."""
function _validate_analysis_type(analysis_type::String)
    valid_types = [
        "mean_amp", "max_peak", "min_peak", "max_peak_lat", "min_peak_lat",
        "rectified_area", "integral", "positive_area", "negative_area",
        "fractional_area_lat", "fractional_peak_lat"
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

"""Extract participant ID from filename, returns Int."""
function _extract_participant_id(filename::String)
    parts = split(replace(filename, ".jld2" => ""), "_")
    participant_str = findfirst(p -> !isempty(p) && all(isdigit, p), parts)

    if participant_str === nothing # No numeric ID found, use hash of filename
        participant = hash(filename) % 10000
        @info "  No numeric participant ID found in filename, using hash: $participant"
        return participant
    else # Numeric ID found, use it
        return parse(Int, parts[participant_str])
    end
end

"""
Find robust peak (peak that is larger than neighbors and local averages).
Returns (peak_value, peak_index) or (nothing, nothing) if no robust peak found.
"""
function _find_robust_peak(
    data::AbstractVector,
    peak_type::Symbol;  # :max or :min
    local_window::Int = 3
)
    n = length(data)
    if local_window >= n
        @minimal_warning "local_window ($local_window) is >= data length ($n), cannot detect robust peak"
        return (nothing, nothing)
    elseif local_window > n ÷ 2
        @minimal_warning "local_window ($local_window) is > 50% of data length ($n), may be too large for robust peak detection"
    end
    
    half_window = div(local_window, 2)
    valid_peaks = Tuple{Float64,Int}[]  # Store (value, index) for all valid peaks
    
    for i in (half_window + 1):(n - half_window)
        val = data[i]
        
        # Use views to avoid allocations
        left_view = @view data[(i - half_window):(i - 1)]
        right_view = @view data[(i + 1):(i + half_window)]
        
        # Check if this is a local peak (unified logic for max/min)
        is_peak = if peak_type == :max
            val > data[i-1] && val > data[i+1] &&
            val > mean(left_view) && val > mean(right_view)
        else  # :min
            val < data[i-1] && val < data[i+1] &&
            val < mean(left_view) && val < mean(right_view)
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
function _fractional_area_latency(
    data::AbstractVector,
    time_col::AbstractVector,
    fraction::Float64
)
    # Compute cumulative area (integral)
    dt = mean(diff(time_col))
    cumulative = cumsum(data) .* dt
    total_area = cumulative[end]
    
    target_area = total_area * fraction
    
    # Find first point where cumulative area >= target
    idx = findfirst(x -> x >= target_area, cumulative)
    isnothing(idx) && return time_col[end]
    
    return time_col[idx]
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
    direction::Symbol  # :onset (before peak) or :offset (after peak)
)
    peak_val = data[peak_idx]
    target_val = peak_val * fraction
    
    if direction == :onset
        # Work backward from peak
        for i in (peak_idx - 1):-1:1
            if data[i] <= target_val
                return time_col[i]
            end
        end
        return time_col[1]
    else  # :offset
        # Work forward from peak
        for i in (peak_idx + 1):length(data)
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
    if local_window >= n_samples
        @minimal_warning "Channel $channel_name: local_window ($local_window) >= analysis window length ($n_samples samples). Cannot detect robust peak, using simple $measurement_name."
    end
    
    peak_val, peak_idx = _find_robust_peak(chan_data, peak_type; local_window=local_window)
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
    if analysis_type == "mean_amp"
        return mean(chan_data)
        
    # Peak measurements (use robust detection with fallback to simple)
    elseif analysis_type in ["max_peak", "min_peak", "max_peak_lat", "min_peak_lat"]
        peak_type = startswith(analysis_type, "max") ? :max : :min
        local_window = measurement_kwargs[:local_window]
        measurement_name = analysis_type == "max_peak" ? "maximum" :
                          analysis_type == "min_peak" ? "minimum" :
                          analysis_type == "max_peak_lat" ? "maximum latency" : "minimum latency"
        
        peak_val, peak_idx = _compute_peak_measurement(chan_data, peak_type, local_window, channel_name, measurement_name)
        
        if analysis_type in ["max_peak", "min_peak"]
            return peak_val
        else  # latency measurements
            return time_col[time_idx[peak_idx]]
        end
        
    # Area/Integral measurements (in µVs) - compute dt only when needed
    elseif analysis_type in ["rectified_area", "integral", "positive_area", "negative_area"]
        dt = length(selected_times) > 1 ? mean(diff(selected_times)) : 0.0
        
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
    elseif analysis_type == "fractional_area_lat"
        isnothing(selected_times) && error("selected_times should not be nothing for fractional_area_lat")
        fraction = measurement_kwargs[:fractional_area_fraction]
        return _fractional_area_latency(chan_data, selected_times, fraction)
    elseif analysis_type == "fractional_peak_lat"
        isnothing(selected_times) && error("selected_times should not be nothing for fractional_peak_lat")
        max_abs_idx = argmax(abs.(chan_data))
        fraction = measurement_kwargs[:fractional_peak_fraction]
        direction = measurement_kwargs[:fractional_peak_direction]
        return _fractional_peak_latency(chan_data, selected_times, max_abs_idx, fraction, direction)
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
    if analysis_type in ["rectified_area", "integral", "positive_area", "negative_area", "fractional_area_lat", "fractional_peak_lat"]
        selected_times = @view time_col[time_idx]
    end

    # Build row as NamedTuple (much faster than DataFrame)
    # Start with metadata
    metadata_pairs = Pair{Symbol,Any}[:participant => participant]
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
        value = _compute_measurement(chan_data, time_col, time_idx, selected_times, analysis_type, measurement_kwargs, chan_symbol)
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

    # Select conditions
    data_var = _condition_select(data_var, condition_selection)
    
    if isempty(data_var)
        return nothing
    end

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
  - **Amplitude**: "mean_amp", "max_peak", "min_peak"
  - **Latency**: "max_peak_lat", "min_peak_lat", "fractional_area_lat", "fractional_peak_lat"
  - **Area/Integral** (in µVs): "rectified_area", "integral", "positive_area", "negative_area"
  
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
erp_measurements("erps", "mean_amp", analysis_window=samples((0.1, 0.2)), baseline_window=samples((-0.2, 0.0)))

# Maximum peak between 300-500 ms for specific participants
erp_measurements("erps", "max_peak", analysis_window=samples((0.3, 0.5)), participant_selection=participants([1, 2, 3]), condition_selection=conditions([1, 2]))

# Exclude specific participants
erp_measurements("erps", "max_peak", analysis_window=samples((0.3, 0.5)), participant_selection=participants_not([10, 11]))

# Minimum peak for specific channels
erp_measurements("erps", "min_peak", analysis_window=samples((0.0, 0.6)), channel_selection=channels([:Fz, :Cz, :Pz]))
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
    kwargs...
)

    # Setup logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)

    try

        @info "ERP measurement analysis started at $(now())"
        @log_call "erp_measurements" (file_pattern, analysis_type)

        # Validation 
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        if (error_msg = _validate_analysis_type(analysis_type)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Merge measurement kwargs with defaults
        measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)

        # Setup directories
        output_dir = something(output_dir, _default_measurements_output_dir(input_dir, analysis_type))
        mkpath(output_dir)

        # Log analysis parameters
        @info "Analysis type: $analysis_type"

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants = participant_selection)
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

        return results_df

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
