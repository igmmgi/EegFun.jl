"""
Batch ERP measurements (amplitude, latency) for EEG/ERP data.
"""

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC VALIDATION
=============================================================================#

"""Validate analysis type, returning error message or nothing."""
function _validate_analysis_type(analysis_type::String)
    valid_types = ["mean_amp", "max_peak", "min_peak", "max_peak_lat", "min_peak_lat"]
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

    if participant_str === nothing
        # No numeric ID found, use hash of filename
        participant = hash(filename) % 10000
        @info "  No numeric participant ID found in filename, using hash: $participant"
        return participant
    else
        return parse(Int, parts[participant_str])
    end
end

"""
Compute measurement for a single channel in a time window.
"""
function _compute_measurement(
    chan_data::AbstractVector,
    time_col::AbstractVector,
    time_idx::AbstractVector,
    analysis_type::String,
)
    if analysis_type == "mean_amp"
        return mean(chan_data)
    elseif analysis_type == "max_peak"
        return maximum(chan_data)
    elseif analysis_type == "min_peak"
        return minimum(chan_data)
    elseif analysis_type == "max_peak_lat"
        max_idx = argmax(chan_data)
        return time_col[time_idx[max_idx]]
    elseif analysis_type == "min_peak_lat"
        min_idx = argmin(chan_data)
        return time_col[time_idx[min_idx]]
    end
end

"""
Process a single DataFrame (epoch or ERP) to extract measurements.
Returns named tuple with metadata and channel measurements.
"""
function _process_dataframe_measurements(
    df::DataFrame,
    selected_channels::Vector{Symbol},
    analysis_window::Function,
    analysis_type::String,
    participant::Int,
)
    # Find time column and apply analysis window predicate
    time_col = df[!, :time]
    sample_mask = analysis_window(df)
    time_idx = findall(sample_mask)

    if isempty(time_idx)
        @minimal_warning "No time points found matching analysis window"
        return nothing
    end

    # Compute measurements for all channels
    channel_values = Vector{Pair{Symbol,Float64}}(undef, length(selected_channels))

    for (i, chan_symbol) in enumerate(selected_channels)
        # Use view to avoid copying (optimization #5)
        chan_data = @view df[time_idx, chan_symbol]

        # Compute measurement
        value = _compute_measurement(chan_data, time_col, time_idx, analysis_type)
        channel_values[i] = chan_symbol => value
    end

    # Get metadata from DataFrame if available
    condition = hasproperty(df, :condition) ? df[1, :condition] : nothing
    condition_name = hasproperty(df, :condition_name) ? df[1, :condition_name] : nothing
    epoch = hasproperty(df, :epoch) ? df[1, :epoch] : nothing

    # Build named tuple with metadata and channel values
    metadata = (participant=participant,)
    
    if condition !== nothing
        metadata = merge(metadata, (condition=condition,))
    end
    if condition_name !== nothing
        metadata = merge(metadata, (condition_name=condition_name,))
    end
    if epoch !== nothing
        metadata = merge(metadata, (epoch=epoch,))
    end

    return merge(metadata, NamedTuple(channel_values))
end

#=============================================================================
    ERP-MEASUREMENTS-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through ERP measurements pipeline.
Returns Vector of named tuples (one per condition/epoch).
"""
function _process_measurements_file(
    filepath::String,
    analysis_window::Function,
    analysis_type::String,
    baseline_window::Function,
    condition_selection::Function,
    channel_selection::Function,
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

    # Process each condition
    results = []

    for (cond_idx, data) in enumerate(data_var)
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

        # Get selected channels from first df (reuse from baseline section above)
        all_channels = setdiff(propertynames(first_df), metadata_cols)
        # Channel selection predicate operates on the entire vector
        selected_channels = all_channels[channel_selection(all_channels)]

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
            )

            if !isnothing(row_data)
                push!(results, row_data)
            end
        end
    end

    return results
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
- `analysis_type::String`: Type of measurement ("mean_amp", "max_peak", "min_peak", "max_peak_lat", "min_peak_lat")
- `analysis_window::Function`: Analysis window sample selection predicate (default: samples() - all samples)
- `baseline_window::Function`: Baseline window sample selection predicate (default: samples() - all samples, baseline skipped)
- `participant_selection::Function`: Participant selection predicate (default: participants() - all)
- `condition_selection::Function`: Condition selection predicate (default: conditions() - all)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `output_file::String`: Base name for output files (default: "erp_measurements")

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
)

    # Setup logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)

    try
        @info "ERP measurement analysis started at $(now())"
        @log_call "erp_measurements" (file_pattern, analysis_type)

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_analysis_type(analysis_type)) !== nothing
            @minimal_error_throw(error_msg)
        end

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

        # Process all files and collect results
        all_results = []
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
                )

                if !isnothing(file_results) && !isempty(file_results)
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

        # Convert to DataFrame
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
