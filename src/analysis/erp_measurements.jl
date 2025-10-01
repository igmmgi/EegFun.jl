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

"""Validate analysis window, returning error message or nothing."""
function _validate_analysis_window(window::Tuple{Real,Real})
    window[1] >= window[2] && return "Analysis window start must be before end. Got: $window"
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
    analysis_window::Tuple{Real,Real},
    analysis_type::String,
    participant::Int,
    condition::Int,
    has_condition_name::Bool,
    has_epoch::Bool,
)
    # Find time column and analysis window indices (optimization #4)
    time_col = df[!, :time]
    time_idx = findall(t -> analysis_window[1] <= t <= analysis_window[2], time_col)

    if isempty(time_idx)
        @minimal_warning "No time points found in analysis window $(analysis_window) for condition $condition"
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

    # Build named tuple with metadata and channel values
    if has_condition_name && has_epoch
        return (;
            participant,
            condition,
            condition_name = df[1, :condition_name],
            epoch = df[1, :epoch],
            channel_values...,
        )
    elseif has_condition_name
        return (; participant, condition, condition_name = df[1, :condition_name], channel_values...)
    elseif has_epoch
        return (; participant, condition, epoch = df[1, :epoch], channel_values...)
    else
        return (; participant, condition, channel_values...)
    end
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
    analysis_window::Tuple{Real,Real},
    analysis_type::String,
    baseline_interval::Union{Tuple{Real,Real},Nothing},
    conditions,
    channel_selection::Function,
)
    filename = basename(filepath)

    # Extract participant ID
    participant = _extract_participant_id(filename)

    # Load data
    file_data = load(filepath)

    # Determine data type
    data_var = nothing
    if haskey(file_data, "erps")
        data_var = file_data["erps"]
    elseif haskey(file_data, "epochs")
        data_var = file_data["epochs"]
    else
        @minimal_warning "No recognized data variable in $filename"
        return nothing
    end

    # Select conditions
    data_var = _select_conditions(data_var, conditions)

    # Process each condition
    results = []

    for (cond_idx, data) in enumerate(data_var)
        condition = cond_idx

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

        # Apply baseline correction if specified
        if baseline_interval !== nothing
            # Convert time interval to index interval
            time_col = dfs_to_process[1][!, :time]
            baseline_start_idx = findfirst(x -> x >= baseline_interval[1], time_col)
            baseline_end_idx = findlast(x -> x <= baseline_interval[2], time_col)

            if isnothing(baseline_start_idx) || isnothing(baseline_end_idx)
                @minimal_warning "Baseline interval $(baseline_interval) is outside data range. Skipping baseline correction."
            else
                # Get channel columns (exclude metadata)
                all_channels = setdiff(propertynames(dfs_to_process[1]), metadata_cols)
                # channels() returns a predicate that operates on a vector, not individual elements
                eeg_channels = all_channels[channels()(all_channels)]

                if !isempty(eeg_channels)
                    # Use existing baseline infrastructure from baseline.jl
                    interval = IntervalIdx(baseline_start_idx, baseline_end_idx)
                    @info "Applying baseline correction to $(length(eeg_channels)) channels over interval: $(baseline_start_idx) to $(baseline_end_idx)"
                    _apply_baseline!(dfs_to_process, eeg_channels, interval)
                else
                    @minimal_warning "No EEG channels found for baseline correction"
                end
            end
        end

        # Get selected channels from first df
        first_df = dfs_to_process[1]
        all_channels = setdiff(propertynames(first_df), metadata_cols)
        # Channel selection predicate operates on the entire vector
        selected_channels = all_channels[channel_selection(all_channels)]

        # Check for condition_name and epoch columns once per condition (optimization #2)
        has_condition_name = hasproperty(first_df, :condition_name)
        has_epoch = hasproperty(first_df, :epoch)

        # Process each dataframe (epoch or single ERP)
        for df in dfs_to_process
            row_data = _process_dataframe_measurements(
                df,
                selected_channels,
                analysis_window,
                analysis_type,
                participant,
                condition,
                has_condition_name,
                has_epoch,
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
    erp_measurements(file_pattern::String, analysis_window::Tuple{Real, Real}, analysis_type::String;
                    input_dir::String = pwd(),
                    baseline_interval::Union{Tuple{Real, Real}, Nothing} = nothing,
                    participants::Union{Int, Vector{Int}, Nothing} = nothing,
                    conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                    channel_selection::Function = channels(),
                    output_dir::Union{String, Nothing} = nothing,
                    output_file::String = "erp_measurements")

Perform standard ERP measurements on averaged or epoched EEG data.

This function computes basic ERP measurements (mean amplitude, peak amplitude, peak latency) 
across specified time windows and saves results to CSV files.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps", "epochs_cleaned")
- `analysis_window::Tuple{Real, Real}`: Time window for analysis in seconds (e.g., (0.1, 0.2))
- `analysis_type::String`: Type of measurement ("mean_amp", "max_peak", "min_peak", "max_peak_lat", "min_peak_lat")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `baseline_interval::Union{Tuple{Real, Real}, Nothing}`: Baseline interval in seconds (e.g., (-0.2, 0.0))
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `output_file::String`: Base name for output files (default: "erp_measurements")

# Examples
```julia
# Mean amplitude between 100-200 ms with baseline
erp_measurements("erps", (0.1, 0.2), "mean_amp", baseline_interval=(-0.2, 0.0))

# Maximum peak between 300-500 ms for specific participants
erp_measurements("erps", (0.3, 0.5), "max_peak", 
                participants=[1, 2, 3], conditions=[1, 2])

# Minimum peak for specific channels
erp_measurements("erps", (0.0, 0.6), "min_peak",
                channel_selection=channels([:Fz, :Cz, :Pz]))
```
"""
function erp_measurements(
    file_pattern::String,
    analysis_window::Tuple{Real,Real},
    analysis_type::String;
    input_dir::String = pwd(),
    baseline_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    conditions::Union{Int,Vector{Int},Nothing} = nothing,
    channel_selection::Function = channels(),
    output_dir::Union{String,Nothing} = nothing,
    output_file::String = "erp_measurements",
)

    # Setup logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)

    try
        @info "ERP measurement analysis started at $(now())"
        @log_call "erp_measurements" (file_pattern, analysis_window, analysis_type)

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_analysis_type(analysis_type)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_analysis_window(analysis_window)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_measurements_output_dir(input_dir, analysis_type))
        mkpath(output_dir)

        # Log analysis parameters
        @info "Analysis type: $analysis_type"
        @info "Analysis window: $(analysis_window[1])s to $(analysis_window[2])s"
        if baseline_interval !== nothing
            @info "Baseline interval: $(baseline_interval[1])s to $(baseline_interval[2])s"
        end

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)

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
                    baseline_interval,
                    conditions,
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

        # Reorder columns: participant, condition, optional metadata, then channels in original order
        metadata_cols = [:participant, :condition]

        # Add optional metadata columns if present
        if hasproperty(results_df, :condition_name)
            push!(metadata_cols, :condition_name)
        end
        if hasproperty(results_df, :epoch)
            push!(metadata_cols, :epoch)
        end

        # Get channel columns (everything else)
        channel_cols = [col for col in names(results_df) if Symbol(col) ∉ metadata_cols]

        # Reorder
        select!(results_df, metadata_cols..., channel_cols...)

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
