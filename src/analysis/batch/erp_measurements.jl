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
- `analysis_type::String`: Type of measurement ("mean_amp", "max_peak", "min_peak")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `baseline_interval::Union{Tuple{Real, Real}, Nothing}`: Baseline interval in seconds (e.g., (-0.2, 0) or (-0.2, 0.0))
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `output_file::String`: Base name for output files (default: "erp_measurements")

# Supported Analysis Types
- `"mean_amp"`: Mean amplitude in the analysis window
- `"max_peak"`: Maximum peak amplitude in the analysis window
- `"min_peak"`: Minimum peak amplitude in the analysis window
- `"max_peak_lat"`: Latency (time) of maximum peak in the analysis window
- `"min_peak_lat"`: Latency (time) of minimum peak in the analysis window

# Returns
- `DataFrame`: Results with columns [participant, condition, channel, value]

# Output Files
- CSV file with measurements for each participant/condition/channel
- Log file documenting the analysis parameters

# Examples
```julia
# Mean amplitude between 100-200 ms with baseline
erp_measurements("erps", (0.1, 0.2), "mean_amp", baseline_interval=(-0.2, 0.0))

# Maximum peak between 300-500 ms for specific participants
erp_measurements("erps", (0.3, 0.5), "max_peak", 
                participants=[1, 2, 3], 
                conditions=[1, 2])

# Minimum peak for specific channels
erp_measurements("erps", (0.0, 0.6), "min_peak",
                channel_selection=channels([:Fz, :Cz, :Pz]))
```

# References
Adapted from analyseMyERPs.m by Ian Greenhouse
Based on methods described in:
- Luck, S. J. (2014). An introduction to the event-related potential technique (2nd ed.). MIT Press.
"""
function erp_measurements(file_pattern::String, 
                         analysis_window::Tuple{Real, Real}, 
                         analysis_type::String;
                         input_dir::String = pwd(),
                         baseline_interval::Union{Tuple{Real, Real}, Nothing} = nothing,
                         participants::Union{Int, Vector{Int}, Nothing} = nothing,
                         conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                         channel_selection::Function = channels(),
                         output_dir::Union{String, Nothing} = nothing,
                         output_file::String = "erp_measurements")
    
    # Set up global logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)
    
    try
        @info "ERP measurement analysis started at $(now())"
        
        # Log the function call
        @log_call "erp_measurements" (file_pattern, analysis_window, analysis_type)
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        valid_types = ["mean_amp", "max_peak", "min_peak", "max_peak_lat", "min_peak_lat"]
        if analysis_type ∉ valid_types
            @minimal_error_throw("Analysis type must be one of: $(join(valid_types, ", ")). Got: $analysis_type")
        end
        
        if analysis_window[1] >= analysis_window[2]
            @minimal_error_throw("Analysis window start must be before end. Got: $analysis_window")
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            output_dir = joinpath(input_dir, "measurements_$(analysis_type)")
        end
        mkpath(output_dir)
        
        @info "Analysis type: $analysis_type"
        @info "Analysis window: $(analysis_window[1])s to $(analysis_window[2])s"
        if baseline_interval !== nothing
            @info "Baseline interval: $(baseline_interval[1])s to $(baseline_interval[2])s"
        end
        
        # Find JLD2 files matching the pattern
        all_files = readdir(input_dir)
        jld2_files = filter(x -> endswith(x, ".jld2") && contains(x, file_pattern), all_files)
        
        # Filter by participant number if specified
        if participants !== nothing
            jld2_files = _filter_files(jld2_files; include = participants)
        end
        
        if isempty(jld2_files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(jld2_files)) JLD2 files"
        
        # Initialize results storage (will be converted to wide format at the end)
        results_long = []
        channel_order = Symbol[] # Store channel order from first file
        
        processed_count = 0
        error_count = 0
        
        for (file_idx, file) in enumerate(jld2_files)
            input_path = joinpath(input_dir, file)
            
            @info "Processing: $file ($(file_idx)/$(length(jld2_files)))"
            
            try
                # Extract participant ID from filename
                # Try to find a numeric participant ID in the filename
                parts = split(replace(file, ".jld2" => ""), "_")
                participant_str = findfirst(p -> !isempty(p) && all(isdigit, p), parts)
                
                if participant_str === nothing
                    # No numeric ID found, use filename as identifier
                    participant = hash(file) % 10000  # Create unique numeric ID from filename
                    @info "  No numeric participant ID found in filename, using hash: $participant"
                else
                    participant = parse(Int, parts[participant_str])
                end
                
                # Load the file
                file_data = load(input_path)
                
                # Determine data type
                data_var = nothing
                if haskey(file_data, "erps")
                    data_var = file_data["erps"]
                elseif haskey(file_data, "epochs")
                    data_var = file_data["epochs"]
                else
                    @minimal_warning "No recognized data variable in $file"
                    error_count += 1
                    continue
                end
                
                # Filter by condition if specified
                if conditions !== nothing
                    condition_nums = conditions isa Int ? [conditions] : conditions
                    data_var = data_var[condition_nums]
                end
                
                # Process each condition
                for (cond_idx, data) in enumerate(data_var)
                    condition = conditions !== nothing ? conditions[cond_idx] : cond_idx
                    
                    @info "  Processing condition $condition"
                    
                    # Apply baseline if specified
                    if baseline_interval !== nothing
                        interval = IntervalTime(baseline_interval[1], baseline_interval[2])
                        data = baseline(data, interval; channel_selection = channel_selection)
                    end
                    
                    # Handle both ErpData and EpochData types
                    is_epoch_data = data isa EpochData
                    
                    # Get list of dataframes to process
                    dfs_to_process = if is_epoch_data
                        data.data  # Vector of DataFrames, one per epoch
                    else
                        [data.data]  # Single DataFrame for ErpData
                    end
                    
                    # Get selected channels once per condition
                    selected_channels = get_selected_channels(data, channel_selection, include_meta = false, include_extra = false)
                    
                    if isempty(selected_channels)
                        @minimal_warning "No channels selected for analysis"
                        continue
                    end
                    
                    # Store channel order from first file/condition
                    if isempty(channel_order)
                        channel_order = selected_channels
                    end
                    
                    # Check for condition_name and epoch once (optimization #4)
                    first_df = dfs_to_process[1]
                    has_condition_name = hasproperty(first_df, :condition_name)
                    has_epoch = is_epoch_data && hasproperty(first_df, :epoch)
                    
                    # Get time vector and compute time indices once per condition (optimization #6)
                    time_col = first_df[!, :time]
                    time_idx = findall(t -> analysis_window[1] <= t <= analysis_window[2], time_col)
                    
                    if isempty(time_idx)
                        @minimal_warning "Analysis window outside data time range for $file condition $condition"
                        continue
                    end
                    
                    # Process each dataframe (epoch or single ERP)
                    for df in dfs_to_process
                        # Compute measurements for all channels
                        channel_values = Vector{Pair{Symbol, Float64}}(undef, length(selected_channels))
                        
                        for (i, chan_symbol) in enumerate(selected_channels)
                            # Use view to avoid copying (optimization #5)
                            chan_data = @view df[time_idx, chan_symbol]
                            
                            # Compute measurement based on analysis type
                            value = if analysis_type == "mean_amp"
                                mean(chan_data)
                            elseif analysis_type == "max_peak"
                                maximum(chan_data)
                            elseif analysis_type == "min_peak"
                                minimum(chan_data)
                            elseif analysis_type == "max_peak_lat"
                                max_idx = argmax(chan_data)
                                time_col[time_idx[max_idx]]
                            elseif analysis_type == "min_peak_lat"
                                min_idx = argmin(chan_data)
                                time_col[time_idx[min_idx]]
                            end
                            
                            channel_values[i] = chan_symbol => value
                        end
                        
                        # Build named tuple with metadata and channel values
                        if has_condition_name && has_epoch
                            row_data = (; participant, condition, 
                                         condition_name = df[1, :condition_name],
                                         epoch = df[1, :epoch],
                                         channel_values...)
                        elseif has_condition_name
                            row_data = (; participant, condition,
                                         condition_name = df[1, :condition_name],
                                         channel_values...)
                        elseif has_epoch
                            row_data = (; participant, condition,
                                         epoch = df[1, :epoch],
                                         channel_values...)
                        else
                            row_data = (; participant, condition, channel_values...)
                        end
                        
                        # Add this row to results
                        push!(results_long, row_data)
                    end
                end
                
                processed_count += 1
                
            catch e
                @error "Error processing $file: $e"
                @error "Stack trace: $(stacktrace(catch_backtrace()))"
                error_count += 1
                continue
            end
        end
        
        # Convert to wide-format DataFrame
        if !isempty(results_long)
            results = DataFrame(results_long)
            
            # Reorder columns: participant, condition, (condition_name if present), (epoch if present), then channels
            col_list = [:participant, :condition]
            if hasproperty(results, :condition_name)
                push!(col_list, :condition_name)
            end
            if hasproperty(results, :epoch)
                push!(col_list, :epoch)
            end
            append!(col_list, channel_order)
            
            select!(results, col_list...)
            
            # Save results to CSV
            csv_path = joinpath(output_dir, "$(output_file).csv")
            CSV.write(csv_path, results)
            @info "Results saved to: $csv_path"
            @info "Total rows: $(nrow(results))"
            @info "Columns: participant, condition, + $(length(channel_order)) channels"
        else
            @minimal_warning "No results to save"
            results = DataFrame()
        end
        
        @info "Analysis complete! Processed $processed_count files successfully, $error_count errors"
        
        return results
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "$(output_file).log"
        log_dest = joinpath(output_dir, "$(output_file).log")
        if log_source != log_dest && isfile(log_source)
            mv(log_source, log_dest, force=true)
        end
    end
end
