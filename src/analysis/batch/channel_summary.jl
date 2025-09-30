"""
Batch channel summary statistics for EEG/ERP data.
"""

#=============================================================================
    CHANNEL-SUMMARY-SPECIFIC HELPERS
=============================================================================#

"""Generate default output directory name for channel summary."""
function _default_channel_summary_output_dir(input_dir::String)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    joinpath(input_dir, "channel_summary_$(timestamp)")
end

#=============================================================================
    CHANNEL-SUMMARY-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through channel summary pipeline.
Returns BatchResult with success/failure info.
Saves individual CSV files for each condition.
"""
function _process_channel_summary_file(filepath::String, output_dir::String,
                                     conditions, sample_selection::Function,
                                     channel_selection::Function, include_extra::Bool,
                                     output_file::String)
    filename = basename(filepath)
    
    # Load data
    data_result = _load_eeg_data(filepath)
    if isnothing(data_result)
        return BatchResult(false, filename, "No recognized data variable")
    end
    
    data_var, var_name = data_result
    
    # Select conditions
    data_var = _select_conditions(data_var, conditions)
    
    # Process each condition
    n_conditions = 0
    for (cond_idx, data) in enumerate(data_var)
        condition = isnothing(conditions) ? cond_idx : 
                   (conditions isa Int ? conditions : conditions[cond_idx])
        
        # Compute channel summary
        summary_df = eegfun.channel_summary(data; 
                                           sample_selection = sample_selection,
                                           channel_selection = channel_selection,
                                           include_extra = include_extra)
        
        # Add metadata columns
        insertcols!(summary_df, 1, :file => splitext(filename)[1])
        insertcols!(summary_df, 2, :condition => condition)
        
        # Save to CSV
        output_filename = "$(splitext(filename)[1])_condition$(condition)_$(output_file).csv"
        output_path = joinpath(output_dir, output_filename)
        CSV.write(output_path, summary_df)
        
        n_conditions += 1
    end
    
    return BatchResult(true, filename, "Processed $n_conditions condition(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    channel_summary(file_pattern::String; 
                    input_dir::String = pwd(), 
                    participants::Union{Int, Vector{Int}, Nothing} = nothing,
                    conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                    sample_selection::Function = samples(),
                    channel_selection::Function = channels(),
                    include_extra::Bool = false,
                    output_dir::Union{String, Nothing} = nothing,
                    output_file::String = "channel_summary")

Batch process EEG/ERP data files to compute channel summary statistics.

This function loads JLD2 files, computes channel summary statistics for each file,
and saves the results to CSV files.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "epochs", "erps", "cleaned")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `sample_selection::Function`: Function for sample filtering (default: all samples)
- `channel_selection::Function`: Function for channel filtering (default: all channels)
- `include_extra::Bool`: Whether to include extra channels (default: false)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `output_file::String`: Base name for output files (default: "channel_summary")

# Examples
```julia
# Compute channel summary for all epoch files
channel_summary("epochs")

# Process specific participants and conditions
channel_summary("erps_cleaned", participants=[1, 2, 3], conditions=[1, 2])

# Compute summary for specific channels only
channel_summary("epochs", channel_selection=channels([:Fp1, :Fp2, :F3, :F4]))

# Include extra channels (EOG, etc.)
channel_summary("epochs", include_extra=true)
```
"""
function channel_summary(file_pattern::String; 
                        input_dir::String = pwd(), 
                        participants::Union{Int, Vector{Int}, Nothing} = nothing,
                        conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                        sample_selection::Function = samples(),
                        channel_selection::Function = channels(),
                        include_extra::Bool = false,
                        output_dir::Union{String, Nothing} = nothing,
                        output_file::String = "channel_summary")
    
    # Setup logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)
    
    try
        @info "Batch channel summary started at $(now())"
        @log_call "channel_summary" (file_pattern,)
        
        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        # Setup directories
        output_dir = something(output_dir, _default_channel_summary_output_dir(input_dir))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files to process"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> 
            _process_channel_summary_file(input_path, output_dir, conditions, 
                                        sample_selection, channel_selection, 
                                        include_extra, output_file)
        
        # Note: output_path not actually used, we save per-condition files
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; 
                                      operation_name="Channel summary")
        
        # Log summary
        _log_batch_summary(results, output_dir)
        
    finally
        # Cleanup logging
        _cleanup_logging(log_file, output_dir)
    end
end