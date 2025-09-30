"""
Batch channel summary statistics for EEG/ERP data.
"""

#=============================================================================
    CHANNEL-SUMMARY-SPECIFIC HELPERS
=============================================================================#

"""Generate default output directory name for channel summary."""
function _default_channel_summary_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "channel_summary_$(pattern)")
end

#=============================================================================
    CHANNEL-SUMMARY-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through channel summary pipeline.
Returns tuple of (BatchResult, Vector{DataFrame}) with all condition results.
"""
function _process_channel_summary_file(filepath::String,
                                     conditions, sample_selection::Function,
                                     channel_selection::Function, include_extra::Bool)
    filename = basename(filepath)
    
    # Load data
    data_result = _load_eeg_data(filepath)
    if isnothing(data_result)
        return (BatchResult(false, filename, "No recognized data variable"), DataFrame[])
    end
    
    data_var, var_name = data_result
    
    # Select conditions
    data_var = _select_conditions(data_var, conditions)
    
    # Process each condition and collect results
    summary_dfs = DataFrame[]
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
        
        push!(summary_dfs, summary_df)
    end
    
    n_conditions = length(summary_dfs)
    return (BatchResult(true, filename, "Processed $n_conditions condition(s)"), summary_dfs)
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
        output_dir = something(output_dir, _default_channel_summary_output_dir(input_dir, file_pattern))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files to process"
        
        # Process all files and collect DataFrames
        all_summaries = DataFrame[]
        n_success = 0
        n_error = 0
        
        for (i, file) in enumerate(files)
            @info "Channel summary: $file ($i/$(length(files)))"
            
            input_path = joinpath(input_dir, file)
            
            result, summary_dfs = try
                _process_channel_summary_file(input_path, conditions, 
                                            sample_selection, channel_selection, 
                                            include_extra)
            catch e
                @error "Error processing $file" exception=(e, catch_backtrace())
                (BatchResult(false, file, "Exception: $(sprint(showerror, e))"), DataFrame[])
            end
            
            # Log result
            if result.success
                @info "  ✓ $(result.message)"
                append!(all_summaries, summary_dfs)
                n_success += 1
            else
                @minimal_warning "  ✗ $(result.message)"
                n_error += 1
            end
        end
        
        # Save combined results to single CSV
        if !isempty(all_summaries)
            combined_df = vcat(all_summaries...)
            output_path = joinpath(output_dir, "$(output_file).csv")
            CSV.write(output_path, combined_df)
            @info "Combined results saved to: $output_path"
        else
            @minimal_warning "No results to save"
        end
        
        @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
        @info "Output saved to: $output_dir"
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end