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

# Returns
- `Nothing`: Saves channel summary statistics to CSV files in output directory

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
    
    # Set up global logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)
    
    try
        @info "Batch channel summary started at $(now())"
        
        # Log the function call
        @log_call "channel_summary" (file_pattern,)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
            output_dir = joinpath(input_dir, "channel_summary_$(timestamp)")
        end
        
        # Create output directory if it doesn't exist
        mkpath(output_dir)
        
        # Find JLD2 files matching the pattern
        all_files = readdir(input_dir)
        jld2_files = Base.filter(x -> endswith(x, ".jld2") && contains(x, file_pattern), all_files)
        
        # Filter by participant number if specified
        if participants !== nothing
            jld2_files = _filter_files(jld2_files; include = participants)
        end
        
        if isempty(jld2_files)
            @minimal_warning "No JLD2 files found matching pattern \"$file_pattern\" in $input_dir"
            return nothing
        end
        
        @info "Found $(length(jld2_files)) JLD2 files to process"
        
        processed_count = 0
        error_count = 0
        
        for (file_idx, file) in enumerate(jld2_files)
            input_path = joinpath(input_dir, file)
            
            @info "Processing: $file ($(file_idx)/$(length(jld2_files)))"
            
            try
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
                    
                    # Compute channel summary
                    summary_df = eegfun.channel_summary(data; 
                                                       sample_selection = sample_selection,
                                                       channel_selection = channel_selection,
                                                       include_extra = include_extra)
                    
                    # Add metadata columns
                    insertcols!(summary_df, 1, :file => splitext(file)[1])
                    insertcols!(summary_df, 2, :condition => condition)
                    
                    # Save to CSV
                    output_filename = "$(splitext(file)[1])_condition$(condition)_$(output_file).csv"
                    output_path = joinpath(output_dir, output_filename)
                    CSV.write(output_path, summary_df)
                    
                    @info "  Saved: $output_filename"
                end
                
                processed_count += 1
                
            catch e
                @error "Error processing $file: $e"
                @error "Stack trace: $(stacktrace(catch_backtrace()))"
                error_count += 1
                continue
            end
        end
        
        @info "Batch channel summary complete!"
        @info "Processed $processed_count files successfully, $error_count errors"
        @info "Results saved to: $output_dir"
        
    finally
        close(global_logger())
    end
    
    return nothing
end
