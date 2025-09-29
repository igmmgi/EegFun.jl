"""
    average_epochs_data(file_pattern::String; 
                       input_dir::String = pwd(), 
                       participants::Union{Int, Vector{Int}, Nothing} = nothing,
                       conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                       output_dir::Union{String, Nothing} = nothing)

Batch process epoch data files to create averaged ERP data.

This function loads JLD2 files containing epoch data, applies `average_epochs` to each condition,
and saves the resulting ERP data to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "epochs_cleaned", "epochs_original")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Returns
- `Nothing`: Saves averaged ERP data to output directory

# Examples
```julia
# Average all epoch files in current directory
average_epochs_data("epochs_cleaned")

# Process specific participants and conditions
average_epochs_data("epochs_cleaned", 
                   input_dir = "/path/to/data", 
                   participants = [1, 2, 3], 
                   conditions = [1, 2])

# Specify custom output directory
average_epochs_data("epochs_cleaned", 
                   input_dir = "/path/to/data", 
                   output_dir = "/path/to/output")
```
"""
function average_epochs_data(file_pattern::String; 
                            input_dir::String = pwd(), 
                            participants::Union{Int, Vector{Int}, Nothing} = nothing,
                            conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                            output_dir::Union{String, Nothing} = nothing)
    
    # Set up global logging
    log_file = "average_epochs_data.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch epoch averaging started at $(now())"
        
        # Log the function call generically
        args = [file_pattern]
        kwargs = [:input_dir => input_dir, :participants => participants, :conditions => conditions, :output_dir => output_dir]
        _log_function_call("average_epochs_data", args, kwargs)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        # Validate that file pattern contains "epochs"
        if !contains(file_pattern, "epochs")
            @minimal_error_throw("average_epochs_data only works with epoch data. File pattern must contain 'epochs', got: '$file_pattern'")
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            output_dir = joinpath(input_dir, "averaged_$(file_pattern)")
        end
        
        # Create output directory if it doesn't exist
        mkpath(output_dir)
        
        # Find JLD2 files matching the pattern
        all_files = readdir(input_dir)
        jld2_files = filter(x -> endswith(x, ".jld2") && contains(x, file_pattern), all_files)
        
        # Filter by participant number if specified
        if participants !== nothing
            jld2_files = _filter_files(jld2_files; include = participants)
        end
        
        if isempty(jld2_files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            @info "Available files: $(filter(x -> endswith(x, ".jld2"), all_files))"
            return nothing
        end
        
        @info "Found $(length(jld2_files)) JLD2 files matching pattern '$file_pattern'"
        
        processed_count = 0
        error_count = 0
        
        for (i, file) in enumerate(jld2_files)
            input_path = joinpath(input_dir, file)
            output_path = joinpath(output_dir, file)
            
            @info("Processing: $file ($(i)/$(length(jld2_files)))")
            
            try
                # Load the file and determine what data we have
                file_data = load(input_path)
                
                # Check what type of data we have
                if haskey(file_data, "epochs")
                    epochs_data = file_data["epochs"]
                    var_name = "epochs"
                else
                    @minimal_warning "No 'epochs' variable found in $file. Available: $(keys(file_data))"
                    continue
                end
                
                # Filter by condition if specified
                if conditions !== nothing
                    condition_nums = conditions isa Int ? [conditions] : conditions
                    # Filter data to only include specified conditions
                    epochs_data = epochs_data[condition_nums]
                    @info "  Filtering conditions: $condition_nums"
                end
                
                # Average epochs for each condition
                erps_data = EegData[]
                for (i, epoch_data) in enumerate(epochs_data)
                    @info "  Averaging epochs $i/$(length(epochs_data))"
                    erp_data = average_epochs(epoch_data)
                    push!(erps_data, erp_data)
                end
                
                # Save averaged data
                save(output_path, "erps", erps_data)
                
                @info "  Saved: $output_path"
                processed_count += 1
                
            catch e
                @error "Error processing $file: $e"
                error_count += 1
                continue
            end
        end
        
        @info "Epoch averaging complete! Processed $processed_count files successfully, $error_count errors, output saved to: $output_dir"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "average_epochs_data.log"
        log_dest = joinpath(output_dir, "average_epochs_data.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end
