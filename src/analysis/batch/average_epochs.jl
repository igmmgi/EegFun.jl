"""
Batch averaging of epoch data to create ERPs.
"""

#=============================================================================
    AVERAGE-SPECIFIC VALIDATION
=============================================================================#

"""Validate that file pattern is for epochs data."""
function _validate_epochs_pattern(pattern::String)
    !contains(pattern, "epochs") && 
        return "average_epochs only works with epoch data. File pattern must contain 'epochs', got: '$pattern'"
    return nothing
end

"""Generate default output directory name for averaging operation."""
function _default_average_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "averaged_$(pattern)")
end

#=============================================================================
    AVERAGE-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single epochs file through averaging pipeline.
Returns BatchResult with success/failure info.
"""
function _process_average_file(filepath::String, output_path::String, conditions)
    filename = basename(filepath)
    
    # Load data
    file_data = load(filepath)
    
    if !haskey(file_data, "epochs")
        return BatchResult(false, filename, "No 'epochs' variable found")
    end
    
    epochs_data = file_data["epochs"]
    
    # Select conditions
    epochs_data = _select_conditions(epochs_data, conditions)
    
    # Average epochs for each condition
    erps_data = average_epochs.(epochs_data)
    
    # Save
    save(output_path, "erps", erps_data)
    
    return BatchResult(true, filename, "Averaged $(length(erps_data)) condition(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    average_epochs(file_pattern::String; 
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

# Examples
```julia
# Average all epoch files in current directory
average_epochs("epochs_cleaned")

# Process specific participants and conditions
average_epochs("epochs_cleaned", 
              input_dir = "/path/to/data", 
              participants = [1, 2, 3], 
              conditions = [1, 2])

# Specify custom output directory
average_epochs("epochs_cleaned", 
                   input_dir = "/path/to/data", 
                   output_dir = "/path/to/output")
```
"""
function average_epochs(file_pattern::String; 
                       input_dir::String = pwd(), 
                       participants::Union{Int, Vector{Int}, Nothing} = nothing,
                       conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                       output_dir::Union{String, Nothing} = nothing)
    
    # Setup logging
    log_file = "average_epochs.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch epoch averaging started at $(now())"
        @log_call "average_epochs" (file_pattern,)
        
        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        if (error_msg = _validate_epochs_pattern(file_pattern)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        # Setup directories
        output_dir = something(output_dir, _default_average_output_dir(input_dir, file_pattern))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> 
            _process_average_file(input_path, output_path, conditions)
        
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name="Averaging")
        
        # Log summary
        _log_batch_summary(results, output_dir)
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end