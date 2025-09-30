"""
Refactored batch filtering with functional, composable style.
All helper functions included in this file for demonstration.
"""

#=============================================================================
    HELPER FUNCTIONS - Data structures for cleaner code
=============================================================================#

"""Immutable struct to hold batch processing results."""
struct BatchResult
    success::Bool
    filename::String
    message::String
end

"""Configuration for batch operations."""
struct BatchConfig
    file_pattern::String
    input_dir::String
    output_dir::String
    participants::Union{Int, Vector{Int}, Nothing}
    conditions::Union{Int, Vector{Int}, Nothing}
end

#=============================================================================
    HELPER FUNCTIONS - Pure functions (no side effects)
=============================================================================#

"""Find JLD2 files matching pattern and participant filter."""
function find_batch_files(pattern::String, dir::String; participants=nothing)
    all_files = readdir(dir)
    
    # Filter by pattern and extension
    files = Base.filter(all_files) do f
        endswith(f, ".jld2") && contains(f, pattern)
    end
    
    # Filter by participant if specified
    if participants !== nothing
        files = _filter_files(files; include=participants)
    end
    
    return files
end

"""Load EEG data from file, returning (data_var, var_name) or nothing."""
function load_eeg_data(filepath::String)
    file_data = load(filepath)
    
    # Try common variable names
    for var_name in ["erps", "epochs"]
        if haskey(file_data, var_name)
            return (file_data[var_name], var_name)
        end
    end
    
    return nothing
end

"""Filter data by conditions."""
function select_conditions(data, conditions)
    isnothing(conditions) && return data
    condition_nums = conditions isa Int ? [conditions] : conditions
    return data[condition_nums]
end

"""Validate filter parameters, returning error message or nothing."""
function validate_filter_params(cutoff_freq::Real, filter_type::String, input_dir::String)
    !isdir(input_dir) && return "Input directory does not exist: $input_dir"
    cutoff_freq <= 0 && return "Cutoff frequency must be positive, got: $cutoff_freq"
    filter_type ∉ ["lp", "hp"] && return "Filter type must be one of: lp, hp, got: $filter_type"
    return nothing
end

"""Generate default output directory name."""
function default_output_dir(input_dir::String, pattern::String, filter_type::String, freq::Real)
    joinpath(input_dir, "filtered_$(pattern)_$(filter_type)_$(freq)hz")
end

#=============================================================================
    HELPER FUNCTIONS - Core processing logic (pure)
=============================================================================#

"""
Process a single file through filtering pipeline.
Returns BatchResult with success/failure info.
"""
function process_filter_file(filepath::String, output_path::String, 
                             filter_type::String, cutoff_freq::Real,
                             conditions)
    filename = basename(filepath)
    
    # Load data
    data_result = load_eeg_data(filepath)
    if isnothing(data_result)
        return BatchResult(false, filename, "No recognized data variable")
    end
    
    data, var_name = data_result
    
    # Select conditions
    data = select_conditions(data, conditions)
    
    # Apply filter (mutates data in-place)
    for item in data
        filter_data!(item, filter_type, cutoff_freq)
    end
    
    # Save
    save(output_path, var_name, data)
    
    return BatchResult(true, filename, "Filtered successfully")
end

#=============================================================================
    HELPER FUNCTIONS - Orchestration (with side effects)
=============================================================================#

"""
Execute batch processing with logging and error handling.
Takes a pure processing function and maps it over files.
"""
function run_batch_operation(process_fn::Function, files::Vector{String}, 
                             input_dir::String, output_dir::String;
                             operation_name::String = "Processing")
    results = BatchResult[]
    
    for (i, file) in enumerate(files)
        @info "$operation_name: $file ($i/$(length(files)))"
        
        input_path = joinpath(input_dir, file)
        output_path = joinpath(output_dir, file)
        
        result = try
            process_fn(input_path, output_path)
        catch e
            @error "Error processing $file" exception=(e, catch_backtrace())
            BatchResult(false, file, "Exception: $(sprint(showerror, e))")
        end
        
        push!(results, result)
        
        # Log individual result
        if result.success
            @info "  ✓ $(result.message)"
        else
            @minimal_warning "  ✗ $(result.message)"
        end
    end
    
    return results
end

"""Log summary statistics from batch results."""
function log_batch_summary(results::Vector{BatchResult}, output_dir::String)
    n_success = count(r -> r.success, results)
    n_error = length(results) - n_success
    
    @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
    @info "Output saved to: $output_dir"
    
    return (success=n_success, errors=n_error)
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    filter(file_pattern::String, cutoff_freq::Real; 
           input_dir::String = pwd(), 
           filter_type::String = "lp", 
           participants::Union{Int, Vector{Int}, Nothing} = nothing,
           conditions::Union{Int, Vector{Int}, Nothing} = nothing,
           output_dir::Union{String, Nothing} = nothing)

Filter EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `cutoff_freq::Real`: Cutoff frequency in Hz
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `filter_type::String`: Type of filter ("lp", "hp")
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on filter settings)

# Example
```julia
# Filter all epochs at 30 Hz
filter("epochs", 30.0)

# Filter specific participant
filter("epochs", 30.0, participants=3)

# Filter specific participants and conditions
filter("epochs", 30.0, participants=[3, 4], conditions=[1, 2])
```
"""
function filter(file_pattern::String, cutoff_freq::Real; 
                input_dir::String = pwd(), 
                filter_type::String = "lp", 
                participants::Union{Int, Vector{Int}, Nothing} = nothing,
                conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                output_dir::Union{String, Nothing} = nothing)
    
    # Setup logging
    log_file = "filter_data.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch filtering started at $(now())"
        @log_call "filter" (file_pattern, cutoff_freq)
        
        # Validation (early return on error)
        if (error_msg = validate_filter_params(cutoff_freq, filter_type, input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        # Setup directories
        output_dir = something(output_dir, default_output_dir(input_dir, file_pattern, filter_type, cutoff_freq))
        mkpath(output_dir)
        
        # Find files
        files = find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Filter settings: $filter_type filter, cutoff: $cutoff_freq Hz"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> 
            process_filter_file(input_path, output_path, filter_type, cutoff_freq, conditions)
        
        # Execute batch operation
        results = run_batch_operation(process_fn, files, input_dir, output_dir; 
                                      operation_name="Filtering")
        
        # Log summary
        log_batch_summary(results, output_dir)
        
    finally
        # Cleanup logging
        close_global_logging()
        log_source = log_file
        log_dest = joinpath(output_dir, log_file)
        if log_source != log_dest
            mv(log_source, log_dest, force=true)
        end
    end
end
