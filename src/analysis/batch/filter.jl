"""
Batch filtering for EEG/ERP data.
"""

#=============================================================================
    FILTER-SPECIFIC VALIDATION
=============================================================================#

"""Validate filter-specific parameters, returning error message or nothing."""
function _validate_filter_params(cutoff_freq::Real, filter_type::String)
    cutoff_freq <= 0 && return "Cutoff frequency must be positive, got: $cutoff_freq"
    filter_type ∉ ["lp", "hp"] && return "Filter type must be one of: lp, hp, got: $filter_type"
    return nothing
end

"""Generate default output directory name for filter operation."""
function _default_filter_output_dir(input_dir::String, pattern::String, filter_type::String, freq::Real)
    joinpath(input_dir, "filtered_$(pattern)_$(filter_type)_$(freq)hz")
end

#=============================================================================
    FILTER-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through filtering pipeline.
Returns BatchResult with success/failure info.
"""
function _process_filter_file(filepath::String, output_path::String, 
                             filter_type::String, cutoff_freq::Real,
                             conditions)
    filename = basename(filepath)
    
    # Load data
    data_result = _load_eeg_data(filepath)
    if isnothing(data_result)
        return BatchResult(false, filename, "No recognized data variable")
    end
    
    data, var_name = data_result
    
    # Select conditions
    data = _select_conditions(data, conditions)
    
    # Apply filter (mutates data in-place)
    filter_data!.(data, filter_type, cutoff_freq)
    
    # Save
    save(output_path, var_name, data)
    
    return BatchResult(true, filename, "Filtered successfully")
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
    log_file = "filter.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch filtering started at $(now())"
        @log_call "filter" (file_pattern, cutoff_freq)
        
        # Validation 
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        if (error_msg = _validate_filter_params(cutoff_freq, filter_type)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        # Setup directories
        output_dir = something(output_dir, _default_filter_output_dir(input_dir, file_pattern, filter_type, cutoff_freq))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Create processing function with captured parameters
        @info "Filter settings: $filter_type filter, cutoff: $cutoff_freq Hz"
        process_fn = (input_path, output_path) -> 
            _process_filter_file(input_path, output_path, filter_type, cutoff_freq, conditions)
        
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name="Filtering")
        
        _log_batch_summary(results, output_dir)
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end
