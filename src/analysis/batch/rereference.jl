"""
Batch rereferencing for EEG/ERP data.
"""

#=============================================================================
    REREFERENCE-SPECIFIC VALIDATION
=============================================================================#

"""Generate default output directory name for rereferencing operation."""
function _default_rereference_output_dir(input_dir::String, pattern::String, ref_selection)
    ref_str = ref_selection isa Symbol ? string(ref_selection) : join(ref_selection, "_")
    joinpath(input_dir, "rereferenced_$(pattern)_$(ref_str)")
end

#=============================================================================
    REREFERENCE-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through rereferencing pipeline.
Returns BatchResult with success/failure info.
"""
function _process_rereference_file(filepath::String, output_path::String, 
                                  reference_selection, conditions)
    filename = basename(filepath)
    
    # Load data
    data_result = _load_eeg_data(filepath)
    if isnothing(data_result)
        return BatchResult(false, filename, "No recognized data variable")
    end
    
    data, var_name = data_result
    
    # Select conditions
    data = _select_conditions(data, conditions)
    
    # Apply rereferencing (mutates data in-place)
    rereference!.(data, reference_selection)
    
    # Save
    save(output_path, var_name, data)
    
    ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
    return BatchResult(true, filename, "Rereferenced to $ref_str")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    rereference(file_pattern::String; 
               input_dir::String = pwd(), 
               reference_selection::Union{Symbol, Vector{Symbol}} = :avg, 
               participants::Union{Int, Vector{Int}, Nothing} = nothing,
               conditions::Union{Int, Vector{Int}, Nothing} = nothing,
               output_dir::Union{String, Nothing} = nothing)

Rereference EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `reference_selection::Union{Symbol, Vector{Symbol}}`: Reference channels, can be:
  - Special symbols: `:avg` (average reference), `:mastoid` (M1+M2)
  - Channel names: `[:Cz]`, `[:M1, :M2]`, etc.
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on reference settings)

# Example
```julia
# Rereference all epochs to average reference
rereference("epochs")

# Rereference specific participant to mastoid reference
rereference("epochs", reference_selection=:mastoid, participants=3)

# Rereference to Cz for specific participants and conditions
rereference("epochs", reference_selection=[:Cz], participants=[3, 4], conditions=[1, 2])
```
"""
function rereference(file_pattern::String; 
                    input_dir::String = pwd(), 
                    reference_selection::Union{Symbol, Vector{Symbol}} = :avg, 
                    participants::Union{Int, Vector{Int}, Nothing} = nothing,
                    conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                    output_dir::Union{String, Nothing} = nothing)
    
    # Setup logging
    log_file = "rereference.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch rereferencing started at $(now())"
        @log_call "rereference" (file_pattern,)
        
        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        # Setup directories
        output_dir = something(output_dir, _default_rereference_output_dir(input_dir, file_pattern, reference_selection))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Reference settings: $ref_str"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> 
            _process_rereference_file(input_path, output_path, reference_selection, conditions)
        
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; 
                                      operation_name="Rereferencing")
        
        # Log summary
        _log_batch_summary(results, output_dir)
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end