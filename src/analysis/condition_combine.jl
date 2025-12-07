"""
Batch combining of conditions for epoch data.
"""

#=============================================================================
    COMBINE-CONDITIONS-SPECIFIC VALIDATION
=============================================================================#

"""Validate that file pattern is for epochs data."""
function _validate_epochs_pattern_combine(pattern::String)
    !contains(pattern, "epochs") &&
        return "condition_combine only works with epoch data. File pattern must contain 'epochs', got: '$pattern'"
    return nothing
end

"""Generate default output directory name for condition combining."""
function _condition_combine_default_output_dir(input_dir::String, pattern::String, groups::Vector{Vector{Int}})
    groups_str = join([join(group, "-") for group in groups], "_")
    joinpath(input_dir, "combined_$(pattern)_$(groups_str)")
end

#=============================================================================
    COMBINE-CONDITIONS-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single epochs file through condition combining pipeline.
Returns BatchResult with success/failure info.
"""
function _condition_combine_process_file(filepath::String, output_path::String, condition_groups::Vector{Vector{Int}})
    filename = basename(filepath)

    # Load data
    data = load_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of EpochData)
    if !(data isa Vector{<:EpochData})
        return BatchResult(false, filename, "Invalid data type: expected Vector{EpochData}")
    end
    max_condition = length(data)

    # Combine conditions for epochs
    combined_data = EpochData[]

    for (group_idx, original_conditions) in enumerate(condition_groups)
        # Validate that all requested conditions exist
        missing_conditions = Base.filter(c -> c > max_condition || c < 1, original_conditions)
        if !isempty(missing_conditions)
            return BatchResult(
                false,
                filename,
                "Condition(s) $missing_conditions not found (only has 1-$max_condition)",
            )
        end

        # Get data for the specified conditions
        condition_data = data[original_conditions]

        # For epochs: combine (concatenate) the conditions
        # Each condition is an EpochData object, we need to concatenate their data fields
        combined_data_frames = reduce(vcat, (epoch_data.data for epoch_data in condition_data))

        # Create new EpochData with concatenated data, using metadata from first condition
        # For combined conditions, use group_idx as condition number and create combined name
        combined_condition_name = join([cond.condition_name for cond in condition_data], "_")
        combined_epochs = EpochData(
            condition_data[1].file,  # All should have same file
            group_idx,  # Use group index as condition number
            combined_condition_name,
            combined_data_frames,
            condition_data[1].layout,
            condition_data[1].sample_rate,
            condition_data[1].analysis_info,
        )

        push!(combined_data, combined_epochs)
    end

    # Save (always use "data" as variable name since load_data finds by type)
    jldsave(output_path; data = combined_data)

    n_groups = length(condition_groups)
    total_epochs = sum(length(cond.data) for cond in combined_data)
    return BatchResult(true, filename, "Combined into $n_groups group(s) with $total_epochs total epochs")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    condition_combine(file_pattern::String, condition_groups::Vector{Vector{Int}}; 
                      input_dir::String = pwd(), 
                      participant_selection::Function = participants(),
                      output_dir::Union{String, Nothing} = nothing)

Combine EEG epoch conditions from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files (should contain "epochs")
- `condition_groups::Vector{Vector{Int}}`: Groups of condition numbers to combine (e.g., [[1,2], [3,4]])
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on groups)

# Example
```julia
# Combine conditions 1,2 into first group and 3,4 into second group
condition_combine("epochs", [[1, 2], [3, 4]])

# Combine specific participant
condition_combine("epochs_cleaned", [[1, 2], [3, 4]], participants=3)

# Combine with custom output directory
condition_combine("epochs", [[1, 2], [3, 4]], output_dir="/path/to/output/")
```

# Note
- Only works with epoch data (not ERPs)
- Conditions are combined (concatenated) into new conditions
- Use `average_epochs()` separately to create ERPs from combined epochs
"""
function condition_combine(
    file_pattern::String,
    condition_groups::Vector{Vector{Int}};
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "condition_combine.log"
    setup_global_logging(log_file)

    try

        @info ""
        @info "Batch condition_combine started at $(now())"
        @log_call "condition_combine"

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_epochs_pattern_combine(file_pattern)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Validate and clean condition groups (modifies in-place)
        if (error_msg = _validate_condition_groups(condition_groups)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir =
            something(output_dir, _condition_combine_default_output_dir(input_dir, file_pattern, condition_groups))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info ""
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Condition groups: $condition_groups\n"

        # Create processing function with captured parameters
        process_fn =
            (input_path, output_path) -> _condition_combine_process_file(input_path, output_path, condition_groups)

        # Execute batch operation
        results =
            _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Combining conditions")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

"""
Batch combining of conditions for epoch data.
"""
