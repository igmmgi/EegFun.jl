"""
Batch combining of conditions for epoch data.
"""

# === COMBINE-CONDITIONS-SPECIFIC VALIDATION ===

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

# === COMBINE-CONDITIONS-SPECIFIC PROCESSING ===

"""
Process a single epochs file through condition combining pipeline.
Returns BatchResult with success/failure info.
"""
function _condition_combine_process_file(filepath::String, output_path::String, condition_groups::Vector{Vector{Int}})
    filename = basename(filepath)

    # Read data
    data = read_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of EpochData)
    if !(data isa Vector{<:EpochData})
        return BatchResult(false, filename, "Invalid data type: expected Vector{EpochData}")
    end

    combined_data = condition_combine(data, condition_groups)

    # Save (always use "data" as variable name since read_data finds by type)
    jldsave(output_path; data = combined_data)

    n_groups = length(condition_groups)
    total_epochs = sum(length(cond.data) for cond in combined_data)
    return BatchResult(true, filename, "Combined into $n_groups group(s) with $total_epochs total epochs")
end

# === IN-MEMORY API FUNCTION ===

"""
    condition_combine(data::Vector{<:EpochData}, condition_groups::Vector{Vector{Int}})::Vector{EpochData}
    condition_combine(file_pattern::String, condition_groups::Vector{Vector{Int}};
    input_dir::String = pwd(), participant_selection::Function = participants(),
    output_dir = nothing)

Combine epoch conditions by concatenating epochs from specified groups.

For each group in `condition_groups`, epochs from those conditions are concatenated into a single
new `EpochData` object. The `file_pattern` method saves results to a new JLD2 directory.

# Arguments
- `data::Vector{<:EpochData}`: Vector of EpochData, one per condition
- `condition_groups::Vector{Vector{Int}}`: Groups of condition indices to combine (e.g., `[[1, 2], [3, 4]]`)

# Example
```julia
epochs = extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
combined = condition_combine(epochs, [[1, 2], [3, 4]])

# Batch
condition_combine("epochs", [[1, 2], [3, 4]])
```
"""
function condition_combine(data::Vector{<:EpochData}, condition_groups::Vector{Vector{Int}})::Vector{EpochData}

    @info "Combining conditions..."

    max_condition = length(data)
    combined_data = EpochData[]

    for (group_idx, original_conditions) in enumerate(condition_groups)
        # Validate that all requested conditions exist
        missing_conditions = filter(c -> c > max_condition || c < 1, original_conditions)
        if !isempty(missing_conditions)
            @minimal_error("Condition(s) $missing_conditions not found (only has 1-$max_condition)")
        end

        # Get data for the specified conditions
        condition_data = data[original_conditions]

        # Concatenate epoch DataFrames
        combined_data_frames = reduce(vcat, (epoch_data.data for epoch_data in condition_data))

        # Create combined condition name
        combined_condition_name = join([cond.condition_name for cond in condition_data], "_")

        combined_epochs = EpochData(
            condition_data[1].file,
            group_idx,
            combined_condition_name,
            combined_data_frames,
            copy(condition_data[1].layout),
            condition_data[1].sample_rate,
            copy(condition_data[1].analysis_info),
        )

        push!(combined_data, combined_epochs)
    end

    @info "Combined into $(length(combined_data)) group(s)"
    return combined_data
end

# === BATCH API FUNCTION ===


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
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        error_msg = _validate_epochs_pattern_combine(file_pattern)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Validate and clean condition groups (modifies in-place)
        error_msg = _validate_condition_groups(condition_groups)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _condition_combine_default_output_dir(input_dir, file_pattern, condition_groups))
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
        process_fn = (input_path, output_path) -> _condition_combine_process_file(input_path, output_path, condition_groups)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Combining conditions")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

"""
Batch combining of conditions for epoch data.
"""
