"""
Batch computation of condition averages for ERP data.
"""

#=============================================================================
    AVERAGE-SPECIFIC VALIDATION
=============================================================================#

"""Generate default output directory name for condition averaging."""
function _condition_average_default_output_dir(input_dir::String, pattern::String, groups::Vector{Vector{Int}})
    groups_str = join([join(group, "-") for group in groups], "_")
    joinpath(input_dir, "averages_$(pattern)_$(groups_str)")
end

#=============================================================================
    AVERAGE-SPECIFIC PROCESSING
=============================================================================#

"""
Create an average wave by averaging ERP data across conditions.
"""
function _create_average_wave(erps::Vector{<:ErpData}, conditions::Vector{Int}, avg_cond::Int)
    # Validate that all ERPs have the same structure
    for i = 2:length(erps)
        _have_same_structure(erps[1], erps[i]) || @minimal_error_throw("ERPs have inconsistent structure")
    end

    # Get EEG channels (exclude metadata columns)
    metadata_cols = meta_labels(erps[1])
    eeg_channels = setdiff(propertynames(erps[1].data), metadata_cols)

    # Create a copy of first ERP's data for the average
    avg_data = copy(erps[1].data)

    # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
    cols_to_remove = [:condition, :condition_name, :n_epochs]
    for col in cols_to_remove
        if hasproperty(avg_data, col)
            select!(avg_data, Not(col))
        end
    end

    # Average EEG channels across conditions
    n = length(erps)
    for ch in eeg_channels
        all_have_ch = all(hasproperty(erp.data, ch) for erp in erps)
        if all_have_ch
            avg_data[!, ch] = sum(erp.data[!, ch] for erp in erps) ./ n
        else
            @minimal_warning "Channel $ch not found in all conditions, keeping values from first condition"
        end
    end

    # Sum n_epochs across all conditions (total epochs contributing to the average)
    total_epochs = sum(erp.n_epochs for erp in erps)
    avg_condition_name = "average_$(join(conditions, "_"))"

    return ErpData(
        erps[1].file,
        avg_cond,
        avg_condition_name,
        avg_data,
        erps[1].layout,
        erps[1].sample_rate,
        erps[1].analysis_info,
        total_epochs,
    )
end

"""
Process a single ERP file through condition averaging pipeline.
Returns BatchResult with success/failure info.
"""
function _condition_average_process_file(filepath::String, output_path::String, condition_groups::Vector{Vector{Int}})
    filename = basename(filepath)

    # Read data (using read_data which finds by type)
    file_data = read_data(filepath)

    if isnothing(file_data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Dispatch to the appropriate in-memory condition_average based on data type
    if file_data isa Vector{<:ErpData} || file_data isa Vector{<:TimeFreqData}
        try
            average_waves = condition_average(file_data, condition_groups)
            jldsave(output_path; data = average_waves)
            return BatchResult(true, filename, "Created $(length(average_waves)) average(s)")
        catch e
            return BatchResult(false, filename, "$(sprint(showerror, e))")
        end
    else
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{TimeFreqData}, got $(typeof(file_data))")
    end
end

#=============================================================================
    IN-MEMORY API FUNCTION
=============================================================================#

"""
    condition_average(data::Vector{<:ErpData}, condition_groups::Vector{Vector{Int}})::Vector{ErpData}

Create condition average waves.

For each group of condition indices, averages the EEG channels across those conditions.

# Arguments
- `data::Vector{<:ErpData}`: Vector of ErpData, one per condition
- `condition_groups::Vector{Vector{Int}}`: Groups of condition indices to average (e.g., `[[1, 2], [3, 4]]`)

# Returns
- `Vector{ErpData}`: New vector with one averaged ERP per group

# Example
```julia
erps = average_epochs(epochs)
avg_waves = condition_average(erps, [[1, 2], [3, 4]])
```
"""
function condition_average(data::Vector{<:ErpData}, condition_groups::Vector{Vector{Int}})::Vector{ErpData}

    @info "Creating condition average waves..."

    average_waves = ErpData[]

    for (group_idx, conditions) in enumerate(condition_groups)
        # Find ERP data for each condition in the group
        group_erps = ErpData[]
        all_found = true

        for cond in conditions
            erp_found = nothing
            for erp in data
                if erp.condition == cond
                    erp_found = erp
                    break
                end
            end

            if isnothing(erp_found)
                @minimal_warning "Condition $cond not found. Skipping group $conditions."
                all_found = false
                break
            end
            push!(group_erps, erp_found)
        end

        if !all_found
            continue
        end

        avg_wave = _create_average_wave(group_erps, conditions, group_idx)
        push!(average_waves, avg_wave)
    end

    if isempty(average_waves)
        @minimal_error_throw("No valid condition groups found")
    end

    @info "Created $(length(average_waves)) average wave(s)"
    return average_waves
end

#=============================================================================
    BATCH API FUNCTION
=============================================================================#

"""
    condition_average(file_pattern::String, condition_groups::Vector{Vector{Int}}; 
                      input_dir::String = pwd(), 
                      participant_selection::Function = participants(),
                      output_dir::Union{String, Nothing} = nothing)

Batch process ERP or time-frequency data files to create condition averages.

This function loads JLD2 files containing ERP or TimeFreqData, computes averages across
specified condition groups, and saves the results to a new directory. Data type is
auto-detected from each file.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "erps_original")
- `condition_groups::Vector{Vector{Int}}`: Groups of conditions to average (e.g., [[1, 2], [3, 4]])
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Average conditions 1 and 2 together, and 3 and 4 together
condition_average("erps_cleaned", [[1, 2], [3, 4]])

# Average all four conditions into one
condition_average("erps_cleaned", [[1, 2, 3, 4]])

# Process specific participants
condition_average("erps_cleaned", [[1, 2]], 
                  input_dir = "/path/to/data", 
                  participant_selection = participants([1, 2, 3]))

# Specify custom output directory
condition_average("erps_cleaned", [[1, 2], [3, 4]], 
                  input_dir = "/path/to/data", 
                  output_dir = "/path/to/output")
```
"""
function condition_average(
    file_pattern::String,
    condition_groups::Vector{Vector{Int}};
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "condition_average.log"
    setup_global_logging(log_file)

    result = (success = 0, errors = 0)  # Default return value
    try
        @info "Batch condition averaging started at $(now())"
        @log_call "condition_average"

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end


        if (error_msg = _validate_condition_groups(condition_groups)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _condition_average_default_output_dir(input_dir, file_pattern, condition_groups))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            result = (success = 0, errors = 0)
        else
            @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
            @info "Condition groups: $condition_groups"

            # Create processing function with captured parameters
            process_fn = (input_path, output_path) -> _condition_average_process_file(input_path, output_path, condition_groups)

            # Execute batch operation
            results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Creating averages")

            result = _log_batch_summary(results, output_dir)
        end

    finally
        _cleanup_logging(log_file, output_dir)
    end

    return result
end
