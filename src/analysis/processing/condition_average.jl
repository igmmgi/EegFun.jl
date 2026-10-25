"""
Batch computation of condition averages for ERP data.
"""

# === AVERAGE-SPECIFIC VALIDATION ===

"""Generate default output directory name for condition averaging."""
function _condition_average_default_output_dir(input_dir::String, pattern::String, groups::Vector{Vector{Int}})
    groups_str = join([join(group, "-") for group in groups], "_")
    joinpath(input_dir, "averages_$(clean_pattern(pattern))_$(groups_str)")
end

# === AVERAGE-SPECIFIC PROCESSING ===

"""
Create an average wave by averaging ERP data across conditions.
"""
function _create_average_wave(erps::AbstractVector{<:ErpData}, conditions::AbstractVector{<:Integer}, avg_cond::Int)
    # Validate that all ERPs have the same structure
    for i = 2:length(erps)
        _have_same_structure(erps[1], erps[i]) || @minimal_error("ERPs have inconsistent structure")
    end

    # Get EEG channels (exclude metadata columns)
    metadata_cols = meta_labels(erps[1])
    eeg_channels = setdiff(propertynames(erps[1].data), metadata_cols)

    # Create a copy of first ERP's data for the average
    avg_data = copy(erps[1].data)

    # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
    cols_to_remove = filter(col -> hasproperty(avg_data, col), [:condition, :condition_name, :n_epochs])
    if !isempty(cols_to_remove)
        select!(avg_data, Not(cols_to_remove))
    end

    # Average EEG channels across conditions
    n = length(erps)
    for ch in eeg_channels
        all_have_ch = all(hasproperty(erp.data, ch) for erp in erps)
        if all_have_ch
            col = avg_data[!, ch]
            for i = 2:n
                col .+= erps[i].data[!, ch]
            end
            col ./= n
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
        copy(erps[1].layout),
        erps[1].sample_rate,
        copy(erps[1].analysis_info),
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
    if file_data isa AbstractVector{<:ErpData} || file_data isa AbstractVector{<:TimeFreqData}
        average_waves = condition_average(file_data, condition_groups)
        jldsave(output_path; data = average_waves)
        return BatchResult(true, filename, "Created $(length(average_waves)) average(s)")
    else
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{TimeFreqData}, got $(typeof(file_data))")
    end
end

# === IN-MEMORY API FUNCTION ===

"""
    condition_average(data::Vector{<:ErpData}, condition_groups::Vector{Vector{Int}}) -> Vector{ErpData}
    condition_average(data::Vector{TimeFreqData}, condition_groups::Vector{Vector{Int}}) -> Vector{TimeFreqData}
    condition_average(file_pattern::String, condition_groups::Vector{Vector{Int}};
    input_dir::String = pwd(), participant_selection::Function = participants(),
    output_dir = nothing)

Create condition average waves (or time-frequency averages), or batch-process JLD2 files.

For each group of condition indices, averages the EEG channels (or TF power/phase) across those
conditions. The `file_pattern` method loads JLD2 files and saves results to a new directory.

# Arguments
- `data`: `Vector{<:ErpData}` or `Vector{TimeFreqData}`, one element per condition
- `condition_groups`: Groups of condition indices to average (e.g., `[[1, 2], [3, 4]]`)

# Example
```julia
erps = average_epochs(epochs)
avg_waves = condition_average(erps, [[1, 2], [3, 4]])

# Batch
condition_average("erps_unrejected", [[1, 2], [3, 4]])
```
"""
function condition_average(data::AbstractVector{<:ErpData}, condition_groups)
    groups = if isempty(condition_groups)
        condition_groups
    elseif first(condition_groups) isa Integer
        [condition_groups]
    else
        condition_groups
    end

    @info "Creating condition average waves..."

    average_waves = ErpData[]

    for (group_idx, conditions) in enumerate(groups)
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
        @minimal_error("No valid condition groups found")
    end

    @info "Created $(length(average_waves)) average wave(s)"
    return average_waves
end

# === BATCH API FUNCTION ===


function condition_average(
    file_pattern::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}};
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    _with_batch_logging("condition averaging", "condition_average.log", input_dir, output_dir) do ctx
        @log_call "condition_average"

        # Validation

        error_msg = _validate_condition_groups(condition_groups)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = _prepare_output_dir!(ctx, _condition_average_default_output_dir(input_dir, file_pattern, condition_groups))

        process_fn = (input_path, output_path) -> _condition_average_process_file(input_path, output_path, condition_groups)
        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Creating averages")
    end
end

condition_average(file_pattern::String, condition_group::AbstractVector{<:Integer}; kwargs...) =
    condition_average(file_pattern, [collect(Int, condition_group)]; kwargs...)
