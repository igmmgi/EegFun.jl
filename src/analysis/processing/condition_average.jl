"""
Condition averages for ERP and time-frequency data (in-memory and batch).
"""

# === AVERAGE-SPECIFIC VALIDATION ===

"""Generate default output directory name for condition averaging."""
function _condition_average_default_output_dir(input_dir::String, pattern::String, groups::Vector{Vector{Int}})
    groups_str = join([join(group, "-") for group in groups], "_")
    joinpath(input_dir, "averages_$(clean_pattern(pattern))_$(groups_str)")
end

# === AVERAGE-SPECIFIC PROCESSING ===

"""
Create an average wave by averaging data across conditions (ERP data, or TF power and phase).
"""
function _create_average_wave(dats::AbstractVector{T}, conditions::AbstractVector{<:Integer}, avg_cond::Int) where {T<:ConditionAverageData}
    # Validate that all datasets have the same structure
    for i = 2:length(dats)
        _have_same_structure(dats[1], dats[i]) || @minimal_error("Data have inconsistent structure")
    end

    # Copy the first dataset's channel data (all containers) and accumulate the others in place
    avg_frames = map(copy, _data_frames(dats[1]))
    n = length(dats)
    inv_n = 1.0 / n
    for ch in _condition_channels(dats[1])
        if all(d -> hasproperty(d.data, ch), dats)
            for (k, avg_df) in enumerate(avg_frames)
                col = avg_df[!, ch]
                for i = 2:n
                    col .+= _data_frames(dats[i])[k][!, ch]
                end
                col .*= inv_n
            end
        else
            @minimal_warning "Channel $ch not found in all conditions, keeping values from first condition"
        end
    end

    # ERP: total number of epochs contributing to the average
    n_epochs = dats[1] isa ErpData ? sum(d.n_epochs for d in dats) : 0
    return _derive_condition(dats[1], avg_cond, "average_$(join(conditions, "_"))", avg_frames; n_epochs = n_epochs)
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
function condition_average(data::AbstractVector{T}, condition_groups) where {T<:ConditionAverageData}
    groups = if isempty(condition_groups)
        condition_groups
    elseif first(condition_groups) isa Integer
        [condition_groups]
    else
        condition_groups
    end

    tag = _kind_tag(T)
    @info "Creating $(tag)condition average waves..."

    average_waves = T[]

    for (group_idx, conditions) in enumerate(groups)
        # Find data for each condition in the group (first match wins)
        group_data = T[]
        all_found = true

        for cond in conditions
            found_idx = findfirst(d -> d.condition == cond, data)
            if isnothing(found_idx)
                @minimal_warning "Condition $cond not found. Skipping group $conditions."
                all_found = false
                break
            end
            push!(group_data, data[found_idx])
        end

        all_found || continue

        avg_wave = _create_average_wave(group_data, conditions, group_idx)
        push!(average_waves, avg_wave)
    end

    if isempty(average_waves)
        @minimal_error("No valid condition groups found")
    end

    @info "Created $(length(average_waves)) $(tag)average wave(s)"
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
