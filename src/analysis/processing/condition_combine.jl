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
function _condition_combine_default_output_dir(input_dir::String, pattern::String, groups::AbstractVector{<:AbstractVector{<:Integer}})
    groups_str = join([join(group, "-") for group in groups], "_")
    joinpath(input_dir, "combined_$(clean_pattern(pattern))_$(groups_str)")
end

# === COMBINE-CONDITIONS-SPECIFIC PROCESSING ===

"""
Process a single epochs file through condition combining pipeline.
Returns BatchResult with success/failure info.
"""
function _condition_combine_process_file(
    filepath::String,
    output_path::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}},
)
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
    jldsave(output_path; data = combined_data, format_version = EEGFUN_FORMAT_VERSION)

    n_groups = length(condition_groups)
    total_epochs = sum(n_epochs(cond) for cond in combined_data)
    return BatchResult(true, filename, "Combined into $n_groups group(s) with $total_epochs total epochs")
end

# === IN-MEMORY API FUNCTION ===

"""
    condition_combine(data::AbstractVector{<:EpochData}, condition_groups::AbstractVector{<:AbstractVector{<:Integer}})::Vector{EpochData}
    condition_combine(file_pattern::String, condition_groups::AbstractVector{<:AbstractVector{<:Integer}};
    input_dir::String = pwd(), participant_selection::Function = participants(),
    output_dir = nothing)

Combine epoch conditions by concatenating epochs from specified groups.

For each group in `condition_groups`, epochs from those conditions are concatenated into a single
new `EpochData` object. Conditions are matched by `.condition` number (falling back to positional
indexing only if `.condition` numbers are non-unique or non-positive), and epoch DataFrames are
copied so the returned `EpochData` objects never alias the input DataFrames.
The `file_pattern` method saves results to a new JLD2 directory.

# Arguments
- `data::AbstractVector{<:EpochData}`: Vector of EpochData, one per condition
- `condition_groups::AbstractVector{<:AbstractVector{<:Integer}}`: Groups of condition numbers to combine (e.g., `[[1, 2], [3, 4]]`)

# Example
```julia
epochs = extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
combined = condition_combine(epochs, [[1, 2], [3, 4]])

# Batch
condition_combine("epochs", [[1, 2], [3, 4]])
```
"""
function condition_combine(
    data::AbstractVector{<:EpochData},
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}},
)::Vector{EpochData}
    isempty(condition_groups) && @minimal_error("Condition groups cannot be empty")
    isempty(data) && @minimal_error("Cannot combine conditions: input data is empty")

    @info "Combining conditions..."

    cond_nums = Int[d.condition for d in data]
    use_cond_field = allunique(cond_nums) && all(>(0), cond_nums)
    max_condition = length(data)
    combined_data = EpochData[]

    for (group_idx, raw_conditions) in enumerate(condition_groups)
        isempty(raw_conditions) && @minimal_error("Condition group $group_idx cannot be empty")

        original_conditions = unique(Int, raw_conditions)
        if length(original_conditions) != length(raw_conditions)
            @minimal_warning "Group $group_idx contains duplicate conditions; only unique conditions will be used."
        end

        invalid_conditions = filter(<(1), original_conditions)
        if !isempty(invalid_conditions)
            @minimal_error("Condition(s) $invalid_conditions must be positive integers (>= 1)")
        end

        condition_data = EpochData[]
        missing_conditions = Int[]
        for c in original_conditions
            if use_cond_field
                idx = findfirst(d -> d.condition == c, data)
                if isnothing(idx)
                    push!(missing_conditions, c)
                else
                    push!(condition_data, data[idx])
                end
            else
                if c < 1 || c > max_condition
                    push!(missing_conditions, c)
                else
                    push!(condition_data, data[c])
                end
            end
        end

        if !isempty(missing_conditions)
            avail_str = use_cond_field ? "available conditions: $cond_nums" : "only has 1-$max_condition"
            @minimal_error("Condition(s) $missing_conditions not found ($avail_str)")
        end

        # Concatenate independent copies of epoch DataFrames
        combined_data_frames = DataFrame[copy(df) for epoch_data in condition_data for df in epoch_data.data]

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

condition_combine(data::AbstractVector{<:EpochData}, condition_group::AbstractVector{<:Integer})::Vector{EpochData} =
    condition_combine(data, [collect(Int, condition_group)])

# === BATCH API FUNCTION ===


function condition_combine(
    file_pattern::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}};
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    _with_batch_logging("condition_combine", "condition_combine.log", input_dir, output_dir) do ctx
        @log_call "condition_combine"

        # Validation
        error_msg = _validate_epochs_pattern_combine(file_pattern)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        for (group_idx, grp) in enumerate(condition_groups)
            isempty(grp) && @minimal_error("Condition group $group_idx cannot be empty")
        end

        # Validate and clean condition groups
        cleaned_groups = Vector{Int}[collect(Int, grp) for grp in condition_groups]
        error_msg = _validate_condition_groups(cleaned_groups)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = _prepare_output_dir!(ctx, _condition_combine_default_output_dir(input_dir, file_pattern, cleaned_groups))

        process_fn = (input_path, output_path) -> _condition_combine_process_file(input_path, output_path, cleaned_groups)
        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Combining conditions")
    end
end

condition_combine(file_pattern::String, condition_group::AbstractVector{<:Integer}; kwargs...) =
    condition_combine(file_pattern, [collect(Int, condition_group)]; kwargs...)
