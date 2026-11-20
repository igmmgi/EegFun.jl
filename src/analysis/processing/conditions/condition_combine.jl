# Pooling and concatenation of epoched trials across experimental conditions.

######################################
########### User Functions ###########
######################################

"""
    condition_combine(data::AbstractVector{<:EpochData}, condition_groups)::Vector{EpochData}
    condition_combine(file_pattern::String, condition_groups; input_dir = pwd(),
                      participant_selection = participants(), output_dir = nothing)

Combine `EpochData` conditions by concatenating independent copies of epochs within each group in
`condition_groups` (matching `.condition` numbers with positional fallback), or batch-process JLD2 files.
"""
function condition_combine(
    data::AbstractVector{<:EpochData},
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}},
)::Vector{EpochData}
    isempty(condition_groups) && @minimal_error("Condition groups cannot be empty")
    isempty(data) && @minimal_error("Cannot combine conditions: input data is empty")

    part_groups = _group_by_participant_file(data)
    !isnothing(part_groups) && return reduce(vcat, [condition_combine(pg, condition_groups) for pg in part_groups])

    @info "Combining conditions..."

    cond_nums = Int[d.condition for d in data]
    use_cond_field = allunique(cond_nums) && all(>(0), cond_nums)
    max_condition = length(data)
    combined_data = Vector{EpochData}(undef, length(condition_groups))

    for (group_idx, raw_conditions) in enumerate(condition_groups)
        isempty(raw_conditions) && @minimal_error("Condition group $group_idx cannot be empty")

        original_conditions = unique(Int, raw_conditions)
        if length(original_conditions) != length(raw_conditions)
            @minimal_warning "Group $group_idx contains duplicate conditions; only unique conditions will be used."
        end

        invalid_conditions = filter(<(1), original_conditions)
        !isempty(invalid_conditions) && @minimal_error("Condition(s) $invalid_conditions must be positive integers (>= 1)")

        condition_data = EpochData[]
        missing_conditions = Int[]
        for c in original_conditions
            if use_cond_field
                idx = findfirst(d -> d.condition == c, data)
                isnothing(idx) ? push!(missing_conditions, c) : push!(condition_data, data[idx])
            else
                (c < 1 || c > max_condition) ? push!(missing_conditions, c) : push!(condition_data, data[c])
            end
        end

        if !isempty(missing_conditions)
            avail_str = use_cond_field ? "available conditions: $cond_nums" : "only has 1-$max_condition"
            @minimal_error("Condition(s) $missing_conditions not found ($avail_str)")
        end

        for ed in @view condition_data[2:end]
            ed.sample_rate != condition_data[1].sample_rate && @minimal_error(
                "Cannot combine EpochData with different sample rates ($(condition_data[1].sample_rate) vs $(ed.sample_rate))"
            )
            channel_labels(ed) != channel_labels(condition_data[1]) &&
                @minimal_error("Cannot combine EpochData with different channel labels")
        end
        nonempty_conds = EpochData[ed for ed in condition_data if !isempty(ed.data)]
        length(nonempty_conds) > 1 &&
            (_have_same_structure(nonempty_conds) || @minimal_error("Cannot combine EpochData with inconsistent structure"))

        total_eps = sum(length(ed.data) for ed in condition_data)
        src_dfs = Vector{DataFrame}(undef, total_eps)
        pos = 0
        for ed in condition_data, df in ed.data
            pos += 1
            src_dfs[pos] = df
        end

        combined_condition_name = join([cond.condition_name for cond in condition_data], "_")
        combined_data_frames = Vector{DataFrame}(undef, total_eps)
        if total_eps >= 16 && Threads.nthreads() > 1
            Threads.@threads :static for i = 1:total_eps
                df_copy = copy(src_dfs[i])
                if hasproperty(df_copy, :epoch)
                    df_copy.epoch isa Vector{Int} ? fill!(df_copy.epoch, i) : (df_copy[!, :epoch] = fill(i, nrow(df_copy)))
                end
                _sync_condition_cols!(df_copy, group_idx, combined_condition_name)
                combined_data_frames[i] = df_copy
            end
        else
            for i = 1:total_eps
                df_copy = copy(src_dfs[i])
                if hasproperty(df_copy, :epoch)
                    df_copy.epoch isa Vector{Int} ? fill!(df_copy.epoch, i) : (df_copy[!, :epoch] = fill(i, nrow(df_copy)))
                end
                _sync_condition_cols!(df_copy, group_idx, combined_condition_name)
                combined_data_frames[i] = df_copy
            end
        end
        combined_data[group_idx] = EpochData(
            condition_data[1].file,
            group_idx,
            combined_condition_name,
            combined_data_frames,
            copy(condition_data[1].layout),
            condition_data[1].sample_rate,
            copy(condition_data[1].analysis_info),
        )
    end

    @info "Combined into $(length(combined_data)) group(s)"
    return combined_data
end

condition_combine(data::AbstractVector{<:EpochData}, condition_group::AbstractVector{<:Integer})::Vector{EpochData} =
    condition_combine(data, [collect(Int, condition_group)])

function condition_combine(
    file_pattern::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}};
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("condition_combine", "condition_combine.log", input_dir, output_dir) do ctx
        @log_call "condition_combine"

        error_msg = _validate_epochs_pattern_combine(file_pattern)
        !isnothing(error_msg) && @minimal_error(error_msg)

        for (group_idx, grp) in enumerate(condition_groups)
            isempty(grp) && @minimal_error("Condition group $group_idx cannot be empty")
        end

        cleaned_groups = Vector{Int}[collect(Int, grp) for grp in condition_groups]
        error_msg = _validate_condition_groups(cleaned_groups)
        !isnothing(error_msg) && @minimal_error(error_msg)

        output_dir = _prepare_output_dir!(ctx, _condition_combine_default_output_dir(input_dir, file_pattern, cleaned_groups))
        process_fn = (input_path, output_path) -> _condition_combine_process_file(input_path, output_path, cleaned_groups)
        return batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Combining conditions")
    end
end

condition_combine(file_pattern::String, condition_group::AbstractVector{<:Integer}; kwargs...) =
    condition_combine(file_pattern, [collect(Int, condition_group)]; kwargs...)

######################################
######### Internal Functions #########
######################################

function _validate_epochs_pattern_combine(pattern::String)
    if !contains(lowercase(pattern), "epoch") && occursin(r"(?:^|[^a-zA-Z])erps?(?:[^a-zA-Z]|$)"i, pattern)
        return "condition_combine only works with epoch data, not ERP files (got pattern: '$pattern'). Use condition_average for ERPs."
    end
    return nothing
end

function _condition_combine_default_output_dir(input_dir::String, pattern::String, groups::AbstractVector{<:AbstractVector{<:Integer}})
    groups_str = join([join(group, "-") for group in groups], "_")
    return joinpath(input_dir, "combined_$(_clean_pattern(pattern))_$(groups_str)")
end

function _condition_combine_process_file(
    filepath::String,
    output_path::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}},
)
    filename = basename(filepath)
    data = read_data(filepath)
    isnothing(data) && return BatchResult(false, filename, "No data variables found")
    !(data isa Vector{<:EpochData}) && return BatchResult(false, filename, "Invalid data type: expected Vector{EpochData}")

    combined_data = condition_combine(data, condition_groups)
    jldsave(output_path; data = combined_data, format_version = EEGFUN_FORMAT_VERSION)

    total_epochs = sum(n_epochs(cond) for cond in combined_data)
    return BatchResult(true, filename, "Combined into $(length(condition_groups)) group(s) with $total_epochs total epochs")
end
