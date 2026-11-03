"""
    prepare_decoding(epochs::AbstractVector{T}; condition_selection = conditions(),
                     channel_selection = channels(), interval_selection = times()) where {T<:MultiDataFrameEeg}
    prepare_decoding(file_pattern::String; input_dir = pwd(), participant_selection = participants(),
                     condition_selection = conditions(), channel_selection = channels(),
                     interval_selection = times())

Organize multi-participant epoch data (`EpochData` or `TimeFreqEpochData`) into per-participant
condition vectors `Vector{Vector{T}}` ready for [`decode_libsvm`](@ref).

Validates consistent participants, channels, sample rates, and time vectors across selected conditions.
"""
function prepare_decoding(
    epochs::AbstractVector{T};
    condition_selection = conditions(),
    channel_selection = channels(),
    interval_selection::Interval = times(),
) where {T<:MultiDataFrameEeg}
    isempty(epochs) && @minimal_error("Cannot prepare decoding with empty epochs vector")

    epochs_by_condition = group_by_condition(epochs)
    all_cond_nums = collect(keys(epochs_by_condition))
    sample_epochs = [epochs_by_condition[c][1] for c in all_cond_nums]
    selected_cond_nums = all_cond_nums[get_selected_conditions(sample_epochs, condition_selection)]

    if length(selected_cond_nums) < 2
        @minimal_error(
            "Decoding requires 2+ conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select 2+ conditions."
        )
    end

    selected_conditions = [epochs_by_condition[c] for c in selected_cond_nums]
    pids_per_cond = [Int[_extract_participant_id(basename(ep.file)) for ep in cond_eps] for cond_eps in selected_conditions]

    sorted_first_pids = sort(pids_per_cond[1])
    for c = 2:length(selected_conditions)
        sort(pids_per_cond[c]) != sorted_first_pids && @minimal_error(
            "Condition $(selected_cond_nums[c]) has different participants than condition $(selected_cond_nums[1]). " *
            "Decoding requires the same participants across all conditions."
        )
    end

    for (c, cond_eps) in enumerate(selected_conditions)
        _have_same_structure(cond_eps) || @minimal_error("Condition $(selected_cond_nums[c]): Epochs have inconsistent structure")
    end
    for c = 2:length(selected_conditions)
        _have_same_structure(selected_conditions[1][1], selected_conditions[c][1]) || @minimal_error(
            "Condition $(selected_cond_nums[1]) vs $(selected_cond_nums[c]): " *
            "Epochs have inconsistent structure (different channels, sample rates, or time vectors)"
        )
    end

    selected_conditions = [
        subset(cond_eps; channel_selection = channel_selection, interval_selection = interval_selection) for cond_eps in selected_conditions
    ]

    for (c, cond_eps) in enumerate(selected_conditions)
        isempty(cond_eps) && @minimal_error("Condition $(selected_cond_nums[c]): No data matched the selection criteria!")
        isempty(channel_labels(cond_eps[1])) && @minimal_error("Channel selection produced no channels")
        isempty(time_vector(cond_eps[1])) && @minimal_error("Interval selection produced no time points")
    end

    unique_pids = unique(sorted_first_pids)
    participant_epochs = Vector{Vector{T}}(undef, length(unique_pids))
    for (p_idx, pid) in enumerate(unique_pids)
        p_data = Vector{T}(undef, length(selected_conditions))
        for (c, cond_eps) in enumerate(selected_conditions)
            ep_idx = findfirst(==(pid), pids_per_cond[c])
            isnothing(ep_idx) && @minimal_error(
                "Participant $pid is missing data for condition. All participants must have data for all selected conditions."
            )
            p_data[c] = cond_eps[ep_idx]
        end
        participant_epochs[p_idx] = p_data
    end

    return participant_epochs
end

function prepare_decoding(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    channel_selection = channels(),
    interval_selection::Interval = times(),
)
    all_data = read_all_data(joinpath(input_dir, file_pattern), participant_selection)
    isempty(all_data) && @minimal_error("No valid data found matching pattern '$file_pattern' in $input_dir")

    epochs = filter(d -> d isa MultiDataFrameEeg, all_data)
    isempty(epochs) && @minimal_error("No epoch data (EpochData or TimeFreqEpochData) found matching pattern '$file_pattern' in $input_dir")

    return prepare_decoding(
        epochs;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        interval_selection = interval_selection,
    )
end
