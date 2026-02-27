"""
    prepare_decoding(epochs::Vector{T}; condition_selection, channel_selection, interval_selection) where {T<:MultiDataFrameEeg}

Prepare epoch data for multivariate pattern analysis (MVPA/decoding).

Organizes epoch data into a structure suitable for decoding analysis by grouping
epochs by participant and selecting specified conditions for classification.
Works with both EpochData and TimeFreqEpochData.

# Arguments
- `epochs::Vector{T}`: Epoch data containing multiple conditions/participants (EpochData or TimeFreqEpochData)
- `condition_selection::Function`: Predicate to select conditions for classification (default: `conditions()` - all conditions)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()` - all channels)
- `interval_selection::Interval`: Time window as tuple (e.g., (0.0, 1.0)) or interval object (default: `times()` - all samples)

# Returns
- `Vector{Vector{T}}`: Vector of participant data, where each element is a vector of data for that participant's conditions

# Examples
```julia
# ERP decoding
all_epochs = read_all_data(EpochData, "epochs_good", input_dir, participants())
participant_epochs = prepare_decoding(
    all_epochs,
    condition_selection = conditions([1, 2]),
    channel_selection = channels(),
    interval_selection = (0.0, 1.0))

# TF decoding
all_tf_epochs = read_all_data(TimeFreqEpochData, "tf_epochs", input_dir, participants())
participant_tf = prepare_decoding(
    all_tf_epochs,
    condition_selection = conditions([1, 2]),
    interval_selection = times(0.0, 1.0))
```
"""
function prepare_decoding(
    epochs::Vector{T};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
) where {T<:MultiDataFrameEeg}

    isempty(epochs) && @minimal_error("Cannot prepare decoding with empty epochs vector")

    # Group all epochs by condition first
    epochs_by_condition = group_by_condition(epochs)

    # Apply condition selection to the sorted condition numbers
    all_cond_nums = collect(keys(epochs_by_condition))  # Already sorted by group_by_condition
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    # Validate at least 2 conditions for classification
    length(selected_cond_nums) >= 2 || @minimal_error(
        "Decoding requires at least 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select at least 2 conditions."
    )

    # Get selected conditions
    selected_conditions = [epochs_by_condition[cond_num] for cond_num in selected_cond_nums]

    # Extract participant IDs from filenames for all conditions
    participant_ids_per_condition =
        [[_extract_participant_id(basename(epoch.file)) for epoch in condition_epochs] for condition_epochs in selected_conditions]

    # Validate all conditions have the same participants
    first_participants = participant_ids_per_condition[1]
    for (cond_idx, participants) in enumerate(participant_ids_per_condition[2:end])
        sort(participants) != sort(first_participants) && @minimal_error(
            "Condition $(selected_cond_nums[cond_idx+1]) has different participants than condition $(selected_cond_nums[1]). " *
            "Decoding requires the same participants across all conditions."
        )
    end

    # Get unique participant list (sorted)
    unique_participants = sort(unique(first_participants))

    # Validate all epochs have same structure within each condition
    for (cond_idx, condition_epochs) in enumerate(selected_conditions)
        _have_same_structure(condition_epochs) ||
            @minimal_error("Condition $(selected_cond_nums[cond_idx]): Epochs have inconsistent structure")
    end

    # Validate structure is consistent across conditions
    for cond_idx = (firstindex(selected_conditions)+1):lastindex(selected_conditions)
        _have_same_structure(selected_conditions[1][1], selected_conditions[cond_idx][1]) || @minimal_error(
            "Condition $(selected_cond_nums[1]) vs $(selected_cond_nums[cond_idx]): " *
            "Epochs have inconsistent structure (different channels, sample rates, or time vectors)"
        )
    end

    # Apply channel and interval selection to all epochs
    selected_conditions = [
        subset(condition_epochs; channel_selection = channel_selection, interval_selection = interval_selection) for
        condition_epochs in selected_conditions
    ]

    # Validate selection produced data
    for (cond_idx, condition_epochs) in enumerate(selected_conditions)
        isempty(condition_epochs) && @minimal_error("Condition $(selected_cond_nums[cond_idx]): No data matched the selection criteria!")
        isempty(channel_labels(condition_epochs[1])) && @minimal_error("Channel selection produced no channels")
        isempty(time(condition_epochs[1])) && @minimal_error("Interval selection produced no time points")
    end

    # Organize by participant: for each participant, collect their data for all selected conditions
    participant_epochs = Vector{Vector{T}}()
    sizehint!(participant_epochs, length(unique_participants))

    for participant_id in unique_participants
        participant_data = Vector{T}()
        sizehint!(participant_data, length(selected_conditions))

        for condition_epochs in selected_conditions
            # Find this participant's epoch for this condition
            participant_epoch = findfirst(epoch -> _extract_participant_id(basename(epoch.file)) == participant_id, condition_epochs)

            if participant_epoch |> !isnothing
                push!(participant_data, condition_epochs[participant_epoch])
            else
                @minimal_error(
                    "Participant $participant_id is missing data for condition. " *
                    "All participants must have data for all selected conditions."
                )
            end
        end

        push!(participant_epochs, participant_data)
    end

    return participant_epochs
end

"""
    prepare_decoding(file_pattern::String; input_dir, participant_selection, condition_selection, channel_selection, interval_selection)

Load epoch data from JLD2 files and prepare for decoding (convenience wrapper).

Loads data from JLD2 files matching the pattern and prepares them for decoding analysis.
Works automatically with both EpochData and TimeFreqEpochData files.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "epochs_good" or "tf_epochs")

# Keyword Arguments
- `input_dir::String`: Directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Predicate to filter participants (default: `participants()`)
- `condition_selection::Function`: Predicate to select conditions (default: `conditions()`)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()`)
- `interval_selection::Interval`: Time window (default: `times()`)

# Returns
- `Vector{Vector{T}}`: Vector of participant data (T is automatically detected from files)

# Examples
```julia
# ERP decoding
participant_epochs = prepare_decoding(
    "epochs_good",
    input_dir = "/path/to/data",
    condition_selection = conditions([1, 2]),
    channel_selection = channels([:Fz, :Cz, :Pz]),
    interval_selection = (0.0, 1.0))

# TF decoding (same syntax, type detected automatically)
participant_tf = prepare_decoding(
    "tf_epochs",
    input_dir = "/path/to/data",
    condition_selection = conditions([1, 2]))
```
"""
function prepare_decoding(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
)
    # Load all data (auto-detects type from files)
    all_data = read_all_data(file_pattern, input_dir, participant_selection)
    isempty(all_data) && @minimal_error("No valid data found matching pattern '$file_pattern' in $input_dir")

    # Filter to epoch-based types suitable for decoding
    epochs = filter(d -> d isa MultiDataFrameEeg, all_data)
    isempty(epochs) && @minimal_error("No epoch data (EpochData or TimeFreqEpochData) found matching pattern '$file_pattern' in $input_dir")

    return prepare_decoding(
        epochs;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        interval_selection = interval_selection,
    )
end
