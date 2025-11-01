"""
    condition_parse_epoch(config::Dict) -> Vector{EpochCondition}

Parse epoch conditions from configuration dictionary.
"""
function condition_parse_epoch(config::Dict)
    defaults = get(config, "epochs", Dict())

    conditions = EpochCondition[]
    condition_configs = get(defaults, "conditions", [])

    for condition_config in condition_configs
        name = condition_config["name"]

        # Parse trigger sequences (unified approach)
        trigger_sequences_raw = get(condition_config, "trigger_sequences", nothing)
        if trigger_sequences_raw === nothing
            @minimal_error_throw("trigger_sequences must be specified for condition '$name'")
        end

        # Parse the unified format
        trigger_sequences = Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}()
        for sequence in trigger_sequences_raw
            if sequence isa Vector
                # Regular sequence - force type conversion
                converted_sequence = Vector{Union{Int,Symbol,UnitRange{Int}}}(sequence)
                push!(trigger_sequences, converted_sequence)
            elseif sequence isa UnitRange
                # Range - convert to single-element sequence
                push!(trigger_sequences, [sequence])
            else
                @minimal_error_throw("Invalid trigger sequence format: $sequence")
            end
        end

        # Parse reference index (single index) - default to 1
        reference_index = get(condition_config, "reference_index", 1)

        # Parse timing pairs (optional - if not specified, no timing constraints)
        timing_pairs_raw = get(condition_config, "timing_pairs", nothing)
        if timing_pairs_raw === nothing
            # No timing constraints
            timing_pairs = nothing
            min_interval = nothing
            max_interval = nothing
        else
            # Parse timing constraints
            timing_pairs = [(pair[1], pair[2]) for pair in timing_pairs_raw]
            min_interval = get(condition_config, "min_interval", nothing)
            max_interval = get(condition_config, "max_interval", nothing)

            # Validate that both min and max intervals are specified if timing_pairs is specified
            if min_interval === nothing || max_interval === nothing
                @minimal_error_throw(
                    "Both min_interval and max_interval must be specified when timing_pairs is specified for condition '$name'",
                )
            end
        end

        # Parse after/before constraints (optional)
        after = get(condition_config, "after", nothing)
        before = get(condition_config, "before", nothing)

        # Validation
        if reference_index < 1 || reference_index > length(trigger_sequences[1])
            @minimal_error_throw(
                "reference_index must be between 1 and $(length(trigger_sequences[1])) for condition '$name'"
            )
        end

        # Only validate timing constraints if they're specified
        if timing_pairs !== nothing
            if min_interval >= max_interval
                error("min_interval must be < max_interval for condition '$name'")
            end

            # Validate timing pairs
            for (start_idx, end_idx) in timing_pairs
                if start_idx < 1 ||
                   start_idx > length(trigger_sequences[1]) ||
                   end_idx < 1 ||
                   end_idx > length(trigger_sequences[1])
                    @minimal_error_throw(
                        "timing_pairs contains invalid indices for sequence of length $(length(trigger_sequences[1])) in condition '$name'",
                    )
                end
                if start_idx >= end_idx
                    @minimal_error_throw("timing_pairs must have start_idx < end_idx for condition '$name'")
                end
            end
        end

        # Validate after/before constraints
        if after !== nothing && before !== nothing
            @minimal_error_throw("Cannot specify both 'after' and 'before' constraints for condition '$name'")
        end

        push!(
            conditions,
            EpochCondition(
                name,
                trigger_sequences,
                reference_index,
                timing_pairs,
                min_interval,
                max_interval,
                after,
                before,
            ),
        )
    end

    return conditions
end

# Helper function to validate epoch window parameters
function _validate_epoch_window_params(dat::ContinuousData, time_window::Vector{<:Real})
    @assert length(time_window) == 2 "Time window must have exactly 2 elements"
    @assert time_window[1] <= time_window[2] "Time window start must be less than or equal to end"
    @assert hasproperty(dat.data, :triggers) "Data must have a triggers column"
    @assert hasproperty(dat.data, :time) "Data must have a time column"
    @assert !isempty(dat.data.time) "Time column cannot be empty"
    @assert !isempty(dat.data.triggers) "Triggers column cannot be empty"
    @assert issorted(dat.data.time) "Time column must be sorted in ascending order"
end

"""
    get_selected_epochs(dat::EpochData, epoch_selection::Function) -> Vector{Int}

Get the indices of epochs that match the epoch selection predicate.

# Arguments
- `dat::EpochData`: The EpochData object containing the epochs
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering

# Returns
- `Vector{Int}`: Indices of selected epochs

# Examples
```julia
# Get all epochs
selected = get_selected_epochs(dat, epochs())

# Get specific epochs
selected = get_selected_epochs(dat, epochs(1:10))

# Get epochs matching a condition
selected = get_selected_epochs(dat, epochs([1, 3, 5]))
```
"""
function get_selected_epochs(dat::EpochData, epoch_selection::Function)
    return findall(epoch_selection(1:length(dat.data)))
end




"""
    _mark_windows_at_indices!(dat::ContinuousData, reference_indices::Vector{Int}, time_window::Vector{<:Real}, channel_out::Symbol)

Internal helper function to mark time windows around specific reference indices.

# Arguments
- `dat`: ContinuousData object containing the EEG data
- `reference_indices`: Vector of sample indices to mark windows around
- `time_window`: Time window in seconds as a vector of two numbers
- `channel_out`: Symbol for the output column name

# Returns
- Number of windows marked
"""
function _mark_windows_at_indices!(
    dat::ContinuousData,
    reference_indices::Vector{Int},
    time_window::Vector{<:Real},
    channel_out::Symbol,
)::Int
    n_marked = 0

    for idx in reference_indices
        # Bounds check
        if idx < 1 || idx > length(dat.data.time)
            @minimal_warning "Reference index $idx is out of bounds, skipping"
            continue
        end

        reference_time = dat.data.time[idx]

        # Calculate window bounds
        window_start = reference_time + time_window[1]
        window_end = reference_time + time_window[2]

        # Mark samples within the window (vectorized for efficiency)
        in_window = (dat.data.time .>= window_start) .& (dat.data.time .<= window_end)
        dat.data[in_window, channel_out] .= true

        n_marked += sum(in_window)
    end

    return n_marked
end

"""
    mark_epoch_windows!(dat::ContinuousData, triggers_of_interest::Vector{Int}, time_window::Vector{<:Real}; 
                         channel_out::Symbol = :epoch_window)

Mark samples that are within a specified time window of triggers of interest.

# Arguments
- `dat`: ContinuousData object containing the EEG data
- `triggers_of_interest`: Vector of trigger values to mark windows around
- `time_window`: Time window in seconds as a vector of two numbers ([-1, 2] for window from -1 to 2 seconds)
- `channel_out`: Symbol for the output column name (default: :epoch_window)

# Returns
- The modified ContinuousData object with a new column indicating samples within trigger windows

# Examples

## Basic Usage
```julia
# Mark windows around trigger 1 (1 second before to 2 seconds after)
mark_epoch_windows!(dat, [1], [-1.0, 2.0])

# Mark windows around multiple (1, 3, and 5) triggers
mark_epoch_windows!(dat, [1, 3, 5], [-1.0, 2.0])

# Custom column name
mark_epoch_windows!(dat, [1, 3], [-0.5, 0.5], channel_out = :near_trigger)
```

## Different Time Windows
```julia
# Pre-stimulus baseline window
mark_epoch_windows!(dat, [1], [-0.2, 0.0], channel_out = :baseline_window)

# Post-stimulus response window
mark_epoch_windows!(dat, [1], [0.0, 0.8], channel_out = :response_window)

# Asymmetric window (more time after trigger)
mark_epoch_windows!(dat, [1], [-0.1, 1.0], channel_out = :asymmetric_window)

# Very short window for fast responses
mark_epoch_windows!(dat, [1], [-0.05, 0.3], channel_out = :fast_response)
```

## Multiple Conditions
```julia
# Different conditions with different windows
mark_epoch_windows!(dat, [1], [-0.2, 0.8], channel_out = :condition_1)
mark_epoch_windows!(dat, [2], [-0.2, 1.0], channel_out = :condition_2)
mark_epoch_windows!(dat, [3], [-0.2, 1.2], channel_out = :condition_3)
```

## Using with Sample Filtering
```julia
# Mark epoch windows
mark_epoch_windows!(dat, [1, 3], [-1.0, 2.0])

# Use in analysis (only samples within epochs)
# summary = channel_summary(dat, samples = samples(:epoch_window))
```
"""
function mark_epoch_windows!(
    dat::ContinuousData,
    triggers_of_interest::Vector{Int},
    time_window::Vector{<:Real};
    channel_out::Symbol = :epoch_window,
)
    # Input validation
    _validate_epoch_window_params(dat, time_window)

    # Initialize result vector with false 
    dat.data[!, channel_out] .= false

    # Collect all relevant trigger indices
    all_reference_indices = Int[]

    for trigger in triggers_of_interest
        trigger_indices = findall(dat.data.triggers .== trigger)
        if isempty(trigger_indices)
            @minimal_warning "Trigger $trigger not found in data"
            continue
        end
        append!(all_reference_indices, trigger_indices)
    end

    # Mark windows around all collected indices
    n_marked = _mark_windows_at_indices!(dat, all_reference_indices, time_window, channel_out)

    return dat
end




"""
    mark_epoch_windows!(dat::ContinuousData, epoch_conditions::Vector{EpochCondition}, time_window::Vector{<:Real}; 
                         channel_out::Symbol = :epoch_window)

Mark samples that are within a specified time window of trigger sequences defined by epoch conditions.

# Arguments
- `dat`: ContinuousData object containing the EEG data
- `epoch_conditions`: Vector of EpochCondition objects defining trigger sequences and reference points
- `time_window`: Time window in seconds as a vector of two numbers ([-1, 2] for window from -1 to 2 seconds)
- `channel_out`: Symbol for the output column name (default: :epoch_window)

# Returns
- The modified ContinuousData object with a new column indicating samples within trigger sequence windows

# Examples
```julia
# Define epoch conditions with wildcards and ranges
condition1 = EpochCondition(name="condition_1", trigger_sequence=[1, :any, 3], reference_index=2)
condition2 = EpochCondition(name="condition_2", trigger_ranges=[1:5, 10:15], reference_index=1)
conditions = [condition1, condition2]

# Mark windows around trigger sequences
mark_epoch_windows!(dat, conditions, [-1.0, 2.0])

# Custom column name
mark_epoch_windows!(dat, conditions, [-0.5, 0.5], channel_out = :near_sequences)
```
"""
function mark_epoch_windows!(
    dat::ContinuousData,
    epoch_conditions::Vector{EpochCondition},
    time_window::Vector{<:Real};
    channel_out::Symbol = :epoch_window,
)
    # Input validation
    _validate_epoch_window_params(dat, time_window)

    # Initialize result vector with false
    dat.data[!, channel_out] .= false

    # Collect all reference indices from all conditions
    all_reference_indices = Int[]

    # For each epoch condition
    for condition in epoch_conditions
        # Find all occurrences of the trigger sequences (unified approach)
        sequence_indices = search_sequence(dat.data.triggers, condition.trigger_sequences)
        if isempty(sequence_indices)
            @minimal_warning "No triggers found for condition '$(condition.name)'"
            continue
        end

        # Apply after/before filtering if specified
        if condition.after !== nothing || condition.before !== nothing
            sequence_indices = Base.filter(sequence_indices) do seq_start_idx
                # Check after constraint
                if condition.after !== nothing
                    found_after = any(dat.data.triggers[1:(seq_start_idx-1)] .== condition.after)
                    if !found_after
                        return false
                    end
                end

                # Check before constraint  
                if condition.before !== nothing
                    sequence_end = seq_start_idx + length(condition.trigger_sequences[1]) - 1
                    found_before = any(dat.data.triggers[(sequence_end+1):end] .== condition.before)
                    if !found_before
                        return false
                    end
                end

                return true
            end
            if isempty(sequence_indices)
                after_msg = condition.after !== nothing ? " after trigger $(condition.after)" : ""
                before_msg = condition.before !== nothing ? " before trigger $(condition.before)" : ""
                @minimal_warning "No trigger sequences found that meet position constraints$(after_msg)$(before_msg) for condition '$(condition.name)'"
                continue
            end
        end

        # Apply timing constraints if specified
        if condition.timing_pairs !== nothing &&
           condition.min_interval !== nothing &&
           condition.max_interval !== nothing

            sequence_indices = Base.filter(sequence_indices) do seq_start_idx
                for (start_idx, end_idx) in condition.timing_pairs
                    # Calculate actual indices in the trigger array
                    actual_start_idx = seq_start_idx + (start_idx - 1)
                    actual_end_idx = seq_start_idx + (end_idx - 1)

                    # Check bounds
                    if actual_start_idx < 1 || actual_end_idx > length(dat.data.triggers)
                        return false
                    end

                    # Calculate time interval between the two triggers
                    start_time_val = dat.data.time[actual_start_idx]
                    end_time_val = dat.data.time[actual_end_idx]
                    interval = end_time_val - start_time_val

                    # Check if interval meets constraints
                    if interval < condition.min_interval || interval > condition.max_interval
                        return false
                    end
                end
                return true
            end
            if isempty(sequence_indices)
                @minimal_warning "No trigger sequences found that meet timing constraints for condition '$(condition.name)'"
                continue
            end
        end

        # Convert sequence start indices to reference indices and collect them
        for seq_start_idx in sequence_indices
            reference_idx = seq_start_idx + (condition.reference_index - 1)
            if reference_idx <= length(dat.data.triggers)
                push!(all_reference_indices, reference_idx)
            else
                @minimal_warning "Reference index $(condition.reference_index) for condition '$(condition.name)' is out of bounds"
            end
        end
    end

    # Mark windows around all collected reference indices
    n_marked = _mark_windows_at_indices!(dat, all_reference_indices, time_window, channel_out)

    return dat
end






"""
    extract_epochs(dat::ContinuousData, condition::Int, epoch_condition::EpochCondition, start_time, end_time)

Extract epochs based on a single EpochCondition object, including timing validation, after/before filtering, 
trigger ranges, wildcard sequences, and multiple sequences.

# Arguments
- `dat::ContinuousData`: The continuous EEG data
- `condition::Int`: Condition number to assign to epochs
- `epoch_condition::EpochCondition`: EpochCondition object defining trigger sequence and timing constraints
- `start_time`: Start time relative to reference point (seconds)
- `end_time`: End time relative to reference point (seconds)

# Returns
- `EpochData`: The extracted epochs

# Examples
```julia
# Single sequence
condition = EpochCondition(name="single", trigger_sequence=[1, 2, 3], reference_index=2)

# Multiple sequences (OR logic)
condition = EpochCondition(
    name="multiple", 
    trigger_sequences=[[1, 2, 1], [1, 3, 1]],  # Match either [1,2,1] OR [1,3,1]
    reference_index=2
)

# Wildcard sequence
condition = EpochCondition(
    name="wildcard", 
    trigger_sequence=[1, :any, 3],  # :any matches any trigger value
    reference_index=2
)

# Trigger ranges
condition = EpochCondition(
    name="ranges", 
    trigger_ranges=[1:5, 10:15],  # Match triggers 1-5 or 10-15
    reference_index=1
)
```
"""
function extract_epochs(dat::ContinuousData, condition::Int, epoch_condition::EpochCondition, start_time, end_time)
    # Find t==0 positions based on trigger_sequences (unified approach)
    zero_idx =
        search_sequence(dat.data.triggers, epoch_condition.trigger_sequences) .+ (epoch_condition.reference_index - 1)
    isempty(zero_idx) && error("None of the trigger sequences $(epoch_condition.trigger_sequences) found!")

    # Apply after/before filtering if specified
    if epoch_condition.after !== nothing || epoch_condition.before !== nothing
        filtered_indices = Int[]

        for zero_pos in zero_idx
            # Find the position of the sequence start
            sequence_start = zero_pos - (epoch_condition.reference_index - 1)

            # Determine sequence length based on the first sequence (all should have same length for timing validation)
            sequence_length = length(epoch_condition.trigger_sequences[1])

            # Check if this sequence meets the after/before constraints
            valid_position = true

            if epoch_condition.after !== nothing
                # Check if there's a trigger with value 'after' before this sequence
                # Look backwards from sequence start to find the trigger
                found_after_trigger = false
                for i = (sequence_start-1):-1:1
                    if dat.data.triggers[i] == epoch_condition.after
                        found_after_trigger = true
                        break
                    end
                end
                if !found_after_trigger
                    valid_position = false
                end
            end

            if epoch_condition.before !== nothing
                # Check if there's a trigger with value 'before' after this sequence
                # Look forwards from sequence end to find the trigger
                sequence_end = sequence_start + sequence_length - 1
                found_before_trigger = false
                for i = (sequence_end+1):length(dat.data.triggers)
                    if dat.data.triggers[i] == epoch_condition.before
                        found_before_trigger = true
                        break
                    end
                end
                if !found_before_trigger
                    valid_position = false
                end
            end

            if valid_position
                push!(filtered_indices, zero_pos)
            end
        end

        zero_idx = filtered_indices
        if isempty(zero_idx)
            after_msg = epoch_condition.after !== nothing ? " after trigger $(epoch_condition.after)" : ""
            before_msg = epoch_condition.before !== nothing ? " before trigger $(epoch_condition.before)" : ""
            error(
                "No trigger sequences found that meet position constraints$(after_msg)$(before_msg) for condition '$(epoch_condition.name)'",
            )
        end
    end

    # Apply timing constraints if specified
    if epoch_condition.timing_pairs !== nothing &&
       epoch_condition.min_interval !== nothing &&
       epoch_condition.max_interval !== nothing

        valid_indices = Int[]

        for zero_pos in zero_idx
            # Check if this sequence meets all timing constraints
            sequence_start = zero_pos - (epoch_condition.reference_index - 1)
            valid_sequence = true

            for (start_idx, end_idx) in epoch_condition.timing_pairs
                # Calculate actual indices in the trigger array
                actual_start_idx = sequence_start + (start_idx - 1)
                actual_end_idx = sequence_start + (end_idx - 1)

                # Check bounds
                if actual_start_idx < 1 || actual_end_idx > length(dat.data.triggers)
                    valid_sequence = false
                    break
                end

                # Calculate time interval between the two triggers
                start_time_val = dat.data.time[actual_start_idx]
                end_time_val = dat.data.time[actual_end_idx]
                interval = end_time_val - start_time_val

                # Check if interval meets constraints
                if interval < epoch_condition.min_interval || interval > epoch_condition.max_interval
                    valid_sequence = false
                    break
                end
            end

            if valid_sequence
                push!(valid_indices, zero_pos)
            end
        end

        zero_idx = valid_indices
        isempty(zero_idx) &&
            error("No trigger sequences found that meet timing constraints for condition '$(epoch_condition.name)'")
    end

    # find number of samples pre/post epoch t = 0 position
    n_pre, n_post = find_idx_start_end(dat.data.time, abs(start_time), abs(end_time))
    pre_idx = zero_idx .- n_pre .+ 1
    post_idx = zero_idx .+ n_post .- 1

    # Extract and create array of dataframes with bounds checking
    epochs = DataFrame[]

    for (epoch, (pre, zero, post)) in enumerate(zip(pre_idx, zero_idx, post_idx))
        # Bounds checking to prevent out-of-bounds errors
        if pre < 1 || post > nrow(dat.data)
            throw(
                BoundsError(
                    dat.data,
                    "Epoch $epoch extends beyond data bounds (pre=$pre, post=$post, data_length=$(nrow(dat.data)))",
                ),
            )
        end

        epoch_df = DataFrame(dat.data[pre:post, :])
        epoch_df.time = epoch_df.time .- dat.data.time[zero]
        # Remove condition, condition_name, file columns if they exist (they're in struct now)
        # Keep epoch column since it represents original epoch number (needed after rejection)
        cols_to_remove = [:condition, :condition_name, :file]
        for col in cols_to_remove
            if hasproperty(epoch_df, col)
                select!(epoch_df, Not(col))
            end
        end
        # Add epoch number (original numbering from extraction)
        insertcols!(epoch_df, 4, :epoch => epoch)
        push!(epochs, epoch_df)
    end

    # Get file name from ContinuousData struct field
    file_name = dat.file
    
    return EpochData(file_name, condition, epoch_condition.name, epochs, dat.layout, dat.sample_rate, dat.analysis_info)
end

function extract_epochs(dat::ContinuousData, epoch_conditions::Vector{EpochCondition}, start_time, end_time)
    epochs = EpochData[]
    for (idx, epoch_condition) in enumerate(epoch_conditions)
        push!(epochs, extract_epochs(dat, idx, epoch_condition, start_time, end_time))
    end
    return epochs
end


"""
    average_epochs(dat::EpochData)

Average epochs to create an ERP. This function:
1. Concatenates all epochs
2. Groups by time point and condition
3. Averages the EEG channels at each time point
4. Adds a count of how many epochs went into each average

# Arguments
- `dat::EpochData`: The epoched data to average

# Returns
- `ErpData`: The averaged ERP data with epoch counts
"""
function average_epochs(dat::EpochData)
    # Input validation
    isempty(dat.data) && @minimal_error_throw("Cannot average empty EpochData")

    # Get all columns from the first epoch
    first_epoch = first(dat.data)
    all_columns = propertynames(first_epoch)
    numeric_columns = Symbol[]

    # Find all numeric columns (excluding Bool)
    for col in all_columns
        col_type = eltype(first_epoch[!, col])
        if col_type <: Number && col_type != Bool
            push!(numeric_columns, col)
        end
    end

    # Define columns that should not be averaged (metadata columns)
    # Keep :time for grouping, but don't average it
    metadata_columns = meta_labels(dat)

    # Get EEG channels to average (numeric columns minus metadata)
    eeg_channels = setdiff(numeric_columns, metadata_columns)

    # Ensure we have some channels to average
    isempty(eeg_channels) && @minimal_error_throw("No EEG channels found to average")

        # Concatenate all epochs with error handling
        try
            all_epochs = reduce(vcat, dat.data)

            # Verify we have time column
            if !hasproperty(all_epochs, :time)
                @minimal_error_throw("Missing required column 'time' for epoch averaging")
            end

            # Group by time only (condition/condition_name are constant in struct)
            # Count epochs by grouping by time
            erp = combine(
                groupby(all_epochs, [:time]),
                eeg_channels .=> mean .=> eeg_channels,  # Average the EEG channels
            )
            
            # Count epochs - each unique time point should have same number of epochs
            # We can check by seeing how many epochs we concatenated
            n_epochs = length(dat.data)

            return ErpData(dat.file, dat.condition, dat.condition_name, erp, dat.layout, dat.sample_rate, dat.analysis_info, n_epochs)
    catch e
        @minimal_error_throw("Failed to average epochs: $(e)")
    end
end


function average_epochs(dat::Vector{EpochData})
    return average_epochs.(dat)
end



"""
    reject_epochs!(dat::EpochData, info::EpochRejectionInfo)::EpochData

Remove epochs identified in `info` from `dat` in-place.

# Arguments
- `dat::EpochData`: Epoched EEG data to modify
- `info::EpochRejectionInfo`: Information about which epochs to reject

# Returns
- `EpochData`: The modified data (same as input)

# Examples
```julia
using eegfun, JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Detect bad epochs
info = detect_bad_epochs_automatic(epochs, z_criterion = 2.0)

# Remove the bad epochs
reject_epochs!(epochs, info)

# Save cleaned data
save("participant_1_cleaned.jld2", "epochs", epochs)
```
"""
function reject_epochs!(dat::EpochData, info::EpochRejectionInfo)::EpochData

    if isempty(info.rejected_epochs)
        @info "No epochs to reject"
        return dat
    end

    n_epochs = length(dat.data)
    rejected_indices = unique([r.epoch for r in info.rejected_epochs])
    epochs_to_keep = setdiff(1:n_epochs, rejected_indices)
    dat.data = dat.data[epochs_to_keep]

    @info "Rejected $(length(rejected_indices)) of $(info.epoch_info.n) epochs."

    return dat
end


"""
    reject_epochs(dat::EpochData, info::EpochRejectionInfo)::EpochData

Non-mutating version. Returns new `EpochData` with rejected epochs removed.

# Arguments
- `dat::EpochData`: Epoched EEG data (not modified)
- `info::EpochRejectionInfo`: Information about which epochs to reject

# Returns
- `EpochData`: New data with rejected epochs removed

# Examples
```julia
# Detect bad epochs
info = detect_bad_epochs_automatic(epochs, z_criterion = 2.0)

# Get cleaned copy (original preserved)
cleaned_epochs = reject_epochs(epochs, info)
```
"""
function reject_epochs(dat::EpochData, info::EpochRejectionInfo)::EpochData
    dat_copy = copy(dat)
    reject_epochs!(dat_copy, info)
    return dat_copy
end

"""
    reject_epochs(dat::Vector{EpochData}, info::Vector{EpochRejectionInfo})::Vector{EpochData}

Remove epochs identified in `info` from each `EpochData` in `dat`.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoched EEG data to filter
- `info::Vector{EpochRejectionInfo}`: Vector of rejection information (one per condition)

# Returns
- `Vector{EpochData}`: Vector of new EpochData objects with rejected epochs removed

# Examples
```julia
# Detect bad epochs for multiple conditions
rejection_info = detect_bad_epochs_automatic(epochs)

# Get cleaned copies (originals preserved)
cleaned_epochs = reject_epochs(epochs, rejection_info)
```
"""
reject_epochs(dat::Vector{EpochData}, info::Vector{EpochRejectionInfo})::Vector{EpochData} = reject_epochs.(dat, info)


"""
    reject_epochs(dat::EpochData, bad_columns::Vector{Symbol})

Remove epochs that contain any true values in the specified boolean columns.

# Arguments
- `dat::EpochData`: The epoched data to filter
- `bad_columns::Vector{Symbol}`: Vector of column names to check for true values

# Returns
- `EpochData`: A new EpochData object with bad epochs removed

"""
function reject_epochs(dat::EpochData, bad_columns::Vector{Symbol})
    # Input validation
    isempty(dat.data) && @minimal_error_throw("Cannot remove bad epochs from empty EpochData")
    isempty(bad_columns) && return dat  # No columns to check

    # Validate that all bad_columns exist in the data
    first_epoch = first(dat.data)
    for col in bad_columns
        if !hasproperty(first_epoch, col)
            @minimal_error_throw("Column '$col' not found in epoch data")
        end

        # Check if column is boolean
        col_type = eltype(first_epoch[!, col])
        if col_type != Bool
            @minimal_warning "Column '$col' is not Boolean type ($(col_type)). Non-zero values will be treated as 'true'"
        end
    end

    # Pre-allocate for better performance
    n_epochs = length(dat.data)
    good_epochs = DataFrame[]
    sizehint!(good_epochs, n_epochs)  # Performance hint

    n_removed = 0

    for epoch_df in dat.data
        # Check if any sample in this epoch has any true value in any bad column
        has_bad_samples = false

        for col in bad_columns
            if any(epoch_df[!, col])
                has_bad_samples = true
                break
            end
        end

        # If no bad samples found, keep this epoch
        if !has_bad_samples
            push!(good_epochs, epoch_df)
        else
            n_removed += 1
        end
    end

    # Log removal statistics
    if n_removed > 0
        @info "Condition $(dat.condition) ($(dat.condition_name)) removed $n_removed of $n_epochs epochs ($(round(100*n_removed/n_epochs, digits=1))%)"
    end

    # Return new EpochData with only good epochs (preserve struct fields)
    return EpochData(dat.file, dat.condition, dat.condition_name, good_epochs, dat.layout, dat.sample_rate, dat.analysis_info)
end

"""
    reject_epochs(dat::EpochData, bad_column::Symbol)

Remove epochs that contain any true values in the specified boolean column.

# Arguments
- `dat::EpochData`: The epoched data to filter
- `bad_column::Symbol`: Column name to check for true values

# Returns
- `EpochData`: A new EpochData object with bad epochs removed

# Examples
```julia
# Remove epochs with extreme values
cleaned_epochs = reject_epochs(epochs, :is_extreme_value_100)

# Remove epochs with EOG artifacts
cleaned_epochs = reject_epochs(epochs, :is_vEOG)
```
"""
reject_epochs(dat::EpochData, bad_column::Symbol) = reject_epochs(dat, [bad_column])

"""
    reject_epochs(dat::Vector{EpochData}, bad_column::Symbol)

Remove epochs that contain any true values in the specified boolean column for each EpochData in the vector.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoched data to filter
- `bad_column::Symbol`: Column name to check for true values

# Returns
- `Vector{EpochData}`: Vector of new EpochData objects with bad epochs removed

# Examples
```julia
# Remove epochs with extreme values from multiple datasets
cleaned_epochs = reject_epochs(epochs_list, :is_extreme_value_100)

# Remove epochs with EOG artifacts from multiple datasets
cleaned_epochs = reject_epochs(epochs_list, :is_vEOG)
```
"""
reject_epochs(dat::Vector{EpochData}, bad_column::Symbol) = reject_epochs.(dat, bad_column)

"""
    reject_epochs(dat::Vector{EpochData}, bad_columns::Vector{Symbol})

Remove epochs that contain any true values in the specified boolean columns for each EpochData in the vector.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoched data to filter
- `bad_columns::Vector{Symbol}`: Vector of column names to check for true values

# Returns
- `Vector{EpochData}`: Vector of new EpochData objects with bad epochs removed

# Examples
```julia
# Remove epochs with multiple artifact types from multiple datasets
cleaned_epochs = reject_epochs(epochs_list, [:is_extreme_value_100, :is_vEOG, :is_hEOG])
```
"""
reject_epochs(dat::Vector{EpochData}, bad_columns::Vector{Symbol}) = reject_epochs.(dat, bad_columns)


# ============================================================================ #
#                           EPOCH TABLE FUNCTIONS                             #
# ============================================================================ #

"""
    epochs_table(epochs::Vector{EpochData})

Display a pretty table showing epoch information to console and return the DataFrame.
"""
function epochs_table(epochs::Vector{EpochData}; print_table::Bool = true, kwargs...)
    isempty(epochs) && throw(ArgumentError("epochs vector cannot be empty"))

    df = _build_base_epochs_df(epochs)
    df.n_epochs = [n_epochs(epoch) for epoch in epochs]

    if print_table
        pretty_table(stdout, df; alignment = [:l, :r, :l, :r], kwargs...)
    end

    return df
end

"""
    epochs_table(epochs_original, epochs_cleaned)

Display comparison table between original and cleaned epochs to console and return DataFrame.
"""
function epochs_table(
    epochs_original::Vector{EpochData},
    epochs_cleaned::Vector{EpochData};
    print_table::Bool = true,
    kwargs...,
)
    length(epochs_original) != length(epochs_cleaned) &&
        throw(ArgumentError("epochs_original and epochs_cleaned must have same length"))

    df = _build_base_epochs_df(epochs_original)
    df.n_epochs_original = [n_epochs(epoch) for epoch in epochs_original]
    df.n_epochs_cleaned = [n_epochs(epoch) for epoch in epochs_cleaned]
    df.percentage = round.((df.n_epochs_cleaned ./ df.n_epochs_original) .* 100; digits = 1)

    if print_table
        pretty_table(stdout, df; alignment = [:l, :r, :l, :r, :r, :r], kwargs...)
    end

    return df
end

# Helper functions to reduce repetition
function _build_base_epochs_df(epochs::Vector{EpochData})::DataFrame
    return DataFrame(
        file = [filename(epoch) for epoch in epochs],
        condition = [condition_number(epoch) for epoch in epochs],
        condition_name = [condition_name(epoch) for epoch in epochs],
    )
end


"""
    log_epochs_table(message::String, epochs...; kwargs...)

Log an epochs table with message and return the DataFrame.
Combines logging and table creation in one clean call.
"""
function log_epochs_table(epochs...; print_table::Bool = false, kwargs...)
    df = epochs_table(epochs...; print_table = print_table);
    table_output = sprint() do output_io
        io_context = IOContext(output_io, :displaysize => displaysize(stdout))
        pretty_table(io_context, df; alignment = [:l, :r, :l, :r, :r, :r], kwargs...)
    end
    @info "\n\n$table_output\n"
    return df
end

"""
Batch averaging of epoch data to create ERPs.
"""

#=============================================================================
    AVERAGE-SPECIFIC VALIDATION
=============================================================================#

"""Validate that file pattern is for epochs data."""
function _validate_epochs_pattern(pattern::String)
    !contains(pattern, "epochs") &&
        return "average_epochs only works with epoch data. File pattern must contain 'epochs', got: '$pattern'"
    return nothing
end

"""Generate default output directory name for averaging operation."""
function _default_average_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "averaged_$(pattern)")
end

#=============================================================================
    AVERAGE-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single epochs file through averaging pipeline.
Returns BatchResult with success/failure info.
"""
function _process_average_file(filepath::String, output_path::String, conditions)
    filename = basename(filepath)

    # Load data
    file_data = load(filepath)

    if !haskey(file_data, "epochs")
        return BatchResult(false, filename, "No 'epochs' variable found")
    end

    epochs_data = file_data["epochs"]

    # Select conditions
    epochs_data = _condition_select(epochs_data, conditions)

    # Average epochs for each condition
    erps_data = average_epochs.(epochs_data)

    # Save
    save(output_path, "erps", erps_data)

    return BatchResult(true, filename, "Averaged $(length(erps_data)) condition(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    average_epochs(file_pattern::String; 
                       input_dir::String = pwd(), 
                       participants::Union{Int, Vector{Int}, Nothing} = nothing,
                       conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                       output_dir::Union{String, Nothing} = nothing)

Batch process epoch data files to create averaged ERP data.

This function loads JLD2 files containing epoch data, applies `average_epochs` to each condition,
and saves the resulting ERP data to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "epochs_cleaned", "epochs_original")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Average all epoch files in current directory
average_epochs("epochs_cleaned")

# Process specific participants and conditions
average_epochs("epochs_cleaned", 
              input_dir = "/path/to/data", 
              participants = [1, 2, 3], 
              conditions = [1, 2])

# Specify custom output directory
average_epochs("epochs_cleaned", 
                   input_dir = "/path/to/data", 
                   output_dir = "/path/to/output")
```
"""
function average_epochs(
    file_pattern::String;
    input_dir::String = pwd(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    conditions::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "average_epochs.log"
    setup_global_logging(log_file)

    try
        @info "Batch epoch averaging started at $(now())"
        @log_call "average_epochs" (file_pattern,)

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_epochs_pattern(file_pattern)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_average_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> _process_average_file(input_path, output_path, conditions)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Averaging")

        # Log summary
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
