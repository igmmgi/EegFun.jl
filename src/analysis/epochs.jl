"""
    parse_epoch_conditions(config::Dict) -> Vector{EpochCondition}

Parse epoch conditions from configuration dictionary.
"""
function parse_epoch_conditions(config::Dict)
    defaults = get(config, "epochs", Dict())

    conditions = EpochCondition[]
    condition_configs = get(defaults, "conditions", [])

    for condition_config in condition_configs
        name = condition_config["name"]

        # Parse trigger sequences (unified approach)
        trigger_sequences_raw = get(condition_config, "trigger_sequences", nothing)
        if trigger_sequences_raw === nothing
            @minimal_error("trigger_sequences must be specified for condition '$name'")
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
                @minimal_error("Invalid trigger sequence format: $sequence")
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
                @minimal_error(
                    "Both min_interval and max_interval must be specified when timing_pairs is specified for condition '$name'",
                )
            end
        end

        # Parse after/before constraints (optional)
        after = get(condition_config, "after", nothing)
        before = get(condition_config, "before", nothing)

        # Validation
        if reference_index < 1 || reference_index > length(trigger_sequences[1])
            @minimal_error(
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
                    @minimal_error(
                        "timing_pairs contains invalid indices for sequence of length $(length(trigger_sequences[1])) in condition '$name'",
                    )
                end
                if start_idx >= end_idx
                    @minimal_error("timing_pairs must have start_idx < end_idx for condition '$name'")
                end
            end
        end

        # Validate after/before constraints
        if after !== nothing && before !== nothing
            @minimal_error("Cannot specify both 'after' and 'before' constraints for condition '$name'")
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

    # For each trigger of interest
    for trigger in triggers_of_interest

        # Find all samples where this trigger occurs
        trigger_indices = findall(dat.data.triggers .== trigger)
        if isempty(trigger_indices)
            @minimal_warning "Trigger $trigger not found in data"
            continue
        end

        # For each trigger occurrence
        for idx in trigger_indices
            trigger_time = dat.data.time[idx]

            # Find all samples within the time window
            window_start = trigger_time + time_window[1]
            window_end = trigger_time + time_window[2]

            # Mark samples within the window
            in_window = (dat.data.time .>= window_start) .& (dat.data.time .<= window_end)
            dat.data[in_window, channel_out] .= true

        end

    end

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

    # For each epoch condition
    for condition in epoch_conditions

        # Find all occurrences of the trigger sequences (unified approach)
        sequence_indices = search_sequences(dat.data.triggers, condition.trigger_sequences)
        if isempty(sequence_indices)
            @minimal_warning "No triggers found for condition '$(condition.name)'"
            continue
        end

        # Apply after/before filtering if specified
        if condition.after !== nothing || condition.before !== nothing

            sequence_indices = filter(sequence_indices) do seq_start_idx
                # Check after constraint
                if condition.after !== nothing
                    found_after = any(dat.data.triggers[1:(seq_start_idx-1)] .== condition.after)
                    if !found_after
                        return false
                    end
                end

                # Check before constraint  
                if condition.before !== nothing
                    sequence_end = seq_start_idx + length(condition.trigger_sequence) - 1
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

            sequence_indices = filter(sequence_indices) do seq_start_idx
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

        # For each valid sequence occurrence, find the reference point (t=0 position)
        for seq_start_idx in sequence_indices

            # Calculate the reference index (t=0 position) within the sequence
            reference_idx = seq_start_idx + (condition.reference_index - 1)
            if reference_idx > length(dat.data.triggers)
                @minimal_warning "Reference index $(condition.reference_index) for condition '$(condition.name)' is out of bounds"
                continue
            end

            reference_time = dat.data.time[reference_idx]

            # Find all samples within the time window
            window_start = reference_time + time_window[1]
            window_end = reference_time + time_window[2]

            # Mark samples within the window
            in_window = (dat.data.time .>= window_start) .& (dat.data.time .<= window_end)
            dat.data[in_window, channel_out] .= true

        end
    end
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
        search_sequences(dat.data.triggers, epoch_condition.trigger_sequences) .+ (epoch_condition.reference_index - 1)
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

    # extract and create array of dataframes
    epochs = []
    for (epoch, (pre, zero, post)) in enumerate(zip(pre_idx, zero_idx, post_idx))
        epoch_df = DataFrame(dat.data[pre:post, :])
        epoch_df.time = epoch_df.time .- dat.data.time[zero]
        insertcols!(epoch_df, 4, :condition => condition)
        insertcols!(epoch_df, 5, :condition_name => epoch_condition.name)
        insertcols!(epoch_df, 6, :epoch => epoch)
        push!(epochs, epoch_df)
    end

    return EpochData(epochs, dat.layout, dat.sample_rate, dat.analysis_info)
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
    # Get all columns from the first epoch
    all_columns = propertynames(first(dat.data))
    numeric_columns = Symbol[]

    # Find all numeric columns (excluding Bool)
    for col in all_columns
        col_type = eltype(first(dat.data)[!, col])
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
    isempty(eeg_channels) && @minimal_error "No EEG channels found to average"

    # Concatenate all epochs
    all_epochs = reduce(vcat, dat.data)

    # Group by time and condition, then average EEG channels and count epochs
    erp = combine(
        groupby(all_epochs, [:time, :condition, :condition_name]),
        :epoch => length => :n_epochs,      # Count how many epochs contributed
        eeg_channels .=> mean .=> eeg_channels,  # Average the EEG channels
    )

    # Get the maximum number of epochs at any time point
    n_epochs = maximum(erp.n_epochs)

    return ErpData(erp, dat.layout, dat.sample_rate, dat.analysis_info, n_epochs)
end

"""
    remove_bad_epochs(dat::EpochData, bad_columns::Vector{Symbol})

Remove epochs that contain any true values in the specified boolean columns.

# Arguments
- `dat::EpochData`: The epoched data to filter
- `bad_columns::Vector{Symbol}`: Vector of column names to check for true values

# Returns
- `EpochData`: A new EpochData object with bad epochs removed

"""
function remove_bad_epochs(dat::EpochData, bad_columns::Vector{Symbol})
    # Validate that all bad_columns exist in the data
    for col in bad_columns
        if !hasproperty(first(dat.data), col)
            error("Column '$col' not found in epoch data")
        end
    end

    # Filter epochs: keep only epochs where ALL bad columns are false for ALL samples
    good_epochs = DataFrame[]

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
        end
    end

    # Return new EpochData with only good epochs
    return EpochData(good_epochs, dat.layout, dat.sample_rate, dat.analysis_info)
end

"""
    remove_bad_epochs(dat::EpochData, bad_column::Symbol)

Remove epochs that contain any true values in the specified boolean column.

# Arguments
- `dat::EpochData`: The epoched data to filter
- `bad_column::Symbol`: Column name to check for true values

# Returns
- `EpochData`: A new EpochData object with bad epochs removed

# Examples
```julia
# Remove epochs with extreme values
cleaned_epochs = remove_bad_epochs(epochs, :is_extreme_value_100)

# Remove epochs with EOG artifacts
cleaned_epochs = remove_bad_epochs(epochs, :is_vEOG)
```
"""
remove_bad_epochs(dat::EpochData, bad_column::Symbol) = remove_bad_epochs(dat, [bad_column])
