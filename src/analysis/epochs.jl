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
            error("trigger_sequences must be specified for condition '$name'")
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
                error("Invalid trigger sequence format: $sequence")
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
                error(
                    "Both min_interval and max_interval must be specified when timing_pairs is specified for condition '$name'",
                )
            end
        end

        # Parse after/before constraints (optional)
        after = get(condition_config, "after", nothing)
        before = get(condition_config, "before", nothing)

        # Validation
        if reference_index < 1 || reference_index > length(trigger_sequences[1])
            error("reference_index must be between 1 and $(length(trigger_sequences[1])) for condition '$name'")
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
                    error(
                        "timing_pairs contains invalid indices for sequence of length $(length(trigger_sequences[1])) in condition '$name'",
                    )
                end
                if start_idx >= end_idx
                    error("timing_pairs must have start_idx < end_idx for condition '$name'")
                end
            end
        end

        # Validate after/before constraints
        if after !== nothing && before !== nothing
            error("Cannot specify both 'after' and 'before' constraints for condition '$name'")
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


"""
    search_sequences(array, sequences)
Return indices of any of the specified sequences within an array, handling wildcards and ranges.
### Examples:
```julia
# Single sequence
idx = search_sequences([1, 2, 3, 4, 5], [[1, 2, 3]])

# Multiple sequences (OR logic)
idx = search_sequences([1, 2, 1, 3, 1, 3, 1], [[1, 2, 1], [1, 3, 1]])

# Wildcards
idx = search_sequences([1, 2, 3, 4, 5], [[1, :any, 3]])

# Ranges
idx = search_sequences([1, 2, 3, 4, 5], [[1:3], [5:5]])

# Mixed sequences
idx = search_sequences([1, 2, 3, 4, 5], [[1, 2:4, 5]])
```
"""
function search_sequences(array, sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}})
    all_indices = Int[]

    for sequence in sequences
        indices = search_single_sequence(array, sequence)
        append!(all_indices, indices)
    end

    return unique(all_indices)
end

"""
    search_single_sequence(array, sequence)
Return indices of a single sequence within an array, handling wildcards and ranges.
"""
function search_single_sequence(array, sequence::Vector{Union{Int,Symbol,UnitRange{Int}}})
    # Handle single trigger case
    if length(sequence) == 1
        if sequence[1] isa Int
            # Use search_sequence from misc.jl to find sequence starts, not all occurrences
            return eegfun.search_sequence(array, sequence[1])
        elseif sequence[1] isa UnitRange
            return search_trigger_ranges(array, [sequence[1]])
        else
            error("Single wildcard sequences not supported")
        end
    end

    # Find starting positions for the first trigger
    first_trigger = sequence[1]
    if first_trigger isa Symbol
        error("First element in sequence cannot be a wildcard")
    elseif first_trigger isa UnitRange
        # For ranges, find all possible starting positions
        idx_start_positions = search_trigger_ranges(array, [first_trigger])
    else
        # Use search_sequence from misc.jl to find sequence starts
        idx_start_positions = eegfun.search_sequence(array, first_trigger)
    end

    # sequence of values
    idx_positions = []
    for idx in idx_start_positions
        good_sequence = true
        for seq = 1:(length(sequence)-1)
            if (idx + seq) > length(array)
                good_sequence = false
                break
            end

            expected_trigger = sequence[seq+1]
            actual_trigger = array[idx+seq]

            if expected_trigger isa Int
                # Exact match required
                if actual_trigger != expected_trigger
                    good_sequence = false
                    break
                end
            elseif expected_trigger == :any
                # Wildcard - matches any value
                continue
            elseif expected_trigger isa UnitRange
                # Range - check if actual trigger is in range
                if !(actual_trigger in expected_trigger)
                    good_sequence = false
                    break
                end
            else
                error("Unsupported trigger type: $expected_trigger")
            end
        end
        if good_sequence
            push!(idx_positions, idx)
        end
    end

    return idx_positions
end

"""
    search_trigger_ranges(array, trigger_ranges)
Return indices where triggers match any of the specified ranges.
### Examples:
```julia
idx = search_trigger_ranges([1, 2, 3, 4, 5, 6], [1:3, 5:6])  # Returns indices of 1,2,3,5,6
```
"""
function search_trigger_ranges(array, trigger_ranges::Vector{UnitRange{Int}})
    idx_positions = Int[]

    for (i, trigger) in enumerate(array)
        for range in trigger_ranges
            if trigger in range
                push!(idx_positions, i)
                break  # Only add each position once
            end
        end
    end

    return idx_positions
end

"""
    search_sequence(array::AbstractVector, sequence::Int) -> Vector{Int}

Find starting indices of a sequence in an array.

# Returns
- `Vector{Int}`: Indices where sequence starts
"""
search_sequence(array::AbstractVector, sequence::Int) =
    intersect(findall(array .== sequence), findall(diff(vcat(0, array)) .>= 1))


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
    # Get all EEG channels (excluding non-numeric columns)
    all_columns = propertynames(first(dat.data))
    numeric_columns = Symbol[]

    for col in all_columns
        # Check if the column contains numeric data (excluding Bool)
        col_type = eltype(first(dat.data)[!, col])
        if col_type <: Number && col_type != Bool
            push!(numeric_columns, col)
        end
    end

    # Remove time column from averaging (keep for grouping)
    eeg_channels = setdiff(numeric_columns, [:time])

    # Concatenate all epochs
    all_epochs = reduce(vcat, dat.data)

    # Group by time and condition, then average EEG channels and count epochs
    erp = combine(
        groupby(all_epochs, [:time, :condition, :condition_name]),
        :epoch => length => :n_epochs,      # Add epoch count
        eeg_channels .=> mean .=> eeg_channels,  # Average EEG channels
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
