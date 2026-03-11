"""
    condition_parse_epoch(config::Dict) -> Vector{EpochCondition}

Parse epoch condition configurations from a TOML-based configuration dictionary.

This function is primarily used by the built-in preprocessing pipelines to parse epoch 
conditions from TOML files. If you're creating a custom pipeline and want to use the 
same TOML configuration format, you can use this function. Otherwise, you can create 
`EpochCondition` objects directly.

# TOML Format
See the pipeline documentation for the expected TOML structure.

# Example
```julia
config = TOML.parsefile("epoch_conditions.toml")
conditions = condition_parse_epoch(config)
```
"""
function condition_parse_epoch(config::Dict)
    epochs_section = get(config, "epochs", Dict())

    conditions = EpochCondition[]
    condition_configs = get(epochs_section, "conditions", [])

    for condition_config in condition_configs
        name = condition_config["name"]

        # Parse trigger sequences (unified approach)
        trigger_sequences_raw = get(condition_config, "trigger_sequences", nothing)
        if isnothing(trigger_sequences_raw)
            @minimal_error("trigger_sequences must be specified for condition '$name'")
        end

        # Parse trigger sequences - convert from TOML arrays to proper types
        trigger_sequences = Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}()
        for sequence in trigger_sequences_raw
            # Normalize: wrap non-vectors (e.g., UnitRange) in a vector
            seq_vec = sequence isa Vector ? sequence : [sequence]
            push!(trigger_sequences, collect(Union{Int,Symbol,UnitRange{Int}}, seq_vec))
        end

        # Parse reference index (single index) - default to 1
        reference_index = get(condition_config, "reference_index", 1)

        # Parse timing pairs (optional - if not specified, no timing constraints)
        timing_pairs_raw = get(condition_config, "timing_pairs", nothing)
        if isnothing(timing_pairs_raw)
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
            if isnothing(min_interval) || isnothing(max_interval)
                @minimal_error("Both min_interval and max_interval must be specified when timing_pairs is specified for condition '$name'",)
            end
        end

        # Parse after/before constraints (optional)
        after = get(condition_config, "after", nothing)
        before = get(condition_config, "before", nothing)

        # Validation - cache sequence length for efficiency
        seq_length = length(trigger_sequences[1])
        if reference_index < 1 || reference_index > seq_length
            @minimal_error("reference_index must be between 1 and $seq_length for condition '$name'")
        end

        # Only validate timing constraints if they're specified
        if !isnothing(timing_pairs)
            if min_interval >= max_interval
                @minimal_error("min_interval must be < max_interval for condition '$name'")
            end

            # Validate timing pairs
            for (start_idx, end_idx) in timing_pairs
                if start_idx < 1 || start_idx > seq_length || end_idx < 1 || end_idx > seq_length
                    @minimal_error("timing_pairs contains invalid indices for sequence of length $seq_length in condition '$name'",)
                end
                if start_idx >= end_idx
                    @minimal_error("timing_pairs must have start_idx < end_idx for condition '$name'")
                end
            end
        end

        # Validate after/before constraints
        if !isnothing(after) && !isnothing(before)
            @minimal_error("Cannot specify both 'after' and 'before' constraints for condition '$name'")
        end

        push!(conditions, EpochCondition(name, trigger_sequences, reference_index, timing_pairs, min_interval, max_interval, after, before))
    end

    return conditions
end

# Helper function to validate epoch interval parameters
function _validate_epoch_interval_params(dat::ContinuousData, epoch_interval::Vector{<:Real})
    length(epoch_interval) == 2 || @minimal_error "Epoch interval must have exactly 2 elements"
    epoch_interval[1] <= epoch_interval[2] || @minimal_error "Epoch interval start must be less than or equal to end"
    hasproperty(dat.data, :trigger) || @minimal_error "Data must have a trigger column"
    hasproperty(dat.data, :time) || @minimal_error "Data must have a time column"
    !isempty(dat.data.time) || @minimal_error "Time column cannot be empty"
    !isempty(dat.data.trigger) || @minimal_error "Triggers column cannot be empty"
    issorted(dat.data.time) || @minimal_error "Time column must be sorted in ascending order"
    return nothing
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
    mark_epoch_intervals!(dat::ContinuousData, triggers_of_interest::Vector{Int}, time_window::Vector{<:Real}; channel_out::Symbol = :epoch_interval)
    mark_epoch_intervals!(dat::ContinuousData, channel_in::Symbol, time_window::Vector{<:Real}; channel_out = Symbol(string(channel_in)*"_window"))
    mark_epoch_intervals!(dat::ContinuousData, channels_in::Vector{Symbol}, time_window::Vector{<:Real}; channel_out::Symbol = :event_window)
    mark_epoch_intervals!(dat::ContinuousData, epoch_conditions::Vector{EpochCondition}, time_window::Vector{<:Real}; channel_out::Symbol = :epoch_interval)

Mark samples within a time window in-place by adding a boolean column to the data.

- **Trigger form**: marks windows around samples matching any trigger in `triggers_of_interest`
- **Column form**: marks windows around each `true` sample in a boolean column (e.g., `:is_vEOG`)
- **Multi-column form**: unions `true` samples across multiple boolean columns and marks windows around all of them
- **EpochCondition form**: marks windows around trigger sequences defined by `epoch_conditions`

# Examples
```julia
mark_epoch_intervals!(dat, [1, 3, 5], [-1.0, 2.0])
mark_epoch_intervals!(dat, :is_vEOG, [-0.1, 0.3])
mark_epoch_intervals!(dat, [:is_vEOG, :is_hEOG], [-0.05, 0.4]; channel_out = :eog_window)
mark_epoch_intervals!(dat, [condition1, condition2], [-1.0, 2.0])
```
"""
function mark_epoch_intervals!(
    dat::ContinuousData,
    triggers_of_interest::Vector{Int},
    time_window::Vector{<:Real};
    channel_out::Symbol = :epoch_interval,
)
    # Input validation
    _validate_epoch_interval_params(dat, time_window)

    # Initialize result vector with false 
    dat.data[!, channel_out] .= false

    # Collect all relevant trigger indices
    all_reference_indices = Int[]

    for trigger in triggers_of_interest
        trigger_indices = findall(dat.data.trigger .== trigger)
        if isempty(trigger_indices)
            @minimal_warning "Trigger $trigger not found in data"
            continue
        end
        append!(all_reference_indices, trigger_indices)
    end

    # Mark windows around all collected indices
    n_marked = _mark_windows_at_indices!(dat, all_reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(all_reference_indices)) trigger locations"

    return dat
end



function mark_epoch_intervals!(
    dat::ContinuousData,
    channel_in::Symbol,
    time_window::Vector{<:Real};
    channel_out::Symbol = Symbol(string(channel_in) * "_window"),
)
    # Input validation
    _validate_epoch_interval_params(dat, time_window)
    channel_in ∉ propertynames(dat.data) && @minimal_error("Column $(channel_in) not found in data")
    eltype(dat.data[!, channel_in]) <: Bool || @minimal_error("Column $(channel_in) must be a Bool column")

    # Initialize result column
    dat.data[!, channel_out] .= false

    # Find all true indices in the boolean column
    reference_indices = findall(dat.data[!, channel_in])

    if isempty(reference_indices)
        @minimal_warning "No true values found in column $(channel_in)"
        return dat
    end

    # Mark windows around all found indices
    n_marked = _mark_windows_at_indices!(dat, reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(reference_indices)) $(channel_in) events"

    return dat
end


function mark_epoch_intervals!(
    dat::ContinuousData,
    channels_in::Vector{Symbol},
    time_window::Vector{<:Real};
    channel_out::Symbol = :event_window,
)
    # Input validation
    _validate_epoch_interval_params(dat, time_window)
    for channel_in in channels_in
        channel_in ∉ propertynames(dat.data) && @minimal_error("Column $(channel_in) not found in data")
        eltype(dat.data[!, channel_in]) <: Bool || @minimal_error("Column $(channel_in) must be a Bool column")
    end

    # Initialize result column
    dat.data[!, channel_out] .= false

    # Union true indices across all boolean columns
    reference_indices = Int[]
    for channel_in in channels_in
        append!(reference_indices, findall(dat.data[!, channel_in]))
    end
    sort!(unique!(reference_indices))

    if isempty(reference_indices)
        @minimal_warning "No true values found in any of the columns $(channels_in)"
        return dat
    end

    # Mark windows around all collected indices
    n_marked = _mark_windows_at_indices!(dat, reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(reference_indices)) events from $(length(channels_in)) columns"

    return dat
end

"""
    _filter_matches(matches, condition, triggers, time_data) -> Vector{Vector{Int}}

Filter sequence matches by after/before position constraints and timing pair constraints.
Used by both `extract_epochs` and `mark_epoch_intervals!`.
"""
function _filter_matches(matches::Vector{Vector{Int}}, condition::EpochCondition, triggers, time_data)
    # Apply after/before filtering
    if !isnothing(condition.after) || !isnothing(condition.before)
        matches = filter(matches) do m
            if !isnothing(condition.after)
                any(triggers[1:(m[1]-1)] .== condition.after) || return false
            end
            if !isnothing(condition.before)
                any(triggers[(m[end]+1):end] .== condition.before) || return false
            end
            return true
        end
    end

    # Apply timing constraints
    if !isnothing(condition.timing_pairs) && !isnothing(condition.min_interval) && !isnothing(condition.max_interval)
        matches = filter(matches) do m
            for (start_idx, end_idx) in condition.timing_pairs
                interval = time_data[m[end_idx]] - time_data[m[start_idx]]
                (interval < condition.min_interval || interval > condition.max_interval) && return false
            end
            return true
        end
    end

    return matches
end



function mark_epoch_intervals!(
    dat::ContinuousData,
    epoch_conditions::Vector{EpochCondition},
    time_window::Vector{<:Real};
    channel_out::Symbol = :epoch_interval,
)
    # Input validation
    _validate_epoch_interval_params(dat, time_window)

    # Initialize result vector with false
    dat.data[!, channel_out] .= false

    # Collect all reference indices from all conditions
    all_reference_indices = Int[]

    # For each epoch condition
    for condition in epoch_conditions

        # Find all occurrences of the trigger sequences — each match is [pos1, pos2, ...]
        matches = search_sequence(dat.data.trigger, condition.trigger_sequences)
        if isempty(matches)
            @minimal_warning "No triggers found for condition '$(condition.name)'"
            continue
        end

        # Apply position and timing constraints
        matches = _filter_matches(matches, condition, dat.data.trigger, dat.data.time)
        if isempty(matches)
            after_msg = !isnothing(condition.after) ? " after trigger $(condition.after)" : ""
            before_msg = !isnothing(condition.before) ? " before trigger $(condition.before)" : ""
            @minimal_warning "No trigger sequences found that meet constraints$(after_msg)$(before_msg) for condition '$(condition.name)'"
            continue
        end

        # Convert matches to reference positions and collect them
        for m in matches
            if condition.reference_index > length(m)
                @minimal_warning "Reference index $(condition.reference_index) exceeds sequence length $(length(m)) for condition '$(condition.name)'"
                continue
            end
            reference_idx = m[condition.reference_index]
            if reference_idx <= length(dat.data.trigger)
                push!(all_reference_indices, reference_idx)
            else
                @minimal_warning "Reference index $(condition.reference_index) for condition '$(condition.name)' is out of bounds"
            end
        end
    end

    # Mark windows around all collected reference indices
    n_marked = _mark_windows_at_indices!(dat, all_reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(all_reference_indices)) reference points"

    return dat
end






"""
    extract_epochs(dat::ContinuousData, condition::Int, epoch_condition::EpochCondition, epoch_interval::Tuple{Real,Real})

Extract epochs based on a single EpochCondition object, including timing validation, after/before filtering, 
trigger ranges, wildcard sequences, and multiple sequences.

# Arguments
- `dat::ContinuousData`: The continuous EEG data
- `condition::Int`: Condition number to assign to epochs
- `epoch_condition::EpochCondition`: EpochCondition object defining trigger sequence and timing constraints
- `epoch_interval::Tuple{Real,Real}`: Time window relative to reference point (start_time, end_time) in seconds

# Returns
- `EpochData`: The extracted epochs

# Examples
```julia
# Single sequence
condition = EpochCondition(name="single", trigger_sequence=[1, 2, 3], reference_index=2)
epochs = extract_epochs(dat, 1, condition, (-0.2, 1.0))  # -200ms to 1000ms

# Multiple sequences (OR logic)
condition = EpochCondition(
    name="multiple", 
    trigger_sequences=[[1, 2, 1], [1, 3, 1]],  # Match either [1,2,1] OR [1,3,1]
    reference_index=2
)
epochs = extract_epochs(dat, 1, condition, (-0.5, 2.0))

# Wildcard sequence
condition = EpochCondition(
    name="wildcard", 
    trigger_sequence=[1, :any, 3],  # :any matches any trigger value
    reference_index=2
)
epochs = extract_epochs(dat, 1, condition, (-0.2, 1.0))

# Trigger ranges
condition = EpochCondition(
    name="ranges", 
    trigger_ranges=[1:5, 10:15],  # Match triggers 1-5 or 10-15
    reference_index=1
)
epochs = extract_epochs(dat, 1, condition, (-0.2, 1.0))
```
"""
function extract_epochs(dat::ContinuousData, condition::Int, epoch_condition::EpochCondition, epoch_interval::Tuple{Real,Real})
    # Extract start and end times from tuple
    start_time, end_time = epoch_interval

    # Find matched positions for trigger sequences — each match is [pos1, pos2, ...]
    matches = search_sequence(dat.data.trigger, epoch_condition.trigger_sequences)
    isempty(matches) && @minimal_error("None of the trigger sequences $(epoch_condition.trigger_sequences) found!")

    # Apply position and timing constraints
    matches = _filter_matches(matches, epoch_condition, dat.data.trigger, dat.data.time)
    if isempty(matches)
        after_msg = !isnothing(epoch_condition.after) ? " after trigger $(epoch_condition.after)" : ""
        before_msg = !isnothing(epoch_condition.before) ? " before trigger $(epoch_condition.before)" : ""
        @minimal_error("No trigger sequences found that meet constraints$(after_msg)$(before_msg) for condition '$(epoch_condition.name)'",)
    end

    # Extract t==0 positions using reference_index (filter out matches where reference_index exceeds match length)
    matches = filter(m -> epoch_condition.reference_index <= length(m), matches)
    if isempty(matches)
        @minimal_error("Reference index $(epoch_condition.reference_index) exceeds sequence length for condition '$(epoch_condition.name)'")
    end
    zero_idx = [m[epoch_condition.reference_index] for m in matches]

    # find number of samples pre/post epoch t = 0 position
    n_pre, n_post = _find_idx_start_end(dat.data.time, abs(start_time), abs(end_time))
    pre_idx = zero_idx .- n_pre .+ 1
    post_idx = zero_idx .+ n_post .- 1

    # Extract and create array of dataframes with bounds checking
    epochs = DataFrame[]

    for (epoch, (pre, zero, post)) in enumerate(zip(pre_idx, zero_idx, post_idx))
        # Bounds checking to prevent out-of-bounds errors
        if pre < 1 || post > nrow(dat.data)
            @minimal_warning "Epoch $epoch extends beyond data bounds (pre=$pre, post=$post, data_length=$(nrow(dat.data))) - skipping"
            continue
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

    return EpochData(dat.file, condition, epoch_condition.name, epochs, dat.layout, dat.sample_rate, dat.analysis_info)
end

function extract_epochs(dat::ContinuousData, epoch_conditions::Vector{EpochCondition}, epoch_interval::Tuple{Real,Real})
    epochs = EpochData[]
    for (idx, epoch_condition) in enumerate(epoch_conditions)
        push!(epochs, extract_epochs(dat, idx, epoch_condition, epoch_interval))
    end
    return epochs
end

function extract_epochs(dat::ContinuousData, epoch_conditions::EpochCondition, epoch_interval::Tuple{Real,Real})
    return extract_epochs(dat, [epoch_conditions], epoch_interval)
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
    isempty(dat.data) && @minimal_error("Cannot average empty EpochData")

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
    isempty(eeg_channels) && @minimal_error("No EEG channels found to average")

    # Average epochs directly by row index (all epochs have same length and time values)
    try
        first_epoch = first(dat.data)
        n_timepoints = nrow(first_epoch)

        # Verify all epochs have the same length
        for (idx, epoch) in enumerate(dat.data)
            if nrow(epoch) != n_timepoints
                @minimal_error("Epoch $idx has $(nrow(epoch)) timepoints, expected $n_timepoints. All epochs must have the same length.")
            end
        end

        # Create result DataFrame with metadata columns from first epoch
        erp = DataFrame()

        # Copy metadata columns (only :time should be kept for an average)
        metadata_cols = filter(c -> c == :time, meta_labels(dat))
        for col in metadata_cols
            if hasproperty(first_epoch, col)
                erp[!, col] = first_epoch[!, col]
            end
        end

        # Average EEG channels across epochs by row index
        for ch in eeg_channels
            # Stack all epochs for this channel: n_epochs × n_timepoints
            channel_matrix = hcat([epoch[!, ch] for epoch in dat.data]...)
            # Average across epochs (mean of each row = each timepoint)
            erp[!, ch] = vec(mean(channel_matrix, dims = 2))
        end

        # Count epochs
        n_epochs = length(dat.data)

        return ErpData(dat.file, dat.condition, dat.condition_name, erp, dat.layout, dat.sample_rate, dat.analysis_info, n_epochs)
    catch e
        @minimal_error("Failed to average epochs: $(e)")
    end
end


average_epochs(dat::Vector{EpochData}) = average_epochs.(dat)



"""
    reject_epochs!(dat::EpochData, info::EpochRejectionInfo) -> EpochData

Remove epochs identified in `info` from `dat` in-place.

# Arguments
- `dat::EpochData`: Epoched EEG data to modify
- `info::EpochRejectionInfo`: Rejection information from `detect_bad_epochs_automatic` etc.

# Examples
```julia
info = detect_bad_epochs_automatic(epochs, z_criterion = 2.0)
reject_epochs!(epochs, info)
```
"""
function reject_epochs!(dat::EpochData, info::EpochRejectionInfo)::EpochData

    if isempty(info.rejected)
        @info "No epochs to reject"
        return dat
    end

    n_epochs = length(dat.data)
    rejected_indices = unique([r.epoch for r in info.rejected])
    epochs_to_keep = setdiff(1:n_epochs, rejected_indices)
    dat.data = dat.data[epochs_to_keep]

    @info "Condition $(dat.condition) ($(dat.condition_name)) - Rejected $(length(rejected_indices)) of $(info.info.n) epochs."

    return dat
end


"""
    reject_epochs(dat::EpochData, info::EpochRejectionInfo) -> EpochData
    reject_epochs(dat::Vector{EpochData}, info::Vector{EpochRejectionInfo}) -> Vector{EpochData}
    reject_epochs(dat::EpochData, bad_column::Symbol) -> EpochData
    reject_epochs(dat::EpochData, bad_columns::Vector{Symbol}) -> EpochData
    reject_epochs(dat::Vector{EpochData}, bad_column::Symbol) -> Vector{EpochData}
    reject_epochs(dat::Vector{EpochData}, bad_columns::Vector{Symbol}) -> Vector{EpochData}

Non-mutating epoch rejection. Returns new `EpochData` with bad epochs removed.

Accepts rejection info from `detect_bad_epochs_automatic`, or one or more boolean column
names — epochs where any sample is `true` in those columns are removed.

# Examples
```julia
# From detection info
info = detect_bad_epochs_automatic(epochs, z_criterion = 2.0)
cleaned = reject_epochs(epochs, info)

# From boolean artifact column
cleaned = reject_epochs(epochs, :is_vEOG)
cleaned = reject_epochs(epochs, [:is_vEOG, :is_extreme_value_100])
```
"""
function reject_epochs(dat::EpochData, info::EpochRejectionInfo)::EpochData
    dat_copy = copy(dat)
    reject_epochs!(dat_copy, info)
    return dat_copy
end

reject_epochs(dat::Vector{EpochData}, info::Vector{EpochRejectionInfo})::Vector{EpochData} = reject_epochs.(dat, info)


function reject_epochs(dat::EpochData, bad_columns::Vector{Symbol})
    # Input validation
    isempty(dat.data) && @minimal_error("Cannot remove bad epochs from empty EpochData")
    isempty(bad_columns) && return dat  # No columns to check

    # Validate that all bad_columns exist in the data
    first_epoch = first(dat.data)
    for col in bad_columns
        if !hasproperty(first_epoch, col)
            @minimal_error("Column '$col' not found in epoch data")
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

reject_epochs(dat::EpochData, bad_column::Symbol) = reject_epochs(dat, [bad_column])

reject_epochs(dat::Vector{EpochData}, bad_column::Symbol) = reject_epochs.(dat, bad_column)
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
function epochs_table(epochs_original::Vector{EpochData}, epochs_cleaned::Vector{EpochData}; print_table::Bool = true, kwargs...)
    length(epochs_original) != length(epochs_cleaned) && throw(ArgumentError("epochs_original and epochs_cleaned must have same length"))

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
    df = epochs_table(epochs...; print_table = print_table)
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
    !contains(pattern, "epochs") && return "average_epochs only works with epoch data. File pattern must contain 'epochs', got: '$pattern'"
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
function _process_average_file(filepath::String, output_path::String, condition_selection::Function)
    filename = basename(filepath)

    # Read data
    epochs_data = read_data(filepath)
    if isnothing(epochs_data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of EpochData)
    if !(epochs_data isa Vector{<:EpochData})
        return BatchResult(false, filename, "Invalid data type: expected Vector{EpochData}")
    end

    # Select conditions
    epochs_data = _condition_select(epochs_data, condition_selection)

    # Average epochs for each condition
    erps_data = average_epochs.(epochs_data)

    # Save (always use "data" as variable name since read_data finds by type)
    jldsave(output_path; data = erps_data)

    return BatchResult(true, filename, "Averaged $(length(erps_data)) condition(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#
function average_epochs(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "average_epochs.log"
    setup_global_logging(log_file)

    try
        @info "Batch epoch averaging started at $(now())"
        @info "average_epochs"
        @log_call "average_epochs"

        # Validation (early return on error)
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        error_msg = _validate_epochs_pattern(file_pattern)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_average_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Create processing function with captured parameters
        # Transform output filenames: replace "epochs" with "erps"
        process_fn = (input_path, output_path) -> begin
            # Transform filename: replace "epochs" with "erps" in the output filename
            output_file = basename(output_path)
            transformed_file = replace(output_file, "epochs" => "erps")
            transformed_output_path = joinpath(dirname(output_path), transformed_file)
            _process_average_file(input_path, transformed_output_path, condition_selection)
        end

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Averaging")

        # Log summary
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
