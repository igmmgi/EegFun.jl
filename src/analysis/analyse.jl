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
    _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}

Cleans trigger data by detecting only the onset (first occurrence) of each trigger value.
Converts sustained trigger signals into single onset events.

# Arguments
- `trigger_data::Vector{<:Integer}`: Raw trigger data vector

# Returns
- `Vector{<:Integer}`: Cleaned trigger data with only onset events

# Examples
```julia
# Input:  [0, 1, 1, 0, 0, 2, 2, 2, 0, 0]
# Output: [0, 1, 0, 0, 0, 2, 0, 0, 0, 0]
```
"""
function _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}

    # find where triggers change (onset detection)
    trigger_changes = diff(vcat(0, trigger_data))
    onset_indices = trigger_changes .> 0

    # insert onsets
    cleaned = zeros(Int, length(trigger_data))
    cleaned[onset_indices] = trigger_data[onset_indices]

    return cleaned

end


"""
    create_eeg_dataframe(data::BioSemiBDF.BioSemiData)::DataFrame

Creates a DataFrame containing EEG data from a BioSemiBDF object.

# Arguments
- `data::BioSemiBDF.BioSemiData`: The BioSemi data structure containing time, triggers, and channel data.

# Returns
A DataFrame with columns: `time`, `sample`, `triggers` (cleaned), and channel data from the BioSemiBDF.

# Note
The trigger data is automatically cleaned to detect only onset events, converting sustained trigger 
signals into single onset events. For example: [0, 1, 1, 0, 0, 2, 2, 2, 0, 0] becomes [0, 1, 0, 0, 0, 2, 0, 0, 0, 0].
"""
function create_eeg_dataframe(data::BioSemiBDF.BioSemiData)::DataFrame
    @info "create_eeg_dataframe: Creating EEG DataFrame"
    df = hcat(
        DataFrame(time = data.time, sample = 1:length(data.time), triggers = _clean_triggers(data.triggers.raw)),
        DataFrame(Float64.(data.data), Symbol.(data.header.channel_labels[1:(end-1)])),  # assumes last channel is trigger
    )
    return df
end


"""
    create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout::DataFrame)::ContinuousData

Creates a ContinuousData object from a BioSemiBDF data structure and a layout DataFrame.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemi data structure containing EEG data.
- `layout::DataFrame`: The DataFrame containing layout information.

# Returns
A ContinuousData object containing the EEG data and layout information.

"""
function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout::Layout)::ContinuousData
    return ContinuousData(create_eeg_dataframe(dat), layout, dat.header.sample_rate[1], AnalysisInfo())
end





# channel_summary
function _channel_summary_impl(
    data::DataFrame,
    sample_selection::Vector{Int},
    channel_selection::Vector{Symbol},
)::DataFrame
    selected_data = @view data[sample_selection, channel_selection]

    # Get base statistics from describe
    stats = describe(selected_data, :min, :max, :std)

    # Add our custom columns
    stats.range = stats.max .- stats.min
    stats.var = var.(eachcol(selected_data))
    stats.zvar = zscore(stats.var)

    # Rename the variable column to channel
    rename!(stats, :variable => :channel)

    return stats
end

# Helper function predicates for easier channel filtering
channels() = x -> fill(true, length(x))  # Default: select all channels given
channels(channel_names::Vector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_name::Symbol) = x -> x .== channel_name
channels(channel_numbers::Union{Vector{Int},UnitRange}) = x -> [i in channel_numbers for i = 1:length(x)]
channels_not(channel_names::Vector{Symbol}) = x -> .!(x .∈ Ref(channel_names))
channels_not(channel_name::Symbol) = x -> .!(x .== channel_name)
channels_not(channel_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in channel_numbers for i = 1:length(x)])

# Helper function predicates for easier component filtering
components() = x -> fill(true, length(x))  # Default: select all components given
components(component_numbers::Union{Vector{Int},UnitRange}) = x -> [i in component_numbers for i = 1:length(x)]
components(component_number::Int) = x -> x .== component_number
components_not(component_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in component_numbers for i = 1:length(x)])
components_not(component_number::Int) = x -> .!(x .== component_number)

# Helper function predicates for easier sample filtering
samples() = x -> fill(true, nrow(x))
samples(column::Symbol) = x -> x[!, column]
samples_or(columns::Vector{Symbol}) = x -> any(x[!, col] for col in columns)
samples_and(columns::Vector{Symbol}) = x -> all(x[!, col] for col in columns)
samples_not(column::Symbol) = x -> .!(x[!, column])
samples_or_not(columns::Vector{Symbol}) = x -> .!(any(x[!, col] for col in columns))
samples_and_not(columns::Vector{Symbol}) = x -> .!(all(x[!, col] for col in columns))

# Helper function predicates for easier epoch filtering
epochs() = x -> fill(true, length(x))  # Default: select all epochs given
epochs(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> [i in epoch_numbers for i in x]
epochs(epoch_number::Int) = x -> x .== epoch_number
epochs_not(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in epoch_numbers for i in x])
epochs_not(epoch_number::Int) = x -> .!(x .== epoch_number)

# Helper to select channels/columns based on a predicate (+ which to include)
function get_selected_channels(dat, channel_selection::Function; include_meta::Bool = true, include_extra::Bool = true)

    # Columns/channels in dataframe to include
    metadata_cols = include_meta ? meta_labels(dat) : Symbol[]
    selectable_cols = include_extra ? vcat(channel_labels(dat), extra_labels(dat)) : channel_labels(dat)

    # Apply channel selection to non-metadata columns
    selected_cols = selectable_cols[channel_selection(selectable_cols)]

    # Return metadata + selected channels
    return vcat(metadata_cols, selected_cols)
end


# Helper to select components based on a predicate
function get_selected_components(ica_result::InfoIca, component_selection::Function)
    # Get all component indices (1 to n_components)
    all_components = 1:length(ica_result.ica_label)
    return all_components[component_selection(all_components)]
end

# Helper to select samples based on a predicate
function get_selected_samples(dat, sample_selection::Function)
    return findall(sample_selection(dat.data))
end

# Helper to select samples from a DataFrame
function get_selected_samples(dat::DataFrame, sample_selection::Function)
    return findall(sample_selection(dat))
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
    channel_summary(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_extr ::Bool = false)::DataFrame

Computes summary statistics for EEG channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_extra::Bool`: Whether to include additional channels (default: false).

# Returns
A DataFrame containing summary statistics for each channel.

# Examples

## Basic Usage
```julia
# Channel summary for layout channels only (default)
summary = channel_summary(dat)

# Channel summary for specific layout channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Channel summary excluding reference channels from layout
summary = channel_summary(dat, channel_selection = channels_not([:M1, :M2]))
```

## Including Additional Channels
```julia
# Channel summary for additional channels (EOG, extreme value flags, etc.)
# The function automatically detects when you specify additional channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Channel Selection
```julia
# Summary for specific channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Summary excluding reference channels
summary = channel_summary(dat, channel_selection = channels_not([:M1, :M2]))

# Summary for frontal channels only (channels 1-10)
summary = channel_summary(dat, channel_selection = channels(1:10))

# Summary for channels starting with "F" (frontal)
summary = channel_summary(dat, channel_selection = x -> startswith.(string.(x), "F"))
```

## Sample Selection
```julia
# Exclude extreme values
summary = channel_summary(dat, sample_selection = samples_not(:is_extreme_value_100))

# Exclude multiple types of bad samples
summary = channel_summary(dat, sample_selection = samples_or_not([:is_extreme_value_100, :is_vEOG, :is_hEOG]))

# Only include samples within epoch windows
summary = channel_summary(dat, sample_selection = samples(:epoch_window))

# Include samples that are both in epoch window AND not extreme
summary = channel_summary(dat, sample_selection = samples_and([:epoch_window, samples_not(:is_extreme_value_100)]))
```

## Combined Selection
```julia
# Exclude reference channels and extreme values
summary = channel_summary(dat, 
    channel_selection = channels_not([:M1, :M2]),
    sample_selection = samples_not(:is_extreme_value_100)
)

# Only frontal channels, exclude bad samples
summary = channel_summary(dat, 
    channel_selection = channels(1:10),
    sample_selection = samples_or_not([:is_extreme_value_100, :is_vEOG])
)

# Complex filtering: frontal channels, good samples, within epochs
summary = channel_summary(dat, 
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
    sample_selection = samples_and([
        :epoch_window, 
        samples_not(:is_extreme_value_100),
        samples_not(:is_vEOG),
        samples_not(:is_hEOG)
    ])
)
```

## Additional Channels (not in layout)
```julia
# Include derived channels like EOG
# The function automatically switches to all available channels when needed
summary = channel_summary(dat, channel_selection = channels([:vEOG, :hEOG]))

# Mix layout channels and additional channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```
"""
function channel_summary(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    selected_channels =
        get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    return _channel_summary_impl(dat.data, selected_samples, selected_channels)
end

"""
    convert(dat::MultiDataFrameEeg, epoch_idx::Int) -> SingleDataFrameEeg

Convert a single epoch from MultiDataFrameEeg to SingleDataFrameEeg.

# Arguments
- `dat::MultiDataFrameEeg`: The multi-DataFrame EEG data
- `epoch_idx::Int`: Index of the epoch to convert (1-based)

# Returns
- `SingleDataFrameEeg`: Single DataFrame containing only the specified epoch

# Examples
```julia
# Convert epoch 3 to single DataFrame
single_dat = convert(dat, 3)

# Now you can use single DataFrame functions
summary = channel_summary(single_dat)
```
"""
function convert(dat::MultiDataFrameEeg, epoch_idx::Int)::SingleDataFrameEeg
    # Validate epoch index
    if epoch_idx < 1 || epoch_idx > length(dat.data)
        @minimal_error "Epoch index $epoch_idx out of range (1:$(length(dat.data)))"
    end
    return ContinuousData(dat.data[epoch_idx], dat.layout, dat.sample_rate, dat.analysis_info)
end

function channel_summary(
    dat::MultiDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    # Process each epoch and collect results
    results = DataFrame[]

    for epoch_df in dat.data
        # Get the original epoch number from the data
        original_epoch_number = epoch_df.epoch[1]  # All rows in an epoch have the same epoch number

        # Create ContinuousData from this epoch DataFrame
        single_dat = ContinuousData(epoch_df, dat.layout, dat.sample_rate, dat.analysis_info)

        # Get summary for this epoch
        epoch_summary = channel_summary(
            single_dat;
            sample_selection = sample_selection,
            channel_selection = channel_selection,
            include_meta = include_meta,
            include_extra = include_extra,
        )

        # Add epoch column as first column with original epoch number
        insertcols!(epoch_summary, 1, :epoch => fill(original_epoch_number, nrow(epoch_summary)))

        push!(results, epoch_summary)
    end

    # Combine all results
    return vcat(results...)
end



"""
    correlation_matrix(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_extra::Bool = false)::Matrix{Float64}

Calculates the correlation matrix for the EEG data.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_extra::Bool`: Whether to include additional channels (default: false).

# Returns
A matrix containing the correlation values between the specified channels.

# Examples

## Basic Usage
```julia
# Correlation matrix for layout channels only (default)
corr_matrix = correlation_matrix(dat)

# Correlation matrix for specific layout channels
corr_matrix = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Correlation matrix excluding reference channels from layout
corr_matrix = correlation_matrix(dat, channel_selection = channels_not([:M1, :M2]))
```

## Including Additional Channels
```julia
# Correlation matrix for additional channels (EOG, extreme value flags, etc.)
corr_matrix = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Sample Selection
```julia
# Correlation matrix only for good samples
corr_matrix = correlation_matrix(dat, sample_selection = samples_not(:is_extreme_value_100))

# Correlation matrix only within epoch windows
corr_matrix = correlation_matrix(dat, sample_selection = samples(:epoch_window))

# Correlation matrix for good samples within epochs
corr_matrix = correlation_matrix(dat, 
    sample_selection = samples_and([
        :epoch_window,
        samples_not(:is_extreme_value_100)
    ])
)
```

## Channel Selection
```julia
# Frontal channels only
corr_matrix = correlation_matrix(dat, channel_selection = channels(1:10))

# Channels starting with "F" (frontal)
corr_matrix = correlation_matrix(dat, channel_selection = x -> startswith.(string.(x), "F"))

# Parietal channels only
corr_matrix = correlation_matrix(dat, channel_selection = x -> startswith.(string.(x), "P"))

# Mix of layout and additional channels
corr_matrix = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Combined Selection
```julia
# Exclude reference channels and bad samples
corr_matrix = correlation_matrix(dat, 
    channel_selection = channels_not([:M1, :M2]),
    sample_selection = samples_not(:is_extreme_value_100)
)

# Only frontal channels, good samples
corr_matrix = correlation_matrix(dat, 
    channel_selection = channels(1:10),
    sample_selection = samples_and([
        :epoch_window, 
        samples_not(:is_extreme_value_100)
    ])
)

# Complex filtering: frontal channels, good samples, within epochs
corr_matrix = correlation_matrix(dat, 
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
    sample_selection = samples_and([
        :epoch_window, 
        samples_not(:is_extreme_value_100),
        samples_not(:is_vEOG),
        samples_not(:is_hEOG)
    ])
)
```

## Quality Control Applications
```julia
# 1. Detect extreme values
is_extreme_value!(dat, 100, channel_out = :is_extreme_100)

# 2. Get correlation matrix for good data
good_corr = correlation_matrix(dat, sample_selection = samples_not(:is_extreme_100))
```
"""
function correlation_matrix(
    dat::ContinuousData;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
)::DataFrame
    selected_channels = get_selected_channels(dat, channel_selection; include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    return _correlation_matrix(dat.data, selected_samples, selected_channels)
end

# Internal function for plain DataFrames with explicit channel specification
function _correlation_matrix(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol},
)::DataFrame
    # Select the specified channels and samples
    selected_data = select(dat, selected_channels)[selected_samples, :]
    # Compute correlation matrix
    df = DataFrame(cor(Matrix(selected_data)), selected_channels)
    insertcols!(df, 1, :row => selected_channels)
    return df
end
