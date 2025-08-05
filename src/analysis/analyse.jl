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


"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Float64, channel_in::Symbol, channel_out::Symbol)

Detects EOG (electrooculogram) onsets in the EEG data based on a specified criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Real`: The threshold for detecting EOG onsets.
- `channel_in::Symbol`: The channel from which to detect EOG onsets.
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event.

# Returns
Nothing. The function modifies the input data in place.

# Examples

## Basic Usage
```julia
# Detect vertical EOG onsets
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)

# Detect horizontal EOG onsets
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# Combine for any EOG artifact
combine_boolean_columns!(dat, [:is_vEOG, :is_hEOG], :or, output_column = :is_any_EOG)

# Create quality flags
combine_boolean_columns!(dat, [:is_extreme_value, :is_any_EOG], :nor, output_column = :is_good_data)

# Complex quality control (good samples = not extreme AND not any EOG)
combine_boolean_columns!(dat, [:is_extreme_value, :is_vEOG, :is_hEOG], :nor, output_column = :is_clean_data)

## Multiple EOG Channels
```julia
# Detect both vertical and horizontal EOG
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# Combine for any EOG artifact
# dat.data[!, :is_any_EOG] = dat.data[!, :is_vEOG] .| dat.data[!, :is_hEOG]
```
"""
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)
    @info "detect_eog_onsets!: Detecting EOG onsets in channel $(channel_in) with stepsize criterion $(criterion)"
    step_size = div(dat.sample_rate, 20)
    eog_signal = dat.data[1:step_size:end, channel_in]
    eog_diff = diff(eog_signal)
    eog_idx = findall(x -> abs(x) >= criterion, eog_diff)
    eog_idx = [idx for (i, idx) in enumerate(eog_idx) if i == 1 || (idx - eog_idx[i-1] > 2)] * step_size
    dat.data[!, channel_out] .= false
    dat.data[eog_idx, channel_out] .= true
    return nothing
end

# Internal function for plain DataFrames with explicit channel specification
function _is_extreme_value(
    dat::DataFrame,
    criterion::Number,
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{Bool}
    # Initialize result vector with false for all samples
    result = fill(false, nrow(dat))

    # Check for extreme values only in selected samples
    if !isempty(selected_samples)
        data_subset = select(dat[selected_samples, :], selected_channels)
        extreme_mask = any(x -> abs.(x) >= criterion, Matrix(data_subset), dims = 2)[:]
        result[selected_samples] = extreme_mask
    end

    return result
end

# Internal function for plain DataFrames with explicit channel specification
function _is_extreme_value!(
    dat::DataFrame,
    criterion::Number,
    selected_channels::Vector{Symbol};
    channel_out::Symbol = :is_extreme_value,
)
    dat[!, channel_out] .= any(x -> abs.(x) >= criterion, Matrix(select(dat, selected_channels)), dims = 2)[:]
end

"""
    is_extreme_value(dat::ContinuousData, criterion::Number; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false)::Vector{Bool}

Checks if any values in the specified channels exceed a given criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).

# Returns
A Boolean vector indicating whether any extreme values were found for each row. Only samples selected by sample_selection are checked for extreme values.

# Examples
```julia
# Check extreme values in layout channels only (default)
is_extreme_value(dat, 100)

# Check extreme values in specific layout channels
is_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Check extreme values only in good samples (exclude already detected extreme values)
is_extreme_value(dat, 100, sample_selection = samples_not(:is_extreme_value_100))

# Check extreme values only within epoch windows
is_extreme_value(dat, 100, sample_selection = samples(:epoch_window))

# Check extreme values in additional channels (automatically detected)
is_extreme_value(dat, 100, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Exclude reference channels from layout
is_extreme_value(dat, 100, channel_selection = channels_not([:M1, :M2]))

# Combined filtering: check specific channels only in good samples
is_extreme_value(dat, 100, 
    sample_selection = samples_and([:epoch_window, samples_not(:is_vEOG)]),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)
```
"""
function is_extreme_value(
    dat::ContinuousData,
    criterion::Number;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_additional_channels::Bool = false,
)::Vector{Bool}

    selected_channels =
        get_selected_channels(dat, channel_selection; include_additional_channels = include_additional_channels)
    selected_samples = get_selected_samples(dat, sample_selection)

    @info "is_extreme_value: Checking for extreme values in channel $(print_vector(selected_channels)) with criterion $(criterion)"
    return _is_extreme_value(dat.data, criterion, selected_channels, selected_samples)
end

"""
    is_extreme_value!(dat::ContinuousData, criterion::Number; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false, channel_out::Symbol = :is_extreme_value)

Checks if any values in the specified channels exceed a given criterion and adds the result as a new column.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).
- `channel_out::Symbol`: Name of the output column (default: :is_extreme_value).

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
# Check extreme values in layout channels only (default)
is_extreme_value!(dat, 100)

# Check extreme values in specific layout channels
is_extreme_value!(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Check extreme values only in good samples (exclude already detected extreme values)
is_extreme_value!(dat, 100, sample_selection = samples_not(:is_extreme_value_100))

# Check extreme values only within epoch windows
is_extreme_value!(dat, 100, sample_selection = samples(:epoch_window))

# Check extreme values in additional channels (automatically detected)
is_extreme_value!(dat, 100, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Exclude reference channels from layout
is_extreme_value!(dat, 100, channel_selection = channels_not([:M1, :M2]))

# Combined filtering: check specific channels only in good samples
is_extreme_value!(dat, 100, 
    sample_selection = samples_and([:epoch_window, samples_not(:is_vEOG)]),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)
```
"""
function is_extreme_value!(
    dat::ContinuousData,
    criterion::Number;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
    channel_out::Symbol = :is_extreme_value,
)
    selected_channels =
        get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)

    @info "is_extreme_value!: Checking for extreme values in channel $(print_vector(selected_channels)) with criterion $(criterion)"
    dat.data[!, channel_out] = _is_extreme_value(dat.data, criterion, selected_channels, selected_samples)
end



"""
    n_extreme_value(dat::ContinuousData, criterion::Number; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false)::Int

Counts the number of extreme values in the specified channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).

# Returns
An integer count of the number of extreme values found in the selected samples and channels.

# Examples
```julia
# Count extreme values in layout channels only (default)
n_extreme_value(dat, 100)

# Count extreme values in specific layout channels
n_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Count extreme values only in good samples (exclude already detected extreme values)
n_extreme_value(dat, 100, sample_selection = samples_not(:is_extreme_value_100))

# Count extreme values only within epoch windows
n_extreme_value(dat, 100, sample_selection = samples(:epoch_window))

# Count extreme values excluding reference channels from layout
n_extreme_value(dat, 100, channel_selection = channels_not([:M1, :M2]))

# Count extreme values in additional channels (EOG, extreme value flags, etc.)
n_extreme_value(dat, 100, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Combined filtering: count extreme values in specific channels only in good samples
n_extreme_value(dat, 100, 
    sample_selection = samples_and([:epoch_window, samples_not(:is_vEOG)]),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)
```
"""
function n_extreme_value(
    dat::ContinuousData,
    criterion::Number;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_additional_channels::Bool = false,
)::Int
    selected_channels =
        get_selected_channels(dat, channel_selection; include_additional_channels = include_additional_channels)
    selected_samples = get_selected_samples(dat, sample_selection)

    return _n_extreme_value(dat.data, criterion, selected_channels, selected_samples)
end

# Internal function for plain DataFrames with explicit channel specification
function _n_extreme_value(
    dat::DataFrame,
    criterion::Number,
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Int
    @info "n_extreme_value: Counting extreme values in channel $(_print_vector(selected_channels)) with criterion $(criterion)"

    # Only count extreme values in selected samples
    if isempty(selected_samples)
        return 0
    end

    data_subset = select(dat[selected_samples, :], selected_channels)
    return sum(sum.(eachcol(abs.(data_subset) .>= criterion)))
end

# Internal function for plain DataFrames with explicit channel specification
function _channel_joint_probability(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol};
    threshold::Float64 = 5.0,
    normval::Int = 2,
)::DataFrame
    @info "channel_joint_probability: Computing probability for channels $(_print_vector(selected_channels))"

    # Select the specified channels and filter by samples
    data = select(dat[selected_samples, :], selected_channels)

    # Convert to matrix and compute joint probability
    jp, indelec = _joint_probability(Matrix(data)', threshold, normval)
    return DataFrame(channel = selected_channels, jp = jp, rejection = indelec)
end

"""
    channel_joint_probability(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false, threshold::Real = 3.0, normval::Real = 2)::DataFrame

Computes joint probability for EEG channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).
- `threshold::Real`: Threshold for joint probability (default: 3.0).
- `normval::Real`: Normalization value (default: 2).

# Returns
A DataFrame containing joint probability values for each channel.

# Examples
```julia
# Basic joint probability for layout channels only (default)
channel_joint_probability(dat)

# Filter samples where epoch_window is true
channel_joint_probability(dat, sample_selection = samples(:epoch_window))

# Filter to specific layout channels
channel_joint_probability(dat, channel_selection = channels([:Fp1, :Fp2]))

# Filter to additional channels (EOG, extreme value flags, etc.)
channel_joint_probability(dat, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Combine both filters
channel_joint_probability(dat, 
    sample_selection = samples(:epoch_window),
    channel_selection = channels_not([:M1, :M2])
)
```
"""
function channel_joint_probability(
    dat::ContinuousData;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_additional_channels::Bool = false,
    threshold::Real = 3.0,
    normval::Real = 2,
)::DataFrame
    selected_channels =
        get_selected_channels(dat, channel_selection; include_additional_channels = include_additional_channels)
    selected_samples = get_selected_samples(dat, sample_selection)

    return _channel_joint_probability(
        dat.data,
        selected_samples,
        selected_channels;
        threshold = threshold,
        normval = normval,
    )
end


function _joint_probability(signal::AbstractMatrix{Float64}, threshold::Float64, normalize::Int, discret::Int = 1000)
    nbchan = size(signal, 1)
    jp = zeros(nbchan)
    dataProba = Vector{Float64}(undef, size(signal, 2)) # Pre-allocate

    @inbounds for rc = 1:nbchan
        compute_probability!(dataProba, view(signal, rc, :), discret)
        jp[rc] = -sum(log, dataProba)
    end

    # Normalize the joint probability
    if normalize != 0
        tmpjp = normalize == 2 ? _trim_extremes(jp) : jp
        jp .= (jp .- mean(tmpjp)) ./ std(tmpjp)
    end

    rej = threshold != 0 ? abs.(jp) .> threshold : falses(nbchan)
    return jp, rej
end

"""
    compute_probability!(probaMap::Vector{Float64}, data::AbstractVector{Float64}, bins::Int)::Vector{Float64}

Computes the probability of each value in the data vector.

# Arguments
- `probaMap::Vector{Float64}`: The vector to store the computed probabilities.
- `data::AbstractVector{Float64}`: The data vector to compute the probabilities for.
- `bins::Int`: The number of bins to use for the probability computation.

# Returns
A vector of probabilities.

"""
function compute_probability!(probaMap::Vector{Float64}, data::AbstractVector{Float64}, bins::Int)::Vector{Float64}

    if bins > 0
        min_val, max_val = extrema(data)
        range_val = max_val - min_val
        sortbox = zeros(Int, bins)

        # Single-pass binning and counting
        @inbounds for x in data
            bin = clamp(floor(Int, (x - min_val) / range_val * (bins - 1)) + 1, 1, bins)
            sortbox[bin] += 1
        end

        # Compute probabilities
        n = length(data)
        @inbounds for (i, x) in enumerate(data)
            bin = clamp(floor(Int, (x - min_val) / range_val * (bins - 1)) + 1, 1, bins)
            probaMap[i] = sortbox[bin] / n
        end
    else
        # Gaussian approximation
        μ, σ = mean(data), std(data)
        inv_sqrt2pi = 1 / (√(2π))
        @inbounds for (i, x) in enumerate(data)
            z = (x - μ) / σ
            probaMap[i] = exp(-0.5 * z * z) * inv_sqrt2pi
        end
        sum_p = sum(probaMap)
        probaMap ./= sum_p
    end

    return probaMap

end


function _trim_extremes(x::Vector{Float64})
    n = length(x)
    trim = round(Int, n * 0.1)
    sorted = sort(x)
    return view(sorted, (trim+1):(n-trim))
end



"""
    get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real, <:Real})

Calculates the mean amplitude for each electrode within a specified time window.

# Arguments
- `erp_data::ErpData`: ERP data structure
- `time_window::Tuple{<:Real, <:Real}`: Time window as (start_time, end_time) in seconds

# Returns
- `DataFrame`: A DataFrame with electrode labels as column names and corresponding mean amplitudes
"""
function get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real,<:Real})
    # Find time indices within the window
    time_indices = findall(x -> x >= time_window[1] && x <= time_window[2], erp_data.time)

    if isempty(time_indices)
        error("No data points found within the specified time window")
    end

    # Calculate mean amplitude for each electrode
    mean_amplitudes = Dict{Symbol,Float64}()
    for electrode in erp_data.layout.label
        if haskey(erp_data.data, electrode)
            mean_amplitudes[electrode] = mean(erp_data.data[time_indices, electrode])
        end
    end

    return DataFrame(mean_amplitudes)
end







"""
    trigger_count(dat::ContinuousData; print_table::Bool = true)::DataFrame

Counts the number of occurrences of each trigger value in the data and optionally prints a formatted table.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `print_table::Bool`: Whether to print the trigger count table (default: true).

# Returns
A DataFrame with columns `trigger` and `count` showing trigger values and their counts, excluding zero values.

# Examples
```julia
# Get trigger counts and print table
trigger_counts = trigger_count(dat)

# Get trigger counts without printing
trigger_counts = trigger_count(dat, print_table = false)
```
"""
function trigger_count(dat::ContinuousData; print_table::Bool = true)::DataFrame
    @assert hasproperty(dat.data, :triggers) "Data must have a triggers column"
    return _trigger_count_impl(dat.data.triggers, print_table = print_table)
end

# Helper function for trigger counting logic
function _trigger_count_impl(
    trigger_data::Vector{<:Integer};
    print_table::Bool = true,
    title::String = "Trigger Count Summary",
)
    # Get unique non-zero trigger values
    unique_triggers = unique(trigger_data)
    non_zero_triggers = filter(x -> x != 0, unique_triggers)

    if isempty(non_zero_triggers)
        if print_table
            println("No non-zero triggers found in the data.")
        end
        return DataFrame()
    end

    # Count occurrences of each trigger value
    trigger_values = Int[]
    trigger_counts = Int[]

    for trigger in sort(non_zero_triggers)
        push!(trigger_values, trigger)
        push!(trigger_counts, count(x -> x == trigger, trigger_data))
    end

    # Create DataFrame
    result_df = DataFrame(trigger = trigger_values, count = trigger_counts)

    # Print table if requested
    if print_table
        pretty_table(result_df, title = title, header = ["Trigger", "Count"], alignment = [:r, :r], crop = :none)
        println()
    end

    return result_df
end

# Helper function for BioSemi data with both raw and cleaned counts
function _trigger_count_biosemi_impl(
    raw_triggers::Vector{<:Integer},
    cleaned_triggers::Vector{<:Integer};
    print_table::Bool = true,
)
    # Get unique non-zero trigger values from both raw and cleaned data
    unique_triggers = unique(vcat(raw_triggers, cleaned_triggers))
    non_zero_triggers = filter(x -> x != 0, unique_triggers)

    if isempty(non_zero_triggers)
        if print_table
            @minimal_warning "No non-zero triggers found in the data."
        end
        return DataFrame()
    end

    # Count occurrences of each trigger value in both raw and cleaned data
    trigger_values = Int[]
    raw_counts = Int[]
    cleaned_counts = Int[]

    for trigger in sort(non_zero_triggers)
        push!(trigger_values, trigger)
        push!(raw_counts, count(x -> x == trigger, raw_triggers))
        push!(cleaned_counts, count(x -> x == trigger, cleaned_triggers))
    end

    # Create DataFrame
    result_df = DataFrame(trigger = trigger_values, raw_count = raw_counts, cleaned_count = cleaned_counts)

    # Print table if requested
    if print_table
        pretty_table(
            result_df,
            title = "Trigger Count Summary (Raw vs Cleaned)",
            header = ["Trigger", "Raw Count", "Cleaned Count"],
            alignment = [:r, :r, :r],
            crop = :none,
        )
        println()
        println("Note: Cleaned counts show only trigger onset events (sustained signals converted to single onsets)")
    end

    return result_df
end


"""
    trigger_count(dat::BioSemiBDF.BioSemiData; print_table::Bool = true)::DataFrame

Counts the number of occurrences of each trigger value in BioSemi data (both raw and cleaned) and optionally prints a formatted table.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemiData object containing EEG data.
- `print_table::Bool`: Whether to print the trigger count table (default: true).

# Returns
A DataFrame with columns `trigger`, `raw_count`, and `cleaned_count` showing trigger values and their counts in both raw and cleaned data, excluding zero values.

# Examples
```julia
# Get trigger counts and print table
trigger_counts = trigger_count(dat)

# Get trigger counts without printing
trigger_counts = trigger_count(dat, print_table = false)
```
"""
function trigger_count(dat::BioSemiBDF.BioSemiData; print_table::Bool = true)::DataFrame
    # Get cleaned trigger data (onset detection only)
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    return _trigger_count_biosemi_impl(dat.triggers.raw, cleaned_triggers, print_table = print_table)
end







"""
    combine_boolean_columns!(dat::ContinuousData, columns::Vector{Symbol}, operation::Symbol; output_column::Symbol = :combined_flags)

Combines multiple boolean columns using a specified logical operation and stores the result in a new column.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing the boolean columns
- `columns::Vector{Symbol}`: Vector of column names to combine
- `operation::Symbol`: Logical operation to apply (:and, :or, :xor, :nand, :nor, :xnor)
- `output_column::Symbol`: Name of the output column (default: :combined_flags)

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
# Combine EOG artifacts (any EOG = vertical OR horizontal)
combine_boolean_columns!(dat, [:is_vEOG, :is_hEOG], :or, output_column = :is_any_EOG)

# Combine quality flags (good data = NOT extreme AND NOT EOG)
combine_boolean_columns!(dat, [:is_extreme_value, :is_any_EOG], :nor, output_column = :is_good_data)

# Combine multiple artifact types (any artifact)
combine_boolean_columns!(dat, [:is_extreme_value, :is_vEOG, :is_hEOG], :or, output_column = :is_any_artifact)

# Create clean data flag (no artifacts)
combine_boolean_columns!(dat, [:is_extreme_value, :is_vEOG, :is_hEOG], :nor, output_column = :is_clean_data)
```

# Available Operations
- `:and` - All columns must be true (logical AND)
- `:or` - At least one column must be true (logical OR)
- `:nand` - Not all columns are true (logical NAND)
- `:nor` - No columns are true (logical NOR)
"""
function combine_boolean_columns!(
    dat::ContinuousData,
    columns::Vector{Symbol},
    operation::Symbol;
    output_column::Symbol = :combined_flags,
)
    # Input validation
    @assert !isempty(columns) "Must specify at least one column to combine"
    @assert all(col -> hasproperty(dat.data, col), columns) "All specified columns must exist in the data"
    @assert operation in [:and, :or, :nand, :nor] "Invalid operation. Must be one of: :and, :or, :nand, :nor"

    # Get the boolean columns
    bool_columns = [dat.data[!, col] for col in columns]

    # Apply the logical operation
    result = if operation == :and
        all.(zip(bool_columns...))
    elseif operation == :or
        any.(zip(bool_columns...))
    elseif operation == :nand
        .!(all.(zip(bool_columns...)))
    elseif operation == :nor
        .!(any.(zip(bool_columns...)))
    end

    # Store the result
    dat.data[!, output_column] = result

    @info "combine_boolean_columns!: Combined $(length(columns)) columns using :$operation operation into column :$output_column"
end
