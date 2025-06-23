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

# Mark windows around multiple triggers
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

# Combine conditions
# dat.data[!, :any_condition] = dat.data[!, :condition_1] .| 
#                              dat.data[!, :condition_2] .| 
#                              dat.data[!, :condition_3]
```

## Using with Sample Filtering
```julia
# Mark epoch windows
mark_epoch_windows!(dat, [1, 3], [-1.0, 2.0])

# Use in analysis (only samples within epochs)
# summary = channel_summary(dat, samples = samples(:epoch_window))

# Use in ICA (only good samples within epochs)
# ica_result = run_ica(dat, 
#     samples = samples_and([
#         :epoch_window,
#         samples_not(:is_extreme_value_100)
#     ])
# )
```

## Quality Control Workflow
```julia
# 1. Mark epoch windows
mark_epoch_windows!(dat, [1, 3, 5], [-1.0, 2.0])

# 2. Detect artifacts
is_extreme_value!(dat, 100, channel_out = :is_extreme_100)
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)

# 3. Get good samples within epochs
# good_epoch_samples = samples_and([
#     :epoch_window,
#     samples_not(:is_extreme_100),
#     samples_not(:is_vEOG)
# ])

# 4. Analyze good data
# summary = channel_summary(dat, samples = good_epoch_samples)
```

## Complex Trigger Sequences
```julia
# Define epoch conditions with wildcards and ranges
# condition1 = EpochCondition(name="condition_1", trigger_sequence=[1, :any, 3], reference_index=2)
# condition2 = EpochCondition(name="condition_2", trigger_ranges=[1:5, 10:15], reference_index=1)
# conditions = [condition1, condition2]

# Mark windows around trigger sequences
# mark_epoch_windows!(dat, conditions, [-1.0, 2.0])

# Custom column name
# mark_epoch_windows!(dat, conditions, [-0.5, 0.5], channel_out = :near_sequences)
```
"""
function mark_epoch_windows!(
    dat::ContinuousData,
    triggers_of_interest::Vector{Int},
    time_window::Vector{<:Real};
    channel_out::Symbol = :epoch_window,
)

    # Input validation
    @assert length(time_window) == 2 "Time window must have exactly 2 elements"
    @assert time_window[1] <= time_window[2] "Time window start must be less than or equal to end"
    @assert !isempty(triggers_of_interest) "Must specify at least one trigger of interest"
    @assert hasproperty(dat.data, :triggers) "Data must have a triggers column"
    @assert hasproperty(dat.data, :time) "Data must have a time column"

    # Initialize result vector with false
    dat.data[!, channel_out] .= false

    # Only need unique trigers
    unique_triggers = unique(dat.data.triggers)

    # For each trigger of interest
    for trigger in triggers_of_interest

        if !(trigger in unique_triggers)
            @minimal_warning "Trigger $trigger not found in data"
            continue
        end

        # Find all samples where this trigger occurs
        trigger_indices = findall(dat.data.triggers .== trigger)

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
    @assert length(time_window) == 2 "Time window must have exactly 2 elements"
    @assert time_window[1] <= time_window[2] "Time window start must be less than or equal to end"
    @assert !isempty(epoch_conditions) "Must specify at least one epoch condition"
    @assert hasproperty(dat.data, :triggers) "Data must have a triggers column"
    @assert hasproperty(dat.data, :time) "Data must have a time column"

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
            filtered_indices = Int[]
            
            for seq_start_idx in sequence_indices
                # Check if this sequence meets the after/before constraints
                valid_position = true
                
                if condition.after !== nothing
                    # Check if there's a trigger with value 'after' before this sequence
                    # Look backwards from sequence start to find the trigger
                    found_after_trigger = false
                    for i in (seq_start_idx-1):-1:1
                        if dat.data.triggers[i] == condition.after
                            found_after_trigger = true
                            break
                        end
                    end
                    if !found_after_trigger
                        valid_position = false
                    end
                end
                
                if condition.before !== nothing
                    # Check if there's a trigger with value 'before' after this sequence
                    # Look forwards from sequence end to find the trigger
                    sequence_end = seq_start_idx + length(condition.trigger_sequence) - 1
                    found_before_trigger = false
                    for i in (sequence_end+1):length(dat.data.triggers)
                        if dat.data.triggers[i] == condition.before
                            found_before_trigger = true
                            break
                        end
                    end
                    if !found_before_trigger
                        valid_position = false
                    end
                end
                
                if valid_position
                    push!(filtered_indices, seq_start_idx)
                end
            end
            
            sequence_indices = filtered_indices
            if isempty(sequence_indices)
                after_msg = condition.after !== nothing ? " after trigger $(condition.after)" : ""
                before_msg = condition.before !== nothing ? " before trigger $(condition.before)" : ""
                @warn "No trigger sequences found that meet position constraints$(after_msg)$(before_msg) for condition '$(condition.name)'"
                continue
            end
        end

        # Apply timing constraints if specified
        if condition.timing_pairs !== nothing && 
           condition.min_interval !== nothing && 
           condition.max_interval !== nothing
            
            valid_indices = Int[]
            
            for seq_start_idx in sequence_indices
                # Check if this sequence meets all timing constraints
                valid_sequence = true
                
                for (start_idx, end_idx) in condition.timing_pairs
                    # Calculate actual indices in the trigger array
                    actual_start_idx = seq_start_idx + (start_idx - 1)
                    actual_end_idx = seq_start_idx + (end_idx - 1)
                    
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
                    if interval < condition.min_interval || interval > condition.max_interval
                        valid_sequence = false
                        break
                    end
                end
                
                if valid_sequence
                    push!(valid_indices, seq_start_idx)
                end
            end
            
            sequence_indices = valid_indices
            if isempty(sequence_indices)
                @warn "No trigger sequences found that meet timing constraints for condition '$(condition.name)'"
                continue
            end
        end

        # For each valid sequence occurrence, find the reference point (t=0 position)
        for seq_start_idx in sequence_indices
            # Calculate the reference index (t=0 position) within the sequence
            reference_idx = seq_start_idx + (condition.reference_index - 1)
            
            # Check if reference index is within bounds
            if reference_idx > length(dat.data.triggers)
                @warn "Reference index $(condition.reference_index) for condition '$(condition.name)' is out of bounds"
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
    create_eeg_dataframe(data::BioSemiBDF.BioSemiData)::DataFrame

Creates a DataFrame containing EEG data from a BioSemiBDF object.

# Arguments
- `data::BioSemiBDF.BioSemiData`: The BioSemi data structure containing time, triggers, and channel data.

# Returns
A DataFrame with two columns: `time` and `triggers`, combined with channel data from the BioSemiBDF.

"""
function create_eeg_dataframe(data::BioSemiBDF.BioSemiData)::DataFrame
    return hcat(
        DataFrame(time = data.time, sample = 1:length(data.time), triggers = data.triggers.raw),
        DataFrame(Float64.(data.data), Symbol.(data.header.channel_labels[1:end-1])),  # assumes last channel is trigger
    )
end



"""
    create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout_file_name::String)::ContinuousData

Creates a ContinuousData object from a BioSemiBDF data structure and a layout file.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemi data structure containing EEG data.
- `layout_file_name::String`: The filename of the layout CSV file.

# Returns
A ContinuousData object containing the EEG data and layout information.

"""
function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout_file_name::String)::ContinuousData
    return ContinuousData(
        create_eeg_dataframe(dat),
        DataFrame(CSV.File(layout_file_name)),
        dat.header.sample_rate[1],
        AnalysisInfo(),
    )
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
function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout::DataFrame)::ContinuousData
    return ContinuousData(create_eeg_dataframe(dat), layout, dat.header.sample_rate[1], AnalysisInfo())
end



# Helper function with common logic
function _channel_summary_impl(selected_data::DataFrame, channel_labels::Vector{Symbol})::DataFrame
    # Initialize a matrix to store summary statistics
    # 6 statistics: min, max, range, std, mad, var
    summary_stats = Matrix{Float64}(undef, length(channel_labels), 6)

    # Compute summary statistics for each column
    for (i, col) in enumerate(channel_labels)
        col_data = @view selected_data[!, col]
        summary_stats[i, :] =
            [minimum(col_data), maximum(col_data), datarange(col_data), std(col_data), mad(col_data), var(col_data)]
    end

    # Create a new DataFrame directly from the matrix and channel labels
    summary_df = DataFrame(
        channel = channel_labels,
        min = summary_stats[:, 1],
        max = summary_stats[:, 2],
        range = summary_stats[:, 3],
        std = summary_stats[:, 4],
        mad = summary_stats[:, 5],
        var = summary_stats[:, 6],
    )

    # Compute z-scores for channel variances
    summary_df[!, :zvar] .= (summary_df[!, :var] .- mean(summary_df[!, :var])) ./ std(summary_df[!, :var])

    return summary_df
end

# Helper functions for easier channel filtering
channels(channel_names::Vector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_name::Symbol) = x -> x .== channel_name
channels(channel_numbers::Union{Vector{Int},UnitRange}) = x -> [i in channel_numbers for i in 1:length(x)]
channels() = x -> fill(true, length(x))
channels_not(channel_names::Vector{Symbol}) = x -> .!(x .∈ Ref(channel_names))
channels_not(channel_name::Symbol) = x -> .!(x .== channel_name)
channels_not(channel_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in channel_numbers for i in 1:length(x)])

# Helper functions for easier sample filtering
samples() = x -> fill(true, nrow(x))
samples(column::Symbol) = x -> x[!, column]
samples(columns::Vector{Symbol}) = x -> any(x[!, col] for col in columns)  # OR logic (backward compatibility)
samples_or(columns::Vector{Symbol}) = x -> any(x[!, col] for col in columns)
samples_and(columns::Vector{Symbol}) = x -> all(x[!, col] for col in columns)
samples_not(column::Symbol) = x -> .!(x[!, column])
samples_not(columns::Vector{Symbol}) = x -> .!(any(x[!, col] for col in columns))  # NOT OR = AND NOT
samples_or_not(columns::Vector{Symbol}) = x -> .!(any(x[!, col] for col in columns))  # NOT OR = AND NOT
samples_and_not(columns::Vector{Symbol}) = x -> .!(all(x[!, col] for col in columns))  # NOT AND = OR NOT

# Helper function to get available channels (layout channels only by default)
function _get_available_channels(dat::EegData)
    return dat.layout.label
end

# Helper function to get all available channels (layout + additional)
function _get_all_available_channels(dat::EegData)
    available_channels = dat.layout.label
    data_channels = propertynames(data(dat))
    additional_channels = setdiff(data_channels, [:time, :sample, :triggers, available_channels...])
    return [available_channels; additional_channels]
end

# Helper function for users to get all available channels
all_channels(dat::EegData) = _get_all_available_channels(dat)



"""
    channel_summary(dat::SingleDataFrameEeg; 
                   samples::Function = samples(),
                   channels::Function = channels())::DataFrame

Computes summary statistics for EEG channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object.
- `samples::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
A DataFrame containing summary statistics for each channel.

# Examples

## Basic Usage
```julia
# Channel summary for layout channels only (default)
summary = channel_summary(dat)

# Channel summary for specific layout channels
summary = channel_summary(dat, channels = channels([:Fp1, :Fp2, :F3, :F4]))

# Channel summary excluding reference channels from layout
summary = channel_summary(dat, channels = channels_not([:M1, :M2]))
```

## Including Additional Channels
```julia
# To include additional channels (EOG, reference, etc.), explicitly specify them
all_channels = _get_all_available_channels(dat)
summary = channel_summary(dat, channels = channels(all_channels))

# Or include specific additional channels
summary = channel_summary(dat, channels = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Channel Filtering
```julia
# Summary for specific channels
# summary = channel_summary(dat, channels = channels([:Fp1, :Fp2, :F3, :F4]))

# Summary excluding reference channels
# summary = channel_summary(dat, channels = channels_not([:M1, :M2]))

# Summary for frontal channels only (channels 1-10)
# summary = channel_summary(dat, channels = channels(1:10))

# Summary for channels starting with "F" (frontal)
# summary = channel_summary(dat, channels = x -> startswith.(string.(x), "F"))
```

## Sample Filtering
```julia
# Exclude extreme values
# summary = channel_summary(dat, samples = samples_not(:is_extreme_value_100))

# Exclude multiple types of bad samples
# summary = channel_summary(dat, samples = samples_or_not([:is_extreme_value_100, :is_vEOG, :is_hEOG]))

# Only include samples within epoch windows
# summary = channel_summary(dat, samples = samples(:epoch_window))

# Include samples that are both in epoch window AND not extreme
# summary = channel_summary(dat, samples = samples_and([:epoch_window, samples_not(:is_extreme_value_100)]))
```

## Combined Filtering
```julia
# Exclude reference channels and extreme values
# summary = channel_summary(dat, 
#     channels = channels_not([:M1, :M2]),
#     samples = samples_not(:is_extreme_value_100)
# )

# Only frontal channels, exclude bad samples
# summary = channel_summary(dat, 
#     channels = channels(1:10),
#     samples = samples_or_not([:is_extreme_value_100, :is_vEOG])
# )

# Complex filtering: frontal channels, good samples, within epochs
# summary = channel_summary(dat, 
#     channels = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
#     samples = samples_and([
#         :epoch_window, 
#         samples_not(:is_extreme_value_100),
#         samples_not(:is_vEOG),
#         samples_not(:is_hEOG)
#     ])
# )
```

## Additional Channels (not in layout)
```julia
# Include derived channels like EOG
# summary = channel_summary(dat, channels = channels([:vEOG, :hEOG]))

# Mix layout channels and additional channels
# summary = channel_summary(dat, channels = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```
"""
function channel_summary(dat::SingleDataFrameEeg; 
                        samples::Function = samples(),
                        channels::Function = channels())::DataFrame
    # Filter samples
    sample_mask = samples(dat.data)
    filtered_data = dat.data[sample_mask, :]
    
    # Always use all available channels as base, let channels() function handle filtering
    all_available_channels = _get_all_available_channels(dat)
    
    channel_mask = channels(all_available_channels)
    selected_channels = all_available_channels[channel_mask]
    
    @info "channel_summary: Selected channels: $(_print_vector(selected_channels))"
    
    return _channel_summary_impl(filtered_data, selected_channels)
end

# For MultiDataFrameEeg
function channel_summary(dat::MultiDataFrameEeg; 
                        samples::Function = samples(),
                        channels::Function = channels())::Vector{DataFrame}
    return [channel_summary(dat.data[trial]; samples = samples, channels = channels) for trial in eachindex(dat.data)]
end







"""
    correlation_matrix(dat::ContinuousData; 
                      samples::Function = samples(),
                      channels::Function = channels(dat))::DataFrame

Calculates the correlation matrix for the EEG data.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `samples::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
A DataFrame containing the correlation matrix of the specified channels.

# Examples

## Basic Usage
```julia
# Correlation matrix for layout channels only (default)
corr_matrix = correlation_matrix(dat)

# Correlation matrix for specific layout channels
corr_matrix = correlation_matrix(dat, channels = channels([:Fp1, :Fp2, :F3, :F4]))

# Correlation matrix excluding reference channels from layout
corr_matrix = correlation_matrix(dat, channels = channels_not([:M1, :M2]))
```

## Including Additional Channels
```julia
# To include additional channels (EOG, reference, etc.), explicitly specify them
all_channels = _get_all_available_channels(dat)
corr_matrix = correlation_matrix(dat, channels = channels(all_channels))

# Or include specific additional channels
corr_matrix = correlation_matrix(dat, channels = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Sample Filtering
```julia
# Correlation matrix only for good samples
# corr_matrix = correlation_matrix(dat, samples = samples_not(:is_extreme_value_100))

# Correlation matrix only within epoch windows
# corr_matrix = correlation_matrix(dat, samples = samples(:epoch_window))

# Correlation matrix for good samples within epochs
# corr_matrix = correlation_matrix(dat, 
#     samples = samples_and([
#         :epoch_window,
#         samples_not(:is_extreme_value_100)
#     ])
# )
```

## Channel Filtering
```julia
# Frontal channels only
# corr_matrix = correlation_matrix(dat, channels = channels(1:10))

# Channels starting with "F" (frontal)
# corr_matrix = correlation_matrix(dat, channels = x -> startswith.(string.(x), "F"))

# Parietal channels only
# corr_matrix = correlation_matrix(dat, channels = x -> startswith.(string.(x), "P"))

# Mix of layout and additional channels
# corr_matrix = correlation_matrix(dat, channels = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Combined Filtering
```julia
# Exclude reference channels and bad samples
# corr_matrix = correlation_matrix(dat, 
#     channels = channels_not([:M1, :M2]),
#     samples = samples_not(:is_extreme_value_100)
# )

# Only frontal channels, good samples
# corr_matrix = correlation_matrix(dat, 
#     channels = channels(1:10),
#     samples = samples_and([
#         :epoch_window, 
#         samples_not(:is_extreme_value_100)
#     ])
# )

# Complex filtering: frontal channels, good samples, within epochs
# corr_matrix = correlation_matrix(dat, 
#     channels = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
#     samples = samples_and([
#         :epoch_window, 
#         samples_not(:is_extreme_value_100),
#         samples_not(:is_vEOG),
#         samples_not(:is_hEOG)
#     ])
# )
```

## Quality Control Applications
```julia
# 1. Detect extreme values
is_extreme_value!(dat, 100, channel_out = :is_extreme_100)

# 2. Get correlation matrix for good data
# good_corr = correlation_matrix(dat, samples = samples_not(:is_extreme_100))

# 3. Check for highly correlated channels (potential duplicates)
# Look for correlation values above 0.95
# Highly correlated channels may indicate duplicates or artifacts
```

## Regional Analysis
```julia
# Frontal correlation matrix
# frontal_corr = correlation_matrix(dat, 
#     channels = channels(x -> startswith.(string.(x), "F")),
#     samples = samples(:epoch_window)
# )

# Parietal correlation matrix
# parietal_corr = correlation_matrix(dat, 
#     channels = channels(x -> startswith.(string.(x), "P")),
#     samples = samples(:epoch_window)
# )

# Compare frontal vs parietal connectivity
# println("Frontal average correlation: ...")
# println("Parietal average correlation: ...")
```

## Time-Based Analysis
```julia
# Correlation matrix for different time periods
# early_corr = correlation_matrix(dat, 
#     samples = samples_and([:epoch_window, x -> x.time .< 0.2])  # First 200ms
# )

# late_corr = correlation_matrix(dat, 
#     samples = samples_and([:epoch_window, x -> x.time .> 0.3])  # After 300ms
# )

# Compare early vs late connectivity
# println("Early connectivity: ...")
# println("Late connectivity: ...")
```

## Visualization
```julia
# Get correlation matrix
# corr_matrix = correlation_matrix(dat, channels = channels_not([:M1, :M2]))

# Convert to matrix for plotting
# corr_values = Matrix(corr_matrix[:, 2:end])
# channel_names = corr_matrix.row

# Plot as heatmap (using your preferred plotting package)
# heatmap(corr_values, xticks=channel_names, yticks=channel_names)
```
"""
function correlation_matrix(dat::ContinuousData; 
                          samples::Function = samples(),
                          channels::Function = channels(dat))::DataFrame
    # Always use all available channels, let channels() function handle filtering
    all_available_channels = propertynames(dat.data)
    all_available_channels = setdiff(all_available_channels, [:time, :sample, :triggers])
    
    channel_mask = channels(all_available_channels)
    selected_channels = all_available_channels[channel_mask]
    
    return _correlation_matrix(dat.data, selected_channels; samples = samples)
end

# Internal function for plain DataFrames with explicit channel specification
function _correlation_matrix(
    dat::DataFrame,
    selected_channels::Vector{Symbol};
    samples::Function = samples(),
)::DataFrame
    # Select the specified channels
    data = select(dat, selected_channels)
    
    # Filter samples
    sample_mask = samples(dat)
    data = data[sample_mask, :]
    
    # Compute correlation matrix
    df = DataFrame(cor(Matrix(data)), selected_channels)
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
```

## Different Thresholds
```julia
# Conservative threshold (fewer detections)
detect_eog_onsets!(dat, 100.0, :vEOG, :is_vEOG_conservative)

# Liberal threshold (more detections)
detect_eog_onsets!(dat, 20.0, :vEOG, :is_vEOG_liberal)

# Standard threshold
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG_standard)
```

## Multiple EOG Channels
```julia
# Detect both vertical and horizontal EOG
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# Combine for any EOG artifact
# dat.data[!, :is_any_EOG] = dat.data[!, :is_vEOG] .| dat.data[!, :is_hEOG]
```

## Using with Sample Filtering
```julia
# Detect EOG artifacts
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# Use in analysis (exclude EOG artifacts)
# summary = channel_summary(dat, 
#     samples = samples_and([
#         samples_not(:is_vEOG),
#         samples_not(:is_hEOG)
#     ])
# )

# Use in ICA (exclude EOG artifacts)
# ica_result = run_ica(dat, 
#     samples = samples_and([
#         :epoch_window,
#         samples_not(:is_vEOG),
#         samples_not(:is_hEOG)
#     ])
# )
```

## Quality Control Workflow
```julia
# 1. Detect extreme values
is_extreme_value!(dat, 100, channel_out = :is_extreme_100)

# 2. Detect EOG artifacts
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# 3. Get summary of good data
# good_samples = samples_and([
#     samples_not(:is_extreme_100),
#     samples_not(:is_vEOG),
#     samples_not(:is_hEOG)
# ])

# summary = channel_summary(dat, samples = good_samples)
```

## Threshold Selection Tips
```julia
# Start with moderate threshold
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG_test)

# Check how many samples were detected
# println("Detected ... EOG samples out of ... total samples")

# Adjust threshold based on results
# if n_detected < 100
#     # Too few detections, lower threshold
#     detect_eog_onsets!(dat, 30.0, :vEOG, :is_vEOG)
# elseif n_detected > 1000
#     # Too many detections, raise threshold
#     detect_eog_onsets!(dat, 80.0, :vEOG, :is_vEOG)
# end
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
function _is_extreme_value(dat::DataFrame, criterion::Number, selected_channels::Vector{Symbol})::Vector{Bool}
    return any(x -> abs.(x) >= criterion, Matrix(select(dat, selected_channels)), dims = 2)[:]
end

# Internal function for plain DataFrames with explicit channel specification
function _is_extreme_value!(dat::DataFrame, criterion::Number, selected_channels::Vector{Symbol}; channel_out::Symbol = :is_extreme_value)
    dat[!, channel_out] .= any(x -> abs.(x) >= criterion, Matrix(select(dat, selected_channels)), dims = 2)[:]
end

"""
    is_extreme_value(dat::ContinuousData, criterion::Number; channels::Function = channels())::Vector{Bool}

Checks if any values in the specified channels exceed a given criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
A Boolean vector indicating whether any extreme values were found for each row.

# Examples
```julia
# Check extreme values in layout channels only (default)
is_extreme_value!(dat, 100)

# Check extreme values in specific layout channels
is_extreme_value!(dat, 100, channels = channels([:Fp1, :Fp2]))

# Exclude reference channels from layout
is_extreme_value!(dat, 100, channels = channels_not([:M1, :M2]))

# To include additional channels (EOG, reference, etc.), explicitly specify them
all_channels = _get_all_available_channels(dat)
is_extreme_value!(dat, 100, channels = channels(all_channels))
```
"""
function is_extreme_value(dat::ContinuousData, criterion::Number; channels::Function = channels())::Vector{Bool}
    # Always use all available channels as base, let channels() function handle filtering
    all_available_channels = _get_all_available_channels(dat)
    
    channel_mask = channels(all_available_channels)
    selected_channels = all_available_channels[channel_mask]
    @info "is_extreme_value!: Checking for extreme values in channel $(print_vector(selected_channels)) with criterion $(criterion)"
    return _is_extreme_value(dat.data, criterion, selected_channels)
end

"""
    is_extreme_value!(dat::ContinuousData, criterion::Number; 
                      channels::Function = channels(),
                      channel_out::Symbol = :is_extreme_value)

Checks if any values in the specified channels exceed a given criterion and adds the result as a new column.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `channel_out::Symbol`: Name of the output column (default: :is_extreme_value).

# Examples
```julia
# Check extreme values in layout channels only (default)
is_extreme_value!(dat, 100)

# Check extreme values in specific layout channels
is_extreme_value!(dat, 100, channels = channels([:Fp1, :Fp2]))

# Exclude reference channels from layout
is_extreme_value!(dat, 100, channels = channels_not([:M1, :M2]))

# To include additional channels (EOG, reference, etc.), explicitly specify them
all_channels = _get_all_available_channels(dat)
is_extreme_value!(dat, 100, channels = channels(all_channels))
```
"""
function is_extreme_value!(dat::ContinuousData, criterion::Number; channels::Function = channels(), channel_out::Symbol = :is_extreme_value)
    # Always use all available channels as base, let channels() function handle filtering
    all_available_channels = _get_all_available_channels(dat)
    
    channel_mask = channels(all_available_channels)
    selected_channels = all_available_channels[channel_mask]
    @info "is_extreme_value!: Checking for extreme values in channel $(print_vector(selected_channels)) with criterion $(criterion)"
    _is_extreme_value!(dat.data, criterion, selected_channels, channel_out = channel_out)
end



"""
    n_extreme_value(dat::ContinuousData, criterion::Number; 
                    channels::Function = channels())::Int

Counts the number of extreme values in the specified channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
An integer count of the number of extreme values found.

# Examples
```julia
# Count extreme values in all channels
n_extreme_value(dat, 100)

# Count extreme values in specific channels
n_extreme_value(dat, 100, channels = channels([:Fp1, :Fp2]))

# Count extreme values excluding reference channels
n_extreme_value(dat, 100, channels = channels_not([:M1, :M2]))
```
"""
function n_extreme_value(dat::ContinuousData, criterion::Number; channels::Function = channels())::Int
    # Start with layout channels only
    layout_channels = dat.layout.label  # Get layout channels directly
    
    # Apply channel filtering to get selected channels
    channel_mask = channels(layout_channels)
    selected_channels = layout_channels[channel_mask]
    
    # Check if user is trying to access channels that aren't in the layout
    # If so, switch to using all available channels
    if !isempty(selected_channels) && !all(ch -> ch in layout_channels, selected_channels)
        # User specified channels not in layout, use all available channels as base
        all_available_channels = _get_all_available_channels(dat)
        channel_mask = channels(all_available_channels)
        selected_channels = all_available_channels[channel_mask]
    end
    
    return _n_extreme_value(dat.data, criterion, selected_channels)
end

# Internal function for plain DataFrames with explicit channel specification
function _n_extreme_value(dat::DataFrame, criterion::Number, selected_channels::Vector{Symbol})::Int
    @info "n_extreme_value: Counting extreme values in channel $(_print_vector(selected_channels)) with criterion $(criterion)"
    return sum(sum.(eachcol(abs.(select(dat, selected_channels)) .>= criterion)))
end

# Internal function for plain DataFrames with explicit channel specification
function _channel_joint_probability(
    dat::DataFrame,
    selected_channels::Vector{Symbol};
    threshold::Float64 = 5.0,
    normval::Int = 2,
    samples::Function = samples(),
)::DataFrame
    @info "channel_joint_probability: Computing probability for channels $(_print_vector(selected_channels))"
    
    # Select the specified channels
    data = select(dat, selected_channels)
    
    # Filter samples
    sample_mask = samples(dat)
    data = data[sample_mask, :]
    
    # Convert to matrix and compute joint probability
    jp, indelec = joint_probability(Matrix(data)', threshold, normval)
    return DataFrame(channel = selected_channels, jp = jp, rejection = indelec)
end

"""
    channel_joint_probability(dat::ContinuousData; 
                             threshold::Float64 = 5.0, 
                             normval::Int = 2,
                             samples::Function = samples(),
                             channels::Function = channels())::DataFrame

Computes joint probability for EEG channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `threshold::Float64`: Threshold for joint probability (default: 5.0).
- `normval::Int`: Normalization value (default: 2).
- `samples::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
A DataFrame containing joint probability values for each channel.

# Examples
```julia
# Basic joint probability
channel_joint_probability(dat)

# Filter samples where epoch_window is true
channel_joint_probability(dat, samples = samples(:epoch_window))

# Filter to specific channels
channel_joint_probability(dat, channels = channels([:Fp1, :Fp2]))

# Combine both filters
channel_joint_probability(dat, 
    samples = samples(:epoch_window),
    channels = channels_not([:M1, :M2])
)
```
"""
function channel_joint_probability(dat::ContinuousData; 
                                 threshold::Float64 = 5.0, 
                                 normval::Int = 2,
                                 samples::Function = samples(),
                                 channels::Function = channels())::DataFrame
    # Always use all available channels as base, let channels() function handle filtering
    all_available_channels = _get_all_available_channels(dat)
    
    channel_mask = channels(all_available_channels)
    selected_channels = all_available_channels[channel_mask]
    
    return _channel_joint_probability(dat.data, selected_channels; 
                                    threshold = threshold, 
                                    normval = normval, 
                                    samples = samples)
end


function joint_probability(signal::AbstractMatrix{Float64}, threshold::Float64, normalize::Int, discret::Int = 1000)
    nbchan = size(signal, 1)
    jp = zeros(nbchan)
    dataProba = Vector{Float64}(undef, size(signal, 2)) # Pre-allocate

    @inbounds for rc = 1:nbchan
        compute_probability!(dataProba, view(signal, rc, :), discret)
        jp[rc] = -sum(log, dataProba)
    end

    # Normalize the joint probability
    if normalize != 0
        tmpjp = normalize == 2 ? trim_extremes(jp) : jp
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


function trim_extremes(x::Vector{Float64})
    n = length(x)
    trim = round(Int, n * 0.1)
    sorted = sort(x)
    return view(sorted, trim+1:n-trim)
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
function get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real, <:Real})
    # Find time indices within the window
    time_indices = findall(x -> x >= time_window[1] && x <= time_window[2], erp_data.time)
    
    if isempty(time_indices)
        error("No data points found within the specified time window")
    end
    
    # Calculate mean amplitude for each electrode
    mean_amplitudes = Dict{Symbol, Float64}()
    for electrode in erp_data.layout.label
        if haskey(erp_data.data, electrode)
            mean_amplitudes[electrode] = mean(erp_data.data[time_indices, electrode])
        end
    end
    
    return DataFrame(mean_amplitudes)
end






