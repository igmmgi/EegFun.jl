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
```julia
# Asymmetric window: 1 second before to 2 seconds after trigger
mark_epoch_windows!(dat, [1, 3], [-1.0, 2.0])

# Custom column name
mark_epoch_windows!(dat, [1, 3], [-0.5, 0.5], channel_out = :near_trigger)
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
samples(column::Symbol) = x -> x[!, column]
samples(column::Symbol, value) = x -> x[!, column] .== value
samples_not(column::Symbol) = x -> .!(x[!, column])
samples() = x -> fill(true, nrow(x))

# Main method with keyword arguments for predicates
function channel_summary(dat::SingleDataFrameEeg; 
                        samples::Function = samples(),
                        channels::Function = channels)::DataFrame
    # Filter samples
    sample_mask = samples(dat.data)
    filtered_data = dat.data[sample_mask, :]
    
    # Filter channels - channels should return boolean vector
    channel_mask = channels(dat.layout.label)
    selected_channels = dat.layout.label[channel_mask]
    
    return _channel_summary_impl(filtered_data, selected_channels)
end

# For MultiDataFrameEeg
function channel_summary(dat::MultiDataFrameEeg; 
                        samples::Function = samples(),
                        channels::Function = channels)::Vector{DataFrame}
    return [channel_summary(dat.data[trial]; samples = samples, channels = channels) for trial in eachindex(dat.data)]
end







"""
    correlation_matrix(dat::ContinuousData; 
                      samples::Function = samples(),
                      channels::Function = channels())::DataFrame

Calculates the correlation matrix for the EEG data.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `samples::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
A DataFrame containing the correlation matrix of the specified channels.

# Examples
```julia
# Basic correlation matrix
correlation_matrix(dat)

# Filter samples where epoch_window is true
correlation_matrix(dat, samples = samples(:epoch_window))

# Filter to specific channels
correlation_matrix(dat, channels = channels([:Fp1, :Fp2]))

# Combine both filters
correlation_matrix(dat, 
    samples = samples(:epoch_window),
    channels = channels_not([:M1, :M2])
)
```
"""
function correlation_matrix(dat::ContinuousData; 
                          samples::Function = samples(),
                          channels::Function = channels())::DataFrame
    # Use layout.label as source of truth for EEG channels
    eeg_channels = dat.layout.label
    channel_mask = channels(eeg_channels)
    selected_channels = eeg_channels[channel_mask]
    
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
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event. This parameter can accept both `Symbol` and `String` types.

# Returns
Nothing. The function modifies the input data in place.

"""
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)
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
    is_extreme_value(dat::ContinuousData, criterion::Number; 
                     channels::Function = channels())::Vector{Bool}

Checks if any values in the specified channels exceed a given criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `channels::Function`: Function that returns boolean vector for channel filtering (default: include all channels).

# Returns
A Boolean vector indicating whether any extreme values were found for each row.

# Examples
```julia
# Check all channels
is_extreme_value(dat, 100)

# Check only specific channels
is_extreme_value(dat, 100, channels = channels([:Fp1, :Fp2]))

# Exclude reference channels
is_extreme_value(dat, 100, channels = channels_not([:M1, :M2]))
```
"""
function is_extreme_value(dat::ContinuousData, criterion::Number; 
                         channels::Function = channels())::Vector{Bool}
    # Use layout.label as source of truth for EEG channels
    eeg_channels = dat.layout.label
    channel_mask = channels(eeg_channels)
    selected_channels = eeg_channels[channel_mask]
    
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
# Check all channels
is_extreme_value!(dat, 100)

# Check only specific channels
is_extreme_value!(dat, 100, channels = channels([:Fp1, :Fp2]))

# Exclude reference channels
is_extreme_value!(dat, 100, channels = channels_not([:M1, :M2]))
```
"""
function is_extreme_value!(dat::ContinuousData, criterion::Number; 
                          channels::Function = channels(),
                          channel_out::Symbol = :is_extreme_value)
    # Use layout.label as source of truth for EEG channels
    eeg_channels = dat.layout.label
    channel_mask = channels(eeg_channels)  # This returns a boolean vector
    selected_channels = eeg_channels[channel_mask]
    
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
    # Use layout.label as source of truth for EEG channels
    eeg_channels = dat.layout.label
    channel_mask = channels(eeg_channels)
    selected_channels = eeg_channels[channel_mask]
    
    return _n_extreme_value(dat.data, criterion, selected_channels)
end

# Internal function for plain DataFrames with explicit channel specification
function _n_extreme_value(dat::DataFrame, criterion::Number, selected_channels::Vector{Symbol})::Int
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
    # Use layout.label as source of truth for EEG channels
    eeg_channels = dat.layout.label
    channel_mask = channels(eeg_channels)
    selected_channels = eeg_channels[channel_mask]
    
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
    # Unpack time window values
    start_time, end_time = time_window
    
    # Get time column name (assuming first column is time)
    time_col = first(names(erp_data.data))
    
    # Find indices corresponding to the time window
    time_indices = findall(t -> start_time <= t <= end_time, erp_data.data[:, time_col])
    
    if isempty(time_indices)
        throw(ArgumentError("No data points found in the specified time window ($start_time, $end_time)"))
    end
    
    # Get electrode names (all columns except the time column)
    electrode_names = filter(col -> col != time_col, names(erp_data.data))
    
    # Calculate mean amplitudes for each electrode
    mean_amplitudes = Dict{String, Float64}()
    for electrode in electrode_names
        mean_amplitudes[electrode] = mean(erp_data.data[time_indices, electrode])
    end
    
    # Return results as DataFrame
    return DataFrame(mean_amplitudes)
end

"""
    get_peak_latency(erp_data::ErpData, time_window::Tuple{<:Real, <:Real}; 
                    peak_type::Symbol=:positive)

Finds the latency (time point) of peaks for each electrode within a time window.

# Arguments
- `erp_data::ErpData`: ERP data structure
- `time_window::Tuple{<:Real, <:Real}`: Time window as (start_time, end_time) in seconds
- `peak_type::Symbol=:positive`: Type of peak to find (:positive for maximum, :negative for minimum)

# Returns
- `DataFrame`: A DataFrame with electrode labels as column names and corresponding peak latencies
"""
function get_peak_latency(
    erp_data::ErpData, 
    time_window::Tuple{<:Real, <:Real}; 
    peak_type::Symbol=:positive
)
    # Get time column name (assuming first column is time)
    time_col = first(names(erp_data.data))
    
    # Find indices corresponding to the time window
    time_values = erp_data.data[:, time_col]
    time_indices = findall(t -> time_window[1] <= t <= time_window[2], time_values)
    
    if isempty(time_indices)
        throw(ArgumentError("No data points found in the specified time window $(time_window)"))
    end
    
    # Get electrode names (all columns except the time column)
    electrode_names = filter(col -> col != time_col, names(erp_data.data))
    
    # Calculate peak latencies for each electrode
    peak_latencies = Dict{String, Float64}()
    
    for electrode in electrode_names
        data_window = erp_data.data[time_indices, electrode]
        
        # Find peak index
        peak_idx = if peak_type == :positive
            argmax(data_window)
        elseif peak_type == :negative
            argmin(data_window)
        else
            throw(ArgumentError("peak_type must be :positive or :negative"))
        end
        
        # Get the time value at the peak
        peak_time = time_values[time_indices[peak_idx]]
        peak_latencies[electrode] = peak_time
    end
    
    # Return results as DataFrame (same format as get_mean_amplitude)
    return DataFrame(peak_latencies)
end

"""
    get_peak_amplitude(erp_data::ErpData, time_window::Tuple{<:Real, <:Real}; 
                      peak_type::Symbol=:positive)

Finds the peak amplitude for each electrode within a time window.

# Arguments
- `erp_data::ErpData`: ERP data structure
- `time_window::Tuple{<:Real, <:Real}`: Time window as (start_time, end_time) in seconds
- `peak_type::Symbol=:positive`: Type of peak to find (:positive for maximum, :negative for minimum)

# Returns
- `DataFrame`: A DataFrame with electrode labels as column names and corresponding peak amplitudes
"""
function get_peak_amplitude(
    erp_data::ErpData, 
    time_window::Tuple{<:Real, <:Real}; 
    peak_type::Symbol=:positive
)
    # Get time column name (assuming first column is time)
    time_col = first(names(erp_data.data))
    
    # Find indices corresponding to the time window
    time_indices = findall(t -> time_window[1] <= t <= time_window[2], erp_data.data[:, time_col])
    
    if isempty(time_indices)
        throw(ArgumentError("No data points found in the specified time window $(time_window)"))
    end
    
    # Get electrode names (all columns except the time column)
    electrode_names = filter(col -> col != time_col, names(erp_data.data))
    
    # Calculate peak amplitudes for each electrode
    peak_amplitudes = Dict{String, Float64}()
    
    for electrode in electrode_names
        data_window = erp_data.data[time_indices, electrode]
        
        # Find peak
        peak_value = if peak_type == :positive
            maximum(data_window)
        elseif peak_type == :negative
            minimum(data_window)
        else
            throw(ArgumentError("peak_type must be :positive or :negative"))
        end
        
        peak_amplitudes[electrode] = peak_value
    end
    
    # Return results as DataFrame (same format as get_mean_amplitude)
    return DataFrame(peak_amplitudes)
end
