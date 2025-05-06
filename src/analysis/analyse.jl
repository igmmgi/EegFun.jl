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
            @warn "Trigger $trigger not found in data"
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



function channel_summary(dat::DataFrame, channel_labels::Vector{Symbol}; filter_samples = nothing)::DataFrame

    # select the specified channels
    selected_data = select(dat, channel_labels)

    # filter samples if requested
    if filter_samples !== nothing
        if filter_samples isa Symbol && hasproperty(dat, filter_samples)
            # # TODO: I want this version only
            # If filter_samples is a column name, use that column
            selected_data = selected_data[dat[!, filter_samples], :]
        else
            # Otherwise assume it's a boolean vector
            selected_data = selected_data[filter_samples, :]
        end
    end

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

function channel_summary(dat::SingleDataFrameEeg; filter_samples = nothing)::DataFrame
    return channel_summary(dat.data, dat.layout.label; filter_samples = filter_samples)
end

function channel_summary(dat::SingleDataFrameEeg, channel_numbers::Union{Vector{Int},UnitRange}; filter_samples = nothing)::DataFrame
    channel_labels = channel_number_to_channel_label(dat.layout.label, channel_numbers)
    return channel_summary(dat.data, channel_labels; filter_samples = filter_samples)
end

function channel_summary(dat::SingleDataFrameEeg, channel_labels::Vector{Symbol}; filter_samples = nothing)::DataFrame
    return channel_summary(dat.data, channel_labels; filter_samples = filter_samples)
end

function channel_summary(dat::MultiDataFrameEeg; filter_samples = nothing)::Vector{DataFrame}
    return [channel_summary(dat.data[trial], dat.layout.label; filter_samples = filter_samples) for trial in eachindex(dat.data)]
end

function channel_summary(dat::MultiDataFrameEeg, channel_numbers::Union{Vector{Int},UnitRange}; filter_samples = nothing)::Vector{DataFrame}
    channel_labels = channel_number_to_channel_label(dat.layout.label, channel_numbers)
    return [channel_summary(dat.data[trial], channel_labels; filter_samples = filter_samples) for trial in eachindex(dat.data)]
end

function channel_summary(dat::MultiDataFrameEeg, channel_labels::Vector{Symbol}; filter_samples = nothing)::Vector{DataFrame}
    return [channel_summary(dat.data[trial], channel_labels; filter_samples = filter_samples) for trial in eachindex(dat.data)]
end







"""
    correlation_matrix(dat::DataFrame, layout::DataFrame)::DataFrame

Calculates the correlation matrix for the EEG data.

# Arguments
- `dat::DataFrame`: The DataFrame containing EEG data.
- `layout::DataFrame`: The DataFrame containing layout information.

# Returns
A DataFrame containing the correlation matrix of the specified channels.

"""
function correlation_matrix(dat::DataFrame, channel_labels::Vector{Symbol}; filter_samples = nothing)::DataFrame
    # Select the specified channels
    data = select(dat, channel_labels)
    
    # Filter samples if requested
    if filter_samples !== nothing
        if filter_samples isa Symbol && hasproperty(dat, filter_samples)
            # If filter_samples is a column name, use that column
            data = data[dat[!, filter_samples], :]
        else
            # Otherwise assume it's a boolean vector
            data = data[filter_samples, :]
        end
    end
    
    # Compute correlation matrix
    df = DataFrame(cor(Matrix(data)), channel_labels)
    insertcols!(df, 1, :row => channel_labels)
    return df
end

function correlation_matrix(dat::ContinuousData; filter_samples = nothing)::DataFrame
    return correlation_matrix(dat.data, dat.layout.label; filter_samples = filter_samples)
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

"""
    is_extreme_value(dat::DataFrame, columns::Vector{Symbol}, criterion::Real)::Bool

Checks if any values in the specified columns exceed a given criterion.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data to check.
- `columns::Vector{Symbol}`: The columns to check for extreme values.
- `criterion::Float64`: The threshold for determining extreme values.

# Returns
A Boolean indicating whether any extreme values were found.

"""
function is_extreme_value(dat::DataFrame, columns::Vector{Symbol}, criterion::Number)::Vector{Bool}
    return any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)), dims = 2)[:]
end

function is_extreme_value!(
    dat::DataFrame,
    columns::Vector{Symbol},
    criterion::Number;
    channel_out::Symbol = :is_extreme_value,
)
    dat[!, channel_out] .= any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)), dims = 2)[:]
end

function is_extreme_value!(
    dat::ContinuousData,
    columns::Vector{Symbol},
    criterion::Number;
    channel_out::Symbol = :is_extreme_value,
)
    is_extreme_value!(dat.data, columns, criterion, channel_out = channel_out)
end

function is_extreme_value!(
    dat::ContinuousData,
    criterion::Number;
    channel_out::Symbol = :is_extreme_value,
)
    is_extreme_value!(dat.data, dat.layout.label, criterion, channel_out = channel_out)
end


"""
    n_extreme_value(dat::DataFrame, columns::Vector{Symbol}, criterion::Number)::Int

Counts the number of extreme values in the specified columns.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data to check.
- `columns::Vector{Symbol}`: The columns to check for extreme values.
- `criterion::Float64`: The threshold for determining extreme values.

# Returns
An integer count of the number of extreme values found.

"""
function n_extreme_value(dat::DataFrame, columns::Vector{Symbol}, criterion::Number)::Int
    return sum(sum.(eachcol(abs.(select(dat, columns)) .>= criterion)))
end

function n_extreme_value(dat::DataFrame, columns::Symbol, criterion::Number)::Int
    return n_extreme_value(dat, [columns], criterion)
end

function n_extreme_value(dat::ContinuousData, columns::Vector{Symbol}, criterion::Number)::Int 
    return n_extreme_value(dat.data, columns, criterion)
end

function n_extreme_value(dat::ContinuousData, columns::Symbol, criterion::Number)::Int 
    return n_extreme_value(dat, [columns], criterion)
end

function n_extreme_value(dat::ContinuousData, criterion::Number)::Int
    return n_extreme_value(dat.data, dat.layout.label, criterion)
end





function channel_joint_probability(
    dat::DataFrame,
    channels::Vector{Symbol};
    threshold::Float64 = 5.0,
    normval::Int = 2,
    filter_samples = nothing,
)::DataFrame
    @info "channel_joint_probability: Computing probability for channels $(print_vector_(channels))"
    
    # Select the specified channels
    data = select(dat, channels)
    
    # Filter samples if requested
    if filter_samples !== nothing
        if filter_samples isa Symbol && hasproperty(dat, filter_samples)
            # If filter_samples is a column name, use that column
            data = data[dat[!, filter_samples], :]
        else
            # Otherwise assume it's a boolean vector
            data = data[filter_samples, :]
        end
    end
    
    # Convert to matrix and compute joint probability
    jp, indelec = joint_probability(Matrix(data)', threshold, normval)
    return DataFrame(channel = channels, jp = jp, rejection = indelec)
end

function channel_joint_probability(dat::ContinuousData; threshold::Float64 = 5.0, normval::Int = 2, filter_samples = nothing)::DataFrame
    return channel_joint_probability(dat.data, dat.layout.label; threshold = threshold, normval = normval, filter_samples = filter_samples)
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
