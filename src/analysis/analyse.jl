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
    return ContinuousData(create_eeg_dataframe(dat), DataFrame(CSV.File(layout_file_name)), dat.header.sample_rate[1], AnalysisInfo())
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
    # Initialize with default AnalysisInfo
    return ContinuousData(create_eeg_dataframe(dat), layout, dat.header.sample_rate[1], AnalysisInfo())
end



function channel_summary(dat::DataFrame, channel_labels::Vector{Symbol})::DataFrame

    # Select the specified channels
    selected_data = select(dat, channel_labels)

    # Initialize a matrix to store summary statistics
    summary_stats = Matrix{Float64}(undef, length(channel_labels), 6)  # 6 statistics: min, max, range, std, mad, var

    # Compute summary statistics for each column
    for (i, col) in enumerate(channel_labels)
        col_data = @view selected_data[!, col]
        summary_stats[i, :] = [minimum(col_data), maximum(col_data), datarange(col_data), std(col_data), mad(col_data), var(col_data)]
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

function channel_summary(dat::SingleDataFrameEeg)::DataFrame
    return channel_summary(dat.data, dat.layout.label)
end

function channel_summary(dat::SingleDataFrameEeg, channel_numbers::Union{Vector{Int}, UnitRange})::DataFrame
    channel_labels = channel_number_to_channel_label(dat.layout.label, channel_numbers)
    return channel_summary(dat.data, channel_labels)
end

function channel_summary(dat::SingleDataFrameEeg, channel_labels::Vector{Symbol})::DataFrame
    return channel_summary(dat.data, channel_labels)
end



function channel_summary(dat::MultiDataFrameEeg)::Vector{DataFrame}
    return [channel_summary(dat.data[trial], dat.layout.label) for trial in eachindex(dat.data)]
end

function channel_summary(dat::MultiDataFrameEeg, channel_numbers::Union{Vector{Int}, UnitRange})::Vector{DataFrame}
    channel_labels = channel_number_to_channel_label(dat.layout.label, channel_numbers)
    return [channel_summary(dat.data[trial], channel_labels) for trial in eachindex(dat.data)]
end

function channel_summary(dat::MultiDataFrameEeg, channel_labels::Vector{Symbol})::Vector{DataFrame}
    return [channel_summary(dat.data[trial], channel_labels) for trial in eachindex(dat.data)]
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
function correlation_matrix(dat::DataFrame, channel_labels::Vector{Symbol})::DataFrame
    df = DataFrame(cor(Matrix(select(dat, channel_labels))), channel_labels)
    insertcols!(df, 1, :row => layout.label)
    return df
end

function correlation_matrix(dat::ContinuousData)::DataFrame
    return correlation_matrix(dat.data, dat.layout.label)
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
function detect_eog_onsets!(
    dat::ContinuousData,
    criterion::Real,
    channel_in::Symbol,
    channel_out::Symbol,
)
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
function is_extreme_value( dat::DataFrame, columns::Vector{Symbol}, criterion::Number,)::Vector{Bool}
    return any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)), dims = 2)[:]
end

function is_extreme_value!(dat::DataFrame, columns::Vector{Symbol}, criterion::Number; channel_out::Symbol = :is_extreme_value)
    dat[!, channel_out] .= any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)), dims = 2)[:]
end

function is_extreme_value!(dat::ContinuousData, columns::Vector{Symbol}, criterion::Number; channel_out::Symbol = :is_extreme_value)
    is_extreme_value!(dat.data, columns, criterion, channel_out = channel_out)
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
    return sum(abs.(select(dat, columns)) .>= criterion)
end



function channel_joint_probability(dat::DataFrame, channels::Vector{Symbol}; threshold::Float64 = 5.0, normval::Int = 2)::DataFrame
    @info "channel_joint_probability: Computing probability for channels $(print_vector_(channels))"
    data = view(Matrix(select(dat, channels)), :, :) # Use view instead of copying
    jp, indelec = joint_probability(data', threshold, normval)
    return DataFrame(channel = channels, jp = jp, rejection = indelec)
end

function channel_joint_probability(dat::ContinuousData; threshold::Float64 = 5.0, normval::Int = 2)::DataFrame
    return channel_joint_probability(dat.data, dat.layout.label; threshold=threshold, normval=normval)
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

