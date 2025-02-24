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
        DataFrame(Float64.(data.data), Symbol.(data.header.channel_labels[1:end-1])),
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
    return ContinuousData(create_eeg_dataframe(dat), DataFrame(CSV.File(layout_file_name)), dat.header.sample_rate[1])
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
    return ContinuousData(create_eeg_dataframe(dat), layout, dat.header.sample_rate[1])
end



function channel_summary(dat::DataFrame, channel_labels::Vector{<:AbstractString})
    # Select the specified channels
    selected_data = select(dat, channel_labels)

    # Initialize a matrix to store summary statistics
    summary_stats = Matrix{Float64}(undef, length(channel_labels), 6)  # 6 statistics: min, max, range, std, mad, var

    # Compute summary statistics for each column
    for (i, col) in enumerate(channel_labels)
        col_data = selected_data[!, col]
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





"""
    correlation_matrix(dat::DataFrame, layout::DataFrame)::DataFrame

Calculates the correlation matrix for the EEG data.

# Arguments
- `dat::DataFrame`: The DataFrame containing EEG data.
- `layout::DataFrame`: The DataFrame containing layout information.

# Returns
A DataFrame containing the correlation matrix of the specified channels.

"""
function correlation_matrix(dat::DataFrame, channel_labels)::DataFrame
    df = DataFrame(cor(Matrix(select(dat, channel_labels))), channel_labels)
    insertcols!(df, 1, :row => layout.label)
    return df
end


"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Float64, channel_in::Symbol, channel_out::Symbol)

Detects EOG (electrooculogram) onsets in the EEG data based on a specified criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Real`: The threshold for detecting EOG onsets.
- `channel_in::Union{Symbol, String}`: The channel from which to detect EOG onsets.
- `channel_out::Union{Symbol, String}`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event. This parameter can accept both `Symbol` and `String` types.
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Float64`: The threshold for detecting EOG onsets.
- `channel_in::Symbol`: The channel from which to detect EOG onsets.
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event. This parameter can accept both `Symbol` and `String` types.

# Returns
Nothing. The function modifies the input data in place.

"""
function detect_eog_onsets!(
    dat::ContinuousData,
    criterion::Real,
    channel_in::Union{Symbol,<:AbstractString},
    channel_out::Union{Symbol,String},
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
    is_extreme_value(dat::DataFrame, columns::Union{Vector{Symbol}, Vector{<:AbstractString}}, criterion::Real)::Bool

Checks if any values in the specified columns exceed a given criterion.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data to check.
- `columns::Vector{Symbol}`: The columns to check for extreme values.
- `criterion::Float64`: The threshold for determining extreme values.

# Returns
A Boolean indicating whether any extreme values were found.

"""
function is_extreme_value(
    dat::DataFrame,
    columns::Union{Vector{Symbol},Vector{<:AbstractString}},
    criterion::Number,
)::Vector{Bool}
    return any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)), dims = 2)[:]
end

function is_extreme_value!(dat::DataFrame, columns::Union{Vector{Symbol},Vector{<:AbstractString}}, criterion::Number)
    dat[!, "is_extreme"] .= any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)), dims = 2)[:]
end

function is_extreme_value!(dat::ContinuousData, columns::Union{Vector{Symbol},Vector{<:AbstractString}}, criterion::Number)
    is_extreme_value!(dat.data, columns, criterion)
end







"""
    n_extreme_value(dat::DataFrame, columns::Union{Vector{Symbol}, Vector{<:AbstractString}}, criterion::Number)::Int

Counts the number of extreme values in the specified columns.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data to check.
- `columns::Vector{Symbol}`: The columns to check for extreme values.
- `criterion::Float64`: The threshold for determining extreme values.

# Returns
An integer count of the number of extreme values found.

"""
function n_extreme_value(
    dat::DataFrame,
    columns::Union{Vector{Symbol},Vector{<:AbstractString}},
    criterion::Number,
)::Int
    return sum(abs.(select(dat, columns)) .>= criterion)
end



function channel_joint_probability(dat::DataFrame, channels, threshold::Float64 = 5.0, normval::Int = 2)

    println("Info (pmd_badChanJointProb): Computing probability for channels...")

    # Call the jointprob function
    data = Float64.(Matrix(select(dat, channels)))

    jp, indelec = joint_probability(permutedims(data), threshold, normval)
    #badElectrodes = findall(indelec)

    summary_df = DataFrame(channel = channels, jp = jp, rejection = indelec)
    return summary_df
end

function joint_probability(signal::Matrix{Float64}, threshold::Float64, normalize::Int, discret::Int = 1000)
    nbchan = size(signal)[1]
    jp = zeros(nbchan)
    for rc = 1:nbchan
        # Compute the density function
        dataProba, _ = realproba(signal[rc, :], discret)
        # Compute joint probability
        jp[rc] = -sum(log.(dataProba))
    end
    # Normalize the joint probability
    if normalize != 0
        tmpjp = jp
        if normalize == 2
            tmpjp = sort(jp)
            tmpjp = tmpjp[round(Int, length(tmpjp) * 0.1)+1:end-round(Int, length(tmpjp) * 0.1)]
        end
        jp = (jp .- mean(tmpjp)) ./ std(tmpjp)
    end
    # Reject channels based on threshold
    if threshold != 0
        rej = abs.(jp) .> threshold
    else
        rej = zeros(Bool, size(jp))
    end
    return jp, rej
end

function realproba(data::Vector{Float64}, bins::Int = 1000)
    if bins > 0
        # Compute the density function
        minimum_value = minimum(data)
        maximum_value = maximum(data)
        data = floor.((data .- minimum_value) ./ (maximum_value - minimum_value) .* (bins - 1)) .+ 1
        if any(isnan.(data))
            @warn "Binning failed - could be due to zeroed out channel"
        end
        # Compute histogram
        sortbox = zeros(bins)
        for d in data
            if !isnan(d)
                sortbox[Int(d)] += 1
            end
        end
        probaMap = sortbox[Int.(data)] / length(data)
        sortbox = sortbox / length(data)
    else
        # Base over error function (Gaussian approximation)
        data = (data .- mean(data)) ./ std(data)
        probaMap = exp.(-0.5 .* (data .* data)) / (2 * π)
        probaMap = probaMap / sum(probaMap)
        sortbox = probaMap / sum(probaMap)
    end
    return probaMap, sortbox
end

