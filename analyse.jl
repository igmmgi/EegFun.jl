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
        DataFrame(time = data.time, triggers = data.triggers.raw),
        DataFrame(data.data, Symbol.(data.header.channel_labels[1:end-1])),
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

"""
    channel_summary(dat::DataFrame, channel_labels::Vector{<:AbstractString})::DataFrame

Calculates summary statistics for specified channels in the EEG data.

# Arguments
- `dat::DataFrame`: The DataFrame containing EEG data.
- `channel_labels::Vector{<:AbstractString}`: A collection of channel labels to summarize.

# Returns
A DataFrame containing summary statistics (min, max, range, std, mad) for each specified channel.

"""
function channel_summary(dat::DataFrame, channel_labels::Vector{<:AbstractString})::DataFrame
    selected_data = select(dat, channel_labels)
    summary_df = combine(selected_data, names(selected_data) .=> [
        minimum, maximum, x -> datarange(x), std, mad
    ] .=> [:min, :max, :range, :std, :mad])
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
function correlation_matrix(dat::DataFrame, layout::DataFrame)::DataFrame
    labels = layout.label
    corr_matrix = cor(Matrix(select(dat, labels)))
    return DataFrame(labels, corr_matrix, :auto)
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
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Union{Symbol, <:AbstractString}, channel_out::Union{Symbol, String})
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
function is_extreme_value(dat::DataFrame, columns::Union{Vector{Symbol}, Vector{<:AbstractString}}, criterion::Number)::Bool
    return any(x -> abs.(x) >= criterion, Matrix(select(dat, columns)))
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
function n_extreme_value(dat::DataFrame, columns::Union{Vector{Symbol}, Vector{<:AbstractString}}, criterion::Number)::Int
    return sum(abs.(select(dat, columns)) .>= criterion)
end

