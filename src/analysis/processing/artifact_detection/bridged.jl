# === BRIDGED CHANNEL DETECTION ===

"""
    find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())
    find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())

Find channels that are highly correlated with each other, indicative of an electrolyte gel bridge.
Returns a `DataFrame` of the bridged pairs and their correlation coefficients.
"""
function find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection::Function = channels())
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)

    pairs = Tuple{Symbol,Symbol}[]
    correlations = Float64[]

    if n_ch < 2
        return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])
    end

    # Calculate pairwise correlation
    for i = 1:(n_ch-1)
        for j = (i+1):n_ch
            ch1 = selected_channels[i]
            ch2 = selected_channels[j]
            corr = cor(dat.data[!, ch1], dat.data[!, ch2])
            if corr > correlation_threshold
                push!(pairs, (ch1, ch2))
                push!(correlations, corr)
            end
        end
    end

    return DataFrame(channel_1 = first.(pairs), channel_2 = last.(pairs), correlation = correlations)
end

function find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection::Function = channels())
    # For epoched data, we flatten the data vertically and compute correlation across all epochs
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)

    if n_ch < 2
        return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])
    end

    # Extract channel data across all epochs
    flattened_data = Dict{Symbol,Vector{Float64}}()
    for ch in selected_channels
        flattened_data[ch] = Float64[]
    end

    for epoch_df in dat
        for ch in selected_channels
            append!(flattened_data[ch], epoch_df[!, ch])
        end
    end

    pairs = Tuple{Symbol,Symbol}[]
    correlations = Float64[]

    for i = 1:(n_ch-1)
        for j = (i+1):n_ch
            ch1 = selected_channels[i]
            ch2 = selected_channels[j]
            corr = cor(flattened_data[ch1], flattened_data[ch2])
            if corr > correlation_threshold
                push!(pairs, (ch1, ch2))
                push!(correlations, corr)
            end
        end
    end

    return DataFrame(channel_1 = first.(pairs), channel_2 = last.(pairs), correlation = correlations)
end

