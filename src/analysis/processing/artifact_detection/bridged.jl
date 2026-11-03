function _bridged_pairs(X::Matrix{Float64}, selected_channels::Vector{Symbol}, correlation_threshold::Real)
    R = cor(X; dims = 1)
    n_ch = length(selected_channels)
    pairs_1 = Symbol[]
    pairs_2 = Symbol[]
    corrs = Float64[]

    for i = 1:(n_ch-1)
        ch1 = selected_channels[i]
        for j = (i+1):n_ch
            r = R[i, j]
            if r > correlation_threshold
                push!(pairs_1, ch1)
                push!(pairs_2, selected_channels[j])
                push!(corrs, r)
            end
        end
    end
    return DataFrame(channel_1 = pairs_1, channel_2 = pairs_2, correlation = corrs)
end

"""
    find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())
    find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())

Find channels that are highly correlated with each other, indicative of an electrolyte gel bridge.
Returns a `DataFrame` of the bridged pairs and their correlation coefficients.
"""
function find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection = channels())
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)
    n_ch < 2 && return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])

    n_samp = n_samples(dat)
    X = Matrix{Float64}(undef, n_samp, n_ch)
    for (c, ch) in enumerate(selected_channels)
        copyto!(view(X, :, c), dat.data[!, ch])
    end
    return _bridged_pairs(X, selected_channels, correlation_threshold)
end

function find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection = channels())
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)
    n_ch < 2 && return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])

    X = Matrix{Float64}(undef, sum(nrow, dat.data), n_ch)
    row_offset = 0
    for epoch_df in dat.data
        n_ep = nrow(epoch_df)
        rows = (row_offset+1):(row_offset+n_ep)
        for (c, ch) in enumerate(selected_channels)
            copyto!(view(X, rows, c), epoch_df[!, ch])
        end
        row_offset += n_ep
    end
    return _bridged_pairs(X, selected_channels, correlation_threshold)
end

function find_bridged_channels(dat::AbstractVector{<:EegData}, correlation_threshold::Real = 0.99; channel_selection = channels())
    return [find_bridged_channels(d, correlation_threshold; channel_selection = channel_selection) for d in dat]
end
