# === BRIDGED CHANNEL DETECTION ===

"""
    find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())
    find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real=0.99; channel_selection=channels())

Find channels that are highly correlated with each other, indicative of an electrolyte gel bridge.
Returns a `DataFrame` of the bridged pairs and their correlation coefficients.
"""
function find_bridged_channels(dat::SingleDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection = channels())
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)

    if n_ch < 2
        return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])
    end

    n_samp = n_samples(dat)
    X = Matrix{Float64}(undef, n_samp, n_ch)
    for (c, ch) in enumerate(selected_channels)
        col = dat.data[!, ch]
        if col isa Vector{Float64}
            col_f = col::Vector{Float64}
            @inbounds @simd for s = 1:n_samp
                X[s, c] = col_f[s]
            end
        else
            @inbounds for s = 1:n_samp
                X[s, c] = Float64(col[s])
            end
        end
    end

    R = cor(X; dims = 1)

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

function find_bridged_channels(dat::MultiDataFrameEeg, correlation_threshold::Real = 0.99; channel_selection = channels())
    # For epoched data, flatten across epochs and compute correlation via single matrix
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    n_ch = length(selected_channels)

    if n_ch < 2
        return DataFrame(channel_1 = Symbol[], channel_2 = Symbol[], correlation = Float64[])
    end

    n_total = sum(nrow, dat.data)
    X = Matrix{Float64}(undef, n_total, n_ch)

    row_offset = 0
    for epoch_df in dat.data
        n_ep = nrow(epoch_df)
        for (c, ch) in enumerate(selected_channels)
            col = epoch_df[!, ch]
            if col isa Vector{Float64}
                col_f = col::Vector{Float64}
                @inbounds @simd for s = 1:n_ep
                    X[row_offset+s, c] = col_f[s]
                end
            else
                @inbounds for s = 1:n_ep
                    X[row_offset+s, c] = Float64(col[s])
                end
            end
        end
        row_offset += n_ep
    end

    R = cor(X; dims = 1)

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

function find_bridged_channels(dat::AbstractVector{<:EegData}, correlation_threshold::Real = 0.99; channel_selection = channels())
    return [find_bridged_channels(d, correlation_threshold; channel_selection = channel_selection) for d in dat]
end
