_first_df(dat::SingleDataFrameEeg) = dat.data
_first_df(dat::MultiDataFrameEeg) = isempty(dat.data) ? nothing : dat.data[1]

function _channel_cols(layout_channels::AbstractVector{Symbol}, df::DataFrame)::Vector{Symbol}
    lookup = DataFrames.index(df).lookup
    out = Vector{Symbol}(undef, length(layout_channels))
    k = 0
    @inbounds for ch in layout_channels
        if haskey(lookup, ch)
            dup = false
            for j = 1:k
                if out[j] === ch
                    dup = true
                    break
                end
            end
            if !dup
                k += 1
                out[k] = ch
            end
        end
    end
    resize!(out, k)
    return out
end

function _meta_cols(layout_channels::AbstractVector{Symbol}, df::DataFrame)::Vector{Symbol}
    idx = DataFrames.index(df)
    names_vec = idx.names
    isempty(layout_channels) && return copy(names_vec)
    lookup = idx.lookup
    first_idx = typemax(Int)
    @inbounds for ch in layout_channels
        pos = get(lookup, ch, typemax(Int))
        pos < first_idx && (first_idx = pos)
    end
    return first_idx == typemax(Int) ? copy(names_vec) : names_vec[1:(first_idx-1)]
end

function _extra_cols(layout_channels::AbstractVector{Symbol}, df::DataFrame)::Vector{Symbol}
    isempty(layout_channels) && return Symbol[]
    idx = DataFrames.index(df)
    lookup = idx.lookup
    last_idx = 0
    @inbounds for ch in layout_channels
        pos = get(lookup, ch, 0)
        pos > last_idx && (last_idx = pos)
    end
    names_vec = idx.names
    return (last_idx == 0 || last_idx >= length(names_vec)) ? Symbol[] : names_vec[(last_idx+1):end]
end

function _count_channels(layout_channels::AbstractVector{Symbol}, df::DataFrame)::Int
    lookup = DataFrames.index(df).lookup
    cnt = 0
    @inbounds for i in eachindex(layout_channels)
        ch = layout_channels[i]
        if haskey(lookup, ch)
            dup = false
            for j = 1:(i-1)
                if layout_channels[j] === ch
                    dup = true
                    break
                end
            end
            !dup && (cnt += 1)
        end
    end
    return cnt
end

function _get_cols_by_group(dat::EegData, group::Symbol)::Vector{Symbol}
    layout_channels = dat.layout.data[!, :label]
    df = _first_df(dat)
    if group === :channels
        return isnothing(df) ? copy(layout_channels) : _channel_cols(layout_channels, df)
    elseif group === :metadata
        return isnothing(df) ? Symbol[] : _meta_cols(layout_channels, df)
    elseif group === :extra
        return isnothing(df) ? Symbol[] : _extra_cols(layout_channels, df)
    end
    @minimal_error "Unknown group type: $group"
end

function _get_cols_by_group(dat::Vector{<:EegData}, group::Symbol)::Vector{Symbol}
    isempty(dat) && return Symbol[]
    return _get_cols_by_group(first(dat), group)
end
