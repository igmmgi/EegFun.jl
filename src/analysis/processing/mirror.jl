@inline _validate_mirror_side(side::Symbol) =
    side in (:pre, :post, :both) || @minimal_error("side must be :pre, :post, or :both, got :$side")

"""
    mirror!(dat::EegData, side::Symbol = :both)
    mirror!(data_vec::Vector{<:EegData}, side::Symbol = :both)

Mirror data in-place by reflecting time-reversed copies at boundaries (`:pre`, `:post`, or `:both`)
to suppress edge artifacts during filtering. Call [`unmirror!`](@ref) with the same `side` afterward.
"""
function mirror!(dat::SingleDataFrameEeg, side::Symbol = :both)
    _validate_mirror_side(side)
    @info "Mirroring data on side: $side"
    dat.data = _mirror_dataframe!(dat.data, side)
    return dat
end

function mirror!(dat::EpochData, side::Symbol = :both)
    _validate_mirror_side(side)
    @info "Mirroring epoched data on side: $side"
    n_ep = length(dat.data)
    if n_ep > 1
        Threads.@threads :static for i = 1:n_ep
            dat.data[i] = _mirror_dataframe!(dat.data[i], side)
        end
    elseif n_ep == 1
        dat.data[1] = _mirror_dataframe!(dat.data[1], side)
    end
    return dat
end

function mirror!(data_vec::Vector{<:EegData}, side::Symbol = :both)
    for dat in data_vec
        mirror!(dat, side)
    end
    return data_vec
end

"""
    mirror(dat::EegData, side::Symbol = :both)
    mirror(data_vec::Vector{<:EegData}, side::Symbol = :both)

Non-mutating version of [`mirror!`](@ref).
"""
mirror(dat::EegData, side::Symbol = :both) = mirror!(copy(dat), side)
mirror(data_vec::Vector{<:EegData}, side::Symbol = :both) = [mirror(dat, side) for dat in data_vec]

"""
    unmirror!(dat::EegData, side::Symbol = :both)
    unmirror!(data_vec::Vector{<:EegData}, side::Symbol = :both)

Remove mirrored boundary segments in-place, restoring original length. `side` must match [`mirror!`](@ref).
"""
function unmirror!(dat::SingleDataFrameEeg, side::Symbol = :both)
    _validate_mirror_side(side)
    @info "Unmirroring data on side: $side"
    dat.data = _unmirror_dataframe!(dat.data, side)
    return dat
end

function unmirror!(dat::EpochData, side::Symbol = :both)
    _validate_mirror_side(side)
    @info "Unmirroring epoched data on side: $side"
    n_ep = length(dat.data)
    if n_ep > 1
        Threads.@threads :static for i = 1:n_ep
            dat.data[i] = _unmirror_dataframe!(dat.data[i], side)
        end
    elseif n_ep == 1
        dat.data[1] = _unmirror_dataframe!(dat.data[1], side)
    end
    return dat
end

function unmirror!(data_vec::Vector{<:EegData}, side::Symbol = :both)
    for dat in data_vec
        unmirror!(dat, side)
    end
    return data_vec
end

"""
    unmirror(dat::EegData, side::Symbol = :both)
    unmirror(data_vec::Vector{<:EegData}, side::Symbol = :both)

Non-mutating version of [`unmirror!`](@ref).
"""
unmirror(dat::EegData, side::Symbol = :both) = unmirror!(copy(dat), side)
unmirror(data_vec::Vector{<:EegData}, side::Symbol = :both) = [unmirror(dat, side) for dat in data_vec]

function _mirror_column(src::AbstractVector{T}, n::Int, pre_len::Int, post_len::Int) where {T}
    n_new = pre_len + n + post_len
    dst = Vector{T}(undef, n_new)
    if pre_len > 0
        @inbounds for i = 1:pre_len
            dst[i] = src[n-i+1]
        end
    end
    @inbounds copyto!(dst, pre_len + 1, src, 1, n)
    if post_len > 0
        base = pre_len + n
        @inbounds for i = 1:post_len
            dst[base+i] = src[n-i]
        end
    end
    return dst
end

function _mirror_time_column(t::AbstractVector, n::Int, pre_len::Int, post_len::Int)
    n_new = pre_len + n + post_len
    t1 = Float64(t[1])
    tn = Float64(t[n])
    dt = (tn - t1) / (n - 1)
    dst = Vector{Float64}(undef, n_new)
    if pre_len > 0
        @inbounds @simd for i = 1:pre_len
            dst[i] = t1 - (pre_len - i + 1) * dt
        end
    end
    @inbounds for i = 1:n
        dst[pre_len+i] = Float64(t[i])
    end
    if post_len > 0
        base = pre_len + n
        @inbounds @simd for i = 1:post_len
            dst[base+i] = tn + i * dt
        end
    end
    return dst
end

function _mirror_dataframe!(df::DataFrame, side::Symbol)
    n = nrow(df)
    n <= 1 && return copy(df)
    pre_len = side == :post ? 0 : (n - 1)
    post_len = side == :pre ? 0 : (n - 1)

    col_names = propertynames(df)
    out_cols = Vector{AbstractVector}(undef, length(col_names))
    for (idx, col) in enumerate(col_names)
        src = df[!, idx]
        out_cols[idx] = col == :time ? _mirror_time_column(src, n, pre_len, post_len) : _mirror_column(src, n, pre_len, post_len)
    end
    return DataFrame(out_cols, col_names; copycols = false)
end

function _unmirror_dataframe!(df::DataFrame, side::Symbol)
    n = nrow(df)
    if side == :pre
        orig_len = div(n + 1, 2)
        rng = orig_len:n
    elseif side == :post
        orig_len = div(n + 1, 2)
        rng = 1:orig_len
    else
        orig_len = div(n + 2, 3)
        rng = orig_len:(2*orig_len-1)
    end
    col_names = propertynames(df)
    out_cols = AbstractVector[df[!, idx][rng] for idx in eachindex(col_names)]
    return DataFrame(out_cols, col_names; copycols = false)
end
