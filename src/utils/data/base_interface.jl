# Partial array/collection interface for EegData:
# - `size`/`axes`/`ndims` report `(samples, channels)` or `(epochs, samples, channels)`.
# - `dat[...]` returns an `EegData` subset; use `Matrix(dat)` / `Array(dat)` for numeric arrays.
# - `MultiDataFrameEeg` iterates over its epoch `DataFrame`s (`length` = `n_epochs`).
# - `dat[:Fz]`, `:Fz in dat`, `haskey(dat, :Fz)`, and `keys(dat)` operate on channels.

Base.axes(dat::SingleDataFrameEeg) = (Base.OneTo(n_samples(dat)), Base.OneTo(n_channels(dat)))
Base.axes(dat::MultiDataFrameEeg) = (Base.OneTo(n_epochs(dat)), Base.OneTo(n_samples(dat)), Base.OneTo(n_channels(dat)))
Base.axes(dat::EegData, d::Int) = d <= length(axes(dat)) ? axes(dat)[d] : Base.OneTo(1)

Base.size(dat::SingleDataFrameEeg) = (n_samples(dat), n_channels(dat))
Base.size(dat::SingleDataFrameEeg, d::Int) = d == 1 ? n_samples(dat) : (d == 2 ? n_channels(dat) : 1)
Base.size(dat::MultiDataFrameEeg) = (n_epochs(dat), n_samples(dat), n_channels(dat))
Base.size(dat::MultiDataFrameEeg, d::Int) = d == 1 ? n_epochs(dat) : (d == 2 ? n_samples(dat) : (d == 3 ? n_channels(dat) : 1))

Base.ndims(::Type{<:SingleDataFrameEeg}) = 2
Base.ndims(::Type{<:MultiDataFrameEeg}) = 3
Base.ndims(dat::EegData) = ndims(typeof(dat))

Base.isempty(dat::SingleDataFrameEeg) = isempty(dat.data)
Base.isempty(dat::MultiDataFrameEeg) = isempty(dat.data)

Base.firstindex(::EegData, d::Int) = 1
Base.lastindex(dat::EegData, d::Int) = Base.size(dat, d)
Base.lastindex(dat::MultiDataFrameEeg) = length(dat.data)

Base.getindex(dat::SingleDataFrameEeg, time_idx, chan_idx) =
    subset(dat; sample_selection = _to_sample_func(time_idx), channel_selection = _to_channel_func(chan_idx))

Base.getindex(dat::MultiDataFrameEeg, epoch_idx, time_idx, chan_idx) = subset(
    dat;
    epoch_selection = _to_epoch_func(epoch_idx),
    sample_selection = _to_sample_func(time_idx),
    channel_selection = _to_channel_func(chan_idx),
)

function Base.setindex!(dat::SingleDataFrameEeg, val, time_idx, chan_idx)
    dat.data[time_idx, channel_labels(dat)[chan_idx]] = val
    return dat
end

function Base.setindex!(dat::MultiDataFrameEeg, val, epoch_idx, time_idx, chan_idx)
    dat.data[epoch_idx][time_idx, channel_labels(dat)[chan_idx]] = val
    return dat
end

Base.setindex!(dat::SingleDataFrameEeg, val, ci::CartesianIndex{2}) = (dat[ci[1], ci[2]] = val; dat)
Base.setindex!(dat::MultiDataFrameEeg, val, ci::CartesianIndex{3}) = (dat[ci[1], ci[2], ci[3]] = val; dat)
Base.getindex(dat::SingleDataFrameEeg, ci::CartesianIndex{2}) = dat[ci[1], ci[2]]
Base.getindex(dat::MultiDataFrameEeg, ci::CartesianIndex{3}) = dat[ci[1], ci[2], ci[3]]

Base.getindex(dat::EegData, chan::Union{Symbol,AbstractString}) = subset(dat; channel_selection = channels(Symbol(chan)))
Base.getindex(dat::EegData, chans::AbstractVector{<:Union{Symbol,AbstractString}}) =
    subset(dat; channel_selection = channels(Symbol.(chans)))

Base.iterate(dat::MultiDataFrameEeg, state = 1) = state > length(dat.data) ? nothing : (dat.data[state], state + 1)
Base.length(dat::MultiDataFrameEeg) = length(dat.data)
Base.eltype(::Type{<:MultiDataFrameEeg}) = DataFrame
Base.IteratorSize(::Type{<:MultiDataFrameEeg}) = Base.HasLength()

function Base.vcat(dats::EpochData...)
    isempty(dats) && @minimal_error("Cannot concatenate empty list of datasets")

    first_dat = dats[1]
    first_chs = channel_labels(first_dat)
    new_data = DataFrame[]
    epoch_counter = 1

    for dat in dats
        if dat.condition != first_dat.condition || dat.sample_rate != first_dat.sample_rate
            @minimal_error("Cannot vcat EpochData with different conditions or sample rates")
        end
        channel_labels(dat) != first_chs && @minimal_error("Cannot vcat EpochData with different channel labels")
        if !isempty(dat.data) && !isempty(first_dat.data) && n_samples(dat) != n_samples(first_dat)
            @minimal_error("Cannot vcat EpochData with different numbers of samples ($(n_samples(dat)) vs $(n_samples(first_dat)))")
        end

        for df in dat.data
            new_df = copy(df)
            hasproperty(new_df, :epoch) && (new_df.epoch .= epoch_counter)
            push!(new_data, new_df)
            epoch_counter += 1
        end
    end

    out_file = length(dats) > 1 ? "$(first_dat.file)_plus_$(length(dats)-1)" : first_dat.file
    return EpochData(
        out_file,
        first_dat.condition,
        first_dat.condition_name,
        new_data,
        copy(first_dat.layout),
        first_dat.sample_rate,
        copy(first_dat.analysis_info),
    )
end

@inline function _fill_col!(M::Matrix{Float64}, c_idx::Int, col::Vector{Float64}, n_samp::Int)
    @inbounds @simd for t = 1:n_samp
        M[t, c_idx] = col[t]
    end
end

@inline function _fill_col!(M::Matrix{Float64}, c_idx::Int, col::AbstractVector, n_samp::Int)
    @inbounds for t = 1:n_samp
        M[t, c_idx] = Float64(col[t])
    end
end

function Base.Matrix(dat::SingleDataFrameEeg)
    cols = channel_labels(dat)
    n_samp = n_samples(dat)
    M = Matrix{Float64}(undef, n_samp, length(cols))
    df = dat.data
    for (c_idx, c) in enumerate(cols)
        _fill_col!(M, c_idx, df[!, c], n_samp)
    end
    return M
end

Base.Array(dat::SingleDataFrameEeg) = Base.Matrix(dat)
Base.Array{Float64}(dat::SingleDataFrameEeg) = Base.Matrix(dat)
Base.Array{Float64}(dat::MultiDataFrameEeg) = Base.Array(dat)
Base.Array{Float64,2}(dat::SingleDataFrameEeg) = Base.Matrix(dat)
Base.Array{Float64,3}(dat::MultiDataFrameEeg) = Base.Array(dat)

@inline function _fill_slice!(A::Array{Float64,3}, e::Int, c_idx::Int, col::Vector{Float64}, n_samp::Int)
    @inbounds @simd for t = 1:n_samp
        A[e, t, c_idx] = col[t]
    end
end

@inline function _fill_slice!(A::Array{Float64,3}, e::Int, c_idx::Int, col::AbstractVector, n_samp::Int)
    @inbounds for t = 1:n_samp
        A[e, t, c_idx] = Float64(col[t])
    end
end

function Base.Array(dat::MultiDataFrameEeg)
    cols = channel_labels(dat)
    n_chan = length(cols)
    isempty(dat.data) && return Array{Float64,3}(undef, 0, 0, n_chan)
    n_ep = n_epochs(dat)
    n_samp = n_samples(dat)
    A = Array{Float64,3}(undef, n_ep, n_samp, n_chan)

    for (c_idx, c) in enumerate(cols)
        for e = 1:n_ep
            _fill_slice!(A, e, c_idx, dat.data[e][!, c], n_samp)
        end
    end
    return A
end

StatsBase.mean(dat::EpochData) = average_epochs(dat)
StatsBase.mean(dats::Vector{EpochData}) = average_epochs(dats)

Base.vec(dat::EegData) = vec(Base.Array(dat))
Base.parent(dat::SingleDataFrameEeg) = dat.data
Base.parent(dat::MultiDataFrameEeg) = dat.data

@inline _symbol_in(sym::Symbol, labels::AbstractVector{Symbol}) = sym in labels

function _has_channel(dat::EegData, sym::Symbol)::Bool
    df = _first_df(dat)
    labels = dat.layout.data[!, :label]
    isnothing(df) && return _symbol_in(sym, labels)
    return haskey(DataFrames.index(df).lookup, sym) && _symbol_in(sym, labels)
end

Base.in(chan::Union{Symbol,AbstractString}, dat::EegData) = _has_channel(dat, Symbol(chan))
Base.haskey(dat::EegData, chan::Union{Symbol,AbstractString}) = _has_channel(dat, Symbol(chan))
Base.keys(dat::EegData) = channel_labels(dat)

@generated function _fields_equal(a::T, b::T) where {T}
    checks = [:((getfield(a, $i) == getfield(b, $i)) === true) for i = 1:fieldcount(T)]
    return isempty(checks) ? :(true) : Expr(:&&, checks...)
end

@generated function _fields_isequal(a::T, b::T) where {T}
    checks = [:(isequal(getfield(a, $i), getfield(b, $i))) for i = 1:fieldcount(T)]
    return isempty(checks) ? :(true) : Expr(:&&, checks...)
end

@generated function _fields_hash(x::T, h::UInt) where {T}
    stmts = [:(h = hash($T, h))]
    for i = 1:fieldcount(T)
        push!(stmts, :(h = hash(getfield(x, $i), h)))
    end
    push!(stmts, :(return h))
    return Expr(:block, stmts...)
end

const _FieldwiseEq = Union{Layout,AnalysisInfo,Neighbors}

Base.:(==)(a::T, b::T) where {T<:_FieldwiseEq} = _fields_equal(a, b)
Base.isequal(a::T, b::T) where {T<:_FieldwiseEq} = _fields_isequal(a, b)
Base.hash(x::_FieldwiseEq, h::UInt) = _fields_hash(x, h)

Base.:(==)(a::T, b::T) where {T<:EegData} = _fields_equal(a, b)
Base.:(==)(::EegData, ::EegData) = false
Base.isequal(a::T, b::T) where {T<:EegData} = _fields_isequal(a, b)
Base.isequal(::EegData, ::EegData) = false
Base.hash(x::EegData, h::UInt) = _fields_hash(x, h)

function Base.isapprox(a::T, b::T; kwargs...) where {T<:SingleDataFrameEeg}
    channel_labels(a) == channel_labels(b) || return false
    size(a) == size(b) || return false
    a.sample_rate == b.sample_rate || return false
    return isapprox(Matrix(a), Matrix(b); kwargs...)
end

function Base.isapprox(a::T, b::T; kwargs...) where {T<:MultiDataFrameEeg}
    channel_labels(a) == channel_labels(b) || return false
    size(a) == size(b) || return false
    a.sample_rate == b.sample_rate || return false
    return isapprox(Array(a), Array(b); kwargs...)
end
Base.isapprox(::EegData, ::EegData; kwargs...) = false

Base.summary(
    dat::SingleDataFrameEeg,
) = "$(n_samples(dat)) timepoints × $(n_channels(dat)) channels $(nameof(typeof(dat))) (\"$(dat.file)\")"
Base.summary(
    dat::MultiDataFrameEeg,
) = "$(n_epochs(dat)) epochs × $(n_samples(dat)) timepoints × $(n_channels(dat)) channels $(nameof(typeof(dat))) (\"$(dat.file)\")"
