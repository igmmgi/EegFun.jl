# === Base Array Interface ===
Base.IndexStyle(::Type{<:EegData}) = IndexCartesian()

# size & axes & dims
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

Base.eltype(::Type{<:EegData}) = Float64
Base.eltype(dat::EegData) = Float64

Base.isempty(dat::SingleDataFrameEeg) = isempty(dat.data)
Base.isempty(dat::MultiDataFrameEeg) = isempty(dat.data)

Base.firstindex(::EegData, d::Int) = 1
Base.firstindex(::EegData) = 1
Base.lastindex(dat::EegData, d::Int) = Base.size(dat, d)
Base.lastindex(dat::EegData) = prod(size(dat))



# getindex implementation
function Base.getindex(dat::SingleDataFrameEeg, time_idx, chan_idx)
    return subset(dat, sample_selection = _to_sample_func(time_idx), channel_selection = _to_channel_func(chan_idx))
end

function Base.getindex(dat::MultiDataFrameEeg, epoch_idx, time_idx, chan_idx)
    return subset(
        dat,
        epoch_selection = _to_epoch_func(epoch_idx),
        sample_selection = _to_sample_func(time_idx),
        channel_selection = _to_channel_func(chan_idx),
    )
end

# setindex! implementation
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

# CartesianIndex getindex
Base.getindex(dat::SingleDataFrameEeg, ci::CartesianIndex{2}) = dat[ci[1], ci[2]]
Base.getindex(dat::MultiDataFrameEeg, ci::CartesianIndex{3}) = dat[ci[1], ci[2], ci[3]]

# Shorthand channel subsetting
Base.getindex(dat::EegData, chan::Union{Symbol,AbstractString}) = subset(dat, channel_selection = channels(Symbol(chan)))
Base.getindex(dat::EegData, chans::AbstractVector{<:Union{Symbol,AbstractString}}) =
    subset(dat, channel_selection = channels(Symbol.(chans)))
# === Base Iteration Interface ===
Base.iterate(dat::MultiDataFrameEeg, state = 1) = state > length(dat.data) ? nothing : (dat.data[state], state + 1)
Base.length(dat::MultiDataFrameEeg) = length(dat.data)

# === Base Concatenation Interface ===
function Base.vcat(dats::EpochData...)
    if isempty(dats)
        throw(ArgumentError("Cannot concatenate empty list of datasets"))
    end

    first_dat = dats[1]
    new_data = DataFrame[]
    epoch_counter = 1

    for dat in dats
        if dat.condition != first_dat.condition || dat.sample_rate != first_dat.sample_rate
            throw(ArgumentError("Cannot vcat EpochData with different conditions or sample rates"))
        end

        for df in dat.data
            new_df = copy(df)
            if hasproperty(new_df, :epoch)
                new_df.epoch .= epoch_counter
            end
            push!(new_data, new_df)
            epoch_counter += 1
        end
    end

    # Just use the first file name with a plus to avoid absurdly long file names
    file_name = length(dats) > 1 ? "$(first_dat.file)_plus_$(length(dats)-1)" : first_dat.file

    return EpochData(
        file_name,
        first_dat.condition,
        first_dat.condition_name,
        new_data,
        copy(first_dat.layout),
        first_dat.sample_rate,
        copy(first_dat.analysis_info),
    )
end

# === Base Tensor Conversion Interface ===
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
    n_chan = length(cols)
    M = Matrix{Float64}(undef, n_samp, n_chan)
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
    n_ep = n_epochs(dat)
    n_samp = n_samples(dat)
    n_chan = length(cols)

    # Preallocate 3D array: (epochs x samples x channels)
    # matching the size(dat) output of (epochs, samples, channels)
    A = Array{Float64,3}(undef, n_ep, n_samp, n_chan)

    for e = 1:n_ep
        df = dat.data[e]
        for (c_idx, c) in enumerate(cols)
            _fill_slice!(A, e, c_idx, df[!, c], n_samp)
        end
    end

    return A
end

# === StatsBase Integration ===
StatsBase.mean(dat::EpochData) = average_epochs(dat)
StatsBase.mean(dats::Vector{EpochData}) = average_epochs(dats)

# === Base Type Conversion Interface ===
Base.convert(::Type{Matrix}, dat::SingleDataFrameEeg) = Base.Matrix(dat)
Base.convert(::Type{Matrix{Float64}}, dat::SingleDataFrameEeg) = Base.Matrix(dat)
Base.convert(::Type{Array}, dat::EegData) = Base.Array(dat)
Base.convert(::Type{DataFrame}, dat::EegData) = DataFrame(dat)
Base.vec(dat::EegData) = vec(Base.Array(dat))
Base.parent(dat::SingleDataFrameEeg) = dat.data
Base.parent(dat::MultiDataFrameEeg) = dat.data


# === Channel Membership & Keys ===
Base.in(chan::Union{Symbol,AbstractString}, dat::EegData) = Symbol(chan) in channel_labels(dat)
Base.haskey(dat::EegData, chan::Union{Symbol,AbstractString}) = Symbol(chan) in channel_labels(dat)
Base.keys(dat::EegData) = channel_labels(dat)

# === Equality and Approximations ===
function Base.:(==)(a::Layout, b::Layout)
    for f in fieldnames(Layout)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

function Base.:(==)(a::AnalysisInfo, b::AnalysisInfo)
    for f in fieldnames(AnalysisInfo)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

function Base.:(==)(a::T, b::T) where {T<:EegData}
    for f in fieldnames(T)
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end
Base.:(==)(::EegData, ::EegData) = false

Base.isequal(a::T, b::T) where {T<:EegData} = a == b

function Base.isapprox(a::T, b::T; kwargs...) where {T<:SingleDataFrameEeg}
    channel_labels(a) == channel_labels(b) || return false
    size(a) == size(b) || return false
    return isapprox(Matrix(a), Matrix(b); kwargs...)
end

function Base.isapprox(a::T, b::T; kwargs...) where {T<:MultiDataFrameEeg}
    channel_labels(a) == channel_labels(b) || return false
    size(a) == size(b) || return false
    return isapprox(Array(a), Array(b); kwargs...)
end
Base.isapprox(::EegData, ::EegData; kwargs...) = false

# === Base.summary ===
Base.summary(
    dat::SingleDataFrameEeg,
) = "$(n_samples(dat)) timepoints × $(n_channels(dat)) channels $(nameof(typeof(dat))) (\"$(dat.file)\")"
Base.summary(
    dat::MultiDataFrameEeg,
) = "$(n_epochs(dat)) epochs × $(n_samples(dat)) timepoints × $(n_channels(dat)) channels $(nameof(typeof(dat))) (\"$(dat.file)\")"
