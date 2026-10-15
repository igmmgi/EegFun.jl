# === Base Array Interface ===
# size
Base.size(dat::SingleDataFrameEeg) = (n_samples(dat), length(channel_labels(dat)))
Base.size(dat::SingleDataFrameEeg, d::Int) = d == 1 ? n_samples(dat) : (d == 2 ? length(channel_labels(dat)) : 1)
Base.size(dat::MultiDataFrameEeg) = (n_epochs(dat), n_samples(dat), length(channel_labels(dat)))
Base.size(dat::MultiDataFrameEeg, d::Int) = d == 1 ? n_epochs(dat) : (d == 2 ? n_samples(dat) : (d == 3 ? length(channel_labels(dat)) : 1))

Base.lastindex(dat::EegData, d::Int) = Base.size(dat, d)

# getindex helpers
_to_sample_func(::Colon) = samples()
_to_sample_func(f::Function) = f
_to_sample_func(idx::Integer) = x -> [i == idx for i = 1:nrow(x)]
_to_sample_func(idx) = x -> [i in idx for i = 1:nrow(x)]

_to_channel_func(::Colon) = channels()
_to_channel_func(f::Function) = f
_to_channel_func(idx::AbstractString) = channels(Symbol(idx))
_to_channel_func(idx::AbstractVector{<:AbstractString}) = channels(Symbol.(idx))
_to_channel_func(idx) = channels(idx)

_to_epoch_func(::Colon) = epochs()
_to_epoch_func(f::Function) = f
_to_epoch_func(idx) = epochs(idx)

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
            if "epoch" in names(new_df)
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
function Base.Matrix(dat::SingleDataFrameEeg)
    cols = channel_labels(dat)
    return Matrix(dat.data[!, cols])
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
        A[e, :, :] .= Matrix(dat.data[e][!, cols])
    end

    return A
end

# === StatsBase Integration ===
StatsBase.mean(dat::EpochData) = average_epochs(dat)
StatsBase.mean(dats::Vector{EpochData}) = average_epochs(dats)
