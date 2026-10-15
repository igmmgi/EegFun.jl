# Automatic epoch rejection based on statistical criteria.
#
# Automatically rejects epochs based on statistical measures 
# (variance, max, min, absolute max, range, kurtosis) using z-score thresholds. 
# This is useful for removing epochs with artifacts without manual inspection.
"""
    Rejection(channel::Symbol, epoch::Int)

Struct representing a rejected epoch for a specific channel.
"""
struct Rejection
    channel::Symbol
    epoch::Int
end

"""Pretty-print a `Rejection`, showing channel and epoch."""
Base.show(io::IO, r::Rejection) = print(io, "Rejection(:$(r.channel), $(r.epoch))")

"""Check whether two `Rejection` objects refer to the same channel and epoch."""
is_equal_rejection(a::Rejection, b::Rejection) = a.channel == b.channel && a.epoch == b.epoch

"""De-duplicate a vector of `Rejection`s by `(channel, epoch)`."""
function unique_rejections(rejections::Vector{Rejection})
    seen = Set{Tuple{Symbol,Int}}()
    out = Rejection[]
    for rejection in rejections
        key = (rejection.channel, rejection.epoch)
        if key ∉ seen
            push!(seen, key)
            push!(out, rejection)
        end
    end
    return out
end

"""Return the distinct channel symbols across a vector of rejections."""
unique_channels(rejections::Vector{Rejection}) = unique(map(x -> x.channel, rejections))
"""Return the distinct epoch indices across a vector of rejections."""
unique_epochs(rejections::Vector{Rejection}) = unique(map(x -> x.epoch, rejections))


"""
    EpochInfo

Stores metadata about the epoch condition.

# Fields
- `number::Int`: Condition number from the epoch data
- `name::String`: Condition name from the epoch data
- `n::Int`: Number of epochs before rejection
"""
struct EpochInfo
    number::Int
    name::String
    n::Int
end

"""
    ZScoreRejectionInfo

Stores z-score based rejection information.

# Fields
- `z_measures::Vector{Symbol}`: Which z-score measures were evaluated
- `z_variance::Vector{Rejection}`: Rejections due to high variance
- `z_max::Vector{Rejection}`: Rejections due to high maximum values
- `z_min::Vector{Rejection}`: Rejections due to low minimum values
- `z_abs::Vector{Rejection}`: Rejections due to high absolute values
- `z_range::Vector{Rejection}`: Rejections due to large range
- `z_kurtosis::Vector{Rejection}`: Rejections due to high kurtosis
"""
struct ZScoreRejectionInfo
    z_measures::Vector{Symbol}
    z_variance::Vector{Rejection}
    z_max::Vector{Rejection}
    z_min::Vector{Rejection}
    z_abs::Vector{Rejection}
    z_range::Vector{Rejection}
    z_kurtosis::Vector{Rejection}
end

"""
    EpochRejectionInfo

Stores information about which epochs were rejected and why, and optionally tracks channel repairs.

# Fields
- `name::String`: Name/identifier for this rejection info (e.g., "rejection_step1", "rejection_step2")
- `info::EpochInfo`: Condition metadata (number, name, n_epochs)
- `n_artifacts::Int`: Total number of artifact detections (channel-epoch pairs)
- `abs_criterion::Real`: Absolute voltage threshold (μV) used for rejection
- `abs_rejections::Union{Vector{Rejection}, Nothing}`: Rejections due to absolute voltage threshold (Nothing if abs_criterion = 0)
- `z_criterion::Real`: Z-score criterion used for rejection
- `z_rejections::Union{ZScoreRejectionInfo, Nothing}`: Z-score based rejection info (Nothing if z_criterion = 0)
- `rejected::Vector{Rejection}`: All rejected epochs
- `repaired::Union{OrderedDict{Int, Vector{Symbol}}, Nothing}`: Channels repaired per epoch, ordered by epoch number (populated during repair, Nothing if no repairs)
- `skipped::Union{OrderedDict{Int, Vector{Symbol}}, Nothing}`: Channels skipped per epoch, ordered by epoch number (populated during repair, Nothing if no repairs)
"""
mutable struct EpochRejectionInfo
    name::String
    info::EpochInfo
    n_artifacts::Int
    abs_criterion::Real
    abs_rejections::Union{Vector{Rejection},Nothing}
    z_criterion::Real
    z_rejections::Union{ZScoreRejectionInfo,Nothing}
    rejected::Vector{Rejection}
    repaired::Union{OrderedDict{Int,Vector{Symbol}},Nothing}
    skipped::Union{OrderedDict{Int,Vector{Symbol}},Nothing}
end

"""Return unique rejections from an `EpochRejectionInfo`."""
unique_rejections(info::EpochRejectionInfo) = unique_rejections(info.rejected)
"""Return unique rejections for each `EpochRejectionInfo` in a vector."""
function unique_rejections(infos::Vector{EpochRejectionInfo})
    results = Vector{Rejection}[]
    for info in infos
        push!(results, unique_rejections(info))
    end
    return results
end

"""Return unique channels from an `EpochRejectionInfo` or vector thereof."""
unique_channels(rejections::EpochRejectionInfo) = unique_channels(rejections.rejected)
unique_channels(info::Vector{EpochRejectionInfo}) = unique_channels.(info)

"""Return unique epochs from an `EpochRejectionInfo` or vector thereof."""
unique_epochs(rejections::EpochRejectionInfo) = unique_epochs(rejections.rejected)
unique_epochs(info::Vector{EpochRejectionInfo}) = unique_epochs.(info)

"""Flatten all rejection data (abs + z-score) into a single DataFrame with reason tags."""
function all_data(info::EpochRejectionInfo)::DataFrame
    rows = NamedTuple{(:condition_number, :condition_name, :epoch, :channel, :reason),Tuple{Int,String,Int,Symbol,Symbol}}[]
    cnum = info.info.number
    cname = info.info.name
    _push!(rs, tag, list) =
        foreach(r -> push!(rs, (condition_number = cnum, condition_name = cname, epoch = r.epoch, channel = r.channel, reason = tag)), list)
    !isnothing(info.abs_rejections) && _push!(rows, :abs_threshold, info.abs_rejections)
    if !isnothing(info.z_rejections)
        z = info.z_rejections
        _push!(rows, :z_variance, z.z_variance)
        _push!(rows, :z_max, z.z_max)
        _push!(rows, :z_min, z.z_min)
        _push!(rows, :z_abs, z.z_abs)
        _push!(rows, :z_range, z.z_range)
        _push!(rows, :z_kurtosis, z.z_kurtosis)
    end
    return DataFrame(rows)
end
all_data(infos::Vector{EpochRejectionInfo})::DataFrame = vcat(all_data.(infos)...)



# === REPORTING FUNCTIONS ===

"""Pretty-print an `EpochRejectionInfo` summary with criterion, counts, and breakdowns."""
function Base.show(io::IO, info::EpochRejectionInfo)
    println(io, "EpochRejectionInfo: $(info.name)")
    println(io, "Condition: $(info.info.number): $(info.info.name)")
    println(io, "  Abs criterion: $(info.abs_criterion > 0 ? string(info.abs_criterion,  " μV") : "disabled")")
    println(io, "  Z-criterion: $(info.z_criterion > 0 ? string(info.z_criterion) : "disabled")")
    println(io, "  Epochs total: $(info.info.n), Epochs rejected: $(length(unique_epochs(info.rejected)))")
    println(io, "  Artifacts total: $(info.n_artifacts)")
    println(io, "  Rejected epochs: $(_print_vector(unique_epochs(info.rejected)))")

    if !isnothing(info.abs_rejections)
        println(io, "  Rejection breakdown (absolute):")
        println(
            io,
            "    Abs threshold: $(length(unique_epochs(info.abs_rejections))) unique epochs, $(length(unique_channels(info.abs_rejections))) unique channels",
        )
    end

    if !isnothing(info.z_rejections)
        z_info = info.z_rejections
        # Map measures to fields and labels
        field_map = Dict(
            :variance => (z_info.z_variance, "Z-Variance"),
            :max => (z_info.z_max, "Z-Maximum"),
            :min => (z_info.z_min, "Z-Minimum"),
            :abs => (z_info.z_abs, "Z-Absolute"),
            :range => (z_info.z_range, "Z-Range"),
            :kurtosis => (z_info.z_kurtosis, "Z-Kurtosis"),
        )

        # Determine which selected measures actually have entries
        nonempty_selected = Tuple{Vector{Rejection},String}[]
        for m in z_info.z_measures
            vec, label = field_map[m]
            if !isempty(vec)
                push!(nonempty_selected, (vec, label))
            end
        end
        if !isempty(nonempty_selected)
            println(io, "  Rejection breakdown (z-score):")
            for (vec, label) in nonempty_selected
                println(io, "    $(label):  $(length(unique_epochs(vec))) unique epochs, $(length(unique_channels(vec))) unique channels")
            end
        end
    end
    println(io, "")

end

"""Show each `EpochRejectionInfo` in a vector."""
Base.show(io::IO, infos::Vector{EpochRejectionInfo}) = Base.show.(Ref(io), infos)
