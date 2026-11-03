"""
    Rejection(channel::Symbol, epoch::Int)

Rejected channel-epoch pair from artifact detection.
"""
struct Rejection
    channel::Symbol
    epoch::Int
end

Base.show(io::IO, r::Rejection) = print(io, "Rejection(:$(r.channel), $(r.epoch))")

is_equal_rejection(a::Rejection, b::Rejection) = a.channel == b.channel && a.epoch == b.epoch

"""
    unique_rejections(rejections)

De-duplicate `Rejection` entries by `(channel, epoch)`, or extract unique rejections from `EpochRejectionInfo`.
"""
function unique_rejections(rejections::Vector{Rejection})
    seen = Set{Tuple{Symbol,Int}}()
    out = Rejection[]
    for r in rejections
        key = (r.channel, r.epoch)
        if key ∉ seen
            push!(seen, key)
            push!(out, r)
        end
    end
    return out
end

unique_channels(rejections::Vector{Rejection}) = unique(map(r -> r.channel, rejections))
unique_epochs(rejections::Vector{Rejection}) = unique(map(r -> r.epoch, rejections))

struct EpochInfo
    number::Int
    name::String
    n::Int
end

"""
    ZScoreRejectionInfo

Per-measure z-score rejection lists (`:variance`, `:max`, `:min`, `:abs`, `:range`, `:kurtosis`).
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

Base.copy(z::ZScoreRejectionInfo) = ZScoreRejectionInfo(
    copy(z.z_measures),
    copy(z.z_variance),
    copy(z.z_max),
    copy(z.z_min),
    copy(z.z_abs),
    copy(z.z_range),
    copy(z.z_kurtosis),
)

"""
    EpochRejectionInfo

Stores epoch rejection criteria, per-channel/epoch rejections, and optional channel repair tracking (`repaired`, `skipped`).
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

function Base.copy(info::EpochRejectionInfo)
    return EpochRejectionInfo(
        info.name,
        info.info,
        info.n_artifacts,
        info.abs_criterion,
        isnothing(info.abs_rejections) ? nothing : copy(info.abs_rejections),
        info.z_criterion,
        isnothing(info.z_rejections) ? nothing : copy(info.z_rejections),
        copy(info.rejected),
        isnothing(info.repaired) ? nothing : OrderedDict{Int,Vector{Symbol}}(k => copy(v) for (k, v) in info.repaired),
        isnothing(info.skipped) ? nothing : OrderedDict{Int,Vector{Symbol}}(k => copy(v) for (k, v) in info.skipped),
    )
end

unique_rejections(info::EpochRejectionInfo) = unique_rejections(info.rejected)
unique_rejections(infos::Vector{EpochRejectionInfo}) = unique_rejections.(infos)

unique_channels(info::EpochRejectionInfo) = unique_channels(info.rejected)
unique_channels(infos::Vector{EpochRejectionInfo}) = unique_channels.(infos)

unique_epochs(info::EpochRejectionInfo) = unique_epochs(info.rejected)
unique_epochs(infos::Vector{EpochRejectionInfo}) = unique_epochs.(infos)

"""Flatten all rejection data (abs + z-score) into a single `DataFrame` with reason tags."""
function all_data(info::EpochRejectionInfo)::DataFrame
    rows = NamedTuple{(:condition_number, :condition_name, :epoch, :channel, :reason),Tuple{Int,String,Int,Symbol,Symbol}}[]
    cnum, cname = info.info.number, info.info.name
    _push!(tag, list) = foreach(
        r -> push!(rows, (condition_number = cnum, condition_name = cname, epoch = r.epoch, channel = r.channel, reason = tag)),
        list,
    )
    !isnothing(info.abs_rejections) && _push!(:abs_threshold, info.abs_rejections)
    if !isnothing(info.z_rejections)
        z = info.z_rejections
        for (tag, list) in (
            (:z_variance, z.z_variance),
            (:z_max, z.z_max),
            (:z_min, z.z_min),
            (:z_abs, z.z_abs),
            (:z_range, z.z_range),
            (:z_kurtosis, z.z_kurtosis),
        )
            _push!(tag, list)
        end
    end
    return DataFrame(rows)
end
all_data(infos::Vector{EpochRejectionInfo})::DataFrame = vcat(all_data.(infos)...)

function Base.show(io::IO, info::EpochRejectionInfo)
    println(io, "EpochRejectionInfo: $(info.name)")
    println(io, "Condition: $(info.info.number): $(info.info.name)")
    println(io, "  Abs criterion: $(info.abs_criterion > 0 ? string(info.abs_criterion, " μV") : "disabled")")
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
        z = info.z_rejections
        entries = (
            (:variance, z.z_variance, "Z-Variance"),
            (:max, z.z_max, "Z-Maximum"),
            (:min, z.z_min, "Z-Minimum"),
            (:abs, z.z_abs, "Z-Absolute"),
            (:range, z.z_range, "Z-Range"),
            (:kurtosis, z.z_kurtosis, "Z-Kurtosis"),
        )
        nonempty = [(vec, label) for (m, vec, label) in entries if m in z.z_measures && !isempty(vec)]
        if !isempty(nonempty)
            println(io, "  Rejection breakdown (z-score):")
            for (vec, label) in nonempty
                println(io, "    $(label):  $(length(unique_epochs(vec))) unique epochs, $(length(unique_channels(vec))) unique channels")
            end
        end
    end
    println(io, "")
end

Base.show(io::IO, infos::Vector{EpochRejectionInfo}) = Base.show.(Ref(io), infos)
