# Tracking and reporting structures for pipeline channel repairs.

######################################
######### Types & Constants ##########
######################################

"""
    ContinuousRepairInfo

Tracks channels repaired and skipped at the continuous data level (`name`, `method`, `repaired`, `skipped`).

# Fields
- `name::String`: Identifier for the continuous repair stage
- `method::Symbol`: Spatial interpolation method used (`:spherical_spline` or `:neighbor_interpolation`)
- `repaired::Vector{Symbol}`: Channels successfully interpolated
- `skipped::Vector{Symbol}`: Bad channels skipped (e.g., insufficient clean neighbors)
"""
mutable struct ContinuousRepairInfo
    name::String
    method::Symbol
    repaired::Vector{Symbol}
    skipped::Vector{Symbol}
end

"""
    EpochRepairInfo

Tracks channels repaired and skipped per epoch for a condition (`condition`, `condition_name`, `repaired`, `method`, `neighbors`, `skipped`).

# Fields
- `condition::Int64`: Numeric condition code
- `condition_name::String`: Condition label
- `repaired::OrderedDict{Int,Vector{Symbol}}`: Repaired channel symbols keyed by epoch index
- `method::Symbol`: Spatial interpolation method used
- `neighbors::Union{Dict{Symbol,Vector{Symbol}},Nothing}`: Neighbor channels used per repaired electrode, if recorded
- `skipped::OrderedDict{Int,Vector{Symbol}}`: Skipped channel symbols keyed by epoch index
"""
struct EpochRepairInfo
    condition::Int64
    condition_name::String
    repaired::OrderedDict{Int,Vector{Symbol}}
    method::Symbol
    neighbors::Union{Dict{Symbol,Vector{Symbol}},Nothing}
    skipped::OrderedDict{Int,Vector{Symbol}}
end

"""
    ChannelRepairInfo

Combined container for `continuous` (`Union{ContinuousRepairInfo,Nothing}`) and `epochs` (`Vector{EpochRepairInfo}`) channel repairs.

# Fields
- `continuous::Union{ContinuousRepairInfo,Nothing}`: Continuous-stage channel repair record, if any
- `epochs::Vector{EpochRepairInfo}`: Per-condition epoch-stage channel repair records
"""
struct ChannelRepairInfo
    continuous::Union{ContinuousRepairInfo,Nothing}
    epochs::Vector{EpochRepairInfo}
end

"""
    ArtifactInfo <: EegFunData

Container for preprocessing artifact records (`continuous_repairs`, `epoch_rejections`, `ica_components`).

# Fields
- `continuous_repairs::Vector{ContinuousRepairInfo}`: Continuous-data channel repair logs
- `epoch_rejections::Vector{EpochRejectionInfo}`: Per-condition epoch rejection and repair logs
- `ica_components::Union{ArtifactComponents,Nothing}`: Identified/removed ICA artifact components, if any
"""
struct ArtifactInfo <: EegFunData
    continuous_repairs::Vector{ContinuousRepairInfo}
    epoch_rejections::Vector{EpochRejectionInfo}
    ica_components::Union{ArtifactComponents,Nothing}
end

######################################
########### User Functions ###########
######################################

ChannelRepairInfo() = ChannelRepairInfo(nothing, EpochRepairInfo[])

ArtifactInfo() = ArtifactInfo(ContinuousRepairInfo[], EpochRejectionInfo[], nothing)

"""
    create_continuous_repair_info(method::Symbol; name::String = "continuous_repair") -> ContinuousRepairInfo

Create an empty `ContinuousRepairInfo` initialized for `method`.
"""
create_continuous_repair_info(method::Symbol; name::String = "continuous_repair") = ContinuousRepairInfo(name, method, Symbol[], Symbol[])

"""
    channel_repairable(repair_info::ContinuousRepairInfo, bad_channels::Vector{Symbol}, layout::Layout) -> ContinuousRepairInfo
    channel_repairable(artifacts::EpochRejectionInfo, layout::Layout) -> EpochRejectionInfo

Partition bad channels into `repaired` and `skipped` based on valid neighbor availability in `layout`,
returning an updated copy without modifying the input.
"""
function channel_repairable(dat, args...; kwargs...)
    dat_copy = copy(dat)
    channel_repairable!(dat_copy, args...; kwargs...)
    return dat_copy
end

"""
    channel_repairable!(repair_info::ContinuousRepairInfo, bad_channels::Vector{Symbol}, layout::Layout) -> ContinuousRepairInfo

In-place mutating version of [`channel_repairable`](@ref).
"""
function channel_repairable!(repair_info::ContinuousRepairInfo, bad_channels::Vector{Symbol}, layout::Layout)
    if isempty(bad_channels)
        repair_info.repaired = Symbol[]
        repair_info.skipped = Symbol[]
        return repair_info
    end

    repairable_channels = check_channel_neighbors(bad_channels, layout)
    repair_info.repaired = repairable_channels
    repair_info.skipped = setdiff(bad_channels, repairable_channels)

    !isempty(repairable_channels) && @info "Can repair $(length(repairable_channels)) channel(s): $repairable_channels"
    !isempty(repair_info.skipped) &&
        @info "Cannot repair $(length(repair_info.skipped)) channel(s) (bad neighbors and/or fewer than 2 neighbors): $(repair_info.skipped)"
    isempty(repairable_channels) &&
        @info "No channels can be repaired (all $(length(bad_channels)) bad channels have bad neighbors and/or fewer than 2 neighbors)"

    return repair_info
end

"""
    create_epoch_repair_info(dat::EpochData, repaired, skipped, method::Symbol, neighbors_dict = nothing) -> EpochRepairInfo

Construct an `EpochRepairInfo` sorted by epoch index from per-epoch `repaired` and `skipped` channel maps.
"""
function create_epoch_repair_info(
    dat::EpochData,
    repaired::Union{OrderedDict{Int,Vector{Symbol}},Dict{Int,Vector{Symbol}}},
    skipped::Union{OrderedDict{Int,Vector{Symbol}},Dict{Int,Vector{Symbol}}},
    method::Symbol,
    neighbors_dict::Union{AbstractDict,Nothing} = nothing,
)
    repairs_ordered = repaired isa OrderedDict ? repaired : OrderedDict(sort!(collect(repaired); by = first))
    skipped_ordered = skipped isa OrderedDict ? skipped : OrderedDict(sort!(collect(skipped); by = first))

    neighbors_info = nothing
    if method == :neighbor_interpolation && !isnothing(neighbors_dict)
        neighbors_info = Dict{Symbol,Vector{Symbol}}()
        for channels in values(repairs_ordered), ch in channels
            if !haskey(neighbors_info, ch) && haskey(neighbors_dict, ch)
                neighbors_info[ch] = neighbors_dict[ch].channels
            end
        end
    end

    return EpochRepairInfo(dat.condition, dat.condition_name, repairs_ordered, method, neighbors_info, skipped_ordered)
end

Base.copy(info::ContinuousRepairInfo) = ContinuousRepairInfo(info.name, info.method, copy(info.repaired), copy(info.skipped))

function Base.show(io::IO, info::ContinuousRepairInfo)
    println(io, "ContinuousRepairInfo: $(info.name)")
    println(io, "  Method: $(info.method)")
    println(io, "  Channels repaired: $(length(info.repaired)) - $(info.repaired)")
    println(io, "  Channels skipped: $(length(info.skipped)) - $(info.skipped)")
end

function Base.show(io::IO, info::EpochRepairInfo)
    println(io, "EpochRepairInfo:")
    println(io, "  Condition: $(info.condition) ($(info.condition_name))")
    println(io, "  Method: $(info.method)")
    total_repairs = sum(length, values(info.repaired); init = 0)
    total_skipped = sum(length, values(info.skipped); init = 0)
    println(io, "  Epochs with repairs: $(length(info.repaired))")
    println(io, "  Total channel repairs: $total_repairs")
    println(io, "  Total channels skipped: $total_skipped")
    if !isnothing(info.neighbors) && !isempty(info.neighbors)
        println(io, "  Neighbors used:")
        for (ch, neighs) in info.neighbors
            println(io, "    $ch → $(neighs)")
        end
    end
end

function Base.show(io::IO, info::ChannelRepairInfo)
    println(io, "ChannelRepairInfo:")
    if !isnothing(info.continuous)
        println(io, "Continuous level:")
        println(io, info.continuous)
    else
        println(io, "Continuous level: No repairs")
    end
    println(io, "Epoch level:")
    if isempty(info.epochs)
        println(io, "  No repairs")
    else
        for info_item in info.epochs
            println(io, "  ", info_item)
        end
    end
end

function Base.show(io::IO, info::ArtifactInfo)
    println(io, "ArtifactInfo:")
    println(io, "  Continuous repairs: $(length(info.continuous_repairs))")
    for repair_info in info.continuous_repairs
        println(io, "    - $(repair_info.name) ($(repair_info.method))")
    end
    println(io, "  Epoch rejections: $(length(info.epoch_rejections))")
    for rejection_info in info.epoch_rejections
        println(io, "    - $(rejection_info.name) (condition $(rejection_info.info.number): $(rejection_info.info.name))")
    end
    if !isnothing(info.ica_components)
        println(io, "  ICA components removed: $(length(get_all_ica_components(info.ica_components)))")
    else
        println(io, "  ICA components: Not applied")
    end
end
