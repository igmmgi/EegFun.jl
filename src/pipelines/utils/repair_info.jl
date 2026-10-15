# === CHANNEL REPAIR TRACKING ===

"""
    ContinuousRepairInfo

Stores information about channels repaired at the continuous data level
(before epoching). These repairs apply to all epochs.

# Fields
- `name::String`: Name/identifier for this repair info (e.g., "continuous_repair_pre_ica")
- `method::Symbol`: Repair method used (:neighbor_interpolation or :spherical_spline)
- `repaired::Vector{Symbol}`: Channels that were successfully repaired (populated during repair)
- `skipped::Vector{Symbol}`: Channels that were identified as bad but couldn't be repaired (populated during repair)
"""
mutable struct ContinuousRepairInfo
    name::String
    method::Symbol
    repaired::Vector{Symbol}
    skipped::Vector{Symbol}
end

Base.copy(info::ContinuousRepairInfo) = ContinuousRepairInfo(info.name, info.method, copy(info.repaired), copy(info.skipped))

"""
    EpochRepairInfo

Stores information about channels repaired at the epoch level for a specific condition.

# Fields
- `condition::Int64`: Condition number
- `condition_name::String`: Condition name
- `repaired::OrderedDict{Int, Vector{Symbol}}`: Maps epoch index to channels repaired in that epoch (ordered by epoch number)
- `method::Symbol`: Repair method used (:neighbor_interpolation or :spherical_spline)
- `neighbors::Union{Dict{Symbol, Vector{Symbol}}, Nothing}`: For neighbor interpolation,
  maps each repaired channel to its neighbor channels (may be the same across epochs)
- `skipped::OrderedDict{Int, Vector{Symbol}}`: Maps epoch index to channels that were
  identified as bad but couldn't be repaired (bad neighbors or too few neighbors) (ordered by epoch number)
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

Comprehensive tracking of all channel repairs during preprocessing.

# Fields
- `continuous::Union{ContinuousRepairInfo, Nothing}`: Repairs at continuous level (Nothing if none)
- `epochs::Vector{EpochRepairInfo}`: Repairs at epoch level, one per condition
"""
struct ChannelRepairInfo
    continuous::Union{ContinuousRepairInfo,Nothing}
    epochs::Vector{EpochRepairInfo}
end

"""Create an empty `ChannelRepairInfo` with no repairs."""
ChannelRepairInfo() = ChannelRepairInfo(nothing, EpochRepairInfo[])

"""Pretty-print continuous-level channel repair info."""
function Base.show(io::IO, info::ContinuousRepairInfo)
    println(io, "ContinuousRepairInfo: $(info.name)")
    println(io, "  Method: $(info.method)")
    println(io, "  Channels repaired: $(length(info.repaired)) - $(info.repaired)")
    println(io, "  Channels skipped: $(length(info.skipped)) - $(info.skipped)")
end

"""Pretty-print epoch-level channel repair info."""
function Base.show(io::IO, info::EpochRepairInfo)
    println(io, "EpochRepairInfo:")
    println(io, "  Condition: $(info.condition) ($(info.condition_name))")
    println(io, "  Method: $(info.method)")
    total_repairs = sum(length(v) for v in values(info.repaired))
    total_skipped = sum(length(v) for v in values(info.skipped))
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

"""Pretty-print combined channel repair info (continuous + epoch)."""
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

"""
    ArtifactInfo

Flexible container for storing all artifact-related information from preprocessing.
Can hold multiple continuous repair infos and epoch rejection infos.

# Fields
- `continuous_repairs::Vector{ContinuousRepairInfo}`: Vector of continuous repair info objects
- `epoch_rejections::Vector{EpochRejectionInfo}`: Vector of epoch rejection info objects
- `ica_components::Union{ArtifactComponents, Nothing}`: ICA artifact components (Nothing if ICA was not applied)
"""
struct ArtifactInfo <: EegFunData
    continuous_repairs::Vector{ContinuousRepairInfo}
    epoch_rejections::Vector{EpochRejectionInfo}
    ica_components::Union{ArtifactComponents,Nothing}
end

"""Create an empty `ArtifactInfo` with no repairs, rejections, or ICA."""
ArtifactInfo() = ArtifactInfo(ContinuousRepairInfo[], EpochRejectionInfo[], nothing)

"""Pretty-print artifact info summary (repairs, rejections, ICA components)."""
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
        all_comps = get_all_ica_components(info.ica_components)
        println(io, "  ICA components removed: $(length(all_comps))")
    else
        println(io, "  ICA components: Not applied")
    end
end

# === HELPER FUNCTIONS FOR CREATING REPAIR INFO ===

"""
    create_continuous_repair_info(method::Symbol; name::String="continuous_repair")

Create a new ContinuousRepairInfo tracking struct for continuous data repairs.
Initializes with empty repaired and skipped vectors.

# Arguments
- `method::Symbol`: Repair method to use
- `name::String`: Optional name/identifier (default: "continuous_repair")
"""
function create_continuous_repair_info(method::Symbol; name::String = "continuous_repair")
    return ContinuousRepairInfo(name, method, Symbol[], Symbol[])
end

"""
    channel_repairable!(repair_info::ContinuousRepairInfo, bad_channels::Vector{Symbol}, layout::Layout)

Analyze which channels can be repaired and which cannot, based on neighbor availability.
Populates `repair_info.repaired` and `repair_info.skipped` with the analysis.

# Arguments
- `repair_info::ContinuousRepairInfo`: Continuous repair tracking struct (mutated to add repair analysis)
- `bad_channels::Vector{Symbol}`: Channels identified as bad that need repair consideration
- `layout::Layout`: Layout object containing neighbor information

# Returns
- `ContinuousRepairInfo`: The same repair_info object (modified in-place)

# Notes
Only analyzes repairability - does not perform any repairs.
Use `repair_channels!` to actually perform the repairs after this analysis.
"""
function channel_repairable!(repair_info::ContinuousRepairInfo, bad_channels::Vector{Symbol}, layout::Layout)
    if isempty(bad_channels)
        repair_info.repaired = Symbol[]
        repair_info.skipped = Symbol[]
        return repair_info
    end

    # Determine which channels can be repaired
    repairable_channels = check_channel_neighbors(bad_channels, layout)

    # Categorize channels
    repair_info.repaired = repairable_channels
    repair_info.skipped = setdiff(bad_channels, repairable_channels)

    # Log the results
    if !isempty(repairable_channels)
        @info "Can repair $(length(repairable_channels)) channel(s): $repairable_channels"
    end

    if !isempty(repair_info.skipped)
        @info "Cannot repair $(length(repair_info.skipped)) channel(s) (bad neighbors and/or fewer than 2 neighbors): $(repair_info.skipped)"
    end

    if isempty(repairable_channels) && !isempty(bad_channels)
        @info "No channels can be repaired (all $(length(bad_channels)) bad channels have bad neighbors and/or fewer than 2 neighbors)"
    end

    return repair_info
end

"""
    channel_repairable(dat, args...; kwargs...)

Non-mutating version of `channel_repairable!`.
"""
function channel_repairable(dat, args...; kwargs...)
    dat_copy = copy(dat)
    channel_repairable!(dat_copy, args...; kwargs...)
    return dat_copy
end

"""
    create_epoch_repair_info(dat::EpochData,
                              repaired::Union{OrderedDict{Int, Vector{Symbol}}, Dict{Int, Vector{Symbol}}},
                              skipped::Union{OrderedDict{Int, Vector{Symbol}}, Dict{Int, Vector{Symbol}}},
                              method::Symbol,
                              neighbors_dict::Union{OrderedDict, Nothing} = nothing)

Create EpochRepairInfo from actual repair tracking data.
The repaired and skipped should be populated during the repair process.
If not already OrderedDict, will convert to OrderedDict sorted by epoch number.
"""
function create_epoch_repair_info(
    dat::EpochData,
    repaired::Union{OrderedDict{Int,Vector{Symbol}},Dict{Int,Vector{Symbol}}},
    skipped::Union{OrderedDict{Int,Vector{Symbol}},Dict{Int,Vector{Symbol}}},
    method::Symbol,
    neighbors_dict::Union{OrderedDict,Nothing} = nothing,
)
    # Convert to OrderedDict if not already, sorted by epoch number
    repairs_ordered = if repaired isa OrderedDict
        repaired
    else
        OrderedDict(sort(pairs(repaired)))
    end

    skipped_ordered = if skipped isa OrderedDict
        skipped
    else
        OrderedDict(sort(pairs(skipped)))
    end

    neighbors_info = nothing
    if method == :neighbor_interpolation && !isnothing(neighbors_dict)
        neighbors_info = Dict{Symbol,Vector{Symbol}}()
        # Collect unique repaired channels and their neighbors
        all_repaired_channels = Set{Symbol}()
        for channels in values(repairs_ordered)
            union!(all_repaired_channels, channels)
        end
        for ch in all_repaired_channels
            if haskey(neighbors_dict, ch)
                neighbors_info[ch] = neighbors_dict[ch].channels
            end
        end
    end

    return EpochRepairInfo(dat.condition, dat.condition_name, repairs_ordered, method, neighbors_info, skipped_ordered)
end
