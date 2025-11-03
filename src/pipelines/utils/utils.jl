# Helper functions used during preprocess to make logging easier/prettier
function _applied_filters(filter_cfg::FilterConfig; filter_sections::Vector{Symbol} = [:highpass, :lowpass])
    applied_filters = []
    for name in filter_sections
        section = getproperty(filter_cfg, name)
        if section.apply
            push!(
                applied_filters,
                "$(section.freq) Hz $(section.type), $(section.method), $(section.func), order $(section.order)",
            )
        end
    end
    return join(applied_filters, "; ")
end

function _eog_config_string(eog_cfg::EogConfig)

    format_channels(channels) = begin
        left = join(channels[1], ", ")
        right = join(channels[2], ", ")
        ref = channels[3][1]
        return "$left - $right → $ref"
    end

    vEOG = format_channels(eog_cfg.vEOG_channels)
    hEOG = format_channels(eog_cfg.hEOG_channels)
return "vEOG: $vEOG ($(eog_cfg.vEOG_criterion) μV); hEOG: $hEOG ($(eog_cfg.hEOG_criterion) μV)"
end

_flag_symbol(base::AbstractString, criterion) = Symbol("$(base)_$(criterion)")

function _center_title(title::String, width::Int)
    total_dashes = width - length(title) - 2  # 2 spaces around title
    left_dashes = div(total_dashes, 2)
    right_dashes = total_dashes - left_dashes
    return "-" ^ left_dashes * " $title " * "-" ^ right_dashes
end

function section(title::String; width::Int = 80)
    dash_line = "-" ^ width
    middle_line = _center_title(title, width)
    return "\n$dash_line\n$middle_line\n$dash_line"
end

subsection(title::String; width::Int = 80) = "\n" * _center_title(title, width)
subsubsection(title::String) = "\n# " * title

# =============================================================================
# CHANNEL REPAIR TRACKING
# =============================================================================

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
    repaired::OrderedDict{Int, Vector{Symbol}}
    method::Symbol
    neighbors::Union{Dict{Symbol, Vector{Symbol}}, Nothing}
    skipped::OrderedDict{Int, Vector{Symbol}}
end

"""
    ChannelRepairInfo

Comprehensive tracking of all channel repairs during preprocessing.

# Fields
- `continuous::Union{ContinuousRepairInfo, Nothing}`: Repairs at continuous level (Nothing if none)
- `epochs::Vector{EpochRepairInfo}`: Repairs at epoch level, one per condition
"""
struct ChannelRepairInfo
    continuous::Union{ContinuousRepairInfo, Nothing}
    epochs::Vector{EpochRepairInfo}
end

# Convenience constructors
ChannelRepairInfo() = ChannelRepairInfo(nothing, EpochRepairInfo[])

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
    total_repairs = sum(length(v) for v in values(info.repaired))
    total_skipped = sum(length(v) for v in values(info.skipped))
    println(io, "  Epochs with repairs: $(length(info.repaired))")
    println(io, "  Total channel repairs: $total_repairs")
    println(io, "  Total channels skipped: $total_skipped")
    if info.neighbors !== nothing && !isempty(info.neighbors)
        println(io, "  Neighbors used:")
        for (ch, neighs) in info.neighbors
            println(io, "    $ch → $(neighs)")
        end
    end
end

function Base.show(io::IO, info::ChannelRepairInfo)
    println(io, "ChannelRepairInfo:")
    if info.continuous !== nothing
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
"""
struct ArtifactInfo
    continuous_repairs::Vector{ContinuousRepairInfo}
    epoch_rejections::Vector{EpochRejectionInfo}
end

ArtifactInfo() = ArtifactInfo(ContinuousRepairInfo[], EpochRejectionInfo[])

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
end

# =============================================================================
# HELPER FUNCTIONS FOR CREATING REPAIR INFO
# =============================================================================

"""
    create_continuous_repair_info(method::Symbol; name::String="continuous_repair")

Create a new ContinuousRepairInfo tracking struct for continuous data repairs.
Initializes with empty repaired and skipped vectors.

# Arguments
- `method::Symbol`: Repair method to use
- `name::String`: Optional name/identifier (default: "continuous_repair")
"""
function create_continuous_repair_info(method::Symbol; name::String="continuous_repair")
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
This function only analyzes repairability - it does not perform any repairs.
Use `repair_channels!` to actually perform the repairs after this analysis.
"""
function channel_repairable!(
    repair_info::ContinuousRepairInfo,
    bad_channels::Vector{Symbol},
    layout::Layout,
)
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

@add_nonmutating channel_repairable!

# =============================================================================
# CONTINUOUS DATA REPAIR FUNCTIONS
# =============================================================================

"""
    repair_channels_neighbor!(data::ContinuousData, repair_info::ContinuousRepairInfo)

Repair channels using weighted neighbor interpolation.
Uses `repair_info.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Updates `repair_info` during repair to track actual repairs vs skips.

# Arguments
- `data::ContinuousData`: Continuous EEG data to repair (modified in-place)
- `repair_info::ContinuousRepairInfo`: Repair tracking struct with `repaired` channels already populated

# Returns
- `Nothing`: All modifications are in-place

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_channels_neighbor!(
    data::ContinuousData,
    repair_info::ContinuousRepairInfo,
)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return nothing
    end

    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using neighbor interpolation"

    repair_channels!(
        data,
        repair_info.repaired;
        method = :neighbor_interpolation,
        repair_info = repair_info,
    )

    return nothing
end

"""
    repair_channels_spherical!(data::ContinuousData, repair_info::ContinuousRepairInfo; m::Int=4, lambda::Float64=1e-5)

Repair channels using spherical spline interpolation.
Uses `repair_info.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Updates `repair_info` during repair to track actual repairs vs skips.

# Arguments
- `data::ContinuousData`: Continuous EEG data to repair (modified in-place)
- `repair_info::ContinuousRepairInfo`: Repair tracking struct with `repaired` channels already populated
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)

# Returns
- `Nothing`: All modifications are in-place

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_channels_spherical!(
    data::ContinuousData,
    repair_info::ContinuousRepairInfo;
    m::Int = 4,
    lambda::Float64 = 1e-5,
)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return nothing
    end

    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using spherical spline interpolation"

    repair_channels!(
        data,
        repair_info.repaired;
        method = :spherical_spline,
        repair_info = repair_info,
        m = m,
        lambda = lambda,
    )

    return nothing
end

"""
    repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Repair channels using the information in `repair_info`.
Orchestrator function similar to `repair_artifacts!` for epoch data.

# Arguments
- `data::ContinuousData`: Continuous EEG data to repair (modified in-place)
- `repair_info::ContinuousRepairInfo`: Repair tracking struct with `repaired` channels already populated
- `method::Symbol`: Repair method to use (default: :neighbor_interpolation)

# Returns
- `Nothing`: All modifications are in-place

# Notes
- `repair_info.repaired` should be populated by `channel_repairable!` before calling this function
"""
function repair_channels!(
    data::ContinuousData,
    repair_info::ContinuousRepairInfo;
    method::Symbol = :neighbor_interpolation,
    kwargs...,
)
    if method == :neighbor_interpolation
        repair_channels_neighbor!(data, repair_info; kwargs...)
    elseif method == :spherical_spline
        repair_channels_spherical!(data, repair_info; kwargs...)
    else
        throw(ArgumentError("Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"))
    end
    return nothing
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
    repaired::Union{OrderedDict{Int, Vector{Symbol}}, Dict{Int, Vector{Symbol}}},
    skipped::Union{OrderedDict{Int, Vector{Symbol}}, Dict{Int, Vector{Symbol}}},
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
    if method == :neighbor_interpolation && neighbors_dict !== nothing
        neighbors_info = Dict{Symbol, Vector{Symbol}}()
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
    
    return EpochRepairInfo(
        dat.condition,
        dat.condition_name,
        repairs_ordered,
        method,
        neighbors_info,
        skipped_ordered,
    )
end
