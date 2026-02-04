# Helper functions used during preprocess to make logging easier/prettier
function _applied_filters(filter_cfg::FilterConfig; filter_sections::Vector{Symbol} = [:highpass, :lowpass])
    applied_filters = []
    for name in filter_sections
        section = getproperty(filter_cfg, name)
        if section.apply
            push!(applied_filters, "$(section.freq) Hz $(section.type), $(section.method), $(section.func), order $(section.order)")
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
- `ica_components::Union{ArtifactComponents, Nothing}`: ICA artifact components (Nothing if ICA was not applied)
"""
struct ArtifactInfo
    continuous_repairs::Vector{ContinuousRepairInfo}
    epoch_rejections::Vector{EpochRejectionInfo}
    ica_components::Union{ArtifactComponents,Nothing}
end

ArtifactInfo() = ArtifactInfo(ContinuousRepairInfo[], EpochRejectionInfo[], nothing)

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
This function only analyzes repairability - it does not perform any repairs.
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
function repair_channels_neighbor!(data::ContinuousData, repair_info::ContinuousRepairInfo)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return nothing
    end

    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using neighbor interpolation"

    repair_channels!(data, repair_info.repaired; method = :neighbor_interpolation, repair_info = repair_info)

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
function repair_channels_spherical!(data::ContinuousData, repair_info::ContinuousRepairInfo; m::Int = 4, lambda::Float64 = 1e-5)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return nothing
    end

    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using spherical spline interpolation"

    repair_channels!(data, repair_info.repaired; method = :spherical_spline, repair_info = repair_info, m = m, lambda = lambda)

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
function repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol = :neighbor_interpolation, kwargs...)
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
    if method == :neighbor_interpolation && neighbors_dict !== nothing
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

"""
    _combine_and_average_epoch_counts(all_epoch_counts::Vector{DataFrame})

Combine epoch counts from all files and create averaged version.

# Arguments
- `all_epoch_counts::Vector{DataFrame}`: Vector of DataFrames with epoch counts per file

# Returns
- `Tuple{DataFrame, DataFrame}`: (combined_counts, avg_counts)
  - `combined_counts`: All epoch counts concatenated
  - `avg_counts`: Averaged percentage per file (grouped by file)
"""
function _epoch_and_file_summary(all_epoch_counts::Vector{DataFrame})
    epoch_summary = vcat(all_epoch_counts...)
    file_summary = combine(groupby(epoch_summary, [:file]), :percentage => mean => :percentage)
    file_summary.percentage = round.(file_summary.percentage; digits = 1)
    return epoch_summary, file_summary
end

"""
    _merge_summaries(new_epoch_summary::DataFrame, new_file_summary::DataFrame, output_directory::String)

Merge new summaries with existing summaries, replacing data for files that already exist.

# Arguments
- `new_epoch_summary::DataFrame`: New epoch summary data
- `new_file_summary::DataFrame`: New file summary data
- `output_directory::String`: Directory where summary files are stored

# Returns
- `Tuple{DataFrame, DataFrame}`: (merged_epoch_summary, merged_file_summary)
"""
function _merge_summaries(new_epoch_summary::DataFrame, new_file_summary::DataFrame, output_directory::String)
    epoch_path = joinpath(output_directory, "epoch_summary.jld2")
    file_path = joinpath(output_directory, "file_summary.jld2")

    # Get list of files in new data (for replacement)
    new_files = unique(new_epoch_summary.file)

    # Load and merge epoch summary
    if isfile(epoch_path)
        existing_epoch = JLD2.load(epoch_path, "data")
        # Remove rows for files that are being updated
        existing_epoch = existing_epoch[.!in.(existing_epoch.file, Ref(new_files)), :]
        merged_epoch = vcat(existing_epoch, new_epoch_summary)
    else
        merged_epoch = new_epoch_summary
    end

    # Load and merge file summary
    if isfile(filepath)
        existing_file = JLD2.load(file_path, "data")
        # Remove rows for files that are being updated
        existing_file = existing_file[.!in.(existing_file.file, Ref(new_files)), :]
        merged_file = vcat(existing_file, new_file_summary)
    else
        merged_file = new_file_summary
    end

    # Sort by filename using natural sorting (file_1, file_2, ..., file_10)
    merged_epoch = sort(merged_epoch, :file, by = natural_sort_key)
    merged_file = sort(merged_file, :file, by = natural_sort_key)

    return merged_epoch, merged_file
end

"""
    compare_rejections(rejection_step1::Vector{EpochRejectionInfo}, rejection_step2::Vector{EpochRejectionInfo})

Compare rejection results before and after repair to show repair effectiveness.

# Arguments
- `rejection_step1::Vector{EpochRejectionInfo}`: Rejection info before repair
- `rejection_step2::Vector{EpochRejectionInfo}`: Rejection info after repair

# Returns
- `DataFrame`: Comparison table with detailed statistics per condition showing:
  - Original number of epochs
  - Rejected before repair (step1)
  - Rejected after repair (step2)
  - Successfully repaired epochs
  - Percentage of step1 rejections that were repaired
  - Final percentage of epochs kept

# Examples
```julia
rejection_step1 = detect_bad_epochs_automatic(epochs, name="step1")
repair_artifacts!(epochs, rejection_step1)
rejection_step2 = detect_bad_epochs_automatic(epochs, name="step2")
comparison = compare_rejections(rejection_step1, rejection_step2)
```
"""
function compare_rejections(rejection1::Vector{EpochRejectionInfo}, rejection2::Vector{EpochRejectionInfo})
    comparison_data = []

    for (info1, info2) in zip(rejection1, rejection2)
        # Get unique epoch indices rejected in each step
        epochs_rejected_step1 = unique_epochs(info1)
        epochs_rejected_step2 = unique_epochs(info2)

        # Calculate statistics
        n_original = info1.info.n
        n_rejected_step1 = length(epochs_rejected_step1)
        n_rejected_step2 = length(epochs_rejected_step2)

        # Epochs that were rejected in step1 but NOT in step2 (successfully repaired)
        epochs_repaired = setdiff(epochs_rejected_step1, epochs_rejected_step2)
        n_repaired = length(epochs_repaired)

        # Epochs rejected in BOTH steps (repair didn't help)
        epochs_still_bad = intersect(epochs_rejected_step1, epochs_rejected_step2)
        n_still_bad = length(epochs_still_bad)

        # New rejections in step2 (epochs that weren't rejected in step1)
        # This can happen if repair introduces artifacts or if thresholds are slightly different
        epochs_new_rejections = setdiff(epochs_rejected_step2, epochs_rejected_step1)
        n_new_rejections = length(epochs_new_rejections)

        # Verify: n_rejected_step2 should equal n_still_bad + n_new_rejections
        # Verify: n_rejected_step1 should equal n_repaired + n_still_bad
        @assert n_rejected_step2 == n_still_bad + n_new_rejections "Logic error: n_rejected_step2 ($n_rejected_step2) != n_still_bad ($n_still_bad) + n_new_rejections ($n_new_rejections)"
        @assert n_rejected_step1 == n_repaired + n_still_bad "Logic error: n_rejected_step1 ($n_rejected_step1) != n_repaired ($n_repaired) + n_still_bad ($n_still_bad)"

        # Calculate key percentages
        pct_repaired = n_rejected_step1 > 0 ? round(100 * n_repaired / n_rejected_step1, digits = 1) : 0.0
        pct_final_kept = round(100 * (n_original - n_rejected_step2) / n_original, digits = 1)

        push!(
            comparison_data,
            (
                condition = info1.info.number,
                condition_name = info1.info.name,
                n_original = n_original,
                rejected_before_repair = n_rejected_step1,
                rejected_after_repair = n_rejected_step2,
                successfully_repaired = n_repaired,
                pct_repaired = pct_repaired,
                pct_final_kept = pct_final_kept,
            ),
        )
    end

    return DataFrame(comparison_data)
end

"""
    summarize_electrode_repairs(continuous_repairs::Vector{ContinuousRepairInfo})

Summarize electrode repairs from a vector of ContinuousRepairInfo objects (continuous level only).
Each ContinuousRepairInfo is treated as coming from a separate participant.

# Arguments
- `continuous_repairs::Vector{ContinuousRepairInfo}`: Vector of ContinuousRepairInfo objects

# Returns
- `DataFrame`: Summary table with columns:
  - `electrode::Symbol`: Electrode name
  - `n_participants::Int`: Number of participants where this electrode was repaired at the continuous level
"""
function summarize_electrode_repairs(continuous_repairs::Vector{ContinuousRepairInfo})
    electrode_participants = Dict{Symbol,Set{Int}}()

    for (participant_idx, repair_info) in enumerate(continuous_repairs)
        # Count continuous repairs (1 per participant)
        for ch in repair_info.repaired
            if !haskey(electrode_participants, ch)
                electrode_participants[ch] = Set{Int}()
            end
            push!(electrode_participants[ch], participant_idx)
        end
    end

    # Create DataFrame sorted by number of participants (descending)
    summary_data = [
        (electrode = ch, n_participants = length(participants)) for
        (ch, participants) in sort(collect(electrode_participants), by = x -> length(x[2]), rev = true)
    ]

    return DataFrame(summary_data)
end

"""
    summarize_electrode_repairs(file_pattern::String; input_dir::String = pwd())

Summarize electrode repairs by loading ArtifactInfo files from a directory (continuous level only).
Uses `_find_batch_files` to match files by pattern.

# Arguments
- `file_pattern::String`: Pattern to match filenames (e.g., "_artifact_info")
- `input_dir::String`: Directory containing `_artifact_info.jld2` files (default: current directory)

# Returns
- `DataFrame`: Summary table with columns:
  - `electrode::Symbol`: Electrode name
  - `n_participants::Int`: Number of participants where this electrode was repaired at the continuous level

# Examples
```julia
# Load all artifact info files from current directory
summary = summarize_electrode_repairs("_artifact_info")

# Load from specific directory
summary = summarize_electrode_repairs("_artifact_info", input_dir="/path/to/output")
```
"""
function summarize_electrode_repairs(file_pattern::String; input_dir::String = pwd())
    # Find matching files using _find_batch_files
    files = _find_batch_files(file_pattern, input_dir)

    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return DataFrame(electrode = Symbol[], n_participants = Int[])
    end

    # Load all ArtifactInfo objects and extract ContinuousRepairInfo
    continuous_repairs = ContinuousRepairInfo[]
    for file in files
        file_path = joinpath(input_dir, file)
        try
            if isfile(filepath)
                artifact_info = read_data(filepath)
                # If read_data returned a Dict, extract "data" key
                if isa(artifact_info, Dict) && haskey(artifact_info, "data")
                    artifact_info = artifact_info["data"]
                end
                for repair_info in artifact_info.continuous_repairs
                    if !isempty(repair_info.repaired)
                        push!(continuous_repairs, repair_info)
                    end
                end
            end
        catch e
            @minimal_warning "Failed to load artifact info from $file_path: $e"
        end
    end

    return summarize_electrode_repairs(continuous_repairs)
end





"""
    summarize_ica_components(ica_components::Vector{ArtifactComponents})

Summarize ICA components removed across all participants.

# Arguments
- `ica_components::Vector{ArtifactComponents}`: Vector of ArtifactComponents objects

# Returns
- `DataFrame`: Summary table with columns:
  - `file::String`: Filename (empty string when called with Vector{ArtifactComponents})
  - `total_components::Int`: Total number of ICA components removed
  - `vEOG::Int`: Number of vEOG components
  - `hEOG::Int`: Number of hEOG components
  - `ECG::Int`: Number of ECG components
  - `line_noise::Int`: Number of line noise components
  - `channel_noise::Int`: Number of channel noise components
- `DataFrame`: Average summary with columns:
  - `component_type::String`: Type of component
  - `avg_per_participant::Float64`: Average number of components per participant
"""
function summarize_ica_components(ica_components::Vector{ArtifactComponents})
    if isempty(ica_components)
        return DataFrame(
            file = String[],
            total_components = Int[],
            vEOG = Int[],
            hEOG = Int[],
            ECG = Int[],
            line_noise = Int[],
            channel_noise = Int[],
        ),
        DataFrame(component_type = String[], avg_per_participant = Float64[])
    end

    # Per-file summary
    per_file_data = []
    for component_artifacts in ica_components
        vEOG_count = length(get(component_artifacts.eog, :vEOG, []))
        hEOG_count = length(get(component_artifacts.eog, :hEOG, []))
        ecg_count = length(component_artifacts.ecg)
        line_noise_count = length(component_artifacts.line_noise)
        channel_noise_count = length(component_artifacts.channel_noise)
        total_count = vEOG_count + hEOG_count + ecg_count + line_noise_count + channel_noise_count

        push!(
            per_file_data,
            (
                file = "",  # No filename when called with Vector{ArtifactComponents}
                total_components = total_count,
                vEOG = vEOG_count,
                hEOG = hEOG_count,
                ECG = ecg_count,
                line_noise = line_noise_count,
                channel_noise = channel_noise_count,
            ),
        )
    end

    per_file_df = DataFrame(per_file_data)

    # Average summary
    n_participants = length(ica_components)
    avg_data = [
        (component_type = "vEOG", avg_per_participant = round(mean(per_file_df.vEOG), digits = 2)),
        (component_type = "hEOG", avg_per_participant = round(mean(per_file_df.hEOG), digits = 2)),
        (component_type = "ECG", avg_per_participant = round(mean(per_file_df.ECG), digits = 2)),
        (component_type = "line_noise", avg_per_participant = round(mean(per_file_df.line_noise), digits = 2)),
        (component_type = "channel_noise", avg_per_participant = round(mean(per_file_df.channel_noise), digits = 2)),
        (component_type = "total", avg_per_participant = round(mean(per_file_df.total_components), digits = 2)),
    ]

    avg_df = DataFrame(avg_data)

    return per_file_df, avg_df
end

"""
    summarize_ica_components(file_pattern::String; input_dir::String = pwd())

Summarize ICA components by loading ArtifactInfo files from a directory.
Uses `_find_batch_files` to match files by pattern.

# Arguments
- `file_pattern::String`: Pattern to match filenames (e.g., "_artifact_info")
- `input_dir::String`: Directory containing `_artifact_info.jld2` files (default: current directory)

# Returns
- `DataFrame`: Per-file summary table with columns:
  - `file::String`: Base filename
  - `total_components::Int`: Total number of ICA components removed
  - `vEOG::Int`: Number of vEOG components
  - `hEOG::Int`: Number of hEOG components
  - `ECG::Int`: Number of ECG components
  - `line_noise::Int`: Number of line noise components
  - `channel_noise::Int`: Number of channel noise components
- `DataFrame`: Average summary with columns:
  - `component_type::String`: Type of component
  - `avg_per_participant::Float64`: Average number of components per participant

# Examples
```julia
# Load all artifact info files from current directory
per_file, avg = summarize_ica_components("_artifact_info")

# Load from specific directory
per_file, avg = summarize_ica_components("_artifact_info", input_dir="/path/to/output")
```
"""
function summarize_ica_components(file_pattern::String; input_dir::String = pwd())
    # Find matching files using _find_batch_files
    files = _find_batch_files(file_pattern, input_dir)

    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return DataFrame(
            file = String[],
            total_components = Int[],
            vEOG = Int[],
            hEOG = Int[],
            ECG = Int[],
            line_noise = Int[],
            channel_noise = Int[],
        ),
        DataFrame(component_type = String[], avg_per_participant = Float64[])
    end

    # Load all ArtifactInfo objects and extract ArtifactComponents with filenames
    ica_components = ArtifactComponents[]
    filenames = String[]

    for file in files
        file_path = joinpath(input_dir, file)
        try
            if isfile(filepath)
                artifact_info = read_data(filepath)
                # If read_data returned a Dict, extract "data" key
                if isa(artifact_info, Dict) && haskey(artifact_info, "data")
                    artifact_info = artifact_info["data"]
                end
                if !isnothing(artifact_info.ica_components)
                    push!(ica_components, artifact_info.ica_components)
                    push!(filenames, file)  # _find_batch_files returns filenames, not full paths
                end
            end
        catch e
            @minimal_warning "Failed to load artifact info from $file_path: $e"
        end
    end

    # Get summary from vector version
    per_file_df, avg_df = summarize_ica_components(ica_components)

    # Update filenames in per_file_df: extract base filename
    for (i, filename) in enumerate(filenames)
        # Extract base filename (remove "_artifact_info.jld2" suffix)
        base_filename = replace(filename, "_artifact_info.jld2" => "")
        per_file_df.file[i] = base_filename
    end

    # Sort by filename using natural sorting (handles numeric parts correctly: Flank_C_3 before Flank_C_12)
    sort!(per_file_df, :file, by = natural_sort_key)

    return per_file_df, avg_df
end
