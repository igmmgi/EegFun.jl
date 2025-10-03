# =============================================================================
# ARTIFACT DETECTION
# =============================================================================

import .==
using Makie

"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)

Detect EOG (eye movement) onsets in the data and mark them as extreme values.

# Arguments
- `dat::ContinuousData`: The continuous EEG data object
- `criterion::Real`: Threshold criterion for EOG detection
- `channel_in::Symbol`: Input channel symbol for EOG detection
- `channel_out::Symbol`: Output channel symbol to mark detected EOG events

# Modifies
- `dat`: Adds extreme value markers to the specified output channel

# Examples
```julia
# Detect EOG onsets using HEOG channel
detect_eog_onsets!(dat, 50.0, :HEOG, :is_eog_onset)
```
"""
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)
    if !(channel_in in propertynames(dat.data))
        @minimal_error_throw("Channel $channel_in not found in data")
    end
    dat.data[!, channel_out] = _is_extreme_value(dat.data[!, channel_in], Float64(criterion))
    return nothing
end

"""
    _is_extreme_value(signal::AbstractVector{Float64}, threshold::Float64)

Detect extreme values in a signal using threshold crossing.

# Arguments
- `signal::AbstractVector{Float64}`: Input signal vector
- `threshold::Float64`: Threshold for extreme value detection

# Returns
- `Vector{Bool}`: Boolean vector indicating extreme values

# Examples
```julia
# Detect values above 50 μV
extreme_mask = _is_extreme_value(signal, 50.0)
```
"""
function _is_extreme_value(signal::AbstractVector{Float64}, threshold::Float64)
    return abs.(signal) .> threshold
end

"""
    _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Float64)

Detect extreme values in a signal and store results in a pre-allocated mask.

# Arguments
- `mask::Vector{Bool}`: Pre-allocated boolean vector to store results
- `signal::AbstractVector{Float64}`: Input signal vector
- `threshold::Float64`: Threshold for extreme value detection

# Modifies
- `mask`: Filled with extreme value detection results

# Examples
```julia
# Pre-allocate mask and detect extreme values
mask = Vector{Bool}(undef, length(signal))
_is_extreme_value!(mask, signal, 50.0)
```
"""
function _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Float64)
    @assert length(mask) == length(signal) "Mask and signal must have the same length"
    @inbounds for i in eachindex(signal)
        mask[i] = abs(signal[i]) > threshold
    end
    return nothing
end

"""
    is_extreme_value!(dat::SingleDataFrameEeg, threshold::Int; 
                     channel_selection::Function = channels(), 
                     sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                     mode::Symbol = :combined,
                     channel_out::Union{Symbol, Nothing} = nothing)

Detect extreme values across selected channels and add results to the data.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
- `mode::Symbol`: Mode of operation - `:combined` (single combined column, default) or `:separate` (separate columns per channel)
- `channel_out::Union{Symbol, Nothing}`: Output channel name for combined mode (default: auto-generated as `is_extreme_value_<threshold>`)

# Modifies
- `dat`: Adds extreme value detection columns to the data

# Examples
```julia
# Detect extreme values and combine into single output channel (default, auto-named)
is_extreme_value!(dat, 100)  # Creates column :is_extreme_value_100

# Detect extreme values with custom output channel name
is_extreme_value!(dat, 100, channel_out = :is_extreme_value)

# Detect extreme values and add separate columns for each channel
is_extreme_value!(dat, 100, mode = :separate)

# Detect extreme values for specific channels (auto-named)
is_extreme_value!(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Detect extreme values only for selected samples (auto-named)
is_extreme_value!(dat, 100, sample_selection = sample_mask)
```
"""
function is_extreme_value!(
    dat::SingleDataFrameEeg,
    threshold::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)

    # Validate mode
    if mode ∉ [:separate, :combined]
        @minimal_error_throw("mode must be :separate or :combined, got :$mode")
    end

    if mode == :combined
        # Combined mode - use all channels (same as original behavior)
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

        # Use provided channel_out or generate default name
        output_channel = channel_out === nothing ? Symbol("is_extreme_value_$(threshold)") : channel_out
        dat.data[!, output_channel] = falses(nrow(dat.data))

        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            dat.data[!, output_channel] .|= extreme_mask
        end

    else  # Separate mode - use specified channel selection
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)
        for (ch, extreme_mask) in results
            column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
            dat.data[!, column_name] = extreme_mask
        end
    end

    return nothing
end

# Helper function to detect extreme values for selected channels
function _detect_extreme_values(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for extreme value detection")
    end

    selected_samples = get_selected_samples(dat, sample_selection)

    results = Dict{Symbol,Vector{Bool}}()

    for ch in selected_channels
        channel_data = dat.data[!, ch]
        extreme_mask = _is_extreme_value(channel_data, Float64(threshold))
        
        # Apply sample selection - only keep extreme values for selected samples
        sample_mask = falses(length(extreme_mask))
        sample_mask[selected_samples] .= true
        extreme_mask = extreme_mask .& sample_mask
        
        results[ch] = extreme_mask
    end

    return results
end

"""
    is_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                    channel_selection::Function = channels(), 
                    sample_selection::Function = samples(),
                    mode::Symbol = :combined)

Detect extreme values across selected channels and return results.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
- `mode::Symbol`: Mode of operation - `:combined` (boolean vector, default) or `:separate` (DataFrame with separate columns per channel)

# Returns
- `DataFrame` or `Vector{Bool}`: Results depending on mode

# Examples
```julia
# Detect extreme values and return combined boolean vector (default)
extreme_mask = is_extreme_value(dat, 100)

# Detect extreme values and return DataFrame with separate columns
results = is_extreme_value(dat, 100, mode = :separate)

# Detect extreme values for specific channels
extreme_mask = is_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Detect extreme values only for selected samples
extreme_mask = is_extreme_value(dat, 100, sample_selection = sample_mask)
```
"""
function is_extreme_value(
    dat::SingleDataFrameEeg,
    threshold::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    mode::Symbol = :combined,
)

    if mode == :combined
        # Combined mode - return boolean vector directly
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

        # Initialize the result with false values
        combined_mask = Vector{Bool}(falses(nrow(dat.data)))

        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            combined_mask .|= extreme_mask
        end

        return combined_mask

    else  # mode == :separate
        # Separate mode - create temporary data object and use mutating version
        temp_dat = deepcopy(dat)
        is_extreme_value!(temp_dat, threshold; channel_selection, sample_selection, mode = :separate)

        # Extract the extreme value columns in the same order as the original channels
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        extreme_cols = [Symbol("is_extreme_value_$(ch)_$(threshold)") for ch in selected_channels]
        return temp_dat.data[:, extreme_cols]
    end
end

"""
    n_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                   channel_selection::Function = channels(), 
                   sample_selection::Function = samples(),
                   mode::Symbol = :combined)

Count the number of extreme values across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
- `mode::Symbol`: Mode for extreme value detection (:separate or :combined, default: :combined)

# Returns
- `DataFrame`: DataFrame with extreme value counts for each channel (separate mode) or total count (combined mode)

# Examples
```julia
# Count total extreme values across all channels (combined mode, default)
total_count = n_extreme_value(dat, 100)

# Count extreme values in specific channels (combined mode)
total_count = n_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Count extreme values only for selected samples (combined mode)
total_count = n_extreme_value(dat, 100, sample_selection = sample_mask)

# Count extreme values in all channels (separate mode)
count_df = n_extreme_value(dat, 100, mode = :separate)
```
"""
function n_extreme_value(
    dat::SingleDataFrameEeg,
    threshold::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    mode::Symbol = :combined,
)
    # Validate mode
    if mode ∉ [:separate, :combined]
        @minimal_error_throw("mode must be :separate or :combined, got :$mode")
    end

    if mode == :combined
        # Combined mode - use all channels (same as original behavior)
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

        # Initialize the result with false values
        combined_mask = Vector{Bool}(falses(nrow(dat.data)))

        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            combined_mask .|= extreme_mask
        end

        # Return total count
        return sum(combined_mask)

    else  # Separate mode - count for each channel
        
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

        # Get channels in original order
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

        # Count extreme values for each channel in original order
        counts = [sum(results[ch]) for ch in selected_channels]

        # Create result DataFrame
        result_df = DataFrame(
            channel = selected_channels,
            n_extreme = counts,
            threshold = fill(threshold, length(selected_channels)),
        )

        return result_df
    end
end

"""
    _n_extreme_value(df::DataFrame, channels::Vector{Symbol}, threshold::Float64)

Internal function to count extreme values for specified channels.

# Arguments
- `df::DataFrame`: DataFrame containing the data
- `channels::Vector{Symbol}`: Vector of channel symbols to analyze
- `threshold::Float64`: Threshold for extreme value detection

# Returns
- `Vector{Int}`: Count of extreme values for each channel
"""
function _n_extreme_value(df::DataFrame, channels::Vector{Symbol}, threshold::Float64)
    counts = Int[]

    for ch in channels
        channel_data = df[!, ch]
        extreme_mask = _is_extreme_value(channel_data, threshold)
        push!(counts, sum(extreme_mask))
    end

    return counts
end


"""
Automatic epoch rejection based on statistical criteria.

This module provides functions to automatically reject epochs based on
statistical measures (variance, max, min, absolute max, range, kurtosis)
using z-score thresholds. This is useful for removing epochs with artifacts
without manual inspection.
"""

#=============================================================================
    REJECTION INFORMATION STRUCTURE
=============================================================================#

struct Rejection
    label::Symbol
    epoch::Int
end

==(a::eegfun.Rejection, b::eegfun.Rejection) = a.label == b.label && a.epoch == b.epoch



"""
    EpochRejectionInfo

Stores information about which epochs were rejected and why.

# Fields
- `n_epochs::Int`: Number of epochs before rejection
- `n_artifacts::Int`: Number of epochs after rejection
- `rejected_epochs::Vector{Int}`: Indices of rejected epochs
- `rejected_by_z_variance::Vector{Int}`: Epochs rejected due to high variance (z-score)
- `rejected_by_z_max::Vector{Int}`: Epochs rejected due to high maximum values (z-score)
- `rejected_by_z_min::Vector{Int}`: Epochs rejected due to low minimum values (z-score)
- `rejected_by_z_abs::Vector{Int}`: Epochs rejected due to high absolute values (z-score)
- `rejected_by_z_range::Vector{Int}`: Epochs rejected due to large range (z-score)
- `rejected_by_z_kurtosis::Vector{Int}`: Epochs rejected due to high kurtosis (z-score)
- `rejected_by_abs_threshold::Vector{Int}`: Epochs rejected due to absolute voltage threshold
- `z_criterion::Float64`: Z-score criterion used for rejection
- `abs_criterion::Union{Float64, Nothing}`: Absolute voltage threshold (μV) used for rejection
"""
struct EpochRejectionInfo
    n_epochs::Int
    n_artifacts::Int
    rejected_epochs::Vector{Rejection}
    z_variance::Vector{Rejection}
    z_max::Vector{Rejection}
    z_min::Vector{Rejection}
    z_abs::Vector{Rejection}
    z_range::Vector{Rejection}
    z_kurtosis::Vector{Rejection}
    absolute_threshold::Vector{Rejection}
    z_criterion::Float64
    abs_criterion::Union{Float64, Nothing}
    channel_metrics::Dict{Symbol, Dict{Symbol, Float64}}
    channel_rejections::Dict{Symbol, Int}
end

#=============================================================================
    CORE REJECTION FUNCTIONS
=============================================================================#

"""
    reject_epochs_automatic!(dat::EpochData, z_criterion::Real;
                             channel_selection::Function = channels())::EpochRejectionInfo

Automatically reject epochs based on statistical criteria using z-score threshold.

This function calculates six statistical measures for each epoch (variance, max, min,
absolute max, range, kurtosis) across selected channels. For each measure, the maximum
across channels is taken for each epoch, then z-scored across epochs. Epochs exceeding
the z-criterion for any measure are rejected.

# Arguments
- `dat::EpochData`: Epoched EEG data to process
- `z_criterion::Real`: Z-score threshold for rejection (typically 2.0 or 3.0)
- `channel_selection::Function`: Channel predicate for selecting channels to analyze (default: all channels)

# Returns
- `EpochRejectionInfo`: Information about which epochs were rejected and why

# Effects
- Modifies the input data in-place by removing rejected epochs

# Examples
```julia
using eegfun, JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Reject epochs with z-score > 2.0
rejection_info = reject_epochs_automatic!(epochs, 2.0)

# Check results
println("Original epochs: \$(rejection_info.n_epochs)")
println("Remaining epochs: \$(rejection_info.n_artifacts)")
println("Rejected epochs: \$(rejection_info.rejected_epochs)")

# Save cleaned data
save("participant_1_epochs_cleaned.jld2", "epochs", epochs)
```

# Notes
- Common z-criteria: 2.0 (more aggressive), 2.5, 3.0 (more conservative)
- Rejection is based on ANY metric exceeding the criterion
- Uses maximum across channels to identify global artifacts
- All six metrics are calculated independently and combined with OR logic
"""
function detect_bad_epochs(
    dat::EpochData,
    z_criterion::Real;
    abs_criterion::Union{Real, Nothing} = nothing,
    channel_selection::Function = channels(),
)::EpochRejectionInfo
    # Validate inputs
    z_criterion <= 0 && @minimal_error_throw("Z-criterion must be positive, got $z_criterion")
    if abs_criterion !== nothing && abs_criterion <= 0
        @minimal_error_throw("Absolute criterion must be positive, got $abs_criterion")
    end

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    println("Selected channels: $(print_vector(selected_channels))")
    isempty(selected_channels) && @minimal_error_throw("No channels selected for epoch rejection")

    # Calculate metrics and identify rejected epochs
    metrics = _calculate_epoch_metrics(dat, selected_channels, Float64(z_criterion), abs_criterion !== nothing ? Float64(abs_criterion) : nothing)

    # Combine all rejected epochs
    rejected_epochs = sort(unique(vcat(
        values(metrics[:z_variance])...,
        values(metrics[:z_max])...,
        values(metrics[:z_min])...,
        values(metrics[:z_abs])...,
        values(metrics[:z_range])...,
        values(metrics[:z_kurtosis])...,
        values(metrics[:absolute_threshold])...
    )))

    # Create Rejection structs for each rejection
    rejected_epochs_info = Rejection[]
    z_variance = Rejection[]
    z_max = Rejection[]
    z_min = Rejection[]
    z_abs = Rejection[]
    z_range = Rejection[]
    z_kurtosis = Rejection[]
    absolute_threshold = Rejection[]

    # Convert metric vectors to Sets for O(1) lookups
    metric_sets = Dict{Symbol, Dict{Symbol, Set{Int}}}()
    for metric_key in [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis, :absolute_threshold]
        metric_sets[metric_key] = Dict{Symbol, Set{Int}}()
        for ch in selected_channels
            metric_sets[metric_key][ch] = Set(metrics[metric_key][ch])
        end
    end

    # Define metric mappings
    metric_mappings = [
        (:z_variance, z_variance, metric_sets[:z_variance]),
        (:z_max, z_max, metric_sets[:z_max]),
        (:z_min, z_min, metric_sets[:z_min]),
        (:z_abs, z_abs, metric_sets[:z_abs]),
        (:z_range, z_range, metric_sets[:z_range]),
        (:z_kurtosis, z_kurtosis, metric_sets[:z_kurtosis])
    ]
    
    # Add absolute threshold if provided
    if abs_criterion !== nothing
        push!(metric_mappings, (:absolute_threshold, absolute_threshold, metric_sets[:absolute_threshold]))
    end

    for epoch_idx in rejected_epochs
        for ch in selected_channels
            rejection = nothing
            for (_, rejection_list, metric_set) in metric_mappings
                if epoch_idx in metric_set[ch]
                    rejection = Rejection(ch, epoch_idx) 
                    if rejection !== nothing
                        push!(rejected_epochs_info, rejection)
                        push!(rejection_list, rejection)
                    end
                end
            end
        end
    end

    # Create rejection info
    rejection_info = EpochRejectionInfo(
        length(dat.data),  # n_epochs
        length(rejected_epochs),
        unique(rejected_epochs_info),
        z_variance,
        z_max,
        z_min,
        z_abs,
        z_range,
        z_kurtosis,
        absolute_threshold,
        Float64(z_criterion),
        abs_criterion !== nothing ? Float64(abs_criterion) : nothing,
        Dict{Symbol, Dict{Symbol, Float64}}(),  # channel_metrics (optional)
        Dict{Symbol, Int}()  # channel_rejections (optional)
    )

    return rejection_info
end

"""
    get_rejected_epochs(state::EpochRejectionState)::Vector{Int}

Get indices of rejected epochs from the rejection state.

# Examples
```julia
state = reject_epochs_interactive(epochs)
# ... after review ...
rejected_indices = get_rejected_epochs(state)
```
"""
function get_rejected_epochs(info::EpochRejectionInfo)::Vector{Int}
    return info.rejected_epochs
end

#=============================================================================
    INTERNAL HELPER FUNCTIONS
=============================================================================#

"""
Validate inputs for epoch rejection.
"""
function _validate_rejection_inputs(dat::EpochData, z_criterion::Real)
    if isempty(dat.data)
        @minimal_error_throw("Cannot reject epochs from empty EpochData")
    end
    if z_criterion <= 0
        @minimal_error_throw("Z-criterion must be positive, got $z_criterion")
    end
    if length(dat.data) < 3
        @minimal_warning "Only $(length(dat.data)) epochs available. Statistical rejection may not be meaningful with so few epochs."
    end
end


"""
Calculate statistical metrics for all epochs.

Returns a Dict with keys: :variance, :max, :min, :abs, :range, :kurtosis
Each value is a vector of length n_epochs containing the maximum of that metric across channels.
"""
function _calculate_epoch_metrics(
    dat::EpochData,
    selected_channels::Vector{Symbol},
    z_criterion::Float64,
    abs_criterion::Union{Float64, Nothing}
)::Dict{Symbol, Dict{Symbol, Vector{Int}}}
    # Initialize metrics dictionary
    metrics = Dict{Symbol, Dict{Symbol, Vector{Int}}}()

    # Preallocate dictionaries for each metric
    metrics[:z_variance] = Dict{Symbol, Vector{Int}}()
    metrics[:z_max] = Dict{Symbol, Vector{Int}}()
    metrics[:z_min] = Dict{Symbol, Vector{Int}}()
    metrics[:z_abs] = Dict{Symbol, Vector{Int}}()
    metrics[:z_range] = Dict{Symbol, Vector{Int}}()
    metrics[:z_kurtosis] = Dict{Symbol, Vector{Int}}()
    metrics[:absolute_threshold] = Dict{Symbol, Vector{Int}}()

    # Calculate metrics for each channel
    for ch in selected_channels
        # Initialize vectors for this channel
        for metric_key in [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis, :absolute_threshold]
            metrics[metric_key][ch] = Int[]
        end

        # Extract all channel data at once
        channel_data_all = [epoch[!, ch] for epoch in dat.data]
        
        # Calculate all metrics vectorized
        variances = var.(channel_data_all)
        max_values = maximum.(channel_data_all)
        min_values = minimum.(channel_data_all)
        abs_values = maximum.([abs.(epoch_data) for epoch_data in channel_data_all])
        ranges = max_values .- min_values
        kurtoses = kurtosis.(channel_data_all)

        # Calculate z-scores and identify exceeding epochs
        z_scores = zscore.([variances, max_values, min_values, abs_values, ranges, kurtoses])
        metric_keys = [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis]

        for (z_score, metric_key) in zip(z_scores, metric_keys)
            exceeding_epochs = findall(abs.(z_score) .> z_criterion)
            append!(metrics[metric_key][ch], exceeding_epochs)
        end

        # Check absolute threshold
        if abs_criterion !== nothing
            abs_threshold_violations = findall(epoch -> any(abs.(epoch[!, ch]) .> abs_criterion), dat.data)
            append!(metrics[:absolute_threshold][ch], abs_threshold_violations)
        end
    end

    return metrics
end





#=============================================================================
    REPORTING FUNCTIONS
=============================================================================#

"""
    Base.show(io::IO, info::EpochRejectionInfo)

Display rejection information in a human-readable format.
"""
function Base.show(io::IO, info::EpochRejectionInfo)
    println(io, "EpochRejectionInfo:")
    println(io, "  Z-criterion: $(info.z_criterion)")
    if info.abs_criterion !== nothing
        println(io, "  Abs criterion: $(info.abs_criterion) μV")
    end
    println(io, "  Number epochs: $(info.n_epochs)")
    println(io, "  Number artifacts: $(info.n_artifacts)")
    println(io, "  Rejected epochs:")
    # for rejection in info.rejected_epochs
    #     println(io, "    $(rejection.label): epoch $(rejection.epoch)")
    # end
    println(io, "")
    println(io, "  Rejection breakdown (z-score):")
    println(io, "    Z-Variance:  $(length(info.z_variance)) epochs")
    println(io, "    Z-Maximum:   $(length(info.z_max)) epochs")
    println(io, "    Z-Minimum:   $(length(info.z_min)) epochs")
    println(io, "    Z-Absolute:  $(length(info.z_abs)) epochs")
    println(io, "    Z-Range:     $(length(info.z_range)) epochs")
    println(io, "    Z-Kurtosis:  $(length(info.z_kurtosis)) epochs")
    if info.abs_criterion !== nothing
        println(io, "")
        println(io, "  Rejection breakdown (absolute):")
        println(io, "    Abs threshold: $(length(info.absolute_threshold)) epochs")
    end
end

#=============================================================================
    UNIFIED ARTIFACT DETECTION & REPAIR API
=============================================================================#


"""
    ArtifactRepairMethod

Enumeration of available artifact repair methods.

# Values
- `:neighbor_interpolation` - Weighted neighbor interpolation (default)
- `:spherical_spline` - Spherical spline interpolation
- `:reject` - Remove bad epochs entirely
"""
@enum ArtifactRepairMethod neighbor_interpolation spherical_spline reject



"""
    repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo, method::Symbol=:neighbor_interpolation; kwargs...)

Repair detected artifacts using the specified method.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use

# Keyword Arguments (for :neighbor_interpolation method)
- `neighbours_dict::Union{OrderedDict, Nothing}`: Neighbor information (default: auto-generate)

# Keyword Arguments (for :spherical_spline method)
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)
"""
function repair_artifacts!(
    dat::EpochData, 
    artifacts::EpochRejectionInfo, 
    method::Symbol=:neighbor_interpolation;
    kwargs...
)
    # Convert symbol to enum
    repair_method = if method == :neighbor_interpolation
        neighbor_interpolation
    elseif method == :spherical_spline
        spherical_spline
    elseif method == :reject
        reject
    else
        throw(ArgumentError("Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline, :reject"))
    end
    
    return repair_artifacts!(dat, artifacts, repair_method; kwargs...)
end

function repair_artifacts!(
    dat::EpochData, 
    artifacts::EpochRejectionInfo, 
    method::ArtifactRepairMethod;
    kwargs...
)
    if method == neighbor_interpolation
        return repair_artifacts_neighbor!(dat, artifacts; kwargs...)
    elseif method == spherical_spline
        return repair_artifacts_spherical_spline!(dat, artifacts; kwargs...)
    elseif method == reject
        return repair_artifacts_reject!(dat, artifacts; kwargs...)
    else
        throw(ArgumentError("Unknown repair method: $method"))
    end
end

"""
    repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo; neighbours_dict::Union{OrderedDict, Nothing}=nothing)

Repair artifacts using weighted neighbor interpolation.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `neighbours_dict::Union{OrderedDict, Nothing}`: Neighbor information (default: auto-generate)

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)
"""
function repair_artifacts_neighbor!(
    dat::EpochData, 
    artifacts::EpochRejectionInfo;
    neighbours_dict::Union{OrderedDict, Nothing}=nothing
)
    # Get neighbor information
    if isnothing(neighbours_dict)
        neighbours_dict = get_electrode_neighbours_xyz(dat.layout)
    end
    
    # Get all rejected epochs
    rejected_epochs = unique([r.epoch for r in artifacts.rejected_epochs])
    
    for epoch_idx in rejected_epochs
        # Get bad channels for this epoch
        bad_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx]
        isempty(bad_channels) && continue
        
        # Get current epoch data
        epoch = dat.data[epoch_idx]
        
        # Repair each bad channel
        for bad_ch in bad_channels
            # Get neighbor information
            neighbours = get(neighbours_dict, bad_ch, nothing)
            if isnothing(neighbours) || isempty(neighbours.electrodes)
                @warn "No neighbors found for channel $bad_ch, skipping"
                continue
            end
            
            # Find good neighbors (not in bad_channels)
            good_neighbors = setdiff(neighbours.electrodes, bad_channels)
            if length(good_neighbors) < 2
                @warn "Not enough good neighbors for channel $bad_ch, skipping"
                continue
            end
            
            # Calculate weights based on distance
            weights = Float64[]
            for neighbor in good_neighbors
                neighbor_idx = findfirst(==(neighbor), neighbours.electrodes)
                if !isnothing(neighbor_idx)
                    weight = 1.0 / neighbours.distances[neighbor_idx]
                    push!(weights, weight)
                end
            end
            
            # Normalize weights
            weights ./= sum(weights)
            
            # Interpolate using weighted average
            interpolated = zeros(nrow(epoch))
            for (weight, neighbor) in zip(weights, good_neighbors)
                if hasproperty(epoch, neighbor)
                    interpolated .+= weight .* epoch[!, neighbor]
                end
            end
            
            # Update the bad channel data
            epoch[!, bad_ch] = interpolated
        end
    end
    
    return dat
end

"""
    repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int=4, lambda::Float64=1e-5)

Repair artifacts using spherical spline interpolation.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)
"""
function repair_artifacts_spherical_spline!(
    dat::EpochData, 
    artifacts::EpochRejectionInfo;
    m::Int=4, 
    lambda::Float64=1e-5
)
    # Get all rejected epochs
    rejected_epochs = unique([r.epoch for r in artifacts.rejected_epochs])
    
    for epoch_idx in rejected_epochs
        # Get bad channels for this epoch
        bad_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx]
        isempty(bad_channels) && continue
        
        # Get current epoch data
        epoch = dat.data[epoch_idx]
        
        # Use the spherical spline repair function from channel_repair.jl
        # We need to convert the epoch data to the format expected by the repair function
        channels = Symbol.(dat.layout.data.label)
        data_matrix = Matrix(epoch[:, channels])'
        
        # Call the spherical spline repair
        repair_channels_spherical_spline!(data_matrix, bad_channels, channels, dat.layout; m=m, lambda=lambda)
        
        # Update the epoch data
        epoch[:, channels] = data_matrix'
    end
    
    return dat
end

"""
    repair_artifacts_reject!(dat::EpochData, artifacts::EpochRejectionInfo)

Remove epochs that contain artifacts entirely.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)
"""
function repair_artifacts_reject!(
    dat::EpochData, 
    artifacts::EpochRejectionInfo
)
    # Get all rejected epochs
    rejected_epochs = unique([r.epoch for r in artifacts.rejected_epochs])
    
    # Remove rejected epochs (in reverse order to maintain indices)
    for epoch_idx in sort(rejected_epochs, rev=true)
        deleteat!(dat.data, epoch_idx)
    end
    
    return dat
end

#=============================================================================
    VISUALIZATION FUNCTIONS
=============================================================================#

"""
    plot_artifact_detection(epochs::EpochData, artifacts::EpochRejectionInfo; epoch_idx::Int=1, channels::Vector{Symbol}=Symbol[])

Plot artifact detection results showing original data with rejected channels/epochs highlighted.

# Arguments
- `epochs::EpochData`: The epoch data
- `artifacts::EpochRejectionInfo`: Artifact detection results
- `epoch_idx::Int`: Epoch index to plot (default: 1)
- `channels::Vector{Symbol}`: Specific channels to plot (default: all channels)

# Returns
- `Figure`: Makie figure showing the plot
"""
function plot_artifact_detection(
    epochs::EpochData, 
    artifacts::EpochRejectionInfo; 
    epoch_idx::Int=1, 
    channels::Vector{Symbol}=Symbol[]
)
    
    # Get channels to plot
    if isempty(channels)
        channels = Symbol.(epochs.layout.data.label)
    end
    
    # Get epoch data
    epoch = epochs.data[epoch_idx]
    time_points = epoch.time
    
    # Find rejected channels for this epoch
    rejected_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx]
    
    # Create figure
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Amplitude (μV)", title="Artifact Detection - Epoch $epoch_idx")
    
    # Plot each channel
    for (i, ch) in enumerate(channels)
        if hasproperty(epoch, ch)
            color = ch in rejected_channels ? :red : :blue
            alpha = ch in rejected_channels ? 0.8 : 0.6
            linewidth = ch in rejected_channels ? 2 : 1
            
            lines!(ax, time_points, epoch[!, ch], 
                   color=color, alpha=alpha, linewidth=linewidth, 
                   label=ch in rejected_channels ? "$ch (rejected)" : "$ch")
        end
    end
    
    # Add legend
    axislegend(ax, position=:topright)
    
    return fig
end

"""
    plot_repair_comparison(epochs_original::EpochData, epochs_repaired::EpochData, artifacts::EpochRejectionInfo; 
                          epoch_idx::Int=1, channels::Vector{Symbol}=Symbol[])

Plot comparison between original and repaired epochs.

# Arguments
- `epochs_original::EpochData`: Original epoch data
- `epochs_repaired::EpochData`: Repaired epoch data  
- `artifacts::EpochRejectionInfo`: Artifact detection results
- `epoch_idx::Int`: Epoch index to plot (default: 1)
- `channels::Vector{Symbol}`: Specific channels to plot (default: all channels)

# Returns
- `Figure`: Makie figure showing before/after comparison
"""
function plot_repair_comparison(
    epochs_original::EpochData, 
    epochs_repaired::EpochData, 
    artifacts::EpochRejectionInfo;
    epoch_idx::Int=1, 
    channels::Vector{Symbol}=Symbol[]
)
    
    # Get channels to plot
    if isempty(channels)
        channels = Symbol.(epochs_original.layout.data.label)
    end
    
    # Get epoch data
    epoch_orig = epochs_original.data[epoch_idx]
    epoch_repaired = epochs_repaired.data[epoch_idx]
    time_points = epoch_orig.time
    
    # Find rejected channels for this epoch
    rejected_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx]
    
    # Create figure with subplots
    fig = Figure(size=(1000, 800))
    
    # Original data
    ax1 = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Amplitude (μV)", title="Original Data - Epoch $epoch_idx")
    
    # Repaired data  
    ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel="Amplitude (μV)", title="Repaired Data - Epoch $epoch_idx")
    
    # Plot each channel
    for (i, ch) in enumerate(channels)
        if hasproperty(epoch_orig, ch) && hasproperty(epoch_repaired, ch)
            # Original plot
            color = ch in rejected_channels ? :red : :blue
            alpha = ch in rejected_channels ? 0.8 : 0.6
            linewidth = ch in rejected_channels ? 2 : 1
            
            lines!(ax1, time_points, epoch_orig[!, ch], 
                   color=color, alpha=alpha, linewidth=linewidth,
                   label=ch in rejected_channels ? "$ch (rejected)" : "$ch")
            
            # Repaired plot
            lines!(ax2, time_points, epoch_repaired[!, ch], 
                   color=:green, alpha=0.8, linewidth=1,
                   label="$ch (repaired)")
        end
    end
    
    # Add legends
    axislegend(ax1, position=:topright)
    axislegend(ax2, position=:topright)
    
    # Link axes
    linkaxes!(ax1, ax2)
    
    return fig
end

# """
#     plot_erp_comparison(epochs_original::EpochData, epochs_repaired::EpochData; 
#                        channels::Vector{Symbol}=Symbol[], condition::Symbol=:all)
# 
# Plot ERP comparison between original and repaired data.
# 
# # Arguments
# - `epochs_original::EpochData`: Original epoch data
# - `epochs_repaired::EpochData`: Repaired epoch data
# - `channels::Vector{Symbol}`: Specific channels to plot (default: all channels)
# - `condition::Symbol`: Condition to plot (default: :all)
# 
# # Returns
# - `Figure`: Makie figure showing ERP comparison
# """
# function plot_erp_comparison(
#     epochs_original::EpochData, 
#     epochs_repaired::EpochData;
#     channels::Vector{Symbol}=Symbol[],
#     condition::Symbol=:all
# )
#     using Makie
#     
#     # Get channels to plot
#     if isempty(channels)
#         channels = Symbol.(epochs_original.layout.data.label)
#     end
#     
#     # Calculate ERPs
#     erp_orig = grandaverage(epochs_original, condition=condition)
#     erp_repaired = grandaverage(epochs_repaired, condition=condition)
#     
#     time_points = erp_orig.time
#     
#     # Create figure
#     fig = Figure(size=(1200, 800))
#     
#     # Plot each channel
#     n_channels = length(channels)
#     n_cols = min(4, n_channels)
#     n_rows = ceil(Int, n_channels / n_cols)
#     
#     for (i, ch) in enumerate(channels)
#         if hasproperty(erp_orig, ch) && hasproperty(erp_repaired, ch)
#             row = ((i-1) ÷ n_cols) + 1
#             col = ((i-1) % n_cols) + 1
#             
#             ax = Axis(fig[row, col], xlabel="Time (s)", ylabel="Amplitude (μV)", title="$ch")
#             
#             # Plot original and repaired ERPs
#             lines!(ax, time_points, erp_orig[!, ch], color=:blue, linewidth=2, label="Original")
#             lines!(ax, time_points, erp_repaired[!, ch], color=:green, linewidth=2, label="Repaired")
#             
#             # Add legend for first subplot
#             if i == 1
#                 axislegend(ax, position=:topright)
#             end
#         end
#     end
#     
#     return fig
# end

# """
#     plot_artifact_summary(artifacts::EpochRejectionInfo)
# 
# Plot summary of artifact detection results.
# 
# # Arguments
# - `artifacts::EpochRejectionInfo`: Artifact detection results
# 
# # Returns
# - `Figure`: Makie figure showing artifact summary
# """
# function plot_artifact_summary(artifacts::EpochRejectionInfo)
#     using Makie
#     
#     # Create figure
#     fig = Figure(size=(800, 600))
#     
#     # Artifact counts by type
#     ax1 = Axis(fig[1, 1], xlabel="Detection Method", ylabel="Number of Rejections", title="Artifact Detection Summary")
#     
#     methods = ["Z-Variance", "Z-Max", "Z-Min", "Z-Abs", "Z-Range", "Z-Kurtosis"]
#     counts = [length(artifacts.z_variance), length(artifacts.z_max), length(artifacts.z_min), 
#               length(artifacts.z_abs), length(artifacts.z_range), length(artifacts.z_kurtosis)]
#     
#     if artifacts.abs_criterion !== nothing
#         push!(methods, "Abs-Threshold")
#         push!(counts, length(artifacts.absolute_threshold))
#     end
#     
#     barplot!(ax1, 1:length(methods), counts, color=:steelblue)
#     ax1.xticks = (1:length(methods), methods)
#     
#     # Overall statistics
#     ax2 = Axis(fig[1, 2], title="Overall Statistics")
#     
#     stats_text = """
#     Total Epochs: $(artifacts.n_epochs)
#     Rejected Epochs: $(artifacts.n_artifacts)
#     Rejection Rate: $(round(100 * artifacts.n_artifacts / artifacts.n_epochs, digits=1))%
#     
#     Z-Criterion: $(artifacts.z_criterion)
#     Abs-Criterion: $(artifacts.abs_criterion === nothing ? "None" : "$(artifacts.abs_criterion) μV")
#     """
#     
#     text!(ax2, 0.1, 0.5, text=stats_text, fontsize=12, align=(:left, :center))
#     xlims!(ax2, 0, 1)
#     ylims!(ax2, 0, 1)
#     hidespines!(ax2)
#     hidedecorations!(ax2)
#     
#     return fig
# end