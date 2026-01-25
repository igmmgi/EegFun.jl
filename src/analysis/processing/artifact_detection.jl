"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Float64, channel_in::Symbol, channel_out::Symbol)

Detects EOG (electrooculogram) onsets in the EEG data based on a specified criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Real`: The threshold for detecting EOG onsets.
- `channel_in::Symbol`: The channel from which to detect EOG onsets.
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event.

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG) # Detect vertical EOG onsets
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG) # Detect horizontal EOG onsets
```
"""

function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol; step_size::Int = 20)
    @info "Detecting EOG onsets in channel $(channel_in) with stepsize criterion $(criterion) μV"
    channel_in ∉ propertynames(dat.data) && @minimal_error_throw("channel $(channel_in) not found in data")

    step_size = div(dat.sample_rate, step_size)
    eog_diff = diff(dat.data[1:step_size:end, channel_in])
    eog_idx = findall(x -> abs(x) >= criterion, eog_diff)
    eog_idx = [idx for (i, idx) in enumerate(eog_idx) if i == 1 || (idx - eog_idx[i-1] > 2)] .* step_size
    dat.data[!, channel_out] .= false
    dat.data[eog_idx, channel_out] .= true
    return nothing
end


"""
    detect_eog_signals!(dat::EegData, eog_cfg::Dict)

Detect EOG onsets for both vertical and horizontal EOG channels based on configuration.

# Arguments
- `dat::EegData`: The EEG data object
- `eog_cfg::Dict`: EOG configuration dictionary containing vEOG and hEOG settings

# Example
```julia
eog_cfg = Dict(
    "vEOG_criterion" => 50.0,
    "hEOG_criterion" => 50.0,
    "vEOG_channels" => [["Fp1"], ["IO1"], ["vEOG"]],
    "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]]
)
detect_eog_signals!(dat, eog_cfg)
```
"""
function detect_eog_signals!(dat::EegData, eog_cfg::Dict)
    vEOG_cfg = eog_cfg["vEOG_channels"]
    detect_eog_onsets!(dat, eog_cfg["vEOG_criterion"], Symbol(vEOG_cfg[3][1]), Symbol("is_" * vEOG_cfg[3][1]))
    hEOG_cfg = eog_cfg["hEOG_channels"]
    detect_eog_onsets!(dat, eog_cfg["hEOG_criterion"], Symbol(hEOG_cfg[3][1]), Symbol("is_" * hEOG_cfg[3][1]))
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
extreme_mask = _is_extreme_value(signal, 50.0) # Detect values ±50 μV
```
"""
_is_extreme_value(signal::AbstractVector{Float64}, threshold::Real) = abs.(signal) .> threshold

"""
    _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Real)

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
function _is_extreme_value!(mask::Vector{Bool}, signal::AbstractVector{Float64}, threshold::Real)
    @assert length(mask) == length(signal) "Mask and signal must have the same length"
    @inbounds for i in eachindex(signal)
        mask[i] = abs(signal[i]) > threshold
    end
    return nothing
end

"""
    is_extreme_value!(dat::SingleDataFrameEeg, threshold::Real; 
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
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)

    mode ∉ [:separate, :combined] && @minimal_error_throw("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error_throw("threshold must be greater than 0")

    results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

    if mode == :combined  # any channel
        channel_out = something(channel_out, Symbol("is_extreme_value_$(threshold)"))
        dat.data[!, channel_out] = falses(nrow(dat.data))
        for (_, extreme_mask) in results
            dat.data[!, channel_out] .|= extreme_mask
        end
    elseif mode == :separate # separate columns for each channel
        for (ch, extreme_mask) in results
            column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
            dat.data[!, column_name] = extreme_mask
        end
    end

    return nothing
end

"""
    is_extreme_value!(dat::MultiDataFrameEeg, threshold::Real; 
                     channel_selection::Function = channels(), 
                     sample_selection::Function = samples(),
                     epoch_selection::Function = epochs(),
                     mode::Symbol = :combined,
                     channel_out::Union{Symbol, Nothing} = nothing)

Detect extreme values across selected channels in multi-DataFrame EEG data (e.g., EpochData).

For each epoch DataFrame in the data, detects extreme values across selected channels and
adds results to that epoch's DataFrame. Works similarly to the SingleDataFrameEeg version
but processes each epoch separately.

# Arguments
- `dat::MultiDataFrameEeg`: The EEG data object (e.g., EpochData)
- `threshold::Real`: Threshold for extreme value detection
- `channel_selection::Function`: Channel predicate for selecting channels (default: all layout channels)
- `sample_selection::Function`: Sample predicate for selecting samples (default: all samples)
- `epoch_selection::Function`: Epoch predicate for selecting epochs (default: all epochs)
- `mode::Symbol`: Mode of operation - `:combined` (single combined column, default) or `:separate` (separate columns per channel)
- `channel_out::Union{Symbol, Nothing}`: Output channel name for combined mode (default: auto-generated as `is_extreme_value_<threshold>`)

# Modifies
- `dat`: Adds extreme value detection columns to each epoch DataFrame

# Examples
```julia
# Detect extreme values in all epochs (default)
is_extreme_value!(epoch_data, 100)

# Detect extreme values only in specific epochs
is_extreme_value!(epoch_data, 100, epoch_selection = epochs([1, 3, 5]))

# Detect with custom output channel name
is_extreme_value!(epoch_data, 100, channel_out = :is_artifact_value_100)
```
"""
function is_extreme_value!(
    dat::MultiDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)

    mode ∉ [:separate, :combined] && @minimal_error_throw("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error_throw("threshold must be greater than 0")

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error_throw("No channels selected for extreme value detection")

    # Get selected epochs
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error_throw("No epochs selected for extreme value detection")

    # Use provided channel_out or generate default name
    channel_out = something(channel_out, Symbol("is_extreme_value_$(threshold)"))

    # Process each selected epoch
    for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]

        # Get selected samples for this epoch
        selected_samples = get_selected_samples(epoch_df, sample_selection)

        # Initialize artifact flag column for this epoch
        epoch_df[!, channel_out] = falses(nrow(epoch_df))

        # Create sample mask once (same for all channels)
        sample_mask = falses(nrow(epoch_df))
        sample_mask[selected_samples] .= true

        if mode == :combined
            for ch in selected_channels
                extreme_mask = _is_extreme_value(epoch_df[!, ch], threshold) .& sample_mask
                epoch_df[!, channel_out] .|= extreme_mask
            end
        else
            for ch in selected_channels
                extreme_mask = _is_extreme_value(epoch_df[!, ch], threshold) .& sample_mask
                column_name = Symbol("is_extreme_value_$(ch)_$(threshold)")
                epoch_df[!, column_name] = extreme_mask
            end
        end
    end

    return nothing
end

"""
    is_extreme_value!(epochs_list::Vector{EpochData}, threshold::Real; kwargs...)

Detect extreme values across multiple EpochData objects (conditions).

Applies extreme value detection to each EpochData object in the vector. This is useful
when processing multiple experimental conditions.

# Arguments
- `epochs_list::Vector{EpochData}`: Vector of EpochData objects (one per condition)
- `threshold::Int`: Threshold for extreme value detection
- `kwargs...`: Same keyword arguments as `is_extreme_value!(dat::MultiDataFrameEeg, ...)`

# Examples
```julia
# Detect extreme values in all conditions
is_extreme_value!(epochs_cleaned, 100)

# Detect with custom output channel name
is_extreme_value!(epochs_cleaned, 100, channel_out = :is_artifact_value_100)
```
"""
is_extreme_value!(dat::Vector{EpochData}, threshold::Real; kwargs...) = is_extreme_value!.(dat, threshold; kwargs...)

# Helper function to detect extreme values for selected channels
function _detect_extreme_values(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error_throw("No channels selected for extreme value detection")

    selected_samples = get_selected_samples(dat, sample_selection)
    isempty(selected_samples) && @minimal_error_throw("No samples selected for extreme value detection")

    results = Dict{Symbol,Vector{Bool}}()

    for ch in selected_channels
        extreme_mask = _is_extreme_value(dat.data[!, ch], Float64(threshold))

        # Apply sample selection - only keep extreme values for selected samples
        sample_mask = falses(length(extreme_mask))
        sample_mask[selected_samples] .= true
        extreme_mask = extreme_mask .& sample_mask

        results[ch] = extreme_mask
    end

    return results
end

"""
    is_extreme_value(dat::SingleDataFrameEeg, threshold::Real; 
                    channel_selection::Function = channels(), 
                    sample_selection::Function = samples(),
                    mode::Symbol = :combined)

Detect extreme values across selected channels and return results.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for extreme value detection
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
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    mode::Symbol = :combined,
)

    mode ∉ [:separate, :combined] && @minimal_error_throw("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error_throw("threshold must be greater than 0")

    results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

    if mode == :combined

        combined_mask = Vector{Bool}(falses(nrow(dat.data)))
        # Combine results from all channels (OR operation)
        for (_, extreme_mask) in results
            combined_mask .|= extreme_mask
        end

        return combined_mask

    elseif mode == :separate
        # Separate mode - create temporary data object and use mutating version
        temp_dat = copy(dat)
        is_extreme_value!(temp_dat, threshold; channel_selection, sample_selection, mode = :separate)

        # Extract the extreme value columns in the same order as the original channels
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        extreme_cols = [Symbol("is_extreme_value_$(ch)_$(threshold)") for ch in selected_channels]
        return temp_dat.data[:, extreme_cols]
    end
end

"""
    n_extreme_value(dat::SingleDataFrameEeg, threshold::Real; 
                   channel_selection::Function = channels(), 
                   sample_selection::Function = samples(),
                   mode::Symbol = :combined)

Count the number of extreme values across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Real`: Threshold for extreme value detection
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
    threshold::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    mode::Symbol = :combined,
)

    mode ∉ [:separate, :combined] && @minimal_error_throw("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error_throw("threshold must be greater than 0")

    results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

    if mode == :combined
        combined_mask = Vector{Bool}(falses(nrow(dat.data)))
        for (_, extreme_mask) in results
            combined_mask .|= extreme_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        counts = [sum(results[ch]) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_extreme = counts)
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
struct Rejection
    label::Symbol
    epoch::Int
end

Base.show(io::IO, r::Rejection) = print(io, "Rejection(:$(r.label), $(r.epoch))")

# Are two Rejections equal?
is_equal_rejection(a::Rejection, b::Rejection) = a.label == b.label && a.epoch == b.epoch

function unique_rejections(rejections::Vector{Rejection})
    seen = Set{Tuple{Symbol,Int}}()
    out = Rejection[]
    for rejection in rejections
        key = (rejection.label, rejection.epoch)
        if key ∉ seen
            push!(seen, key)
            push!(out, rejection)
        end
    end
    return out
end

unique_channels(rejections::Vector{Rejection}) = unique(map(x -> x.label, rejections))
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
- `rejected::Vector{Rejection}`: All rejected epochs
- `z_criterion::Real`: Z-score criterion used for rejection
- `z_rejections::Union{ZScoreRejectionInfo, Nothing}`: Z-score based rejection info (Nothing if z_criterion = 0)
- `abs_rejections::Union{Vector{Rejection}, Nothing}`: Rejections due to absolute voltage threshold (Nothing if abs_criterion = 0)
- `abs_criterion::Real`: Absolute voltage threshold (μV) used for rejection
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

# Methods depending on EpochRejectionInfo must be defined after the struct
unique_rejections(info::EpochRejectionInfo) = unique_rejections(info.rejected)
function unique_rejections(infos::Vector{EpochRejectionInfo})
    results = Vector{Rejection}[]
    for info in infos
        push!(results, unique_rejections(info))
    end
    return results
end

unique_channels(rejections::EpochRejectionInfo) = unique_channels(rejections.rejected)
unique_channels(info::Vector{EpochRejectionInfo}) = unique_channels.(info)

unique_epochs(rejections::EpochRejectionInfo) = unique_epochs(rejections.rejected)
unique_epochs(info::Vector{EpochRejectionInfo}) = unique_epochs.(info)


"""
    detect_bad_epochs_automatic(dat::EpochData; z_criterion::Real = 3,
                     abs_criterion::Real = 100,
                     channel_selection::Function = channels(),
                     z_measures::Vector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis])::EpochRejectionInfo

Detect bad epochs using statistical criteria and/or absolute voltage thresholds.

This function can use two types of criteria for epoch rejection:
1. **Z-score criteria**: Calculates six statistical measures for each epoch (variance, max, min,
   absolute max, range, kurtosis) across selected channels. For each measure, the maximum
   across channels is taken for each epoch, then z-scored across epochs. Epochs exceeding
   the z-criterion for any selected measure are rejected. You can choose which measures to
   apply using `measures` (default: `[:variance, :max, :min, :abs, :range, :kurtosis]`).
2. **Absolute voltage criteria**: Epochs with any channel exceeding the absolute voltage
   threshold (in μV) are rejected.

# Arguments
- `dat::EpochData`: Epoched EEG data to process
- `z_criterion::Real`: Z-score threshold for rejection (default: 3.0). Set to 0 to disable z-score based rejection.
- `abs_criterion::Real`: Absolute voltage threshold in μV for rejection (default: 100.0). Set to 0 to disable absolute threshold rejection.
- `channel_selection::Function`: Channel predicate for selecting channels to analyze (default: all channels)
- `z_measures::Vector{Symbol}`: Which z-score measures to apply (default: all: `[:variance, :max, :min, :abs, :range, :kurtosis]`)

# Returns
- `EpochRejectionInfo`: Information about which epochs were rejected and why

# Requirements
- At least one of `z_criterion` or `abs_criterion` must be greater than 0

# Examples
```julia
using EegFun, JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2")

# Use default criteria (z_criterion=3.0, abs_criterion=100.0)
rejection_info = detect_bad_epochs_automatic(epochs)

# Customize z-score criterion only
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 2.0)

# Customize absolute voltage threshold only
rejection_info = detect_bad_epochs_automatic(epochs, abs_criterion = 80.0)  # 80 μV

# Use both criteria with custom values
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 2.5, abs_criterion = 80.0)

# Disable z-score rejection (use only absolute threshold)
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100.0)

# Disable absolute threshold rejection (use only z-score)
rejection_info = detect_bad_epochs_automatic(epochs, z_criterion = 3.0, abs_criterion = 0)

# Check results
println("Original epochs: \$(rejection_info.n_epochs)")
println("Rejected epochs: \$(rejection_info.n_artifacts)")
println("Rejected epochs: \$(rejection_info.rejected)")
```

# Notes
- Default z-criterion: 3.0 
- Default abs-criterion: 100 μV 
- Common z-criteria: 2.0 (more aggressive), 2.5, 3.0 (more conservative)
- Common absolute criteria: 50-100 μV 
- Set either criterion to 0 to disable that type of rejection
- Rejection is based on ANY metric exceeding the criterion
- All metrics are calculated independently and combined with OR logic
"""
function detect_bad_epochs_automatic(
    dat::EpochData;
    z_criterion::Real = 3,
    abs_criterion::Real = 100,
    channel_selection::Function = channels(),
    z_measures::Vector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis],
    name::String = "rejection_info",
)::EpochRejectionInfo

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Detecting bad epochs"

    # Validate inputs
    z_criterion < 0 && @minimal_error_throw("Z-criterion must be non-negative")
    abs_criterion < 0 && @minimal_error_throw("Absolute criterion must be non-negative")
    z_criterion == 0 && abs_criterion == 0 && @minimal_error_throw("One of z_criterion or abs_criterion must be > 0")

    # Validate measures
    allowed_measures = Set([:variance, :max, :min, :abs, :range, :kurtosis])
    selected_measures = Set(z_measures)
    if !issubset(selected_measures, allowed_measures)
        invalid = collect(setdiff(selected_measures, allowed_measures))
        @minimal_error_throw("Invalid measures: $(invalid).")
    end

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    @info "Selected channels: $(print_vector(selected_channels))"
    isempty(selected_channels) && @minimal_error_throw("No channels selected for epoch rejection")

    # Calculate metrics and identify rejected epochs
    metrics = _calculate_epoch_metrics(dat, selected_channels, Float64(z_criterion), Float64(abs_criterion))

    # Initialize rejection lists (all needed for EpochRejectionInfo struct)
    rejected_info = Rejection[]
    z_variance = Rejection[]
    z_max = Rejection[]
    z_min = Rejection[]
    z_abs = Rejection[]
    z_range = Rejection[]
    z_kurtosis = Rejection[]

    # Map measure directly to (rejection_list, metrics_key)
    measure_to_list_and_key = Dict(
        :variance => (z_variance, :z_variance),
        :max => (z_max, :z_max),
        :min => (z_min, :z_min),
        :abs => (z_abs, :z_abs),
        :range => (z_range, :z_range),
        :kurtosis => (z_kurtosis, :z_kurtosis),
    )

    # Build Rejection objects directly from metric results
    # Handle z-score metrics
    if z_criterion > 0
        for measure in z_measures
            rejection_list, metrics_key = measure_to_list_and_key[measure]
            for channel in selected_channels
                for epoch_idx in metrics[metrics_key][channel]
                    rejection = Rejection(channel, epoch_idx)
                    push!(rejected_info, rejection)
                    push!(rejection_list, rejection)
                end
            end
        end
    end

    # Handle absolute threshold
    abs_rejections = abs_criterion > 0 ? Rejection[] : nothing
    if abs_criterion > 0
        for channel in selected_channels
            for epoch_idx in metrics[:absolute_threshold][channel]
                rejection = Rejection(channel, epoch_idx)
                push!(rejected_info, rejection)
                push!(abs_rejections, rejection)
            end
        end
    end

    # Create z-score rejection info if z_criterion > 0
    z_rejections = z_criterion > 0 ? ZScoreRejectionInfo(z_measures, z_variance, z_max, z_min, z_abs, z_range, z_kurtosis) : nothing

    # Create rejection info
    info = EpochInfo(dat.condition, dat.condition_name, length(dat.data))

    rejection_info = EpochRejectionInfo(
        name,
        info,
        length(rejected_info),
        abs_criterion,
        abs_rejections,
        z_criterion,
        z_rejections,
        unique_rejections(rejected_info),
        nothing,  # repaired - populated during repair
        nothing,  # skipped - populated during repair
    )

    return rejection_info
end

detect_bad_epochs_automatic(dat::Vector{EpochData}; kwargs...) = detect_bad_epochs_automatic.(dat; kwargs...)


"""
    get_rejected(state::EpochRejectionState)::Vector{Int}

Get indices of rejected epochs from the rejection state.

# Examples
```julia
state = detect_bad_epochs_interactive(epochs)
rejected_indices = get_rejected(state)
```
"""
get_rejected(info::EpochRejectionInfo)::Vector{Rejection} = info.rejected

"""
    get_rejected(state::Vector{EpochRejectionState})::Vector{Int}

Get indices of rejected epochs from the rejection state.

# Examples
```julia
state = detect_bad_epochs_interactive(epochs)
rejected_indices = get_rejected(state)
```
"""
get_rejected(info::Vector{EpochRejectionInfo})::Vector{Vector{Rejection}} = get_rejected.(info)

"""
Validate inputs for epoch rejection.
"""
function _validate_rejection_inputs(dat::EpochData, z_criterion::Real)
    isempty(dat.data) && @minimal_error_throw("Cannot reject epochs from empty EpochData")
    z_criterion <= 0 && @minimal_error_throw("Z-criterion must be positive")
    length(dat.data) < 3 && @minimal_warning "Only $(length(dat.data)) epochs available."
end


"""
Calculate statistical metrics for all epochs.

Returns a Dict with keys: :z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis, :absolute_threshold
Each value is a Dict mapping channel symbols to vectors of epoch indices that exceeded the criteria.
"""
function _calculate_epoch_metrics(
    dat::EpochData,
    selected_channels::Vector{Symbol},
    z_criterion::Real,
    abs_criterion::Real,
)::Dict{Symbol,Dict{Symbol,Vector{Int}}}

    metric_keys = [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis, :absolute_threshold]
    metrics = Dict(k => Dict(ch => Int[] for ch in selected_channels) for k in metric_keys)

    num_epochs = n_epochs(dat)
    for ch in selected_channels
        channel_data_all = Vector{Vector{Float64}}(undef, num_epochs)
        for (i, epoch) in enumerate(dat.data)
            channel_data_all[i] = epoch[!, ch]
        end

        if abs_criterion > 0
            abs_threshold_violations = findall(epoch_data -> maximum(abs.(epoch_data)) > abs_criterion, channel_data_all)
            append!(metrics[:absolute_threshold][ch], abs_threshold_violations)
        end

        if z_criterion > 0

            variances = var.(channel_data_all)
            max_values = maximum.(channel_data_all)
            min_values = minimum.(channel_data_all)
            abs_values = [maximum(abs, epoch_data) for epoch_data in channel_data_all]
            ranges = max_values .- min_values
            kurtoses = kurtosis.(channel_data_all)

            z_scores = zscore.([variances, max_values, min_values, abs_values, ranges, kurtoses])
            z_metric_keys = [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis]

            for (z_score, metric_key) in zip(z_scores, z_metric_keys)
                bad_epochs = findall(abs.(z_score) .> z_criterion)
                append!(metrics[metric_key][ch], bad_epochs)
            end

        end

    end

    return metrics
end


# ===================
# REPORTING FUNCTIONS
# ===================

"""
    Base.show(io::IO, info::EpochRejectionInfo)

Display rejection information in a human-readable format.
"""
function Base.show(io::IO, info::EpochRejectionInfo)
    println(io, "EpochRejectionInfo: $(info.name)")
    println(io, "Condition: $(info.info.number): $(info.info.name)")
    println(io, "  Abs criterion: $(info.abs_criterion > 0 ? string(info.abs_criterion,  " μV") : "disabled")")
    println(io, "  Z-criterion: $(info.z_criterion > 0 ? string(info.z_criterion) : "disabled")")
    println(io, "  Epochs total: $(info.info.n), Epochs rejected: $(length(unique_epochs(info.rejected)))")
    println(io, "  Artifacts total: $(info.n_artifacts)")
    println(io, "  Rejected epochs: $(print_vector(unique_epochs(info.rejected)))")

    if info.abs_rejections !== nothing
        println(io, "  Rejection breakdown (absolute):")
        println(
            io,
            "    Abs threshold: $(length(unique_epochs(info.abs_rejections))) unique epochs, $(length(unique_channels(info.abs_rejections))) unique channels",
        )
    end

    if info.z_rejections !== nothing
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

Base.show(io::IO, infos::Vector{EpochRejectionInfo}) = Base.show.(Ref(io), infos)




"""
    repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Repair detected artifacts using the specified method.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use

# Available Methods
- `:neighbor_interpolation` - Weighted neighbor interpolation (default). Uses `dat.layout.neighbours` for neighbor information.
- `:spherical_spline` - Spherical spline interpolation

# Keyword Arguments (for :spherical_spline method)
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)

# Returns
Nothing (mutates dat in-place)
"""
function repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)

    if method ∉ [:neighbor_interpolation, :spherical_spline]
        throw(ArgumentError("Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"))
    end

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Repairing artifacts using method: $method"
    if method == :neighbor_interpolation
        # Determine which channels can be repaired if not already done
        if isnothing(artifacts.repaired)
            channel_repairable!(artifacts, dat.layout)
        end
        repair_artifacts_neighbor!(dat, artifacts; kwargs...)
    elseif method == :spherical_spline
        repair_artifacts_spherical_spline!(dat, artifacts; kwargs...)
    end
    return nothing
end


"""
    repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Non-mutating version of repair_artifacts!. Creates a copy of the data and repairs artifacts without modifying the original.

# Arguments
- `dat::EpochData`: The epoch data to repair (NOT modified)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use (default: :neighbor_interpolation)

# Keyword Arguments
Same as repair_artifacts! for the respective methods.

# Returns
- `EpochData`: A new EpochData object with repaired artifacts

# Examples
```julia
# Basic artifact repair (creates new object)
repaired_epochs = repair_artifacts(epochs, artifacts)

# Using spherical spline method
repaired_epochs = repair_artifacts(epochs, artifacts, :spherical_spline)

# Rejecting bad epochs entirely (use reject_epochs instead)
clean_epochs = reject_epochs(epochs, artifacts)
```
"""
function repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)
    dat_copy = copy(dat)
    repair_artifacts!(dat_copy, artifacts; method, kwargs...)
    return dat_copy
end

function repair_artifacts(dat::Vector{EpochData}, artifacts::Vector{EpochRejectionInfo}; kwargs...)
    return repair_artifacts.(dat, artifacts; kwargs...)
end

function repair_artifacts!(dat::Vector{EpochData}, artifacts::Vector{EpochRejectionInfo}; kwargs...)
    repair_artifacts!.(dat, artifacts; kwargs...)
    return nothing
end



"""
    channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)

Analyze which channels can be repaired and which cannot, based on neighbor availability.
Populates `artifacts.repaired` and `artifacts.skipped` with the analysis.

# Arguments
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts (mutated to add repair analysis)
- `layout::Layout`: Layout object containing neighbor information

# Returns
- `EpochRejectionInfo`: The same artifacts object (modified in-place)

# Notes
This function only analyzes repairability - it does not perform any repairs.
Use `repair_artifacts_neighbor!` to actually perform the repairs after this analysis.
"""
function channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)
    rejected = unique([r.epoch for r in artifacts.rejected])

    # Initialize tracking dictionaries in artifacts struct (OrderedDict to maintain sorted order)
    artifacts.repaired = OrderedDict{Int,Vector{Symbol}}()
    artifacts.skipped = OrderedDict{Int,Vector{Symbol}}()

    # Process epochs in sorted order to maintain ordering in OrderedDict
    for epoch_idx in sort(rejected)
        bad_channels = [artifact.label for artifact in artifacts.rejected if artifact.epoch == epoch_idx]
        isempty(bad_channels) && continue

        repairable_channels = check_channel_neighbors(bad_channels, layout)

        if isempty(repairable_channels)
            if length(bad_channels) == 1
                @info "Epoch $epoch_idx: Cannot repair channel $(bad_channels[1]) (fewer than 2 neighbors)"
            else
                @info "Epoch $epoch_idx: Cannot repair channels $(bad_channels) (bad neighbors and/or fewer than 2 neighbors)"
            end
            # Track that all were skipped
            artifacts.skipped[epoch_idx] = bad_channels
        else
            skipped_channels = setdiff(bad_channels, repairable_channels)
            artifacts.repaired[epoch_idx] = repairable_channels
            if !isempty(skipped_channels)
                @info "Epoch $epoch_idx: Skipping repair of $(length(skipped_channels)) channel(s) with bad neighbors: $skipped_channels"
                artifacts.skipped[epoch_idx] = skipped_channels
            end
        end
    end

    return artifacts
end

channel_repairable!(artifacts::Vector{EpochRejectionInfo}, layout::Layout) = channel_repairable!.(artifacts, Ref(layout))


"""
    repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)

Repair artifacts using weighted neighbor interpolation.
Uses `artifacts.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Uses `dat.layout.neighbours` for neighbor information.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information with `repaired` already populated

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)
    # Check if repaired has been populated
    if isnothing(artifacts.repaired)
        throw(ArgumentError("repaired not populated. Call channel_repairable! first."))
    end

    if isempty(artifacts.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return dat
    end

    # Process epochs in sorted order (already sorted in OrderedDict)
    for (epoch_idx, repairable_channels) in artifacts.repaired
        @info "Epoch $epoch_idx: Repairing channels $(repairable_channels) using neighbor interpolation"

        # Use unified channel repair function with epoch selection
        repair_channels!(
            dat,
            repairable_channels;
            method = :neighbor_interpolation,
            epoch_selection = epochs([epoch_idx]),
            neighbours_dict = dat.layout.neighbours,
        )
        @info "" # formatting
    end

    return nothing
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
function repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int = 4, lambda::Float64 = 1e-5)

    _ensure_coordinates_3d!(dat.layout)

    # Get all rejected epochs with their bad channels
    rejected = unique([r.epoch for r in artifacts.rejected])

    for epoch_idx in rejected
        # Get bad channels for this epoch
        bad_channels = [artifact.label for artifact in artifacts.rejected if artifact.epoch == epoch_idx]
        isempty(bad_channels) && continue

        @info "Repairing epoch $epoch_idx channels $(bad_channels) using spherical spline interpolation"

        # Use unified channel repair function with epoch selection
        repair_channels!(dat, bad_channels; method = :spherical_spline, epoch_selection = epochs([epoch_idx]), m = m, lambda = lambda)
    end

    return dat
end



"""
    subset_bad_data(data_path::String, threshold::Float64; 
                    subset_directory::String = "excluded")

Identify and move files from participants with low data retention to a separate directory.

Participants are considered "bad" if they have less than the threshold percentage
of data remaining in ANY condition. All files associated with bad participants
are moved to a subdirectory (default: "excluded").

The function searches for "epoch_summary.jld2" in the specified directory.

# Arguments
- `data_path::String`: Path to the directory containing epoch_summary.jld2 and preprocessed files
- `threshold::Float64`: Minimum percentage threshold (e.g., 75.0 means 75% retention required)

# Keyword Arguments
- `subset_directory::String`: Name of subdirectory for excluded participants (default: "excluded")

# Examples
```julia
# Move participants with < 75% data in any condition to "excluded" subdirectory
subset_bad_data("preprocessed_files", 75.0)

# Specify custom subset directory name
subset_bad_data("/path/to/preprocessed", 80.0, subset_directory="excluded")
```
"""
function subset_bad_data(data_path::String, threshold::Float64; subset_directory::String = "excluded")

    # Validate inputs/outputs
    (threshold < 0.0 || threshold > 100.0) && @minimal_error_throw("threshold must be 0 < threshold < 100, got $threshold")
    !isdir(data_path) && @minimal_error_throw("data_path must be a directory: $data_path")

    epoch_summary_path = joinpath(data_path, "epoch_summary.jld2")
    !isfile(epoch_summary_path) && @minimal_error_throw("epoch_summary.jld2 not found: $data_path")

    # Load epoch summary and check required columns
    epoch_summary = load_data(epoch_summary_path)
    missing_cols = [col for col in [:file, :percentage] if !(col in propertynames(epoch_summary))]
    !isempty(missing_cols) && @minimal_error_throw("epoch_summary missing required columns: $(missing_cols)")

    # Use data_path as output directory
    output_directory = abspath(data_path)
    subset_dir_path = joinpath(output_directory, subset_directory)
    !isdir(subset_dir_path) && mkpath(subset_dir_path)

    # Find participants with any condition below threshold
    bad_participants = unique(epoch_summary.file[epoch_summary.percentage.<threshold])
    println("Subsetting data: $(length(bad_participants))")
    println("   N remaining: $(length(unique(epoch_summary.file)) - length(bad_participants))")
    println("   N removed: $(length(bad_participants))")

    # Create subset summary files with only non-excluded participants
    # Filter epoch_summary to exclude bad participants
    epoch_summary_subset = epoch_summary[.!in.(epoch_summary.file, Ref(bad_participants)), :]

    # Load and filter file_summary_subset
    file_summary_path = joinpath(data_path, "file_summary.jld2")
    !isfile(file_summary_path) && @minimal_error_throw("file_summary.jld2 not found: $data_path")

    file_summary = load_data(file_summary_path)
    file_summary_subset = file_summary[.!in.(file_summary.file, Ref(bad_participants)), :]

    # Save subset summary files
    epoch_summary_subset_path = joinpath(output_directory, "epoch_summary_subset.jld2")
    file_summary_subset_path = joinpath(output_directory, "file_summary_subset.jld2")

    jldsave(epoch_summary_subset_path; data = epoch_summary_subset)
    jldsave(file_summary_subset_path; data = file_summary_subset)

    # Find and move all files for bad participants
    all_files = readdir(output_directory)
    for participant_id in bad_participants
        # Find all files that contain participant_id
        matching_files = filter(all_files) do filename
            occursin(participant_id, filename)
        end
        # Move each matching file
        for filename in matching_files
            src_path = joinpath(output_directory, filename)
            dst_path = joinpath(subset_dir_path, filename)
            if isfile(src_path)
                mv(src_path, dst_path, force = true)
            end
        end
    end
end
