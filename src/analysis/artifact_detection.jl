# =============================================================================
# ARTIFACT DETECTION
# =============================================================================

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

    # Get the input channel data
    input_data = dat.data[!, channel_in]

    # Detect onsets (positive-going threshold crossings)
    onsets = _is_extreme_value(input_data, Float64(criterion))

    # Add the output channel to the data
    dat.data[!, channel_out] = onsets

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
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Union{Vector{Bool}, Nothing}`: Boolean vector for sample selection (default: nothing - all samples)
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
    sample_selection::Union{Vector{Bool},Nothing} = nothing,
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

        # Initialize the output channel with false values
        dat.data[!, output_channel] = falses(nrow(dat.data))

        # Combine results from all channels (OR operation)
        for (ch, extreme_mask) in results
            dat.data[!, output_channel] .|= extreme_mask
        end

    else  # mode == :separate
        # Separate mode - use specified channel selection
        results = _detect_extreme_values(dat, threshold; channel_selection, sample_selection)

        # Add results to data object
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
    sample_selection::Union{Vector{Bool},Nothing} = nothing,
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for extreme value detection")
    end

    results = Dict{Symbol,Vector{Bool}}()

    for ch in selected_channels
        channel_data = dat.data[!, ch]
        extreme_mask = _is_extreme_value(channel_data, Float64(threshold))

        # Apply sample selection if provided
        if sample_selection !== nothing
            extreme_mask = extreme_mask .& sample_selection
        end

        results[ch] = extreme_mask
    end

    return results
end

"""
    is_extreme_value(dat::SingleDataFrameEeg, threshold::Int; 
                    channel_selection::Function = channels(), 
                    sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                    mode::Symbol = :combined)

Detect extreme values across selected channels and return results.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Union{Vector{Bool}, Nothing}`: Boolean vector for sample selection (default: nothing - all samples)
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
    sample_selection::Union{Vector{Bool},Nothing} = nothing,
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
                   sample_selection::Union{Vector{Bool}, Nothing} = nothing,
                   mode::Symbol = :combined)

Count the number of extreme values across selected channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `threshold::Int`: Threshold for extreme value detection
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Union{Vector{Bool}, Nothing}`: Boolean vector for sample selection (default: nothing - all samples)
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
    sample_selection::Union{Vector{Bool},Nothing} = nothing,
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

    else  # mode == :separate
        # Separate mode - count for each channel
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

"""
    EpochRejectionInfo

Stores information about which epochs were rejected and why.

# Fields
- `n_original::Int`: Number of epochs before rejection
- `n_remaining::Int`: Number of epochs after rejection
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
    n_original::Int
    n_remaining::Int
    rejected_epochs::Vector{Int}
    rejected_by_z_variance::Vector{Int}
    rejected_by_z_max::Vector{Int}
    rejected_by_z_min::Vector{Int}
    rejected_by_z_abs::Vector{Int}
    rejected_by_z_range::Vector{Int}
    rejected_by_z_kurtosis::Vector{Int}
    rejected_by_abs_threshold::Vector{Int}
    z_criterion::Float64
    abs_criterion::Union{Float64, Nothing}
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
println("Original epochs: \$(rejection_info.n_original)")
println("Remaining epochs: \$(rejection_info.n_remaining)")
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

    @info "Starting epoch detection with z-criterion = $z_criterion" * 
          (abs_criterion === nothing ? "" : ", abs_criterion = $abs_criterion μV")
    
    # Validate inputs
    _validate_rejection_inputs(dat, z_criterion)
    
    if abs_criterion !== nothing && abs_criterion <= 0
        @minimal_error_throw("Absolute criterion must be positive, got $abs_criterion")
    end
    
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for epoch rejection")
    end
    
    @info "Using $(length(selected_channels)) channels for rejection analysis"
    
    # Calculate metrics for all epochs
    n_epochs = length(dat.data)
    metrics = _calculate_epoch_metrics(dat, selected_channels)
    
    # Identify rejected epochs for each z-score metric
    rejection_results = _identify_rejected_epochs(metrics, z_criterion)
    
    # Identify epochs rejected by absolute threshold
    abs_rejected = Int[]
    if abs_criterion !== nothing
        abs_rejected = _identify_abs_threshold_rejected(dat, selected_channels, Float64(abs_criterion))
    end
    
    # Combine all rejected epochs
    rejected_epochs = sort(unique(vcat(
        rejection_results[:variance],
        rejection_results[:max],
        rejection_results[:min],
        rejection_results[:abs],
        rejection_results[:range],
        rejection_results[:kurtosis],
        abs_rejected
    )))
    
    n_rejected = length(rejected_epochs)
    n_remaining = n_epochs - n_rejected

    # Create rejection info before modifying data
    rejection_info = EpochRejectionInfo(
        n_epochs,
        n_remaining,
        rejected_epochs,
        rejection_results[:variance],
        rejection_results[:max],
        rejection_results[:min],
        rejection_results[:abs],
        rejection_results[:range],
        rejection_results[:kurtosis],
        abs_rejected,
        Float64(z_criterion),
        abs_criterion === nothing ? nothing : Float64(abs_criterion)
    )
    @info rejection_info
    
    # # Remove rejected epochs
    # if !isempty(rejected_epochs)
    #     epochs_to_keep = setdiff(1:n_epochs, rejected_epochs)
    #     dat.data = dat.data[epochs_to_keep]
    #     @info "Removed $n_rejected epochs. $(length(dat.data)) epochs remaining."
    # else
    #     @info "No epochs rejected."
    # end
    
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
function _calculate_epoch_metrics(dat::EpochData, selected_channels::Vector{Symbol})::Dict{Symbol, Vector{Float64}}
    n_epochs = length(dat.data)
    
    # Pre-allocate arrays for metrics
    variance_max = zeros(Float64, n_epochs)
    max_vals = zeros(Float64, n_epochs)
    min_vals = zeros(Float64, n_epochs)
    abs_max = zeros(Float64, n_epochs)
    range_vals = zeros(Float64, n_epochs)
    kurtosis_max = zeros(Float64, n_epochs)
    
    # Calculate metrics for each epoch
    for (epoch_idx, epoch) in enumerate(dat.data)
        # For each channel, calculate metric, then take max across channels
        channel_variance = Float64[]
        channel_max = Float64[]
        channel_min = Float64[]
        channel_abs = Float64[]
        channel_range = Float64[]
        channel_kurtosis = Float64[]
        
        for ch in selected_channels
            channel_data = epoch[!, ch]
            
            push!(channel_variance, var(channel_data))
            push!(channel_max, maximum(channel_data))
            push!(channel_min, minimum(channel_data))
            push!(channel_abs, maximum(abs.(channel_data)))
            push!(channel_range, maximum(channel_data) - minimum(channel_data))
            push!(channel_kurtosis, kurtosis(channel_data))
        end
        
        # Take maximum across channels for each metric
        variance_max[epoch_idx] = maximum(channel_variance)
        max_vals[epoch_idx] = maximum(channel_max)
        min_vals[epoch_idx] = minimum(channel_min)  # Note: minimum of minimums!
        abs_max[epoch_idx] = maximum(channel_abs)
        range_vals[epoch_idx] = maximum(channel_range)
        kurtosis_max[epoch_idx] = maximum(channel_kurtosis)
    end
    
    return Dict{Symbol, Vector{Float64}}(
        :variance => variance_max,
        :max => max_vals,
        :min => min_vals,
        :abs => abs_max,
        :range => range_vals,
        :kurtosis => kurtosis_max
    )
end


"""
Identify epochs that exceed the z-criterion for each metric.

Returns a Dict with keys: :variance, :max, :min, :abs, :range, :kurtosis
Each value is a vector of epoch indices that were rejected by that criterion.
"""
function _identify_rejected_epochs(metrics::Dict{Symbol, Vector{Float64}}, z_criterion::Real)::Dict{Symbol, Vector{Int}}
    results = Dict{Symbol, Vector{Int}}()
    
    for (metric_name, values) in metrics
        # Calculate z-scores
        z_scores = zscore(values)
        
        # Find epochs exceeding criterion
        rejected = findall(z_scores .> z_criterion)
        
        results[metric_name] = rejected
    end
    
    return results
end


"""
Identify epochs that exceed the absolute voltage threshold.

Returns a vector of epoch indices that contain any sample exceeding the threshold.
"""
function _identify_abs_threshold_rejected(dat::EpochData, selected_channels::Vector{Symbol}, abs_threshold::Float64)::Vector{Int}
    rejected = Int[]
    
    for (epoch_idx, epoch) in enumerate(dat.data)
        # Check if any sample in any selected channel exceeds threshold
        for ch in selected_channels
            channel_data = epoch[!, ch]
            if any(abs.(channel_data) .> abs_threshold)
                push!(rejected, epoch_idx)
                break  # Don't need to check other channels once one exceeds
            end
        end
    end
    
    return rejected
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
    println(io, "  Original epochs: $(info.n_original)")
    println(io, "  Remaining epochs: $(info.n_remaining)")
    println(io, "  Rejected epochs: $(length(info.rejected_epochs))")
    println(io, "")
    println(io, "  Rejection breakdown (z-score):")
    println(io, "    Z-Variance:  $(length(info.rejected_by_z_variance)) epochs")
    println(io, "    Z-Maximum:   $(length(info.rejected_by_z_max)) epochs")
    println(io, "    Z-Minimum:   $(length(info.rejected_by_z_min)) epochs")
    println(io, "    Z-Absolute:  $(length(info.rejected_by_z_abs)) epochs")
    println(io, "    Z-Range:     $(length(info.rejected_by_z_range)) epochs")
    println(io, "    Z-Kurtosis:  $(length(info.rejected_by_z_kurtosis)) epochs")
    if info.abs_criterion !== nothing
        println(io, "")
        println(io, "  Rejection breakdown (absolute):")
        println(io, "    Abs threshold: $(length(info.rejected_by_abs_threshold)) epochs")
    end
    
    if length(info.rejected_epochs) <= 20
        println(io, "")
        println(io, "  Rejected epoch indices: $(info.rejected_epochs)")
    end
end


"""
    save_rejection_report(info::EpochRejectionInfo, filename::String)

Save rejection information to a text file.

# Arguments
- `info::EpochRejectionInfo`: Rejection information to save
- `filename::String`: Output filename (typically ends in .txt)

# Examples
```julia
save_rejection_report(rejection_info, "participant_1_rejection_report.txt")
```
"""
function save_rejection_report(info::EpochRejectionInfo, filename::String)
    open(filename, "w") do io
        println(io, repeat("=", 70))
        println(io, "AUTOMATIC EPOCH REJECTION REPORT")
        println(io, repeat("=", 70))
        println(io, "")
        println(io, "Z-criterion: $(info.z_criterion)")
        if info.abs_criterion !== nothing
            println(io, "Abs criterion: $(info.abs_criterion) μV")
        end
        println(io, "Original epochs: $(info.n_original)")
        println(io, "Remaining epochs: $(info.n_remaining)")
        println(io, "Rejected epochs: $(length(info.rejected_epochs)) ($(round(100 * length(info.rejected_epochs) / info.n_original, digits=1))%)")
        println(io, "")
        println(io, "REJECTION BREAKDOWN (Z-SCORE):")
        println(io, repeat("-", 70))
        println(io, "Z-Variance:  $(length(info.rejected_by_z_variance)) epochs - $(info.rejected_by_z_variance)")
        println(io, "Z-Maximum:   $(length(info.rejected_by_z_max)) epochs - $(info.rejected_by_z_max)")
        println(io, "Z-Minimum:   $(length(info.rejected_by_z_min)) epochs - $(info.rejected_by_z_min)")
        println(io, "Z-Absolute:  $(length(info.rejected_by_z_abs)) epochs - $(info.rejected_by_z_abs)")
        println(io, "Z-Range:     $(length(info.rejected_by_z_range)) epochs - $(info.rejected_by_z_range)")
        println(io, "Z-Kurtosis:  $(length(info.rejected_by_z_kurtosis)) epochs - $(info.rejected_by_z_kurtosis)")
        if info.abs_criterion !== nothing
            println(io, "")
            println(io, "REJECTION BREAKDOWN (ABSOLUTE):")
            println(io, repeat("-", 70))
            println(io, "Abs threshold: $(length(info.rejected_by_abs_threshold)) epochs - $(info.rejected_by_abs_threshold)")
        end
        println(io, "")
        println(io, "ALL REJECTED EPOCHS:")
        println(io, repeat("-", 70))
        println(io, info.rejected_epochs)
        println(io, "")
        println(io, repeat("=", 70))
    end
    
    @info "Rejection report saved to: $filename"
end


# #=============================================================================
#     BATCH PROCESSING FUNCTIONS
# =============================================================================#
# 
# """Generate default output directory name for rejection operation."""
# function _default_rejection_output_dir(input_dir::String, pattern::String, z_criterion::Real)
#     z_str = replace(string(z_criterion), "." => "p")  # Replace . with p for filename
#     joinpath(input_dir, "rejected_z$(z_str)_$(pattern)")
# end
# 
# 
# """
# Process a single epoch file through automatic rejection pipeline.
# Returns BatchResult with success/failure info.
# """
# function _process_rejection_file(
#     filepath::String,
#     output_path::String,
#     z_criterion::Real,
#     channel_selection::Function,
# )
#     filename = basename(filepath)
#     
#     # Load data
#     file_data = load(filepath)
#     
#     # Try common variable names for epoched data
#     epoch_var_names = ["epochs", "epoch_data", "data"]
#     epochs_data = nothing
#     
#     for var_name in epoch_var_names
#         if haskey(file_data, var_name)
#             epochs_data = file_data[var_name]
#             break
#         end
#     end
#     
#     if isnothing(epochs_data)
#         return BatchResult(false, filename, "No epoched data variable found (tried: $(epoch_var_names))")
#     end
#     
#     if !(epochs_data isa EpochData)
#         return BatchResult(false, filename, "Data is not EpochData type")
#     end
#     
#     # Reject epochs
#     try
#         rejection_info = reject_epochs_automatic!(epochs_data, z_criterion; channel_selection = channel_selection)
#         
#         # Save cleaned epochs
#         save(output_path, "epochs", epochs_data)
#         
#         # Save rejection report in same directory as output file
#         output_dir = dirname(output_path)
#         base_filename = splitext(basename(filename))[1]
#         report_path = joinpath(output_dir, "$(base_filename)_rejection_report.txt")
#         save_rejection_report(rejection_info, report_path)
#         
#         message = "Rejected $(length(rejection_info.rejected_epochs)) of $(rejection_info.n_original) epochs"
#         return BatchResult(true, filename, message)
#     catch e
#         return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
#     end
# end


# """
#     reject_epochs_automatic(file_pattern::String, z_criterion::Real;
#                            input_dir::String = pwd(),
#                            channel_selection::Function = channels(),
#                            participants::Union{Int, Vector{Int}, Nothing} = nothing,
#                            output_dir::Union{String, Nothing} = nothing)
# 
# Batch process epoch files with automatic rejection and save to a new directory.
# 
# This function processes multiple participant files at once, automatically rejecting
# epochs based on statistical criteria.
# 
# # Arguments
# - `file_pattern::String`: Pattern to match files (e.g., "epochs", "epochs_cleaned")
# - `z_criterion::Real`: Z-score threshold for rejection (typically 2.0 or 3.0)
# - `input_dir::String`: Input directory containing JLD2 files (default: current directory)
# - `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
# - `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
# - `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)
# 
# # Examples
# ```julia
# # Reject epochs with z > 2.0 for all participants
# reject_epochs_automatic("epochs", 2.0)
# 
# # Specific participants only
# reject_epochs_automatic("epochs", 2.5, participants = [1, 2, 3])
# 
# # Custom output directory
# reject_epochs_automatic("epochs", 3.0, output_dir = "/path/to/output")
# 
# # Only use specific channels for rejection analysis
# reject_epochs_automatic("epochs", 2.0,
#                        channel_selection = channels(x -> startswith.(string.(x), "F")))
# 
# # Full example
# reject_epochs_automatic("epochs", 2.5,
#                        input_dir = "/data/study1",
#                        participants = 1:20)
# ```
# 
# # Output
# - Creates new directory with cleaned epoch data files
# - Each output file contains "epochs" variable with cleaned EpochData
# - Rejection reports saved as *_rejection_report.txt files
# - Log file saved to output directory
# 
# # Notes
# - Common z-criteria: 2.0 (aggressive), 2.5 (moderate), 3.0 (conservative)
# - Higher z-criterion = fewer epochs rejected = more lenient
# - Lower z-criterion = more epochs rejected = more strict
# """
# function reject_epochs_automatic(
#     file_pattern::String,
#     z_criterion::Real;
#     input_dir::String = pwd(),
#     channel_selection::Function = channels(),
#     participants::Union{Int,Vector{Int},Nothing} = nothing,
#     output_dir::Union{String,Nothing} = nothing,
# )
#     
#     # Setup logging
#     log_file = "reject_epochs_automatic.log"
#     setup_global_logging(log_file)
#     
#     try
#         @info "Batch automatic epoch rejection started at $(now())"
#         @log_call "reject_epochs_automatic" (file_pattern, z_criterion)
#         
#         # Validation
#         if (error_msg = _validate_input_dir(input_dir)) !== nothing
#             @minimal_error_throw(error_msg)
#         end
#         
#         if z_criterion <= 0
#             @minimal_error_throw("Z-criterion must be positive, got $z_criterion")
#         end
#         
#         # Setup directories
#         output_dir = something(output_dir, _default_rejection_output_dir(input_dir, file_pattern, z_criterion))
#         mkpath(output_dir)
#         
#         # Find files
#         files = _find_batch_files(file_pattern, input_dir; participants)
#         
#         if isempty(files)
#             @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
#             return nothing
#         end
#         
#         @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
#         @info "Z-criterion: $z_criterion"
#         @info "Channel selection: $(channel_selection == channels() ? "all channels" : "custom")"
#         
#         # Create processing function with captured parameters
#         process_fn = (input_path, output_path) ->
#             _process_rejection_file(input_path, output_path, z_criterion, channel_selection)
#         
#         # Execute batch operation
#         results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Rejecting epochs")
#         
#         _log_batch_summary(results, output_dir)
#         
#     finally
#         _cleanup_logging(log_file, output_dir)
#     end
# end

