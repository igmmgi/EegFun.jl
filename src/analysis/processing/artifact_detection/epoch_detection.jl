"""
    detect_bad_epochs_automatic(dat::EpochData; z_criterion::Real = 3,
                     abs_criterion::Real = 100,
                     channel_selection::Function = channels(),
                     z_measures::Vector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis])::EpochRejectionInfo

Detect bad epochs using statistical criteria and/or absolute voltage thresholds.

Supports two types of criteria for epoch rejection:
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
    interval_selection::Interval = times(),
    z_measures::AbstractVector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis],
    name::String = "rejection_info",
)::EpochRejectionInfo

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Detecting bad epochs"

    # Validate inputs
    z_criterion < 0 && @minimal_error("Z-criterion must be non-negative")
    abs_criterion < 0 && @minimal_error("Absolute criterion must be non-negative")
    z_criterion == 0 && abs_criterion == 0 && @minimal_error("One of z_criterion or abs_criterion must be > 0")

    # Validate measures
    allowed_measures = Set([:variance, :max, :min, :abs, :range, :kurtosis])
    selected_measures = Set(z_measures)
    if !issubset(selected_measures, allowed_measures)
        invalid = collect(setdiff(selected_measures, allowed_measures))
        @minimal_error("Invalid measures: $(invalid).")
    end

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    @info "Selected channels: $(_print_vector(selected_channels))"
    isempty(selected_channels) && @minimal_error("No channels selected for epoch rejection")

    # Get sample indices for the selected interval
    selected_samples = get_selected_samples(dat, interval_selection)
    @info "Selected interval: $(length(selected_samples)) samples"
    isempty(selected_samples) && @minimal_error("No samples selected for epoch rejection")

    # Calculate metrics and identify rejected epochs
    metrics = _calculate_epoch_metrics(dat, selected_channels, selected_samples, Float64(z_criterion), Float64(abs_criterion))

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
    info = EpochInfo(dat.condition, dat.condition_name, n_epochs(dat))

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
    get_rejected(info::EpochRejectionInfo)::Vector{Rejection}

Get the list of rejected channel-epoch pairs from the rejection info.

# Examples
```julia
info = detect_bad_epochs_automatic(epochs)
rejected = get_rejected(info)
```
"""
get_rejected(info::EpochRejectionInfo)::Vector{Rejection} = info.rejected
get_rejected(info::Vector{EpochRejectionInfo})::Vector{Vector{Rejection}} = get_rejected.(info)



"""
Calculate statistical metrics for all epochs.

Returns a Dict with keys: `:z_variance`, `:z_max`, `:z_min`, `:z_abs`, `:z_range`, 
`:z_kurtosis`, `:absolute_threshold`. Each value is a Dict mapping channel symbols 
to vectors of epoch indices that exceeded the criteria.
"""
function _calculate_epoch_metrics(
    dat::EpochData,
    selected_channels::AbstractVector{Symbol},
    selected_samples::AbstractVector{<:Integer},
    z_criterion::Real,
    abs_criterion::Real,
)::Dict{Symbol,Dict{Symbol,Vector{Int}}}

    metric_keys = [:z_variance, :z_max, :z_min, :z_abs, :z_range, :z_kurtosis, :absolute_threshold]
    metrics = Dict(k => Dict(ch => Int[] for ch in selected_channels) for k in metric_keys)
    isempty(dat.data) && return metrics

    first_epoch = first(dat.data)
    col_map = Dict(ch => columnindex(first_epoch, ch) for ch in selected_channels)

    # Use a UnitRange when sample indices are contiguous 
    samples_idx = if selected_samples isa AbstractUnitRange
        selected_samples
    elseif length(selected_samples) > 1 &&
           issorted(selected_samples) &&
           (selected_samples[end] - selected_samples[1] + 1 == length(selected_samples))
        selected_samples[1]:selected_samples[end]
    else
        selected_samples
    end

    Threads.@threads for ch in selected_channels
        col_idx = col_map[ch]
        # Zero-allocation views into each epoch's column
        channel_data_all = [@views(epoch[!, col_idx][samples_idx]) for epoch in dat]

        if abs_criterion > 0
            abs_threshold_violations = findall(epoch_data -> maximum(abs, epoch_data) > abs_criterion, channel_data_all)
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
