"""
Representational Similarity Analysis (RSA) functions.

This module provides functions for computing Representational Dissimilarity
Matrices (RDMs) from EEG data and comparing them to model RDMs.
"""

# ==============================================================================
#   DATA PREPARATION (reuse from decoding)
# ==============================================================================

"""
    _prepare_rsa_data(epochs::Vector{EpochData})

Prepare epoch data for RSA analysis.

Extracts data from multiple EpochData conditions (already subsetted).
Returns arrays ready for RSA computation.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData, one per condition (already subsetted)

# Returns
- `data_arrays::Vector{Array{Float64, 3}}`: Vector of [channels × time × trials] arrays, one per condition
- `times::Vector{Float64}`: Time points in seconds
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
"""
function _prepare_rsa_data(epochs::Vector{EpochData})
    channels = channel_labels(epochs)
    times = time(epochs)

    # Prepare data arrays for each condition
    data_arrays = Vector{Array{Float64,3}}()
    n_trials_per_condition = Int[]

    for epoch_data in epochs
        n_trials = length(epoch_data.data)
        push!(n_trials_per_condition, n_trials)

        # Preallocate: [channels × time × trials]
        condition_data = zeros(Float64, length(channels), length(times), n_trials)

        for (trial_idx, trial_df) in enumerate(epoch_data.data)
            # Extract channel data (data is already subsetted)
            for (ch_idx, ch_name) in enumerate(channels)
                ch_data = trial_df[!, ch_name]
                condition_data[ch_idx, :, trial_idx] = ch_data
            end
        end

        push!(data_arrays, condition_data)
    end

    return data_arrays, times, n_trials_per_condition
end

"""
    _compute_pooled_covariance(data_arrays::Vector{Array{Float64, 3}}, time_idx::Int)

Compute pooled covariance matrix across all conditions at a specific time point.

The pooled covariance is the weighted average of within-condition covariances,
which provides a robust estimate of the noise covariance structure.

# Arguments
- `data_arrays::Vector{Array{Float64, 3}}`: Vector of [channels × time × trials] arrays, one per condition
- `time_idx::Int`: Time point index to compute covariance for

# Returns
- `pooled_cov::Matrix{Float64}`: Pooled covariance matrix [features × features]

# Notes
- Uses unbiased estimator (dividing by n-1)
- Regularization is applied if matrix is near-singular
"""
function _compute_pooled_covariance(data_arrays::Vector{Array{Float64,3}}, time_idx::Int)
    n_conditions = length(data_arrays)
    n_features = size(data_arrays[1], 1)  # Number of channels

    # Collect all trials across conditions at this time point
    all_trials = Matrix{Float64}[]

    for cond_data in data_arrays
        n_trials = size(cond_data, 3)
        # Extract [features × trials] at time_idx
        timepoint_data = cond_data[:, time_idx, :]  # [channels × trials]
        push!(all_trials, timepoint_data)
    end

    # Compute pooled covariance
    # Pool all trials together for robust covariance estimation
    all_data = hcat(all_trials...)  # [features × total_trials]

    # Compute covariance matrix
    pooled_cov = cov(all_data', corrected = true)  # Transpose to [trials × features]

    # Add small regularization for numerical stability
    # This helps with near-singular matrices
    λ = 1e-6 * tr(pooled_cov) / n_features
    pooled_cov += λ * I(n_features)

    return pooled_cov
end

# ==============================================================================
#   DISSIMILARITY MEASURES
# ==============================================================================

"""
    _compute_dissimilarity(
        pattern1::Vector{Float64},
        pattern2::Vector{Float64},
        measure::Symbol;
        covariance_matrix::Union{Matrix{Float64}, Nothing} = nothing
    )

Compute dissimilarity between two neural patterns.

# Arguments
- `pattern1::Vector{Float64}`: First pattern (e.g., average across trials for condition 1)
- `pattern2::Vector{Float64}`: Second pattern (e.g., average across trials for condition 2)
- `measure::Symbol`: Dissimilarity measure (:correlation, :euclidean, :mahalanobis)
- `covariance_matrix::Union{Matrix{Float64}, Nothing}`: Pooled covariance matrix for Mahalanobis distance (optional)

# Returns
- `dissimilarity::Float64`: Dissimilarity value (higher = more dissimilar)
"""
function _compute_dissimilarity(
    pattern1::Vector{Float64},
    pattern2::Vector{Float64},
    measure::Symbol;
    covariance_matrix::Union{Matrix{Float64},Nothing} = nothing,
)
    if measure == :correlation || measure == :pearson
        # 1 - Pearson correlation (dissimilarity)
        corr = cor(pattern1, pattern2)
        return 1.0 - corr
    elseif measure == :spearman
        # 1 - Spearman correlation (dissimilarity)
        corr = StatsBase.corspearman(pattern1, pattern2)
        return 1.0 - corr
    elseif measure == :euclidean
        # Euclidean distance
        return sqrt(sum((pattern1 .- pattern2) .^ 2))
    elseif measure == :mahalanobis
        # Mahalanobis distance with pooled covariance matrix
        # Note: The pooled covariance is automatically computed in the rsa() function
        # when dissimilarity_measure == :mahalanobis, so this error should never occur
        # in normal usage. This is a safety check for direct function calls.
        if isnothing(covariance_matrix)
            @minimal_error_throw(
                "Mahalanobis distance requires a covariance matrix. " *
                "This error should not occur when using rsa() or rsa_crossvalidated() functions, " *
                "as they automatically compute pooled covariance. " *
                "If calling _compute_dissimilarity() directly, provide covariance_matrix parameter."
            )
        end

        # Compute Mahalanobis distance: sqrt((x-y)' * inv(Σ) * (x-y))
        diff = pattern1 .- pattern2

        # Use pseudo-inverse for numerical stability
        try
            inv_cov = pinv(covariance_matrix)
            mahal_dist = sqrt(max(0.0, dot(diff, inv_cov * diff)))
            return mahal_dist
        catch e
            @minimal_warning "Mahalanobis distance computation failed (singular covariance?), using Euclidean instead"
            return sqrt(sum(diff .^ 2))
        end
    else
        @minimal_error_throw("Unknown dissimilarity measure: $measure. Use :correlation, :spearman, :euclidean, or :mahalanobis")
    end
end

"""
    _compute_rdm(
        condition_patterns::Vector{Matrix{Float64}},
        measure::Symbol;
        covariance_matrix::Union{Matrix{Float64}, Nothing} = nothing
    )

Compute Representational Dissimilarity Matrix (RDM) from condition patterns.

# Arguments
- `condition_patterns::Vector{Matrix{Float64}}`: Vector of [channels × time] matrices, one per condition
- `measure::Symbol`: Dissimilarity measure
- `covariance_matrix::Union{Matrix{Float64}, Nothing}`: Pooled covariance matrix for Mahalanobis distance (optional)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]
"""
function _compute_rdm(
    condition_patterns::Vector{Matrix{Float64}},
    measure::Symbol;
    covariance_matrix::Union{Matrix{Float64},Nothing} = nothing,
)
    n_conditions = length(condition_patterns)
    rdm = zeros(Float64, n_conditions, n_conditions)

    # Only compute upper triangle (excluding diagonal) for efficiency
    # RDM is symmetric, so we can fill lower triangle by copying
    for i = 1:n_conditions
        rdm[i, i] = 0.0  # Diagonal is always 0 (self-dissimilarity)
        for j = (i+1):n_conditions
            # Flatten patterns to vectors for comparison
            pattern_i = vec(condition_patterns[i])
            pattern_j = vec(condition_patterns[j])
            dissim = _compute_dissimilarity(pattern_i, pattern_j, measure; covariance_matrix = covariance_matrix)
            rdm[i, j] = dissim
            rdm[j, i] = dissim  # Symmetric
        end
    end

    return rdm
end

"""
    normalize_rdm(rdm::Matrix{Float64}; method::Symbol = :none)

Normalize a Representational Dissimilarity Matrix (RDM).

Different normalization schemes can affect correlation results and interpretability.
Choose based on your analysis goals and the properties of your dissimilarity measure.

# Arguments
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]
- `method::Symbol`: Normalization method
  - `:none` - No normalization (default)
  - `:zscore` - Z-score normalization (mean=0, std=1)
  - `:rank` - Rank transformation (converts to ordinal ranks)
  - `:minmax` - Min-max normalization (scales to [0, 1])

# Returns
- `normalized_rdm::Matrix{Float64}`: Normalized RDM matrix

# Examples
```julia
# Z-score normalization (useful for comparing RDMs with different scales)
rdm_z = normalize_rdm(rdm, method=:zscore)

# Rank normalization (robust to outliers, used for Spearman-like comparisons)
rdm_rank = normalize_rdm(rdm, method=:rank)

# Min-max normalization (scales to [0, 1])
rdm_minmax = normalize_rdm(rdm, method=:minmax)
```

# Notes
- Normalization is applied to the upper triangle (excluding diagonal)
- Diagonal remains zero after normalization
- For `:rank`, ties are handled using average ranks
"""
function normalize_rdm(rdm::Matrix{Float64}; method::Symbol = :none)
    if method == :none
        return copy(rdm)
    end

    n = size(rdm, 1)
    normalized_rdm = copy(rdm)

    # Extract upper triangle (excluding diagonal) for normalization
    upper_indices = [CartesianIndex(i, j) for i = 1:n for j = (i+1):n]
    upper_values = [rdm[idx] for idx in upper_indices]

    if method == :zscore
        # Z-score normalization: (x - mean) / std
        μ = mean(upper_values)
        σ = std(upper_values)
        if σ > 0
            normalized_values = (upper_values .- μ) ./ σ
        else
            # All values are the same, return zeros
            normalized_values = zeros(Float64, length(upper_values))
        end
    elseif method == :rank
        # Rank transformation (ordinal ranks)
        normalized_values = StatsBase.ordinalrank(upper_values)
    elseif method == :minmax
        # Min-max normalization: (x - min) / (max - min)
        min_val = minimum(upper_values)
        max_val = maximum(upper_values)
        if max_val > min_val
            normalized_values = (upper_values .- min_val) ./ (max_val - min_val)
        else
            # All values are the same, return zeros
            normalized_values = zeros(Float64, length(upper_values))
        end
    else
        @minimal_error_throw("Unknown normalization method: $method. Use :none, :zscore, :rank, or :minmax")
    end

    # Fill normalized values back into matrix (both upper and lower triangles)
    for (idx_pos, idx) in enumerate(upper_indices)
        i, j = idx.I
        normalized_rdm[i, j] = normalized_values[idx_pos]
        normalized_rdm[j, i] = normalized_values[idx_pos]  # Symmetric
    end

    # Ensure diagonal is zero
    for i = 1:n
        normalized_rdm[i, i] = 0.0
    end

    return normalized_rdm
end

# ==============================================================================
#   RSA ANALYSIS
# ==============================================================================

"""
    rsa(
        epochs::Vector{EpochData};
        channel_selection::Function = channels(),
        sample_selection::Function = samples(),
        dissimilarity_measure::Symbol = :correlation,
        average_trials::Bool = true,
        normalize_method::Symbol = :none,
    )

Perform Representational Similarity Analysis (RSA) on epoch data.

Computes Representational Dissimilarity Matrices (RDMs) at each time point,
showing how dissimilar neural patterns are between different conditions.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData, one per condition, from a SINGLE participant
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels([:Fz, :Cz, :Pz])` for specific channels
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels()` for all channels (default)
- `sample_selection::Function=samples()`: Sample selection predicate. See `samples()` for options.
  - Example: `sample_selection=samples((-0.2, 0.8))` for time window from -0.2 to 0.8 seconds
  - Example: `sample_selection=samples()` for all time points (default)
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)
  - `:correlation` or `:pearson` - 1 - Pearson correlation (default, most common)
  - `:spearman` - 1 - Spearman rank correlation (robust to outliers)
  - `:euclidean` - Euclidean distance
  - `:mahalanobis` - Mahalanobis distance (accounts for covariance structure, automatically computed)
- `average_trials::Bool`: Whether to average across trials before computing RDM (default: true)
- `normalize_method::Symbol`: RDM normalization method (default: :none)
  - `:none` - No normalization
  - `:zscore` - Z-score normalization (mean=0, std=1)
  - `:rank` - Rank transformation (ordinal ranks)
  - `:minmax` - Min-max normalization (scales to [0, 1])

# Returns
- `RsaData`: Object containing RSA results

# Examples
```julia
# Basic RSA with correlation-based dissimilarity
epochs = [epoch_condition1, epoch_condition2, epoch_condition3]
rsa_result = rsa(epochs, channel_selection=channels([:Fz, :Cz, :Pz]), sample_selection=samples((-0.2, 0.8)))

# Using Euclidean distance
rsa_result = rsa(epochs, dissimilarity_measure=:euclidean)

# Using Mahalanobis distance (accounts for covariance structure)
rsa_result = rsa(epochs, dissimilarity_measure=:mahalanobis)

# Without averaging trials (computes RDM for each trial, then averages RDMs)
rsa_result = rsa(epochs, average_trials=false)

# With z-score normalization
rsa_result = rsa(epochs, normalize_method=:zscore)

# Use all channels and all time points (default)
rsa_result = rsa(epochs)
```
"""
function rsa(
    epochs::Vector{EpochData};
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    dissimilarity_measure::Symbol = :correlation,
    average_trials::Bool = true,
    normalize_method::Symbol = :none,
)
    # Prepare and validate data using shared helper
    data_arrays, times, n_trials_per_condition, condition_names, selected_channels, first_epoch =
        _prepare_and_validate_rsa(epochs, channel_selection, sample_selection)

    n_conditions = length(condition_names)
    n_timepoints = length(times)

    # Preallocate RDM array: [time × condition × condition]
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    # Compute RDM at each time point
    for t = 1:n_timepoints
        if average_trials
            # Use shared helper to compute RDMs from data
            rdms = _compute_rdms_from_data(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
            break # _compute_rdms_from_data handles the loop
        else
            # Compute pooled covariance if using Mahalanobis distance
            pooled_cov = nothing
            if dissimilarity_measure == :mahalanobis
                pooled_cov = _compute_pooled_covariance(data_arrays, t)
            end

            # Compute RDM for each trial, then average RDMs
            # This approach is more robust to noise and outliers
            n_trials = minimum(n_trials_per_condition)
            trial_rdms = Vector{Matrix{Float64}}()

            for trial_idx = 1:n_trials
                condition_patterns = Matrix{Float64}[]
                for cond_data in data_arrays
                    if trial_idx <= size(cond_data, 3)
                        # Extract [channels] at time t, trial trial_idx
                        pattern = vec(cond_data[:, t, trial_idx])
                        push!(condition_patterns, reshape(pattern, length(selected_channels), 1))
                    end
                end
                if length(condition_patterns) == n_conditions
                    trial_rdm = _compute_rdm(condition_patterns, dissimilarity_measure; covariance_matrix = pooled_cov)
                    push!(trial_rdms, trial_rdm)
                end
            end

            # Average RDMs across trials (proper matrix averaging)
            if !isempty(trial_rdms)
                # Sum all RDM matrices
                avg_rdm = zeros(Float64, n_conditions, n_conditions)
                for trial_rdm in trial_rdms
                    avg_rdm .+= trial_rdm
                end
                # Divide by number of trials to get average
                avg_rdm ./= length(trial_rdms)
                rdms[t, :, :] = avg_rdm
            end
        end
    end

    # Apply normalization if requested using shared helper
    _normalize_rdms!(rdms, normalize_method)

    # Create RsaData object
    rsa_result = RsaData(
        first_epoch.file,
        condition_names,
        times,
        rdms,
        dissimilarity_measure,
        selected_channels,
        first_epoch.layout,
        first_epoch.sample_rate,
        # Default analysis_info if not present
        analysis_info = first_epoch.analysis_info,
    )

    return rsa_result
end

# ==============================================================================
#   MODEL COMPARISON
# ==============================================================================

"""
    compare_models(
        rsa_data::RsaData,
        model_rdms::Union{Vector{Matrix{Float64}}, Vector{Array{Float64, 3}}, Vector{Union{Matrix{Float64}, Array{Float64, 3}}}};
        model_names::Union{Vector{String}, Nothing} = nothing,
        correlation_type::Symbol = :spearman,
        n_permutations::Int = 1000,
        rng::AbstractRNG = Random.GLOBAL_RNG,
    )

Compare neural RDMs to model RDMs (static or temporal).

Computes correlations between neural RDMs and model RDMs at each time point,
with optional permutation-based significance testing.

# Input Format

**Model RDMs must be in one of these formats:**

1. **Static RDM**: `Matrix{Float64}` [n_conditions × n_conditions]
   - Single RDM used at all time points
   - Example: Reaction time RDM, similarity ratings RDM

2. **Temporal RDM**: `Array{Float64, 3}` [n_timepoints × n_conditions × n_conditions]
   - RDM at each time point
   - Must have same number of timepoints as `rsa_data.times`
   - Example: Eye tracking RDM, EDA RDM

**Note**: Convert your raw data to RDMs first using helper functions:
- Static: `create_rdm_from_reaction_times()`, `create_rdm_from_vectors()`, etc.
- Temporal: `create_temporal_rdm()` (handles temporal alignment automatically)

# Arguments
- `rsa_data::RsaData`: RSA results from neural data (temporal)
- `model_rdms`: Vector of model RDMs
  - Static: `Matrix{Float64}` [condition × condition]
  - Temporal: `Array{Float64, 3}` [time × condition × condition]
  - Mixed: Can contain both static and temporal models
- `model_names::Union{Vector{String}, Nothing}`: Names for models (default: Model1, Model2, ...)
- `correlation_type::Symbol`: Type of correlation (:spearman, :pearson)
- `n_permutations::Int`: Number of permutations for significance testing (0 = no testing)
- `rng::AbstractRNG`: Random number generator

# Returns
- `RsaData`: Updated RsaData with model correlations and p-values

# Examples

## Static Model (Reaction Times)

```julia
# Step 1: Convert your data to RDM
rts = [0.3, 0.5, 0.4]  # Your reaction time data
rt_rdm = create_rdm_from_reaction_times(rts)

# Step 2: Compare
neural_rsa = rsa(epochs)
rsa_with_model = compare_models(neural_rsa, [rt_rdm], model_names=["RTs"])
```

## Temporal Model (Eye Tracking)

```julia
# Step 1: Convert your data to temporal RDM (with automatic alignment)
eye_data = Array{Float64, 3}  # [conditions × features × time] - your format
eye_times = Vector{Float64}     # Your time vector
eye_rdms = create_temporal_rdm(eye_data, eye_times; align_to=neural_rsa)

# Step 2: Compare
neural_rsa = rsa(epochs)
rsa_with_model = compare_models(neural_rsa, [eye_rdms], model_names=["Eye Tracking"])
```

## Multiple Models (Mixed Static and Temporal)

```julia
# Convert all your model data to RDMs
rt_rdm = create_rdm_from_reaction_times(rts)
eye_rdms = create_temporal_rdm(eye_data, eye_times; align_to=neural_rsa)
eda_rdms = create_temporal_rdm(eda_data, eda_times; align_to=neural_rsa)

# Compare all at once
rsa_with_models = compare_models(
    neural_rsa,
    [rt_rdm, eye_rdms, eda_rdms],
    model_names=["RTs", "Eye Tracking", "EDA"]
)
```

## Pre-computed RDMs

If you already have RDMs computed (from other analyses):

```julia
# Your pre-computed RDM
my_rdm = Matrix{Float64}  # [conditions × conditions]

# Use directly - no conversion needed
rsa_with_model = compare_models(neural_rsa, [my_rdm], model_names=["My Model"])
```
"""
function compare_models(
    rsa_data::RsaData,
    model_rdms::Union{Vector{Matrix{Float64}},Vector{Array{Float64,3}},Vector{Union{Matrix{Float64},Array{Float64,3}}}};
    model_names::Union{Vector{String},Nothing} = nothing,
    correlation_type::Symbol = :spearman,
    n_permutations::Int = 1000,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    n_models = length(model_rdms)
    n_conditions = length(rsa_data.condition_names)
    n_times = length(rsa_data.times)

    # Validate model RDMs and determine which are temporal
    is_temporal = _validate_model_rdms(model_rdms, n_conditions, n_times)

    # Generate model names if not provided
    if isnothing(model_names)
        model_names = ["Model$(i)" for i = 1:n_models]
    elseif length(model_names) != n_models
        @minimal_error_throw("Number of model names ($(length(model_names))) doesn't match number of models ($n_models)")
    end

    # Compute correlations at each time point
    correlations, p_values = _compute_model_correlations(rsa_data, model_rdms, is_temporal, correlation_type, n_permutations, rng)

    # Update RsaData
    rsa_data.model_correlations = correlations
    rsa_data.model_names = model_names
    rsa_data.p_values = p_values

    return rsa_data
end

# ==============================================================================
#   GRAND AVERAGE
# ==============================================================================

"""
    grand_average(rsa_data_list::Vector{RsaData}; compute_noise_ceiling::Bool = true)

Compute grand average RSA results across participants.

Averages RDMs across participants at each time point. Optionally computes
noise ceiling using leave-one-out cross-validation.

# Arguments
- `rsa_data_list::Vector{RsaData}`: Vector of RsaData objects, one per participant
- `compute_noise_ceiling::Bool`: Whether to compute noise ceiling (default: true)

# Returns
- `RsaData`: Grand-averaged RSA results with noise ceiling (if computed)

# Examples
```julia
# Average across participants with noise ceiling
all_rsa = [rsa_p1, rsa_p2, rsa_p3]
grand_avg_rsa = grand_average(all_rsa)

# Without noise ceiling
grand_avg_rsa = grand_average(all_rsa, compute_noise_ceiling=false)
```
"""
function grand_average(rsa_data_list::Vector{RsaData}; compute_noise_ceiling::Bool = true)
    if isempty(rsa_data_list)
        @minimal_error_throw("Cannot compute grand average with empty RSA data list")
    end

    if length(rsa_data_list) == 1
        return rsa_data_list[1]
    end

    # Validate all have same structure
    first_rsa = rsa_data_list[1]
    n_conditions = length(first_rsa.condition_names)
    n_times = length(first_rsa.times)

    for (idx, rsa_data) in enumerate(rsa_data_list)
        if length(rsa_data.condition_names) != n_conditions
            @minimal_error_throw("RSA data $idx has $(length(rsa_data.condition_names)) conditions, expected $n_conditions")
        end
        if length(rsa_data.times) != n_times
            @minimal_error_throw("RSA data $idx has $(length(rsa_data.times)) time points, expected $n_times")
        end
        if rsa_data.dissimilarity_measure != first_rsa.dissimilarity_measure
            @minimal_warning "RSA data $idx uses dissimilarity measure $(rsa_data.dissimilarity_measure), expected $(first_rsa.dissimilarity_measure)"
        end
    end

    # Average RDMs across participants
    grand_rdm = zeros(Float64, n_times, n_conditions, n_conditions)
    for t = 1:n_times
        for i = 1:n_conditions
            for j = 1:n_conditions
                grand_rdm[t, i, j] = mean([rsa_data.rdm[t, i, j] for rsa_data in rsa_data_list])
            end
        end
    end

    # Compute noise ceiling if requested
    noise_ceiling = nothing
    if compute_noise_ceiling && length(rsa_data_list) >= 2
        noise_ceiling = EegFun.compute_noise_ceiling(rsa_data_list)
    end

    # Create grand-averaged RsaData
    grand_rsa = RsaData(
        "Grand Average",
        first_rsa.condition_names,
        first_rsa.times,
        grand_rdm,
        first_rsa.dissimilarity_measure,
        first_rsa.channels,
        first_rsa.layout,
        first_rsa.sample_rate;
        noise_ceiling = noise_ceiling,
        analysis_info = first_rsa.analysis_info,
    )

    return grand_rsa
end

# ==============================================================================
#   INTERNAL HELPERS
# ==============================================================================

"""
    _validate_model_rdms(model_rdms, n_conditions, n_times)

Validate dimensions and types of model RDMs. Returns a vector of booleans
indicating whether each model is temporal.
"""
function _validate_model_rdms(model_rdms, n_conditions, n_times)
    is_temporal = Bool[]
    for (idx, model_rdm) in enumerate(model_rdms)
        if isa(model_rdm, Array{Float64,3})
            # Temporal model: [time × condition × condition]
            if size(model_rdm, 1) != n_times
                @minimal_error_throw(
                    "Temporal model RDM $idx has $(size(model_rdm, 1)) timepoints, " * "expected $n_times to match neural data"
                )
            end
            if size(model_rdm, 2) != n_conditions || size(model_rdm, 3) != n_conditions
                @minimal_error_throw(
                    "Temporal model RDM $idx has condition dimensions $(size(model_rdm, 2))×$(size(model_rdm, 3)), " *
                    "expected $n_conditions×$n_conditions"
                )
            end
            push!(is_temporal, true)
        elseif isa(model_rdm, Matrix{Float64})
            # Static model: [condition × condition]
            if size(model_rdm) != (n_conditions, n_conditions)
                @minimal_error_throw("Static model RDM $idx has size $(size(model_rdm)), " * "expected ($n_conditions, $n_conditions)")
            end
            if !issymmetric(model_rdm)
                @minimal_warning "Static model RDM $idx is not symmetric, symmetrizing"
                model_rdms[idx] = (model_rdm + model_rdm') / 2
            end
            push!(is_temporal, false)
        else
            @minimal_error_throw(
                "Model RDM $idx has unsupported type $(typeof(model_rdm)). " *
                "Expected Matrix{Float64} (static) or Array{Float64, 3} (temporal)"
            )
        end
    end
    return is_temporal
end

"""
    _compute_model_correlations(rsa_data, model_rdms, is_temporal, correlation_type, n_permutations, rng)

Compute correlations between neural RDMs and model RDMs across time.
"""
function _compute_model_correlations(
    rsa_data::RsaData,
    model_rdms,
    is_temporal::Vector{Bool},
    correlation_type::Symbol,
    n_permutations::Int,
    rng::AbstractRNG,
)
    n_times = length(rsa_data.times)
    n_models = length(model_rdms)

    correlations = zeros(Float64, n_times, n_models)
    p_values = n_permutations > 0 ? zeros(Float64, n_times, n_models) : nothing

    for t = 1:n_times
        neural_rdm = rsa_data.rdm[t, :, :]
        neural_vec = _extract_upper_triangular(neural_rdm)

        for (model_idx, model_rdm) in enumerate(model_rdms)
            # Get model RDM for this time point
            if is_temporal[model_idx]
                model_rdm_t = model_rdm[t, :, :]
                if !issymmetric(model_rdm_t)
                    model_rdm_t = (model_rdm_t + model_rdm_t') / 2
                end
            else
                model_rdm_t = model_rdm
            end

            model_vec = _extract_upper_triangular(model_rdm_t)

            # Compute correlation
            corr = _correlate_vectors(neural_vec, model_vec, correlation_type)
            correlations[t, model_idx] = corr

            # Permutation test if requested
            if n_permutations > 0
                p_values[t, model_idx] = _run_model_permutations(neural_vec, model_vec, corr, correlation_type, n_permutations, rng)
            end
        end
    end

    return correlations, p_values
end

"""
    _run_model_permutations(neural_vec, model_vec, observed_corr, correlation_type, n_permutations, rng)

Perform permutation testing for model correlations.
"""
function _run_model_permutations(neural_vec, model_vec, observed_corr, correlation_type, n_permutations, rng)
    permuted_corrs = zeros(Float64, n_permutations)
    for perm_idx = 1:n_permutations
        shuffled_model = shuffle(rng, model_vec)
        permuted_corrs[perm_idx] = _correlate_vectors(neural_vec, shuffled_model, correlation_type)
    end
    # Two-tailed p-value
    return sum(abs.(permuted_corrs) .>= abs(observed_corr)) / n_permutations
end

"""
    _correlate_vectors(vec1, vec2, correlation_type)

Helper to compute correlation between two vectors.
"""
function _correlate_vectors(vec1, vec2, correlation_type)
    if correlation_type == :spearman
        return StatsBase.corspearman(vec1, vec2)
    elseif correlation_type == :pearson
        return cor(vec1, vec2)
    else
        @minimal_error_throw("Unknown correlation type: $correlation_type. Use :spearman or :pearson")
    end
end

"""
    _extract_upper_triangular(rdm::Matrix{Float64})

Extract upper triangular values (excluding diagonal) from an RDM.
"""
function _extract_upper_triangular(rdm::Matrix{Float64})
    n = size(rdm, 1)
    triu_indices = [CartesianIndex(i, j) for i = 1:n for j = (i+1):n]
    return [rdm[idx] for idx in triu_indices]
end

"""
    _prepare_and_validate_rsa(epochs, channel_selection, sample_selection)

Shared data preparation and validation for RSA modules.
"""
function _prepare_and_validate_rsa(epochs, channel_selection, sample_selection)
    # Input validations
    if isempty(epochs)
        @minimal_error_throw("Cannot perform RSA with empty epochs vector")
    end

    if length(epochs) < 2
        @minimal_error_throw("Need at least 2 conditions for RSA, got $(length(epochs))")
    end

    # Subset epochs by channel and sample selection
    epochs_subset = subset(epochs; channel_selection = channel_selection, sample_selection = sample_selection, include_extra = false)

    if isempty(epochs_subset) || isempty(channel_labels(epochs_subset[1]))
        @minimal_error_throw("Channel selection produced no channels")
    end

    if isempty(epochs_subset[1].data) || isempty(epochs_subset[1].data[1][!, :time])
        @minimal_error_throw("Sample selection produced no time points")
    end

    # Prepare data from subsetted epochs
    data_arrays, times, n_trials_per_condition = _prepare_rsa_data(epochs_subset)

    # Metadata extraction
    first_epoch = epochs_subset[1]
    condition_names = [e.condition_name for e in epochs_subset]
    selected_channels = channel_labels(epochs_subset)

    return data_arrays, times, n_trials_per_condition, condition_names, selected_channels, first_epoch
end

"""
    _normalize_rdms!(rdms, method)

In-place normalization of an RDM array [time × condition × condition].
"""
function _normalize_rdms!(rdms::Array{Float64,3}, method::Symbol)
    if method == :none
        return rdms
    end

    n_times = size(rdms, 1)
    for t = 1:n_times
        rdms[t, :, :] = normalize_rdm(rdms[t, :, :]; method = method)
    end
    return rdms
end

"""
    _compute_rdms_from_data(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)

Compute RDMs at each time point from data arrays [channels × time × trials].

This is a shared helper used by both `rsa()` and `rsa_crossvalidated()`.
"""
function _compute_rdms_from_data(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
)
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    for t = 1:n_timepoints
        # Compute pooled covariance if using Mahalanobis distance
        pooled_cov = nothing
        if dissimilarity_measure == :mahalanobis
            pooled_cov = _compute_pooled_covariance(data_arrays, t)
        end

        # Average across trials first, then compute RDM
        condition_patterns = Matrix{Float64}[]
        for cond_data in data_arrays
            # Extract [channels × trials] at time t
            timepoint_data = cond_data[:, t, :]  # [channels × trials]
            # Average across trials: [channels]
            avg_pattern = vec(mean(timepoint_data, dims = 2))
            push!(condition_patterns, reshape(avg_pattern, length(selected_channels), 1))
        end
        rdms[t, :, :] = _compute_rdm(condition_patterns, dissimilarity_measure; covariance_matrix = pooled_cov)
    end

    return rdms
end

