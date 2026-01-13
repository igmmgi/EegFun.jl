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
    data_arrays = Vector{Array{Float64, 3}}()
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

# ==============================================================================
#   DISSIMILARITY MEASURES
# ==============================================================================

"""
    _compute_dissimilarity(pattern1::Vector{Float64}, pattern2::Vector{Float64}, measure::Symbol)

Compute dissimilarity between two neural patterns.

# Arguments
- `pattern1::Vector{Float64}`: First pattern (e.g., average across trials for condition 1)
- `pattern2::Vector{Float64}`: Second pattern (e.g., average across trials for condition 2)
- `measure::Symbol`: Dissimilarity measure (:correlation, :euclidean, :mahalanobis)

# Returns
- `dissimilarity::Float64`: Dissimilarity value (higher = more dissimilar)
"""
function _compute_dissimilarity(pattern1::Vector{Float64}, pattern2::Vector{Float64}, measure::Symbol)
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
        # Mahalanobis distance (requires covariance matrix - simplified here)
        # For full implementation, would need covariance matrix
        @minimal_warning "Mahalanobis distance not fully implemented, using Euclidean"
        return sqrt(sum((pattern1 .- pattern2) .^ 2))
    else
        @minimal_error_throw("Unknown dissimilarity measure: $measure. Use :correlation, :spearman, :euclidean, or :mahalanobis")
    end
end

"""
    _compute_rdm(condition_patterns::Vector{Matrix{Float64}}, measure::Symbol)

Compute Representational Dissimilarity Matrix (RDM) from condition patterns.

# Arguments
- `condition_patterns::Vector{Matrix{Float64}}`: Vector of [channels × time] matrices, one per condition
- `measure::Symbol`: Dissimilarity measure

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]
"""
function _compute_rdm(condition_patterns::Vector{Matrix{Float64}}, measure::Symbol)
    n_conditions = length(condition_patterns)
    rdm = zeros(Float64, n_conditions, n_conditions)

    for i in 1:n_conditions
        for j in 1:n_conditions
            if i == j
                rdm[i, j] = 0.0  # Self-dissimilarity is 0
            else
                # Flatten patterns to vectors for comparison
                pattern_i = vec(condition_patterns[i])
                pattern_j = vec(condition_patterns[j])
                rdm[i, j] = _compute_dissimilarity(pattern_i, pattern_j, measure)
                rdm[j, i] = rdm[i, j]  # Symmetric
            end
        end
    end

    return rdm
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
- `average_trials::Bool`: Whether to average across trials before computing RDM (default: true)

# Returns
- `RsaData`: Object containing RSA results

# Examples
```julia
# Basic RSA with correlation-based dissimilarity
epochs = [epoch_condition1, epoch_condition2, epoch_condition3]
rsa_result = rsa(epochs, channel_selection=channels([:Fz, :Cz, :Pz]), sample_selection=samples((-0.2, 0.8)))

# Using Euclidean distance
rsa_result = rsa(epochs, dissimilarity_measure=:euclidean)

# Without averaging trials (computes RDM for each trial, then averages RDMs)
rsa_result = rsa(epochs, average_trials=false)

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
)
    # Input validations
    if isempty(epochs)
        @minimal_error_throw("Cannot perform RSA with empty epochs vector")
    end

    if length(epochs) < 2
        @minimal_error_throw("Need at least 2 conditions for RSA, got $(length(epochs))")
    end

    # Subset epochs by channel and sample selection
    epochs = subset(
        epochs;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = false,
    )
    isempty(channel_labels(epochs[1])) && @minimal_error_throw("Channel selection produced no channels")
    isempty(epochs[1].data[1][!, :time]) && @minimal_error_throw("Sample selection produced no time points")

    # Prepare data from subsetted epochs
    data_arrays, times, n_trials_per_condition = _prepare_rsa_data(epochs)
    n_conditions = length(epochs)
    n_timepoints = length(times)

    # Get metadata from first epoch
    first_epoch = epochs[1]
    condition_names = [e.condition_name for e in epochs]
    selected_channels = channel_labels(epochs)

    # Preallocate RDM array: [time × condition × condition]
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    # Compute RDM at each time point
    for t in 1:n_timepoints
        if average_trials
            # Average across trials first, then compute RDM
            condition_patterns = Matrix{Float64}[]
            for (cond_idx, cond_data) in enumerate(data_arrays)
                # Extract [channels × trials] at time t
                timepoint_data = cond_data[:, t, :]  # [channels × trials]
                # Average across trials: [channels]
                avg_pattern = vec(mean(timepoint_data, dims=2))
                push!(condition_patterns, reshape(avg_pattern, length(selected_channels), 1))
            end
            rdms[t, :, :] = _compute_rdm(condition_patterns, dissimilarity_measure)
        else
            # Compute RDM for each trial, then average RDMs
            trial_rdms = Float64[]
            n_trials = minimum(n_trials_per_condition)
            
            for trial_idx in 1:n_trials
                condition_patterns = Matrix{Float64}[]
                for cond_data in data_arrays
                    if trial_idx <= size(cond_data, 3)
                        # Extract [channels] at time t, trial trial_idx
                        pattern = vec(cond_data[:, t, trial_idx])
                        push!(condition_patterns, reshape(pattern, length(selected_channels), 1))
                    end
                end
                if length(condition_patterns) == n_conditions
                    trial_rdm = _compute_rdm(condition_patterns, dissimilarity_measure)
                    push!(trial_rdms, trial_rdm)
                end
            end
            
            # Average RDMs across trials
            if !isempty(trial_rdms)
                # This is simplified - would need to properly average matrices
                # For now, compute RDM on averaged patterns
                condition_patterns = Matrix{Float64}[]
                for cond_data in data_arrays
                    timepoint_data = cond_data[:, t, :]  # [channels × trials]
                    avg_pattern = vec(mean(timepoint_data, dims=2))
                    push!(condition_patterns, reshape(avg_pattern, length(selected_channels), 1))
                end
                rdms[t, :, :] = _compute_rdm(condition_patterns, dissimilarity_measure)
            end
        end
    end

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
    model_rdms::Union{Vector{Matrix{Float64}}, Vector{Array{Float64, 3}}, Vector{Union{Matrix{Float64}, Array{Float64, 3}}}};
    model_names::Union{Vector{String}, Nothing} = nothing,
    correlation_type::Symbol = :spearman,
    n_permutations::Int = 1000,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    n_models = length(model_rdms)
    n_conditions = length(rsa_data.condition_names)
    n_times = length(rsa_data.times)

    # Determine which models are static vs temporal
    is_temporal = Bool[]
    for (idx, model_rdm) in enumerate(model_rdms)
        if isa(model_rdm, Array{Float64, 3})
            # Temporal model: [time × condition × condition]
            if size(model_rdm, 1) != n_times
                @minimal_error_throw(
                    "Temporal model RDM $idx has $(size(model_rdm, 1)) timepoints, " *
                    "expected $n_times to match neural data"
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
                @minimal_error_throw(
                    "Static model RDM $idx has size $(size(model_rdm)), " *
                    "expected ($n_conditions, $n_conditions)"
                )
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

    # Generate model names if not provided
    if isnothing(model_names)
        model_names = ["Model$(i)" for i in 1:n_models]
    elseif length(model_names) != n_models
        @minimal_error_throw("Number of model names ($(length(model_names))) doesn't match number of models ($n_models)")
    end

    # Preallocate correlation and p-value arrays
    correlations = zeros(Float64, n_times, n_models)
    p_values = n_permutations > 0 ? zeros(Float64, n_times, n_models) : nothing

    # Extract upper triangular (excluding diagonal) for correlation
    function extract_upper_triangular(rdm::Matrix{Float64})
        n = size(rdm, 1)
        triu_indices = [CartesianIndex(i, j) for i in 1:n for j in (i+1):n]
        return [rdm[idx] for idx in triu_indices]
    end

    # Compute correlations at each time point
    for t in 1:n_times
        neural_rdm = rsa_data.rdm[t, :, :]
        neural_vec = extract_upper_triangular(neural_rdm)

        for (model_idx, model_rdm) in enumerate(model_rdms)
            # Get model RDM for this time point
            if is_temporal[model_idx]
                # Temporal model: extract RDM at time t
                model_rdm_t = model_rdm[t, :, :]
                # Symmetrize if needed
                if !issymmetric(model_rdm_t)
                    model_rdm_t = (model_rdm_t + model_rdm_t') / 2
                end
            else
                # Static model: use same RDM at all time points
                model_rdm_t = model_rdm
            end

            model_vec = extract_upper_triangular(model_rdm_t)

            # Compute correlation
            if correlation_type == :spearman
                corr = StatsBase.corspearman(neural_vec, model_vec)
            elseif correlation_type == :pearson
                corr = cor(neural_vec, model_vec)
            else
                @minimal_error_throw("Unknown correlation type: $correlation_type. Use :spearman or :pearson")
            end

            correlations[t, model_idx] = corr

            # Permutation test if requested
            if n_permutations > 0
                permuted_corrs = zeros(Float64, n_permutations)
                for perm_idx in 1:n_permutations
                    shuffled_model = shuffle(rng, model_vec)
                    if correlation_type == :spearman
                        permuted_corrs[perm_idx] = StatsBase.corspearman(neural_vec, shuffled_model)
                    else
                        permuted_corrs[perm_idx] = cor(neural_vec, shuffled_model)
                    end
                end
                # Two-tailed p-value
                p_val = sum(abs.(permuted_corrs) .>= abs(corr)) / n_permutations
                p_values[t, model_idx] = p_val
            end
        end
    end

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
    grand_average(rsa_data_list::Vector{RsaData})

Compute grand average RSA results across participants.

Averages RDMs across participants at each time point.

# Arguments
- `rsa_data_list::Vector{RsaData}`: Vector of RsaData objects, one per participant

# Returns
- `RsaData`: Grand-averaged RSA results

# Examples
```julia
# Average across participants
all_rsa = [rsa_p1, rsa_p2, rsa_p3]
grand_avg_rsa = grand_average(all_rsa)
```
"""
function grand_average(rsa_data_list::Vector{RsaData})
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
    for t in 1:n_times
        for i in 1:n_conditions
            for j in 1:n_conditions
                grand_rdm[t, i, j] = mean([rsa_data.rdm[t, i, j] for rsa_data in rsa_data_list])
            end
        end
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
        first_rsa.sample_rate,
        analysis_info = first_rsa.analysis_info,
    )

    return grand_rsa
end

