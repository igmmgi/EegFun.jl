"""
Cross-validated RDM computation for Representational Similarity Analysis (RSA).

This module implements cross-validated RDM computation methods to assess
reliability and robustness of representational dissimilarity estimates.
"""

# ==============================================================================
#   CROSS-VALIDATED RDM COMPUTATION
# ==============================================================================

"""
    rsa_crossvalidated(
        epochs::Vector{EpochData};
        channel_selection::Function = channels(),
        sample_selection::Function = samples(),
        dissimilarity_measure::Symbol = :correlation,
        cv_method::Symbol = :splithalf,
        n_folds::Int = 5,
        n_iterations::Int = 100,
        normalize_method::Symbol = :none,
        rng::AbstractRNG = Random.GLOBAL_RNG,
    )

Compute cross-validated RDMs to assess reliability of representational structure.

Cross-validation provides a robust estimate of RDM structure by computing RDMs
on independent subsets of trials and averaging them. This reduces the impact of
noise and outliers compared to computing a single RDM on all trials.

# Cross-Validation Methods

## Split-Half (`:splithalf`)
- Randomly split trials into two halves
- Compute RDM on each half independently
- Average the two RDMs
- Repeat for `n_iterations` and average results
- **Use when**: You have many trials (20+) and want robust estimates

## Leave-One-Out (`:leaveoneout`)
- For each trial, compute RDM using all other trials
- Average all leave-one-out RDMs
- **Use when**: You have few trials (10-20) and want maximum data usage
- **Note**: Computationally expensive for many trials

## K-Fold (`:kfold`)
- Split trials into `n_folds` groups
- For each fold, compute RDM using all other folds
- Average all fold RDMs
- **Use when**: You want a balance between split-half and leave-one-out

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData, one per condition
- `channel_selection::Function`: Channel selection predicate (default: all channels)
- `sample_selection::Function`: Sample selection predicate (default: all time points)
- `dissimilarity_measure::Symbol`: Dissimilarity measure (:correlation, :spearman, :euclidean)
- `cv_method::Symbol`: Cross-validation method (:splithalf, :leaveoneout, :kfold)
- `n_folds::Int`: Number of folds for k-fold CV (default: 5)
- `n_iterations::Int`: Number of iterations for split-half (default: 100)
- `normalize_method::Symbol`: RDM normalization method (default: :none)
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
- `RsaData`: Cross-validated RSA results

# Examples
```julia
# Split-half cross-validation (default)
rsa_cv = rsa_crossvalidated(epochs, cv_method=:splithalf, n_iterations=100)

# Leave-one-out cross-validation
rsa_cv = rsa_crossvalidated(epochs, cv_method=:leaveoneout)

# K-fold cross-validation
rsa_cv = rsa_crossvalidated(epochs, cv_method=:kfold, n_folds=5)

# With normalization
rsa_cv = rsa_crossvalidated(epochs, cv_method=:splithalf, normalize_method=:zscore)
```

# Notes
- Cross-validation requires sufficient trials per condition (minimum 10, preferably 20+)
- Split-half is fastest and recommended for most cases
- Leave-one-out is most thorough but computationally expensive
- Results are more reliable than single RDM computation
"""
function rsa_crossvalidated(
    epochs::Vector{EpochData};
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    dissimilarity_measure::Symbol = :correlation,
    cv_method::Symbol = :splithalf,
    n_folds::Int = 5,
    n_iterations::Int = 100,
    normalize_method::Symbol = :none,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    # Input validations
    if isempty(epochs)
        @minimal_error_throw("Cannot perform cross-validated RSA with empty epochs vector")
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

    # Check minimum trials
    min_trials = minimum(n_trials_per_condition)
    if cv_method == :splithalf && min_trials < 10
        @minimal_warning "Split-half CV works best with 20+ trials per condition, got $min_trials"
    elseif cv_method == :leaveoneout && min_trials < 5
        @minimal_warning "Leave-one-out CV requires at least 5 trials per condition, got $min_trials"
    elseif cv_method == :kfold && min_trials < n_folds
        @minimal_error_throw("K-fold CV requires at least $n_folds trials per condition, got $min_trials")
    end

    # Get metadata
    first_epoch = epochs[1]
    condition_names = [e.condition_name for e in epochs]
    selected_channels = channel_labels(epochs)

    # Preallocate RDM array
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    # Compute cross-validated RDMs based on method
    if cv_method == :splithalf
        rdms = _cv_splithalf(
            data_arrays,
            n_timepoints,
            n_conditions,
            selected_channels,
            dissimilarity_measure,
            n_iterations,
            rng,
        )
    elseif cv_method == :leaveoneout
        rdms = _cv_leaveoneout(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
    elseif cv_method == :kfold
        rdms =
            _cv_kfold(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure, n_folds, rng)
    else
        @minimal_error_throw("Unknown CV method: $cv_method. Use :splithalf, :leaveoneout, or :kfold")
    end

    # Apply normalization if requested
    if normalize_method != :none
        for t = 1:n_timepoints
            rdms[t, :, :] = normalize_rdm(rdms[t, :, :]; method = normalize_method)
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
#   CROSS-VALIDATION METHODS
# ==============================================================================

"""
    _cv_splithalf(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure, n_iterations, rng)

Compute split-half cross-validated RDMs.
"""
function _cv_splithalf(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
    n_iterations::Int,
    rng::AbstractRNG,
)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    for iter = 1:n_iterations
        # For each condition, randomly split trials into two halves
        half1_data = Vector{Array{Float64,3}}()
        half2_data = Vector{Array{Float64,3}}()

        for cond_data in data_arrays
            n_trials = size(cond_data, 3)
            n_half = div(n_trials, 2)

            # Random permutation
            perm = randperm(rng, n_trials)
            half1_indices = perm[1:n_half]
            half2_indices = perm[(n_half+1):end]

            push!(half1_data, cond_data[:, :, half1_indices])
            push!(half2_data, cond_data[:, :, half2_indices])
        end

        # Compute RDMs for each half
        rdm_half1 =
            _compute_rdms_from_data(half1_data, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
        rdm_half2 =
            _compute_rdms_from_data(half2_data, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)

        # Average the two halves
        rdms_sum .+= (rdm_half1 .+ rdm_half2) ./ 2
    end

    # Average across iterations
    return rdms_sum ./ n_iterations
end

"""
    _cv_leaveoneout(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)

Compute leave-one-out cross-validated RDMs.
"""
function _cv_leaveoneout(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
)
    min_trials = minimum([size(d, 3) for d in data_arrays])
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    for trial_idx = 1:min_trials
        # Create data arrays excluding this trial
        loo_data = Vector{Array{Float64,3}}()

        for cond_data in data_arrays
            n_trials = size(cond_data, 3)
            if trial_idx <= n_trials
                # Exclude trial_idx
                indices = [i for i = 1:n_trials if i != trial_idx]
                push!(loo_data, cond_data[:, :, indices])
            else
                # Use all trials if this condition has fewer trials
                push!(loo_data, cond_data)
            end
        end

        # Compute RDM without this trial
        rdm_loo =
            _compute_rdms_from_data(loo_data, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
        rdms_sum .+= rdm_loo
    end

    # Average across all leave-one-out iterations
    return rdms_sum ./ min_trials
end

"""
    _cv_kfold(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure, n_folds, rng)

Compute k-fold cross-validated RDMs.
"""
function _cv_kfold(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
    n_folds::Int,
    rng::AbstractRNG,
)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    for fold = 1:n_folds
        # For each condition, select trials for this fold
        fold_data = Vector{Array{Float64,3}}()

        for cond_data in data_arrays
            n_trials = size(cond_data, 3)

            # Random permutation
            perm = randperm(rng, n_trials)

            # Select trials NOT in this fold
            fold_size = div(n_trials, n_folds)
            fold_start = (fold - 1) * fold_size + 1
            fold_end = fold == n_folds ? n_trials : fold * fold_size

            # Use all trials except those in current fold
            indices = [perm[i] for i = 1:n_trials if i < fold_start || i > fold_end]

            push!(fold_data, cond_data[:, :, indices])
        end

        # Compute RDM for this fold
        rdm_fold =
            _compute_rdms_from_data(fold_data, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
        rdms_sum .+= rdm_fold
    end

    # Average across folds
    return rdms_sum ./ n_folds
end

"""
    _compute_rdms_from_data(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)

Helper function to compute RDMs from data arrays.
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
