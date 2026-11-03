# Cross-validated RDM computation for Representational Similarity Analysis (RSA).
#
# Implements unbiased cross-fold RDM computation (Crossnobis, cross-validated squared
# Euclidean, and cross-fold correlation distances; Walther et al., 2016).

# === CROSS-VALIDATED RDM COMPUTATION ===

"""
    rsa_crossvalidated(
        epochs::AbstractVector{<:EpochData};
        channel_selection = channels(),
        sample_selection = samples(),
        interval_selection::Interval = times(),
        dissimilarity_measure::Symbol = :correlation,
        cv_method::Symbol = :splithalf,
        n_folds::Int = 5,
        n_iterations::Int = 100,
        normalize_method::Symbol = :none,
        rng::AbstractRNG = Random.default_rng(),
    )

Compute unbiased cross-validated RDMs across disjoint trial splits (Walther et al., 2016).

For `:euclidean`, evaluates cross-validated squared Euclidean distance
`(1 / P) * (x̄_i^(A) - x̄_j^(A))' * (x̄_i^(B) - x̄_j^(B))`, which has zero expectation under `H_0`.
For `:mahalanobis` (or `:crossnobis`), evaluates the symmetrized Crossnobis distance using
within-condition demeaned Ledoit–Wolf noise covariance estimated strictly within each fold.

# Arguments
- `epochs::AbstractVector{<:EpochData}`: Vector of `EpochData`, one per condition
- `channel_selection`: Channel selection predicate (default: all channels)
- `sample_selection`: Sample selection predicate (default: all time points)
- `interval_selection::Interval`: Time interval filter
- `dissimilarity_measure::Symbol`: `:correlation` (`:pearson`), `:spearman`, `:euclidean`, or `:mahalanobis` (`:crossnobis`)
- `cv_method::Symbol`: `:splithalf` (default), `:leaveoneout`, or `:kfold`
- `n_folds::Int`: Number of folds for `:kfold` (default: `5`)
- `n_iterations::Int`: Number of random split-half iterations for `:splithalf` (default: `100`)
- `normalize_method::Symbol`: RDM normalization (`:none`, `:zscore`, `:rank`, `:minmax`)
- `rng::AbstractRNG`: Random number generator for reproducible splits

# Returns
- `RsaData`: Cross-validated RSA results

# Examples
```julia
rsa_cv = EegFun.rsa_crossvalidated(epochs; cv_method = :splithalf, n_iterations = 100)
rsa_cn = EegFun.rsa_crossvalidated(epochs; dissimilarity_measure = :crossnobis, cv_method = :kfold, n_folds = 5)
```
"""
function rsa_crossvalidated(
    epochs::AbstractVector{<:EpochData};
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    dissimilarity_measure::Symbol = :correlation,
    cv_method::Symbol = :splithalf,
    n_folds::Int = 5,
    n_iterations::Int = 100,
    normalize_method::Symbol = :none,
    rng::AbstractRNG = Random.default_rng(),
)
    if !(dissimilarity_measure in (:correlation, :pearson, :spearman, :euclidean, :mahalanobis, :crossnobis))
        @minimal_error(
            "Unknown dissimilarity measure: $dissimilarity_measure. " *
            "Use :correlation, :spearman, :euclidean, :mahalanobis, or :crossnobis"
        )
    end

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    data_arrays, times, n_trials_per_condition, condition_names, selected_channels, first_epoch =
        _prepare_and_validate_rsa(epochs, channel_selection, combined_sel)

    n_conditions = length(condition_names)
    n_timepoints = length(times)
    min_trials = minimum(n_trials_per_condition)

    if min_trials < 2
        @minimal_error("Cross-validated RSA requires at least 2 trials per condition, got $min_trials")
    end

    if cv_method == :splithalf
        if n_iterations < 1
            @minimal_error("n_iterations must be >= 1, got $n_iterations")
        end
        if min_trials < 10
            @minimal_warning "Split-half CV works best with 20+ trials per condition, got $min_trials"
        end
        rdms = _cv_splithalf(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure, n_iterations; rng = rng)
    elseif cv_method == :leaveoneout
        if min_trials < 5
            @minimal_warning "Leave-one-out CV works best with at least 5 trials per condition, got $min_trials"
        end
        rdms = _cv_leaveoneout(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
    elseif cv_method == :kfold
        if n_folds < 2
            @minimal_error("K-fold CV requires n_folds >= 2, got $n_folds")
        end
        if min_trials < n_folds
            @minimal_error("K-fold CV requires at least $n_folds trials per condition, got $min_trials")
        end
        rdms = _cv_kfold(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure, n_folds; rng = rng)
    else
        @minimal_error("Unknown CV method: $cv_method. Use :splithalf, :leaveoneout, or :kfold")
    end

    _normalize_rdms!(rdms, normalize_method)

    return RsaData(
        first_epoch.file,
        condition_names,
        times,
        rdms,
        dissimilarity_measure,
        selected_channels,
        first_epoch.layout,
        first_epoch.sample_rate;
        analysis_info = first_epoch.analysis_info,
    )
end

# === CROSS-FOLD RDM EVALUATION ===

function _balanced_fold_ranges(n_trials::Int, n_folds::Int)
    base, rem = fldmod(n_trials, n_folds)
    ranges = Vector{UnitRange{Int}}(undef, n_folds)
    start_idx = 1
    for k = 1:n_folds
        len = base + (k <= rem ? 1 : 0)
        stop_idx = start_idx + len - 1
        ranges[k] = start_idx:stop_idx
        start_idx = stop_idx + 1
    end
    return ranges
end

function _compute_subset_means!(
    means_out::Vector{Matrix{Float64}},
    data_arrays::Vector{Array{Float64,3}},
    trial_indices::Vector{Vector{Int}},
)
    n_channels = size(data_arrays[1], 1)
    n_timepoints = size(data_arrays[1], 2)
    for (c, cond_data) in enumerate(data_arrays)
        idxs = trial_indices[c]
        cm = means_out[c]
        fill!(cm, 0.0)
        n_tr = length(idxs)
        n_tr == 0 && continue
        @inbounds for tr in idxs
            for t = 1:n_timepoints
                @simd for ch = 1:n_channels
                    cm[ch, t] += cond_data[ch, t, tr]
                end
            end
        end
        cm ./= n_tr
    end
    return means_out
end

function _accumulate_crossfold_rdms!(
    rdms_sum::Array{Float64,3},
    data_arrays::Vector{Array{Float64,3}},
    foldA_indices::Vector{Vector{Int}},
    foldB_indices::Vector{Vector{Int}},
    meansA::Vector{Matrix{Float64}},
    meansB::Vector{Matrix{Float64}},
    dissimilarity_measure::Symbol,
)
    n_channels = size(data_arrays[1], 1)
    n_timepoints = size(data_arrays[1], 2)
    n_conditions = length(data_arrays)

    _compute_subset_means!(meansA, data_arrays, foldA_indices)
    _compute_subset_means!(meansB, data_arrays, foldB_indices)

    diffA = Vector{Float64}(undef, n_channels)
    diffB = Vector{Float64}(undef, n_channels)
    use_mahal = dissimilarity_measure in (:mahalanobis, :crossnobis)
    foldB_has_cov = use_mahal && any(idxs -> length(idxs) >= 2, foldB_indices)
    inv_P = 1.0 / n_channels

    for t = 1:n_timepoints
        inv_cov = if use_mahal
            inv_cov_A = _invert_spd_covariance(_compute_pooled_covariance_subset(data_arrays, foldA_indices, t))
            if foldB_has_cov
                inv_cov_B = _invert_spd_covariance(_compute_pooled_covariance_subset(data_arrays, foldB_indices, t))
                0.5 .* (inv_cov_A .+ inv_cov_B)
            else
                inv_cov_A
            end
        else
            nothing
        end

        for i = 1:n_conditions
            mAi = meansA[i]
            mBi = meansB[i]
            for j = (i+1):n_conditions
                mAj = meansA[j]
                mBj = meansB[j]

                d_ij = if dissimilarity_measure == :euclidean
                    ip = 0.0
                    @inbounds @simd for ch = 1:n_channels
                        ip += (mAi[ch, t] - mAj[ch, t]) * (mBi[ch, t] - mBj[ch, t])
                    end
                    ip * inv_P
                elseif use_mahal
                    @inbounds @simd for ch = 1:n_channels
                        diffA[ch] = mAi[ch, t] - mAj[ch, t]
                        diffB[ch] = mBi[ch, t] - mBj[ch, t]
                    end
                    dot(diffA, inv_cov::Matrix{Float64}, diffB)
                else
                    pAi = @view mAi[:, t]
                    pAj = @view mAj[:, t]
                    pBi = @view mBi[:, t]
                    pBj = @view mBj[:, t]
                    0.5 *
                    (_compute_dissimilarity(pAi, pBj, dissimilarity_measure) + _compute_dissimilarity(pBi, pAj, dissimilarity_measure))
                end

                rdms_sum[t, i, j] += d_ij
                rdms_sum[t, j, i] += d_ij
            end
        end
    end

    return rdms_sum
end

# === CROSS-VALIDATION METHODS ===

function _cv_splithalf(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
    n_iterations::Int;
    rng::AbstractRNG = Random.default_rng(),
)
    n_channels = length(selected_channels)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    foldA = Vector{Vector{Int}}(undef, n_conditions)
    foldB = Vector{Vector{Int}}(undef, n_conditions)
    meansA = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]
    meansB = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]

    for _ = 1:n_iterations
        for (c, cond_data) in enumerate(data_arrays)
            n_trials = size(cond_data, 3)
            n_half = div(n_trials, 2)
            perm = randperm(rng, n_trials)
            foldA[c] = perm[1:n_half]
            foldB[c] = perm[(n_half+1):n_trials]
        end
        _accumulate_crossfold_rdms!(rdms_sum, data_arrays, foldA, foldB, meansA, meansB, dissimilarity_measure)
    end

    return rdms_sum ./ n_iterations
end

function _cv_leaveoneout(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
)
    n_channels = length(selected_channels)
    n_loo = minimum(size(d, 3) for d in data_arrays)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    foldA = Vector{Vector{Int}}(undef, n_conditions)
    foldB = Vector{Vector{Int}}(undef, n_conditions)
    meansA = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]
    meansB = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]

    cond_ranges = [_balanced_fold_ranges(size(d, 3), n_loo) for d in data_arrays]

    for fold = 1:n_loo
        for (c, cond_data) in enumerate(data_arrays)
            n_trials = size(cond_data, 3)
            r = cond_ranges[c][fold]
            foldB[c] = collect(r)
            foldA[c] = [i for i = 1:n_trials if !(i in r)]
        end
        _accumulate_crossfold_rdms!(rdms_sum, data_arrays, foldA, foldB, meansA, meansB, dissimilarity_measure)
    end

    return rdms_sum ./ n_loo
end

function _cv_kfold(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
    n_folds::Int;
    rng::AbstractRNG = Random.default_rng(),
)
    n_channels = length(selected_channels)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    foldA = Vector{Vector{Int}}(undef, n_conditions)
    foldB = Vector{Vector{Int}}(undef, n_conditions)
    meansA = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]
    meansB = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]

    cond_perms = [randperm(rng, size(d, 3)) for d in data_arrays]
    cond_ranges = [_balanced_fold_ranges(size(d, 3), n_folds) for d in data_arrays]

    for fold = 1:n_folds
        for (c, cond_data) in enumerate(data_arrays)
            n_trials = size(cond_data, 3)
            perm = cond_perms[c]
            r = cond_ranges[c][fold]
            foldB[c] = perm[r]
            foldA[c] = [perm[i] for i = 1:n_trials if !(i in r)]
        end
        _accumulate_crossfold_rdms!(rdms_sum, data_arrays, foldA, foldB, meansA, meansB, dissimilarity_measure)
    end

    return rdms_sum ./ n_folds
end
