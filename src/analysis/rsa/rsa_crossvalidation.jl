# Cross-validated unbiased RDM computation for RSA (Walther et al., 2016).

struct _CrossFoldWorkspace
    meansA::Vector{Matrix{Float64}}
    meansB::Vector{Matrix{Float64}}
    diffA::Vector{Float64}
    diffB::Vector{Float64}
    cov_ws::Union{Nothing,_PooledCovWorkspace}
    inv_cov_A::Union{Nothing,Matrix{Float64}}
    inv_cov_B::Union{Nothing,Matrix{Float64}}
    unitA::Union{Nothing,Matrix{Float64}}
    unitB::Union{Nothing,Matrix{Float64}}
    varA::Vector{Bool}
    varB::Vector{Bool}
end

function _CrossFoldWorkspace(n_channels::Int, n_timepoints::Int, n_conditions::Int, max_total_trials::Int, dissimilarity_measure::Symbol)
    use_mahal = dissimilarity_measure in (:mahalanobis, :crossnobis)
    use_pearson = dissimilarity_measure in (:correlation, :pearson)
    return _CrossFoldWorkspace(
        [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions],
        [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions],
        Vector{Float64}(undef, n_channels),
        Vector{Float64}(undef, n_channels),
        use_mahal ? _PooledCovWorkspace(n_channels, max_total_trials) : nothing,
        use_mahal ? Matrix{Float64}(undef, n_channels, n_channels) : nothing,
        use_mahal ? Matrix{Float64}(undef, n_channels, n_channels) : nothing,
        use_pearson ? Matrix{Float64}(undef, n_channels, n_conditions) : nothing,
        use_pearson ? Matrix{Float64}(undef, n_channels, n_conditions) : nothing,
        Vector{Bool}(undef, n_conditions),
        Vector{Bool}(undef, n_conditions),
    )
end

function _balanced_fold_ranges(n_trials::Int, n_folds::Int)
    base, rem = fldmod(n_trials, n_folds)
    ranges = Vector{UnitRange{Int}}(undef, n_folds)
    start_idx = 1
    for k = 1:n_folds
        stop_idx = start_idx + base + (k <= rem ? 1 : 0) - 1
        ranges[k] = start_idx:stop_idx
        start_idx = stop_idx + 1
    end
    return ranges
end

function _fill_complement_range!(foldA::Vector{Int}, foldB::Vector{Int}, perm::Union{Nothing,Vector{Int}}, n_trials::Int, r::UnitRange{Int})
    n_b = length(r)
    resize!(foldB, n_b)
    resize!(foldA, n_trials - n_b)
    if isnothing(perm)
        @inbounds for (k, idx) in enumerate(r)
            foldB[k] = idx
        end
        pos = 1
        @inbounds for i = 1:(first(r)-1)
            foldA[pos] = i
            pos += 1
        end
        @inbounds for i = (last(r)+1):n_trials
            foldA[pos] = i
            pos += 1
        end
    else
        @inbounds for (k, idx) in enumerate(r)
            foldB[k] = perm[idx]
        end
        pos = 1
        @inbounds for i = 1:(first(r)-1)
            foldA[pos] = perm[i]
            pos += 1
        end
        @inbounds for i = (last(r)+1):n_trials
            foldA[pos] = perm[i]
            pos += 1
        end
    end
    return nothing
end

function _compute_subset_means!(
    means_out::Vector{Matrix{Float64}},
    data_arrays::Vector{Array{Float64,3}},
    trial_indices::Vector{Vector{Int}},
)
    n_channels, n_timepoints, _ = size(data_arrays[1])
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
    ws::_CrossFoldWorkspace,
    data_arrays::Vector{Array{Float64,3}},
    foldA_indices::Vector{Vector{Int}},
    foldB_indices::Vector{Vector{Int}},
    dissimilarity_measure::Symbol,
)
    n_channels, n_timepoints, _ = size(data_arrays[1])
    n_conditions = length(data_arrays)

    _compute_subset_means!(ws.meansA, data_arrays, foldA_indices)
    _compute_subset_means!(ws.meansB, data_arrays, foldB_indices)

    meansA, meansB = ws.meansA, ws.meansB
    diffA, diffB = ws.diffA, ws.diffB
    use_mahal = dissimilarity_measure in (:mahalanobis, :crossnobis)
    use_pearson = dissimilarity_measure in (:correlation, :pearson)
    foldB_has_cov = use_mahal && any(idxs -> length(idxs) >= 2, foldB_indices)
    inv_P = 1.0 / n_channels

    for t = 1:n_timepoints
        if use_mahal
            inv_cov = ws.inv_cov_A::Matrix{Float64}
            cov_ws = ws.cov_ws::_PooledCovWorkspace
            _compute_inv_pooled_covariance!(inv_cov, cov_ws, data_arrays, foldA_indices, t)
            if foldB_has_cov
                inv_cov_B = ws.inv_cov_B::Matrix{Float64}
                _compute_inv_pooled_covariance!(inv_cov_B, cov_ws, data_arrays, foldB_indices, t)
                @. inv_cov = 0.5 * (inv_cov + inv_cov_B)
            end
        elseif use_pearson
            uA = ws.unitA::Matrix{Float64}
            uB = ws.unitB::Matrix{Float64}
            for c = 1:n_conditions
                ws.varA[c] = _zscore_unit_norm!(@view(uA[:, c]), @view(meansA[c][:, t]))
                ws.varB[c] = _zscore_unit_norm!(@view(uB[:, c]), @view(meansB[c][:, t]))
            end
        end

        @inbounds for i = 1:n_conditions
            mAi = meansA[i]
            mBi = meansB[i]
            for j = (i+1):n_conditions
                mAj = meansA[j]
                mBj = meansB[j]

                d_ij = if dissimilarity_measure == :euclidean
                    ip = 0.0
                    @simd for ch = 1:n_channels
                        ip += (mAi[ch, t] - mAj[ch, t]) * (mBi[ch, t] - mBj[ch, t])
                    end
                    ip * inv_P
                elseif use_mahal
                    @simd for ch = 1:n_channels
                        diffA[ch] = mAi[ch, t] - mAj[ch, t]
                        diffB[ch] = mBi[ch, t] - mBj[ch, t]
                    end
                    dot(diffA, ws.inv_cov_A::Matrix{Float64}, diffB)
                elseif use_pearson
                    uA = ws.unitA::Matrix{Float64}
                    uB = ws.unitB::Matrix{Float64}
                    cAB = (ws.varA[i] && ws.varB[j]) ? dot(@view(uA[:, i]), @view(uB[:, j])) : NaN
                    cBA = (ws.varB[i] && ws.varA[j]) ? dot(@view(uB[:, i]), @view(uA[:, j])) : NaN
                    1.0 - 0.5 * (cAB + cBA)
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

function _cv_splithalf(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    n_channels::Int,
    dissimilarity_measure::Symbol,
    n_iterations::Int;
    rng::AbstractRNG = Random.default_rng(),
)
    max_total = sum(size(d, 3) for d in data_arrays)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    ws = _CrossFoldWorkspace(n_channels, n_timepoints, n_conditions, max_total, dissimilarity_measure)
    perms = [Vector{Int}(undef, size(d, 3)) for d in data_arrays]
    foldA = [Vector{Int}(undef, div(size(d, 3), 2)) for d in data_arrays]
    foldB = [Vector{Int}(undef, size(d, 3) - div(size(d, 3), 2)) for d in data_arrays]

    for _ = 1:n_iterations
        for (c, cond_data) in enumerate(data_arrays)
            n_trials = size(cond_data, 3)
            n_half = div(n_trials, 2)
            perm = Random.randperm!(rng, perms[c])
            copyto!(foldA[c], 1, perm, 1, n_half)
            copyto!(foldB[c], 1, perm, n_half + 1, n_trials - n_half)
        end
        _accumulate_crossfold_rdms!(rdms_sum, ws, data_arrays, foldA, foldB, dissimilarity_measure)
    end

    rdms_sum ./= n_iterations
    return rdms_sum
end

function _cv_folded(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    n_channels::Int,
    dissimilarity_measure::Symbol,
    n_folds::Int;
    shuffle_trials::Bool,
    rng::AbstractRNG = Random.default_rng(),
)
    max_total = sum(size(d, 3) for d in data_arrays)
    rdms_sum = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    ws = _CrossFoldWorkspace(n_channels, n_timepoints, n_conditions, max_total, dissimilarity_measure)
    foldA = [sizehint!(Int[], size(d, 3)) for d in data_arrays]
    foldB = [sizehint!(Int[], size(d, 3)) for d in data_arrays]
    cond_perms = shuffle_trials ? [Random.randperm(rng, size(d, 3)) for d in data_arrays] : nothing
    cond_ranges = [_balanced_fold_ranges(size(d, 3), n_folds) for d in data_arrays]

    for fold = 1:n_folds
        for (c, cond_data) in enumerate(data_arrays)
            perm = isnothing(cond_perms) ? nothing : cond_perms[c]
            _fill_complement_range!(foldA[c], foldB[c], perm, size(cond_data, 3), cond_ranges[c][fold])
        end
        _accumulate_crossfold_rdms!(rdms_sum, ws, data_arrays, foldA, foldB, dissimilarity_measure)
    end

    rdms_sum ./= n_folds
    return rdms_sum
end

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
    ) -> RsaData

Compute unbiased cross-validated RDMs across disjoint trial splits (Walther et al., 2016),
supporting `:splithalf`, `:leaveoneout`, and `:kfold` cross-validation and `:correlation`
(`:pearson`), `:spearman`, `:euclidean`, and `:mahalanobis` (`:crossnobis`) measures.
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
    dissimilarity_measure in (:correlation, :pearson, :spearman, :euclidean, :mahalanobis, :crossnobis) || @minimal_error(
        "Unknown dissimilarity measure: $dissimilarity_measure. Use :correlation, :spearman, :euclidean, :mahalanobis, or :crossnobis"
    )

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    data_arrays, times, n_trials_per_condition, condition_names, selected_channels, layout_sub, first_epoch =
        _prepare_and_validate_rsa(epochs, channel_selection, combined_sel)

    n_conditions = length(condition_names)
    n_timepoints = length(times)
    n_ch = length(selected_channels)
    min_trials = minimum(n_trials_per_condition)
    min_trials < 2 && @minimal_error("Cross-validated RSA requires at least 2 trials per condition, got $min_trials")

    rdms = if cv_method == :splithalf
        n_iterations < 1 && @minimal_error("n_iterations must be >= 1, got $n_iterations")
        min_trials < 10 && @minimal_warning "Split-half CV works best with 20+ trials per condition, got $min_trials"
        _cv_splithalf(data_arrays, n_timepoints, n_conditions, n_ch, dissimilarity_measure, n_iterations; rng = rng)
    elseif cv_method == :leaveoneout
        min_trials < 5 && @minimal_warning "Leave-one-out CV works best with at least 5 trials per condition, got $min_trials"
        _cv_folded(data_arrays, n_timepoints, n_conditions, n_ch, dissimilarity_measure, min_trials; shuffle_trials = false)
    elseif cv_method == :kfold
        n_folds < 2 && @minimal_error("K-fold CV requires n_folds >= 2, got $n_folds")
        min_trials < n_folds && @minimal_error("K-fold CV requires at least $n_folds trials per condition, got $min_trials")
        _cv_folded(data_arrays, n_timepoints, n_conditions, n_ch, dissimilarity_measure, n_folds; shuffle_trials = true, rng = rng)
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
        layout_sub,
        first_epoch.sample_rate;
        analysis_info = copy(first_epoch.analysis_info),
    )
end
