"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; kwargs...)
    fit_mass_lmm(eeg_data::AbstractArray{<:Real,3}, meta_df::DataFrame, f::FormulaTerm; kwargs...)

Fit a linear mixed model at every channel × time point and (optionally) control the family-wise
error rate with a whole-block sign-flip permutation test.

# Inference
Each point is fitted once with MixedModels.jl (Wald z = β̂ / SE, uncorrected p in `p_parametric`).
For permutation inference, the Freedman–Lane residuals of the reduced model (without the tested
coefficient) are multiplied by a random sign per *block* (by default the outermost random-effects
grouping factor, e.g. subject). This is valid for the population null (the average effect across
subjects is zero) for both within- and between-subject predictors and unequal trial counts, assuming
independent blocks and symmetric residuals. Because flipping a whole block leaves its covariance
unchanged, the variance parameters estimated on the observed data are used for every flip, and the
null statistics are computed in closed form (no refitting). The permutation null retains all
between-subject variability, which also protects against an under-specified random-effects structure.

All random-effects grouping factors must be nested within the flipping block (e.g. `(1 | subject)`,
`(1 + condition | subject)`, `(1 | subject & session)`). Crossed designs (e.g. subject × item) have no
exact permutation scheme and raise an error when `n_perms > 0`; use `n_perms = 0` for Wald-only results.

# Keyword Arguments
- `n_perms::Int=1000`: Number of random sign flips (`0` = no permutation inference). If all
  `2^(n_blocks-1)` distinct patterns fit within `n_perms`, they are enumerated exactly instead.
- `tested_coefs=nothing`: Coefficients to permutation-test.
  - `nothing` or `:effects`: all non-intercept coefficients (default); `:all`: including the intercept.
  - `Symbol`, `String`, `Integer` or a `Vector` of these. An exact coefficient name selects that column;
    a term name (e.g. `:condition`) selects all of its contrast columns but not interactions.
- `test_intercept::Bool=false`: Also test the intercept.
- `flip_by=nothing`: Column defining the exchangeability blocks. Default: the outermost grouping factor.
- `correction::Symbol=:maxt`: FWER correction: `:maxt` (max-|z|), `:cluster` (spatio-temporal
  cluster mass) or `:tfce`.
- `cluster_threshold::Real=2.0`: |z| threshold for cluster formation (`:cluster`).
- `tfce_E=0.5`, `tfce_H=2.0`, `tfce_dh=0.1`: TFCE parameters (`:tfce`).
- `spatial_connectivity=nothing`: Channel adjacency for clustering/TFCE (built from the layout for EpochData).
- `REML::Bool=false`: Use REML instead of ML.
- `rng=Random.default_rng()`: Random number generator for the sign flips.
- Array method only: `channel_names`, `time_points`, `dim_order=(:channels, :timepoints, :epochs)`.

# Parallelism
Uses Distributed worker processes if any are available (load the packages on them first:
`@everywhere using EegFun, MixedModels`), otherwise Julia threads (`julia -t auto`).

# Rank-deficient designs
If the fixed-effects design matrix is rank deficient (e.g. a predictor that is constant in the
data, an empty design cell combined with an interaction, or exactly collinear predictors), the
aliased coefficient(s) cannot be estimated. A single warning names them; their `beta`, `se`,
`z_values`, `p_parametric` and `p_corrected` values are `NaN`, and they are excluded from
permutation testing. All estimable coefficients are fitted and tested as usual.

Returns an [`LmmStatsResult`](@ref).

**Requires `MixedModels.jl` and `StatsModels.jl` to be loaded.**
"""
function fit_mass_lmm(args...; kwargs...)
    @minimal_error "To use Mass-Univariate Mixed Models, you must first load the packages: `using MixedModels, StatsModels`"
end

# Permutation p-value of `obs` given the sorted null. With exhaustive enumeration the identity pattern is part
# of the null set (p = #{null ≥ obs} / N); with random flips the observed statistic is added (Phipson & Smyth, 2010).
function _lmm_perm_pvalue(null_sorted::AbstractVector, obs::Real, exhaustive::Bool)
    n = length(null_sorted)
    n_ge = n - searchsortedlast(null_sorted, obs * (1 - 1e-8))   # tolerance for floating-point ties
    return exhaustive ? n_ge / n : (n_ge + 1) / (n + 1)
end

function _update_clusters(clusters, stats, null_sorted, alpha, exhaustive, dims)
    updated = Cluster[]
    mask = zeros(Bool, dims)
    for (i, c) in enumerate(clusters)
        p_val = _lmm_perm_pvalue(null_sorted, abs(stats[i]), exhaustive)
        sig = p_val <= alpha
        push!(updated, Cluster(c.id, c.electrodes, c.time_indices, c.time_range, stats[i], p_val, sig, c.polarity, c.members))
        if sig
            for (e_idx, t_idx) in c.members
                mask[e_idx, t_idx] = true
            end
        end
    end
    return updated, mask
end

"""
    extract_predictor_stats(result::LmmStatsResult, coef_name; alpha = 0.05)

Format the permutation-corrected results for one predictor as a `PermutationResult`, so they can be
plotted with `plot_erp_stats`, `plot_topography_stats`, etc.

The correction method and its settings (cluster threshold, TFCE parameters, connectivity) are taken
from `result.inference`, i.e. exactly those used to build the null distribution in `fit_mass_lmm`.
Requires a result from the `EpochData` method of `fit_mass_lmm`.
"""
function extract_predictor_stats(result::LmmStatsResult, coef_name; alpha = 0.05)
    coef_str = String(coef_name)
    coef_idx = findfirst(==(coef_str), result.coefficients)
    isnothing(coef_idx) && @minimal_error "Coefficient '$coef_str' not found in model. Available: $(result.coefficients)"
    isnothing(result.p_corrected) && @minimal_error "No permutation inference in this result: run fit_mass_lmm with n_perms > 0."
    coef_idx in result.tested || @minimal_error "Coefficient '$coef_str' was not permutation-tested. Include it in `tested_coefs` " *
                                                "(or set `test_intercept = true`) when calling fit_mass_lmm."
    isnothing(result.epochs) && @minimal_error "extract_predictor_stats requires a result from the EpochData method of fit_mass_lmm."

    inf = result.inference
    n_electrodes = length(result.channels)
    n_time_points = length(result.time_points)
    z_map = result.z_values[:, :, coef_idx]
    null_dist = result.null_max[:, coef_idx]
    n_perms = length(null_dist)

    pos_clusters, neg_clusters = Cluster[], Cluster[]
    if inf.correction === :cluster
        spatial_connectivity = isnothing(inf.spatial_connectivity) ?
            _build_connectivity_matrix(result.channels, result.epochs.layout, :spatiotemporal) : inf.spatial_connectivity
        thr = inf.cluster_threshold
        raw_pos, raw_neg = _find_clusters(z_map .> thr, z_map .< -thr, result.channels, result.time_points,
            spatial_connectivity, :spatiotemporal)
        pos_stats = _compute_cluster_statistics(raw_pos, z_map; return_clusters = false)
        neg_stats = _compute_cluster_statistics(raw_neg, z_map; return_clusters = false)
        null_sorted = sort(null_dist)
        dims = (n_electrodes, n_time_points)
        pos_clusters, final_mask_pos = _update_clusters(raw_pos, pos_stats, null_sorted, alpha, inf.exhaustive, dims)
        neg_clusters, final_mask_neg = _update_clusters(raw_neg, neg_stats, null_sorted, alpha, inf.exhaustive, dims)
    else
        # :maxt and :tfce give point-wise FWER-corrected p-values
        sig = result.p_corrected[:, :, coef_idx] .<= alpha
        final_mask_pos = sig .& (z_map .> 0)
        final_mask_neg = sig .& (z_map .< 0)
    end

    se_diff = result.se[:, :, coef_idx]
    correction_symbol = inf.correction === :cluster ? :cluster_permutation : inf.correction
    test_info = TestInfo(:LMM, 0.0, Float64(alpha), :both, correction_symbol,
        ClusterInfo(:parametric, :spatiotemporal, n_perms))

    # ERP data: the average of all epochs as placeholder data for plotting
    erp_dummy = average_epochs(result.epochs)
    thr = inf.cluster_threshold
    return PermutationResult(
        test_info,
        [erp_dummy, erp_dummy], # Dummy data for Cond1/Cond2
        StatMatrix(z_map, nothing),
        Masks(final_mask_pos, final_mask_neg),
        Clusters(pos_clusters, neg_clusters),
        PermutationDistribution(null_dist, null_dist),
        result.channels,
        result.time_points,
        (thr, -thr),
        se_diff,
        zeros(n_electrodes, n_time_points),
        zeros(n_electrodes, n_time_points),
        se_diff
    )
end
