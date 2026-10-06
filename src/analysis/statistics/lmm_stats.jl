"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; kwargs...)
    fit_mass_lmm(eeg_data::AbstractArray, meta_df::DataFrame, f::FormulaTerm; kwargs...)

Fit a mass-univariate linear mixed model across all channels and timepoints.

# Keyword Arguments
- `n_perms::Int=1000`: Number of permutations for null hypothesis cluster testing.
- `tested_coefs=nothing`: Which fixed effect coefficients to permute for cluster null distributions.
  - `nothing` or `:effects`: Tests all non-intercept predictors (default).
  - `:all`: Tests all predictors including the intercept.
  - `Symbol`, `String`, or `Vector`: Test specific predictor(s) by name, substring, or index (e.g., `:condition`, `"condition"`, or `[:condition, :rt]`).
- `test_intercept::Bool=false`: If `true`, tests the intercept coefficient. Defaults to `false` because in ERP research, testing the baseline against zero is rarely meaningful, and skipping it yields an immediate ~2× speedup.

- `use_clusters::Bool=false`: Use spatio-temporal cluster-mass testing.
- `cluster_threshold::Float64=2.0`: T-statistic threshold for clustering.
- `permute_block`: Factor for independent shuffling within blocks (e.g., `:Subject`).
- `permute_crossed`: Tuple `(Subject, Item)` for synchronized crossed permutations.
- `perm_matrix`: Explicit precomputed permutation matrix (`n_epochs × n_perms`).

# Rank-deficient designs
If the fixed-effects design matrix is rank deficient (e.g. a predictor that is constant in the
data, an empty design cell combined with an interaction, or exactly collinear predictors), the
aliased coefficient(s) cannot be estimated. A single warning names them; their `beta`, `se`,
`t`, `p`, `p_uncorrected` and `p_corrected` values are `NaN`, and they are excluded from
permutation testing. All estimable coefficients are fitted and tested as usual.

**Requires `MixedModels.jl` and `StatsModels.jl` to be loaded.**
"""
function fit_mass_lmm(args...; kwargs...)
    error("To use Mass-Univariate Mixed Models, you must first load the packages: `using MixedModels, StatsModels`")
end

function _update_clusters(clusters, stats, null_dist, alpha, n_perms, dims)
    updated = Cluster[]
    mask = zeros(Bool, dims)
    for (i, c) in enumerate(clusters)
        stat_mag = abs(stats[i])
        p_val = (count(>=(stat_mag), null_dist) + 1) / (n_perms + 1)
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
    extract_predictor_stats(result::LmmStatsResult, coef_name; alpha=0.05, cluster_threshold=2.0, kwargs...)

Extracts the cluster-corrected statistics for a specific LMM predictor and formats them as a standard `PermutationResult`.
This allows the result to be plotted directly using `plot_erp_stats`, `plot_topography_stats`, etc.
"""
function extract_predictor_stats(result::LmmStatsResult, coef_name; 
    alpha = 0.05, 
    cluster_threshold = 2.0,
    use_tfce = false,
    tfce_E = 0.5,
    tfce_H = 2.0,
    tfce_dh = 0.1
)
    coef_str = String(coef_name)
    coef_idx = findfirst(==(coef_str), result.coefficients)
    if isnothing(coef_idx)
        error("Coefficient '$coef_str' not found in model. Available: $(result.coefficients)")
    end

    n_electrodes = length(result.channels)
    n_time_points = length(result.time_points)

    t_map = result.t_values[:, :, coef_idx]
    
    null_dist = result.max_cluster_mass_null[:, coef_idx]
    n_perms = length(null_dist)
    if n_perms > 0 && all(==(0), null_dist)
        error("Coefficient '$coef_str' was not permuted during fit_mass_lmm (e.g., intercept was skipped). To test this coefficient, run fit_mass_lmm with `test_intercept=true` or include it in `tested_coefs`.")
    end
    
    spatial_connectivity = EegFun._build_connectivity_matrix(result.channels, result.epochs.layout, :spatiotemporal)
    
    if use_tfce
        tfce_map = EegFun._compute_tfce(t_map, result.channels, Float64.(result.time_points), spatial_connectivity, :spatiotemporal; E=tfce_E, H=tfce_H, dh=tfce_dh)
        final_mask_pos = zeros(Bool, n_electrodes, n_time_points)
        final_mask_neg = zeros(Bool, n_electrodes, n_time_points)
        
        for t_idx in 1:n_time_points, e_idx in 1:n_electrodes
            val = tfce_map[e_idx, t_idx]
            iszero(val) && continue
            p_val = (count(>=(abs(val)), null_dist) + 1) / (n_perms + 1)
            p_val <= alpha || continue
            val > 0 ? (final_mask_pos[e_idx, t_idx] = true) : (final_mask_neg[e_idx, t_idx] = true)
        end
        updated_pos_clusters = Cluster[]
        updated_neg_clusters = Cluster[]
    else
        pos_clusters, neg_clusters = EegFun._find_clusters(
            t_map .> cluster_threshold, t_map .< -cluster_threshold,
            result.channels, result.time_points, spatial_connectivity, :spatiotemporal
        )
        
        electrode_to_idx = Dict(e => i for (i, e) in enumerate(result.channels))
        pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
        neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
        
        dims = (n_electrodes, n_time_points)
        updated_pos_clusters, final_mask_pos = _update_clusters(pos_clusters, pos_stats, null_dist, alpha, n_perms, dims)
        updated_neg_clusters, final_mask_neg = _update_clusters(neg_clusters, neg_stats, null_dist, alpha, n_perms, dims)
    end
    
    # Use stored standard errors directly from the LmmStatsResult
    se_diff = result.se[:, :, coef_idx]
    
    # Construct PermutationResult
    test_info = TestInfo(:LMM, 0.0, Float64(alpha), :both, :cluster_permutation, 
        ClusterInfo(:parametric, :spatiotemporal, n_perms))
        
    stat_matrix = StatMatrix(t_map, nothing)
    masks = Masks(final_mask_pos, final_mask_neg)
    clusters = Clusters(updated_pos_clusters, updated_neg_clusters)
    perm_dist = PermutationDistribution(null_dist, null_dist)
    
    # ERP Data: we can just provide the average of epochs as dummy data for plotting
    erp_dummy = EegFun.average_epochs(result.epochs)
    
    c_thresh_f64 = Float64(cluster_threshold)
    return PermutationResult(
        test_info,
        [erp_dummy, erp_dummy], # Dummy data for Cond1/Cond2
        stat_matrix,
        masks,
        clusters,
        perm_dist,
        result.channels,
        result.time_points,
        (c_thresh_f64, -c_thresh_f64),
        se_diff,
        zeros(n_electrodes, n_time_points),
        zeros(n_electrodes, n_time_points),
        se_diff
    )
end

