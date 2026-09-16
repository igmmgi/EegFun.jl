"""
    fit_mass_lmm(args...; kwargs...)

Fit a mass-univariate linear mixed model across all channels and timepoints.

**Requires `MixedModels.jl` and `StatsModels.jl` to be loaded.**

To use this function, you must first import the required packages:
```julia
using EegFun
using MixedModels, StatsModels
```
"""
function fit_mass_lmm(args...; kwargs...)
    error("To use Mass-Univariate Mixed Models, you must first load the packages: `using MixedModels, StatsModels`")
end

"""
    extract_predictor_stats(result::LmmStatsResult, coef_name::String; alpha::Float64=0.05, cluster_threshold::Float64=2.0)

Extracts the cluster-corrected statistics for a specific LMM predictor and formats them as a standard `PermutationResult`.
This allows the result to be plotted directly using `plot_erp_stats`, `plot_topography_stats`, etc.
"""
function extract_predictor_stats(result::LmmStatsResult, coef_name::String; alpha::Float64=0.05, cluster_threshold::Float64=2.0)
    coef_idx = findfirst(==(coef_name), result.coefficients)
    if isnothing(coef_idx)
        error("Coefficient '\$coef_name' not found in model.")
    end

    n_electrodes = length(result.channels)
    n_time_points = length(result.time_points)

    t_map = result.t_values[:, :, coef_idx]
    
    # 1. Compute significance masks
    mask_pos = t_map .> cluster_threshold
    mask_neg = t_map .< -cluster_threshold
    
    # Need to find clusters in the TRUE data
    spatial_connectivity = EegFun._build_connectivity_matrix(result.channels, result.epochs.layout, :spatiotemporal)
    pos_clusters, neg_clusters = EegFun._find_clusters(
        mask_pos, mask_neg, result.channels, result.time_points, spatial_connectivity, :spatiotemporal
    )
    
    electrode_to_idx = Dict(e => i for (i, e) in enumerate(result.channels))
    
    # Calculate cluster stats
    pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
    neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
    
    # Compute p-values based on null distribution
    null_dist = result.max_cluster_mass_null[:, coef_idx]
    
    for (i, c) in enumerate(pos_clusters)
        p_val = count(>=(pos_stats[i]), null_dist) / length(null_dist)
        c.p_value = p_val
        c.is_significant = p_val <= alpha
    end
    
    for (i, c) in enumerate(neg_clusters)
        p_val = count(>=(abs(neg_stats[i])), null_dist) / length(null_dist)
        c.p_value = p_val
        c.is_significant = p_val <= alpha
    end
    
    # Create significance masks based on cluster p-values
    final_mask_pos = zeros(Bool, n_electrodes, n_time_points)
    final_mask_neg = zeros(Bool, n_electrodes, n_time_points)
    
    for c in pos_clusters
        if c.is_significant
            for pt in c.points
                ch_idx = electrode_to_idx[pt.electrode]
                t_idx = findfirst(==(pt.time), result.time_points)
                final_mask_pos[ch_idx, t_idx] = true
            end
        end
    end
    
    for c in neg_clusters
        if c.is_significant
            for pt in c.points
                ch_idx = electrode_to_idx[pt.electrode]
                t_idx = findfirst(==(pt.time), result.time_points)
                final_mask_neg[ch_idx, t_idx] = true
            end
        end
    end
    
    # Fake standard error (we just use t-values to infer se_diff if we really need it, but LMM doesn't export se directly in the struct right now)
    # se_diff = beta / t_value
    beta_map = result.beta[:, :, coef_idx]
    se_diff = abs.(beta_map ./ t_map)
    se_diff[isnan.(se_diff)] .= 0.0
    se_diff[isinf.(se_diff)] .= 0.0
    
    # Construct PermutationResult
    test_info = TestInfo(:LMM, 0.0, alpha, :both, :cluster_permutation, 
        ClusterInfo(:parametric, :spatiotemporal, length(null_dist)))
        
    stat_matrix = StatMatrix(t_map, nothing)
    masks = Masks(final_mask_pos, final_mask_neg)
    clusters = Clusters(pos_clusters, neg_clusters)
    perm_dist = PermutationDistribution(null_dist, null_dist)
    
    # ERP Data: we can just provide the average of epochs as dummy data for plotting
    erp_dummy = EegFun.average_epochs(result.epochs)
    
    return PermutationResult(
        test_info,
        [erp_dummy, erp_dummy], # Dummy data for Cond1/Cond2
        stat_matrix,
        masks,
        clusters,
        perm_dist,
        result.channels,
        result.time_points,
        (cluster_threshold, -cluster_threshold),
        se_diff,
        zeros(n_electrodes, n_time_points),
        zeros(n_electrodes, n_time_points),
        se_diff
    )
end
