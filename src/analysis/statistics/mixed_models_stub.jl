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
function extract_predictor_stats(result::LmmStatsResult, coef_name::String; 
    alpha::Float64=0.05, 
    cluster_threshold::Float64=2.0,
    use_tfce::Bool=false,
    tfce_E::Float64=0.5,
    tfce_H::Float64=2.0,
    tfce_dh::Float64=0.1
)
    coef_idx = findfirst(==(coef_name), result.coefficients)
    if isnothing(coef_idx)
        error("Coefficient '\$coef_name' not found in model.")
    end

    n_electrodes = length(result.channels)
    n_time_points = length(result.time_points)

    t_map = result.t_values[:, :, coef_idx]
    
    null_dist = result.max_cluster_mass_null[:, coef_idx]
    n_perms = length(null_dist)
    
    spatial_connectivity = EegFun._build_connectivity_matrix(result.channels, result.epochs.layout, :spatiotemporal)
    
    updated_pos_clusters = Cluster[]
    updated_neg_clusters = Cluster[]
    final_mask_pos = zeros(Bool, n_electrodes, n_time_points)
    final_mask_neg = zeros(Bool, n_electrodes, n_time_points)
    tfce_map = zeros(Float64, n_electrodes, n_time_points)
    
    if use_tfce
        tfce_map = EegFun._compute_tfce(t_map, result.channels, Float64.(result.time_points), spatial_connectivity, :spatiotemporal; E=tfce_E, H=tfce_H, dh=tfce_dh)
        
        for e_idx in 1:n_electrodes
            for t_idx in 1:n_time_points
                val = tfce_map[e_idx, t_idx]
                if val != 0
                    p_val = (count(>=(abs(val)), null_dist) + 1) / (n_perms + 1)
                    if p_val <= alpha
                        if val > 0
                            final_mask_pos[e_idx, t_idx] = true
                        else
                            final_mask_neg[e_idx, t_idx] = true
                        end
                    end
                end
            end
        end
    else
        mask_pos = t_map .> cluster_threshold
        mask_neg = t_map .< -cluster_threshold
        
        pos_clusters, neg_clusters = EegFun._find_clusters(
            mask_pos, mask_neg, result.channels, result.time_points, spatial_connectivity, :spatiotemporal
        )
        
        electrode_to_idx = Dict(e => i for (i, e) in enumerate(result.channels))
        pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
        neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
        
        for (i, c) in enumerate(pos_clusters)
            p_val = (count(>=(pos_stats[i]), null_dist) + 1) / (n_perms + 1)
            push!(updated_pos_clusters, Cluster(
                c.id, c.electrodes, c.time_indices, c.time_range,
                pos_stats[i], p_val, p_val <= alpha, c.polarity, c.members
            ))
            if p_val <= alpha
                for (e_idx, t_idx) in c.members
                    final_mask_pos[e_idx, t_idx] = true
                end
            end
        end
        
        for (i, c) in enumerate(neg_clusters)
            p_val = (count(>=(abs(neg_stats[i])), null_dist) + 1) / (n_perms + 1)
            push!(updated_neg_clusters, Cluster(
                c.id, c.electrodes, c.time_indices, c.time_range,
                neg_stats[i], p_val, p_val <= alpha, c.polarity, c.members
            ))
            if p_val <= alpha
                for (e_idx, t_idx) in c.members
                    final_mask_neg[e_idx, t_idx] = true
                end
            end
        end
    end
    
    # Use stored standard errors directly from the LmmStatsResult
    se_diff = result.se[:, :, coef_idx]
    
    # Construct PermutationResult
    test_info = TestInfo(:LMM, 0.0, alpha, :both, :cluster_permutation, 
        ClusterInfo(:parametric, :spatiotemporal, n_perms))
        
    stat_matrix = StatMatrix(t_map, nothing)
    masks = Masks(final_mask_pos, final_mask_neg)
    clusters = Clusters(updated_pos_clusters, updated_neg_clusters)
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

