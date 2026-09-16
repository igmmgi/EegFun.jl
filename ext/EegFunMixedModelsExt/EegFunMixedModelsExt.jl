module EegFunMixedModelsExt

using EegFun
using MixedModels
using StatsModels
using DataFrames
using Base.Threads
using Random
using ProgressMeter

import EegFun: fit_mass_lmm

"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; n_perms=0, use_clusters=false, cluster_threshold=2.0)

Fit a linear mixed model natively across the spatial-temporal grid.
"""
function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; n_perms=0, use_clusters=false, cluster_threshold=2.0)
    dfs = epochs.data
    n_epochs = length(dfs)
    
    if n_epochs == 0
        error("EpochData is empty.")
    end
    
    # 1. Identify channels vs metadata
    first_df = dfs[1]
    all_channels = intersect(propertynames(first_df), epochs.layout.data.label)
    meta_cols = setdiff(propertynames(first_df), all_channels)
    
    n_timepoints = size(first_df, 1)
    n_channels = length(all_channels)
    
    # 2. Extract metadata into a master DataFrame
    meta_df = reduce(vcat, [df[1:1, meta_cols] for df in dfs])
    meta_df.amplitude = randn(n_epochs)
    
    # 3. Parse formula
    @info "Compiling MixedModel design matrix..."
    m_initial = fit(MixedModel, f, meta_df)
    m_initial.optsum.ftol_rel = 1e-5
    m_initial.optsum.xtol_rel = 1e-5
    
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
    beta_matrix = zeros(n_channels, n_timepoints, n_coefs)
    t_matrix = zeros(n_channels, n_timepoints, n_coefs)
    p_matrix = zeros(n_channels, n_timepoints, n_coefs)
    
    @info "Fitting $(n_channels * n_timepoints) Mixed Models across $n_epochs epochs (n_perms=$n_perms)..."
    
    # Generate sign flips for Residual Permutation
    signs = [rand([-1, 1], n_epochs) for _ in 1:n_perms]
    
    n_threads_max = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    
    # If not using clusters, track max-t on the fly. If using clusters, store full tensor.
    if use_clusters
        t_perm_matrix = zeros(n_channels, n_timepoints, n_perms, n_coefs)
    else
        thread_max_t = [zeros(n_perms, n_coefs) for _ in 1:n_threads_max]
    end
    
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    
    @info "Pre-extracting EEG data tensor..."
    eeg_tensor = zeros(n_channels, n_timepoints, n_epochs)
    for (c_idx, c_name) in enumerate(all_channels)
        for (i, df) in enumerate(dfs)
            eeg_tensor[c_idx, :, i] .= df[:, c_name]
        end
    end
    
    p = Progress(length(tasks), 1, "Fitting models...")
    
    Threads.@threads for (c_idx, t_idx) in tasks
        m_thread = get!(task_local_storage(), :mme_thread_model) do
            deepcopy(m_initial)
        end::typeof(m_initial)
        
        y_true = get!(task_local_storage(), :mme_y_true) do
            zeros(n_epochs)
        end::Vector{Float64}
        
        y_true .= @view eeg_tensor[c_idx, t_idx, :]
        
        refit!(m_thread, y_true)
        
        beta_matrix[c_idx, t_idx, :] .= coef(m_thread)
        ct = coeftable(m_thread)
        t_matrix[c_idx, t_idx, :] .= ct.cols[3]
        p_matrix[c_idx, t_idx, :] .= ct.cols[4]
        
        if n_perms > 0
            # Extract fitted values and residuals for valid Freedman-Lane style permutations
            # MixedModels fitted() and residuals() return the vectors
            y_hat = fitted(m_thread)
            res = residuals(m_thread)
            
            tid = Threads.threadid()
            y_perm = get!(task_local_storage(), :mme_y_perm) do
                zeros(n_epochs)
            end::Vector{Float64}
            
            for p in 1:n_perms
                @. y_perm = y_hat + (signs[p] * res)
                
                # Ultra-fast GLS update: Keep variance components (theta) fixed to the true model.
                # This mathematically solves the mixed model equations without NLopt.
                copyto!(m_thread.y, y_perm)
                MixedModels.unfit!(m_thread)
                MixedModels.updateL!(m_thread)
                
                t_perm = coef(m_thread) ./ stderror(m_thread)
                
                if use_clusters
                    t_perm_matrix[c_idx, t_idx, p, :] .= t_perm
                else
                    for coef_idx in 1:n_coefs
                        val = abs(t_perm[coef_idx])
                        if val > thread_max_t[tid][p, coef_idx]
                            thread_max_t[tid][p, coef_idx] = val
                        end
                    end
                end
            end
        end
        next!(p)
    end
    
    max_t_null = zeros(n_perms, n_coefs)
    max_cluster_mass_null = zeros(n_perms, n_coefs)
    
    if n_perms > 0
        if use_clusters
            @info "Extracting Spatio-Temporal Clusters (threshold = $cluster_threshold)..."
            spatial_connectivity = EegFun._build_connectivity_matrix(all_channels, epochs.layout, :spatiotemporal)
            electrode_to_idx = Dict(e => i for (i, e) in enumerate(all_channels))
            
            for p in 1:n_perms
                for coef_idx in 1:n_coefs
                    t_map = t_perm_matrix[:, :, p, coef_idx]
                    
                    mask_pos = t_map .> cluster_threshold
                    mask_neg = t_map .< -cluster_threshold
                    
                    pos_clusters, neg_clusters = EegFun._find_clusters(
                        mask_pos, mask_neg, all_channels, first_df.time, spatial_connectivity, :spatiotemporal
                    )
                    
                    pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
                    neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
                    
                    max_pos = maximum(pos_stats, init=0.0)
                    max_neg = maximum(abs, neg_stats, init=0.0)
                    
                    max_cluster_mass_null[p, coef_idx] = max(max_pos, max_neg)
                end
            end
        else
            for tid in 1:n_threads_max
                for p in 1:n_perms
                    for coef_idx in 1:n_coefs
                        if thread_max_t[tid][p, coef_idx] > max_t_null[p, coef_idx]
                            max_t_null[p, coef_idx] = thread_max_t[tid][p, coef_idx]
                        end
                    end
                end
            end
        end
    end
    
    return EegFun.LmmStatsResult(
        coef_names,
        all_channels,
        first_df.time,
        beta_matrix,
        t_matrix,
        p_matrix,
        max_t_null,
        max_cluster_mass_null,
        epochs
    )
end

end
