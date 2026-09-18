module EegFunMixedModelsExt

using EegFun
using MixedModels
using StatsModels
using DataFrames
using Base.Threads
using Random
using ProgressMeter
using Distributions

import EegFun: fit_mass_lmm

"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; n_perms=0, use_clusters=false, cluster_threshold=2.0, rng=Random.GLOBAL_RNG)

Fit a linear mixed model natively across the spatial-temporal grid.

# Arguments
- `epochs::EpochData`: The epoch data to fit.
- `f::FormulaTerm`: The model formula (response is replaced internally).
- `n_perms::Int=0`: Number of Freedman-Lane permutations. 0 = no permutation test.
- `use_clusters::Bool=false`: If true, store full permutation t-maps for cluster-based inference.
- `cluster_threshold::Float64=2.0`: t-value threshold for cluster formation.
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator for reproducible permutations.
"""
function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; n_perms=0, use_clusters=false, cluster_threshold=2.0, rng::AbstractRNG=Random.GLOBAL_RNG)
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
    
    # 2. Extract metadata into a master DataFrame — use zeros for deterministic initial fit
    meta_df = reduce(vcat, [df[1:1, meta_cols] for df in dfs])
    meta_df.amplitude = randn(Random.MersenneTwister(0), n_epochs)
    
    # 3. Parse formula and compile design matrix
    @info "Compiling MixedModel design matrix..."
    m_initial = fit(MixedModel, f, meta_df)
    
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
    beta_matrix = zeros(n_channels, n_timepoints, n_coefs)
    se_matrix = zeros(n_channels, n_timepoints, n_coefs)
    t_matrix = zeros(n_channels, n_timepoints, n_coefs)
    p_matrix = zeros(n_channels, n_timepoints, n_coefs)
    
    @info "Fitting $(n_channels * n_timepoints) Mixed Models across $n_epochs epochs (n_perms=$n_perms)..."
    
    # Generate sign flips for Freedman-Lane permutation with controlled RNG
    signs = [rand(rng, [-1, 1], n_epochs) for _ in 1:n_perms]
    
    n_threads_max = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    
    # If not using clusters, track max-t on the fly. If using clusters, store full tensor.
    if use_clusters
        t_perm_matrix = zeros(n_channels, n_timepoints, n_perms, n_coefs)
    else
        thread_max_t = [zeros(n_perms, n_coefs) for _ in 1:n_threads_max]
    end
    
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    
    # Pre-extract EEG data tensor — use df[!, col] (no-copy view) instead of df[:, col]
    @info "Pre-extracting EEG data tensor..."
    eeg_tensor = zeros(n_channels, n_timepoints, n_epochs)
    for (c_idx, c_name) in enumerate(all_channels)
        for (i, df) in enumerate(dfs)
            eeg_tensor[c_idx, :, i] .= df[!, c_name]
        end
    end
    
    # Pre-allocate per-thread buffers to avoid allocations in the hot loop
    m_threads = [deepcopy(m_initial) for _ in 1:n_threads_max]
    y_true_threads = [zeros(n_epochs) for _ in 1:n_threads_max]
    y_perm_threads = [zeros(n_epochs) for _ in 1:n_threads_max]
    coef_buf_threads = [zeros(n_coefs) for _ in 1:n_threads_max]
    se_buf_threads = [zeros(n_coefs) for _ in 1:n_threads_max]
    y_hat_threads = [zeros(n_epochs) for _ in 1:n_threads_max]
    res_threads = [zeros(n_epochs) for _ in 1:n_threads_max]
    t_perm_buf_threads = [zeros(n_coefs) for _ in 1:n_threads_max]
    
    prog = Progress(length(tasks), 1, "Fitting models...")
    
    Threads.@threads :static for (c_idx, t_idx) in tasks
        tid = Threads.threadid()
        m_thread = m_threads[tid]
        y_true = y_true_threads[tid]
        y_perm = y_perm_threads[tid]
        coef_buf = coef_buf_threads[tid]
        se_buf = se_buf_threads[tid]
        y_hat = y_hat_threads[tid]
        res_buf = res_threads[tid]
        t_perm_buf = t_perm_buf_threads[tid]
        
        y_true .= @view eeg_tensor[c_idx, t_idx, :]
        
        # Fit the true model
        refit!(m_thread, y_true)
        
        # Extract coefficients and standard errors into pre-allocated buffers
        copyto!(coef_buf, coef(m_thread))
        copyto!(se_buf, stderror(m_thread))
        
        beta_matrix[c_idx, t_idx, :] .= coef_buf
        se_matrix[c_idx, t_idx, :] .= se_buf
        
        # Compute Wald z-tests directly (bypasses expensive coeftable string allocations)
        for coef_idx in 1:n_coefs
            z = coef_buf[coef_idx] / se_buf[coef_idx]
            t_matrix[c_idx, t_idx, coef_idx] = z
            p_matrix[c_idx, t_idx, coef_idx] = 2.0 * ccdf(Normal(), abs(z))
        end
        
        if n_perms > 0
            # Extract fitted values and residuals for Freedman-Lane permutations
            copyto!(y_hat, fitted(m_thread))
            copyto!(res_buf, residuals(m_thread))
            
            for perm_idx in 1:n_perms
                @. y_perm = y_hat + (signs[perm_idx] * res_buf)
                
                # Full NLopt refit — re-estimates variance components (rigorous MLE)
                refit!(m_thread, y_perm)
                
                copyto!(coef_buf, coef(m_thread))
                copyto!(se_buf, stderror(m_thread))
                
                for coef_idx in 1:n_coefs
                    t_perm_buf[coef_idx] = coef_buf[coef_idx] / se_buf[coef_idx]
                end
                
                if use_clusters
                    t_perm_matrix[c_idx, t_idx, perm_idx, :] .= t_perm_buf
                else
                    for coef_idx in 1:n_coefs
                        val = abs(t_perm_buf[coef_idx])
                        if val > thread_max_t[tid][perm_idx, coef_idx]
                            thread_max_t[tid][perm_idx, coef_idx] = val
                        end
                    end
                end
            end
        end
        next!(prog)
    end
    
    max_t_null = zeros(n_perms, n_coefs)
    max_cluster_mass_null = zeros(n_perms, n_coefs)
    
    if n_perms > 0
        if use_clusters
            @info "Extracting Spatio-Temporal Clusters (threshold = $cluster_threshold)..."
            spatial_connectivity = EegFun._build_connectivity_matrix(all_channels, epochs.layout, :spatiotemporal)
            electrode_to_idx = Dict(e => i for (i, e) in enumerate(all_channels))
            
            for perm_idx in 1:n_perms
                for coef_idx in 1:n_coefs
                    t_map = t_perm_matrix[:, :, perm_idx, coef_idx]
                    
                    mask_pos = t_map .> cluster_threshold
                    mask_neg = t_map .< -cluster_threshold
                    
                    pos_clusters, neg_clusters = EegFun._find_clusters(
                        mask_pos, mask_neg, all_channels, Float64.(first_df.time), spatial_connectivity, :spatiotemporal
                    )
                    
                    pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
                    neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
                    
                    max_pos = maximum(pos_stats, init=0.0)
                    max_neg = maximum(abs, neg_stats, init=0.0)
                    
                    max_cluster_mass_null[perm_idx, coef_idx] = max(max_pos, max_neg)
                end
            end
        else
            for tid in 1:n_threads_max
                for perm_idx in 1:n_perms
                    for coef_idx in 1:n_coefs
                        if thread_max_t[tid][perm_idx, coef_idx] > max_t_null[perm_idx, coef_idx]
                            max_t_null[perm_idx, coef_idx] = thread_max_t[tid][perm_idx, coef_idx]
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
        se_matrix,
        t_matrix,
        p_matrix,
        max_t_null,
        max_cluster_mass_null,
        epochs
    )
end

end
