module EegFunMixedModelsExt

using EegFun
using MixedModels
using StatsModels
using DataFrames
using Base.Threads
using Random
using ProgressMeter
using Distributions
using LinearAlgebra
using Logging
using Distributed
using SharedArrays

import EegFun: fit_mass_lmm

"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; n_perms=0, use_clusters=false, cluster_threshold=2.0, rng=Random.GLOBAL_RNG)

Fit a linear mixed model natively across the spatial-temporal grid using Distributed processing.
"""
function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; n_perms=0, use_clusters=false, cluster_threshold=2.0, rng::AbstractRNG=Random.GLOBAL_RNG, fast::Bool=false)
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
    meta_df.amplitude = randn(Random.MersenneTwister(0), n_epochs)
    
    # 3. Parse formula and compile design matrix
    @info "Compiling MixedModel design matrix..."
    m_initial = LinearMixedModel(f, meta_df)
    m_initial.optsum.maxfeval = 1
    fit!(m_initial; progress=false)
    
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
    # Pre-extract EEG data tensor into a SharedArray so workers don't copy it
    @info "Pre-extracting EEG data tensor to SharedArray..."
    eeg_tensor = SharedArray{Float64}((n_channels, n_timepoints, n_epochs))
    for (c_idx, c_name) in enumerate(all_channels)
        for (i, df) in enumerate(dfs)
            eeg_tensor[c_idx, :, i] .= df[!, c_name]
        end
    end
    
    beta_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    se_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    t_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    p_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    
    @info "Fitting $(n_channels * n_timepoints) Mixed Models across $n_epochs epochs (n_perms=$n_perms) on $(nprocs()) processes..."
    
    if isempty(m_initial.reterms)
        error("Model has no random effects. fit_mass_lmm requires a mixed model.")
    end
    # Extract the primary grouping factor (e.g., subject) to respect exchangeability blocks
    group_refs = m_initial.reterms[1].refs
    n_groups = length(m_initial.reterms[1].levels)
    signs = [rand(rng, [-1, 1], n_groups) for _ in 1:n_perms]    
    if use_clusters
        t_perm_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_perms, n_coefs))
    else
        worker_max_t = SharedArray{Float64}((maximum(workers()), n_perms, n_coefs))
    end
    
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    
    prog = Progress(length(tasks), 1, "Fitting models...")
    prog_channel = RemoteChannel(()->Channel{Bool}(n_channels * n_timepoints), 1)
    
    # Extract fixed model matrix once to avoid allocations
    X_fixed = modelmatrix(m_initial)
    
    # Run Distributed loop
    @sync begin
    @async while take!(prog_channel)
        next!(prog)
    end
    @distributed for (c_idx, t_idx) in tasks
        # Cache worker-local allocations in task_local_storage to reuse them across iterations
        buffers = get!(task_local_storage(), :lmm_buffers) do
            # Disable BLAS threads on this worker process
            LinearAlgebra.BLAS.set_num_threads(1)
            
            # Allocate buffers
            m_thread = LinearMixedModel(f, meta_df)
            y_true = zeros(n_epochs)
            y_perm = zeros(n_epochs)
            coef_buf = zeros(n_coefs)
            se_buf = zeros(n_coefs)
            res_buf = zeros(n_epochs)
            t_perm_buf = zeros(n_coefs)
            
            (m_thread, y_true, y_perm, coef_buf, se_buf, res_buf, t_perm_buf)
        end
        (m_thread, y_true, y_perm, coef_buf, se_buf, res_buf, t_perm_buf) = buffers
        
        y_true .= @view eeg_tensor[c_idx, t_idx, :]
        
        Logging.with_logger(Logging.NullLogger()) do
            m_thread.optsum.maxfeval = 1000
            m_thread.optsum.maxtime = 1.0
            refit!(m_thread, y_true; progress=false)
            
            # Extract optimal theta to use as frozen rulebook or warm start
            true_theta = copy(m_thread.theta)
            
            copyto!(coef_buf, m_thread.beta)
            copyto!(se_buf, stderror(m_thread))
            
            beta_matrix[c_idx, t_idx, :] .= coef_buf
            se_matrix[c_idx, t_idx, :] .= se_buf
            
            for coef_idx in 1:n_coefs
                z = coef_buf[coef_idx] / se_buf[coef_idx]
                t_matrix[c_idx, t_idx, coef_idx] = z
                p_matrix[c_idx, t_idx, coef_idx] = 2.0 * ccdf(Normal(), abs(z))
            end
            
            if n_perms > 0
                fitted_vals = fitted(m_thread)
                for i in 1:n_epochs
                    res_buf[i] = y_true[i] - fitted_vals[i]
                end
                
                X = X_fixed
                
                for perm_idx in 1:n_perms
                    curr_signs = signs[perm_idx]
                    
                    # ter Braak method: Permute for each partial effect separately
                    for test_coef in 1:n_coefs
                        for i in 1:n_epochs
                            # y* = y_hat_{reduced} + permuted_residuals
                            # where y_hat_{reduced} = y_hat - X_k * beta_k
                            y_perm[i] = (fitted_vals[i] - X[i, test_coef] * coef_buf[test_coef]) + curr_signs[group_refs[i]] * res_buf[i]
                        end
                        
                        m_thread.optsum.initial .= true_theta
                        m_thread.optsum.ftol_rel = 1e-5
                        m_thread.optsum.maxfeval = fast ? 0 : 1000
                        m_thread.optsum.maxtime = 1.0
                        refit!(m_thread, y_perm; progress=false)
                        
                        # Extract permuted t-value for this specific coefficient
                        t_perm_buf[test_coef] = m_thread.beta[test_coef] / stderror(m_thread)[test_coef]
                    end
                    
                    if use_clusters
                        t_perm_matrix[c_idx, t_idx, perm_idx, :] .= t_perm_buf
                    else
                        wid = myid()
                        for coef_idx in 1:n_coefs
                            val = abs(t_perm_buf[coef_idx])
                            if val > worker_max_t[wid, perm_idx, coef_idx]
                                worker_max_t[wid, perm_idx, coef_idx] = val
                            end
                        end
                    end
                end
            end
        end # with_logger
        put!(prog_channel, true)
    end # distributed
    put!(prog_channel, false)
    end # sync
    
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
            for wid in 1:size(worker_max_t, 1)
                for perm_idx in 1:n_perms
                    for coef_idx in 1:n_coefs
                        if worker_max_t[wid, perm_idx, coef_idx] > max_t_null[perm_idx, coef_idx]
                            max_t_null[perm_idx, coef_idx] = worker_max_t[wid, perm_idx, coef_idx]
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
        Array(beta_matrix),
        Array(se_matrix),
        Array(t_matrix),
        Array(p_matrix),
        max_t_null,
        max_cluster_mass_null,
        epochs
    )
end

end
