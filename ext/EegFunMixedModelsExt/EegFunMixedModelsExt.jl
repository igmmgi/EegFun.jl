module EegFunMixedModelsExt

using EegFun
using MixedModels
using StatsModels
using DataFrames
using Base.Threads
using Random

import EegFun: fit_mass_lmm

function generate_permutation_indices(meta_df::DataFrame, group::Symbol, n_perms::Int)
    n = nrow(meta_df)
    perm_indices = [collect(1:n) for _ in 1:n_perms]
    
    if group !== :none && string(group) in names(meta_df)
        group_col = meta_df[!, group]
        groups = unique(group_col)
        group_indices = [findall(x -> x == g, group_col) for g in groups]
        
        for p in 1:n_perms
            for idxs in group_indices
                # Shuffle the indices within this group
                perm_indices[p][idxs] = Random.shuffle(idxs)
            end
        end
    else
        for p in 1:n_perms
            Random.shuffle!(perm_indices[p])
        end
    end
    return perm_indices
end

"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; correction=:fdr, n_perms=0, group=:subject)

Fit a linear mixed model natively across the spatial-temporal grid.
"""
function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; correction=:fdr, n_perms=0, group=:subject)
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
    
    # 2. Extract metadata into a master DataFrame (one row per epoch)
    # Assuming metadata is constant across time within an epoch
    meta_df = DataFrame()
    for df in dfs
        append!(meta_df, DataFrame(df[1:1, meta_cols]))
    end
    
    # Add dummy amplitude column
    meta_df.amplitude = randn(n_epochs)
    
    # 3. Parse the formula and build the initial model
    # We fit the model once to build the design matrices X and Z
    @info "Compiling MixedModel design matrix..."
    m_initial = fit(MixedModel, f, meta_df)
    
    # Extract coefficient names to build result matrices
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
    # Pre-allocate 3D arrays for coefficients, t-values, and p-values
    # Dimensions: (Channels × Timepoints × Coefficients)
    beta_matrix = zeros(n_channels, n_timepoints, n_coefs)
    t_matrix = zeros(n_channels, n_timepoints, n_coefs)
    p_matrix = zeros(n_channels, n_timepoints, n_coefs)
    
    @info "Fitting $(n_channels * n_timepoints) Mixed Models across $n_epochs epochs (n_perms=$n_perms)..."
    
    perm_indices = n_perms > 0 ? generate_permutation_indices(meta_df, group, n_perms) : []
    
    n_threads_max = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    thread_max_t = [zeros(n_perms, n_coefs) for _ in 1:n_threads_max]
    
    # 4. Mass-Univariate Loop
    # We will use task_local_storage to securely maintain a deepcopy per task
    
    # Create a vector of tasks (channel, timepoint)
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    
    Threads.@threads for (c_idx, t_idx) in tasks
        m_thread = get!(task_local_storage(), :mme_thread_model) do
            deepcopy(m_initial)
        end::typeof(m_initial)
        
        # Extract true response vector (amplitude) for this specific channel/time across all epochs
        y_true = zeros(n_epochs)
        c_name = all_channels[c_idx]
        for (i, df) in enumerate(dfs)
            y_true[i] = df[t_idx, c_name]
        end
        
        # Lightning fast refit (avoids parsing formula and rebuilding matrices)
        refit!(m_thread, y_true)
        
        # Extract statistics for observed data
        beta_matrix[c_idx, t_idx, :] .= coef(m_thread)
        
        ct = coeftable(m_thread)
        t_matrix[c_idx, t_idx, :] .= ct.cols[3]
        p_matrix[c_idx, t_idx, :] .= ct.cols[4]
        
        # Process Permutations
        if n_perms > 0
            tid = Threads.threadid()
            y_perm = zeros(n_epochs)
            for p in 1:n_perms
                for i in 1:n_epochs
                    y_perm[i] = y_true[perm_indices[p][i]]
                end
                
                refit!(m_thread, y_perm)
                
                # Update local max absolute t-value for FWER control
                t_perm = coeftable(m_thread).cols[3]
                for coef_idx in 1:n_coefs
                    val = abs(t_perm[coef_idx])
                    if val > thread_max_t[tid][p, coef_idx]
                        thread_max_t[tid][p, coef_idx] = val
                    end
                end
            end
        end
    end
    
    # Reduce permutation distributions across threads
    max_t_null = zeros(n_perms, n_coefs)
    if n_perms > 0
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
    
    # In a full implementation, we would wrap this in a StatisticalData object
    # For now, we return a NamedTuple of the results
    return (
        coefficients = coef_names,
        channels = all_channels,
        timepoints = first_df.time,
        beta = beta_matrix,
        t_values = t_matrix,
        p_values = p_matrix,
        max_t_null = max_t_null
    )
end

end
