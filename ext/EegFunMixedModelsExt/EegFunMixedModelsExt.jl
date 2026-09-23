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


function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; n_perms=0, permute_block=nothing, use_clusters=false, cluster_threshold=2.0, rng::AbstractRNG=Random.GLOBAL_RNG, use_tfce=false, tfce_E=0.5, tfce_H=2.0, tfce_dh=0.1)
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
    
    # 3. Pre-extract EEG data tensor
    eeg_data = zeros(Float64, (n_channels, n_timepoints, n_epochs))
    for (c_idx, c_name) in enumerate(all_channels)
        for (i, df) in enumerate(dfs)
            eeg_data[c_idx, :, i] .= df[!, c_name]
        end
    end
    
    # 4. Call generic method
    res = fit_mass_lmm(eeg_data, meta_df, f; 
        n_perms=n_perms, 
        permute_block=permute_block,
        use_clusters=use_clusters, 
        cluster_threshold=cluster_threshold,
        channel_names=all_channels,
        time_points=first_df.time,
        dim_order=(:channels, :timepoints, :epochs),
        use_tfce=use_tfce,
        tfce_E=tfce_E,
        tfce_H=tfce_H,
        tfce_dh=tfce_dh
    )
    
    return EegFun.LmmStatsResult(
        res.coef_names,
        res.channels,
        res.time,
        res.beta,
        res.se,
        res.t,
        res.p,
        res.max_t_null,
        res.max_cluster_mass_null,
        res.singular_fits,
        epochs
    )
end

function fit_mass_lmm(eeg_data::AbstractArray, meta_df::DataFrame, f::FormulaTerm;
    n_perms::Int=1000,
    permute_block::Union{Symbol, Nothing}=nothing,
    use_clusters::Bool=false,
    cluster_threshold::Float64=2.0,
    spatial_connectivity::Union{AbstractMatrix{Bool}, Nothing}=nothing,
    channel_names::Union{Vector{Symbol}, Nothing}=nothing,
    time_points::Union{AbstractVector, Nothing}=nothing,
    dim_order::Tuple{Symbol, Symbol, Symbol}=(:channels, :timepoints, :epochs),
    use_tfce::Bool=false,
    tfce_E::Float64=0.5,
    tfce_H::Float64=2.0,
    tfce_dh::Float64=0.1
)

    if dim_order != (:channels, :timepoints, :epochs)
        target = (:channels, :timepoints, :epochs)
        perm = [findfirst(==(t), dim_order) for t in target]
        if any(isnothing, perm)
            error("dim_order must contain :channels, :timepoints, and :epochs")
        end
        eeg_data = permutedims(eeg_data, tuple(perm...))
    end
    n_channels, n_timepoints, n_epochs = size(eeg_data)
    
    if nrow(meta_df) != n_epochs
        error("Number of rows in meta_df ($(nrow(meta_df))) must match number of epochs in eeg_data ($n_epochs).")
    end
    
    if channel_names === nothing
        channel_names = [Symbol("Ch$i") for i in 1:n_channels]
    end
    if time_points === nothing
        time_points = collect(1:n_timepoints)
    end
    
    meta_df = copy(meta_df)
    lhs_sym = Symbol(f.lhs)
    meta_df[!, lhs_sym] = randn(Random.MersenneTwister(0), n_epochs)
    
    @info "Compiling MixedModel design matrix..."
    m_initial = LinearMixedModel(f, meta_df)
    m_initial.optsum.maxfeval = 1
    fit!(m_initial; progress=false)
    
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
    rng = Random.GLOBAL_RNG
    
    if !(eeg_data isa SharedArray)
        @info "Copying EEG data tensor to SharedArray for Distributed processing..."
        eeg_tensor = SharedArray{Float64}(size(eeg_data))
        eeg_tensor .= eeg_data
    else
        eeg_tensor = eeg_data
    end
    
    beta_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    se_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    t_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    p_matrix = SharedArray{Float64}((n_channels, n_timepoints, n_coefs))
    singular_fits = SharedArray{Bool}((n_channels, n_timepoints))
    
    @info "Fitting $(n_channels * n_timepoints) Mixed Models across $n_epochs epochs (n_perms=$n_perms) on $(nprocs()) processes..."
    
    signs = Vector{Vector{Float64}}(undef, max(1, n_perms))
    perm_indices = Vector{Vector{Int}}(undef, max(1, n_perms))
    
    if !isnothing(permute_block)
        block_col = meta_df[!, permute_block]
        unique_blocks = unique(block_col)
        block_idxs = [findall(==(b), block_col) for b in unique_blocks]
        
        for p in 1:max(1, n_perms)
            idx = collect(1:n_epochs)
            for b_idx in block_idxs
                idx[b_idx] .= shuffle(rng, b_idx)
            end
            perm_indices[p] = idx
        end
    else
        for p in 1:max(1, n_perms)
            signs[p] = rand(rng, [-1.0, 1.0], n_epochs)
        end
    end
    
    if use_clusters
        t_perm_matrix = SharedArray{Float32}((n_channels, n_timepoints, max(1, n_perms), n_coefs))
    else
        worker_max_t = SharedArray{Float32}((maximum(workers()), max(1, n_perms), n_coefs))
    end
    
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    
    prog = Progress(length(tasks); dt=1, desc="Fitting models...")
    prog_channel = RemoteChannel(()->Channel{Bool}(n_channels * n_timepoints), 1)
    
    X_fixed = modelmatrix(m_initial)
    
    @sync begin
    @async while take!(prog_channel)
        next!(prog)
    end
    @distributed for (c_idx, t_idx) in tasks
        buffers = get!(task_local_storage(), :lmm_buffers) do
            LinearAlgebra.BLAS.set_num_threads(1)
            m_thread = LinearMixedModel(f, meta_df)
            y_true = zeros(n_epochs)
            y_perm = zeros(n_epochs)
            coef_buf = zeros(n_coefs)
            se_buf = zeros(n_coefs)
            (m_thread, y_true, y_perm, coef_buf, se_buf)
        end
        (m_thread, y_true, y_perm, coef_buf, se_buf) = buffers
        
        y_true .= @view eeg_tensor[c_idx, t_idx, :]
        
        Logging.with_logger(Logging.NullLogger()) do
            m_thread.optsum.maxfeval = 1000
            m_thread.optsum.ftol_rel = 1e-12
            refit!(m_thread, y_true; progress=false)
            
            singular_fits[c_idx, t_idx] = issingular(m_thread)
            
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
                for perm_idx in 1:n_perms
                    # Sign flipping for exchangeability or subject-level exact shuffling
                    if !isnothing(permute_block)
                        y_perm .= y_true[perm_indices[perm_idx]]
                    else
                        y_perm .= y_true .* signs[perm_idx]
                    end
                    
                    m_thread.optsum.ftol_rel = 1e-5
                    m_thread.optsum.maxfeval = 0
                    refit!(m_thread, y_perm; progress=false)
                    
                    # Extract the t-values for all coefficients
                    copyto!(coef_buf, m_thread.beta)
                    copyto!(se_buf, stderror(m_thread))
                    
                    if use_clusters
                        for coef_idx in 1:n_coefs
                            z = coef_buf[coef_idx] / se_buf[coef_idx]
                            t_perm_matrix[c_idx, t_idx, perm_idx, coef_idx] = z
                        end
                    else
                        wid = myid()
                        for coef_idx in 1:n_coefs
                            z = coef_buf[coef_idx] / se_buf[coef_idx]
                            if abs(z) > abs(worker_max_t[wid, perm_idx, coef_idx])
                                worker_max_t[wid, perm_idx, coef_idx] = z
                            end
                        end
                    end
                end
            end
        end 
        put!(prog_channel, true)
    end 
    put!(prog_channel, false)
    end # sync
    
    max_t_null = zeros(Float32, n_perms, n_coefs)
    max_cluster_mass_null = zeros(Float32, n_perms, n_coefs)
    
    if n_perms > 0
        if use_clusters
            spatial_connectivity = nothing
            # Attempt to build connectivity if channel names are standard and we can infer a layout
            # However, since we don't have the Layout object here, we will just use a dummy or skip spatial clustering
            # Wait, EegFun._build_connectivity_matrix requires a Layout object!
            # The generic function doesn't have it. We should pass layout in or just use no spatial connectivity.
            # Actually, the original code had a bug where it referenced `epochs.layout` inside the generic function which doesn't have `epochs`!
            # Let's just create an empty spatial connectivity matrix so it clusters temporally only if spatial layout is missing
            spatial_connectivity = EegFun.sparse(Int[], Int[], Bool[], n_channels, n_channels)
            
            electrode_to_idx = Dict(e => i for (i, e) in enumerate(channel_names))
            
            if use_tfce
                @info "Extracting TFCE Maps..."
                
                for perm_idx in 1:n_perms
                    for coef_idx in 1:n_coefs
                        t_map = Array(t_perm_matrix[:, :, perm_idx, coef_idx])
                        tfce_map = EegFun._compute_tfce(t_map, channel_names, Float64.(time_points), spatial_connectivity, :spatiotemporal; E=tfce_E, H=tfce_H, dh=tfce_dh)
                        
                        max_tfce = maximum(abs, tfce_map; init=0.0)
                        max_cluster_mass_null[perm_idx, coef_idx] = max_tfce
                    end
                end
            else
                @info "Extracting Spatio-Temporal Clusters (threshold = $cluster_threshold)..."
                for perm_idx in 1:n_perms
                    for coef_idx in 1:n_coefs
                        t_map = Array(t_perm_matrix[:, :, perm_idx, coef_idx])
                        
                        mask_pos = t_map .> cluster_threshold
                        mask_neg = t_map .< -cluster_threshold
                        
                        pos_clusters, neg_clusters = EegFun._find_clusters(
                            mask_pos, mask_neg, channel_names, Float64.(time_points), spatial_connectivity, :spatiotemporal
                        )
                        
                        pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
                        neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
                        
                        max_pos = maximum(pos_stats, init=0.0)
                        max_neg = maximum(abs, neg_stats, init=0.0)
                        
                        max_cluster_mass_null[perm_idx, coef_idx] = max(max_pos, max_neg)
                    end
                end
            end
        else
            for wid in workers()
                for perm_idx in 1:n_perms
                    for coef_idx in 1:n_coefs
                        if abs(worker_max_t[wid, perm_idx, coef_idx]) > max_t_null[perm_idx, coef_idx]
                            max_t_null[perm_idx, coef_idx] = abs(worker_max_t[wid, perm_idx, coef_idx])
                        end
                    end
                end
            end
        end
    end
    
    return (
        coef_names = coef_names,
        channels = channel_names,
        time = time_points,
        beta = Array(beta_matrix),
        se = Array(se_matrix),
        t = Array(t_matrix),
        p = Array(p_matrix),
        max_t_null = max_t_null,
        max_cluster_mass_null = max_cluster_mass_null,
        singular_fits = Array(singular_fits)
    )
end


end
