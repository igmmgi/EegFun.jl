module EegFunMixedModelsExt

using EegFun
using MixedModels
using StatsModels
using DataFrames
using Random
using ProgressMeter
using Distributions
using LinearAlgebra
using Logging
using Distributed
using SharedArrays
using SparseArrays

import EegFun: fit_mass_lmm, generate_permutation_matrix

# Relative objective tolerance for the optimizer in permutation refits (the observed-data fit uses 1e-12).
# Roughly 2x faster than 1e-12, with negligible change in the permutation t-statistics.
const PERMUTATION_FTOL_REL = 1e-8

function _resolve_tested_coefs(coef_names::Vector{String}, tested_coefs, test_intercept::Bool)
    n_coefs = length(coef_names)
    intercept_idx = findfirst(x -> x in ("(Intercept)", "Intercept", ":(Intercept)"), coef_names)

    if tested_coefs === nothing || tested_coefs == :effects || tested_coefs == :auto
        if test_intercept
            return collect(1:n_coefs)
        else
            effects = isnothing(intercept_idx) ? collect(1:n_coefs) : filter(!=(intercept_idx), 1:n_coefs)
            return isempty(effects) ? [intercept_idx] : effects
        end
    elseif tested_coefs == :all
        return collect(1:n_coefs)
    elseif tested_coefs isa Integer
        if !(1 <= tested_coefs <= n_coefs)
            error("tested_coefs index $tested_coefs out of bounds (1:$n_coefs)")
        end
        idxs = [Int(tested_coefs)]
        if test_intercept && !isnothing(intercept_idx) && !(intercept_idx in idxs)
            push!(idxs, intercept_idx)
        end
        return sort(idxs)
    elseif tested_coefs isa Union{Symbol, AbstractString}
        target_str = String(tested_coefs)
        matching = findall(x -> occursin(lowercase(target_str), lowercase(x)), coef_names)
        if isempty(matching)
            error("No model coefficient matches '$target_str'. Available coefficients: $(coef_names)")
        end
        if test_intercept && !isnothing(intercept_idx) && !(intercept_idx in matching)
            push!(matching, intercept_idx)
        end
        return sort(unique(matching))
    elseif tested_coefs isa Union{AbstractVector, Tuple}
        idxs = Int[]
        for item in tested_coefs
            if item isa Integer
                if !(1 <= item <= n_coefs)
                    error("tested_coefs index $item out of bounds (1:$n_coefs)")
                end
                push!(idxs, Int(item))
            elseif item isa Union{Symbol, AbstractString}
                item_str = String(item)
                matching = findall(x -> occursin(lowercase(item_str), lowercase(x)), coef_names)
                if isempty(matching)
                    error("No model coefficient matches '$item_str'. Available coefficients: $(coef_names)")
                end
                append!(idxs, matching)
            else
                error("Unsupported item in tested_coefs: $item. Must be String, Symbol, or Integer.")
            end
        end
        if test_intercept && !isnothing(intercept_idx) && !(intercept_idx in idxs)
            push!(idxs, intercept_idx)
        end
        return sort(unique(idxs))
    else
        error("Unsupported type for tested_coefs: $(typeof(tested_coefs)). Expected Symbol, String, Integer, or Vector.")
    end
end



function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; 
    n_perms=0, 
    tested_coefs=nothing,
    test_intercept::Bool=false,
    permute_block=nothing, 
    permute_crossed=nothing,
    perm_matrix=nothing, 
    use_clusters=false, 
    cluster_threshold=2.0, 
    spatial_connectivity=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG, 
    use_tfce=false, 
    tfce_E=0.5, 
    tfce_H=2.0, 
    tfce_dh=0.1
)
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
    
    spatial_conn = if !isnothing(spatial_connectivity)
        spatial_connectivity
    elseif use_clusters && !isnothing(epochs.layout)
        EegFun._build_connectivity_matrix(all_channels, epochs.layout, :spatiotemporal)
    else
        nothing
    end
    
    # 4. Call generic method
    res = fit_mass_lmm(eeg_data, meta_df, f; 
        n_perms=n_perms, 
        tested_coefs=tested_coefs,
        test_intercept=test_intercept,
        permute_block=permute_block,
        permute_crossed=permute_crossed,
        perm_matrix=perm_matrix,
        use_clusters=use_clusters, 
        cluster_threshold=cluster_threshold,
        spatial_connectivity=spatial_conn,
        channel_names=all_channels,
        time_points=first_df.time,
        dim_order=(:channels, :timepoints, :epochs),
        use_tfce=use_tfce,
        tfce_E=tfce_E,
        tfce_H=tfce_H,
        tfce_dh=tfce_dh,
        rng=rng
    )
    
    return EegFun.LmmStatsResult(
        res.coef_names,
        res.channels,
        res.time,
        res.beta,
        res.se,
        res.t,
        res.p,
        res.p_uncorrected,
        res.p_corrected,
        res.max_t_null,
        res.max_cluster_mass_null,
        res.singular_fits,
        epochs
    )
end

function fit_mass_lmm(eeg_data::AbstractArray, meta_df::DataFrame, f::FormulaTerm;
    n_perms::Int=1000,
    tested_coefs=nothing,
    test_intercept::Bool=false,
    permute_block::Union{Symbol, Nothing}=nothing,
    permute_crossed::Union{Tuple{Symbol, Symbol}, Vector{Symbol}, Nothing}=nothing,
    perm_matrix::Union{AbstractMatrix{Int}, Nothing}=nothing,
    use_clusters::Bool=false,
    cluster_threshold::Float64=2.0,
    spatial_connectivity::Union{AbstractMatrix{Bool}, Nothing}=nothing,
    channel_names::Union{Vector{Symbol}, Nothing}=nothing,
    time_points::Union{AbstractVector, Nothing}=nothing,
    dim_order::Tuple{Symbol, Symbol, Symbol}=(:channels, :timepoints, :epochs),
    use_tfce::Bool=false,
    tfce_E::Float64=0.5,
    tfce_H::Float64=2.0,
    tfce_dh::Float64=0.1,
    rng::AbstractRNG=Random.GLOBAL_RNG
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
    # MixedModels' own rank-deficiency warning is silenced here; a single, more
    # informative EegFun warning is issued below instead.
    m_initial = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df)
    end
    m_initial.optsum.maxfeval = 1
    fit!(m_initial; progress=false)
    
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
    # Rank deficiency: the fixed-effects design is identical for every pixel, so it
    # can be detected once. MixedModels pivots aliased columns to the end and drops
    # them; `pivot[(rank+1):end]` gives their indices in the original coefficient order.
    rank_X = m_initial.feterm.rank
    dropped_coefs = sort(MixedModels.pivot(m_initial)[(rank_X + 1):end])
    if !isempty(dropped_coefs)
        @warn "Fixed-effects design matrix is rank deficient (rank $rank_X of $n_coefs). " *
              "Aliased coefficient(s) $(join(coef_names[dropped_coefs], ", ")) cannot be estimated: " *
              "their beta, se, t and p values are NaN and they are excluded from permutation testing. " *
              "Common causes: a predictor that is constant in the data, an empty design cell combined " *
              "with an interaction, or exactly collinear predictors."
    end
    
    tested_coef_indices = _resolve_tested_coefs(coef_names, tested_coefs, test_intercept)
    untestable = intersect(tested_coef_indices, dropped_coefs)
    if !isempty(untestable)
        tested_coef_indices = setdiff(tested_coef_indices, dropped_coefs)
        if n_perms > 0 || !isnothing(perm_matrix)
            @info "Skipping permutation testing for aliased coefficient(s): $(join(coef_names[untestable], ", "))"
        end
    end
    if n_perms > 0
        @info "Permutation testing will be computed for $(length(tested_coef_indices)) of $n_coefs coefficients: $(join(coef_names[tested_coef_indices], ", ")) (test_intercept=$test_intercept)"
    end
    
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
    
    is_permutation = !isnothing(perm_matrix) || !isnothing(permute_block) || !isnothing(permute_crossed)
    
    if is_permutation
        perm_indices = SharedArray{Int}((n_epochs, max(1, n_perms)))
        signs = SharedArray{Float32}((0, 0))
    else
        perm_indices = SharedArray{Int}((0, 0))
        signs = SharedArray{Float32}((n_epochs, max(1, n_perms)))
    end
    
    if !isnothing(perm_matrix)
        if size(perm_matrix, 1) != n_epochs
            error("perm_matrix must have exactly $n_epochs rows (one for each epoch).")
        end
        n_perms = size(perm_matrix, 2)
        @info "Using user-provided permutation matrix with $n_perms permutations."
        perm_indices .= perm_matrix
    elseif !isnothing(permute_crossed)
        sub_sym = permute_crossed[1]
        item_sym = permute_crossed[2]
        @info "Generating $n_perms synchronized crossed permutations for $(sub_sym) and $(item_sym)..."
        perm_indices .= generate_permutation_matrix(meta_df, item_sym; sync_col=sub_sym, n_perms=n_perms, type=:synchronized, rng=rng)
    elseif !isnothing(permute_block)
        @info "Generating $n_perms within-block permutations for $(permute_block)..."
        perm_indices .= generate_permutation_matrix(meta_df, permute_block; n_perms=n_perms, type=:within, rng=rng)
    else
        for p in 1:max(1, n_perms)
            signs[:, p] .= rand(rng, Float32[-1.0, 1.0], n_epochs)
        end
    end
    
    if use_clusters
        t_perm_matrix = SharedArray{Float32}((n_channels, n_timepoints, max(1, n_perms), n_coefs))
    else
        worker_max_t = SharedArray{Float32}((maximum(workers()), max(1, n_perms), n_coefs))
    end
    
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    perm_ftol_rel = PERMUTATION_FTOL_REL
    
    prog = Progress(length(tasks); dt=1, desc="Fitting models...")
    prog_channel = RemoteChannel(()->Channel{Bool}(n_channels * n_timepoints), 1)
    
    @sync begin
    @async while take!(prog_channel)
        next!(prog)
    end
    @sync @distributed for (c_idx, t_idx) in tasks
        Logging.with_logger(Logging.NullLogger()) do
        buffers = get!(task_local_storage(), :lmm_buffers) do
            LinearAlgebra.BLAS.set_num_threads(1)
            m_thread = Logging.with_logger(Logging.NullLogger()) do
                LinearMixedModel(f, meta_df)
            end
            # modelmatrix(m) is column-PIVOTED in rank-deficient designs, while m.beta,
            # stderror! and coefnames use the ORIGINAL column order. Align X with them.
            X_mat = modelmatrix(m_thread)[:, invperm(MixedModels.pivot(m_thread))]
            y_true = zeros(n_epochs)
            coef_buf = zeros(n_coefs)
            se_buf = zeros(n_coefs)

            fitted_buf = zeros(n_epochs)
            nuisance_fitted = zeros(n_epochs, n_coefs)
            partial_resid_fl = zeros(n_epochs, n_coefs)
            (m_thread, X_mat, y_true, coef_buf, se_buf, fitted_buf, nuisance_fitted, partial_resid_fl)
        end
        (m_thread, X_mat, y_true, coef_buf, se_buf, fitted_buf, nuisance_fitted, partial_resid_fl) = buffers
        
        y_true .= @view eeg_tensor[c_idx, t_idx, :]
        
        # The observed-data fit uses full precision; ftol_rel is relaxed to PERMUTATION_FTOL_REL
        # only for the permutation refits below.
        m_thread.optsum.ftol_rel = 1e-12
        m_thread.optsum.maxfeval = 1000
        
        refit!(m_thread, y_true; progress=false)
        

        
        singular_fits[c_idx, t_idx] = issingular(m_thread)
        
        MixedModels.fixef!(coef_buf, m_thread)
        MixedModels.stderror!(se_buf, m_thread)
        
        beta_matrix[c_idx, t_idx, :] .= coef_buf
        se_matrix[c_idx, t_idx, :] .= se_buf
        
        for coef_idx in 1:n_coefs
            z = coef_buf[coef_idx] / se_buf[coef_idx]
            t_matrix[c_idx, t_idx, coef_idx] = z
            p_matrix[c_idx, t_idx, coef_idx] = 2.0 * ccdf(Normal(), abs(z))
        end
        
        # Freedman-Lane: compute partial residuals for each coefficient.
        # For coefficient j: nuisance_j = ŷ - X[:,j]*β[j]  (everything EXCEPT tested effect)
        #                    partial_resid_j = y - nuisance_j = ε̂ + X[:,j]*β[j]
        # During permutation, only partial_resid_j is shuffled; nuisance stays fixed.
        # This correctly tests each predictor while holding others constant.
        if n_perms > 0
            # Compute marginal fitted values (fixed effects only). 
            # This ensures random effects are left in the permutable residuals.
            mul!(fitted_buf, X_mat, coef_buf)
            for j in tested_coef_indices
                @inbounds @simd for i in 1:n_epochs
                    nuisance_fitted[i, j] = fitted_buf[i] - X_mat[i, j] * coef_buf[j]
                    partial_resid_fl[i, j] = y_true[i] - nuisance_fitted[i, j]
                end
            end
            
            m_y = m_thread.y

            # Permutation refits use the relaxed ftol_rel (see PERMUTATION_FTOL_REL).
            m_thread.optsum.ftol_rel = perm_ftol_rel

            
            for coef_idx in tested_coef_indices
                nuisance_col = @view nuisance_fitted[:, coef_idx]
                resid_col = @view partial_resid_fl[:, coef_idx]

                    
                    for perm_idx in 1:n_perms
                        # Freedman-Lane: write permuted partial residuals directly into m_y
                        if is_permutation
                            perm_vec = @view perm_indices[:, perm_idx]
                            @inbounds for i in 1:n_epochs
                                m_y[i] = nuisance_col[i] + resid_col[perm_vec[i]]
                            end
                        else
                            sign_vec = @view signs[:, perm_idx]
                            @inbounds @simd for i in 1:n_epochs
                                m_y[i] = nuisance_col[i] + resid_col[i] * sign_vec[i]
                            end
                        end
                        
                        # Every permutation is a full refit: all variance components are re-estimated.
                        refit!(m_thread, m_y; progress=false)
                        MixedModels.stderror!(se_buf, m_thread)
                        MixedModels.fixef!(coef_buf, m_thread)
                        z = coef_buf[coef_idx] / se_buf[coef_idx]
                        
                        # Store only the tested coefficient's t-value for its null distribution
                        if use_clusters
                            t_perm_matrix[c_idx, t_idx, perm_idx, coef_idx] = z
                        else
                            wid = myid()
                            if abs(z) > abs(worker_max_t[wid, perm_idx, coef_idx])
                                worker_max_t[wid, perm_idx, coef_idx] = z
                            end
                        end
                    end
                end
            end
        put!(prog_channel, true)
        end # logger
    end 
    put!(prog_channel, false)
    end # sync
    
    max_t_null = zeros(Float32, n_perms, n_coefs)
    max_cluster_mass_null = zeros(Float32, n_perms, n_coefs)
    
    # Aliased coefficients: MixedModels reports beta = -0.0 (and NaN se); make every
    # output NaN so they cannot be mistaken for estimated zero effects.
    if !isempty(dropped_coefs)
        beta_matrix[:, :, dropped_coefs] .= NaN
        se_matrix[:, :, dropped_coefs] .= NaN
        t_matrix[:, :, dropped_coefs] .= NaN
        p_matrix[:, :, dropped_coefs] .= NaN
    end
    
    if n_perms > 0
        if use_clusters
            spatial_connectivity = isnothing(spatial_connectivity) ? EegFun.sparse(Int[], Int[], Bool[], n_channels, n_channels) : SparseMatrixCSC{Bool}(spatial_connectivity)
            electrode_to_idx = Dict(e => i for (i, e) in enumerate(channel_names))
            
            if use_tfce
                @info "Extracting TFCE Maps..."
                
                for perm_idx in 1:n_perms
                    for coef_idx in tested_coef_indices
                        t_map = Array(t_perm_matrix[:, :, perm_idx, coef_idx])
                        tfce_map = EegFun._compute_tfce(t_map, channel_names, Float64.(time_points), spatial_connectivity, :spatiotemporal; E=tfce_E, H=tfce_H, dh=tfce_dh)
                        
                        max_tfce = maximum(abs, tfce_map; init=0.0)
                        max_cluster_mass_null[perm_idx, coef_idx] = max_tfce
                    end
                end
            else
                @info "Extracting Spatio-Temporal Clusters (threshold = $cluster_threshold)..."
                for perm_idx in 1:n_perms
                    for coef_idx in tested_coef_indices
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
                    for coef_idx in tested_coef_indices
                        if abs(worker_max_t[wid, perm_idx, coef_idx]) > max_t_null[perm_idx, coef_idx]
                            max_t_null[perm_idx, coef_idx] = abs(worker_max_t[wid, perm_idx, coef_idx])
                        end
                    end
                end
            end
        end
    end
    
    p_uncorr = Array(p_matrix)
    p_corrected = if n_perms > 0 && !use_clusters
        p_corr = similar(p_uncorr)
        fill!(p_corr, 1.0)
        t_arr = Array(t_matrix)
        for c in tested_coef_indices
            null_col = max_t_null[:, c]
            if any(!=(0), null_col)
                null_sorted = sort(null_col)
                for t_idx in 1:n_timepoints, c_idx in 1:n_channels
                    obs_t = abs(t_arr[c_idx, t_idx, c])
                    n_exceeding = n_perms - searchsortedlast(null_sorted, obs_t - eps(obs_t))
                    p_corr[c_idx, t_idx, c] = (n_exceeding + 1) / (n_perms + 1)
                end
            end
        end
        p_corr[:, :, dropped_coefs] .= NaN
        p_corr
    else
        nothing
    end
    p_main = !isnothing(p_corrected) ? p_corrected : p_uncorr

    return (
        coef_names = coef_names,
        channels = channel_names,
        time = time_points,
        beta = Array(beta_matrix),
        se = Array(se_matrix),
        t = Array(t_matrix),
        p = p_main,
        p_uncorrected = p_uncorr,
        p_corrected = p_corrected,
        max_t_null = max_t_null,
        max_cluster_mass_null = max_cluster_mass_null,
        singular_fits = Array(singular_fits)
    )
end

"""
    generate_permutation_matrix(df::DataFrame, block_col::Symbol; sync_col::Union{Symbol, Nothing}=nothing, n_perms::Int=1000, type::Symbol=:within, rng::AbstractRNG=Random.GLOBAL_RNG)

Helper tool for end-users to generate mathematically valid permutation matrices for `fit_mass_lmm`.
Returns an `N × n_perms` matrix of integer indices.

Types:
- `:within`: Independent shuffling within blocks (e.g. within Subjects). Standard for uncrossed designs.
- `:synchronized`: Shuffles the `block_col` (e.g. Items) and perfectly synchronizes that shuffle across `sync_col` (e.g. Subjects). 
  Requires perfectly balanced crossed designs where every subject sees every item.
"""
function generate_permutation_matrix(df::DataFrame, block_col::Symbol; sync_col::Union{Symbol, Nothing}=nothing, n_perms::Int=1000, type::Symbol=:within, rng::AbstractRNG=Random.GLOBAL_RNG)
    N = nrow(df)
    actual_perms = max(1, n_perms)
    perm_matrix = zeros(Int, N, actual_perms)
    blocks = df[!, block_col]
    unique_blocks = unique(blocks)
    
    if type == :within
        block_idxs = [findall(==(b), blocks) for b in unique_blocks]
        for p in 1:actual_perms
            idx = collect(1:N)
            for b_idx in block_idxs
                idx[b_idx] .= shuffle(rng, b_idx)
            end
            perm_matrix[:, p] = idx
        end
    elseif type == :synchronized
        if isnothing(sync_col)
            error("Synchronized shuffling requires a `sync_col` to synchronize across (e.g. sync_col=:Subject, block_col=:Item).")
        end
        
        syncs = df[!, sync_col]
        unique_blocks = unique(blocks)
        n_blocks = length(unique_blocks)
        
        # Build lookup table: (Sync, Block) -> Vector of RowIndices
        lookup = Dict{Tuple{Any, Any}, Vector{Int}}()
        for i in 1:N
            key = (syncs[i], blocks[i])
            if !haskey(lookup, key)
                lookup[key] = Int[]
            end
            push!(lookup[key], i)
        end
        
        for p in 1:actual_perms
            shuffled_blocks = shuffle(rng, unique_blocks)
            block_map = Dict(unique_blocks[i] => shuffled_blocks[i] for i in 1:n_blocks)
            
            cell_counters = Dict{Tuple{Any, Any}, Int}()
            idx = zeros(Int, N)
            for i in 1:N
                sync_val = syncs[i]
                old_block = blocks[i]
                new_block = block_map[old_block]
                
                key = (sync_val, new_block)
                if haskey(lookup, key)
                    count = get(cell_counters, key, 0) + 1
                    if count <= length(lookup[key])
                        idx[i] = lookup[key][count]
                        cell_counters[key] = count
                    else
                        error("Synchronized shuffling failed! Subject '$sync_val' has fewer trials for Block '$new_block' than required. The design is unbalanced.")
                    end
                else
                    error("Synchronized shuffling failed! Synchronizing `$block_col` across `$sync_col` requires a fully balanced crossed design. Row $i (Sync: $sync_val, Block: $old_block) mapped to new Block '$new_block', but Subject '$sync_val' never saw Item '$new_block'. For unbalanced designs, please supply a custom `perm_matrix`.")
                end
            end
            perm_matrix[:, p] = idx
        end
    else
        error("Unknown permutation type: $type")
    end
    
    return perm_matrix
end

end
