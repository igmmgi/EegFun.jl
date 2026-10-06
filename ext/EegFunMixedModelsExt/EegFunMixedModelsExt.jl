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

import EegFun: fit_mass_lmm

# Relative objective tolerance for the optimizer in permutation refits 
const PERMUTATION_FTOL_REL = 1e-8

# ============================================================================
# 1. COEFFICIENT & INTERCEPT RESOLUTION HELPERS
# ============================================================================

function _add_intercept!(idxs, coef_names, test_intercept)
    if test_intercept
        intercept_idx = findfirst(==("(Intercept)"), coef_names)
        if !isnothing(intercept_idx)
            push!(idxs, intercept_idx)
        end
    end
    return sort!(unique!(idxs))
end

function _resolve_tested_coefs(coef_names, tested_coefs, test_intercept)
    n = length(coef_names)
    
    if tested_coefs === nothing || tested_coefs === :effects || tested_coefs === :auto
        if test_intercept
            return Vector(1:n)
        end
        intercept_idx = findfirst(==("(Intercept)"), coef_names)
        effects = isnothing(intercept_idx) ? Vector(1:n) : filter(!=(intercept_idx), 1:n)
        return isempty(effects) ? [intercept_idx] : effects
    elseif tested_coefs === :all
        return Vector(1:n)
    end
    
    targets = tested_coefs isa Union{AbstractVector, Tuple} ? tested_coefs : (tested_coefs,)
    idxs = Int[]
    
    for t in targets
        if t isa Integer
            (1 <= t <= n) || error("Index $t out of bounds (1:$n)")
            push!(idxs, t)
        elseif t isa Union{AbstractString, Symbol}
            target_str = lowercase(String(t))
            matching = findall(x -> occursin(target_str, lowercase(x)), coef_names)
            isempty(matching) && error("No model coefficient matches '$t'. Available: $coef_names")
            append!(idxs, matching)
        else
            error("Unsupported type for tested_coefs: $(typeof(t))")
        end
    end
    
    return _add_intercept!(idxs, coef_names, test_intercept)
end

# ============================================================================
# 2. HIGH-LEVEL WRAPPER FOR EpochData
# ============================================================================

function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; n_perms=0, kwargs...)
    dfs = epochs.data
    n_epochs = length(dfs)
    n_epochs == 0 && error("EpochData is empty.")
    
    first_df = dfs[1]
    all_channels = intersect(propertynames(first_df), epochs.layout.data.label)
    meta_cols = setdiff(propertynames(first_df), all_channels)
    
    n_timepoints = size(first_df, 1)
    n_channels = length(all_channels)
    
    meta_df = reduce(vcat, [df[1:1, meta_cols] for df in dfs])
    
    eeg_data = zeros(Float64, (n_channels, n_timepoints, n_epochs))
    for (c_idx, c_name) in enumerate(all_channels)
        for (i, df) in enumerate(dfs)
            eeg_data[c_idx, :, i] .= df[!, c_name]
        end
    end
    
    kwdict = Dict{Symbol, Any}(kwargs)
    spatial_conn = get(kwdict, :spatial_connectivity, nothing)
    use_clusters = get(kwdict, :use_clusters, false)
    if isnothing(spatial_conn) && use_clusters && !isnothing(epochs.layout)
        spatial_conn = EegFun._build_connectivity_matrix(all_channels, epochs.layout, :spatiotemporal)
    end
    
    res = fit_mass_lmm(eeg_data, meta_df, f; 
        channel_names = all_channels,
        time_points = first_df.time,
        dim_order = (:channels, :timepoints, :epochs),
        n_perms = n_perms,
        kwargs...,
        spatial_connectivity = spatial_conn
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

# ============================================================================
# 3. LOW-LEVEL PERFORMANCE KERNELS & IN-PLACE MATH
# ============================================================================

@inline _alloc(T, dims...) = nprocs() > 1 ? SharedArray{T}(dims) : zeros(T, dims)

@inline function _apply_permutation!(m_y, y_src, perm_vec)
    @inbounds for i in eachindex(m_y)
        m_y[i] = y_src[perm_vec[i]]
    end
end

@inline function _apply_sign_flip!(m_y, y_src, sign_vec)
    @inbounds @simd for i in eachindex(m_y)
        m_y[i] = y_src[i] * sign_vec[i]
    end
end

@inline function _apply_fl_perm!(m_y, nuisance, resid, perm_vec)
    @inbounds for i in eachindex(m_y)
        m_y[i] = nuisance[i] + resid[perm_vec[i]]
    end
end

@inline function _apply_fl_flip!(m_y, nuisance, resid, sign_vec)
    @inbounds @simd for i in eachindex(m_y)
        m_y[i] = nuisance[i] + resid[i] * sign_vec[i]
    end
end

@inline function _approx_ols_z(m_y, X_pinv, X_mat, XtX_inv_diag, coef_idx, beta_buf)
    n_epochs, n_coefs = size(X_mat)
    mul!(beta_buf, X_pinv, m_y)
    sigma2 = 0.0
    @inbounds @simd for i in 1:n_epochs
        pred = 0.0
        for j in 1:n_coefs
            pred += X_mat[i, j] * beta_buf[j]
        end
        sigma2 += abs2(m_y[i] - pred)
    end
    sigma2 /= (n_epochs - n_coefs)
    se_perm_coef = sqrt(XtX_inv_diag[coef_idx] * sigma2)
    return beta_buf[coef_idx] / se_perm_coef
end

# ============================================================================
# 4. WORKER SCRATCHPAD CONTAINER
# ============================================================================

struct LmmWorkerScratchpad{M<:LinearMixedModel}
    m::M
    X_mat::Matrix{Float64}
    y_true::Vector{Float64}
    coef_buf::Vector{Float64}
    se_buf::Vector{Float64}
    fitted_buf::Vector{Float64}
    nuisance::Vector{Float64}
    resid::Vector{Float64}
    X_pinv::Matrix{Float64}
    XtX_inv_diag::Vector{Float64}
    Zb_buf::Vector{Float64}
    y_marg::Vector{Float64}
end

function _create_scratchpad(f, meta_df, n_epochs, n_coefs, approximate_ols)
    LinearAlgebra.BLAS.set_num_threads(1)
    m_thread = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df)
    end
    # Align X with original column order of m.beta and coefnames
    X_mat = modelmatrix(m_thread)[:, invperm(MixedModels.pivot(m_thread))]
    
    return LmmWorkerScratchpad(
        m_thread,
        X_mat,
        zeros(Float64, n_epochs),
        zeros(Float64, n_coefs),
        zeros(Float64, n_coefs),
        zeros(Float64, n_epochs),
        zeros(Float64, n_epochs),
        zeros(Float64, n_epochs),
        approximate_ols ? pinv(X_mat) : zeros(Float64, 0, 0),
        approximate_ols ? diag(inv(X_mat' * X_mat)) : zeros(Float64, 0),
        approximate_ols ? zeros(Float64, n_epochs) : zeros(Float64, 0),
        approximate_ols ? zeros(Float64, n_epochs) : zeros(Float64, 0)
    )
end

@inline function _fit_pixel!(
    scratch, eeg_tensor, c_idx, t_idx,
    beta_matrix, se_matrix, t_matrix, p_matrix, singular_fits, n_coefs
)
    scratch.y_true .= @view eeg_tensor[c_idx, t_idx, :]
    scratch.m.optsum.ftol_rel = 1e-12
    scratch.m.optsum.maxfeval = 1000
    
    refit!(scratch.m, scratch.y_true; progress=false)
    singular_fits[c_idx, t_idx] = issingular(scratch.m)
    
    scratch.coef_buf .= scratch.m.beta
    MixedModels.stderror!(scratch.se_buf, scratch.m)
    
    @inbounds for coef_idx in 1:n_coefs
        b = scratch.coef_buf[coef_idx]
        se = scratch.se_buf[coef_idx]
        z = b / se
        beta_matrix[c_idx, t_idx, coef_idx] = b
        se_matrix[c_idx, t_idx, coef_idx] = se
        t_matrix[c_idx, t_idx, coef_idx] = z
        p_matrix[c_idx, t_idx, coef_idx] = 2.0 * ccdf(Normal(), abs(z))
    end
end

# ============================================================================
# 5. MODEL DESIGN & PERMUTATION PIPELINE SETUP
# ============================================================================

function _setup_lmm_design(f, meta_df, n_epochs, tested_coefs, test_intercept, n_perms, has_perm_matrix)
    meta_df_init = copy(meta_df)
    lhs_sym = Symbol(f.lhs)
    meta_df_init[!, lhs_sym] = randn(Random.MersenneTwister(0), n_epochs)
    
    @info "Compiling MixedModel design matrix..."
    m_initial = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df_init)
    end
    m_initial.optsum.maxfeval = 1
    fit!(m_initial; progress=false)
    
    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)
    
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
        if n_perms > 0 || has_perm_matrix
            @info "Skipping permutation testing for aliased coefficient(s): $(join(coef_names[untestable], ", "))"
        end
    end
    if n_perms > 0
        @info "Permutation testing will be computed for $(length(tested_coef_indices)) of $n_coefs coefficients: $(join(coef_names[tested_coef_indices], ", ")) (test_intercept=$test_intercept)"
    end
    
    return coef_names, n_coefs, dropped_coefs, tested_coef_indices, meta_df_init
end

function _setup_permutations(meta_df, n_epochs, n_perms, perm_matrix, permute_block, permute_crossed, rng)
    is_permutation = !isnothing(perm_matrix) || !isnothing(permute_block) || !isnothing(permute_crossed)
    actual_perms = max(1, n_perms)
    
    if is_permutation
        perm_indices = _alloc(Int, n_epochs, actual_perms)
        signs = _alloc(Float32, 0, 0)
        
        if !isnothing(perm_matrix)
            size(perm_matrix, 1) == n_epochs || error("perm_matrix must have exactly $n_epochs rows (one for each epoch).")
            actual_perms = size(perm_matrix, 2)
            @info "Using user-provided permutation matrix with $actual_perms permutations."
            perm_indices .= perm_matrix
        elseif !isnothing(permute_crossed)
            sub_sym, item_sym = permute_crossed[1], permute_crossed[2]
            @info "Generating $n_perms synchronized crossed permutations for $sub_sym and $item_sym..."
            perm_indices .= generate_permutation_matrix(meta_df, item_sym; sync_col=sub_sym, n_perms=n_perms, type=:synchronized, rng=rng)
        elseif !isnothing(permute_block)
            @info "Generating $n_perms within-block permutations for $permute_block..."
            perm_indices .= generate_permutation_matrix(meta_df, permute_block; n_perms=n_perms, type=:within, rng=rng)
        end
    else
        perm_indices = _alloc(Int, 0, 0)
        signs = _alloc(Float32, n_epochs, actual_perms)
        for p in 1:actual_perms
            signs[:, p] .= rand(rng, Float32[-1.0, 1.0], n_epochs)
        end
    end
    
    return perm_indices, signs, is_permutation, actual_perms
end

# ============================================================================
# 6. NULL DISTRIBUTION EXTRACTION & STATISTICAL INFERENCE
# ============================================================================

function _build_cluster_stat_map(t_map, channel_names, time_points, spatial_conn, electrode_to_idx, cluster_threshold)
    n_channels, n_timepoints = size(t_map)
    pos_clusters, neg_clusters = EegFun._find_clusters(
        t_map .> cluster_threshold, t_map .< -cluster_threshold,
        channel_names, Float64.(time_points), spatial_conn, :spatiotemporal
    )
    pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
    neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
    
    smap = zeros(Float64, n_channels, n_timepoints)
    for (i, cluster) in enumerate(pos_clusters)
        stat = pos_stats[i]
        for (e_idx, t_idx) in cluster.members
            smap[e_idx, t_idx] = max(smap[e_idx, t_idx], stat)
        end
    end
    for (i, cluster) in enumerate(neg_clusters)
        stat = abs(neg_stats[i])
        for (e_idx, t_idx) in cluster.members
            smap[e_idx, t_idx] = max(smap[e_idx, t_idx], stat)
        end
    end
    return smap
end

function _extract_null_distributions!(
    max_cluster_mass_null, max_t_null, t_perm_matrix, worker_max_t,
    n_perms, tested_coef_indices, channel_names, time_points, spatial_conn, electrode_to_idx,
    use_clusters, use_tfce, tfce_E, tfce_H, tfce_dh, cluster_threshold
)
    if !use_clusters
        w_list = workers()
        for coef_idx in tested_coef_indices
            @inbounds for perm_idx in 1:n_perms
                m_val = 0.0f0
                for wid in w_list
                    val = abs(worker_max_t[wid, perm_idx, coef_idx])
                    if val > m_val
                        m_val = val
                    end
                end
                max_t_null[perm_idx, coef_idx] = m_val
            end
        end
        return
    end

    if use_tfce
        @info "Extracting TFCE Maps..."
        for perm_idx in 1:n_perms
            for coef_idx in tested_coef_indices
                t_map = Array(t_perm_matrix[:, :, perm_idx, coef_idx])
                tfce_map = EegFun._compute_tfce(t_map, channel_names, Float64.(time_points), spatial_conn, :spatiotemporal; E=tfce_E, H=tfce_H, dh=tfce_dh)
                max_cluster_mass_null[perm_idx, coef_idx] = maximum(abs, tfce_map; init=0.0)
            end
        end
    else
        @info "Extracting Spatio-Temporal Clusters (threshold = $cluster_threshold)..."
        for perm_idx in 1:n_perms
            for coef_idx in tested_coef_indices
                t_map = Array(t_perm_matrix[:, :, perm_idx, coef_idx])
                pos_clusters, neg_clusters = EegFun._find_clusters(
                    t_map .> cluster_threshold, t_map .< -cluster_threshold,
                    channel_names, Float64.(time_points), spatial_conn, :spatiotemporal
                )
                pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map, electrode_to_idx; return_clusters=false)
                neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map, electrode_to_idx; return_clusters=false)
                
                max_pos = maximum(pos_stats, init=0.0)
                max_neg = maximum(abs, neg_stats, init=0.0)
                max_cluster_mass_null[perm_idx, coef_idx] = max(max_pos, max_neg)
            end
        end
    end
end

function _compute_p_corrected(
    t_matrix, max_t_null, max_cluster_mass_null,
    n_perms, tested_coef_indices, channel_names, time_points, spatial_conn, electrode_to_idx,
    use_clusters, use_tfce, tfce_E, tfce_H, tfce_dh, cluster_threshold, dropped_coefs
)
    n_channels = length(channel_names)
    n_timepoints = length(time_points)
    p_corr = ones(Float64, n_channels, n_timepoints, size(t_matrix, 3))
    t_arr = Array(t_matrix)

    for c in tested_coef_indices
        null_col = use_clusters ? @view(max_cluster_mass_null[:, c]) : @view(max_t_null[:, c])
        if !any(!=(0), null_col)
            continue
        end
        null_sorted = sort(null_col)
        
        stat_map = if !use_clusters
            nothing
        elseif use_tfce
            abs.(EegFun._compute_tfce(t_arr[:, :, c], channel_names, Float64.(time_points), spatial_conn, :spatiotemporal; E=tfce_E, H=tfce_H, dh=tfce_dh))
        else
            _build_cluster_stat_map(t_arr[:, :, c], channel_names, time_points, spatial_conn, electrode_to_idx, cluster_threshold)
        end

        for t_idx in 1:n_timepoints, c_idx in 1:n_channels
            obs_val = isnothing(stat_map) ? abs(t_arr[c_idx, t_idx, c]) : stat_map[c_idx, t_idx]
            n_exceeding = n_perms - searchsortedlast(null_sorted, obs_val - eps(obs_val))
            p_corr[c_idx, t_idx, c] = (n_exceeding + 1) / (n_perms + 1)
        end
    end
    
    p_corr[:, :, dropped_coefs] .= NaN
    return p_corr
end

# ============================================================================
# 7. PRIMARY MASS-UNIVARIATE LMM FITTER
# ============================================================================

function fit_mass_lmm(eeg_data, meta_df, f;
    n_perms = 1000,
    tested_coefs = nothing,
    test_intercept = false,
    permute_block = nothing,
    permute_crossed = nothing,
    perm_matrix = nothing,
    use_clusters = false,
    cluster_threshold = 2.0,
    spatial_connectivity = nothing,
    channel_names = nothing,
    time_points = nothing,
    dim_order = (:channels, :timepoints, :epochs),
    use_tfce = false,
    approximate_ols = false,
    tfce_E = 0.5,
    tfce_H = 2.0,
    tfce_dh = 0.1,
    rng = Random.GLOBAL_RNG
)
    # 1. Dimension ordering and validation
    if dim_order != (:channels, :timepoints, :epochs)
        target = (:channels, :timepoints, :epochs)
        perm = [findfirst(==(t), dim_order) for t in target]
        any(isnothing, perm) && error("dim_order must contain :channels, :timepoints, and :epochs")
        eeg_data = permutedims(eeg_data, tuple(perm...))
    end
    n_channels, n_timepoints, n_epochs = size(eeg_data)
    nrow(meta_df) == n_epochs || error("Number of rows in meta_df ($(nrow(meta_df))) must match number of epochs in eeg_data ($n_epochs).")
    
    channel_names = isnothing(channel_names) ? [Symbol("Ch$i") for i in 1:n_channels] : channel_names
    time_points = isnothing(time_points) ? Vector(1:n_timepoints) : time_points

    # 2. Setup model design, rank deficiency & tested coefficients
    has_perm_matrix = !isnothing(perm_matrix)
    coef_names, n_coefs, dropped_coefs, tested_coef_indices, meta_df_init = _setup_lmm_design(
        f, meta_df, n_epochs, tested_coefs, test_intercept, n_perms, has_perm_matrix
    )

    # 3. Setup permutation matrices / signs
    perm_indices, signs, is_permutation, n_perms = _setup_permutations(
        meta_df, n_epochs, n_perms, perm_matrix, permute_block, permute_crossed, rng
    )

    # 4. Allocate tensors & null distribution buffers
    eeg_tensor = (nprocs() > 1 && !(eeg_data isa SharedArray)) ? (SharedArray{Float64}(size(eeg_data)) .= eeg_data) : eeg_data
    
    beta_matrix = _alloc(Float64, n_channels, n_timepoints, n_coefs)
    se_matrix = _alloc(Float64, n_channels, n_timepoints, n_coefs)
    t_matrix = _alloc(Float64, n_channels, n_timepoints, n_coefs)
    p_matrix = _alloc(Float64, n_channels, n_timepoints, n_coefs)
    singular_fits = _alloc(Bool, n_channels, n_timepoints)
    
    t_perm_matrix = use_clusters ? _alloc(Float32, n_channels, n_timepoints, max(1, n_perms), n_coefs) : _alloc(Float32, 0, 0, 0, 0)
    worker_max_t = !use_clusters ? _alloc(Float32, maximum(workers()), max(1, n_perms), n_coefs) : _alloc(Float32, 0, 0, 0)

    @info "Fitting $(n_channels * n_timepoints) Mixed Models across $n_epochs epochs (n_perms=$n_perms) on $(nprocs()) processes..."

    # 5. Distributed mass-univariate execution
    tasks = [(c_idx, t_idx) for c_idx in 1:n_channels for t_idx in 1:n_timepoints]
    perm_ftol_rel = PERMUTATION_FTOL_REL
    
    prog = Progress(length(tasks); dt=1, desc="Fitting models...")
    prog_channel = RemoteChannel(() -> Channel{Bool}(n_channels * n_timepoints), 1)

    @sync begin
        @async while take!(prog_channel)
            next!(prog)
        end
        
        @sync @distributed for (c_idx, t_idx) in tasks
            Logging.with_logger(Logging.NullLogger()) do
                scratch = get!(task_local_storage(), :lmm_buffers) do
                    _create_scratchpad(f, meta_df_init, n_epochs, n_coefs, approximate_ols)
                end
                
                # Fit observed data
                _fit_pixel!(scratch, eeg_tensor, c_idx, t_idx, beta_matrix, se_matrix, t_matrix, p_matrix, singular_fits, n_coefs)
                
                # Permutation refits (Freedman-Lane)
                if n_perms > 0
                    mul!(scratch.fitted_buf, scratch.X_mat, scratch.coef_buf)
                    if approximate_ols
                        scratch.Zb_buf .= fitted(scratch.m) .- scratch.fitted_buf
                        scratch.y_marg .= scratch.y_true .- scratch.Zb_buf
                    end
                    
                    scratch.m.optsum.ftol_rel = perm_ftol_rel
                    m_y = scratch.m.y
                    wid = !use_clusters ? myid() : 0
                    
                    for coef_idx in tested_coef_indices
                        if !approximate_ols
                            @inbounds @simd for i in 1:n_epochs
                                scratch.nuisance[i] = scratch.fitted_buf[i] - scratch.X_mat[i, coef_idx] * scratch.coef_buf[coef_idx]
                                scratch.resid[i] = scratch.y_true[i] - scratch.nuisance[i]
                            end
                        end
                        
                        for perm_idx in 1:n_perms
                            if approximate_ols
                                if is_permutation
                                    _apply_permutation!(m_y, scratch.y_marg, @view perm_indices[:, perm_idx])
                                else
                                    _apply_sign_flip!(m_y, scratch.y_marg, @view signs[:, perm_idx])
                                end
                                z = _approx_ols_z(m_y, scratch.X_pinv, scratch.X_mat, scratch.XtX_inv_diag, coef_idx, scratch.coef_buf)
                            else
                                if is_permutation
                                    _apply_fl_perm!(m_y, scratch.nuisance, scratch.resid, @view perm_indices[:, perm_idx])
                                else
                                    _apply_fl_flip!(m_y, scratch.nuisance, scratch.resid, @view signs[:, perm_idx])
                                end
                                
                                refit!(scratch.m, m_y; progress=false)
                                scratch.coef_buf .= scratch.m.beta
                                MixedModels.stderror!(scratch.se_buf, scratch.m)
                                z = scratch.coef_buf[coef_idx] / scratch.se_buf[coef_idx]
                            end
                            
                            if use_clusters
                                t_perm_matrix[c_idx, t_idx, perm_idx, coef_idx] = z
                            else
                                if abs(z) > abs(worker_max_t[wid, perm_idx, coef_idx])
                                    worker_max_t[wid, perm_idx, coef_idx] = z
                                end
                            end
                        end
                    end
                end
                put!(prog_channel, true)
            end
        end
        put!(prog_channel, false)
    end

    # 6. Post-hoc masking for aliased coefficients
    if !isempty(dropped_coefs)
        beta_matrix[:, :, dropped_coefs] .= NaN
        se_matrix[:, :, dropped_coefs] .= NaN
        t_matrix[:, :, dropped_coefs] .= NaN
        p_matrix[:, :, dropped_coefs] .= NaN
    end

    max_t_null = zeros(Float32, n_perms, n_coefs)
    max_cluster_mass_null = zeros(Float32, n_perms, n_coefs)
    p_uncorr = Array(p_matrix)
    p_corrected = nothing

    # 7. Extract null distribution & compute corrected p-values
    if n_perms > 0
        spatial_conn = isnothing(spatial_connectivity) ? EegFun.sparse(Int[], Int[], Bool[], n_channels, n_channels) : SparseMatrixCSC{Bool}(spatial_connectivity)
        electrode_to_idx = Dict(e => i for (i, e) in enumerate(channel_names))
        
        _extract_null_distributions!(
            max_cluster_mass_null, max_t_null, t_perm_matrix, worker_max_t,
            n_perms, tested_coef_indices, channel_names, time_points, spatial_conn, electrode_to_idx,
            use_clusters, use_tfce, tfce_E, tfce_H, tfce_dh, cluster_threshold
        )
        p_corrected = _compute_p_corrected(
            t_matrix, max_t_null, max_cluster_mass_null,
            n_perms, tested_coef_indices, channel_names, time_points, spatial_conn, electrode_to_idx,
            use_clusters, use_tfce, tfce_E, tfce_H, tfce_dh, cluster_threshold, dropped_coefs
        )
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

end # module
