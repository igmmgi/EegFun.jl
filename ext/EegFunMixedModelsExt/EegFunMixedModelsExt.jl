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

import EegFun: fit_mass_lmm, @minimal_error

# Relative objective tolerance for the (single) observed fit per pixel
const OBSERVED_FTOL_REL = 1e-12

# Minimum number of distinct sign-flip patterns (2^(S-1)); below this the smallest attainable
# p-value exceeds 0.05 and the test cannot reject.
const MIN_SIGN_FLIPS = 20

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

# Match a user-supplied name to coefficient names (case-insensitive). An exact match wins.
# Otherwise the name selects all contrast columns of that term (e.g. `:condition` matches
# "condition: B" and "condition: C"), but never interactions involving it, unless the
# requested name is itself an interaction (contains '&').
function _match_coef_names(coef_names, target::AbstractString)
    t = lowercase(strip(target))
    names_lc = lowercase.(coef_names)
    exact = findall(==(t), names_lc)
    isempty(exact) || return exact
    is_interaction = occursin('&', t)
    return findall(names_lc) do x
        startswith(x, t * ": ") && (is_interaction || !occursin('&', x))
    end
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
            (1 <= t <= n) || @minimal_error "Index $t out of bounds (1:$n)"
            push!(idxs, t)
        elseif t isa Union{AbstractString, Symbol}
            matching = _match_coef_names(coef_names, String(t))
            isempty(matching) && @minimal_error "No model coefficient matches '$t'. Available: $coef_names"
            append!(idxs, matching)
        else
            @minimal_error "Unsupported type for tested_coefs: $(typeof(t))"
        end
    end

    return _add_intercept!(idxs, coef_names, test_intercept)
end

# ============================================================================
# 2. HIGH-LEVEL WRAPPER FOR EpochData
# ============================================================================

@inline function _fill_lmm_epoch_slice!(eeg_data::Array{Float64,3}, c_idx::Int, ep_idx::Int, col::Vector{Float64}, n_samp::Int)
    @inbounds @simd for t = 1:n_samp
        eeg_data[c_idx, t, ep_idx] = col[t]
    end
end
@inline function _fill_lmm_epoch_slice!(eeg_data::Array{Float64,3}, c_idx::Int, ep_idx::Int, col::AbstractVector, n_samp::Int)
    @inbounds for t = 1:n_samp
        eeg_data[c_idx, t, ep_idx] = Float64(col[t])
    end
end

function fit_mass_lmm(epochs::EegFun.EpochData, f::FormulaTerm; correction::Symbol = :maxt, spatial_connectivity = nothing, kwargs...)
    isempty(epochs) && @minimal_error "EpochData is empty."

    dfs = epochs.data
    n_epochs = length(dfs)

    all_channels = EegFun.channel_labels(epochs)
    n_channels = length(all_channels)
    n_timepoints = EegFun.n_samples(epochs)
    meta_cols = setdiff(EegFun.all_labels(epochs), all_channels)

    meta_df = DataFrame([col => [df[1, col] for df in dfs] for col in meta_cols])

    eeg_data = Array{Float64, 3}(undef, n_channels, n_timepoints, n_epochs)
    for (i, df) in enumerate(dfs)
        for (c_idx, c_name) in enumerate(all_channels)
            _fill_lmm_epoch_slice!(eeg_data, c_idx, i, df[!, c_name], n_timepoints)
        end
    end

    if isnothing(spatial_connectivity) && correction in (:cluster, :tfce) && !isnothing(epochs.layout)
        spatial_connectivity = EegFun._build_connectivity_matrix(all_channels, epochs.layout, :spatiotemporal)
    end

    tvec = EegFun.time_vector(epochs)
    time_pts = isempty(tvec) ? Vector(1:n_timepoints) : tvec

    return _fit_mass_lmm(eeg_data, meta_df, f, epochs;
        channel_names = all_channels,
        time_points = time_pts,
        correction = correction,
        spatial_connectivity = spatial_connectivity,
        kwargs...
    )
end

# ============================================================================
# 3. EXCHANGEABILITY BLOCKS (whole-block sign flipping)
# ============================================================================
#
# Under H0, whole-block sign flipping of the reduced-model (Freedman–Lane) residuals preserves
# their joint distribution when blocks are independent and residuals are symmetric. This targets
# the population null (fixed effect = 0 averaged over blocks), handles within- and between-block
# predictors, and unequal block sizes. It requires every random-effects grouping factor to be
# nested within the flipping block.

_block_ids(col) = (lv = unique(col); idx = Dict(v => i for (i, v) in enumerate(lv)); [idx[v] for v in col])

# Every level of every random-effects term must occur in exactly one block
function _is_nested(m::LinearMixedModel, block_ids::AbstractVector{<:Integer})
    for A in m.reterms
        owner = zeros(Int, length(A.levels))
        for (obs, l) in enumerate(A.refs)
            b = block_ids[obs]
            owner[l] == 0 && (owner[l] = b)
            owner[l] == b || return false
        end
    end
    return true
end

function _resolve_flip_by(m::LinearMixedModel, meta_df::DataFrame, flip_by)
    cols = propertynames(meta_df)
    if isnothing(flip_by)
        # Outermost grouping factor: the coarsest plain-column factor within which all others nest
        candidates = Symbol[MixedModels.fname(A) for A in m.reterms if MixedModels.fname(A) in cols]
        sort!(candidates; by = c -> length(unique(meta_df[!, c])))
        for c in candidates
            ids = _block_ids(meta_df[!, c])
            _is_nested(m, ids) && return c, ids
        end
        @minimal_error "The random-effects grouping factors ($(join(MixedModels.fname.(m.reterms), ", "))) are crossed " *
                       "(e.g. subject × item): no single factor contains all others, so no exact permutation scheme exists. " *
                       "Use `n_perms = 0` for Wald-only inference."
    end
    flip_by = Symbol(flip_by)
    flip_by in cols || @minimal_error "flip_by = :$flip_by is not a column of the metadata. Available: $(cols)"
    ids = _block_ids(meta_df[!, flip_by])
    _is_nested(m, ids) || @minimal_error "Random-effects grouping factors ($(join(MixedModels.fname.(m.reterms), ", "))) " *
                                         "are not all nested within flip_by = :$flip_by. Whole-block sign flipping requires that every " *
                                         "random-effects level belongs to a single block."
    return flip_by, ids
end

# Sign-flip matrix (n_blocks × n_flips). Because |z| is invariant to η → -η, only 2^(S-1) patterns are
# distinct: enumerate them all (including the identity) when that is no more than n_perms.
function _sign_flip_matrix(n_blocks::Int, n_perms::Int, rng)
    if n_blocks - 1 < 62 && 2^(n_blocks - 1) <= n_perms
        n_flips = 2^(n_blocks - 1)
        H = ones(Float64, n_blocks, n_flips)
        for k in 0:(n_flips - 1), s in 2:n_blocks
            ((k >> (s - 2)) & 1) == 1 && (H[s, k + 1] = -1.0)
        end
        return H, true
    end
    return rand(rng, [-1.0, 1.0], n_blocks, n_perms), false
end

# Pixel-independent per-block design pieces
struct BlockDesign
    rows::Vector{Vector{Int}}
    X::Vector{Matrix{Float64}}                    # n_s × p (estimable fixed-effects columns)
    Z::Vector{Matrix{Float64}}                    # n_s × q_s (unscaled random-effects columns)
    XtX::Vector{Matrix{Float64}}
    ZtX::Vector{Matrix{Float64}}
    ZtZ::Vector{Matrix{Float64}}
    lambda_blocks::Vector{Vector{Tuple{Int,Int}}} # per block: (term index, column offset) for each level
end

function BlockDesign(m::LinearMixedModel, X::Matrix{Float64}, block_ids::Vector{Int})
    n_blocks = maximum(block_ids)
    rows = [Int[] for _ in 1:n_blocks]
    for (i, b) in enumerate(block_ids)
        push!(rows[b], i)
    end
    Xs, Zs, lam = Matrix{Float64}[], Matrix{Float64}[], Vector{Tuple{Int,Int}}[]
    for r in rows
        cols, blocks, off = Matrix{Float64}[], Tuple{Int,Int}[], 0
        for (k, A) in enumerate(m.reterms)
            q = size(A.z, 1)
            levs = sort!(unique(view(A.refs, r)))
            Zk = zeros(length(r), q * length(levs))
            for (i, obs) in enumerate(r)
                j = searchsortedfirst(levs, A.refs[obs])
                Zk[i, ((j - 1) * q + 1):(j * q)] .= view(A.z, :, obs)
            end
            for j in eachindex(levs)
                push!(blocks, (k, off + (j - 1) * q))
            end
            off += q * length(levs)
            push!(cols, Zk)
        end
        push!(Xs, X[r, :])
        push!(Zs, reduce(hcat, cols))
        push!(lam, blocks)
    end
    return BlockDesign(rows, Xs, Zs, [x' * x for x in Xs], [z' * x for (z, x) in zip(Zs, Xs)], [z' * z for z in Zs], lam)
end

# Per-block GLS sufficient statistics at the model's current θ, with V_s = I + Z_s Λ_s Λ_s' Z_s':
#   XtVX_s = X_s' V_s⁻¹ X_s,  XtVy_s = X_s' V_s⁻¹ y_s,  ytVy_s = y_s' V_s⁻¹ y_s   (Woodbury)
function _block_gls_stats(design::BlockDesign, lams::Vector{Matrix{Float64}}, y::AbstractVector{Float64})
    n_blocks = length(design.rows)
    XtVX = Vector{Matrix{Float64}}(undef, n_blocks)
    XtVy = Vector{Vector{Float64}}(undef, n_blocks)
    ytVy = Vector{Float64}(undef, n_blocks)
    for s in 1:n_blocks
        ys = view(y, design.rows[s])
        q = size(design.Z[s], 2)
        Λ = zeros(q, q)
        for (k, off) in design.lambda_blocks[s]
            λ = lams[k]
            qk = size(λ, 1)
            Λ[(off + 1):(off + qk), (off + 1):(off + qk)] .= λ
        end
        Xty = design.X[s]' * ys
        Gty = Λ' * (design.Z[s]' * ys)
        GtX = Λ' * design.ZtX[s]
        C = cholesky!(Symmetric(Λ' * design.ZtZ[s] * Λ + I))
        FGX = C \ GtX
        XtVX[s] = design.XtX[s] - GtX' * FGX
        XtVy[s] = Xty - FGX' * Gty
        ytVy[s] = dot(ys, ys) - dot(Gty, C \ Gty)
    end
    return XtVX, XtVy, ytVy
end

# Frozen-θ Freedman–Lane z-statistics for coefficient position `pc` under every sign pattern in H.
# With reduced-model GLS nuisance γ̂, per-block w_s = X_s'V_s⁻¹ r_s and d = Σ r_s'V_s⁻¹ r_s:
#   β̂*_c = (M W η)_c,   pwrss* = d − (Wη)' M (Wη)   (using η_s² = 1),   z* = β̂*_c / sqrt(pwrss*/ν · M_cc)
# At η = 1 this reproduces the LMM's own Wald z exactly.
function _null_z!(out::AbstractVector, XtVX_s, XtVy_s, ytVy_s, pc::Int, H::Matrix{Float64}, ν::Real)
    p, n_blocks = size(first(XtVX_s), 1), length(XtVy_s)
    XtVX, XtVy = sum(XtVX_s), sum(XtVy_s)
    M = inv(cholesky(Symmetric(XtVX)))
    idx = [1:(pc - 1); (pc + 1):p]
    γ = isempty(idx) ? Float64[] : cholesky(Symmetric(XtVX[idx, idx])) \ XtVy[idx]
    W = Matrix{Float64}(undef, p, n_blocks)
    d = 0.0
    for s in 1:n_blocks
        XtVXr = view(XtVX_s[s], :, idx)
        W[:, s] .= XtVy_s[s] .- XtVXr * γ
        d += ytVy_s[s] - 2 * dot(γ, view(XtVy_s[s], idx)) + dot(γ, view(XtVX_s[s], idx, idx) * γ)
    end
    U = W * H
    MU = M * U
    scale = M[pc, pc] / ν
    @inbounds for j in axes(H, 2)
        quad = 0.0
        for i in 1:p
            quad += U[i, j] * MU[i, j]
        end
        out[j] = MU[pc, j] / sqrt(max(d - quad, eps(d)) * scale)
    end
    return out
end

# ============================================================================
# 4. MODEL DESIGN SETUP
# ============================================================================

function _setup_lmm_design(f, meta_df, n_epochs, tested_coefs, test_intercept)
    meta_df_init = copy(meta_df)
    meta_df_init[!, Symbol(f.lhs)] = randn(Random.MersenneTwister(0), n_epochs)

    @info "Compiling MixedModel design matrix..."
    # Coefficient names, rank and pivot are available from the unfitted model
    m_initial = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df_init)
    end
    length(m_initial.y) == n_epochs || @minimal_error "The model uses $(length(m_initial.y)) of $n_epochs epochs " *
                                                      "(missing values in model variables?). Remove incomplete rows first."

    coef_names = coefnames(m_initial)
    n_coefs = length(coef_names)

    rank_X = m_initial.feterm.rank
    dropped_coefs = sort(MixedModels.pivot(m_initial)[(rank_X + 1):end])
    if !isempty(dropped_coefs)
        @warn "Fixed-effects design matrix is rank deficient (rank $rank_X of $n_coefs). " *
              "Aliased coefficient(s) $(join(coef_names[dropped_coefs], ", ")) cannot be estimated: " *
              "their beta, se, z and p values are NaN and they are excluded from permutation testing. " *
              "Common causes: a predictor that is constant in the data, an empty design cell combined " *
              "with an interaction, or exactly collinear predictors."
    end

    tested = setdiff(_resolve_tested_coefs(coef_names, tested_coefs, test_intercept), dropped_coefs)
    return m_initial, coef_names, dropped_coefs, tested, meta_df_init
end

# ============================================================================
# 5. PER-PIXEL WORK (one LMM fit + closed-form sign-flip null)
# ============================================================================

struct LmmScratch{M<:LinearMixedModel}
    m::M
    X::Matrix{Float64}          # estimable fixed-effects columns, coefnames order
    design::Union{BlockDesign,Nothing}
    y::Vector{Float64}
    se::Vector{Float64}
    znull::Vector{Float64}
end

function _create_scratch(f, meta_df_init, est_cols, block_ids, n_flips)
    m = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df_init)
    end
    m.optsum.ftol_rel = OBSERVED_FTOL_REL
    m.optsum.maxfeval = 1000
    X = Matrix{Float64}(modelmatrix(m)[:, invperm(MixedModels.pivot(m))][:, est_cols])
    design = isnothing(block_ids) ? nothing : BlockDesign(m, X, block_ids)
    n_obs = length(m.y)
    return LmmScratch(m, X, design, zeros(n_obs), zeros(length(coefnames(m))), zeros(n_flips))
end

function _process_pixel!(scr::LmmScratch, eeg, c_idx, t_idx, out, ctx)
    scr.y .= @view eeg[c_idx, t_idx, :]
    m = scr.m
    refit!(m, scr.y; REML = ctx.REML, progress = false)
    out.singular[c_idx, t_idx] = issingular(m)
    MixedModels.stderror!(scr.se, m)
    β = m.beta
    @inbounds for k in eachindex(scr.se)
        z = β[k] / scr.se[k]
        out.beta[c_idx, t_idx, k] = β[k]
        out.se[c_idx, t_idx, k] = scr.se[k]
        out.z[c_idx, t_idx, k] = z
        out.p[c_idx, t_idx, k] = 2.0 * ccdf(Normal(), abs(z))
    end

    ctx.n_flips == 0 && return
    lams = Matrix{Float64}[Matrix(A.λ) for A in m.reterms]
    stats = _block_gls_stats(scr.design, lams, scr.y)
    ν = MixedModels.ssqdenom(m)
    for (j, pc) in enumerate(ctx.tested_pos)
        _null_z!(scr.znull, stats..., pc, ctx.H, ν)
        if ctx.store_maps
            @inbounds for p in 1:ctx.n_flips
                out.zperm[c_idx, t_idx, p, j] = scr.znull[p]
            end
        else
            slot = ctx.slot
            @inbounds for p in 1:ctx.n_flips
                a = abs(scr.znull[p])
                a > out.maxabs[slot, p, j] && (out.maxabs[slot, p, j] = a)
            end
        end
    end
    return
end

# ============================================================================
# 6. NULL DISTRIBUTIONS & CORRECTED P-VALUES
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

# Max cluster mass / max TFCE per sign pattern from the stored null z-maps
function _map_null!(null_max, zperm, tested, ctx)
    n_channels, n_timepoints, n_flips, _ = size(zperm)
    chunks = Iterators.partition(1:n_flips, cld(n_flips, Threads.nthreads()))
    tasks = map(chunks) do flip_chunk
        Threads.@spawn begin
            # One buffer set per task (never index buffers by threadid(): tasks may migrate)
            labels_buf = zeros(Int, n_channels, n_timepoints)
            queue_buf = Vector{Tuple{Int,Int}}(undef, n_channels * n_timepoints)
            pos_mask = BitArray{2}(undef, n_channels, n_timepoints)
            neg_mask = BitArray{2}(undef, n_channels, n_timepoints)
            t_map = Array{Float64,2}(undef, n_channels, n_timepoints)
            for p in flip_chunk, (j, c) in enumerate(tested)
                copyto!(t_map, view(zperm, :, :, p, j))
                if ctx.correction === :tfce
                    tfce_map = EegFun._compute_tfce(t_map, ctx.channel_names, Float64.(ctx.time_points), ctx.spatial_conn,
                        :spatiotemporal; E = ctx.tfce_E, H = ctx.tfce_H, dh = ctx.tfce_dh)
                    null_max[p, c] = maximum(abs, tfce_map; init = 0.0)
                else
                    @inbounds for i in eachindex(t_map)
                        pos_mask[i] = t_map[i] > ctx.cluster_threshold
                        neg_mask[i] = t_map[i] < -ctx.cluster_threshold
                    end
                    max_pos = EegFun._max_cluster_stat(pos_mask, t_map, ctx.spatial_conn, :spatiotemporal, labels_buf, queue_buf; is_negative = false)
                    min_neg = EegFun._max_cluster_stat(neg_mask, t_map, ctx.spatial_conn, :spatiotemporal, labels_buf, queue_buf; is_negative = true)
                    null_max[p, c] = max(max_pos, abs(min_neg))
                end
            end
        end
    end
    foreach(wait, tasks)
    return null_max
end

function _compute_p_corrected(z, null_max, tested, ctx)
    n_channels, n_timepoints, n_coefs = size(z)
    # Untested coefficients have no permutation p-value: NaN (not 1.0, which reads as "tested, not significant")
    p_corr = fill(NaN, n_channels, n_timepoints, n_coefs)
    electrode_to_idx = Dict(e => i for (i, e) in enumerate(ctx.channel_names))
    for c in tested
        null_sorted = sort(null_max[:, c])
        zc = z[:, :, c]
        stat_map = if ctx.correction === :maxt
            abs.(zc)
        elseif ctx.correction === :tfce
            abs.(EegFun._compute_tfce(zc, ctx.channel_names, Float64.(ctx.time_points), ctx.spatial_conn, :spatiotemporal;
                E = ctx.tfce_E, H = ctx.tfce_H, dh = ctx.tfce_dh))
        else
            _build_cluster_stat_map(zc, ctx.channel_names, ctx.time_points, ctx.spatial_conn, electrode_to_idx, ctx.cluster_threshold)
        end
        for i in eachindex(stat_map)
            p_corr[i + (c - 1) * n_channels * n_timepoints] = EegFun._lmm_perm_pvalue(null_sorted, stat_map[i], ctx.exhaustive)
        end
    end
    return p_corr
end

# ============================================================================
# 7. PRIMARY MASS-UNIVARIATE LMM FITTER
# ============================================================================

function fit_mass_lmm(eeg_data::AbstractArray{<:Real,3}, meta_df::DataFrame, f::FormulaTerm; dim_order = (:channels, :timepoints, :epochs), kwargs...)
    if dim_order != (:channels, :timepoints, :epochs)
        target = (:channels, :timepoints, :epochs)
        perm = [findfirst(==(t), dim_order) for t in target]
        any(isnothing, perm) && @minimal_error "dim_order must contain :channels, :timepoints, and :epochs"
        eeg_data = permutedims(eeg_data, tuple(perm...))
    end
    return _fit_mass_lmm(eeg_data, meta_df, f, nothing; kwargs...)
end

function _fit_mass_lmm(eeg_data, meta_df, f, epochs;
    n_perms::Integer = 1000,
    tested_coefs = nothing,
    test_intercept::Bool = false,
    flip_by = nothing,
    correction::Symbol = :maxt,
    cluster_threshold::Real = 2.0,
    tfce_E::Real = 0.5,
    tfce_H::Real = 2.0,
    tfce_dh::Real = 0.1,
    spatial_connectivity = nothing,
    channel_names = nothing,
    time_points = nothing,
    REML::Bool = false,
    rng::AbstractRNG = Random.default_rng()
)
    correction in (:maxt, :cluster, :tfce) || @minimal_error "correction must be :maxt, :cluster or :tfce (got :$correction)"
    n_perms >= 0 || @minimal_error "n_perms must be ≥ 0"
    n_channels, n_timepoints, n_epochs = size(eeg_data)
    nrow(meta_df) == n_epochs || @minimal_error "Number of rows in meta_df ($(nrow(meta_df))) must match number of epochs in eeg_data ($n_epochs)."
    channel_names = isnothing(channel_names) ? [Symbol("Ch$i") for i in 1:n_channels] : Symbol.(channel_names)
    time_points = Float64.(isnothing(time_points) ? (1:n_timepoints) : time_points)

    # Many tiny BLAS calls per model: single-threaded BLAS is faster. Restore the user's setting afterwards.
    blas_threads = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    try
        # 1. Design, tested coefficients, exchangeability blocks
        m0, coef_names, dropped, tested, meta_df_init = _setup_lmm_design(f, meta_df, n_epochs, tested_coefs, test_intercept)
        n_coefs = length(coef_names)
        est_cols = setdiff(1:n_coefs, dropped)
        do_perm = n_perms > 0 && !isempty(tested)

        flip_col, block_ids, H, exhaustive = nothing, nothing, zeros(0, 0), false
        if do_perm
            flip_col, block_ids = _resolve_flip_by(m0, meta_df, flip_by)
            n_blocks = maximum(block_ids)
            2.0^(n_blocks - 1) >= MIN_SIGN_FLIPS || @minimal_error "Only $n_blocks levels of :$flip_col: whole-block sign flipping " *
                "has $(2^(n_blocks - 1)) distinct patterns, too few to reach p < 0.05. At least 6 blocks are required."
            n_blocks < 8 && @warn "Only $n_blocks levels of :$flip_col: the permutation distribution is coarse (≤ $(2^(n_blocks - 1)) distinct values)."
            H, exhaustive = _sign_flip_matrix(n_blocks, n_perms, rng)
            @info "Permutation inference: whole-block sign flipping by :$flip_col ($n_blocks blocks, " *
                  (exhaustive ? "all $(size(H, 2)) distinct patterns" : "$(size(H, 2)) random patterns") *
                  "), correction = :$correction, testing $(join(coef_names[tested], ", "))"
        end
        n_flips = size(H, 2)
        store_maps = do_perm && correction !== :maxt
        tested_pos = [findfirst(==(c), est_cols) for c in tested]

        # 2. Outputs (SharedArrays when worker processes exist)
        use_workers = nprocs() > 1
        alloc(T, dims...) = use_workers ? SharedArray{T}(dims) : zeros(T, dims)
        tasks = [(c, t) for c in 1:n_channels for t in 1:n_timepoints]
        n_slots = use_workers ? maximum(workers()) : max(1, min(length(tasks), 4 * Threads.nthreads()))
        out = (
            beta = alloc(Float64, n_channels, n_timepoints, n_coefs),
            se = alloc(Float64, n_channels, n_timepoints, n_coefs),
            z = alloc(Float64, n_channels, n_timepoints, n_coefs),
            p = alloc(Float64, n_channels, n_timepoints, n_coefs),
            singular = alloc(Bool, n_channels, n_timepoints),
            zperm = store_maps ? alloc(Float32, n_channels, n_timepoints, n_flips, length(tested)) : zeros(Float32, 0, 0, 0, 0),
            maxabs = (do_perm && !store_maps) ? alloc(Float64, n_slots, n_flips, length(tested)) : zeros(0, 0, 0),
        )
        eeg = (use_workers && !(eeg_data isa SharedArray)) ? (SharedArray{Float64}(size(eeg_data)) .= eeg_data) : eeg_data
        base_ctx = (REML = REML, n_flips = do_perm ? n_flips : 0, tested_pos = tested_pos, H = H, store_maps = store_maps)

        # 3. Mass-univariate fits: one LMM per pixel; the sign-flip null is closed-form at the fitted θ
        par_desc = use_workers ? "$(nworkers()) worker processes" : "$(Threads.nthreads()) threads"
        @info "Fitting $(length(tasks)) mixed models across $n_epochs epochs on $par_desc..."
        prog = Progress(length(tasks); dt = 1, desc = "Fitting models...")
        Logging.with_logger(Logging.NullLogger()) do
            if use_workers
                prog_channel = RemoteChannel(() -> Channel{Bool}(length(tasks)), 1)
                @sync begin
                    @async while take!(prog_channel)
                        next!(prog)
                    end
                    @sync @distributed for (c_idx, t_idx) in tasks
                        scr = get!(task_local_storage(), :lmm_scratch) do
                            _create_scratch(f, meta_df_init, est_cols, block_ids, n_flips)
                        end
                        _process_pixel!(scr, eeg, c_idx, t_idx, out, merge(base_ctx, (slot = myid(),)))
                        put!(prog_channel, true)
                    end
                    put!(prog_channel, false)
                end
            else
                chunks = collect(Iterators.partition(tasks, cld(length(tasks), n_slots)))
                prog_lock = ReentrantLock()
                @sync for (slot, chunk) in enumerate(chunks)
                    Threads.@spawn begin
                        scr = _create_scratch(f, meta_df_init, est_cols, block_ids, n_flips)
                        ctx = merge(base_ctx, (slot = slot,))
                        for (c_idx, t_idx) in chunk
                            _process_pixel!(scr, eeg, c_idx, t_idx, out, ctx)
                            lock(() -> next!(prog), prog_lock)
                        end
                    end
                end
            end
        end
        finish!(prog)

        beta, se, z, p_param = Array(out.beta), Array(out.se), Array(out.z), Array(out.p)
        for arr in (beta, se, z, p_param)
            arr[:, :, dropped] .= NaN
        end

        # 4. Null distributions & corrected p-values
        null_max = fill(NaN, n_flips, n_coefs)
        p_corrected = nothing
        spatial_conn = nothing
        if do_perm
            spatial_conn = isnothing(spatial_connectivity) ? sparse(Int[], Int[], Bool[], n_channels, n_channels) : SparseMatrixCSC{Bool}(spatial_connectivity)
            spatial_conn = issymmetric(spatial_conn) ? spatial_conn : (spatial_conn .| copy(transpose(spatial_conn)))
            ctx = (correction = correction, cluster_threshold = Float64(cluster_threshold), tfce_E = tfce_E, tfce_H = tfce_H,
                tfce_dh = tfce_dh, channel_names = channel_names, time_points = time_points, spatial_conn = spatial_conn,
                exhaustive = exhaustive)
            if store_maps
                @info correction === :tfce ? "Computing TFCE null distribution..." : "Computing cluster-mass null distribution (threshold = $cluster_threshold)..."
                _map_null!(null_max, out.zperm, tested, ctx)
            else
                maxabs = Array(out.maxabs)
                for (j, c) in enumerate(tested)
                    null_max[:, c] .= vec(maximum(view(maxabs, :, :, j); dims = 1))
                end
            end
            p_corrected = _compute_p_corrected(z, null_max, tested, ctx)
        end

        inference = (
            correction = do_perm ? correction : :none,
            n_perms = n_flips,
            exhaustive = exhaustive,
            flip_by = flip_col,
            n_blocks = isnothing(block_ids) ? 0 : maximum(block_ids),
            cluster_threshold = Float64(cluster_threshold),
            tfce = (E = Float64(tfce_E), H = Float64(tfce_H), dh = Float64(tfce_dh)),
            spatial_connectivity = spatial_conn,
            REML = REML,
            n_obs = n_epochs,
            formula = f,
        )
        return EegFun.LmmStatsResult(coef_names, channel_names, time_points, beta, se, z, p_param, p_corrected,
            null_max, tested, Array(out.singular), inference, epochs)
    finally
        LinearAlgebra.BLAS.set_num_threads(blas_threads)
    end
end

end # module
