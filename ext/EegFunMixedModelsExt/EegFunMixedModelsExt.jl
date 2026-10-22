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

# Maximum optimizer evaluations per fit (a fit hitting it is reported as :max_evals)
const MAX_EVALS = 1000

# Fit outcome codes: zero-based indices into `EegFun.LMM_FIT_STATUS`
const STATUS_CONVERGED = 0x00
const STATUS_ROUNDOFF = 0x01
const STATUS_MAXEVAL = 0x02
const STATUS_OPTFAIL = 0x03
const STATUS_FLAT = 0x04    # constant / non-finite data: not fitted
const STATUS_ERROR = 0x05   # fit threw, or produced non-finite estimates

# Named predicate for statistical inference eligibility.
# Only successfully converged fits and roundoff-limited fits (which in MixedModels/NLopt
# indicate convergence at machine precision or on a boundary θ_k -> 0) are eligible to
# enter the observed statistic map, the permutation null distribution, and cluster/TFCE formation.
# Fits that reached evaluation limits (STATUS_MAXEVAL), optimizer failures (STATUS_OPTFAIL),
# flat/constant data (STATUS_FLAT), or threw errors (STATUS_ERROR) are ineligible.
is_inference_eligible(status::UInt8) = status <= STATUS_ROUNDOFF

# Stored null maps (:cluster / :tfce) are warned about / refused above these fractions of free memory
const ZPERM_WARN_FRACTION = 0.5
const ZPERM_ERROR_FRACTION = 0.9

# Below this many random sign flips the permutation p-values are very coarse
const SMALL_N_PERMS = 100

# Map the optimizer's return code to a fit outcome. The names are NLopt's (MixedModels' default backend);
# anything unrecognised is conservatively reported as an optimizer failure.
function _classify_returnvalue(rv::Symbol)
    rv in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :SMALL_TR_RADIUS, :FTARGET_ACHIEVED) && return STATUS_CONVERGED
    rv === :ROUNDOFF_LIMITED && return STATUS_ROUNDOFF
    rv in (:MAXEVAL_REACHED, :MAXTIME_REACHED, :MAXFUN_REACHED, :MAXTR_REACHED) && return STATUS_MAXEVAL
    return STATUS_OPTFAIL
end

# Refuse / warn when the stored null z-maps (channels × time × flips × tested, Float32) would not fit in memory
function _check_zperm_memory(n_channels, n_timepoints, n_flips, n_tested; free::Real = Sys.free_memory())
    bytes = 4.0 * n_channels * n_timepoints * n_flips * n_tested
    gb(x) = round(x / 2^30; digits = 1)
    detail = "The cluster/TFCE null needs the z-map of every permutation ($n_channels channels × $n_timepoints time points × " *
             "$n_flips permutations × $n_tested tested coefficients = $(gb(bytes)) GB, $(gb(free)) GB currently free). "
    advice = "Reduce `n_perms`, test fewer coefficients (`tested_coefs`), use a shorter time window or fewer channels, " *
             "or use `correction = :maxt`, which does not store the null maps."
    if bytes > ZPERM_ERROR_FRACTION * free
        @minimal_error detail * "This exceeds the available memory. " * advice
    elseif bytes > ZPERM_WARN_FRACTION * free
        @warn detail * "This is a large share of the available memory. " * advice
    end
    return bytes
end

# Random sign flips give p-values with resolution 1/(n+1); say so when that is too coarse
function _warn_small_n_perms(n_flips::Int, exhaustive::Bool)
    exhaustive && return
    pmin = 1 / (n_flips + 1)
    if pmin >= 0.05
        @warn "Only $n_flips random sign flips: the smallest attainable corrected p-value is $(round(pmin; digits = 3)), " *
              "so no point can reach p < 0.05. Use at least 1000 permutations."
    elseif n_flips < SMALL_N_PERMS
        @warn "Only $n_flips random sign flips: corrected p-values have a resolution of $(round(pmin; digits = 3)) and " *
              "are noisy. Use at least 1000 permutations for reportable results."
    end
    return
end

# One consolidated warning for points that could not be fitted or did not converge
function _report_fit_problems(status::AbstractMatrix{UInt8}, channel_names, first_error)
    n_total = length(status)
    n_flat, n_err = count(==(STATUS_FLAT), status), count(==(STATUS_ERROR), status)
    n_max, n_opt = count(==(STATUS_MAXEVAL), status), count(==(STATUS_OPTFAIL), status)
    n_unfit, n_unconv = n_flat + n_err, n_max + n_opt
    n_unfit + n_unconv == 0 && return
    msgs = String[]
    if n_unfit > 0
        bad = status .>= STATUS_FLAT
        by_ch = [(channel_names[c], count(@view bad[c, :])) for c in axes(bad, 1) if any(@view bad[c, :])]
        shown = join(["$ch ($k/$(size(bad, 2)) time points)" for (ch, k) in first(by_ch, 5)], ", ")
        more = length(by_ch) > 5 ? ", ... ($(length(by_ch)) channels in total)" : ""
        push!(msgs, "$n_unfit of $n_total channel × time points could not be fitted (flat or non-finite data: $n_flat; " *
                    "fit error: $n_err). Their beta, se, z and p are NaN, and they are excluded from the permutation null and " *
                    "from cluster/TFCE formation. Affected: $shown$more.")
        isnothing(first_error) || push!(msgs, "First fit error: $first_error")
    end
    if n_unconv > 0
        push!(msgs, "$n_unconv of $n_total fits did not converge (evaluation limit reached: $n_max; other optimizer " *
                    "failure: $n_opt). Their estimates are kept but may be unreliable; see `fit_status` in the result.")
    end
    @warn join(msgs, "\n")
    return
end

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
        return (isempty(effects) && !isnothing(intercept_idx)) ? [intercept_idx] : effects
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
    all(nrow(df) == n_timepoints for df in dfs) || @minimal_error "All epochs must have exactly $n_timepoints time points."
    :channel_names in keys(kwargs) && @minimal_error "Cannot pass `channel_names` when calling fit_mass_lmm with EpochData."
    :time_points in keys(kwargs) && @minimal_error "Cannot pass `time_points` when calling fit_mass_lmm with EpochData."
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
    # Monte Carlo sampling: draw independent sign vectors from {-1, 1}^S.
    # To strictly avoid double-counting the observed identity pattern under
    # the Phipson & Smyth (2010) formulation p = (n_ge + 1) / (B + 1),
    # ensure that none of the sampled columns is all +1.0 or all -1.0
    # (which under two-sided |z| invariance yields the observed statistic).
    H = rand(rng, [-1.0, 1.0], n_blocks, n_perms)
    for col in 1:n_perms
        while all(==(1.0), @view H[:, col]) || all(==(-1.0), @view H[:, col])
            H[:, col] .= rand(rng, [-1.0, 1.0], n_blocks)
        end
    end
    return H, false
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
# Under Freedman–Lane GLS on pseudo-data y*(η) = X_(-c) γ̂ + r*(η), the nuisance projection satisfies
# (I - P_X) X_(-c) γ̂ = 0 because X_(-c) ⊆ col(X). Hence the full GLS residual is e*(η) = (I - P_X) r*(η),
# and because η_s² = 1, r*(η)' V⁻¹ r*(η) = Σ r_s' V_s⁻¹ r_s = d identically for all η.
# The target coefficient is β̂*_c = (M W η)_c and the penalized residual quadratic form is
# pwrss*(η) = d − (Wη)' M (Wη) identically for ALL sign patterns η ∈ {-1, +1}^S (machine precision match to dense GLS).
# At η = 1 this reproduces the LMM's own Wald z exactly.
function _null_z!(out::AbstractVector, XtVX_s, XtVy_s, ytVy_s, pc::Int, H::Matrix{Float64}, ν::Real)
    p, n_blocks = size(first(XtVX_s), 1), length(XtVy_s)
    XtVX, XtVy = sum(XtVX_s), sum(XtVy_s)
    M = inv(cholesky(Symmetric(XtVX)))
    (ν > 0 && M[pc, pc] > 0) || return fill!(out, NaN)

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
    tol = 1e-11 * max(d, 1.0)
    @inbounds for j in axes(H, 2)
        quad = 0.0
        for i in 1:p
            quad += U[i, j] * MU[i, j]
        end
        diff = d - quad
        if diff < -tol
            out[j] = NaN  # Materially negative residual variance indicates numerical breakdown
        else
            diff_clamped = max(diff, 0.0)
            denom = sqrt(diff_clamped * scale)
            out[j] = (denom > 0 && isfinite(denom)) ? (MU[pc, j] / denom) : NaN
        end
    end
    return out
end

# ============================================================================
# 4. MODEL DESIGN SETUP
# ============================================================================

function _setup_lmm_design(f, meta_df, n_epochs, tested_coefs, test_intercept, contrasts)
    meta_df_init = copy(meta_df)
    meta_df_init[!, Symbol(f.lhs)] = randn(Random.MersenneTwister(0), n_epochs)

    @info "Compiling MixedModel design matrix..."
    # Coefficient names, rank and pivot are available from the unfitted model
    m_initial = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df_init; contrasts = contrasts)
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
    design::Union{BlockDesign,Nothing}
    y::Vector{Float64}
    beta::Vector{Float64}       # results of the current pixel, committed to the shared outputs only on success
    se::Vector{Float64}
    z::Vector{Float64}
    p::Vector{Float64}
    znull::Matrix{Float64}      # n_flips × n_tested
    singular::Base.RefValue{Bool}
end

function _create_scratch(f, meta_df_init, est_cols, block_ids, n_flips, n_tested, contrasts)
    LinearAlgebra.BLAS.set_num_threads(1)
    m = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df_init; contrasts = contrasts)
    end
    m.optsum.ftol_rel = OBSERVED_FTOL_REL
    m.optsum.maxfeval = MAX_EVALS
    X = Matrix{Float64}(modelmatrix(m)[:, invperm(MixedModels.pivot(m))][:, est_cols])
    design = isnothing(block_ids) ? nothing : BlockDesign(m, X, block_ids)
    n_obs = length(m.y)
    n_coefs = length(coefnames(m))
    return LmmScratch(m, design, zeros(n_obs), zeros(n_coefs), zeros(n_coefs), zeros(n_coefs), zeros(n_coefs),
        zeros(n_flips, n_tested), Ref(false))
end

# Fit one pixel and compute its sign-flip null into the scratch buffers (nothing is written to shared outputs).
# Returns the fit outcome code; may throw. Flat or non-finite data are not fitted (MixedModels throws for constant
# data and silently returns NaN for non-finite data).
function _compute_pixel!(scr::LmmScratch, ctx)
    y = scr.y
    (all(isfinite, y) && any(!=(y[1]), y)) || return STATUS_FLAT
    m = scr.m
    refit!(m, y; REML = ctx.REML, progress = false)
    status = _classify_returnvalue(m.optsum.returnvalue)
    scr.singular[] = issingular(m)
    MixedModels.stderror!(scr.se, m)
    copyto!(scr.beta, coef(m))
    @inbounds for k in eachindex(scr.se)
        z = scr.beta[k] / scr.se[k]
        scr.z[k] = z
        scr.p[k] = 2.0 * ccdf(Normal(), abs(z))
    end
    # Aliased (dropped) coefficients are NaN by design; every estimable one must be finite
    for k in ctx.est_cols
        (isfinite(scr.beta[k]) && isfinite(scr.z[k])) || return STATUS_ERROR
    end

    ctx.n_flips == 0 && return status
    lams = Matrix{Float64}[Matrix(A.λ) for A in m.reterms]
    stats = _block_gls_stats(scr.design, lams, scr.y)
    ν = MixedModels.ssqdenom(m)
    for (j, pc) in enumerate(ctx.tested_pos)
        zn = view(scr.znull, :, j)
        _null_z!(zn, stats..., pc, ctx.H, ν)
        all(isfinite, zn) || return STATUS_ERROR
    end
    return status
end

# Copy a successfully computed pixel into the shared outputs
function _commit_pixel!(scr::LmmScratch, out, c_idx, t_idx, ctx, eligible::Bool)
    out.singular[c_idx, t_idx] = scr.singular[]
    @inbounds for k in eachindex(scr.beta)
        out.beta[c_idx, t_idx, k] = scr.beta[k]
        out.se[c_idx, t_idx, k] = scr.se[k]
        out.z[c_idx, t_idx, k] = scr.z[k]
        out.p[c_idx, t_idx, k] = scr.p[k]
    end
    (ctx.n_flips == 0 || !eligible) && return
    for j in eachindex(ctx.tested_pos)
        if ctx.store_maps
            @inbounds for p in 1:ctx.n_flips
                out.zperm[c_idx, t_idx, p, j] = scr.znull[p, j]
            end
        else
            slot = myid()
            @inbounds for p in 1:ctx.n_flips
                a = abs(scr.znull[p, j])
                a > out.maxabs[slot, p, j] && (out.maxabs[slot, p, j] = a)
            end
        end
    end
    return
end

function _process_pixel!(scr::LmmScratch, eeg, c_idx, t_idx, out, ctx)
    scr.y .= @view eeg[c_idx, t_idx, :]
    status = try
        _compute_pixel!(scr, ctx)
    catch e
        e isa Union{InterruptException,UndefVarError,BoundsError,TypeError} && rethrow()
        STATUS_ERROR
    end
    eligible = is_inference_eligible(status)
    status < STATUS_FLAT && _commit_pixel!(scr, out, c_idx, t_idx, ctx, eligible)
    out.status[c_idx, t_idx] = status
    return
end

# Re-run the first pixel that errored on this process to report the reason (workers cannot return messages)
function _first_fit_error(status, eeg, make_scratch, ctx)
    idx = findfirst(==(STATUS_ERROR), status)
    isnothing(idx) && return nothing
    c, t = Tuple(idx)
    scr = make_scratch()
    scr.y .= @view eeg[c, t, :]
    reason = try
        code = Logging.with_logger(Logging.NullLogger()) do
            _compute_pixel!(scr, ctx)
        end
        code == STATUS_ERROR ? "non-finite estimates or null statistics (e.g. zero standard error)" : "could not be reproduced"
    catch e
        sprint(showerror, e)
    end
    return "channel index $c, time index $t: $reason"
end

# ============================================================================
# 6. NULL DISTRIBUTIONS & CORRECTED P-VALUES
# ============================================================================

function _build_cluster_stat_map(t_map, channel_names, time_points, spatial_conn, cluster_threshold)
    n_channels, n_timepoints = size(t_map)
    pos_clusters, neg_clusters = EegFun._find_clusters(
        t_map .> cluster_threshold, t_map .< -cluster_threshold,
        channel_names, time_points, spatial_conn, :spatiotemporal
    )
    pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map; return_clusters=false)
    neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map; return_clusters=false)

    smap = zeros(Float64, n_channels, n_timepoints)
    for (clusters, stats) in ((pos_clusters, pos_stats), (neg_clusters, neg_stats))
        for (i, cluster) in enumerate(clusters)
            stat = abs(stats[i])
            for (e_idx, t_idx) in cluster.members
                smap[e_idx, t_idx] = max(smap[e_idx, t_idx], stat)
            end
        end
    end
    return smap
end

# Max cluster mass / max TFCE per sign pattern from the stored null z-maps
function _map_null!(null_max, zperm, tested, ctx)
    n_channels, n_timepoints, n_flips, _ = size(zperm)
    chunks = Iterators.partition(1:n_flips, cld(n_flips, Threads.nthreads()))
    @sync for flip_chunk in chunks
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
                    tfce_map = EegFun._compute_tfce(t_map, ctx.channel_names, ctx.time_points, ctx.spatial_conn,
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
    return null_max
end

function _compute_p_corrected(z, null_max, tested, ctx, bad::AbstractMatrix{Bool})
    n_channels, n_timepoints, n_coefs = size(z)
    # Untested coefficients have no permutation p-value: NaN (not 1.0, which reads as "tested, not significant")
    p_corr = fill(NaN, n_channels, n_timepoints, n_coefs)
    for c in tested
        null_sorted = sort(null_max[:, c])
        zc = copy(z[:, :, c])
        zc[bad] .= 0.0   # unfitted points (NaN) carry no evidence and never form clusters
        stat_map = if ctx.correction === :maxt
            abs.(zc)
        elseif ctx.correction === :tfce
            abs.(EegFun._compute_tfce(zc, ctx.channel_names, ctx.time_points, ctx.spatial_conn, :spatiotemporal;
                E = ctx.tfce_E, H = ctx.tfce_H, dh = ctx.tfce_dh))
        else
            _build_cluster_stat_map(zc, ctx.channel_names, ctx.time_points, ctx.spatial_conn, ctx.cluster_threshold)
        end
        p_c = view(p_corr, :, :, c)
        for i in eachindex(stat_map)
            p_c[i] = EegFun._lmm_perm_pvalue(null_sorted, stat_map[i], ctx.exhaustive)
        end
        p_c[bad] .= NaN
    end
    return p_corr
end

# ============================================================================
# 7. PRIMARY MASS-UNIVARIATE LMM FITTER
# ============================================================================

function fit_mass_lmm(eeg_data::AbstractArray{<:Real,3}, meta_df::DataFrame, f::FormulaTerm; dim_order = (:channels, :timepoints, :epochs), kwargs...)
    target = (:channels, :timepoints, :epochs)
    (length(dim_order) == 3 && Set(dim_order) == Set(target)) ||
        @minimal_error "dim_order must be a 3-element permutation of (:channels, :timepoints, :epochs), got: $dim_order"
    if dim_order != target
        perm = [findfirst(==(t), dim_order) for t in target]
        eeg_data = permutedims(eeg_data, tuple(perm...))
    end
    return _fit_mass_lmm(eeg_data, meta_df, f, nothing; kwargs...)
end

# Only what the workers need in Main: the extension's trigger packages (see [extensions] in EegFun's
# Project.toml: MixedModels and StatsModels) plus SharedArrays. Packages that are merely dependencies
# of EegFun are resolved by the extension itself.
const _WORKER_USING = :(using EegFun, MixedModels, StatsModels, SharedArrays)

function _ensure_workers_loaded()
    ws = filter(!=(1), workers())
    isempty(ws) && return
    proj = Base.active_project()
    try
        @sync for p in ws
            @async remotecall_eval(Main, p, quote
                using Pkg
                Pkg.activate($proj; io = devnull)
                using Logging
                Logging.global_logger(Logging.NullLogger())
                using EegFun, MixedModels, StatsModels, SharedArrays
            end)
        end
    catch e
        @minimal_error "Could not load the required packages on the worker processes: $(sprint(showerror, e)). " *
                       "Start workers with `addprocs(n; exeflags = \"--project=\$(Base.active_project())\")`."
    end
    for p in ws
        ok = remotecall_fetch(Core.eval, p, Main, :(!isnothing(Base.get_extension(EegFun, :EegFunMixedModelsExt))))
        ok || @minimal_error "EegFunMixedModelsExt did not load on worker $p. Start workers with " *
                             "`addprocs(n; exeflags = \"--project=\$(Base.active_project())\")`."
    end
    return
end

function _set_worker_blas_threads!(n::Int)
    old = Dict{Int,Int}()
    for p in workers()
        p == 1 && continue
        old[p] = remotecall_fetch(Core.eval, p, Main,
            :(using LinearAlgebra; let o = BLAS.get_num_threads(); BLAS.set_num_threads($n); o end))
    end
    return old
end

function _restore_worker_blas!(old::Dict{Int,Int})
    for (p, n) in old
        try
            remotecall_fetch(Core.eval, p, Main, :(using LinearAlgebra; BLAS.set_num_threads($n)))
        catch
        end
    end
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
    contrasts::AbstractDict = Dict{Symbol,Any}(),
    progress::Bool = true,
    rng::AbstractRNG = Random.default_rng()
)
    correction in (:maxt, :cluster, :tfce) || @minimal_error "correction must be :maxt, :cluster or :tfce (got :$correction)"
    n_perms >= 0 || @minimal_error "n_perms must be ≥ 0"
    n_channels, n_timepoints, n_epochs = size(eeg_data)
    nrow(meta_df) == n_epochs || @minimal_error "Number of rows in meta_df ($(nrow(meta_df))) must match number of epochs in eeg_data ($n_epochs)."
    if !isnothing(channel_names)
        length(channel_names) == n_channels || @minimal_error "length(channel_names) ($(length(channel_names))) must match n_channels ($n_channels)."
    end
    if !isnothing(time_points)
        length(time_points) == n_timepoints || @minimal_error "length(time_points) ($(length(time_points))) must match n_timepoints ($n_timepoints)."
        all(isfinite, time_points) || @minimal_error "time_points must contain finite values."
    end
    channel_names = isnothing(channel_names) ? [Symbol("Ch$i") for i in 1:n_channels] : Symbol.(channel_names)
    time_points = Float64.(isnothing(time_points) ? (1:n_timepoints) : time_points)

    # Many tiny BLAS calls per model: single-threaded BLAS is faster. Restore the user's setting afterwards.
    blas_threads = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    worker_blas = Dict{Int,Int}()
    try
        # 1. Design, tested coefficients, exchangeability blocks
        m0, coef_names, dropped, tested, meta_df_init = _setup_lmm_design(f, meta_df, n_epochs, tested_coefs, test_intercept, contrasts)
        n_coefs = length(coef_names)
        est_cols = setdiff(1:n_coefs, dropped)
        do_perm = n_perms > 0 && !isempty(tested)
        do_perm && correction !== :maxt && isnothing(spatial_connectivity) && n_channels > 1 &&
            @warn "No spatial connectivity or layout available: clusters/TFCE will form across time only."

        flip_col, block_ids, H, exhaustive = nothing, nothing, zeros(0, 0), false
        if do_perm
            flip_col, block_ids = _resolve_flip_by(m0, meta_df, flip_by)
            n_blocks = maximum(block_ids)
            2.0^(n_blocks - 1) >= MIN_SIGN_FLIPS || @minimal_error "Only $n_blocks levels of :$flip_col: whole-block sign flipping " *
                "has $(2^(n_blocks - 1)) distinct patterns, too few to reach p < 0.05. At least 6 blocks are required."
            n_blocks < 8 && @warn "Only $n_blocks levels of :$flip_col: the permutation distribution is coarse (≤ $(2^(n_blocks - 1)) distinct values)."
            H, exhaustive = _sign_flip_matrix(n_blocks, n_perms, rng)
            _warn_small_n_perms(size(H, 2), exhaustive)
            @info "Permutation inference: whole-block sign flipping by :$flip_col ($n_blocks blocks, " *
                  (exhaustive ? "all $(size(H, 2)) distinct patterns" : "$(size(H, 2)) random patterns") *
                  "), correction = :$correction, testing $(join(coef_names[tested], ", "))"
        end
        n_flips = size(H, 2)
        store_maps = do_perm && correction !== :maxt
        store_maps && _check_zperm_memory(n_channels, n_timepoints, n_flips, length(tested))
        tested_pos = [findfirst(==(c), est_cols) for c in tested]

        # 2. Worker initialization & Outputs (SharedArrays across workers, or regular Arrays if single process)
        is_distributed = nprocs() > 1
        if is_distributed
            _ensure_workers_loaded()
            worker_blas = _set_worker_blas_threads!(1)
        end

        tasks = [(c, t) for c in 1:n_channels for t in 1:n_timepoints]
        n_slots = is_distributed ? maximum(workers()) : 1
        out = if is_distributed
            (
                beta = SharedArray{Float64}(n_channels, n_timepoints, n_coefs),
                se = SharedArray{Float64}(n_channels, n_timepoints, n_coefs),
                z = SharedArray{Float64}(n_channels, n_timepoints, n_coefs),
                p = SharedArray{Float64}(n_channels, n_timepoints, n_coefs),
                singular = SharedArray{Bool}(n_channels, n_timepoints),
                status = SharedArray{UInt8}(n_channels, n_timepoints),
                zperm = store_maps ? SharedArray{Float32}(n_channels, n_timepoints, n_flips, length(tested)) : zeros(Float32, 0, 0, 0, 0),
                maxabs = (do_perm && !store_maps) ? SharedArray{Float64}(n_slots, n_flips, length(tested)) : zeros(0, 0, 0),
            )
        else
            (
                beta = fill(NaN, n_channels, n_timepoints, n_coefs),
                se = fill(NaN, n_channels, n_timepoints, n_coefs),
                z = fill(NaN, n_channels, n_timepoints, n_coefs),
                p = fill(NaN, n_channels, n_timepoints, n_coefs),
                singular = fill(false, n_channels, n_timepoints),
                status = fill(STATUS_FLAT, n_channels, n_timepoints),
                zperm = store_maps ? zeros(Float32, n_channels, n_timepoints, n_flips, length(tested)) : zeros(Float32, 0, 0, 0, 0),
                maxabs = (do_perm && !store_maps) ? zeros(Float64, n_slots, n_flips, length(tested)) : zeros(0, 0, 0),
            )
        end
        if is_distributed
            fill!(out.beta, NaN)
            fill!(out.se, NaN)
            fill!(out.z, NaN)
            fill!(out.p, NaN)
            fill!(out.singular, false)
            do_perm && !store_maps && fill!(out.maxabs, 0.0)
            store_maps && fill!(out.zperm, 0.0f0)
            fill!(out.status, STATUS_FLAT)
        end
        eeg = is_distributed ? ((eeg_data isa SharedArray) ? eeg_data : (SharedArray{Float64}(size(eeg_data)) .= eeg_data)) : eeg_data
        base_ctx = (REML = REML, n_flips = do_perm ? n_flips : 0, tested_pos = tested_pos, H = H, store_maps = store_maps,
            est_cols = est_cols)

        # 3. Mass-univariate fits: one LMM per pixel; the sign-flip null is closed-form at the fitted θ
        if is_distributed
            progress && @info "Fitting $(length(tasks)) mixed models across $n_epochs epochs on $(nworkers()) worker processes..."
        elseif progress
            @info "Fitting $(length(tasks)) mixed models across $n_epochs epochs on 1 worker. " *
                  "To parallelize across multiple CPU cores with up to 6x faster throughput, " *
                  "launch Julia with `julia -p auto` or call `addprocs()` before `fit_mass_lmm`."
        end

        make_scratch() = _create_scratch(f, meta_df_init, est_cols, block_ids, n_flips, length(tested), contrasts)
        prog = Progress(length(tasks); dt = 1, desc = "Fitting models...", enabled = progress)
        if is_distributed
            prog_channel = RemoteChannel(() -> Channel{Bool}(length(tasks)), 1)
            Logging.with_logger(Logging.NullLogger()) do
                @sync begin
                    @async while take!(prog_channel)
                        next!(prog)
                    end
                    @sync @distributed for (c_idx, t_idx) in tasks
                        scr = get!(make_scratch, task_local_storage(), :lmm_scratch)
                        _process_pixel!(scr, eeg, c_idx, t_idx, out, base_ctx)
                        put!(prog_channel, true)
                    end
                    put!(prog_channel, false)
                end
            end
        else
            scr = make_scratch()
            Logging.with_logger(Logging.NullLogger()) do
                for (c_idx, t_idx) in tasks
                    _process_pixel!(scr, eeg, c_idx, t_idx, out, base_ctx)
                    progress && next!(prog)
                end
            end
        end
        finish!(prog)

        # Points that could not be fitted / did not converge: one consolidated warning, NaN statistics already set
        fit_status = Array(out.status)
        first_error = _first_fit_error(fit_status, eeg, make_scratch, base_ctx)
        _report_fit_problems(fit_status, channel_names, first_error)
        unfitted = fit_status .>= STATUS_FLAT
        ineligible = .!is_inference_eligible.(fit_status)

        beta, se, z, p_param = Array(out.beta), Array(out.se), Array(out.z), Array(out.p)
        for arr in (beta, se, z, p_param)
            for c in 1:size(arr, 3)
                view(arr, :, :, c)[unfitted] .= NaN
            end
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
            p_corrected = _compute_p_corrected(z, null_max, tested, ctx, ineligible)
        end

        inference = (
            correction = do_perm ? correction : :none,
            n_perms = n_flips,
            exhaustive = exhaustive,
            p_value_convention = exhaustive ? :exhaustive_hypercube : :monte_carlo_phipson_smyth,
            n_unique_classes = isnothing(block_ids) ? 0 : (maximum(block_ids) - 1 < 62 ? 2^(maximum(block_ids) - 1) : -1),
            flip_by = flip_col,
            n_blocks = isnothing(block_ids) ? 0 : maximum(block_ids),
            cluster_threshold = Float64(cluster_threshold),
            tfce = (E = Float64(tfce_E), H = Float64(tfce_H), dh = Float64(tfce_dh)),
            spatial_connectivity = spatial_conn,
            REML = REML,
            contrasts = contrasts,
            n_obs = n_epochs,
            formula = f,
            prespecified_field_size = (n_channels, n_timepoints),
            field_mask_policy = :intersection_of_valid_pixels,
            n_valid_pixels = count(is_inference_eligible, fit_status),
            n_ineligible_pixels = count(!is_inference_eligible, fit_status),
            n_total_pixels = length(fit_status),
        )
        return EegFun.LmmStatsResult(coef_names, channel_names, time_points, beta, se, z, p_param, p_corrected,
            null_max, tested, Array(out.singular), fit_status, inference, epochs)
    finally
        _restore_worker_blas!(worker_blas)
        LinearAlgebra.BLAS.set_num_threads(blas_threads)
    end
end

end # module
