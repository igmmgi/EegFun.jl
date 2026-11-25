# Mass-univariate linear mixed-effects modeling (MixedModels.jl extension).

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

using SparseArrays

import EegFun: fit_mass_lmm, @minimal_error, @minimal_warning

######################################
######### Types & Constants ##########
######################################

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

# Stored null maps (:cluster / :tfce) are warned about / refused above these fractions of free memory
const ZPERM_WARN_FRACTION = 0.5

const ZPERM_ERROR_FRACTION = 0.9

# Below this many random sign flips the permutation p-values are very coarse
const SMALL_N_PERMS = 100

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

# Preallocated per-pixel CR2 & Woodbury workspaces
struct Cr2Workspace
    lams::Vector{Matrix{Float64}}
    XtVX_s::Vector{Matrix{Float64}}
    XtVy_s::Vector{Vector{Float64}}
    Λ_s::Vector{Matrix{Float64}}
    Zty_s::Vector{Vector{Float64}}
    Gty_s::Vector{Vector{Float64}}
    GtX_s::Vector{Matrix{Float64}}
    FGX_s::Vector{Matrix{Float64}}
    ZΛ_s::Vector{Matrix{Float64}}
    C_s::Vector{Matrix{Float64}}
    Xty::Vector{Float64}
    XtVX::Matrix{Float64}
    XtVy::Vector{Float64}
    chol_buf::Matrix{Float64}
    M::Matrix{Float64}
    R::Matrix{Float64}
    Hs_pre::Vector{Matrix{Float64}}
    As_pre::Vector{Matrix{Float64}}
    tmp_pp::Matrix{Float64}
    tmp_pp2::Matrix{Float64}
    W::Matrix{Float64}
    B::Matrix{Float64}
    F::Matrix{Float64}
    FtF::Matrix{Float64}
    Q::Matrix{Float64}
    Q_eta::Matrix{Float64}
    G::Matrix{Float64}
    psi::Vector{Float64}
    d::Vector{Float64}
    cs::Vector{Float64}
    fs::Vector{Float64}
    vs::Vector{Float64}
    q_obs::Vector{Float64}
end

# 5. PER-PIXEL WORK (one LMM fit + closed-form pivotal CR2 sign-flip null)
struct LmmScratch{M<:LinearMixedModel}
    m::M
    design::Union{BlockDesign,Nothing}
    cr2_ws::Union{Cr2Workspace,Nothing}
    y::Vector{Float64}
    beta::Vector{Float64}
    se::Vector{Float64}
    z::Vector{Float64}
    p::Vector{Float64}
    z_perm_obs::Vector{Float64}
    znull::Matrix{Float64}      # n_flips × n_tested
    singular::Base.RefValue{Bool}
end

######################################
########### User Functions ###########
######################################

function fit_mass_lmm(
    epochs::EegFun.EpochData,
    f::FormulaTerm;
    correction_method::Symbol = :maxt,
    correction::Union{Nothing,Symbol} = nothing,
    spatial_connectivity = nothing,
    kwargs...,
)
    eff_corr = isnothing(correction) ? correction_method : correction
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

    eeg_data = Array{Float64,3}(undef, n_channels, n_timepoints, n_epochs)
    for (i, df) in enumerate(dfs)
        for (c_idx, c_name) in enumerate(all_channels)
            _fill_lmm_epoch_slice!(eeg_data, c_idx, i, df[!, c_name], n_timepoints)
        end
    end

    if isnothing(spatial_connectivity) &&
       eff_corr in (:cluster, :tfce) &&
       !isnothing(epochs.layout) &&
       (!isnothing(epochs.layout.neighbors) || EegFun.has_valid_coordinates(epochs.layout))
        spatial_connectivity = EegFun._build_connectivity_matrix(all_channels, epochs.layout, :spatiotemporal)
    end

    tvec = EegFun.time_vector(epochs)
    time_pts = isempty(tvec) ? Vector(1:n_timepoints) : tvec

    return _fit_mass_lmm(
        eeg_data,
        meta_df,
        f,
        epochs;
        channel_names = all_channels,
        time_points = time_pts,
        correction_method = eff_corr,
        spatial_connectivity = spatial_connectivity,
        kwargs...,
    )
end

function fit_mass_lmm(
    epochs_vec::AbstractVector{<:EegFun.EpochData},
    f::FormulaTerm;
    condition_selection = EegFun.conditions(),
    correction_method::Symbol = :maxt,
    correction::Union{Nothing,Symbol} = nothing,
    spatial_connectivity = nothing,
    kwargs...,
)
    eff_corr = isnothing(correction) ? correction_method : correction
    isempty(epochs_vec) && @minimal_error "EpochData vector is empty."
    :channel_names in keys(kwargs) && @minimal_error "Cannot pass `channel_names` when calling fit_mass_lmm with EpochData."
    :time_points in keys(kwargs) && @minimal_error "Cannot pass `time_points` when calling fit_mass_lmm with EpochData."

    grouped = EegFun.group_by_condition(epochs_vec)
    reps = EegFun.EpochData[first(v) for v in values(grouped)]
    sel_idx = EegFun.get_selected_conditions(reps, condition_selection)
    isempty(sel_idx) && @minimal_error "No EpochData conditions matched `condition_selection`."
    keep_conds = Set{Int}(reps[i].condition for i in sel_idx)
    selected_sets = EegFun.EpochData[ep for ep in epochs_vec if ep.condition in keep_conds]

    nonempty_sets = EegFun.EpochData[ep for ep in selected_sets if !isempty(ep.data)]
    isempty(nonempty_sets) && @minimal_error "All selected EpochData objects are empty (0 trials)."

    ref_ep = first(nonempty_sets)
    all_channels = EegFun.channel_labels(ref_ep)
    n_channels = length(all_channels)
    n_timepoints = EegFun.n_samples(ref_ep)

    for (idx, ep) in enumerate(nonempty_sets)
        EegFun._have_same_structure(ref_ep, ep) ||
            @minimal_error "EpochData at index $idx (file='$(ep.file)', condition=$(ep.condition)) has inconsistent channel, sample-rate, or time structure."
        all(nrow(df) == n_timepoints for df in ep.data) ||
            @minimal_error "All epochs in '$(ep.file)' (condition $(ep.condition)) must have $n_timepoints time points."
    end

    all_dfs = DataFrame[df for ep in nonempty_sets for df in ep.data]
    total_epochs = length(all_dfs)

    common_meta = foldl(intersect, [setdiff(EegFun.all_labels(ep), all_channels) for ep in nonempty_sets])
    meta_df = DataFrame([col => [df[1, col] for df in all_dfs] for col in common_meta])

    distinct_conds = unique(ep.condition for ep in nonempty_sets)
    raw_names = [isempty(ep.condition_name) ? string(ep.condition) : ep.condition_name for ep in nonempty_sets]
    use_names = length(unique(raw_names)) >= length(distinct_conds)

    if !hasproperty(meta_df, :participant)
        meta_df[!, :participant] = String[_resolve_epoch_participant_label(ep) for ep in nonempty_sets for _ in ep.data]
    end
    if !hasproperty(meta_df, :condition)
        meta_df[!, :condition] = String[
            (use_names ? (isempty(ep.condition_name) ? string(ep.condition) : ep.condition_name) : string(ep.condition)) for
            ep in nonempty_sets for _ in ep.data
        ]
    end
    if !hasproperty(meta_df, :condition_num)
        meta_df[!, :condition_num] = Int[ep.condition for ep in nonempty_sets for _ in ep.data]
    end
    if !hasproperty(meta_df, :condition_name)
        meta_df[!, :condition_name] = String[ep.condition_name for ep in nonempty_sets for _ in ep.data]
    end
    if !hasproperty(meta_df, :file)
        meta_df[!, :file] = String[ep.file for ep in nonempty_sets for _ in ep.data]
    end

    eeg_data = Array{Float64,3}(undef, n_channels, n_timepoints, total_epochs)
    for (i, df) in enumerate(all_dfs)
        for (c_idx, c_name) in enumerate(all_channels)
            _fill_lmm_epoch_slice!(eeg_data, c_idx, i, df[!, c_name], n_timepoints)
        end
    end

    if isnothing(spatial_connectivity) &&
       eff_corr in (:cluster, :tfce) &&
       !isnothing(ref_ep.layout) &&
       (!isnothing(ref_ep.layout.neighbors) || EegFun.has_valid_coordinates(ref_ep.layout))
        spatial_connectivity = EegFun._build_connectivity_matrix(all_channels, ref_ep.layout, :spatiotemporal)
    end

    tvec = EegFun.time_vector(ref_ep)
    time_pts = isempty(tvec) ? Vector(1:n_timepoints) : tvec
    combined_ep = EegFun.EpochData(
        ref_ep.file,
        ref_ep.condition,
        ref_ep.condition_name,
        all_dfs,
        copy(ref_ep.layout),
        ref_ep.sample_rate,
        copy(ref_ep.analysis_info),
    )

    return _fit_mass_lmm(
        eeg_data,
        meta_df,
        f,
        combined_ep;
        channel_names = all_channels,
        time_points = time_pts,
        correction_method = eff_corr,
        spatial_connectivity = spatial_connectivity,
        kwargs...,
    )
end

# 7. PRIMARY MASS-UNIVARIATE LMM FITTER
function fit_mass_lmm(
    eeg_data::AbstractArray{<:Real,3},
    meta_df::DataFrame,
    f::FormulaTerm;
    dim_order = (:channels, :timepoints, :epochs),
    kwargs...,
)
    target = (:channels, :timepoints, :epochs)
    (length(dim_order) == 3 && Set(dim_order) == Set(target)) ||
        @minimal_error "dim_order must be a 3-element permutation of (:channels, :timepoints, :epochs), got: $dim_order"
    if dim_order != target
        perm = [findfirst(==(t), dim_order) for t in target]
        eeg_data = permutedims(eeg_data, tuple(perm...))
    end
    return _fit_mass_lmm(eeg_data, meta_df, f, nothing; kwargs...)
end

function BlockDesign(m::LinearMixedModel, X::Matrix{Float64}, block_ids::Vector{Int})
    n_blocks = maximum(block_ids)
    rows = [Int[] for _ = 1:n_blocks]
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
                Zk[i, ((j-1)*q+1):(j*q)] .= view(A.z, :, obs)
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

######################################
######### Internal Functions #########
######################################

# Named predicate for statistical inference eligibility.
# Only successfully converged fits and roundoff-limited fits (which in MixedModels/NLopt
# indicate convergence at machine precision or on a boundary θ_k -> 0) are eligible to
# enter the observed statistic map, the permutation null distribution, and cluster/TFCE formation.
# Fits that reached evaluation limits (STATUS_MAXEVAL), optimizer failures (STATUS_OPTFAIL),
# flat/constant data (STATUS_FLAT), or threw errors (STATUS_ERROR) are ineligible.
_is_inference_eligible(status::UInt8) = status <= STATUS_ROUNDOFF

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
    detail =
        "The cluster/TFCE null needs the z-map of every permutation ($n_channels channels × $n_timepoints time points × " *
        "$n_flips permutations × $n_tested tested coefficients = $(gb(bytes)) GB, $(gb(free)) GB currently free). "
    advice =
        "Reduce `n_permutations`, test fewer coefficients (`tested_coefs`), use a shorter time window or fewer channels, " *
        "or use `correction_method = :maxt`, which does not store the null maps."
    if bytes > ZPERM_ERROR_FRACTION * free
        @minimal_error detail * "This exceeds the available memory. " * advice
    elseif bytes > ZPERM_WARN_FRACTION * free
        @minimal_warning detail * "This is a large share of the available memory. " * advice
    end
    return bytes
end

# Random sign flips give p-values with resolution 1/(n+1); say so when that is too coarse
function _warn_small_n_perms(n_flips::Int, exhaustive::Bool)
    exhaustive && return
    pmin = 1 / (n_flips + 1)
    if pmin >= 0.05
        @minimal_warning "Only $n_flips random sign flips: the smallest attainable corrected p-value is $(round(pmin; digits = 3)), " *
                         "so no point can reach p < 0.05. Use at least 1000 permutations."
    elseif n_flips < SMALL_N_PERMS
        @minimal_warning "Only $n_flips random sign flips: corrected p-values have a resolution of $(round(pmin; digits = 3)) and " *
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
        push!(
            msgs,
            "$n_unfit of $n_total channel × time points could not be fitted (flat or non-finite data: $n_flat; " *
            "fit error: $n_err). Their beta, se, z and p are NaN, and they are excluded from the permutation null and " *
            "from cluster/TFCE formation. Affected: $shown$more.",
        )
        isnothing(first_error) || push!(msgs, "First fit error: $first_error")
    end
    if n_unconv > 0
        push!(
            msgs,
            "$n_unconv of $n_total fits did not converge (evaluation limit reached: $n_max; other optimizer " *
            "failure: $n_opt). Their estimates are kept in beta/se/z_values, and they are excluded from permutation " *
            "inference; see `fit_status` in the result.",
        )
    end
    @minimal_warning join(msgs, "\n")
    return
end

# 1. COEFFICIENT & INTERCEPT RESOLUTION HELPERS
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

    targets = tested_coefs isa Union{AbstractVector,Tuple} ? tested_coefs : (tested_coefs,)
    idxs = Int[]

    for t in targets
        if t isa Integer
            (1 <= t <= n) || @minimal_error "Index $t out of bounds (1:$n)"
            push!(idxs, t)
        elseif t isa Union{AbstractString,Symbol}
            matching = _match_coef_names(coef_names, String(t))
            isempty(matching) && @minimal_error "No model coefficient matches '$t'. Available: $coef_names"
            append!(idxs, matching)
        else
            @minimal_error "Unsupported type for tested_coefs: $(typeof(t))"
        end
    end

    return _add_intercept!(idxs, coef_names, test_intercept)
end

# 2. HIGH-LEVEL WRAPPER FOR EpochData
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

function _resolve_epoch_participant_label(ep::EegFun.EpochData)::String
    if !isempty(ep.data) && !isempty(ep.data[1]) && hasproperty(ep.data[1], :participant)
        return string(ep.data[1][1, :participant])
    end
    stem, _ = splitext(basename(ep.file))
    any(isdigit, stem) && return string(EegFun._extract_participant_id(basename(ep.file)))
    return isempty(stem) ? "1" : stem
end

# 3. EXCHANGEABILITY BLOCKS & PIVOTAL CR2 SIGN-FLIP INFERENCE
#
# Under H0, whole-block sign flipping of the reduced-model (Freedman–Lane) residuals preserves
# their joint distribution when blocks are independent and block-level errors are symmetric.
# For multi-pixel FWER control under random slopes, boundary fits (θ̂ = 0), between-subject
# predictors, and unequal trial counts without refitting, each flip η ∈ {-1, +1}^S (and the
# observed identity flip η = 1) is studentized by its Bell–McCaffrey / Pustejovsky–Tipton CR2
# cluster-robust standard error and mapped to a Normal-equivalent z-score via Satterthwaite df.
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
                       "Use `n_permutations = 0` for Wald-only inference."
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
        for k = 0:(n_flips-1), s = 2:n_blocks
            ((k >> (s - 2)) & 1) == 1 && (H[s, k+1] = -1.0)
        end
        return H, true
    end
    # Monte Carlo sampling: draw independent sign vectors from {-1, 1}^S, excluding ±ones(S)
    # to avoid double-counting the observed identity equivalence class under Phipson & Smyth (2010).
    H = rand(rng, [-1.0, 1.0], n_blocks, n_perms)
    for col = 1:n_perms
        while all(==(1.0), @view H[:, col]) || all(==(-1.0), @view H[:, col])
            H[:, col] .= rand(rng, [-1.0, 1.0], n_blocks)
        end
    end
    return H, false
end

# Between-Within degrees of freedom per fixed-effect column for parametric t-inference (`p_parametric`).
# Columns that are constant within levels of a grouping factor A (between-group predictors, intercept)
# receive G_A - p_between(A); columns that vary within levels of A and lie in A's random-effects span
# (random slopes) receive G_A - max(1, p_between(A)); pure fixed within-group columns receive N - G_max - p_within.
function _compute_bw_dfs(m::LinearMixedModel, X_est::Matrix{Float64}, est_cols::Vector{Int}, n_coefs::Int)
    dfs = fill(NaN, n_coefs)
    p = size(X_est, 2)
    n_obs = size(X_est, 1)
    p == 0 && return dfs

    if isempty(m.reterms)
        df_ols = max(1.0, Float64(n_obs - p))
        for (j, col_idx) in enumerate(est_cols)
            dfs[col_idx] = df_ols
        end
        return dfs
    end

    col_group_df = fill(Inf, p)
    for A in m.reterms
        G_A = length(A.levels)
        level_rows = [Int[] for _ = 1:G_A]
        for (obs, l) in enumerate(A.refs)
            push!(level_rows[l], obs)
        end

        is_const = falses(p)
        in_re_span = falses(p)
        for j = 1:p
            xj = view(X_est, :, j)
            tol = 1e-10 * max(maximum(abs, xj), 1.0)
            c_flag = true
            for rows_l in level_rows
                isempty(rows_l) && continue
                v0 = xj[first(rows_l)]
                if any(r -> abs(xj[r] - v0) > tol, rows_l)
                    c_flag = false
                    break
                end
            end
            is_const[j] = c_flag
            if !c_flag
                span_flag = true
                for rows_l in level_rows
                    isempty(rows_l) && continue
                    yl = xj[rows_l]
                    Zl = Matrix(view(A.z, :, rows_l)')
                    resid = yl - Zl * (Zl \ yl)
                    if norm(resid) > 1e-8 * max(norm(yl), 1.0)
                        span_flag = false
                        break
                    end
                end
                in_re_span[j] = span_flag
            end
        end

        p_between_A = count(is_const)
        for j = 1:p
            if is_const[j]
                col_group_df[j] = min(col_group_df[j], max(1.0, Float64(G_A - p_between_A)))
            elseif in_re_span[j]
                col_group_df[j] = min(col_group_df[j], max(1.0, Float64(G_A - max(1, p_between_A))))
            end
        end
    end

    G_max = maximum(length(A.levels) for A in m.reterms)
    p_within = count(isinf, col_group_df)
    df_within = max(1.0, Float64(n_obs - G_max - p_within))

    for (j, col_idx) in enumerate(est_cols)
        dfs[col_idx] = isfinite(col_group_df[j]) ? col_group_df[j] : df_within
    end
    return dfs
end

@inline function _t_to_z(t::Float64, df::Float64)
    (t == 0.0 || !isfinite(t)) && return t
    lp = logccdf(TDist(df), abs(t))
    lp == -Inf && return copysign(38.0, t)
    lp >= -0.6931471805599453 && return 0.0
    return copysign(invlogccdf(Normal(), lp), t)
end

function _make_cr2_workspace(m::LinearMixedModel, design::BlockDesign, p::Int, n_flips::Int)
    S = length(design.rows)
    lams = Matrix{Float64}[Matrix{Float64}(undef, size(A.λ)) for A in m.reterms]
    XtVX_s = [Matrix{Float64}(undef, p, p) for _ = 1:S]
    XtVy_s = [Vector{Float64}(undef, p) for _ = 1:S]
    Λ_s = [zeros(Float64, size(design.Z[s], 2), size(design.Z[s], 2)) for s = 1:S]
    Zty_s = [Vector{Float64}(undef, size(design.Z[s], 2)) for s = 1:S]
    Gty_s = [Vector{Float64}(undef, size(design.Z[s], 2)) for s = 1:S]
    GtX_s = [Matrix{Float64}(undef, size(design.Z[s], 2), p) for s = 1:S]
    FGX_s = [Matrix{Float64}(undef, size(design.Z[s], 2), p) for s = 1:S]
    ZΛ_s = [Matrix{Float64}(undef, size(design.Z[s], 2), size(design.Z[s], 2)) for s = 1:S]
    C_s = [Matrix{Float64}(undef, size(design.Z[s], 2), size(design.Z[s], 2)) for s = 1:S]
    return Cr2Workspace(
        lams,
        XtVX_s,
        XtVy_s,
        Λ_s,
        Zty_s,
        Gty_s,
        GtX_s,
        FGX_s,
        ZΛ_s,
        C_s,
        Vector{Float64}(undef, p),
        Matrix{Float64}(undef, p, p),
        Vector{Float64}(undef, p),
        Matrix{Float64}(undef, p, p),
        Matrix{Float64}(undef, p, p),
        Matrix{Float64}(undef, p, p),
        [Matrix{Float64}(undef, p, p) for _ = 1:S],
        [Matrix{Float64}(undef, p, p) for _ = 1:S],
        Matrix{Float64}(undef, p, p),
        Matrix{Float64}(undef, p, p),
        Matrix{Float64}(undef, p, S),
        Matrix{Float64}(undef, S, p),
        Matrix{Float64}(undef, S, p),
        Matrix{Float64}(undef, p, p),
        Matrix{Float64}(undef, p, S),
        Matrix{Float64}(undef, p, n_flips),
        Matrix{Float64}(undef, S, n_flips),
        Vector{Float64}(undef, S),
        Vector{Float64}(undef, S),
        Vector{Float64}(undef, p),
        Vector{Float64}(undef, p),
        Vector{Float64}(undef, p),
        Vector{Float64}(undef, p),
    )
end

# Per-block GLS sufficient statistics at the model's current θ, with V_s = I + Z_s Λ_s Λ_s' Z_s':
#   XtVX_s = X_s' V_s⁻¹ X_s,  XtVy_s = X_s' V_s⁻¹ y_s   (Woodbury)
function _block_gls_stats!(ws::Cr2Workspace, design::BlockDesign, lams::Vector{Matrix{Float64}}, y::AbstractVector{Float64})
    n_blocks = length(design.rows)
    for s = 1:n_blocks
        ys = view(y, design.rows[s])
        Λ = ws.Λ_s[s]
        fill!(Λ, 0.0)
        for (k, off) in design.lambda_blocks[s]
            λ = lams[k]
            qk = size(λ, 1)
            @inbounds for j = 1:qk, i = 1:qk
                Λ[off+i, off+j] = λ[i, j]
            end
        end
        mul!(ws.Xty, design.X[s]', ys)
        mul!(ws.Zty_s[s], design.Z[s]', ys)
        mul!(ws.Gty_s[s], Λ', ws.Zty_s[s])
        GtX = ws.GtX_s[s]
        mul!(GtX, Λ', design.ZtX[s])
        mul!(ws.ZΛ_s[s], design.ZtZ[s], Λ)
        C_mat = ws.C_s[s]
        mul!(C_mat, Λ', ws.ZΛ_s[s])
        @inbounds for i = 1:size(C_mat, 1)
            C_mat[i, i] += 1.0
        end
        C = cholesky!(Symmetric(C_mat))
        FGX = ws.FGX_s[s]
        copyto!(FGX, GtX)
        ldiv!(C, FGX)
        copyto!(ws.XtVX_s[s], design.XtX[s])
        mul!(ws.XtVX_s[s], GtX', FGX, -1.0, 1.0)
        copyto!(ws.XtVy_s[s], ws.Xty)
        mul!(ws.XtVy_s[s], FGX', ws.Gty_s[s], -1.0, 1.0)
    end
    return ws.XtVX_s, ws.XtVy_s
end

function _block_gls_stats(design::BlockDesign, lams::Vector{Matrix{Float64}}, y::AbstractVector{Float64})
    n_blocks = length(design.rows)
    XtVX = Vector{Matrix{Float64}}(undef, n_blocks)
    XtVy = Vector{Vector{Float64}}(undef, n_blocks)
    for s = 1:n_blocks
        ys = view(y, design.rows[s])
        q = size(design.Z[s], 2)
        Λ = zeros(q, q)
        for (k, off) in design.lambda_blocks[s]
            λ = lams[k]
            qk = size(λ, 1)
            Λ[(off+1):(off+qk), (off+1):(off+qk)] .= λ
        end
        Xty = design.X[s]' * ys
        Gty = Λ' * (design.Z[s]' * ys)
        GtX = Λ' * design.ZtX[s]
        C = cholesky!(Symmetric(Λ' * design.ZtZ[s] * Λ + I))
        FGX = C \ GtX
        XtVX[s] = design.XtX[s] - GtX' * FGX
        XtVy[s] = Xty - FGX' * Gty
    end
    return XtVX, XtVy
end

# Closed-form Bell–McCaffrey / Pustejovsky–Tipton CR2 cluster-robust pivotal z-statistics
# for coefficient position `pc` at the observed identity flip (returned as scalar `z_obs`) and
# under every sign pattern in `H` (written to `out_null`).
function _cr2_pivotal_z!(
    out_null::AbstractVector{Float64},
    XtVX_s::Vector{Matrix{Float64}},
    XtVy_s::Vector{Vector{Float64}},
    pc::Int,
    H::Matrix{Float64};
    Hs_pre::Union{Nothing,Vector{Matrix{Float64}}} = nothing,
    As_pre::Union{Nothing,Vector{Matrix{Float64}}} = nothing,
    M_pre::Union{Nothing,Matrix{Float64}} = nothing,
    R_pre::Union{Nothing,Matrix{Float64}} = nothing,
    XtVX_pre::Union{Nothing,Matrix{Float64}} = nothing,
    XtVy_pre::Union{Nothing,Vector{Float64}} = nothing,
    cr2_ws::Union{Nothing,Cr2Workspace} = nothing,
)
    p = size(first(XtVX_s), 1)
    S = length(XtVy_s)
    n_flips = size(H, 2)
    XtVX = isnothing(XtVX_pre) ? sum(XtVX_s) : XtVX_pre
    XtVy = isnothing(XtVy_pre) ? sum(XtVy_s) : XtVy_pre

    if isnothing(M_pre) || isnothing(R_pre)
        chol_X = cholesky(Symmetric(XtVX); check = false)
        issuccess(chol_X) || (fill!(out_null, NaN); return NaN)
        M = inv(chol_X)
        R = Matrix(chol_X.U) \ Matrix{Float64}(I, p, p)
    else
        M = M_pre
        R = R_pre
    end
    M[pc, pc] > 0 || (fill!(out_null, NaN); return NaN)

    has_ws = !isnothing(cr2_ws)
    W = has_ws ? cr2_ws.W : Matrix{Float64}(undef, p, S)

    if p == 1
        @inbounds for s = 1:S
            W[1, s] = XtVy_s[s][1]
        end
    else
        idx = [1:(pc-1); (pc+1):p]
        chol_r = cholesky(Symmetric(XtVX[idx, idx]); check = false)
        issuccess(chol_r) || (fill!(out_null, NaN); return NaN)
        γ = chol_r \ XtVy[idx]
        @inbounds for s = 1:S
            Xs = XtVX_s[s]
            ys_v = XtVy_s[s]
            for i = 1:p
                acc = ys_v[i]
                for (k, col_k) in enumerate(idx)
                    acc -= Xs[i, col_k] * γ[k]
                end
                W[i, s] = acc
            end
        end
    end

    a = @view R[pc, :]
    psi = has_ws ? cr2_ws.psi : Vector{Float64}(undef, S)
    B = has_ws ? cr2_ws.B : Matrix{Float64}(undef, S, p)
    d = has_ws ? cr2_ws.d : Vector{Float64}(undef, S)
    F = has_ws ? cr2_ws.F : Matrix{Float64}(undef, S, p)
    FtF = has_ws ? cr2_ws.FtF : Matrix{Float64}(undef, p, p)
    cs = has_ws ? cr2_ws.cs : Vector{Float64}(undef, p)
    fs = has_ws ? cr2_ws.fs : Vector{Float64}(undef, p)
    vs = has_ws ? cr2_ws.vs : Vector{Float64}(undef, p)

    for s = 1:S
        if isnothing(Hs_pre) || isnothing(As_pre)
            Hs = Matrix(Symmetric(R' * XtVX_s[s] * R))
            eig = eigen(Symmetric(Hs))
            # Moore–Penrose pseudoinverse (I_p - H_s)^(+1/2): zero out unit-leverage directions (ev ≈ 1)
            # rather than amplifying floating-point roundoff on single-block predictors.
            inv_sqrt = [ev < 1.0 - 1e-8 ? 1.0 / sqrt(1.0 - ev) : 0.0 for ev in eig.values]
            As = eig.vectors * Diagonal(inv_sqrt) * eig.vectors'
        else
            Hs = Hs_pre[s]
            As = As_pre[s]
        end
        mul!(cs, As, a)
        mul!(fs, Hs, cs)
        d[s] = dot(cs, fs)
        @inbounds for j = 1:p
            F[s, j] = fs[j]
        end
        mul!(vs, R, cs)
        psi[s] = dot(vs, view(W, :, s))
        mul!(cs, XtVX_s[s], vs)
        @inbounds for j = 1:p
            B[s, j] = cs[j]
        end
    end

    tr1 = 0.0
    @inbounds for s = 1:S
        f_ss = 0.0
        for j = 1:p
            v = F[s, j]
            f_ss += v * v
        end
        tr1 += max(0.0, d[s] - f_ss)
    end
    mul!(FtF, F', F)
    cross_sum = 0.0
    @inbounds for s = 1:S
        f_ss = 0.0
        for j = 1:p
            v = F[s, j]
            f_ss += v * v
        end
        cross_sum += d[s] * f_ss
    end
    tr2 = sum(abs2, d) - 2.0 * cross_sum + sum(abs2, FtF)
    (tr1 >= (1.0 - 1e-5) * M[pc, pc] && tr2 > 0 && isfinite(tr1) && isfinite(tr2)) || (fill!(out_null, NaN); return NaN)
    raw_df = tr1^2 / tr2
    raw_df >= 1.0 - 1e-6 || (fill!(out_null, NaN); return NaN)
    df_bm = clamp(raw_df, 1.0, Float64(max(1, S - 1)))

    Q = has_ws ? cr2_ws.Q : Matrix{Float64}(undef, p, S)
    mul!(Q, M, W)
    q_obs = has_ws ? cr2_ws.q_obs : Vector{Float64}(undef, p)
    fill!(q_obs, 0.0)
    @inbounds for s = 1:S, i = 1:p
        q_obs[i] += Q[i, s]
    end
    ss_obs = 0.0
    @inbounds for s = 1:S
        gs = psi[s] - dot(view(B, s, :), q_obs)
        ss_obs += gs * gs
    end
    se_obs = sqrt(ss_obs)
    (se_obs > 0 && isfinite(se_obs)) || (fill!(out_null, NaN); return NaN)
    z_obs = _t_to_z(q_obs[pc] / se_obs, df_bm)
    isfinite(z_obs) || (fill!(out_null, NaN); return NaN)

    if n_flips > 0
        Q_eta = (has_ws && size(cr2_ws.Q_eta, 2) == n_flips) ? cr2_ws.Q_eta : Matrix{Float64}(undef, p, n_flips)
        G = (has_ws && size(cr2_ws.G, 2) == n_flips) ? cr2_ws.G : Matrix{Float64}(undef, S, n_flips)
        mul!(Q_eta, Q, H)
        mul!(G, B, Q_eta)
        @inbounds for col = 1:n_flips
            ss = 0.0
            for s = 1:S
                g = psi[s] * H[s, col] - G[s, col]
                ss += g * g
            end
            se_col = sqrt(ss)
            if se_col > 0 && isfinite(se_col)
                out_null[col] = _t_to_z(Q_eta[pc, col] / se_col, df_bm)
            else
                out_null[col] = NaN
            end
        end
    end
    return z_obs
end

# 4. MODEL DESIGN SETUP
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
    dropped_coefs = sort(MixedModels.pivot(m_initial)[(rank_X+1):end])
    if !isempty(dropped_coefs)
        @minimal_warning "Fixed-effects design matrix is rank deficient (rank $rank_X of $n_coefs). " *
                         "Aliased coefficient(s) $(join(coef_names[dropped_coefs], ", ")) cannot be estimated: " *
                         "their beta, se, z and p values are NaN and they are excluded from permutation testing. " *
                         "Common causes: a predictor that is constant in the data, an empty design cell combined " *
                         "with an interaction, or exactly collinear predictors."
    end

    tested = setdiff(_resolve_tested_coefs(coef_names, tested_coefs, test_intercept), dropped_coefs)
    return m_initial, coef_names, dropped_coefs, tested, meta_df_init
end

function _create_scratch(f, meta_df_init, est_cols, block_ids, n_flips, n_tested, contrasts)
    m = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f, meta_df_init; contrasts = contrasts)
    end
    m.optsum.ftol_rel = OBSERVED_FTOL_REL
    m.optsum.maxfeval = MAX_EVALS
    X = Matrix{Float64}(modelmatrix(m)[:, invperm(MixedModels.pivot(m))][:, est_cols])
    design = isnothing(block_ids) ? nothing : BlockDesign(m, X, block_ids)
    cr2_ws = (!isnothing(design) && n_flips > 0) ? _make_cr2_workspace(m, design, length(est_cols), n_flips) : nothing
    n_obs = length(m.y)
    n_coefs = length(coefnames(m))
    return LmmScratch(
        m,
        design,
        cr2_ws,
        zeros(n_obs),
        zeros(n_coefs),
        zeros(n_coefs),
        zeros(n_coefs),
        zeros(n_coefs),
        fill(NaN, n_coefs),
        zeros(n_flips, n_tested),
        Ref(false),
    )
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
    has_df = hasproperty(ctx, :df_bw)
    @inbounds for k in eachindex(scr.se)
        z = scr.beta[k] / scr.se[k]
        scr.z[k] = z
        if has_df
            df_k = ctx.df_bw[k]
            scr.p[k] = isfinite(df_k) ? 2.0 * ccdf(TDist(df_k), abs(z)) : NaN
        else
            scr.p[k] = 2.0 * ccdf(Normal(), abs(z))
        end
    end
    # Aliased (dropped) coefficients are NaN by design; every estimable one must be finite
    for k in ctx.est_cols
        (isfinite(scr.beta[k]) && isfinite(scr.z[k])) || return STATUS_ERROR
    end

    ctx.n_flips == 0 && return status

    ws = scr.cr2_ws::Cr2Workspace
    for (k, A) in enumerate(m.reterms)
        copyto!(ws.lams[k], A.λ)
    end
    XtVX_s, XtVy_s = _block_gls_stats!(ws, scr.design::BlockDesign, ws.lams, scr.y)
    p = length(ctx.est_cols)
    S = length(XtVy_s)
    XtVX = ws.XtVX
    XtVy = ws.XtVy
    fill!(XtVX, 0.0)
    fill!(XtVy, 0.0)
    @inbounds for s = 1:S
        Xs = XtVX_s[s]
        ys = XtVy_s[s]
        for j = 1:p
            XtVy[j] += ys[j]
            for i = 1:p
                XtVX[i, j] += Xs[i, j]
            end
        end
    end

    copyto!(ws.chol_buf, XtVX)
    chol_X = cholesky!(Symmetric(ws.chol_buf); check = false)
    issuccess(chol_X) || return STATUS_ERROR
    M = ws.M
    R = ws.R
    fill!(M, 0.0)
    fill!(R, 0.0)
    @inbounds for i = 1:p
        M[i, i] = 1.0
        R[i, i] = 1.0
    end
    ldiv!(chol_X, M)
    ldiv!(chol_X.U, R)

    Hs_pre = ws.Hs_pre
    As_pre = ws.As_pre
    tmp_pp = ws.tmp_pp
    tmp_pp2 = ws.tmp_pp2
    for s = 1:S
        Hs = Hs_pre[s]
        mul!(tmp_pp, XtVX_s[s], R)
        mul!(Hs, R', tmp_pp)
        @inbounds for j = 1:p, i = 1:(j-1)
            avg_ij = 0.5 * (Hs[i, j] + Hs[j, i])
            Hs[i, j] = avg_ij
            Hs[j, i] = avg_ij
        end
        copyto!(tmp_pp, Hs)
        eig = eigen!(Symmetric(tmp_pp))
        V = eig.vectors
        ev = eig.values
        @inbounds for k = 1:p
            w_k = ev[k] < 1.0 - 1e-8 ? (1.0 / sqrt(1.0 - ev[k])) : 0.0
            for i = 1:p
                tmp_pp2[i, k] = V[i, k] * w_k
            end
        end
        mul!(As_pre[s], tmp_pp2, V')
    end

    fill!(scr.z_perm_obs, NaN)
    for (j, pc) in enumerate(ctx.tested_pos)
        zn = view(scr.znull, :, j)
        z_obs = _cr2_pivotal_z!(
            zn,
            XtVX_s,
            XtVy_s,
            pc,
            ctx.H;
            Hs_pre = Hs_pre,
            As_pre = As_pre,
            M_pre = M,
            R_pre = R,
            XtVX_pre = XtVX,
            XtVy_pre = XtVy,
            cr2_ws = ws,
        )
        (isfinite(z_obs) && all(isfinite, zn)) || return STATUS_ERROR
        scr.z_perm_obs[ctx.tested[j]] = z_obs
    end
    return status
end

# Copy a successfully computed pixel into the shared outputs
function _commit_pixel!(scr::LmmScratch, out, c_idx::Int, t_idx::Int, slot::Int, ctx, eligible::Bool)
    out.singular[c_idx, t_idx] = scr.singular[]
    @inbounds for k in eachindex(scr.beta)
        out.beta[c_idx, t_idx, k] = scr.beta[k]
        out.se[c_idx, t_idx, k] = scr.se[k]
        out.z[c_idx, t_idx, k] = scr.z[k]
        out.p[c_idx, t_idx, k] = scr.p[k]
    end
    (ctx.n_flips == 0 || !eligible) && return
    @inbounds for (j, c) in enumerate(ctx.tested)
        out.z_perm_obs[c_idx, t_idx, c] = scr.z_perm_obs[c]
        if ctx.store_maps
            for p = 1:ctx.n_flips
                out.zperm[c_idx, t_idx, p, j] = Float32(scr.znull[p, j])
            end
        else
            for p = 1:ctx.n_flips
                a = abs(scr.znull[p, j])
                a > out.maxabs[slot, p, j] && (out.maxabs[slot, p, j] = a)
            end
        end
    end
    return
end

function _process_pixel!(scr::LmmScratch, eeg, c_idx::Int, t_idx::Int, slot::Int, out, ctx)
    scr.y .= @view eeg[c_idx, t_idx, :]
    status = try
        _compute_pixel!(scr, ctx)
    catch e
        e isa Union{InterruptException,UndefVarError,BoundsError,TypeError} && rethrow()
        STATUS_ERROR
    end
    eligible = _is_inference_eligible(status)
    status < STATUS_FLAT && _commit_pixel!(scr, out, c_idx, t_idx, slot, ctx, eligible)
    out.status[c_idx, t_idx] = status
    return
end

function _run_pixel_chunks!(scratches, task_chunks, eeg::Array{Float64,3}, out, base_ctx, prog, progress::Bool)
    if length(task_chunks) == 1
        scr1 = scratches[1]
        for (c_idx, t_idx) in task_chunks[1]
            _process_pixel!(scr1, eeg, c_idx, t_idx, 1, out, base_ctx)
            progress && next!(prog)
        end
    else
        Threads.@threads :static for slot in eachindex(task_chunks)
            local scr = scratches[slot]
            for (c_idx, t_idx) in task_chunks[slot]
                _process_pixel!(scr, eeg, c_idx, t_idx, slot, out, base_ctx)
                progress && next!(prog)
            end
        end
    end
    return
end

# Re-run the first pixel that errored to report the diagnostic reason
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

# 6. NULL DISTRIBUTIONS & CORRECTED P-VALUES
function _build_cluster_stat_map(t_map, channel_names, time_points, spatial_conn, cluster_threshold)
    n_channels, n_timepoints = size(t_map)
    pos_clusters, neg_clusters = EegFun._find_clusters(
        t_map .> cluster_threshold,
        t_map .< -cluster_threshold,
        channel_names,
        time_points,
        spatial_conn,
        :spatiotemporal,
    )
    pos_stats = EegFun._compute_cluster_statistics(pos_clusters, t_map; return_clusters = false)
    neg_stats = EegFun._compute_cluster_statistics(neg_clusters, t_map; return_clusters = false)

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
    chunks = collect(Iterators.partition(1:n_flips, cld(n_flips, max(1, Threads.nthreads()))))
    Threads.@threads :static for chunk_idx in eachindex(chunks)
        flip_chunk = chunks[chunk_idx]
        labels_buf = zeros(Int, n_channels, n_timepoints)
        queue_buf = Vector{Tuple{Int,Int}}(undef, n_channels * n_timepoints)
        pos_mask = BitArray{2}(undef, n_channels, n_timepoints)
        neg_mask = BitArray{2}(undef, n_channels, n_timepoints)
        t_map = Array{Float64,2}(undef, n_channels, n_timepoints)
        tfce_buf = ctx.correction === :tfce ? zeros(Float64, n_channels, n_timepoints) : Array{Float64,2}(undef, 0, 0)
        for p in flip_chunk, (j, c) in enumerate(tested)
            copyto!(t_map, view(zperm, :, :, p, j))
            if ctx.correction === :tfce
                EegFun._compute_tfce!(
                    tfce_buf,
                    t_map,
                    ctx.spatial_conn,
                    :spatiotemporal,
                    labels_buf,
                    queue_buf;
                    E = ctx.tfce_E,
                    H = ctx.tfce_H,
                    dh = ctx.tfce_dh,
                )
                null_max[p, c] = maximum(abs, tfce_buf; init = 0.0)
            else
                @inbounds for i in eachindex(t_map)
                    pos_mask[i] = t_map[i] > ctx.cluster_threshold
                    neg_mask[i] = t_map[i] < -ctx.cluster_threshold
                end
                max_pos =
                    EegFun._max_cluster_stat(pos_mask, t_map, ctx.spatial_conn, :spatiotemporal, labels_buf, queue_buf; is_negative = false)
                min_neg =
                    EegFun._max_cluster_stat(neg_mask, t_map, ctx.spatial_conn, :spatiotemporal, labels_buf, queue_buf; is_negative = true)
                null_max[p, c] = max(max_pos, abs(min_neg))
            end
        end
    end
    return null_max
end

function _compute_p_corrected(z_obs_map, null_max, tested, ctx, bad::AbstractMatrix{Bool})
    n_channels, n_timepoints, n_coefs = size(z_obs_map)
    # Untested coefficients have no permutation p-value: NaN (not 1.0, which reads as "tested, not significant")
    p_corr = fill(NaN, n_channels, n_timepoints, n_coefs)
    for c in tested
        null_sorted = sort(null_max[:, c])
        zc = copy(z_obs_map[:, :, c])
        zc[bad] .= 0.0   # ineligible/unfitted points carry no evidence and never form clusters
        stat_map = if ctx.correction === :maxt
            abs.(zc)
        elseif ctx.correction === :tfce
            abs.(
                EegFun._compute_tfce(
                    zc,
                    ctx.channel_names,
                    ctx.time_points,
                    ctx.spatial_conn,
                    :spatiotemporal;
                    E = ctx.tfce_E,
                    H = ctx.tfce_H,
                    dh = ctx.tfce_dh,
                ),
            )
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

function _fit_mass_lmm(
    eeg_data,
    meta_df,
    f,
    epochs;
    n_permutations::Integer = 1000,
    n_perms::Union{Nothing,Integer} = nothing,
    tested_coefs = nothing,
    test_intercept::Bool = false,
    flip_by = nothing,
    correction_method::Symbol = :maxt,
    correction::Union{Nothing,Symbol} = nothing,
    cluster_threshold::Real = 2.0,
    tfce_E::Real = 0.5,
    tfce_H::Real = 2.0,
    tfce_dh::Real = 0.1,
    spatial_connectivity = nothing,
    channel_names = nothing,
    time_points = nothing,
    REML::Bool = false,
    contrasts::AbstractDict = Dict{Symbol,Any}(),
    show_progress::Bool = true,
    progress::Union{Nothing,Bool} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    eff_n_perms = isnothing(n_perms) ? n_permutations : n_perms
    eff_corr = isnothing(correction) ? correction_method : correction
    eff_progress = isnothing(progress) ? show_progress : progress
    eff_corr in (:maxt, :cluster, :tfce) || @minimal_error "correction_method must be :maxt, :cluster or :tfce (got :$eff_corr)"
    eff_n_perms >= 0 || @minimal_error "n_permutations must be ≥ 0"
    n_channels, n_timepoints, n_epochs = size(eeg_data)
    nrow(meta_df) == n_epochs ||
        @minimal_error "Number of rows in meta_df ($(nrow(meta_df))) must match number of epochs in eeg_data ($n_epochs)."
    if !isnothing(channel_names)
        length(channel_names) == n_channels ||
            @minimal_error "length(channel_names) ($(length(channel_names))) must match n_channels ($n_channels)."
    end
    if !isnothing(time_points)
        length(time_points) == n_timepoints ||
            @minimal_error "length(time_points) ($(length(time_points))) must match n_timepoints ($n_timepoints)."
        all(isfinite, time_points) || @minimal_error "time_points must contain finite values."
    end
    channel_names = isnothing(channel_names) ? [Symbol("Ch$i") for i = 1:n_channels] : Symbol.(channel_names)
    time_points = Float64.(isnothing(time_points) ? (1:n_timepoints) : time_points)

    # Many tiny BLAS calls per model: single-threaded BLAS avoids oversubscription across Julia threads.
    blas_threads = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    try
        # 1. Design, tested coefficients, Between-Within parametric df, exchangeability blocks
        m0, coef_names, dropped, tested, meta_df_init = _setup_lmm_design(f, meta_df, n_epochs, tested_coefs, test_intercept, contrasts)
        n_coefs = length(coef_names)
        est_cols = setdiff(1:n_coefs, dropped)
        X0_est = Matrix{Float64}(modelmatrix(m0)[:, invperm(MixedModels.pivot(m0))][:, est_cols])
        df_bw = _compute_bw_dfs(m0, X0_est, est_cols, n_coefs)

        do_perm = eff_n_perms > 0 && !isempty(tested)
        do_perm &&
            eff_corr !== :maxt &&
            isnothing(spatial_connectivity) &&
            n_channels > 1 &&
            @minimal_warning "No spatial connectivity or layout available: clusters/TFCE will form across time only."

        flip_col, block_ids, H, exhaustive = nothing, nothing, zeros(0, 0), false
        if do_perm
            flip_col, block_ids = _resolve_flip_by(m0, meta_df, flip_by)
            n_blocks = maximum(block_ids)
            2.0^(n_blocks - 1) >= MIN_SIGN_FLIPS || @minimal_error "Only $n_blocks levels of :$flip_col: whole-block sign flipping " *
                           "has $(2^(n_blocks - 1)) distinct patterns, too few to reach p < 0.05. At least 6 blocks are required."
            n_blocks < 8 &&
                @minimal_warning "Only $n_blocks levels of :$flip_col: the permutation distribution is coarse (≤ $(2^(n_blocks - 1)) distinct values)."

            # Check OLS block leverage up front: any coefficient identified from a single block
            # (or with < 1 effective CR2 degree of freedom across blocks) cannot be block-permuted.
            p_est = length(est_cols)
            rows_b = [findall(==(b), block_ids) for b = 1:n_blocks]
            XtX_s = [X0_est[r, :]' * X0_est[r, :] for r in rows_b]
            chol_0 = cholesky(Symmetric(sum(XtX_s)); check = false)
            if issuccess(chol_0)
                M0 = inv(chol_0)
                R0 = Matrix(chol_0.U) \ Matrix{Float64}(I, p_est, p_est)
                Hs0 = [Matrix(Symmetric(R0' * XtX_s[s] * R0)) for s = 1:n_blocks]
                As0 = map(Hs0) do Hs
                    eig = eigen(Symmetric(Hs))
                    w = [ev < 1.0 - 1e-8 ? 1.0 / sqrt(1.0 - ev) : 0.0 for ev in eig.values]
                    eig.vectors * Diagonal(w) * eig.vectors'
                end
                untestable = Int[]
                for c in tested
                    pc = findfirst(==(c), est_cols)::Int
                    a = R0[pc, :]
                    d0 = zeros(n_blocks)
                    F0 = zeros(n_blocks, p_est)
                    for s = 1:n_blocks
                        cs = As0[s] * a
                        fs = Hs0[s] * cs
                        d0[s] = dot(cs, fs)
                        F0[s, :] .= fs
                    end
                    tr1 = sum(max(0.0, d0[s] - sum(abs2, view(F0, s, :))) for s = 1:n_blocks)
                    FtF = F0' * F0
                    tr2 = sum(abs2, d0) - 2.0 * sum(d0[s] * sum(abs2, view(F0, s, :)) for s = 1:n_blocks) + sum(abs2, FtF)
                    if !(tr1 >= (1.0 - 1e-5) * M0[pc, pc] && tr2 > 0 && tr1^2 / tr2 >= 2.0 - 1e-6)
                        push!(untestable, c)
                    end
                end
                if !isempty(untestable)
                    @minimal_warning "Coefficient(s) $(join(coef_names[untestable], ", ")) have fewer than 2 effective across-block " *
                                     "degrees of freedom (e.g. identified from only one or two levels of :$flip_col) and are excluded from " *
                                     "permutation testing; their p_corrected is NaN."
                    tested = setdiff(tested, untestable)
                end
            end

            if isempty(tested)
                do_perm = false
            else
                H, exhaustive = _sign_flip_matrix(n_blocks, eff_n_perms, rng)
                _warn_small_n_perms(size(H, 2), exhaustive)
                @info "Permutation inference: whole-block sign flipping by :$flip_col ($n_blocks blocks, " *
                      (exhaustive ? "all $(size(H, 2)) distinct patterns" : "$(size(H, 2)) random patterns") *
                      "), correction = :$eff_corr, testing $(join(coef_names[tested], ", "))"
            end
        end
        n_flips = size(H, 2)
        store_maps = do_perm && eff_corr !== :maxt
        store_maps && _check_zperm_memory(n_channels, n_timepoints, n_flips, length(tested))
        tested_pos = Int[findfirst(==(c), est_cols)::Int for c in tested]

        # 2. Shared-memory outputs and per-thread-chunk scratch models
        tasks = [(c, t) for c = 1:n_channels for t = 1:n_timepoints]
        n_threads = max(1, Threads.nthreads())
        task_chunks = collect(Iterators.partition(tasks, cld(max(1, length(tasks)), n_threads)))
        n_slots = max(1, length(task_chunks))

        out = (
            beta = fill(NaN, n_channels, n_timepoints, n_coefs),
            se = fill(NaN, n_channels, n_timepoints, n_coefs),
            z = fill(NaN, n_channels, n_timepoints, n_coefs),
            p = fill(NaN, n_channels, n_timepoints, n_coefs),
            z_perm_obs = do_perm ? fill(NaN, n_channels, n_timepoints, n_coefs) : zeros(0, 0, 0),
            singular = fill(false, n_channels, n_timepoints),
            status = fill(STATUS_FLAT, n_channels, n_timepoints),
            zperm = store_maps ? zeros(Float32, n_channels, n_timepoints, n_flips, length(tested)) : zeros(Float32, 0, 0, 0, 0),
            maxabs = (do_perm && !store_maps) ? zeros(Float64, n_slots, n_flips, length(tested)) : zeros(0, 0, 0),
        )

        eeg = eeg_data isa Array{Float64,3} ? eeg_data : Array{Float64,3}(eeg_data)
        base_ctx = (
            REML = REML,
            n_flips = do_perm ? n_flips : 0,
            tested = tested,
            tested_pos = tested_pos,
            H = H,
            store_maps = store_maps,
            est_cols = est_cols,
            df_bw = df_bw,
        )

        # 3. Mass-univariate fits across threads
        if eff_progress
            if n_threads > 1
                @info "Fitting $(length(tasks)) mixed models across $n_epochs epochs on $n_threads threads..."
            else
                @info "Fitting $(length(tasks)) mixed models across $n_epochs epochs on 1 thread (start Julia with `julia -t auto` to use multiple threads)..."
            end
        end

        make_scratch() = _create_scratch(f, copy(meta_df_init; copycols = true), est_cols, block_ids, n_flips, length(tested), contrasts)
        scratches = [make_scratch() for _ = 1:n_slots]
        prog = Progress(length(tasks); dt = 1, desc = "Fitting models...", enabled = eff_progress)

        Logging.with_logger(Logging.NullLogger()) do
            _run_pixel_chunks!(scratches, task_chunks, eeg, out, base_ctx, prog, eff_progress)
        end
        finish!(prog)

        # Points that could not be fitted / did not converge: one consolidated warning, NaN statistics already set
        fit_status = out.status
        first_error = _first_fit_error(fit_status, eeg, make_scratch, base_ctx)
        _report_fit_problems(fit_status, channel_names, first_error)
        unfitted = fit_status .>= STATUS_FLAT
        ineligible = .!_is_inference_eligible.(fit_status)

        beta, se, z, p_param = out.beta, out.se, out.z, out.p
        for arr in (beta, se, z, p_param)
            for c = 1:size(arr, 3)
                view(arr, :, :, c)[unfitted] .= NaN
            end
            arr[:, :, dropped] .= NaN
        end

        # 4. Null distributions & corrected p-values
        null_max = fill(NaN, n_flips, n_coefs)
        p_corrected = nothing
        spatial_conn = nothing
        z_perm_obs_out = nothing
        if do_perm
            z_perm_obs_out = out.z_perm_obs
            for c = 1:size(z_perm_obs_out, 3)
                view(z_perm_obs_out, :, :, c)[ineligible] .= NaN
            end
            z_perm_obs_out[:, :, dropped] .= NaN

            spatial_conn =
                isnothing(spatial_connectivity) ? sparse(Int[], Int[], Bool[], n_channels, n_channels) :
                SparseMatrixCSC{Bool}(spatial_connectivity)
            spatial_conn = issymmetric(spatial_conn) ? spatial_conn : (spatial_conn .| copy(transpose(spatial_conn)))
            ctx = (
                correction = eff_corr,
                cluster_threshold = Float64(cluster_threshold),
                tfce_E = tfce_E,
                tfce_H = tfce_H,
                tfce_dh = tfce_dh,
                channel_names = channel_names,
                time_points = time_points,
                spatial_conn = spatial_conn,
                exhaustive = exhaustive,
            )
            if store_maps
                @info eff_corr === :tfce ? "Computing TFCE null distribution..." :
                      "Computing cluster-mass null distribution (threshold = $cluster_threshold)..."
                _map_null!(null_max, out.zperm, tested, ctx)
            else
                for (j, c) in enumerate(tested)
                    null_max[:, c] .= vec(maximum(view(out.maxabs, :, :, j); dims = 1))
                end
            end
            p_corrected = _compute_p_corrected(z_perm_obs_out, null_max, tested, ctx, ineligible)
        end

        inference = (
            correction = do_perm ? eff_corr : :none,
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
            df_parametric = df_bw,
            z_perm_obs = z_perm_obs_out,
            prespecified_field_size = (n_channels, n_timepoints),
            field_mask_policy = :intersection_of_valid_pixels,
            n_valid_pixels = count(_is_inference_eligible, fit_status),
            n_ineligible_pixels = count(!_is_inference_eligible, fit_status),
            n_total_pixels = length(fit_status),
        )
        return EegFun.LmmStatsResult(
            coef_names,
            channel_names,
            time_points,
            beta,
            se,
            z,
            p_param,
            p_corrected,
            null_max,
            tested,
            out.singular,
            fit_status,
            inference,
            epochs,
        )
    finally
        LinearAlgebra.BLAS.set_num_threads(blas_threads)
    end
end

end # module
