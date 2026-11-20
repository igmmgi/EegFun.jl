# Stubs and predictor extraction utilities for mass-univariate linear mixed models.

######################################
########### User Functions ###########
######################################

"""
    fit_mass_lmm(epochs::EpochData, f::FormulaTerm; kwargs...)
    fit_mass_lmm(epochs::AbstractVector{<:EpochData}, f::FormulaTerm; condition_selection = conditions(), kwargs...)
    fit_mass_lmm(eeg_data::AbstractArray{<:Real,3}, meta_df::DataFrame, f::FormulaTerm; kwargs...)

Fit a linear mixed model at every channel × time point and (optionally) control the family-wise
error rate with a whole-block sign-flip permutation test.

When given an `AbstractVector{<:EpochData}` (e.g. from `extract_epochs` or `read_all_data("epochs")`),
trials are stacked across datasets after applying `condition_selection`, automatically populating
`:participant`, `:condition`, `:condition_num`, `:condition_name`, and `:file` metadata columns when
not already present in the trial `DataFrame`s.

**Requires `MixedModels.jl` and `StatsModels.jl` to be loaded.**

# Keyword Arguments
- `n_permutations::Int=1000`: Number of random sign flips (`0` = no permutation inference). If all
  `2^(n_blocks-1)` distinct patterns fit within `n_permutations`, they are enumerated exactly instead.
- `tested_coefs=nothing`: Coefficients to permutation-test.
  - `nothing` or `:effects`: all non-intercept coefficients (default); `:all`: including the intercept.
  - `Symbol`, `String`, `Integer` or a `Vector` of these. An exact coefficient name selects that column;
    a term name (e.g. `:condition`) selects all of its contrast columns but not interactions.
- `test_intercept::Bool=false`: Also test the intercept.
- `flip_by=nothing`: Column defining the exchangeability blocks. Default: the outermost grouping factor.
- `correction_method::Symbol=:maxt`: FWER correction: `:maxt` (max-|z|), `:cluster` (spatio-temporal
  cluster mass) or `:tfce`.
- `cluster_threshold::Real=2.0`: |z| threshold for cluster formation (`:cluster`).
- `tfce_E=0.5`, `tfce_H=2.0`, `tfce_dh=0.1`: TFCE parameters (`:tfce`).
- `spatial_connectivity=nothing`: Channel adjacency for clustering/TFCE (built from the layout for EpochData).
- `REML::Bool=false`: Use REML instead of ML.
- `contrasts=Dict{Symbol,Any}()`: Contrast coding of categorical predictors, passed to MixedModels as in
  `fit(MixedModel, f, df; contrasts = Dict(:condition => EffectsCoding(base = "A")))`. Default: treatment (dummy)
  coding. Coefficient names in `tested_coefs` refer to the resulting columns.
- `show_progress::Bool=true`: Show progress bar and summary log messages.
- `rng=Random.default_rng()`: Random number generator for the sign flips.
- Array method only: `channel_names`, `time_points`, `dim_order=(:channels, :timepoints, :epochs)`.

# Returns
- [`LmmStatsResult`](@ref): Fitted mass-univariate mixed-model estimates, diagnostics, and optional permutation inference.

# Details

## Inference
Each point is fitted once with `MixedModels.jl`, storing the model-based fixed-effect estimates `beta`,
standard errors `se`, Wald statistics `z_values = β̂ / SE`, and uncorrected two-sided parametric p-values
`p_parametric` (evaluated against a t-distribution with Between-Within degrees of freedom stored in
`inference.df_parametric`).

For permutation inference, Freedman–Lane residuals of the reduced model (without the tested coefficient)
are multiplied by a random sign per *block* (by default the outermost random-effects grouping factor,
e.g. subject). At each flip, the test statistic is studentized using a closed-form Bell–McCaffrey (CR2)
cluster-robust standard error at the fitted variance ratios (no refitting) and mapped to a
Normal-equivalent z-score (`inference.z_perm_obs`) via Satterthwaite degrees of freedom. This tests the
population null (zero fixed effect across blocks) for within- and between-block predictors, random slopes,
unequal trial counts, and under-specified random-effects structures, assuming independent blocks and
symmetric block-level errors. Because nuisance regression and Satterthwaite degrees of freedom are
finite-sample approximations, calibration improves with the number of blocks (at least 10–12 blocks are
recommended in practice, especially when trial counts vary strongly across blocks).

All random-effects grouping factors must be nested within the flipping block (e.g. `(1 | subject)`,
`(1 + condition | subject)`, `(1 | subject & session)`). Crossed designs (e.g. subject × item) have no
single-factor block permutation scheme and raise an error when `n_permutations > 0`; use `n_permutations = 0` for
parametric results.

## Parallelism
Fits run across Julia threads (`Threads.@threads`). Start Julia with `julia -t auto` (or set
`JULIA_NUM_THREADS`) to use multiple threads.

## Rank-deficient designs
If the fixed-effects design matrix is rank deficient (e.g. a predictor that is constant in the
data, an empty design cell combined with an interaction, or exactly collinear predictors), the
aliased coefficient(s) cannot be estimated. A single warning names them; their `beta`, `se`,
`z_values`, `p_parametric` and `p_corrected` values are `NaN`, and they are excluded from
permutation testing. All estimable coefficients are fitted and tested as usual.

## Fit failures and convergence
Each channel × time point is fitted independently, and one bad point never aborts the run:
- Points whose data are constant (e.g. a flat or dead channel) or contain non-finite values are not fitted;
  points where the fit throws or yields non-finite estimates are treated the same. Their `beta`, `se`,
  `z_values`, `p_parametric` and `p_corrected` are `NaN`, and they are excluded from the permutation null
  distribution and from cluster/TFCE formation. A single warning lists the affected channels and the first
  error message.
- The optimizer outcome of every fit is stored in `fit_status` (see [`LmmStatsResult`](@ref)). Fits that hit the
  evaluation limit or otherwise failed to converge keep their estimates in `beta`/`se`/`z_values`, are
  excluded from permutation FWER inference, and are reported in one consolidated warning;
  `show(result)` summarises all outcomes.

## Memory and permutation count
`correction_method = :cluster` and `:tfce` store the z-map of every permutation (`n_channels × n_timepoints × n_permutations ×
n_tested` single-precision values). A warning is issued when this exceeds half of the free memory, and an
error before it would exceed ~90%; `:maxt` does not store the maps. With few random sign flips (`n_permutations < 100`)
the p-value resolution is `1/(n_permutations+1)`, and a warning is issued; below 19 no point can reach p < 0.05.
"""
function fit_mass_lmm(args...; kwargs...)
    if isnothing(Base.get_extension(@__MODULE__, :EegFunMixedModelsExt))
        @minimal_error "To use Mass-Univariate Mixed Models, you must first load the packages: `using MixedModels, StatsModels`"
    end
    throw(MethodError(fit_mass_lmm, args))
end

"""
    extract_predictor_stats(result::LmmStatsResult, coef; alpha = 0.05)

Format the permutation-corrected results for one predictor (specified by name, `Symbol`, or integer index)
as a [`PermutationResult`](@ref), so they can be plotted with `plot_erp_stats`, `plot_topography_stats`, etc.

The correction method and its settings (cluster threshold, TFCE parameters, connectivity) are taken
from `result.inference`, i.e. exactly those used to build the null distribution in `fit_mass_lmm`.
The returned `PermutationResult.data` contains `[erp_beta, erp_zero]`, where `erp_beta` holds the
fitted `β̂(t)` waveform per channel (with `se_cond1 = se_diff = SE(β̂)`) and `erp_zero` is zero,
so `plot_erp_stats` displays the `β̂` effect waveform and `plot_topography_stats(...; topo_data = :difference)`
displays `β̂` topographies.

Requires a result from the `EpochData` method of `fit_mass_lmm`.
"""
function extract_predictor_stats(result::LmmStatsResult, coef; alpha = 0.05)
    coef_idx = _resolve_lmm_coef_idx(result, coef)
    coef_str = result.coefficients[coef_idx]
    isnothing(result.p_corrected) && @minimal_error "No permutation inference in this result: run fit_mass_lmm with n_permutations > 0."
    coef_idx in result.tested || @minimal_error "Coefficient '$coef_str' was not permutation-tested. Include it in `tested_coefs` " *
                   "(or set `test_intercept = true`) when calling fit_mass_lmm."
    isnothing(result.epochs) && @minimal_error "extract_predictor_stats requires a result from the EpochData method of fit_mass_lmm."

    inf = result.inference
    n_electrodes = length(result.channels)
    n_time_points = length(result.time_points)
    z_raw = Matrix{Float64}(_lmm_stat_slice(result, coef_idx))
    bad = result.fit_status .> 0x01
    z_clean = copy(z_raw)
    z_clean[bad] .= 0.0
    null_dist = result.null_max[:, coef_idx]
    n_perms = length(null_dist)
    thr = Float64(inf.cluster_threshold)

    pos_clusters, neg_clusters = Cluster[], Cluster[]
    if inf.correction === :cluster
        spatial_connectivity =
            isnothing(inf.spatial_connectivity) ? _build_connectivity_matrix(result.channels, result.epochs.layout, :spatiotemporal) :
            inf.spatial_connectivity
        raw_pos, raw_neg =
            _find_clusters(z_clean .> thr, z_clean .< -thr, result.channels, result.time_points, spatial_connectivity, :spatiotemporal)
        pos_stats = _compute_cluster_statistics(raw_pos, z_clean; return_clusters = false)
        neg_stats = _compute_cluster_statistics(raw_neg, z_clean; return_clusters = false)
        null_sorted = sort(null_dist)
        dims = (n_electrodes, n_time_points)
        pos_clusters, final_mask_pos = _update_clusters(raw_pos, pos_stats, null_sorted, alpha, inf.exhaustive, dims)
        neg_clusters, final_mask_neg = _update_clusters(raw_neg, neg_stats, null_sorted, alpha, inf.exhaustive, dims)
    else
        # :maxt and :tfce give point-wise FWER-corrected p-values (NaN on ineligible pixels -> false)
        sig = BitArray(result.p_corrected[:, :, coef_idx] .<= alpha)
        final_mask_pos = sig .& (z_clean .> 0)
        final_mask_neg = sig .& (z_clean .< 0)
    end

    se_diff = copy(result.se[:, :, coef_idx])
    correction_symbol = inf.correction === :cluster ? :cluster_permutation : inf.correction
    df_coef = (haskey(inf, :df_parametric) && coef_idx <= length(inf.df_parametric)) ? Float64(inf.df_parametric[coef_idx]) : 0.0
    cluster_alpha = 2.0 * Distributions.ccdf(Distributions.Normal(), thr)
    test_info =
        TestInfo(:LMM, df_coef, Float64(alpha), :both, correction_symbol, ClusterInfo(:parametric, :spatiotemporal, n_perms, cluster_alpha))

    df_beta = DataFrame(:time => copy(result.time_points))
    df_zero = DataFrame(:time => copy(result.time_points))
    for (ch_i, ch) in enumerate(result.channels)
        df_beta[!, ch] = Float64[result.beta[ch_i, t_i, coef_idx] for t_i = 1:n_time_points]
        df_zero[!, ch] = zeros(Float64, n_time_points)
    end
    n_ep = length(result.epochs.data)
    erp_beta = ErpData(
        result.epochs.file,
        coef_idx,
        coef_str,
        df_beta,
        copy(result.epochs.layout),
        result.epochs.sample_rate,
        copy(result.epochs.analysis_info),
        n_ep,
    )
    erp_zero = ErpData(
        result.epochs.file,
        0,
        "0",
        df_zero,
        copy(result.epochs.layout),
        result.epochs.sample_rate,
        copy(result.epochs.analysis_info),
        n_ep,
    )
    p_mat = inf.correction === :cluster ? nothing : copy(result.p_corrected[:, :, coef_idx])

    return PermutationResult(
        test_info,
        [erp_beta, erp_zero],
        StatMatrix(z_raw, p_mat),
        Masks(final_mask_pos, final_mask_neg),
        Clusters(pos_clusters, neg_clusters),
        PermutationDistribution(copy(null_dist), copy(null_dist)),
        copy(result.channels),
        copy(result.time_points),
        (thr, -thr),
        se_diff,
        copy(se_diff),
        zeros(n_electrodes, n_time_points),
        copy(se_diff),
    )
end

######################################
######### Internal Functions #########
######################################

function _resolve_lmm_coef_idx(result::LmmStatsResult, coef)
    if coef isa Integer
        (1 <= coef <= length(result.coefficients)) || @minimal_error("coef index $coef out of bounds (1:$(length(result.coefficients)))")
        return Int(coef)
    end
    coef_str = String(coef)
    idx = findfirst(==(coef_str), result.coefficients)
    isnothing(idx) && @minimal_error("Coefficient '$coef_str' not found in model. Available: $(result.coefficients)")
    return idx
end

function _lmm_stat_slice(result::LmmStatsResult, c_idx::Int)
    if !isnothing(result.p_corrected) && c_idx in result.tested && haskey(result.inference, :z_perm_obs)
        return @view result.inference.z_perm_obs[:, :, c_idx]
    end
    return @view result.z_values[:, :, c_idx]
end

# Permutation p-value of `obs` given the sorted null. With exhaustive enumeration the identity pattern is part
# of the null set (p = #{null ≥ obs} / N); with random flips the observed statistic is added (Phipson & Smyth, 2010).
function _lmm_perm_pvalue(null_sorted::AbstractVector, obs::Real, exhaustive::Bool)
    n = length(null_sorted)
    # Tolerance for floating-point ties, ensuring exact matches (such as the identity flip
    # in exhaustive enumeration) are included in the null count even with slight roundoff.
    tol = max(abs(obs) * 1e-8, 1e-12)
    threshold = obs - tol
    idx = searchsortedfirst(null_sorted, threshold)
    n_ge = n - idx + 1
    return exhaustive ? n_ge / n : (n_ge + 1) / (n + 1)
end

function _update_clusters(clusters, stats, null_sorted, alpha, exhaustive, dims)
    updated = Cluster[]
    mask = falses(dims)
    for (i, c) in enumerate(clusters)
        p_val = _lmm_perm_pvalue(null_sorted, abs(stats[i]), exhaustive)
        sig = p_val <= alpha
        push!(updated, Cluster(c.id, c.electrodes, c.time_indices, c.time_range, stats[i], p_val, sig, c.polarity, c.members))
        if sig
            for (e_idx, t_idx) in c.members
                mask[e_idx, t_idx] = true
            end
        end
    end
    return updated, mask
end
