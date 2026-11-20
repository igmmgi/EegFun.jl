# StatsAPI.jl interface implementations for EegFun statistical result types.

######################################
########### User Functions ###########
######################################

"""
    StatsAPI.coef(res::LmmStatsResult)

Returns the fixed-effect coefficients (`beta`) for the fitted linear mixed models as a 3D array
of dimensions `[channels × time × coefficients]`.
"""
StatsAPI.coef(res::LmmStatsResult) = res.beta

"""
    StatsAPI.stderror(res::LmmStatsResult)

Returns the standard errors (`se`) for the fitted linear mixed models as a 3D array
of dimensions `[channels × time × coefficients]`.
"""
StatsAPI.stderror(res::LmmStatsResult) = res.se

"""
    StatsAPI.nobs(res::LmmStatsResult)

Returns the number of observations (epochs) used in the mass-univariate fit.
"""
StatsAPI.nobs(res::LmmStatsResult) = res.inference.n_obs

"""
    StatsAPI.pvalue(res::LmmStatsResult; corrected::Bool = !isnothing(res.p_corrected))

Returns the p-values for the fitted linear mixed models as a 3D array
of dimensions `[channels × time × coefficients]`.

If permutations were performed, `corrected` defaults to `true` and returns the FWER-corrected
permutation p-values (method in `res.inference.correction`; `NaN` for untested coefficients).
Otherwise, or with `corrected = false`, returns the uncorrected parametric p-values (`res.p_parametric`,
evaluated against a t-distribution with Between-Within degrees of freedom).
"""
function StatsAPI.pvalue(res::LmmStatsResult; corrected::Bool = !isnothing(res.p_corrected))
    if corrected
        isnothing(res.p_corrected) &&
            @minimal_error "No permutation-corrected p-values: run fit_mass_lmm with n_permutations > 0, or use corrected = false."
        return res.p_corrected
    end
    return res.p_parametric
end

"""
    StatsAPI.confint(res::LmmStatsResult; level::Real = 0.95)

Returns parametric confidence intervals for the fitted linear mixed models as a tuple
of `(lower, upper)` bound 3D arrays of dimensions `[channels × time × coefficients]`.

Uses the per-coefficient Between-Within t-distribution degrees of freedom stored in
`res.inference.df_parametric` when finite, or the asymptotic Normal quantile otherwise.
"""
function StatsAPI.confint(res::LmmStatsResult; level::Real = 0.95)
    0 < level < 1 || @minimal_error "Confidence level must be in (0, 1), got $level."
    tail_prob = 1 - (1 - Float64(level)) / 2
    n_ch, n_tp, n_coef = size(res.beta)
    lower = similar(res.beta)
    upper = similar(res.beta)
    dfs = haskey(res.inference, :df_parametric) ? res.inference.df_parametric : nothing
    for k = 1:n_coef
        df_k = (!isnothing(dfs) && k <= length(dfs)) ? dfs[k] : Inf
        crit =
            (isfinite(df_k) && df_k > 0) ? Distributions.quantile(Distributions.TDist(df_k), tail_prob) :
            Distributions.quantile(Distributions.Normal(), tail_prob)
        @views begin
            lower[:, :, k] .= res.beta[:, :, k] .- crit .* res.se[:, :, k]
            upper[:, :, k] .= res.beta[:, :, k] .+ crit .* res.se[:, :, k]
        end
    end
    return (lower, upper)
end

# === Support for ERP and Time-Frequency Permutation/Analytic Tests ===
"""
    StatsAPI.coef(res::StatsResult)
    StatsAPI.coef(res::TFStatsResult)

Returns the mean effect (`condition 1 - condition 2`, or `condition 1` for single-condition tests)
aligned with `res.stat_matrix.t` (`[electrodes × time]` for ERP stats, `[electrodes × frequencies × time]`
for time-frequency stats).
"""
function StatsAPI.coef(res::StatsResult)
    n_cond = length(res.data)
    (n_cond == 1 || n_cond == 2) || @minimal_error("Unsupported number of conditions ($n_cond) in StatsResult.")
    df1 = res.data[1].data
    t_col1 = df1[!, :time]
    t_idx1 = Int[argmin(abs.(t_col1 .- t)) for t in res.time_points]
    out = Matrix{Float64}(undef, length(res.electrodes), length(res.time_points))
    if n_cond == 2
        df2 = res.data[2].data
        t_col2 = df2[!, :time]
        t_idx2 = Int[argmin(abs.(t_col2 .- t)) for t in res.time_points]
        for (c_i, ch) in enumerate(res.electrodes)
            col1 = df1[!, ch]
            col2 = df2[!, ch]
            @inbounds for (k, (i1, i2)) in enumerate(zip(t_idx1, t_idx2))
                out[c_i, k] = Float64(col1[i1]) - Float64(col2[i2])
            end
        end
    else
        for (c_i, ch) in enumerate(res.electrodes)
            col1 = df1[!, ch]
            @inbounds for (k, i1) in enumerate(t_idx1)
                out[c_i, k] = Float64(col1[i1])
            end
        end
    end
    return out
end

function StatsAPI.coef(res::TFStatsResult)
    n_cond = length(res.data)
    (n_cond == 1 || n_cond == 2) || @minimal_error("Unsupported number of conditions ($n_cond) in TFStatsResult.")
    n_elec = length(res.electrodes)
    n_freq = length(res.frequencies)
    n_time = length(res.time_points)
    out = fill(NaN, n_elec, n_freq, n_time)

    df1 = res.data[1].data_power
    f_map = Dict{Float64,Int}(f => i for (i, f) in enumerate(res.frequencies))
    t_map = Dict{Float64,Int}(t => i for (i, t) in enumerate(res.time_points))

    for row = 1:nrow(df1)
        fi = get(f_map, Float64(df1[row, :freq]), 0)
        ti = get(t_map, Float64(df1[row, :time]), 0)
        (fi == 0 || ti == 0) && continue
        for (ci, ch) in enumerate(res.electrodes)
            out[ci, fi, ti] = Float64(df1[row, ch])
        end
    end
    if n_cond == 2
        df2 = res.data[2].data_power
        for row = 1:nrow(df2)
            fi = get(f_map, Float64(df2[row, :freq]), 0)
            ti = get(t_map, Float64(df2[row, :time]), 0)
            (fi == 0 || ti == 0) && continue
            for (ci, ch) in enumerate(res.electrodes)
                out[ci, fi, ti] -= Float64(df2[row, ch])
            end
        end
    end
    return out
end

"""
    StatsAPI.stderror(res::StatsResult)

Returns the standard error of the mean difference over the analysis interval for ERP stats.
"""
StatsAPI.stderror(res::StatsResult) =
    hasfield(typeof(res), :se) ? res.se :
    hasfield(typeof(res), :se_diff) ? res.se_diff : @minimal_error("Standard error not stored for this test type.")

"""
    StatsAPI.pvalue(res::StatsResult)

Returns the p-values for the statistical test if available.
Note: For cluster-based permutation tests, point-wise p-values are mathematically undefined,
so this will return `nothing`. You should inspect `res.clusters` or `res.masks` directly.
"""
function StatsAPI.pvalue(res::StatsResult)
    if isnothing(res.stat_matrix.p)
        @minimal_warning "Cluster-based permutation tests do not produce true point-wise p-values (FWER is controlled at the cluster level). Use `res.masks` or `res.clusters` to determine significance."
    end
    return res.stat_matrix.p
end
