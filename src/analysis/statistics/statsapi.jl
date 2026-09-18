# === StatsAPI Integration ===

"""
    StatsAPI.coef(res::LmmStatsResult)

Returns the coefficients (betas) for the fitted linear mixed models as a 3D array 
of dimensions [channels × time × coefficients].
"""
StatsAPI.coef(res::LmmStatsResult) = res.beta

"""
    StatsAPI.stderror(res::LmmStatsResult)

Returns the standard errors for the fitted linear mixed models as a 3D array 
of dimensions [channels × time × coefficients].
"""
StatsAPI.stderror(res::LmmStatsResult) = res.se

"""
    StatsAPI.nobs(res::LmmStatsResult)

Returns the number of observations (epochs) used in the mass-univariate fit.
"""
StatsAPI.nobs(res::LmmStatsResult) = length(res.epochs.data)

"""
    StatsAPI.pvalue(res::LmmStatsResult; corrected::Bool=false)

Returns the p-values for the fitted linear mixed models as a 3D array 
of dimensions [channels × time × coefficients].

If `corrected=true`, computes permutation-corrected max-t p-values by 
comparing the observed t-values against the permuted max-t null distribution.
Uses the formula `(count + 1) / (n_perms + 1)` (Phipson & Smyth, 2010).
"""
function StatsAPI.pvalue(res::LmmStatsResult; corrected::Bool=false)
    if !corrected
        return res.p_values
    else
        if isempty(res.max_t_null)
            error("Cannot compute permutation-corrected p-values: max_t_null is empty. Run fit_mass_lmm with n_perms > 0 and use_clusters=false.")
        end
        n_perms, n_coefs = size(res.max_t_null)
        p_corrected = similar(res.p_values)
        
        for c in 1:n_coefs
            # Sort the null distribution once for O(log n) lookups
            null_sorted = sort(res.max_t_null[:, c])
            for t in 1:size(res.t_values, 2)
                for ch in 1:size(res.t_values, 1)
                    obs_t = abs(res.t_values[ch, t, c])
                    # Number of null values >= obs_t, using sorted search
                    n_exceeding = n_perms - searchsortedlast(null_sorted, obs_t - eps(obs_t))
                    # Phipson & Smyth (2010) corrected formula
                    p_corrected[ch, t, c] = (n_exceeding + 1) / (n_perms + 1)
                end
            end
        end
        return p_corrected
    end
end

"""
    StatsAPI.confint(res::LmmStatsResult; level::Real=0.95)

Returns Wald-type confidence intervals for the fitted linear mixed models as a tuple 
of (lower, upper) bound 3D arrays of dimensions [channels × time × coefficients].

Uses the asymptotic normal approximation (consistent with MixedModels.jl's approach).
"""
function StatsAPI.confint(res::LmmStatsResult; level::Real=0.95)
    z = Distributions.quantile(Distributions.Normal(), 1 - (1-level)/2)
    se = res.se
    lower = res.beta .- z .* se
    upper = res.beta .+ z .* se
    return (lower, upper)
end

# === Support for Original Permutation/Analytic Tests ===

"""
    StatsAPI.coef(res::StatsResult)
    StatsAPI.coef(res::TFStatsResult)

Returns the mean difference (condition 1 - condition 2) for the statistical result.
"""
function StatsAPI.coef(res::Union{StatsResult, TFStatsResult})
    if length(res.data) == 2
        return res.data[1].data .- res.data[2].data
    elseif length(res.data) == 1
        return res.data[1].data
    else
        error("Unsupported number of conditions in data for coef.")
    end
end

"""
    StatsAPI.stderror(res::StatsResult)

Returns the standard error of the mean difference for standard ERP stats.
"""
StatsAPI.stderror(res::StatsResult) = hasfield(typeof(res), :se_diff) ? res.se_diff : error("Standard error not stored for this test type.")

"""
    StatsAPI.pvalue(res::Union{StatsResult, TFStatsResult})

Returns the p-values for the statistical test if available.
Note: For cluster-based permutation tests, point-wise p-values are mathematically undefined, 
so this will return `nothing`. You should inspect `res.clusters` or `res.masks` directly.
"""
function StatsAPI.pvalue(res::Union{StatsResult, TFStatsResult})
    if isnothing(res.stat_matrix.p)
        @warn "Cluster-based permutation tests do not produce true point-wise p-values (FWER is controlled at the cluster level). Use `res.masks` or `res.clusters` to determine significance."
    end
    return res.stat_matrix.p
end
