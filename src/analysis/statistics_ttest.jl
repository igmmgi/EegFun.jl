"""
Basic t-test functions for EEG/ERP data analysis.

This module provides fundamental t-test functions for computing t-statistics
and p-values for paired and independent sample designs.
"""

# ==============
# PAIRED T-TEST
# ==============
"""
    paired_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

Compute paired t-test with degrees of freedom, t-statistic, and p-value.

# Arguments
- `x::AbstractVector`: First condition data (must have same length as y)
- `y::AbstractVector`: Second condition data (must have same length as x)
- `tail::Symbol`: Type of test - `:both` (two-tailed, default), `:left` (one-tailed, A < B), or `:right` (one-tailed, A > B)

# Returns
- `TTestResult`: Struct containing `df` (degrees of freedom), `t` (t-statistic), and `p` (p-value).
  Returns `NaN` for p-value if t-value is `NaN` or `Inf`.
"""
function paired_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

    # validate equal lengths
    length(x) == length(y) || error("Paired t-test requires equal sample sizes")
    
    n = length(x)
    n < 2 && return (df = NaN, t = NaN, p = NaN)  # Need at least 2 observations
    
    df = n - 1  # degrees of freedom for paired t-test
    
    # Compute mean and std of differences 
    diff = x .- y
    mean_diff = mean(diff)
    std_diff = std(diff, corrected = true)  # Sample standard deviation (n-1)
    
    # Compute t-value
    if std_diff == 0.0
        if mean_diff == 0.0
            t = NaN
        else
            t = Inf * sign(mean_diff)
        end
    else
        t = mean_diff / (std_diff / sqrt(n))
    end
    
    # Handle edge cases
    (isnan(t) || isinf(t)) && return (df = df, t = t, p = NaN)
    
    # Compute p-value
    dist = TDist(df)
    if tail == :both # Two-tailed test
        p = 2 * (1 - cdf(dist, abs(t)))
    elseif tail == :left # One-tailed: H0: A >= B, H1: A < B
        p = cdf(dist, t)
    elseif tail == :right # One-tailed: H0: A <= B, H1: A > B
        p = 1 - cdf(dist, t)
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return (df = df, t = t, p = p)
end

# ==================
# INDEPENDENT T-TEST
# ==================
"""
    independent_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

Compute independent sample t-test with degrees of freedom, t-statistic, and p-value.
Assumes equal variances (standard independent t-test).

# Arguments
- `x::AbstractVector`: First group data (must have at least 2 observations)
- `y::AbstractVector`: Second group data (must have at least 2 observations)
- `tail::Symbol`: Type of test - `:both` (two-tailed, default), `:left` (one-tailed, A < B), or `:right` (one-tailed, A > B)

# Returns
- `TTestResult`: Struct containing `df` (degrees of freedom), `t` (t-statistic), and `p` (p-value).
  Returns `NaN` for p-value if t-value is `NaN` or `Inf`.
"""
function independent_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

    # Validate input lengths
    n_A, n_B = length(x), length(y)
    (n_A < 2 || n_B < 2) && error("Independent t-test requires at least 2 observations per group")
    
    df = n_A + n_B - 2  # degrees of freedom for independent t-test
    
    # Compute means and variances (Statistics.jl is already optimized with SIMD)
    mean_x, mean_y = mean(x), mean(y)
    
    # Pooled variance (assuming equal variances)
    var_x = var(x, corrected = true)
    var_y = var(y, corrected = true)
    pooled_var = ((n_A - 1) * var_x + (n_B - 1) * var_y) / df
    
    # Compute t-value
    if pooled_var == 0.0
        if mean_x == mean_y
            t = NaN  # Both groups identical
        else
            t = Inf * sign(mean_x - mean_y)
        end
    else
        t = (mean_x - mean_y) / sqrt(pooled_var * (1/n_A + 1/n_B))
    end
    
    # Handle edge cases
    (isnan(t) || isinf(t)) && return (df = df, t = t, p = NaN)
    
    # Compute p-value
    dist = TDist(df)
    if tail == :both # Two-tailed test
        p = 2 * (1 - cdf(dist, abs(t)))
    elseif tail == :left # One-tailed: H0: A >= B, H1: A < B
        p = cdf(dist, t)
    elseif tail == :right # One-tailed: H0: A <= B, H1: A > B
        p = 1 - cdf(dist, t)
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return (df = df, t = t, p = p)
end

