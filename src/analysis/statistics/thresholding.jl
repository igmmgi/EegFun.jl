# This file contains all thresholding logic for statistical testing,
# including both parametric (t-distribution) and non-parametric (permutation-based)
# thresholding methods.
"""
    threshold_t_matrix_parametric!(mask_positive, mask_negative, t_matrix, critical_t_values, tail)

In-place version: Threshold t-matrix using parametric critical values.

Fills pre-allocated mask arrays with significance results.

# Arguments
- `mask_positive::BitArray{2}`: Output mask for positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Output mask for negative significant points [electrodes × time]
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time]
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right`
"""
function _threshold_t_matrix_parametric!(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    t_matrix::Array{Float64,2},
    critical_t_values::Array{Float64,2},
    tail::Symbol = :both,
)
    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]

            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = t_val > crit_t
                mask_negative[i] = t_val < -crit_t
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]

            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
            else
                mask_positive[i] = t_val > crit_t
            end
        end
    elseif tail == :left
        fill!(mask_positive, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]

            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_negative[i] = false
            else
                mask_negative[i] = t_val < crit_t
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

"""
    threshold_t_matrix_parametric(t_matrix, critical_t_values, tail)

Threshold t-matrix using parametric critical values.

Creates new mask arrays and calls the in-place version.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time]
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `mask_positive::BitArray{2}`: Mask for positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Mask for negative significant points [electrodes × time]

# Examples
```julia
mask_pos, mask_neg = threshold_t_matrix_parametric(t_matrix, critical_t, :both)
```
"""
function _threshold_t_matrix_parametric(t_matrix::Array{Float64,2}, critical_t_values::Array{Float64,2}, tail::Symbol = :both)
    n_electrodes, n_time = size(t_matrix)
    mask_positive = BitArray{2}(undef, n_electrodes, n_time)
    mask_negative = BitArray{2}(undef, n_electrodes, n_time)

    _threshold_t_matrix_parametric!(mask_positive, mask_negative, t_matrix, critical_t_values, tail)

    return mask_positive, mask_negative
end

# ===================
# NON-PARAMETRIC THRESHOLDING
# ===================

"""
    _collect_valid_t_values(permutation_t_matrices, predicate, transform)

Helper function to collect valid t-values from permutation distribution.

Reduces code duplication in nonparametric threshold computation.

# Arguments
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations
- `predicate::Function`: Function to filter t-values (e.g., t -> t > 0)
- `transform::Function`: Function to transform t-values (e.g., abs)

# Returns
- `all_t_values::Vector{Float64}`: Collected valid t-values
"""
function _collect_valid_t_values(permutation_t_matrices::Array{Float64,3}, predicate::Function, transform::Function)
    n_electrodes, n_time, n_permutations = size(permutation_t_matrices)
    all_t_values = Float64[]
    sizehint!(all_t_values, n_electrodes * n_time * n_permutations)

    for perm_idx = 1:n_permutations
        for i = 1:n_electrodes
            for j = 1:n_time
                t_val = permutation_t_matrices[i, j, perm_idx]
                if !isnan(t_val) && !isinf(t_val) && predicate(t_val)
                    push!(all_t_values, transform(t_val))
                end
            end
        end
    end

    return all_t_values
end

"""
    compute_nonparametric_threshold_common(permutation_t_matrices, alpha, tail)

Compute a single non-parametric threshold from pooled permutation distribution (common threshold).

# Arguments
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]
- `alpha::Float64`: Significance level (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `threshold_positive::Float64`: Threshold for positive t-values
- `threshold_negative::Float64`: Threshold for negative t-values (absolute value)

# Examples
```julia
thresh_pos, thresh_neg = compute_nonparametric_threshold_common(perm_t_matrices, 0.05, :both)
```
"""
function _compute_nonparametric_threshold_common(permutation_t_matrices::Array{Float64,3}, alpha::Float64 = 0.05, tail::Symbol = :both)
    if tail == :both
        # Two-tailed: collect all absolute t-values
        all_t_values = _collect_valid_t_values(permutation_t_matrices, t -> true, abs)

        if isempty(all_t_values)
            error("No valid t-values found in permutation distribution")
        end

        # Compute (1 - alpha) percentile for two-tailed on absolute values
        # Since |t| is always positive, the (1-alpha) percentile of |t|
        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)

        return threshold, threshold

    elseif tail == :right
        # One-tailed right: collect all positive t-values
        all_t_values = _collect_valid_t_values(permutation_t_matrices, t -> t > 0, identity)

        if isempty(all_t_values)
            error("No valid positive t-values found in permutation distribution")
        end

        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)

        return threshold, NaN

    elseif tail == :left
        # One-tailed left: collect all negative t-values (absolute values)
        all_t_values = _collect_valid_t_values(permutation_t_matrices, t -> t < 0, abs)

        if isempty(all_t_values)
            error("No valid negative t-values found in permutation distribution")
        end

        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)

        return NaN, threshold

    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

"""
    compute_nonparametric_threshold_individual(permutation_t_matrices, alpha, tail)

Compute point-specific non-parametric thresholds from permutation distribution (individual thresholds).

# Arguments
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]
- `alpha::Float64`: Significance level (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `thresholds_positive::Array{Float64, 2}`: Point-specific thresholds for positive t-values [electrodes × time]
- `thresholds_negative::Array{Float64, 2}`: Point-specific thresholds for negative t-values [electrodes × time]

# Examples
```julia
thresh_pos, thresh_neg = compute_nonparametric_threshold_individual(perm_t_matrices, 0.05, :both)
```
"""
function _compute_nonparametric_threshold_individual(permutation_t_matrices::Array{Float64,3}, alpha::Float64 = 0.05, tail::Symbol = :both)
    n_electrodes, n_time, n_permutations = size(permutation_t_matrices)

    thresholds_positive = Array{Float64,2}(undef, n_electrodes, n_time)
    thresholds_negative = Array{Float64,2}(undef, n_electrodes, n_time)

    if tail == :both
        # Two-tailed: for each point, compute (1 - alpha) percentile of |t|
        # Since |t| is always positive, 1-alpha of |t| = 1-alpha/2 of raw t
        percentile_level = 1.0 - alpha

        for i = 1:n_electrodes
            for j = 1:n_time
                # Collect t-values at this point across all permutations
                t_values = Float64[]
                sizehint!(t_values, n_permutations)

                for perm_idx = 1:n_permutations
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val)
                        push!(t_values, abs(t_val))
                    end
                end

                if isempty(t_values)
                    thresholds_positive[i, j] = NaN
                    thresholds_negative[i, j] = NaN
                else
                    threshold = quantile(t_values, percentile_level)
                    thresholds_positive[i, j] = threshold
                    thresholds_negative[i, j] = threshold
                end
            end
        end

    elseif tail == :right
        # One-tailed right: for each point, compute (1 - alpha) percentile of positive values
        percentile_level = 1.0 - alpha

        for i = 1:n_electrodes
            for j = 1:n_time
                # Collect positive t-values at this point
                t_values = Float64[]
                sizehint!(t_values, n_permutations)

                for perm_idx = 1:n_permutations
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val > 0
                        push!(t_values, t_val)
                    end
                end

                if isempty(t_values)
                    thresholds_positive[i, j] = NaN
                else
                    thresholds_positive[i, j] = quantile(t_values, percentile_level)
                end
                thresholds_negative[i, j] = NaN
            end
        end

    elseif tail == :left
        # One-tailed left: for each point, compute (1 - alpha) percentile of negative values (absolute)
        percentile_level = 1.0 - alpha

        for i = 1:n_electrodes
            for j = 1:n_time
                # Collect negative t-values at this point (absolute values)
                t_values = Float64[]
                sizehint!(t_values, n_permutations)

                for perm_idx = 1:n_permutations
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val < 0
                        push!(t_values, abs(t_val))
                    end
                end

                thresholds_positive[i, j] = NaN
                if isempty(t_values)
                    thresholds_negative[i, j] = NaN
                else
                    thresholds_negative[i, j] = quantile(t_values, percentile_level)
                end
            end
        end

    else
        error("tail must be :both, :left, or :right, got :$tail")
    end

    return thresholds_positive, thresholds_negative
end

"""
    threshold_t_matrix_nonparametric!(mask_positive, mask_negative, t_matrix, thresholds_positive, thresholds_negative, tail)

In-place version: Threshold t-matrix using non-parametric thresholds.

Fills pre-allocated mask arrays with significance results.

# Arguments
- `mask_positive::BitArray{2}`: Output mask for positive significant points
- `mask_negative::BitArray{2}`: Output mask for negative significant points
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `thresholds_positive::Union{Float64, Array{Float64, 2}}`: Threshold(s) for positive t-values
- `thresholds_negative::Union{Float64, Array{Float64, 2}}`: Threshold(s) for negative t-values
- `tail::Symbol`: Test tail - `:both`, `:left`, or `:right`
"""
function _threshold_t_matrix_nonparametric!(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    t_matrix::Array{Float64,2},
    thresholds_positive::Union{Float64,Array{Float64,2}},
    thresholds_negative::Union{Float64,Array{Float64,2}},
    tail::Symbol = :both,
)
    # Determine if thresholds are scalar (common) or matrix (individual)
    is_common = isa(thresholds_positive, Float64)

    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]

            if isnan(t_val) || isinf(t_val)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                # Get threshold for this point
                thresh_pos = is_common ? thresholds_positive::Float64 : thresholds_positive[i]
                thresh_neg = is_common ? thresholds_negative::Float64 : thresholds_negative[i]

                if isnan(thresh_pos) || isnan(thresh_neg)
                    mask_positive[i] = false
                    mask_negative[i] = false
                else
                    mask_positive[i] = t_val > thresh_pos
                    mask_negative[i] = t_val < -thresh_neg
                end
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]

            if isnan(t_val) || isinf(t_val)
                mask_positive[i] = false
            else
                thresh_pos = is_common ? thresholds_positive::Float64 : thresholds_positive[i]
                if isnan(thresh_pos)
                    mask_positive[i] = false
                else
                    mask_positive[i] = t_val > thresh_pos
                end
            end
        end
    elseif tail == :left
        fill!(mask_positive, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]

            if isnan(t_val) || isinf(t_val)
                mask_negative[i] = false
            else
                thresh_neg = is_common ? thresholds_negative::Float64 : thresholds_negative[i]
                if isnan(thresh_neg)
                    mask_negative[i] = false
                else
                    mask_negative[i] = t_val < -thresh_neg
                end
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

"""
    threshold_t_matrix_nonparametric(t_matrix, thresholds_positive, thresholds_negative, tail)

Threshold t-matrix using non-parametric thresholds.

Creates new mask arrays and calls the in-place version.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `thresholds_positive::Union{Float64, Array{Float64, 2}}`: Threshold(s) for positive t-values (scalar for common, matrix for individual)
- `thresholds_negative::Union{Float64, Array{Float64, 2}}`: Threshold(s) for negative t-values (scalar for common, matrix for individual)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `mask_positive::BitArray{2}`: Mask for positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Mask for negative significant points [electrodes × time]

# Examples
```julia
# Common threshold (scalar)
mask_pos, mask_neg = threshold_t_matrix_nonparametric(t_matrix, thresh_pos, thresh_neg, :both)

# Individual thresholds (matrices)
mask_pos, mask_neg = threshold_t_matrix_nonparametric(t_matrix, thresh_pos_mat, thresh_neg_mat, :both)
```
"""
function _threshold_t_matrix_nonparametric(
    t_matrix::Array{Float64,2},
    thresholds_positive::Union{Float64,Array{Float64,2}},
    thresholds_negative::Union{Float64,Array{Float64,2}},
    tail::Symbol = :both,
)
    n_electrodes, n_time = size(t_matrix)
    mask_positive = BitArray{2}(undef, n_electrodes, n_time)
    mask_negative = BitArray{2}(undef, n_electrodes, n_time)

    _threshold_t_matrix_nonparametric!(mask_positive, mask_negative, t_matrix, thresholds_positive, thresholds_negative, tail)

    return mask_positive, mask_negative
end


# ===================
# TF PARAMETRIC THRESHOLDING (3D: electrodes × frequencies × time)
# ===================

"""
    _threshold_t_matrix_parametric_tf!(mask_positive, mask_negative, t_matrix, critical_t_values, tail)

In-place parametric thresholding of 3D t-matrix for TF data.
"""
function _threshold_t_matrix_parametric_tf!(
    mask_positive::BitArray{3},
    mask_negative::BitArray{3},
    t_matrix::Array{Float64,3},
    critical_t_values::Array{Float64,3},
    tail::Symbol = :both,
)
    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]

            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = t_val > crit_t
                mask_negative[i] = t_val < -crit_t
            end
        end
    elseif tail == :right
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]

            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = t_val > crit_t
                mask_negative[i] = false
            end
        end
    elseif tail == :left
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]

            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = false
                mask_negative[i] = t_val < crit_t
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

"""
    _threshold_t_matrix_parametric_tf(t_matrix, critical_t_values, tail)

Parametric thresholding of 3D t-matrix (allocating version).
"""
function _threshold_t_matrix_parametric_tf(t_matrix::Array{Float64,3}, critical_t_values::Array{Float64,3}, tail::Symbol = :both)
    mask_positive = falses(size(t_matrix))
    mask_negative = falses(size(t_matrix))
    _threshold_t_matrix_parametric_tf!(mask_positive, mask_negative, t_matrix, critical_t_values, tail)
    return mask_positive, mask_negative
end


# ===================
# TF NON-PARAMETRIC THRESHOLDING (3D)
# ===================

"""
    _collect_valid_t_values_tf(permutation_t_matrices, predicate, transform)

Collect valid t-values from 4D TF permutation distribution.
"""
function _collect_valid_t_values_tf(permutation_t_matrices::Array{Float64,4}, predicate::Function, transform::Function)
    all_t_values = Float64[]
    sizehint!(all_t_values, length(permutation_t_matrices) ÷ 2)
    @inbounds for i in eachindex(permutation_t_matrices)
        t_val = permutation_t_matrices[i]
        if !isnan(t_val) && !isinf(t_val) && predicate(t_val)
            push!(all_t_values, transform(t_val))
        end
    end
    return all_t_values
end


"""
    _compute_nonparametric_threshold_common_tf(permutation_t_matrices, alpha, tail)

Compute common non-parametric thresholds from 4D TF permutation distribution.
"""
function _compute_nonparametric_threshold_common_tf(permutation_t_matrices::Array{Float64,4}, alpha::Float64 = 0.05, tail::Symbol = :both)
    if tail == :both
        all_t_values = _collect_valid_t_values_tf(permutation_t_matrices, t -> true, abs)
        isempty(all_t_values) && return NaN, NaN
        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)
        return threshold, threshold
    elseif tail == :right
        all_t_values = _collect_valid_t_values_tf(permutation_t_matrices, t -> t > 0, identity)
        isempty(all_t_values) && return NaN, NaN
        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)
        return threshold, NaN
    elseif tail == :left
        all_t_values = _collect_valid_t_values_tf(permutation_t_matrices, t -> t < 0, abs)
        isempty(all_t_values) && return NaN, NaN
        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)
        return NaN, threshold
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end


"""
    _threshold_t_matrix_nonparametric_tf!(mask_positive, mask_negative, t_matrix, thresh_pos, thresh_neg, tail)

In-place non-parametric thresholding of 3D t-matrix for TF data.
"""
function _threshold_t_matrix_nonparametric_tf!(
    mask_positive::BitArray{3},
    mask_negative::BitArray{3},
    t_matrix::Array{Float64,3},
    thresholds_positive::Union{Float64,Array{Float64,3}},
    thresholds_negative::Union{Float64,Array{Float64,3}},
    tail::Symbol = :both,
)
    is_common = isa(thresholds_positive, Float64)

    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            thresh_pos = is_common ? thresholds_positive : thresholds_positive[i]
            thresh_neg = is_common ? thresholds_negative : thresholds_negative[i]

            if isnan(t_val) || isinf(t_val) || isnan(thresh_pos) || isnan(thresh_neg)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = t_val > thresh_pos
                mask_negative[i] = t_val < -thresh_neg
            end
        end
    elseif tail == :right
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            thresh_pos = is_common ? thresholds_positive : thresholds_positive[i]

            if isnan(t_val) || isinf(t_val) || isnan(thresh_pos)
                mask_positive[i] = false
            else
                mask_positive[i] = t_val > thresh_pos
            end
            mask_negative[i] = false
        end
    elseif tail == :left
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            thresh_neg = is_common ? thresholds_negative : thresholds_negative[i]

            mask_positive[i] = false
            if isnan(t_val) || isinf(t_val) || isnan(thresh_neg)
                mask_negative[i] = false
            else
                mask_negative[i] = t_val < -thresh_neg
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end
