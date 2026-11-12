# Parametric (t-distribution) and non-parametric (permutation-quantile) thresholding.

function _threshold_t_matrix_parametric!(
    mask_positive::BitArray{N},
    mask_negative::BitArray{N},
    t_matrix::AbstractArray{Float64,N},
    critical_t_values::AbstractArray{Float64,N},
    tail::Symbol = :both,
) where {N}
    if tail == :both
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            if isfinite(t_val) && isfinite(crit_t)
                mask_positive[i] = t_val > crit_t
                mask_negative[i] = t_val < -crit_t
            else
                mask_positive[i] = false
                mask_negative[i] = false
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            mask_positive[i] = isfinite(t_val) && isfinite(crit_t) && t_val > crit_t
        end
    elseif tail == :left
        fill!(mask_positive, false)
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            mask_negative[i] = isfinite(t_val) && isfinite(crit_t) && t_val < crit_t
        end
    else
        @minimal_error("tail must be :both, :left, or :right, got :$tail")
    end
    return nothing
end

function _threshold_t_matrix_parametric(
    t_matrix::AbstractArray{Float64,N},
    critical_t_values::AbstractArray{Float64,N},
    tail::Symbol = :both,
) where {N}
    mask_positive = BitArray{N}(undef, size(t_matrix))
    mask_negative = BitArray{N}(undef, size(t_matrix))
    _threshold_t_matrix_parametric!(mask_positive, mask_negative, t_matrix, critical_t_values, tail)
    return mask_positive, mask_negative
end

function _compute_nonparametric_threshold_common(
    permutation_t_matrices::AbstractArray{Float64,N},
    alpha::Real = 0.05,
    tail::Symbol = :both,
) where {N}
    tail ∈ (:both, :left, :right) || @minimal_error("tail must be :both, :left, or :right, got :$tail")

    buf = Vector{Float64}(undef, length(permutation_t_matrices))
    cnt = 0
    if tail == :both
        @inbounds for i in eachindex(permutation_t_matrices)
            t_val = permutation_t_matrices[i]
            if isfinite(t_val)
                cnt += 1
                buf[cnt] = abs(t_val)
            end
        end
        cnt == 0 && @minimal_error("No valid t-values found in permutation distribution")
        thresh = quantile!(view(buf, 1:cnt), 1.0 - alpha)
        return thresh, thresh
    else
        @inbounds for i in eachindex(permutation_t_matrices)
            t_val = permutation_t_matrices[i]
            if isfinite(t_val)
                cnt += 1
                buf[cnt] = t_val
            end
        end
        cnt == 0 && @minimal_error("No valid t-values found in permutation distribution")
        if tail == :right
            return quantile!(view(buf, 1:cnt), 1.0 - alpha), NaN
        else
            return NaN, -quantile!(view(buf, 1:cnt), alpha)
        end
    end
end

function _compute_nonparametric_threshold_individual(
    permutation_t_matrices::Array{Float64,N},
    alpha::Real = 0.05,
    tail::Symbol = :both,
) where {N}
    tail ∈ (:both, :left, :right) || @minimal_error("tail must be :both, :left, or :right, got :$tail")

    grid_dims = ntuple(d -> size(permutation_t_matrices, d), Val(N - 1))
    n_voxels = prod(grid_dims)
    n_permutations = size(permutation_t_matrices, N)

    thresholds_positive = Array{Float64,N - 1}(undef, grid_dims)
    thresholds_negative = Array{Float64,N - 1}(undef, grid_dims)
    tail == :right && fill!(thresholds_negative, NaN)
    tail == :left && fill!(thresholds_positive, NaN)

    perm_flat = reshape(permutation_t_matrices, n_voxels, n_permutations)
    alpha_f = Float64(alpha)
    percentile_level = 1.0 - alpha_f

    _foreach_chunk(n_voxels) do v_range
        scratch = Vector{Float64}(undef, n_permutations)
        if tail == :both
            @inbounds for v in v_range
                cnt = 0
                for perm_idx = 1:n_permutations
                    t_val = perm_flat[v, perm_idx]
                    if isfinite(t_val)
                        cnt += 1
                        scratch[cnt] = abs(t_val)
                    end
                end
                thresh = cnt == 0 ? NaN : quantile!(view(scratch, 1:cnt), percentile_level)
                thresholds_positive[v] = thresh
                thresholds_negative[v] = thresh
            end
        elseif tail == :right
            @inbounds for v in v_range
                cnt = 0
                for perm_idx = 1:n_permutations
                    t_val = perm_flat[v, perm_idx]
                    if isfinite(t_val)
                        cnt += 1
                        scratch[cnt] = t_val
                    end
                end
                thresholds_positive[v] = cnt == 0 ? NaN : quantile!(view(scratch, 1:cnt), percentile_level)
            end
        else
            @inbounds for v in v_range
                cnt = 0
                for perm_idx = 1:n_permutations
                    t_val = perm_flat[v, perm_idx]
                    if isfinite(t_val)
                        cnt += 1
                        scratch[cnt] = t_val
                    end
                end
                thresholds_negative[v] = cnt == 0 ? NaN : -quantile!(view(scratch, 1:cnt), alpha_f)
            end
        end
    end

    return thresholds_positive, thresholds_negative
end

function _threshold_t_matrix_nonparametric!(
    mask_positive::BitArray{N},
    mask_negative::BitArray{N},
    t_matrix::AbstractArray{Float64,N},
    thresh_pos::Float64,
    thresh_neg::Float64,
    tail::Symbol = :both,
) where {N}
    if tail == :both
        if isnan(thresh_pos) || isnan(thresh_neg)
            fill!(mask_positive, false)
            fill!(mask_negative, false)
            return nothing
        end
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            if isfinite(t_val)
                mask_positive[i] = t_val > thresh_pos
                mask_negative[i] = t_val < -thresh_neg
            else
                mask_positive[i] = false
                mask_negative[i] = false
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        if isnan(thresh_pos)
            fill!(mask_positive, false)
            return nothing
        end
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            mask_positive[i] = isfinite(t_val) && t_val > thresh_pos
        end
    elseif tail == :left
        fill!(mask_positive, false)
        if isnan(thresh_neg)
            fill!(mask_negative, false)
            return nothing
        end
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            mask_negative[i] = isfinite(t_val) && t_val < -thresh_neg
        end
    else
        @minimal_error("tail must be :both, :left, or :right, got :$tail")
    end
    return nothing
end

function _threshold_t_matrix_nonparametric!(
    mask_positive::BitArray{N},
    mask_negative::BitArray{N},
    t_matrix::AbstractArray{Float64,N},
    thresholds_positive::AbstractArray{Float64,N},
    thresholds_negative::AbstractArray{Float64,N},
    tail::Symbol = :both,
) where {N}
    if tail == :both
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            tp = thresholds_positive[i]
            tn = thresholds_negative[i]
            if isfinite(t_val) && !isnan(tp) && !isnan(tn)
                mask_positive[i] = t_val > tp
                mask_negative[i] = t_val < -tn
            else
                mask_positive[i] = false
                mask_negative[i] = false
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            tp = thresholds_positive[i]
            mask_positive[i] = isfinite(t_val) && !isnan(tp) && t_val > tp
        end
    elseif tail == :left
        fill!(mask_positive, false)
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            tn = thresholds_negative[i]
            mask_negative[i] = isfinite(t_val) && !isnan(tn) && t_val < -tn
        end
    else
        @minimal_error("tail must be :both, :left, or :right, got :$tail")
    end
    return nothing
end

function _threshold_t_matrix_nonparametric(
    t_matrix::AbstractArray{Float64,N},
    thresholds_positive::Union{Float64,AbstractArray{Float64,N}},
    thresholds_negative::Union{Float64,AbstractArray{Float64,N}},
    tail::Symbol = :both,
) where {N}
    mask_positive = BitArray{N}(undef, size(t_matrix))
    mask_negative = BitArray{N}(undef, size(t_matrix))
    _threshold_t_matrix_nonparametric!(mask_positive, mask_negative, t_matrix, thresholds_positive, thresholds_negative, tail)
    return mask_positive, mask_negative
end

const _threshold_t_matrix_parametric_tf! = _threshold_t_matrix_parametric!
const _threshold_t_matrix_parametric_tf = _threshold_t_matrix_parametric
const _compute_nonparametric_threshold_common_tf = _compute_nonparametric_threshold_common
const _compute_nonparametric_threshold_individual_tf = _compute_nonparametric_threshold_individual
const _threshold_t_matrix_nonparametric_tf! = _threshold_t_matrix_nonparametric!
