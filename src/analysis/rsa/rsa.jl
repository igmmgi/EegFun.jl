# Representational Similarity Analysis (RSA): RDM computation, normalization, model comparison, and grand averaging.

function _extract_upper_triangular!(vec::Vector{Float64}, rdm::AbstractMatrix{<:Real})
    n = size(rdm, 1)
    k = 1
    @inbounds for i = 1:n
        for j = (i+1):n
            vec[k] = Float64(rdm[i, j])
            k += 1
        end
    end
    return vec
end

function _extract_upper_triangular(rdm::AbstractMatrix{Float64})
    n = size(rdm, 1)
    return _extract_upper_triangular!(Vector{Float64}(undef, (n * (n - 1)) ÷ 2), rdm)
end

function _symmetrized_upper_triangular!(vec::AbstractVector{Float64}, rdm::AbstractMatrix{<:Real})
    n = size(rdm, 1)
    k = 1
    @inbounds for i = 1:n
        for j = (i+1):n
            vec[k] = 0.5 * (Float64(rdm[i, j]) + Float64(rdm[j, i]))
            k += 1
        end
    end
    return vec
end

function _symmetrized_upper_triangular(rdm::AbstractMatrix{<:Real})
    n = size(rdm, 1)
    return _symmetrized_upper_triangular!(Vector{Float64}(undef, (n * (n - 1)) ÷ 2), rdm)
end

function _tiedrank!(ranks::Vector{Float64}, perm::Vector{Int}, x::AbstractVector{Float64})
    n = length(x)
    any(isnan, x) && (fill!(ranks, NaN); return ranks)
    sortperm!(perm, x)
    i = 1
    @inbounds while i <= n
        v = x[perm[i]]
        j = i
        while j < n && x[perm[j+1]] == v
            j += 1
        end
        avg_rank = 0.5 * (i + j)
        for k = i:j
            ranks[perm[k]] = avg_rank
        end
        i = j + 1
    end
    return ranks
end

function _zscore_unit_norm!(dest::AbstractVector{Float64}, src::AbstractVector{Float64})
    n = length(src)
    if n < 2 || any(!isfinite, src)
        fill!(dest, 0.0)
        return false
    end
    μ = mean(src)
    norm_sq = 0.0
    @inbounds @simd for k = 1:n
        d = src[k] - μ
        dest[k] = d
        norm_sq += d * d
    end
    if norm_sq <= 1e-30
        fill!(dest, 0.0)
        return false
    end
    inv_norm = 1.0 / sqrt(norm_sq)
    @inbounds @simd for k = 1:n
        dest[k] *= inv_norm
    end
    return true
end

# Pre-allocated buffers for within-condition demeaned Ledoit-Wolf shrinkage covariance
struct _PooledCovWorkspace
    residuals::Matrix{Float64}
    cond_mean::Vector{Float64}
    S::Matrix{Float64}
end

function _PooledCovWorkspace(n_features::Int, max_total_trials::Int)
    return _PooledCovWorkspace(
        Matrix{Float64}(undef, n_features, max(max_total_trials, 1)),
        Vector{Float64}(undef, n_features),
        Matrix{Float64}(undef, n_features, n_features),
    )
end

function _accumulate_cond_residuals!(
    residuals::Matrix{Float64},
    cond_mean::Vector{Float64},
    cond_data::Array{Float64,3},
    time_idx::Int,
    idxs::AbstractVector{Int},
    col_offset::Int,
)
    n_tr = length(idxs)
    n_tr < 2 && return col_offset
    n_features = size(residuals, 1)
    fill!(cond_mean, 0.0)
    @inbounds for tr_idx in idxs
        @simd for f = 1:n_features
            cond_mean[f] += cond_data[f, time_idx, tr_idx]
        end
    end
    inv_n = 1.0 / n_tr
    @inbounds @simd for f = 1:n_features
        cond_mean[f] *= inv_n
    end
    @inbounds for tr_idx in idxs
        col_offset += 1
        @simd for f = 1:n_features
            residuals[f, col_offset] = cond_data[f, time_idx, tr_idx] - cond_mean[f]
        end
    end
    return col_offset
end

function _compute_pooled_covariance!(
    out_cov::Matrix{Float64},
    ws::_PooledCovWorkspace,
    data_arrays::Vector{Array{Float64,3}},
    trial_indices::Union{Nothing,Vector{Vector{Int}}},
    time_idx::Int,
)
    n_features = size(out_cov, 1)
    total_trials = 0
    n_dof = 0
    for (c, cond_data) in enumerate(data_arrays)
        n_tr = isnothing(trial_indices) ? size(cond_data, 3) : length(trial_indices[c])
        if n_tr >= 2
            total_trials += n_tr
            n_dof += n_tr - 1
        end
    end

    if n_dof < 1 || total_trials == 0
        fill!(out_cov, 0.0)
        @inbounds for i = 1:n_features
            out_cov[i, i] = 1.0
        end
        return out_cov
    end

    col_offset = 0
    for (c, cond_data) in enumerate(data_arrays)
        if isnothing(trial_indices)
            col_offset = _accumulate_cond_residuals!(ws.residuals, ws.cond_mean, cond_data, time_idx, 1:size(cond_data, 3), col_offset)
        else
            col_offset = _accumulate_cond_residuals!(ws.residuals, ws.cond_mean, cond_data, time_idx, trial_indices[c], col_offset)
        end
    end

    res_view = @view ws.residuals[:, 1:total_trials]
    S = ws.S
    mul!(S, res_view, transpose(res_view), 1.0 / total_trials, 0.0)

    mu = tr(S) / n_features
    normS_sq = 0.0
    d2 = 0.0
    @inbounds for j = 1:n_features, i = 1:n_features
        s_ij = S[i, j]
        normS_sq += s_ij * s_ij
        diff = s_ij - (i == j ? mu : 0.0)
        d2 += diff * diff
    end
    d2 /= n_features

    delta = if d2 <= 1e-30
        1.0
    else
        sum_norm4 = 0.0
        @inbounds for col = 1:total_trials
            n2 = 0.0
            @simd for f = 1:n_features
                v = ws.residuals[f, col]
                n2 += v * v
            end
            sum_norm4 += n2 * n2
        end
        b2 = min(d2, max(0.0, sum_norm4 - total_trials * normS_sq) / (Float64(total_trials)^2 * n_features))
        clamp(b2 / d2, 0.0, 1.0)
    end

    scale = total_trials / n_dof
    alpha_s = scale * (1.0 - delta)
    diag_add = scale * delta * mu + (mu <= 1e-15 ? 1e-12 : 0.0)
    @inbounds for j = 1:n_features, i = 1:n_features
        out_cov[i, j] = alpha_s * S[i, j] + (i == j ? diag_add : 0.0)
    end

    return out_cov
end

function _compute_inv_pooled_covariance!(
    out_inv::Matrix{Float64},
    ws::_PooledCovWorkspace,
    data_arrays::Vector{Array{Float64,3}},
    trial_indices::Union{Nothing,Vector{Vector{Int}}},
    time_idx::Int,
)
    _compute_pooled_covariance!(out_inv, ws, data_arrays, trial_indices, time_idx)
    chol = cholesky!(Symmetric(out_inv); check = false)
    if issuccess(chol)
        LinearAlgebra.inv!(chol)
    else
        _compute_pooled_covariance!(ws.S, ws, data_arrays, trial_indices, time_idx)
        copyto!(out_inv, pinv(ws.S))
    end
    return out_inv
end

function _compute_pooled_covariance_subset(
    data_arrays::Vector{Array{Float64,3}},
    trial_indices::Union{Nothing,Vector{Vector{Int}}},
    time_idx::Int,
)
    n_features = size(data_arrays[1], 1)
    max_trials = sum(size(d, 3) for d in data_arrays)
    ws = _PooledCovWorkspace(n_features, max_trials)
    out_cov = Matrix{Float64}(undef, n_features, n_features)
    return _compute_pooled_covariance!(out_cov, ws, data_arrays, trial_indices, time_idx)
end

function _compute_pooled_covariance(data_arrays::Vector{Array{Float64,3}}, time_idx::Int)
    return _compute_pooled_covariance_subset(data_arrays, nothing, time_idx)
end

function _compute_dissimilarity(
    pattern1::AbstractVector{Float64},
    pattern2::AbstractVector{Float64},
    measure::Symbol;
    inv_covariance_matrix::Union{Matrix{Float64},Nothing} = nothing,
)
    if measure == :correlation || measure == :pearson
        return 1.0 - cor(pattern1, pattern2)
    elseif measure == :spearman
        return 1.0 - StatsBase.corspearman(pattern1, pattern2)
    elseif measure == :euclidean
        d2 = 0.0
        @inbounds @simd for i in eachindex(pattern1)
            diff = pattern1[i] - pattern2[i]
            d2 += diff * diff
        end
        return sqrt(d2)
    elseif measure == :mahalanobis
        isnothing(inv_covariance_matrix) &&
            @minimal_error("Mahalanobis distance requires an inverse covariance matrix (provide inv_covariance_matrix).")
        n = length(pattern1)
        q = 0.0
        @inbounds for j = 1:n
            dj = pattern1[j] - pattern2[j]
            col_sum = 0.0
            @simd for i = 1:n
                col_sum += inv_covariance_matrix[i, j] * (pattern1[i] - pattern2[i])
            end
            q += dj * col_sum
        end
        return sqrt(max(0.0, q))
    else
        @minimal_error("Unknown dissimilarity measure: $measure. Use :correlation, :spearman, :euclidean, or :mahalanobis")
    end
end

function _compute_rdm!(
    rdm::Matrix{Float64},
    condition_patterns::Vector{<:AbstractVector{Float64}},
    measure::Symbol;
    inv_covariance_matrix::Union{Matrix{Float64},Nothing} = nothing,
)
    n_conditions = length(condition_patterns)
    @inbounds for i = 1:n_conditions
        rdm[i, i] = 0.0
        for j = (i+1):n_conditions
            dissim =
                _compute_dissimilarity(condition_patterns[i], condition_patterns[j], measure; inv_covariance_matrix = inv_covariance_matrix)
            rdm[i, j] = dissim
            rdm[j, i] = dissim
        end
    end
    return rdm
end

function _compute_rdm(
    condition_patterns::Vector{<:AbstractVector{Float64}},
    measure::Symbol;
    inv_covariance_matrix::Union{Matrix{Float64},Nothing} = nothing,
)
    n_conditions = length(condition_patterns)
    return _compute_rdm!(
        zeros(Float64, n_conditions, n_conditions),
        condition_patterns,
        measure;
        inv_covariance_matrix = inv_covariance_matrix,
    )
end

"""
    normalize_rdm(rdm::AbstractMatrix{Float64}; method::Symbol = :none) -> Matrix{Float64}

Normalize upper/lower off-diagonal entries of a symmetric RDM (preserving zero diagonal)
using `:none`, `:zscore`, `:rank` (average tied ranks), or `:minmax` (`[0, 1]`).
"""
function normalize_rdm(rdm::AbstractMatrix{Float64}; method::Symbol = :none)
    method == :none && return Matrix(rdm)

    n = size(rdm, 1)
    normalized_rdm = Matrix(rdm)
    upper_values = _extract_upper_triangular(rdm)

    normalized_values = if method == :zscore
        μ = mean(upper_values)
        σ = std(upper_values)
        σ > 0 ? (upper_values .- μ) ./ σ : zeros(Float64, length(upper_values))
    elseif method == :rank
        StatsBase.tiedrank(upper_values)
    elseif method == :minmax
        min_val, max_val = extrema(upper_values)
        max_val > min_val ? (upper_values .- min_val) ./ (max_val - min_val) : zeros(Float64, length(upper_values))
    else
        @minimal_error("Unknown normalization method: $method. Use :none, :zscore, :rank, or :minmax")
    end

    k = 1
    @inbounds for i = 1:n
        normalized_rdm[i, i] = 0.0
        for j = (i+1):n
            v = normalized_values[k]
            normalized_rdm[i, j] = v
            normalized_rdm[j, i] = v
            k += 1
        end
    end

    return normalized_rdm
end

function _normalize_rdms!(rdms::Array{Float64,3}, method::Symbol)
    method == :none && return rdms
    for t = 1:size(rdms, 1)
        rdms[t, :, :] = normalize_rdm(@view(rdms[t, :, :]); method = method)
    end
    return rdms
end

function _prepare_and_validate_rsa(epochs::AbstractVector{<:EpochData}, channel_selection, sample_selection)
    isempty(epochs) && @minimal_error("Cannot perform RSA with empty epochs vector")
    length(epochs) < 2 && @minimal_error("Need at least 2 conditions for RSA, got $(length(epochs))")

    first_ep = epochs[1]
    selected_channels = get_selected_channels(first_ep, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("Channel selection produced no channels")

    selected_samples = get_selected_samples(first_ep, sample_selection)
    isempty(selected_samples) && @minimal_error("Sample selection produced no time points")

    times = _extract_selected_f64(first_ep.data[1][!, :time], selected_samples)
    isempty(times) && @minimal_error("Sample selection produced no time points")

    data_arrays, n_trials_per_condition = _extract_epoch_arrays(epochs, selected_channels, selected_samples)
    condition_names = String[e.condition_name for e in epochs]
    layout_sub = subset_layout(first_ep.layout; channel_selection = channel_selection)

    return data_arrays, times, n_trials_per_condition, condition_names, selected_channels, layout_sub, first_ep
end

function _compute_rdms_from_data(
    data_arrays::Vector{Array{Float64,3}},
    n_timepoints::Int,
    n_conditions::Int,
    selected_channels::Vector{Symbol},
    dissimilarity_measure::Symbol,
)
    n_channels = length(selected_channels)
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    rdm_buffer = zeros(Float64, n_conditions, n_conditions)

    cond_means = [zeros(Float64, n_channels, n_timepoints) for _ = 1:n_conditions]
    for (i, cond_data) in enumerate(data_arrays)
        n_trials = size(cond_data, 3)
        cm = cond_means[i]
        if n_trials > 0
            @inbounds for trial = 1:n_trials, t = 1:n_timepoints
                @simd for ch = 1:n_channels
                    cm[ch, t] += cond_data[ch, t, trial]
                end
            end
            cm ./= n_trials
        end
    end

    use_mahal = dissimilarity_measure == :mahalanobis
    cov_ws = use_mahal ? _PooledCovWorkspace(n_channels, sum(size(d, 3) for d in data_arrays)) : nothing
    inv_pooled_cov = use_mahal ? Matrix{Float64}(undef, n_channels, n_channels) : nothing
    condition_patterns = [Vector{Float64}(undef, n_channels) for _ = 1:n_conditions]

    for t = 1:n_timepoints
        if use_mahal
            _compute_inv_pooled_covariance!(inv_pooled_cov::Matrix{Float64}, cov_ws::_PooledCovWorkspace, data_arrays, nothing, t)
        end

        @inbounds for i = 1:n_conditions
            cp = condition_patterns[i]
            cm = cond_means[i]
            @simd for ch = 1:n_channels
                cp[ch] = cm[ch, t]
            end
        end

        _compute_rdm!(rdm_buffer, condition_patterns, dissimilarity_measure; inv_covariance_matrix = inv_pooled_cov)
        @inbounds for j = 1:n_conditions, i = 1:n_conditions
            rdms[t, i, j] = rdm_buffer[i, j]
        end
    end

    return rdms
end

"""
    rsa(
        epochs::AbstractVector{<:EpochData};
        channel_selection = channels(),
        sample_selection = samples(),
        interval_selection::Interval = times(),
        dissimilarity_measure::Symbol = :correlation,
        average_trials::Bool = true,
        normalize_method::Symbol = :none,
    ) -> RsaData
    rsa(all_participant_epochs::Vector{Vector{EpochData}}; kwargs...) -> Vector{RsaData}

Compute time-resolved Representational Dissimilarity Matrices (RDMs) across conditions for a
participant (or vector of participants).

# Keyword Arguments
- `channel_selection`: Channel selection predicate (default: `channels()`)
- `sample_selection`: Sample selection predicate (default: `samples()`)
- `interval_selection::Interval`: Time interval filter (default: `times()`)
- `dissimilarity_measure::Symbol`: `:correlation` (`:pearson`), `:spearman`, `:euclidean`, or `:mahalanobis`
- `average_trials::Bool`: Average trials per condition before computing RDMs (default: `true`)
- `normalize_method::Symbol`: `:none` (default), `:zscore`, `:rank`, or `:minmax`
"""
function rsa(
    epochs::AbstractVector{<:EpochData};
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    dissimilarity_measure::Symbol = :correlation,
    average_trials::Bool = true,
    normalize_method::Symbol = :none,
)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    data_arrays, times, n_trials_per_condition, condition_names, selected_channels, layout_sub, first_epoch =
        _prepare_and_validate_rsa(epochs, channel_selection, combined_sel)

    n_conditions = length(condition_names)
    n_timepoints = length(times)

    if any(==(0), n_trials_per_condition)
        empty_conditions = condition_names[n_trials_per_condition .== 0]
        @minimal_error("One or more conditions have zero trials: $(join(empty_conditions, ", "))")
    end

    if average_trials
        rdms = _compute_rdms_from_data(data_arrays, n_timepoints, n_conditions, selected_channels, dissimilarity_measure)
    else
        n_ch = length(selected_channels)
        rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)
        avg_rdm = zeros(Float64, n_conditions, n_conditions)
        trial_rdm_buffer = zeros(Float64, n_conditions, n_conditions)
        condition_patterns = [Vector{Float64}(undef, n_ch) for _ = 1:n_conditions]
        n_trials = minimum(n_trials_per_condition)
        use_mahal = dissimilarity_measure == :mahalanobis
        cov_ws = use_mahal ? _PooledCovWorkspace(n_ch, sum(size(d, 3) for d in data_arrays)) : nothing
        inv_pooled_cov = use_mahal ? Matrix{Float64}(undef, n_ch, n_ch) : nothing

        for t = 1:n_timepoints
            if use_mahal
                _compute_inv_pooled_covariance!(inv_pooled_cov::Matrix{Float64}, cov_ws::_PooledCovWorkspace, data_arrays, nothing, t)
            end
            fill!(avg_rdm, 0.0)
            for trial_idx = 1:n_trials
                @inbounds for cond_idx = 1:n_conditions
                    cp = condition_patterns[cond_idx]
                    cd = data_arrays[cond_idx]
                    @simd for ch = 1:n_ch
                        cp[ch] = cd[ch, t, trial_idx]
                    end
                end
                _compute_rdm!(trial_rdm_buffer, condition_patterns, dissimilarity_measure; inv_covariance_matrix = inv_pooled_cov)
                avg_rdm .+= trial_rdm_buffer
            end
            if n_trials > 0
                inv_nt = 1.0 / n_trials
                @inbounds for j = 1:n_conditions, i = 1:n_conditions
                    rdms[t, i, j] = avg_rdm[i, j] * inv_nt
                end
            end
        end
    end

    _normalize_rdms!(rdms, normalize_method)

    return RsaData(
        first_epoch.file,
        condition_names,
        times,
        rdms,
        dissimilarity_measure,
        selected_channels,
        layout_sub,
        first_epoch.sample_rate;
        analysis_info = copy(first_epoch.analysis_info),
    )
end

function rsa(all_participant_epochs::Vector{Vector{EpochData}}; kwargs...)
    return [rsa(participant_epochs; kwargs...) for participant_epochs in all_participant_epochs]
end

# === MODEL COMPARISON & MANTEL PERMUTATIONS ===

function _validate_model_rdms(model_rdms, n_conditions::Int, n_times::Int)
    is_temporal = Vector{Bool}(undef, length(model_rdms))
    for (idx, model_rdm) in enumerate(model_rdms)
        if model_rdm isa AbstractArray{<:Real,3}
            size(model_rdm, 1) != n_times &&
                @minimal_error("Temporal model RDM $idx has $(size(model_rdm, 1)) timepoints, expected $n_times to match neural data")
            (size(model_rdm, 2) != n_conditions || size(model_rdm, 3) != n_conditions) && @minimal_error(
                "Temporal model RDM $idx has condition dimensions $(size(model_rdm, 2))×$(size(model_rdm, 3)), expected $n_conditions×$n_conditions"
            )
            is_temporal[idx] = true
        elseif model_rdm isa AbstractMatrix{<:Real}
            size(model_rdm) != (n_conditions, n_conditions) &&
                @minimal_error("Static model RDM $idx has size $(size(model_rdm)), expected ($n_conditions, $n_conditions)")
            issymmetric(model_rdm) || @minimal_warning "Static model RDM $idx is not symmetric, symmetrizing copy"
            is_temporal[idx] = false
        else
            @minimal_error(
                "Model RDM $idx has unsupported type $(typeof(model_rdm)). Expected Matrix{Float64} (static) or Array{Float64, 3} (temporal)"
            )
        end
    end
    return is_temporal
end

function _all_permutations(n::Int)
    perms = sizehint!(Vector{Vector{Int}}(), factorial(n))
    cur = collect(1:n)
    function _gen!(k::Int)
        if k == 1
            push!(perms, copy(cur))
            return
        end
        _gen!(k - 1)
        for i = 1:(k-1)
            swap_idx = iseven(k) ? i : 1
            cur[swap_idx], cur[k] = cur[k], cur[swap_idx]
            _gen!(k - 1)
        end
    end
    _gen!(n)
    return perms
end

function _correlate_vectors(vec1, vec2, correlation_type::Symbol)
    correlation_type == :spearman && return StatsBase.corspearman(vec1, vec2)
    correlation_type == :pearson && return cor(vec1, vec2)
    @minimal_error("Unknown correlation type: $correlation_type. Use :spearman or :pearson")
end

function _fill_mantel_unit_matrix!(
    U::Matrix{Float64},
    u_buf::Vector{Float64},
    rank_buf::Vector{Float64},
    rank_perm::Vector{Int},
    model_vec::AbstractVector{Float64},
    correlation_type::Symbol,
)
    src = correlation_type == :spearman ? _tiedrank!(rank_buf, rank_perm, model_vec) : model_vec
    has_var = _zscore_unit_norm!(u_buf, src)
    N = size(U, 1)
    k = 1
    @inbounds for i = 1:N
        U[i, i] = 0.0
        for j = (i+1):N
            v = u_buf[k]
            U[i, j] = v
            U[j, i] = v
            k += 1
        end
    end
    return has_var
end

@inline function _eval_mantel_dot(u_neural::Vector{Float64}, U_model::Matrix{Float64}, perm::Vector{Int})
    N = length(perm)
    s = 0.0
    idx = 1
    @inbounds for i = 1:N
        pi = perm[i]
        for j = (i+1):N
            s += u_neural[idx] * U_model[pi, perm[j]]
            idx += 1
        end
    end
    return clamp(s, -1.0, 1.0)
end

function _mantel_pvalue(
    u_neural::Vector{Float64},
    U_model::Matrix{Float64},
    observed_corr::Float64,
    n_permutations::Int,
    tail::Symbol,
    exact_perms::Union{Nothing,Vector{Vector{Int}}},
    perm_buf::Vector{Int},
    rng::AbstractRNG,
)
    tol = 1e-12
    abs_obs = abs(observed_corr)
    @inline is_ext(c::Float64) =
        tail == :right ? (c >= observed_corr - tol) : tail == :left ? (c <= observed_corr + tol) : (abs(c) >= abs_obs - tol)

    if !isnothing(exact_perms)
        extreme_count = 0
        for perm in exact_perms
            extreme_count += is_ext(_eval_mantel_dot(u_neural, U_model, perm))
        end
        return max(extreme_count, 1) / length(exact_perms)
    end

    extreme_count = 0
    for _ = 1:n_permutations
        Random.randperm!(rng, perm_buf)
        extreme_count += is_ext(_eval_mantel_dot(u_neural, U_model, perm_buf))
    end
    return (extreme_count + 1) / (n_permutations + 1)
end

function _compute_model_correlations(
    rsa_data::RsaData,
    model_rdms,
    is_temporal::Vector{Bool},
    correlation_type::Symbol,
    n_permutations::Int;
    tail::Symbol = :right,
    rng::AbstractRNG = Random.default_rng(),
)
    n_times = length(rsa_data.times)
    n_models = length(model_rdms)
    n_conditions = length(rsa_data.condition_names)
    n_upper = (n_conditions * (n_conditions - 1)) ÷ 2

    correlations = zeros(Float64, n_times, n_models)
    p_values = n_permutations > 0 ? zeros(Float64, n_times, n_models) : nothing

    static_model_vecs = Vector{Union{Nothing,Vector{Float64}}}(nothing, n_models)
    temporal_model_vecs = Vector{Union{Nothing,Matrix{Float64}}}(nothing, n_models)
    static_U_models = Vector{Union{Nothing,Matrix{Float64}}}(nothing, n_models)
    static_has_var = falses(n_models)

    rank_buf = Vector{Float64}(undef, n_upper)
    rank_perm = Vector{Int}(undef, n_upper)
    u_buf = Vector{Float64}(undef, n_upper)

    for model_idx = 1:n_models
        model_rdm = model_rdms[model_idx]
        if is_temporal[model_idx]
            t_vecs = Matrix{Float64}(undef, n_upper, n_times)
            for t = 1:n_times
                _symmetrized_upper_triangular!(@view(t_vecs[:, t]), @view(model_rdm[t, :, :]))
            end
            temporal_model_vecs[model_idx] = t_vecs
        else
            m_vec = _symmetrized_upper_triangular(model_rdm)
            static_model_vecs[model_idx] = m_vec
            if n_permutations > 0 && n_conditions >= 3 && n_upper >= 2
                U = zeros(Float64, n_conditions, n_conditions)
                static_has_var[model_idx] = _fill_mantel_unit_matrix!(U, u_buf, rank_buf, rank_perm, m_vec, correlation_type)
                static_U_models[model_idx] = U
            end
        end
    end

    exact_perms = (n_permutations > 0 && 3 <= n_conditions <= 4) ? _all_permutations(n_conditions) : nothing
    perm_buf = Vector{Int}(undef, n_conditions)
    neural_vec = Vector{Float64}(undef, n_upper)
    u_neural = Vector{Float64}(undef, n_upper)
    U_temp = (n_permutations > 0 && any(is_temporal)) ? zeros(Float64, n_conditions, n_conditions) : nothing

    for t = 1:n_times
        _extract_upper_triangular!(neural_vec, @view(rsa_data.rdm[t, :, :]))
        neural_has_var = false
        if n_permutations > 0 && n_conditions >= 3
            n_src = correlation_type == :spearman ? _tiedrank!(rank_buf, rank_perm, neural_vec) : neural_vec
            neural_has_var = _zscore_unit_norm!(u_neural, n_src)
        end

        for model_idx = 1:n_models
            is_temp = is_temporal[model_idx]
            model_vec =
                is_temp ? @view((temporal_model_vecs[model_idx]::Matrix{Float64})[:, t]) : static_model_vecs[model_idx]::Vector{Float64}

            # Observe correlation using pre-computed unit vectors when valid, or fallback to _correlate_vectors for exact NaN behavior
            corr = _correlate_vectors(neural_vec, model_vec, correlation_type)
            correlations[t, model_idx] = corr

            if n_permutations > 0
                p_val = if isnan(corr) || n_conditions < 3 || !neural_has_var
                    1.0
                elseif is_temp
                    U_t = U_temp::Matrix{Float64}
                    _fill_mantel_unit_matrix!(U_t, u_buf, rank_buf, rank_perm, model_vec, correlation_type) ?
                    _mantel_pvalue(u_neural, U_t, corr, n_permutations, tail, exact_perms, perm_buf, rng) : 1.0
                else
                    static_has_var[model_idx] ?
                    _mantel_pvalue(
                        u_neural,
                        static_U_models[model_idx]::Matrix{Float64},
                        corr,
                        n_permutations,
                        tail,
                        exact_perms,
                        perm_buf,
                        rng,
                    ) : 1.0
                end
                (p_values::Matrix{Float64})[t, model_idx] = p_val
            end
        end
    end

    return correlations, p_values
end

"""
    compare_models!(
        rsa_data::RsaData,
        model_rdms::AbstractVector;
        model_names::Union{Vector{String},Nothing} = nothing,
        correlation_type::Symbol = :spearman,
        n_permutations::Int = 1000,
        tail::Symbol = :right,
        rng::AbstractRNG = Random.default_rng(),
    ) -> RsaData

Compare neural RDMs in `rsa_data` to `model_rdms` in-place (without mutating `model_rdms`).
See [`compare_models`](@ref) for full parameter details.
"""
function compare_models!(
    rsa_data::RsaData,
    model_rdms::AbstractVector;
    model_names::Union{Vector{String},Nothing} = nothing,
    correlation_type::Symbol = :spearman,
    n_permutations::Int = 1000,
    tail::Symbol = :right,
    rng::AbstractRNG = Random.default_rng(),
)
    correlation_type in (:spearman, :pearson) || @minimal_error("Unknown correlation type: $correlation_type. Use :spearman or :pearson")
    tail in (:right, :left, :both) || @minimal_error("Unknown tail: $tail. Use :right, :left, or :both")
    n_permutations < 0 && @minimal_error("n_permutations must be >= 0, got $n_permutations")

    n_models = length(model_rdms)
    n_conditions = length(rsa_data.condition_names)
    n_times = length(rsa_data.times)
    is_temporal = _validate_model_rdms(model_rdms, n_conditions, n_times)

    resolved_names = if isnothing(model_names)
        ["Model$(i)" for i = 1:n_models]
    elseif length(model_names) != n_models
        @minimal_error("Number of model names ($(length(model_names))) doesn't match number of models ($n_models)")
    else
        copy(model_names)
    end

    correlations, p_values =
        _compute_model_correlations(rsa_data, model_rdms, is_temporal, correlation_type, n_permutations; tail = tail, rng = rng)

    rsa_data.model_correlations = correlations
    rsa_data.model_names = resolved_names
    rsa_data.p_values = p_values
    return rsa_data
end

"""
    compare_models(
        rsa_data::RsaData,
        model_rdms::AbstractVector;
        model_names::Union{Vector{String},Nothing} = nothing,
        correlation_type::Symbol = :spearman,
        n_permutations::Int = 1000,
        tail::Symbol = :right,
        rng::AbstractRNG = Random.default_rng(),
    ) -> RsaData

Compare neural RDMs to static (`Matrix{Float64}`) or temporal (`Array{Float64,3}`) model RDMs,
returning a new `RsaData` copy with Mantel condition-permutation `p_values` (exact `N!` for
`3 <= N <= 4`, Monte Carlo for `N >= 5`).
"""
function compare_models(
    rsa_data::RsaData,
    model_rdms::AbstractVector;
    model_names::Union{Vector{String},Nothing} = nothing,
    correlation_type::Symbol = :spearman,
    n_permutations::Int = 1000,
    tail::Symbol = :right,
    rng::AbstractRNG = Random.default_rng(),
)
    return compare_models!(
        copy(rsa_data),
        model_rdms;
        model_names = model_names,
        correlation_type = correlation_type,
        n_permutations = n_permutations,
        tail = tail,
        rng = rng,
    )
end

function grand_average(rsa_data_list::AbstractVector{<:RsaData}; compute_noise_ceiling::Bool = true)
    isempty(rsa_data_list) && @minimal_error("Cannot compute grand average with empty RSA data list")
    length(rsa_data_list) == 1 && return copy(rsa_data_list[1])

    first_rsa = rsa_data_list[1]
    n_conditions = length(first_rsa.condition_names)
    n_times = length(first_rsa.times)

    for (idx, rsa_data) in enumerate(rsa_data_list)
        length(rsa_data.condition_names) != n_conditions &&
            @minimal_error("RSA data $idx has $(length(rsa_data.condition_names)) conditions, expected $n_conditions")
        length(rsa_data.times) != n_times && @minimal_error("RSA data $idx has $(length(rsa_data.times)) time points, expected $n_times")
        rsa_data.dissimilarity_measure != first_rsa.dissimilarity_measure &&
            @minimal_warning "RSA data $idx uses dissimilarity measure $(rsa_data.dissimilarity_measure), expected $(first_rsa.dissimilarity_measure)"
    end

    grand_rdm = zeros(Float64, n_times, n_conditions, n_conditions)
    for rsa_data in rsa_data_list
        grand_rdm .+= rsa_data.rdm
    end
    grand_rdm ./= length(rsa_data_list)

    noise_ceiling = (compute_noise_ceiling && length(rsa_data_list) >= 2) ? EegFun.compute_noise_ceiling(rsa_data_list) : nothing

    return RsaData(
        "Grand Average",
        copy(first_rsa.condition_names),
        copy(first_rsa.times),
        grand_rdm,
        first_rsa.dissimilarity_measure,
        copy(first_rsa.channels),
        copy(first_rsa.layout),
        first_rsa.sample_rate;
        noise_ceiling = noise_ceiling,
        analysis_info = copy(first_rsa.analysis_info),
    )
end
