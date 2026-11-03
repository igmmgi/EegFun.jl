# Leave-one-out noise ceiling estimation for RSA (Nili et al., 2014).

function _safe_pearson_centered(unit_x::Vector{Float64}, has_var_x::Bool, y::Vector{Float64})
    (!has_var_x || length(unit_x) < 2) && return 0.0
    μ_y = mean(y)
    num = 0.0
    norm_y_sq = 0.0
    @inbounds @simd for k in eachindex(unit_x)
        dy = y[k] - μ_y
        num += unit_x[k] * dy
        norm_y_sq += dy * dy
    end
    norm_y_sq <= 1e-30 && return 0.0
    return clamp(num / sqrt(norm_y_sq), -1.0, 1.0)
end

"""
    compute_noise_ceiling(
        rsa_data_list::AbstractVector{<:RsaData};
        correlation_method::Symbol = :spearman,
    ) -> NoiseCeiling

Compute leave-one-out lower and upper RSA noise ceiling bounds across `>= 2` participants
(Nili et al., 2014), rank-transforming (`:spearman`) and z-scoring each participant's
upper-triangular RDM vector before group averaging so `lower_bound <= upper_bound` holds strictly.
"""
function compute_noise_ceiling(rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)
    length(rsa_data_list) < 2 && @minimal_error("Noise ceiling requires at least 2 participants, got $(length(rsa_data_list))")
    correlation_method in (:spearman, :pearson) ||
        @minimal_error("Unknown correlation method: $correlation_method. Use :spearman or :pearson")

    first_rsa = rsa_data_list[1]
    n_conditions = length(first_rsa.condition_names)
    n_times = length(first_rsa.times)

    for (idx, rsa_data) in enumerate(rsa_data_list)
        length(rsa_data.condition_names) != n_conditions &&
            @minimal_error("Participant $idx has $(length(rsa_data.condition_names)) conditions, expected $n_conditions")
        length(rsa_data.times) != n_times && @minimal_error("Participant $idx has $(length(rsa_data.times)) time points, expected $n_times")
    end

    n_participants = length(rsa_data_list)
    lower_bound = zeros(Float64, n_times)
    upper_bound = zeros(Float64, n_times)

    n_upper = (n_conditions * (n_conditions - 1)) ÷ 2
    n_upper < 2 && return NoiseCeiling(lower_bound, upper_bound, n_participants)

    raw_vec = Vector{Float64}(undef, n_upper)
    rank_buf = Vector{Float64}(undef, n_upper)
    rank_perm = Vector{Int}(undef, n_upper)
    participant_vecs = [Vector{Float64}(undef, n_upper) for _ = 1:n_participants]
    has_var = Vector{Bool}(undef, n_participants)
    sum_vec = Vector{Float64}(undef, n_upper)
    avg_vec_all = Vector{Float64}(undef, n_upper)
    avg_vec_others = Vector{Float64}(undef, n_upper)
    inv_n = 1.0 / n_participants
    inv_n_minus_1 = 1.0 / (n_participants - 1)

    for t = 1:n_times
        fill!(sum_vec, 0.0)

        for p = 1:n_participants
            _extract_upper_triangular!(raw_vec, @view(rsa_data_list[p].rdm[t, :, :]))
            transformed = correlation_method == :spearman ? _tiedrank!(rank_buf, rank_perm, raw_vec) : raw_vec
            has_var[p] = _zscore_unit_norm!(participant_vecs[p], transformed)

            v = participant_vecs[p]
            @inbounds @simd for idx = 1:n_upper
                sum_vec[idx] += v[idx]
            end
        end

        @inbounds @simd for k = 1:n_upper
            avg_vec_all[k] = sum_vec[k] * inv_n
        end

        sum_lower = 0.0
        sum_upper = 0.0
        for left_out_idx = 1:n_participants
            left_vec = participant_vecs[left_out_idx]
            @inbounds @simd for k = 1:n_upper
                avg_vec_others[k] = (sum_vec[k] - left_vec[k]) * inv_n_minus_1
            end
            sum_lower += _safe_pearson_centered(left_vec, has_var[left_out_idx], avg_vec_others)
            sum_upper += _safe_pearson_centered(left_vec, has_var[left_out_idx], avg_vec_all)
        end

        lower_bound[t] = sum_lower * inv_n
        upper_bound[t] = sum_upper * inv_n
    end

    return NoiseCeiling(lower_bound, upper_bound, n_participants)
end

"""
    add_noise_ceiling!(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman) -> RsaData

Compute and attach a [`NoiseCeiling`](@ref) to `rsa_data` in-place.
"""
function add_noise_ceiling!(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)
    rsa_data.noise_ceiling = compute_noise_ceiling(rsa_data_list; correlation_method = correlation_method)
    return rsa_data
end

"""
    add_noise_ceiling(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman) -> RsaData

Return a copy of `rsa_data` with the computed [`NoiseCeiling`](@ref) attached.
"""
function add_noise_ceiling(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)
    return add_noise_ceiling!(copy(rsa_data), rsa_data_list; correlation_method = correlation_method)
end
