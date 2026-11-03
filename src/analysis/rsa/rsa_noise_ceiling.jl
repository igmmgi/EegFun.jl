# Noise ceiling estimation for Representational Similarity Analysis (RSA).
#
# Implements noise ceiling computation using leave-one-out cross-validation
# across participants, following the methodology of Nili et al. (2014).
#
# References
# Nili, H., Wingfield, C., Walther, A., Su, L., Marslen-Wilson, W., & Kriegeskorte, N. (2014).
# A toolbox for representational similarity analysis. PLoS computational biology, 10(4), e1003553.

# === HELPER FUNCTIONS ===

"""
    _correlate_rdms(rdm1::AbstractMatrix{Float64}, rdm2::AbstractMatrix{Float64}; method::Symbol = :spearman)

Compute correlation between two RDMs using upper triangle values.
"""
function _correlate_rdms(rdm1::AbstractMatrix{Float64}, rdm2::AbstractMatrix{Float64}; method::Symbol = :spearman)
    if size(rdm1) != size(rdm2)
        @minimal_error("RDMs must have the same dimensions, got $(size(rdm1)) and $(size(rdm2))")
    end
    return _correlate_vectors(_extract_upper_triangular(rdm1), _extract_upper_triangular(rdm2), method)
end

function _zscore_unit_norm!(dest::Vector{Float64}, src::AbstractVector{Float64})
    n = length(src)
    if n < 2
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

# === NOISE CEILING COMPUTATION ===

"""
    compute_noise_ceiling(
        rsa_data_list::AbstractVector{<:RsaData};
        correlation_method::Symbol = :spearman,
    )

Compute noise ceiling for RSA analysis using leave-one-out cross-validation (Nili et al., 2014).

Before averaging across participants at each time point, each participant's upper-triangular RDM
vector is rank-transformed via `StatsBase.tiedrank` (for `:spearman`) and z-scored (for both
`:spearman` and `:pearson`) so that the group average is the exact closed-form maximizer of mean
correlation and `lower_bound <= upper_bound` holds strictly.

# Arguments
- `rsa_data_list::AbstractVector{<:RsaData}`: Vector of `RsaData` objects, one per participant (`>= 2`)
- `correlation_method::Symbol`: `:spearman` (default) or `:pearson`

# Returns
- `noise_ceiling::NoiseCeiling`: Noise ceiling estimates with `lower_bound` and `upper_bound`

# Examples
```julia
all_rsa = [EegFun.rsa(epochs_p1), EegFun.rsa(epochs_p2), EegFun.rsa(epochs_p3)]
nc = EegFun.compute_noise_ceiling(all_rsa)
grand_avg = EegFun.grand_average(all_rsa)
```
"""
function compute_noise_ceiling(rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)
    if length(rsa_data_list) < 2
        @minimal_error("Noise ceiling requires at least 2 participants, got $(length(rsa_data_list))")
    end
    if !(correlation_method in (:spearman, :pearson))
        @minimal_error("Unknown correlation method: $correlation_method. Use :spearman or :pearson")
    end

    first_rsa = rsa_data_list[1]
    n_conditions = length(first_rsa.condition_names)
    n_times = length(first_rsa.times)

    for (idx, rsa_data) in enumerate(rsa_data_list)
        if length(rsa_data.condition_names) != n_conditions
            @minimal_error("Participant $idx has $(length(rsa_data.condition_names)) conditions, expected $n_conditions")
        end
        if length(rsa_data.times) != n_times
            @minimal_error("Participant $idx has $(length(rsa_data.times)) time points, expected $n_times")
        end
    end

    n_participants = length(rsa_data_list)
    lower_bound = zeros(Float64, n_times)
    upper_bound = zeros(Float64, n_times)

    n_upper = (n_conditions * (n_conditions - 1)) ÷ 2
    if n_upper < 2
        # With only 2 conditions (1 upper-triangular element), correlation across pairs is undefined (0.0)
        return NoiseCeiling(lower_bound, upper_bound, n_participants)
    end

    raw_vec = Vector{Float64}(undef, n_upper)
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
            rdm_p = @view rsa_data_list[p].rdm[t, :, :]
            k = 1
            @inbounds for i = 1:n_conditions
                for j = (i+1):n_conditions
                    raw_vec[k] = rdm_p[i, j]
                    k += 1
                end
            end

            transformed = correlation_method == :spearman ? StatsBase.tiedrank(raw_vec) : raw_vec
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
    add_noise_ceiling!(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)

Compute and attach a `NoiseCeiling` to `rsa_data` in-place.
"""
function add_noise_ceiling!(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)
    rsa_data.noise_ceiling = compute_noise_ceiling(rsa_data_list; correlation_method = correlation_method)
    return rsa_data
end

"""
    add_noise_ceiling(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)

Return a copy of `rsa_data` with the computed `NoiseCeiling` attached (leaving `rsa_data` unmodified).
"""
function add_noise_ceiling(rsa_data::RsaData, rsa_data_list::AbstractVector{<:RsaData}; correlation_method::Symbol = :spearman)
    return add_noise_ceiling!(copy(rsa_data), rsa_data_list; correlation_method = correlation_method)
end
