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
    _correlate_rdms(rdm1::Matrix{Float64}, rdm2::Matrix{Float64}; method::Symbol = :spearman)

Compute correlation between two RDMs using upper triangle values.

# Arguments
- `rdm1::Matrix{Float64}`: First RDM [condition × condition]
- `rdm2::Matrix{Float64}`: Second RDM [condition × condition]
- `method::Symbol`: Correlation method (:spearman or :pearson)

# Returns
- `correlation::Float64`: Correlation coefficient between the two RDMs
"""
function _correlate_rdms(rdm1::AbstractMatrix{Float64}, rdm2::AbstractMatrix{Float64}; method::Symbol = :spearman)
    # Validate dimensions
    if size(rdm1) != size(rdm2)
        @minimal_error("RDMs must have the same dimensions, got $(size(rdm1)) and $(size(rdm2))")
    end

    return _correlate_vectors(_extract_upper_triangular(rdm1), _extract_upper_triangular(rdm2), method)
end

# === NOISE CEILING COMPUTATION ===

"""
    compute_noise_ceiling(
        rsa_data_list::Vector{RsaData};
        correlation_method::Symbol = :spearman,
    )

Compute noise ceiling for RSA analysis using leave-one-out cross-validation.

The noise ceiling quantifies the maximum correlation that an ideal model could
achieve with the observed neural data, given the noise and variability in the
measurements. It provides a benchmark for evaluating model performance.

# Methodology (Nili et al., 2014)

**Lower Bound** (conservative estimate):
- For each participant, correlate their RDM with the average RDM of all OTHER participants
- Average these correlations across participants
- This underestimates the true ceiling because it doesn't use all available data

**Upper Bound** (liberal estimate):
- For each participant, correlate their RDM with the average RDM of ALL participants (including themselves)
- Average these correlations across participants
- This overestimates the true ceiling because it includes the participant in their own reference

The true noise ceiling likely lies between these bounds.

# Arguments
- `rsa_data_list::Vector{RsaData}`: Vector of RsaData objects, one per participant
  - All must have the same number of conditions and time points
  - All must use the same dissimilarity measure
- `correlation_method::Symbol`: Method for correlating RDMs (:spearman or :pearson, default: :spearman)

# Returns
- `noise_ceiling::NoiseCeiling`: Noise ceiling estimates with lower and upper bounds

# Examples
```julia
# Compute RSA for each participant
all_rsa = [rsa(epochs_p1), rsa(epochs_p2), rsa(epochs_p3)]

# Compute noise ceiling
nc = compute_noise_ceiling(all_rsa)

# Noise ceiling is automatically added to grand average
grand_avg = grand_average(all_rsa)
println(grand_avg.noise_ceiling)
```

# Notes
- Requires at least 2 participants (preferably 10+)
- More participants = more reliable noise ceiling estimate
- The noise ceiling is time-varying (computed at each time point)
- Use this to evaluate whether your model performance is limited by the model or by data quality

# References
Nili, H., Wingfield, C., Walther, A., Su, L., Marslen-Wilson, W., & Kriegeskorte, N. (2014).
A toolbox for representational similarity analysis. PLoS computational biology, 10(4), e1003553.
"""
function compute_noise_ceiling(rsa_data_list::Vector{RsaData}; correlation_method::Symbol = :spearman)
    # Validate input
    if length(rsa_data_list) < 2
        @minimal_error("Noise ceiling requires at least 2 participants, got $(length(rsa_data_list))")
    end

    # Validate all participants have same structure
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

    # Preallocate arrays for bounds and upper-triangle vectors
    lower_bound = zeros(Float64, n_times)
    upper_bound = zeros(Float64, n_times)

    n_upper = (n_conditions * (n_conditions - 1)) ÷ 2
    participant_vecs = [Vector{Float64}(undef, n_upper) for _ = 1:n_participants]
    sum_vec = Vector{Float64}(undef, n_upper)
    avg_vec_all = Vector{Float64}(undef, n_upper)
    avg_vec_others = Vector{Float64}(undef, n_upper)
    inv_n_minus_1 = 1.0 / (n_participants - 1)

    # Compute noise ceiling at each time point
    for t = 1:n_times
        fill!(sum_vec, 0.0)

        # Extract upper-triangular elements for all participants and accumulate sum
        for p = 1:n_participants
            v = participant_vecs[p]
            rdm_p = @view rsa_data_list[p].rdm[t, :, :]
            k = 1
            @inbounds for i = 1:n_conditions
                for j = (i+1):n_conditions
                    val = rdm_p[i, j]
                    v[k] = val
                    sum_vec[k] += val
                    k += 1
                end
            end
        end

        # Grand average upper triangle across all participants
        @inbounds @simd for k = 1:n_upper
            avg_vec_all[k] = sum_vec[k] / n_participants
        end

        # Leave-one-out cross-validation
        sum_lower = 0.0
        sum_upper = 0.0

        for left_out_idx = 1:n_participants
            left_vec = participant_vecs[left_out_idx]

            # Leave-one-out mean: (sum_all - left_out) / (n - 1)
            @inbounds @simd for k = 1:n_upper
                avg_vec_others[k] = (sum_vec[k] - left_vec[k]) * inv_n_minus_1
            end

            sum_lower += _correlate_vectors(left_vec, avg_vec_others, correlation_method)
            sum_upper += _correlate_vectors(left_vec, avg_vec_all, correlation_method)
        end

        lower_bound[t] = sum_lower / n_participants
        upper_bound[t] = sum_upper / n_participants
    end

    return NoiseCeiling(lower_bound, upper_bound, n_participants)
end

"""
    add_noise_ceiling!(rsa_data::RsaData, rsa_data_list::Vector{RsaData}; correlation_method::Symbol = :spearman)

Compute and add noise ceiling to an RsaData object (typically grand average).

This is a convenience function that computes the noise ceiling and adds it to
the provided RsaData object.

# Arguments
- `rsa_data::RsaData`: RsaData object to add noise ceiling to (modified in-place)
- `rsa_data_list::Vector{RsaData}`: Vector of individual participant RsaData objects
- `correlation_method::Symbol`: Correlation method (:spearman or :pearson)

# Returns
- `rsa_data::RsaData`: The same RsaData object with noise_ceiling field populated

# Examples
```julia
# Compute RSA for each participant
all_rsa = [rsa(epochs_p1), rsa(epochs_p2), rsa(epochs_p3)]

# Create grand average
grand_avg = grand_average(all_rsa)

# Add noise ceiling
add_noise_ceiling!(grand_avg, all_rsa)
```
"""
function add_noise_ceiling!(rsa_data::RsaData, rsa_data_list::Vector{RsaData}; correlation_method::Symbol = :spearman)
    rsa_data.noise_ceiling = compute_noise_ceiling(rsa_data_list; correlation_method = correlation_method)
    return rsa_data
end
