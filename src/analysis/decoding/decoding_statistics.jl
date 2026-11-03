"""
    TemporalCluster

Stores information about a contiguous 1D temporal cluster in decoding results (`id`, `time_indices`,
`time_range`, `cluster_stat`, `p_value`, `is_significant`).
"""
struct TemporalCluster
    id::Int
    time_indices::Vector{Int}
    time_range::Tuple{Float64,Float64}
    cluster_stat::Float64
    p_value::Float64
    is_significant::Bool
end

"""
    DecodingStatisticsResult

Stores statistical test results on multi-participant decoding accuracy (`times`, `t_statistics`,
`p_values`, `significant_mask`, `df`, `alpha`, `correction_method`, `clusters`).
"""
struct DecodingStatisticsResult
    times::Vector{Float64}
    t_statistics::Vector{Float64}
    p_values::Vector{Float64}
    significant_mask::BitVector
    df::Float64
    alpha::Float64
    correction_method::Symbol
    clusters::Union{Vector{TemporalCluster},Nothing}
end

function _validate_and_extract_diffs(decoded_list::AbstractVector{<:DecodedData}, fn_name::AbstractString, alpha::Real)
    isempty(decoded_list) && @minimal_error("Cannot test empty decoded data list")
    length(decoded_list) < 2 && @minimal_error("$fn_name requires at least 2 participants, got $(length(decoded_list))")
    alpha_f = Float64(alpha)
    (alpha_f <= 0.0 || alpha_f >= 1.0) && @minimal_error("alpha must be between 0 and 1, got $alpha_f")

    first_decoded = decoded_list[1]
    first_times = first_decoded.times
    chance_level = first_decoded.parameters.chance_level

    for decoded in @view decoded_list[2:end]
        decoded.times != first_times && @minimal_error("DecodedData objects have inconsistent time vectors")
        decoded.parameters.chance_level != chance_level && @minimal_warning "DecodedData objects have different chance_level"
    end

    n_participants = length(decoded_list)
    n_timepoints = length(first_times)
    diffs = Matrix{Float64}(undef, n_participants, n_timepoints)
    for (p_idx, decoded) in enumerate(decoded_list)
        scores = decoded.average_score
        @inbounds @simd for t_idx = 1:n_timepoints
            diffs[p_idx, t_idx] = scores[t_idx] - chance_level
        end
    end

    return diffs, first_times, alpha_f
end

function _onesample_ttest_right!(t_statistics::Vector{Float64}, p_values::Vector{Float64}, diffs::Matrix{Float64})
    n_participants, n_timepoints = size(diffs)
    df = Float64(n_participants - 1)
    dist = TDist(df)

    @inbounds for t_idx = 1:n_timepoints
        sum_d = 0.0
        @simd for p_idx = 1:n_participants
            sum_d += diffs[p_idx, t_idx]
        end
        mean_d = sum_d / n_participants

        ss_d = 0.0
        @simd for p_idx = 1:n_participants
            delta = diffs[p_idx, t_idx] - mean_d
            ss_d += delta * delta
        end
        se = sqrt(max(ss_d / df, 0.0) / n_participants)
        t_val = se == 0.0 ? (mean_d == 0.0 ? NaN : copysign(Inf, mean_d)) : (mean_d / se)
        t_statistics[t_idx] = t_val
        p_values[t_idx] = isnan(t_val) ? NaN : ccdf(dist, t_val)
    end

    return df
end

function _apply_correction(p_values::Vector{Float64}, alpha::Float64, method::Symbol, t_statistics::Vector{Float64} = Float64[])
    method == :none && return .!isnan.(p_values) .& (p_values .<= alpha)
    method == :bonferroni && return _apply_bonferroni_correction(p_values, alpha, t_statistics)
    @minimal_error("Unknown correction method: $method. Must be :none or :bonferroni.")
end

"""
    test_against_chance(decoded_list::Vector{DecodedData}; alpha::Real = 0.05, correction_method::Symbol = :none)

Perform a right-tailed one-sample t-test against chance level across participants at each time point.

# Keyword Arguments
- `alpha::Real`: Significance threshold in `(0, 1)` (default: `0.05`)
- `correction_method::Symbol`: `:none` (default) or `:bonferroni`

# Returns
- `DecodingStatisticsResult`
"""
function test_against_chance(decoded_list::Vector{DecodedData}; alpha::Real = 0.05, correction_method::Symbol = :none)
    correction_method == :no && (correction_method = :none)
    correction_method ∈ (:none, :bonferroni) || @minimal_error("correction_method must be :none or :bonferroni, got :$correction_method")

    diffs, first_times, alpha_f = _validate_and_extract_diffs(decoded_list, "test_against_chance", alpha)
    n_timepoints = length(first_times)
    t_statistics = Vector{Float64}(undef, n_timepoints)
    p_values = Vector{Float64}(undef, n_timepoints)
    df = _onesample_ttest_right!(t_statistics, p_values, diffs)
    significant_mask = _apply_correction(p_values, alpha_f, correction_method, t_statistics)

    return DecodingStatisticsResult(first_times, t_statistics, p_values, significant_mask, df, alpha_f, correction_method, nothing)
end

function _find_temporal_clusters(mask::AbstractVector{Bool}, times::AbstractVector{Float64})
    clusters = TemporalCluster[]
    n = length(mask)
    t = 1
    id = 0
    @inbounds while t <= n
        if mask[t]
            s = t
            while t < n && mask[t+1]
                t += 1
            end
            id += 1
            push!(clusters, TemporalCluster(id, collect(s:t), (times[s], times[t]), 0.0, 1.0, false))
        end
        t += 1
    end
    return clusters
end

function _compute_cluster_statistics(clusters::Vector{TemporalCluster}, t_statistics::Vector{Float64}, statistic_type::Symbol = :sum)
    isempty(clusters) && return Float64[]
    cluster_stats = Vector{Float64}(undef, length(clusters))
    @inbounds for (i, cluster) in enumerate(clusters)
        if statistic_type == :sum
            s = 0.0
            for t_idx in cluster.time_indices
                s += t_statistics[t_idx]
            end
            cluster_stats[i] = s
        elseif statistic_type == :max
            m = -Inf
            for t_idx in cluster.time_indices
                a = abs(t_statistics[t_idx])
                a > m && (m = a)
            end
            cluster_stats[i] = m
        else
            @minimal_error("cluster_statistic must be :sum or :max, got :$statistic_type")
        end
    end
    return cluster_stats
end

"""
    test_against_chance_cluster(
        decoded_list::Vector{DecodedData};
        alpha::Real = 0.05,
        n_permutations::Int = 1000,
        cluster_statistic::Symbol = :sum,
        show_progress::Bool = true,
        rng::AbstractRNG = Random.default_rng(),
    )

Perform a right-tailed 1D temporal cluster-based sign-flip permutation test against chance level
across participants.

# Keyword Arguments
- `alpha::Real`: Cluster-forming and cluster-level significance threshold (default: `0.05`)
- `n_permutations::Int`: Number of random sign-flip permutations (default: `1000`)
- `cluster_statistic::Symbol`: `:sum` (default) or `:max`
- `show_progress::Bool`: Show progress bar (default: `true`)
- `rng::AbstractRNG`: Random number generator for reproducible sign flips

# Returns
- `DecodingStatisticsResult`
"""
function test_against_chance_cluster(
    decoded_list::Vector{DecodedData};
    alpha::Real = 0.05,
    n_permutations::Int = 1000,
    cluster_statistic::Symbol = :sum,
    show_progress::Bool = true,
    rng::AbstractRNG = Random.default_rng(),
)
    cluster_statistic ∈ (:sum, :max) || @minimal_error("cluster_statistic must be :sum or :max, got :$cluster_statistic")
    n_permutations < 1 && @minimal_error("n_permutations must be >= 1, got $n_permutations")

    diffs, first_times, alpha_f = _validate_and_extract_diffs(decoded_list, "test_against_chance_cluster", alpha)
    n_participants, n_timepoints = size(diffs)

    t_statistics = Vector{Float64}(undef, n_timepoints)
    p_values = Vector{Float64}(undef, n_timepoints)
    df = _onesample_ttest_right!(t_statistics, p_values, diffs)

    critical_t = quantile(TDist(df), 1.0 - alpha_f)
    mask_observed = .!isnan.(t_statistics) .& (t_statistics .> critical_t)
    observed_clusters = _find_temporal_clusters(mask_observed, first_times)
    cluster_stats_observed = _compute_cluster_statistics(observed_clusters, t_statistics, cluster_statistic)

    flips = rand(rng, Bool, n_participants, n_permutations)
    permutation_max = zeros(Float64, n_permutations)
    use_sum_stat = cluster_statistic == :sum
    completed = _start_permutation_progress(show_progress, n_permutations, "Permutations: ")

    _foreach_chunk(n_permutations) do perm_range
        @inbounds for perm_idx in perm_range
            best_stat = 0.0
            in_cluster = false
            cur_stat = 0.0

            for t_idx = 1:n_timepoints
                sum_d = 0.0
                @simd for p_idx = 1:n_participants
                    d = diffs[p_idx, t_idx]
                    sum_d += flips[p_idx, perm_idx] ? -d : d
                end
                mean_d = sum_d / n_participants

                ss_d = 0.0
                @simd for p_idx = 1:n_participants
                    d = diffs[p_idx, t_idx]
                    delta = (flips[p_idx, perm_idx] ? -d : d) - mean_d
                    ss_d += delta * delta
                end
                se = sqrt(max(ss_d / df, 0.0) / n_participants)
                t_val = se == 0.0 ? (mean_d == 0.0 ? NaN : copysign(Inf, mean_d)) : (mean_d / se)

                if !isnan(t_val) && t_val > critical_t
                    in_cluster = true
                    if use_sum_stat
                        cur_stat += t_val
                    else
                        abs_t = abs(t_val)
                        abs_t > cur_stat && (cur_stat = abs_t)
                    end
                elseif in_cluster
                    cur_stat > best_stat && (best_stat = cur_stat)
                    in_cluster = false
                    cur_stat = 0.0
                end
            end

            in_cluster && cur_stat > best_stat && (best_stat = cur_stat)
            permutation_max[perm_idx] = best_stat
            _tick!(completed)
        end
    end

    if !isempty(observed_clusters)
        inv_denom = 1.0 / (n_permutations + 1)
        updated_clusters = Vector{TemporalCluster}(undef, length(observed_clusters))
        for (i, cluster) in enumerate(observed_clusters)
            cluster_stat = cluster_stats_observed[i]
            p_value = (count(>=(cluster_stat), permutation_max) + 1) * inv_denom
            updated_clusters[i] =
                TemporalCluster(cluster.id, cluster.time_indices, cluster.time_range, cluster_stat, p_value, p_value <= alpha_f)
        end
        sort!(updated_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
        observed_clusters = updated_clusters
    end

    significant_mask = falses(n_timepoints)
    for cluster in observed_clusters
        cluster.is_significant && (significant_mask[cluster.time_indices] .= true)
    end

    return DecodingStatisticsResult(
        first_times,
        t_statistics,
        p_values,
        significant_mask,
        df,
        alpha_f,
        :cluster_permutation,
        observed_clusters,
    )
end
