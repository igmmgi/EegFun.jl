"""
    TemporalCluster

Stores information about a temporal cluster in decoding results.

# Fields
- `id::Int`: Unique cluster identifier
- `time_indices::Vector{Int}`: Time point indices in this cluster
- `time_range::Tuple{Float64, Float64}`: Time range in seconds (start, end)
- `cluster_stat::Float64`: Cluster-level statistic (e.g., sum of t-values)
- `p_value::Float64`: P-value from permutation test
- `is_significant::Bool`: Whether cluster is significant (p < alpha)
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

Stores results from statistical tests on decoding data.

# Fields
- `times::Vector{Float64}`: Time points in seconds
- `t_statistics::Vector{Float64}`: T-statistics at each time point
- `p_values::Vector{Float64}`: P-values at each time point
- `significant_mask::BitVector`: Boolean mask indicating significant time points
- `df::Float64`: Degrees of freedom
- `alpha::Float64`: Significance threshold used
- `correction_method::Symbol`: Multiple comparison correction method used
- `clusters::Union{Vector{TemporalCluster}, Nothing}`: Temporal clusters (if cluster-based test was performed)
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

# Centered two-pass one-sample right-tailed t-test against `chance_level`
function _onesample_ttest_right!(
    t_statistics::Vector{Float64},
    p_values::Vector{Float64},
    accuracies::Matrix{Float64},
    chance_level::Float64,
)
    n_participants, n_timepoints = size(accuracies)
    df = Float64(n_participants - 1)
    dist = TDist(df)

    @inbounds for t_idx = 1:n_timepoints
        sum_d = 0.0
        @simd for p_idx = 1:n_participants
            sum_d += accuracies[p_idx, t_idx] - chance_level
        end
        mean_d = sum_d / n_participants

        ss_d = 0.0
        @simd for p_idx = 1:n_participants
            delta = (accuracies[p_idx, t_idx] - chance_level) - mean_d
            ss_d += delta * delta
        end
        variance = max(ss_d / df, 0.0)
        se = sqrt(variance / n_participants)

        t_val = se == 0.0 ? (mean_d == 0.0 ? NaN : copysign(Inf, mean_d)) : (mean_d / se)
        t_statistics[t_idx] = t_val
        p_values[t_idx] = isnan(t_val) ? NaN : ccdf(dist, t_val)
    end

    return df
end

"""
    test_against_chance(decoded_list::Vector{DecodedData}; alpha::Real = 0.05, correction_method::Symbol = :none)

Perform one-sample t-test against chance level across multiple participants.

Tests whether decoding accuracy is significantly above chance at each time point,
using the across-participant variance. This is a one-tailed test (right tail) testing
if accuracy > chance level.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects, one per participant

# Keyword Arguments
- `alpha::Real`: Significance threshold (default: 0.05)
- `correction_method::Symbol`: Multiple comparison correction - `:none` (default) or `:bonferroni`

# Returns
- `DecodingStatisticsResult`: Results with t-statistics, p-values, and significance masks

# Examples
```julia
# Test across participants
decoded_list = [decoded_p1, decoded_p2, decoded_p3]
stats = EegFun.test_against_chance(decoded_list; alpha = 0.05, correction_method = :bonferroni)
```
"""
function test_against_chance(decoded_list::Vector{DecodedData}; alpha::Real = 0.05, correction_method::Symbol = :none)

    # validate some inputs
    if isempty(decoded_list)
        @minimal_error("Cannot test empty decoded data list")
    end
    length(decoded_list) < 2 && @minimal_error("test_against_chance requires at least 2 participants, got $(length(decoded_list))")
    alpha_f = Float64(alpha)
    (alpha_f <= 0.0 || alpha_f >= 1.0) && @minimal_error("alpha must be between 0 and 1, got $alpha_f")
    correction_method == :no && (correction_method = :none)
    if correction_method ∉ (:none, :bonferroni)
        @minimal_error("correction_method must be :none or :bonferroni, got :$correction_method")
    end

    # Validate all inputs have same structure
    first_decoded = decoded_list[1]
    first_times = first_decoded.times
    first_params = first_decoded.parameters
    chance_level = first_params.chance_level
    for decoded in decoded_list[2:end]
        if decoded.times != first_times
            @minimal_error("DecodedData inputs have inconsistent time vectors")
        end
        if decoded.parameters.chance_level != chance_level
            @minimal_warning "DecodedData objects have different chance_level"
        end
    end

    # Extract accuracies: [participants × timepoints]
    n_participants = length(decoded_list)
    n_timepoints = length(first_times)
    accuracies = Matrix{Float64}(undef, n_participants, n_timepoints)
    for (p_idx, decoded) in enumerate(decoded_list)
        accuracies[p_idx, :] = decoded.average_score
    end

    t_statistics = Vector{Float64}(undef, n_timepoints)
    p_values = Vector{Float64}(undef, n_timepoints)
    df = _onesample_ttest_right!(t_statistics, p_values, accuracies, chance_level)

    # Apply multiple comparison correction
    significant_mask = _apply_correction(p_values, alpha_f, correction_method, t_statistics)

    return DecodingStatisticsResult(
        first_times,
        t_statistics,
        p_values,
        significant_mask,
        df,
        alpha_f,
        correction_method,
        nothing,  # No clusters for simple t-test
    )
end

"""
    _apply_correction(p_values::Vector{Float64}, alpha::Float64, method::Symbol, t_statistics::Vector{Float64})

Apply multiple comparison correction to p-values.

# Arguments
- `p_values::Vector{Float64}`: Uncorrected p-values
- `alpha::Float64`: Significance threshold
- `method::Symbol`: Correction method - `:none` or `:bonferroni`
- `t_statistics::Vector{Float64}`: T-statistics corresponding to p-values

# Returns
- `BitVector`: Boolean mask indicating significant time points after correction
"""
function _apply_correction(p_values::Vector{Float64}, alpha::Float64, method::Symbol, t_statistics::Vector{Float64} = Float64[])
    if method == :none
        return .!isnan.(p_values) .& (p_values .<= alpha)
    elseif method == :bonferroni
        return _apply_bonferroni_correction(p_values, alpha, t_statistics)
    else
        @minimal_error("Unknown correction method: $method. Must be :none or :bonferroni.")
    end
end

# === CLUSTER-BASED PERMUTATION TESTING ===

"""
    _find_temporal_clusters(mask::BitVector, times::Vector{Float64})

Find temporal clusters in 1D significant mask.

Adapter function that converts 1D temporal mask to 2D format, calls the general
clustering function with cluster_type=:temporal, and converts results back to
TemporalCluster format.

# Arguments
- `mask::BitVector`: Boolean mask indicating significant time points (already thresholded)
- `times::Vector{Float64}`: Time points in seconds

# Returns
- `clusters::Vector{TemporalCluster}`: Found temporal clusters
"""
function _find_temporal_clusters(mask::BitVector, times::Vector{Float64})

    if !any(mask)
        return TemporalCluster[]
    end

    # Convert 1D mask to 2D (1 electrode × n_timepoints)
    mask_2d = BitArray{2}(undef, 1, length(mask))
    mask_2d[1, :] = mask

    # Create dummy electrode and connectivity (not used for temporal clustering)
    dummy_electrode = [:TimePoint]
    dummy_connectivity = sparse([1], [1], Bool[true], 1, 1)

    # Use general clustering function with :temporal type
    general_clusters = _find_clusters_connected_components(mask_2d, dummy_electrode, times, dummy_connectivity, :temporal)

    # Convert general Cluster objects to TemporalCluster objects
    temporal_clusters =
        [TemporalCluster(c.id, c.time_indices, c.time_range, c.cluster_stat, c.p_value, c.is_significant) for c in general_clusters]

    return temporal_clusters
end


"""
    _compute_cluster_statistics(clusters::Vector{TemporalCluster}, t_statistics::Vector{Float64}, statistic_type::Symbol = :sum)

Compute cluster-level statistics for temporal clusters.

Adapter function that uses the general cluster statistics computation.
For temporal-only clustering (1D), we use the :sum statistic (maxsum).

# Arguments
- `clusters::Vector{TemporalCluster}`: Temporal clusters
- `t_statistics::Vector{Float64}`: T-statistics at each time point
- `statistic_type::Symbol`: Statistic type - `:sum` (default) or `:max`

# Returns
- `cluster_stats::Vector{Float64}`: Cluster statistics
"""
function _compute_cluster_statistics(clusters::Vector{TemporalCluster}, t_statistics::Vector{Float64}, statistic_type::Symbol = :sum)

    if isempty(clusters)
        return Float64[]
    end

    cluster_stats = Vector{Float64}(undef, length(clusters))

    for (i, cluster) in enumerate(clusters)
        if statistic_type == :sum
            cluster_stats[i] = sum(t_statistics[t_idx] for t_idx in cluster.time_indices)
        elseif statistic_type == :max
            cluster_stats[i] = maximum(abs(t_statistics[t_idx]) for t_idx in cluster.time_indices)
        end
    end

    return cluster_stats
end


"""
    test_against_chance_cluster(decoded_list::Vector{DecodedData};
                               alpha::Float64 = 0.05,
                               n_permutations::Int = 1000,
                               cluster_statistic::Symbol = :sum,
                               show_progress::Bool = true,
                               rng::AbstractRNG = Random.default_rng())

Perform cluster-based permutation test against chance level for decoding results.

This method uses cluster-based permutation testing to control for multiple comparisons
across time points. Clusters of contiguous significant time points are identified, and
their cluster-level statistics are compared to a permutation distribution.
This is a one-tailed test (right tail) testing if accuracy > chance level.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects, one per participant

# Keyword Arguments
- `alpha::Float64`: Significance threshold (default: 0.05)
- `n_permutations::Int`: Number of permutations (default: 1000)
- `cluster_statistic::Symbol`: Cluster statistic - `:sum` (default) or `:max`
- `show_progress::Bool`: Show progress bar (default: true)
- `rng::AbstractRNG`: Random number generator used for the sign-flip permutations (default: `Random.default_rng()`)

# Notes
For reproducible results pass a seeded RNG, e.g. `rng = Random.Xoshiro(42)`.

# Returns
- `DecodingStatisticsResult`: Results with t-statistics, p-values, significance masks, and clusters

# Examples
```julia
# Test across participants with cluster-based correction
decoded_list = [decoded_p1, decoded_p2, decoded_p3]
stats = test_against_chance_cluster(decoded_list, alpha=0.05, n_permutations=1000)
```
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
    isempty(decoded_list) && @minimal_error("Cannot test empty decoded data list")
    length(decoded_list) < 2 && @minimal_error("test_against_chance_cluster requires at least 2 participants, got $(length(decoded_list))")
    n_permutations < 1 && @minimal_error("n_permutations must be >= 1, got $n_permutations")
    alpha_f = Float64(alpha)
    (alpha_f <= 0.0 || alpha_f >= 1.0) && @minimal_error("alpha must be between 0 and 1, got $alpha_f")

    # Validate all have same structure
    first_decoded = decoded_list[1]
    first_times = first_decoded.times
    first_params = first_decoded.parameters
    chance_level = first_params.chance_level

    for decoded in decoded_list[2:end]
        decoded.times != first_times && @minimal_error("DecodedData objects have inconsistent time vectors")
        decoded.parameters.chance_level != chance_level && @minimal_warning "DecodedData objects have different chance_level"
    end

    # Extract accuracies and centered differences: [participants × timepoints]
    n_participants = length(decoded_list)
    n_timepoints = length(first_times)
    accuracies = Matrix{Float64}(undef, n_participants, n_timepoints)
    diffs = Matrix{Float64}(undef, n_participants, n_timepoints)

    for (p_idx, decoded) in enumerate(decoded_list)
        @inbounds for t_idx = 1:n_timepoints
            acc = decoded.average_score[t_idx]
            accuracies[p_idx, t_idx] = acc
            diffs[p_idx, t_idx] = acc - chance_level
        end
    end

    # Compute observed t-statistics and p-values (right-tailed one-sample test)
    t_statistics = Vector{Float64}(undef, n_timepoints)
    p_values = Vector{Float64}(undef, n_timepoints)
    df = _onesample_ttest_right!(t_statistics, p_values, accuracies, chance_level)

    # Compute critical t-value for thresholding (one-tailed right)
    dist = TDist(df)
    critical_t = quantile(dist, 1.0 - alpha_f)

    # Threshold observed data (right tail: accuracy > chance)
    mask_observed = .!isnan.(t_statistics) .& (t_statistics .> critical_t)

    # Find observed clusters
    observed_clusters = _find_temporal_clusters(mask_observed, first_times)

    # Compute observed cluster statistics
    cluster_stats_observed = _compute_cluster_statistics(observed_clusters, t_statistics, cluster_statistic)

    # Pre-draw all sign-flip permutations up front for thread-count-independent reproducibility
    flips = rand(rng, Bool, n_participants, n_permutations)
    permutation_max = zeros(Float64, n_permutations)
    use_sum_stat = cluster_statistic == :sum
    completed = _start_permutation_progress(show_progress, n_permutations, "Permutations: ")

    # Fast zero-allocation threaded 1D contiguous cluster scan
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
                    val = flips[p_idx, perm_idx] ? -d : d
                    delta = val - mean_d
                    ss_d += delta * delta
                end
                variance = max(ss_d / df, 0.0)
                se = sqrt(variance / n_participants)
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

            if in_cluster && cur_stat > best_stat
                best_stat = cur_stat
            end
            permutation_max[perm_idx] = best_stat
            _tick!(completed)
        end
    end

    # Compute p-values for observed clusters
    if !isempty(observed_clusters)
        updated_clusters = TemporalCluster[]
        for (i, cluster) in enumerate(observed_clusters)
            cluster_stat = cluster_stats_observed[i]  # Don't use abs() for one-tailed test
            # Count permutations with max >= observed
            count_exceed = sum(permutation_max .>= cluster_stat) + 1
            p_value = count_exceed / (n_permutations + 1)
            is_significant = p_value <= alpha_f

            updated_cluster = TemporalCluster(cluster.id, cluster.time_indices, cluster.time_range, cluster_stat, p_value, is_significant)
            push!(updated_clusters, updated_cluster)
        end

        # Sort: significant first, then by cluster statistic (largest first)
        sort!(updated_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
        observed_clusters = updated_clusters
    end

    # Create significance mask: only points from significant clusters
    significant_mask = falses(n_timepoints)
    for cluster in observed_clusters
        if cluster.is_significant
            significant_mask[cluster.time_indices] .= true
        end
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
