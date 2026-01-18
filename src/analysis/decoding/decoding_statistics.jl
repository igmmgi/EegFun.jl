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
    time_range::Tuple{Float64, Float64}
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
mutable struct DecodingStatisticsResult
    times::Vector{Float64}
    t_statistics::Vector{Float64}
    p_values::Vector{Float64}
    significant_mask::BitVector
    df::Float64
    alpha::Float64
    correction_method::Symbol
    clusters::Union{Vector{TemporalCluster}, Nothing}
end

"""
    test_against_chance(decoded_list::Vector{DecodedData}; alpha::Float64 = 0.05, correction_method::Symbol = :none)

Perform one-sample t-test against chance level across multiple participants.

Tests whether decoding accuracy is significantly above chance at each time point,
using the across-participant variance. This is a one-tailed test (right tail) testing
if accuracy > chance level.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects, one per participant

# Keyword Arguments
- `alpha::Float64`: Significance threshold (default: 0.05)
- `correction_method::Symbol`: Multiple comparison correction - `:none` (default) or `:bonferroni`

# Returns
- `DecodingStatisticsResult`: Results with t-statistics, p-values, and significance masks

# Examples
```julia
# Test across participants
decoded_list = [decoded_p1, decoded_p2, decoded_p3]
stats = test_against_chance(decoded_list, alpha=0.05, correction_method=:bonferroni)
```
"""
function test_against_chance(
    decoded_list::Vector{DecodedData};
    alpha::Float64 = 0.05,
    correction_method::Symbol = :none,
)
    correction_method ∈ (:none, :bonferroni) || @minimal_error_throw("correction_method must be :none or :bonferroni, got :$correction_method")
    
    isempty(decoded_list) && @minimal_error_throw("Cannot test empty decoded data list")
    
    # Validate all have same structure
    first_decoded = decoded_list[1]
    first_times = first_decoded.times
    first_params = first_decoded.parameters
    chance_level = first_params.chance_level
    
    for decoded in decoded_list[2:end]
        decoded.times != first_times && @minimal_error_throw("DecodedData objects have inconsistent time vectors")
        decoded.parameters.chance_level != chance_level && @minimal_warning "DecodedData objects have different chance_level"
    end
    
    # Extract accuracies: [participants × timepoints]
    n_participants = length(decoded_list)
    n_timepoints = length(first_times)
    accuracies = Matrix{Float64}(undef, n_participants, n_timepoints)
    
    for (p_idx, decoded) in enumerate(decoded_list)
        accuracies[p_idx, :] = decoded.average_score
    end
    
    # Compute t-statistics using AnovaFun's paired_ttest
    # For one-sample t-test against chance, we compare accuracies against a constant vector
    # of chance level values. This tests if mean(accuracies - chance) is significantly different from 0.
    t_statistics = Vector{Float64}(undef, n_timepoints)
    p_values = Vector{Float64}(undef, n_timepoints)
    chance_vector = fill(chance_level, n_participants)
    df = Float64(n_participants - 1)  # df for one-sample t-test
    
    for t_idx in 1:n_timepoints
        # Extract accuracies for this timepoint across all participants
        accuracies_at_t = accuracies[:, t_idx]
        
        # Use paired_ttest to compare accuracies against chance level
        # This is equivalent to a one-sample t-test: H0: mean(accuracies) = chance_level
        # One-tailed test (right tail): tests if accuracy > chance
        result = paired_ttest(accuracies_at_t, chance_vector, tail = :right)
        
        t_statistics[t_idx] = result.t
        p_values[t_idx] = result.p
    end
    
    # Apply multiple comparison correction
    significant_mask = _apply_correction(p_values, alpha, correction_method)
    
    return DecodingStatisticsResult(
        first_times,
        t_statistics,
        p_values,
        significant_mask,
        df,
        alpha,
        correction_method,
        nothing,  # No clusters for simple t-test
    )
end

"""
    _apply_correction(p_values::Vector{Float64}, alpha::Float64, method::Symbol)

Apply multiple comparison correction to p-values.

# Arguments
- `p_values::Vector{Float64}`: Uncorrected p-values
- `alpha::Float64`: Significance threshold
- `method::Symbol`: Correction method - `:none` or `:bonferroni`

# Returns
- `BitVector`: Boolean mask indicating significant time points after correction
"""
function _apply_correction(p_values::Vector{Float64}, alpha::Float64, method::Symbol)
    if method == :none
        return p_values .<= alpha
    elseif method == :bonferroni
        n_comparisons = count(!isnan, p_values)
        corrected_alpha = n_comparisons > 0 ? alpha / n_comparisons : 0.0
        return p_values .<= corrected_alpha
    else
        @minimal_error_throw("Unknown correction method: $method. Must be :none or :bonferroni.")
    end
end

# ===================
# CLUSTER-BASED PERMUTATION TESTING
# ===================

"""
    _find_temporal_clusters(mask::BitVector, times::Vector{Float64})

Find temporal clusters in 1D significant mask.

# Arguments
- `mask::BitVector`: Boolean mask indicating significant time points (already thresholded)
- `times::Vector{Float64}`: Time points in seconds

# Returns
- `clusters::Vector{TemporalCluster}`: Found temporal clusters
"""
function _find_temporal_clusters(
    mask::BitVector,
    times::Vector{Float64},
)
    clusters = TemporalCluster[]
    
    if !any(mask)
        return clusters
    end
    
    n_timepoints = length(mask)
    current_cluster_id = 0
    in_cluster = false
    cluster_start = 0
    
    for t_idx in 1:n_timepoints
        if mask[t_idx] && !in_cluster
            # Start new cluster
            in_cluster = true
            cluster_start = t_idx
            current_cluster_id += 1
        elseif !mask[t_idx] && in_cluster
            # End current cluster
            in_cluster = false
            cluster_end = t_idx - 1
            time_indices = collect(cluster_start:cluster_end)
            time_range = (times[cluster_start], times[cluster_end])
            
            cluster = TemporalCluster(
                current_cluster_id,
                time_indices,
                time_range,
                0.0,  # cluster_stat - will be computed later
                1.0,  # p_value - will be computed later
                false,  # is_significant - will be set later
            )
            push!(clusters, cluster)
        end
    end
    
    # Handle cluster that extends to end
    if in_cluster
        cluster_end = n_timepoints
        time_indices = collect(cluster_start:cluster_end)
        time_range = (times[cluster_start], times[cluster_end])
        
        cluster = TemporalCluster(
            current_cluster_id,
            time_indices,
            time_range,
            0.0,
            1.0,
            false,
        )
        push!(clusters, cluster)
    end
    
    return clusters
end

"""
    _compute_cluster_statistics(clusters::Vector{TemporalCluster}, t_statistics::Vector{Float64}, statistic_type::Symbol = :sum)

Compute cluster-level statistics.

# Arguments
- `clusters::Vector{TemporalCluster}`: Temporal clusters
- `t_statistics::Vector{Float64}`: T-statistics at each time point
- `statistic_type::Symbol`: Statistic type - `:sum` (default) or `:max`

# Returns
- `cluster_stats::Vector{Float64}`: Cluster statistics
"""
function _compute_cluster_statistics(
    clusters::Vector{TemporalCluster},
    t_statistics::Vector{Float64},
    statistic_type::Symbol = :sum,
)
    cluster_stats = Float64[]
    
    for cluster in clusters
        if statistic_type == :sum # Sum of t-values in cluster
            cluster_stat = sum(t_statistics[t_idx] for t_idx in cluster.time_indices)
        elseif statistic_type == :max # Maximum absolute t-value in cluster
            cluster_stat = maximum(abs(t_statistics[t_idx]) for t_idx in cluster.time_indices)
        else
            @minimal_error_throw("statistic_type must be :sum or :max, got :$statistic_type")
        end
        push!(cluster_stats, cluster_stat)
    end
    
    return cluster_stats
end

"""
    test_against_chance_cluster(decoded_list::Vector{DecodedData};
                               alpha::Float64 = 0.05,
                               n_permutations::Int = 1000,
                               cluster_statistic::Symbol = :sum,
                               random_seed::Union{Int, Nothing} = nothing,
                               show_progress::Bool = true)

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
- `random_seed::Union{Int, Nothing}`: Random seed for reproducibility (default: nothing)
- `show_progress::Bool`: Show progress bar (default: true)

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
    alpha::Float64 = 0.05,
    n_permutations::Int = 1000,
    cluster_statistic::Symbol = :sum,
    random_seed::Union{Int, Nothing} = nothing,
    show_progress::Bool = true,
)
    cluster_statistic ∈ (:sum, :max) || @minimal_error_throw("cluster_statistic must be :sum or :max, got :$cluster_statistic")
    
    isempty(decoded_list) && @minimal_error_throw("Cannot test empty decoded data list")
    
    # Validate all have same structure
    first_decoded = decoded_list[1]
    first_times = first_decoded.times
    first_params = first_decoded.parameters
    chance_level = first_params.chance_level
    
    for decoded in decoded_list[2:end]
        decoded.times != first_times && @minimal_error_throw("DecodedData objects have inconsistent time vectors")
        decoded.parameters.chance_level != chance_level && @minimal_warning "DecodedData objects have different chance_level"
    end
    
    # Extract accuracies: [participants × timepoints]
    n_participants = length(decoded_list)
    n_timepoints = length(first_times)
    accuracies = Matrix{Float64}(undef, n_participants, n_timepoints)
    
    for (p_idx, decoded) in enumerate(decoded_list)
        accuracies[p_idx, :] = decoded.average_score
    end
    
    # Set random seed if provided
    rng = random_seed !== nothing ? MersenneTwister(random_seed) : Random.GLOBAL_RNG
    
    # Compute observed t-statistics
    t_statistics = Vector{Float64}(undef, n_timepoints)
    p_values = Vector{Float64}(undef, n_timepoints)
    chance_vector = fill(chance_level, n_participants)
    df = Float64(n_participants - 1)
    
    for t_idx in 1:n_timepoints
        accuracies_at_t = accuracies[:, t_idx]
        # One-tailed test (right tail): tests if accuracy > chance
        result = paired_ttest(accuracies_at_t, chance_vector, tail = :right)
        t_statistics[t_idx] = result.t
        p_values[t_idx] = result.p
    end
    
    # Compute critical t-value for thresholding (one-tailed right)
    dist = TDist(df)
    critical_t = quantile(dist, 1.0 - alpha)
    
    # Threshold observed data (right tail: accuracy > chance)
    mask_observed = t_statistics .> critical_t
    
    # Find observed clusters
    observed_clusters = _find_temporal_clusters(mask_observed, first_times)
    
    # Compute observed cluster statistics
    if !isempty(observed_clusters)
        cluster_stats_observed = _compute_cluster_statistics(observed_clusters, t_statistics, cluster_statistic)
    else
        cluster_stats_observed = Float64[]
    end
    
    # Run permutations
    permutation_max = Float64[]
    sizehint!(permutation_max, n_permutations)
    
    if show_progress
        progress = Progress(n_permutations, desc = "Permutations: ", showspeed = true)
    end
    
    for perm_idx in 1:n_permutations
        # Shuffle participant labels (sign-flip for one-sample test)
        # For one-sample t-test, we randomly flip signs of accuracies
        shuffled_accuracies = copy(accuracies)
        for p_idx in 1:n_participants
            if rand(rng, Bool)
                # Flip sign: subtract from chance and negate
                shuffled_accuracies[p_idx, :] = 2 * chance_level .- shuffled_accuracies[p_idx, :]
            end
        end
        
        # Compute t-statistics for permuted data
        t_statistics_perm = Vector{Float64}(undef, n_timepoints)
        for t_idx in 1:n_timepoints
            accuracies_at_t = shuffled_accuracies[:, t_idx]
            # One-tailed test (right tail): tests if accuracy > chance
            result = paired_ttest(accuracies_at_t, chance_vector, tail = :right)
            t_statistics_perm[t_idx] = result.t
        end
        
        # Threshold permuted data (right tail: accuracy > chance)
        mask_perm = t_statistics_perm .> critical_t
        
        # Find clusters in permuted data
        clusters_perm = _find_temporal_clusters(mask_perm, first_times)
        
        # Compute max cluster statistic
        if !isempty(clusters_perm)
            cluster_stats_perm = _compute_cluster_statistics(clusters_perm, t_statistics_perm, cluster_statistic)
            max_stat = maximum(abs.(cluster_stats_perm))
        else
            max_stat = 0.0
        end
        
        push!(permutation_max, max_stat)
        
        if show_progress
            next!(progress)
        end
    end
    
    # Compute p-values for observed clusters
    if !isempty(observed_clusters)
        updated_clusters = TemporalCluster[]
        for (i, cluster) in enumerate(observed_clusters)
            cluster_stat = abs(cluster_stats_observed[i])
            # Count permutations with max >= observed
            count_exceed = sum(permutation_max .>= cluster_stat) + 1
            p_value = count_exceed / (n_permutations + 1)
            is_significant = p_value < alpha
            
            updated_cluster = TemporalCluster(
                cluster.id,
                cluster.time_indices,
                cluster.time_range,
                cluster_stat,
                p_value,
                is_significant,
            )
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
        alpha,
        :cluster_permutation,
        observed_clusters,
    )
end
