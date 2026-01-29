# This file contains statistical inference logic for computing p-values
# and significance masks from statistical test results.
"""
    _compute_cluster_pvalues(clusters, cluster_stats, permutation_max, n_permutations, alpha)

Compute p-values for clusters by comparing to permutation distribution.

# Arguments
- `clusters::Vector{Cluster}`: Observed clusters
- `cluster_stats::Vector{Float64}`: Observed cluster statistics
- `permutation_max::Vector{Float64}`: Maximum cluster stats from permutations
- `n_permutations::Int`: Number of permutations
- `alpha::Float64`: Significance level (default: 0.05)

# Returns
- `updated_clusters::Vector{Cluster}`: Clusters with p-values and significance flags

# Examples
```julia
clusters = _compute_cluster_pvalues(clusters, stats, perm_max, 1000, 0.05)
```
"""
function _compute_cluster_pvalues(
    clusters::Vector{Cluster},
    cluster_stats::Vector{Float64},
    permutation_max::Vector{Float64},
    n_permutations::Int,
    alpha::Float64 = 0.05,
)
    if isempty(clusters)
        return Cluster[]
    end

    updated_clusters = Cluster[]

    for (i, cluster) in enumerate(clusters)
        cluster_stat = cluster_stats[i]

        # For positive clusters: count permutations with max >= observed
        # For negative clusters: count permutations with max <= observed (more extreme negative)
        # FieldTrip approach: compare observed to permutation distribution
        if cluster.polarity == :positive
            # Positive clusters: larger is more extreme
            count_exceed = sum(permutation_max .>= cluster_stat) + 1
        else
            # Negative clusters: more negative (smaller) is more extreme
            # permutation_max contains minimum (most negative) values from each permutation
            count_exceed = sum(permutation_max .<= cluster_stat) + 1
        end
        p_value = count_exceed / (n_permutations + 1)

        is_significant = p_value < alpha

        updated_cluster = Cluster(
            cluster.id,
            cluster.electrodes,
            cluster.time_indices,
            cluster.time_range,
            cluster_stat,
            p_value,
            is_significant,
            cluster.polarity,
        )
        push!(updated_clusters, updated_cluster)
    end

    return updated_clusters
end

"""
    _create_significance_mask(p_matrix, alpha, method)

Create significance mask from p-values.

# Arguments
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]
- `alpha::Float64`: Significance threshold
- `method::Symbol`: Correction method (`:no` or `:bonferroni`)

# Returns
- `mask::BitArray{2}`: Significance mask [electrodes × time]
"""
function _create_significance_mask(p_matrix::Array{Float64,2}, alpha::Float64, method::Symbol)
    method == :bonferroni && (alpha = alpha / count(!isnan, p_matrix))
    return .!isnan.(p_matrix) .& (p_matrix .<= alpha)
end

"""
    _apply_bonferroni_correction(p_values::Vector{Float64}, alpha::Float64, t_statistics::Vector{Float64})

Apply Bonferroni correction to p-values and return significance mask.

Shared utility function used by both general statistics and decoding statistics modules.
Handles NaN p-values that can occur with zero variance by treating them as significant
when the corresponding t-statistic is positive.

# Arguments
- `p_values::Vector{Float64}`: Uncorrected p-values
- `alpha::Float64`: Significance threshold
- `t_statistics::Vector{Float64}`: T-statistics corresponding to p-values

# Returns
- `BitVector`: Boolean mask indicating significant values after Bonferroni correction

# Examples
```julia
significant_mask = _apply_bonferroni_correction(p_values, 0.05, t_statistics)
```
"""
function _apply_bonferroni_correction(p_values::Vector{Float64}, alpha::Float64, t_statistics::Vector{Float64})
    n_comparisons = count(!isnan, p_values)
    corrected_alpha = n_comparisons > 0 ? alpha / n_comparisons : 0.0
    # Handle NaN p-values: if p is NaN but t > 0, treat as significant
    return (p_values .<= corrected_alpha) .| (isnan.(p_values) .& (t_statistics .> 0))
end

