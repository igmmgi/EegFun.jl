# Cluster-level p-value inference and multiple-comparison significance masks.

@inline function _cluster_pval_and_sig(
    polarity::Symbol,
    cluster_stat::Float64,
    permutation_max::Vector{Float64},
    n_permutations::Int,
    alpha::Real,
    tail::Symbol,
)
    count_exceed = polarity == :positive ? count(>=(cluster_stat), permutation_max) + 1 : count(<=(cluster_stat), permutation_max) + 1
    p_one_sided = count_exceed / (n_permutations + 1)
    p_value = tail == :both ? min(1.0, 2.0 * p_one_sided) : p_one_sided
    return p_value, p_value <= alpha
end

function _compute_cluster_pvalues(
    clusters::Vector{Cluster},
    cluster_stats::Vector{Float64},
    permutation_max::Vector{Float64},
    n_permutations::Int,
    alpha::Real = 0.05,
    tail::Symbol = :both,
)
    isempty(clusters) && return Cluster[]
    updated = Vector{Cluster}(undef, length(clusters))
    for (i, c) in enumerate(clusters)
        stat = cluster_stats[i]
        p_val, is_sig = _cluster_pval_and_sig(c.polarity, stat, permutation_max, n_permutations, alpha, tail)
        updated[i] = Cluster(c.id, c.electrodes, c.time_indices, c.time_range, stat, p_val, is_sig, c.polarity, c.members)
    end
    return updated
end

function _compute_cluster_pvalues(
    clusters::Vector{TFCluster},
    cluster_stats::Vector{Float64},
    permutation_max::Vector{Float64},
    n_permutations::Int,
    alpha::Real = 0.05,
    tail::Symbol = :both,
)
    isempty(clusters) && return TFCluster[]
    updated = Vector{TFCluster}(undef, length(clusters))
    for (i, c) in enumerate(clusters)
        stat = cluster_stats[i]
        p_val, is_sig = _cluster_pval_and_sig(c.polarity, stat, permutation_max, n_permutations, alpha, tail)
        updated[i] = TFCluster(
            c.id,
            c.electrodes,
            c.freq_indices,
            c.time_indices,
            c.freq_range,
            c.time_range,
            stat,
            p_val,
            is_sig,
            c.polarity,
            c.pixels,
        )
    end
    return updated
end

function _create_significance_mask(p_matrix::AbstractArray{Float64,N}, alpha::Float64, method::Symbol) where {N}
    thresh = if method == :bonferroni
        n_comp = count(!isnan, p_matrix)
        n_comp > 0 ? alpha / n_comp : 0.0
    else
        alpha
    end
    mask = BitArray{N}(undef, size(p_matrix))
    @inbounds for i in eachindex(p_matrix)
        p = p_matrix[i]
        mask[i] = !isnan(p) && p <= thresh
    end
    return mask
end

_apply_bonferroni_correction(p_values::Vector{Float64}, alpha::Float64, ::Vector{Float64} = Float64[]) =
    _create_significance_mask(p_values, alpha, :bonferroni)
