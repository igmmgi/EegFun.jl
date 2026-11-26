# Cluster-level p-value inference and multiple-comparison significance masks.

######################################
######### Internal Functions #########
######################################

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
    clusters::Vector{C},
    cluster_stats::Vector{Float64},
    permutation_max::Vector{Float64},
    n_permutations::Int,
    alpha::Real = 0.05,
    tail::Symbol = :both,
) where {C<:Union{Cluster,TFCluster}}
    isempty(clusters) && return C[]
    updated = Vector{C}(undef, length(clusters))
    for (i, c) in enumerate(clusters)
        stat = cluster_stats[i]
        p_val, is_sig = _cluster_pval_and_sig(c.polarity, stat, permutation_max, n_permutations, alpha, tail)
        updated[i] = _rebuild_cluster(c, stat, p_val, is_sig)
    end
    return updated
end

function _significant_cluster_mask(clusters::Vector{<:Union{Cluster,TFCluster}}, dims::NTuple{N,Int}) where {N}
    mask = falses(dims)
    for c in clusters
        c.is_significant || continue
        @inbounds for pt in _cluster_points(c)
            mask[_to_ci(pt)] = true
        end
    end
    return mask
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
