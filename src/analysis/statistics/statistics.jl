# Public statistical testing entry points: `permutation_test` and `analytic_test` (ERP and TF).

######################################
######### Types & Constants ##########
######################################

const _VALID_THRESHOLD_METHODS = (:parametric, :nonparametric_common, :nonparametric_individual)

const _compute_tf_cluster_pvalues = _compute_cluster_pvalues

######################################
########### User Functions ###########
######################################

"""
    permutation_test(prepared::StatisticalData; kwargs...)
    permutation_test(prepared::TFStatisticalData; kwargs...)

Perform a cluster-based permutation test (maxsum) on prepared ERP (`StatisticalData`) or
time-frequency (`TFStatisticalData`) data.

# Keyword Arguments
- `n_permutations::Int = 1000`: Number of Monte Carlo permutations
- `alpha::Real = 0.05`: Cluster-level family-wise error rate significance threshold
- `cluster_alpha::Real = 0.05`: Point-wise cluster-forming p-value threshold
- `threshold_method::Symbol = :parametric`: `:parametric`, `:nonparametric_common`, or `:nonparametric_individual`
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal` for ERP (default `:spatiotemporal`); additionally `:spectral`, `:spectrotemporal`, or `:full` for TF (default `:full`)
- `min_num_neighbors::Int = 2`: Minimum neighboring significant channels required to retain a point
- `tail::Symbol = :both`: `:both`, `:left`, or `:right`
- `show_progress::Bool = true`: Show progress bar
- `use_gpu::Bool = false`: (ERP only) Compute permutation t-matrices on an active GPU backend
- `rng::AbstractRNG = Random.default_rng()`: Random number generator (thread-count independent)
"""
function permutation_test(
    prepared::StatisticalData;
    n_permutations::Int = 1000,
    alpha::Real = 0.05,
    cluster_alpha::Real = 0.05,
    threshold_method::Symbol = :parametric,
    cluster_type::Symbol = :spatiotemporal,
    min_num_neighbors::Int = 2,
    tail::Symbol = :both,
    show_progress::Bool = true,
    use_gpu::Bool = false,
    rng::AbstractRNG = Random.default_rng(),
)
    alpha_f, cluster_alpha_f = Float64(alpha), Float64(cluster_alpha)
    _validate_perm_params(
        n_permutations,
        alpha_f,
        cluster_alpha_f,
        tail,
        threshold_method,
        cluster_type,
        (:spatial, :temporal, :spatiotemporal),
    )

    gpu_active = false
    if use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Computing $n_permutations permutations on $(gpu_device_name())..."
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true) but no functional GPU package is loaded. Falling back to CPU."
        end
    end

    plan = _draw_permutation_plan(rng, prepared, n_permutations)
    permutation_t_matrices = gpu_active ? _collect_permutation_t_matrices_gpu(prepared, plan) : nothing
    t_matrix, df, _, se_matrix = _compute_t_matrix(prepared; compute_p_values = false)

    critical_t_values, threshold_for_permutations, mask_positive, mask_negative, permutation_t_matrices = _resolve_perm_thresholds(
        prepared,
        t_matrix,
        df,
        plan,
        threshold_method,
        cluster_alpha_f,
        tail,
        show_progress;
        permutation_t_matrices = permutation_t_matrices,
    )

    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points
    spatial_connectivity =
        _prepare_clusters_masks!(mask_positive, mask_negative, electrodes, prepared.data[1].layout, cluster_type, min_num_neighbors)

    pos_clusters, neg_clusters = _find_clusters(mask_positive, mask_negative, electrodes, time_points, spatial_connectivity, cluster_type)
    pos_clusters, pos_stats = _score_clusters(pos_clusters, t_matrix)
    neg_clusters, neg_stats = _score_clusters(neg_clusters, t_matrix)

    perm_max_pos, perm_max_neg = _run_permutations(
        prepared,
        plan,
        threshold_for_permutations,
        spatial_connectivity,
        cluster_type,
        tail,
        min_num_neighbors,
        show_progress;
        permutation_t_matrices = permutation_t_matrices,
    )

    pos_clusters, sig_mask_pos =
        _finalize_polarity_clusters(pos_clusters, pos_stats, perm_max_pos, n_permutations, alpha_f, tail, size(mask_positive))
    neg_clusters, sig_mask_neg =
        _finalize_polarity_clusters(neg_clusters, neg_stats, perm_max_neg, n_permutations, alpha_f, tail, size(mask_negative))

    result = PermutationResult(
        TestInfo(
            prepared.analysis.design,
            df,
            alpha_f,
            tail,
            :cluster_permutation,
            ClusterInfo(threshold_method, cluster_type, n_permutations, cluster_alpha_f),
        ),
        prepared.data,
        StatMatrix(t_matrix, nothing),
        Masks(sig_mask_pos, sig_mask_neg),
        Clusters(pos_clusters, neg_clusters),
        PermutationDistribution(perm_max_pos, perm_max_neg),
        electrodes,
        time_points,
        critical_t_values,
        se_matrix,
        prepared.se_cond1,
        prepared.se_cond2,
        prepared.se_diff,
    )

    @info "Permutation test complete. Found $(length(pos_clusters)) positive and $(length(neg_clusters)) negative clusters."
    return result
end

function permutation_test(
    prepared::TFStatisticalData;
    n_permutations::Int = 1000,
    alpha::Real = 0.05,
    cluster_alpha::Real = 0.05,
    threshold_method::Symbol = :parametric,
    cluster_type::Symbol = :full,
    min_num_neighbors::Int = 2,
    tail::Symbol = :both,
    show_progress::Bool = true,
    rng::AbstractRNG = Random.default_rng(),
)
    alpha_f, cluster_alpha_f = Float64(alpha), Float64(cluster_alpha)
    _validate_perm_params(
        n_permutations,
        alpha_f,
        cluster_alpha_f,
        tail,
        threshold_method,
        cluster_type,
        (:temporal, :spatial, :spectral, :spatiotemporal, :spectrotemporal, :full),
    )

    plan = _draw_permutation_plan(rng, prepared, n_permutations)
    t_matrix, df, _ = _compute_t_matrix_tf(prepared; tail = tail, compute_p_values = false)

    electrodes = channel_labels(prepared.data[1])
    frequencies = prepared.analysis.frequencies
    time_points = prepared.analysis.time_points

    _, threshold_for_permutations, mask_positive, mask_negative, permutation_t_matrices =
        _resolve_perm_thresholds(prepared, t_matrix, df, plan, threshold_method, cluster_alpha_f, tail, show_progress; stream_large = true)

    spatial_connectivity =
        _prepare_clusters_masks!(mask_positive, mask_negative, electrodes, prepared.data[1].layout, cluster_type, min_num_neighbors)

    pos_clusters, neg_clusters =
        _find_clusters_tf(mask_positive, mask_negative, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)
    pos_clusters, pos_stats = _score_clusters(pos_clusters, t_matrix)
    neg_clusters, neg_stats = _score_clusters(neg_clusters, t_matrix)

    perm_max_pos, perm_max_neg = _run_permutations_tf(
        prepared,
        plan,
        threshold_for_permutations,
        spatial_connectivity,
        cluster_type,
        tail,
        min_num_neighbors,
        electrodes,
        frequencies,
        time_points,
        show_progress;
        permutation_t_matrices = permutation_t_matrices,
    )

    pos_clusters, sig_mask_pos =
        _finalize_polarity_clusters(pos_clusters, pos_stats, perm_max_pos, n_permutations, alpha_f, tail, size(mask_positive))
    neg_clusters, sig_mask_neg =
        _finalize_polarity_clusters(neg_clusters, neg_stats, perm_max_neg, n_permutations, alpha_f, tail, size(mask_negative))

    result = TFClusterPermutationResult(
        TestInfo(
            prepared.analysis.design,
            df,
            alpha_f,
            tail,
            :cluster_permutation,
            ClusterInfo(threshold_method, cluster_type, n_permutations, cluster_alpha_f),
        ),
        prepared.data,
        TFStatMatrix(t_matrix, nothing),
        TFMasks(sig_mask_pos, sig_mask_neg),
        TFClusters(pos_clusters, neg_clusters),
        PermutationDistribution(perm_max_pos, perm_max_neg),
        electrodes,
        frequencies,
        time_points,
    )

    @info "TF Permutation test complete. Found $(length(pos_clusters)) positive and $(length(neg_clusters)) negative clusters."
    return result
end

"""
    analytic_test(prepared::StatisticalData; alpha = 0.05, tail = :both, correction_method = :none)
    analytic_test(prepared::TFStatisticalData; alpha = 0.05, tail = :both, correction_method = :none)

Perform an analytic (parametric) mass-univariate t-test on prepared ERP (`StatisticalData`)
or time-frequency (`TFStatisticalData`) data.

# Keyword Arguments
- `alpha::Real = 0.05`: Significance threshold
- `tail::Symbol = :both`: Test tail (`:both`, `:left`, or `:right`)
- `correction_method::Symbol = :none`: Multiple comparison correction (`:none` / `:no` or `:bonferroni`)
"""
function analytic_test(prepared::StatisticalData; alpha::Real = 0.05, tail::Symbol = :both, correction_method::Symbol = :none)
    alpha_f, correction_method = _validate_analytic_params(alpha, tail, correction_method)
    t_matrix, df, p_matrix, se_matrix = _compute_t_matrix(prepared, tail = tail)
    mask_positive, mask_negative, critical_t = _validate_and_mask_analytic(t_matrix, p_matrix, df, alpha_f, tail, correction_method)

    result = AnalyticResult(
        TestInfo(prepared.analysis.design, df, alpha_f, tail, correction_method, nothing),
        prepared.data,
        StatMatrix(t_matrix, p_matrix),
        Masks(mask_positive, mask_negative),
        channel_labels(prepared.data[1]),
        prepared.analysis.time_points,
        critical_t,
        se_matrix,
        prepared.se_cond1,
        prepared.se_cond2,
        prepared.se_diff,
    )

    @info "Analytic test complete: $(count(mask_positive)) +ve, $(count(mask_negative)) -ve significant points."
    return result
end

function analytic_test(prepared::TFStatisticalData; alpha::Real = 0.05, tail::Symbol = :both, correction_method::Symbol = :none)
    alpha_f, correction_method = _validate_analytic_params(alpha, tail, correction_method)
    t_matrix, df, p_matrix = _compute_t_matrix_tf(prepared, tail = tail)
    mask_positive, mask_negative, critical_t = _validate_and_mask_analytic(t_matrix, p_matrix, df, alpha_f, tail, correction_method)

    result = TFAnalyticResult(
        TestInfo(prepared.analysis.design, df, alpha_f, tail, correction_method, nothing),
        prepared.data,
        TFStatMatrix(t_matrix, p_matrix),
        TFMasks(mask_positive, mask_negative),
        channel_labels(prepared.data[1]),
        prepared.analysis.frequencies,
        prepared.analysis.time_points,
        critical_t,
    )

    @info "TF Analytic test complete: $(count(mask_positive)) +ve, $(count(mask_negative)) -ve significant points."
    return result
end

######################################
######### Internal Functions #########
######################################

function _validate_perm_params(
    n_permutations::Int,
    alpha_f::Float64,
    cluster_alpha_f::Float64,
    tail::Symbol,
    threshold_method::Symbol,
    cluster_type::Symbol,
    valid_cluster_types::Tuple,
)
    n_permutations < 1 && @minimal_error "n_permutations must be >= 1, got $n_permutations"
    (0.0 < alpha_f < 1.0) || @minimal_error "alpha must be between 0 and 1, got $alpha_f"
    (0.0 < cluster_alpha_f < 1.0) || @minimal_error "cluster_alpha must be between 0 and 1, got $cluster_alpha_f"
    cluster_type ∈ valid_cluster_types || @minimal_error "cluster_type must be one of $valid_cluster_types, got :$cluster_type"
    tail ∈ (:both, :left, :right) || @minimal_error "tail must be :both, :left, or :right, got :$tail"
    threshold_method ∈ _VALID_THRESHOLD_METHODS ||
        @minimal_error "threshold_method must be :parametric, :nonparametric_common, or :nonparametric_individual, got :$threshold_method"
    return nothing
end

function _prepare_clusters_masks!(
    mask_positive::BitArray{N},
    mask_negative::BitArray{N},
    electrodes::Vector{Symbol},
    layout::Layout,
    cluster_type::Symbol,
    min_num_neighbors::Int,
) where {N}
    spatial_connectivity = _build_connectivity_matrix(electrodes, layout, cluster_type)
    if min_num_neighbors > 0 && cluster_type ∈ (:spatial, :spatiotemporal, :full)
        _prefilter_mask_by_neighbors!(mask_positive, spatial_connectivity, min_num_neighbors)
        _prefilter_mask_by_neighbors!(mask_negative, spatial_connectivity, min_num_neighbors)
    end
    return spatial_connectivity
end

_score_clusters(clusters::Vector{C}, t_matrix::AbstractArray{<:AbstractFloat}) where {C<:Union{Cluster,TFCluster}} =
    isempty(clusters) ? (clusters, Float64[]) : _compute_cluster_statistics(clusters, t_matrix)

function _finalize_polarity_clusters(
    clusters::Vector{C},
    stats::Vector{Float64},
    perm_max::Vector{Float64},
    n_permutations::Int,
    alpha_f::Float64,
    tail::Symbol,
    dims::NTuple{N,Int},
) where {C<:Union{Cluster,TFCluster},N}
    if !isempty(clusters)
        clusters = _compute_cluster_pvalues(clusters, stats, perm_max, n_permutations, alpha_f, tail)
        sort!(clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end
    return clusters, _significant_cluster_mask(clusters, dims)
end

function _resolve_perm_thresholds(
    prepared::Union{StatisticalData,TFStatisticalData},
    t_matrix::Array{Float64,N},
    df::Float64,
    plan::AbstractMatrix,
    threshold_method::Symbol,
    cluster_alpha_f::Float64,
    tail::Symbol,
    show_progress::Bool;
    permutation_t_matrices::Union{Nothing,Array{Float64}} = nothing,
    stream_large::Bool = false,
) where {N}
    if threshold_method == :parametric
        crit = _compute_critical_t_values(df, size(t_matrix), cluster_alpha_f, tail)
        mask_pos, mask_neg = _threshold_t_matrix_parametric(t_matrix, crit, tail)
        return crit, crit, mask_pos, mask_neg, permutation_t_matrices
    end

    use_streaming = stream_large && length(t_matrix) * size(plan, 2) > _TF_PERM_CACHE_MAX_FLOATS
    if !use_streaming && isnothing(permutation_t_matrices)
        permutation_t_matrices = _collect_perm_t_generic(
            prepared.analysis.data[1],
            prepared.analysis.data[2],
            prepared.analysis.design,
            plan,
            show_progress,
            prepared isa TFStatisticalData ? "Collecting TF t-matrices: " : "Collecting permutation t-matrices: ",
        )
    end

    if threshold_method == :nonparametric_common
        thresh_pos, thresh_neg = if use_streaming
            _stream_nonparametric_threshold_common_tf(prepared::TFStatisticalData, plan, cluster_alpha_f, tail, show_progress)
        else
            _compute_nonparametric_threshold_common(permutation_t_matrices, cluster_alpha_f, tail)
        end
        mask_pos, mask_neg = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos, thresh_neg, tail)
        return (thresh_pos, -thresh_neg), (thresh_pos, thresh_neg), mask_pos, mask_neg, permutation_t_matrices
    else
        thresh_pos_mat, thresh_neg_mat = if use_streaming
            _stream_nonparametric_threshold_individual_tf(prepared::TFStatisticalData, plan, cluster_alpha_f, tail, show_progress)
        else
            _compute_nonparametric_threshold_individual(permutation_t_matrices, cluster_alpha_f, tail)
        end
        mask_pos, mask_neg = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos_mat, thresh_neg_mat, tail)
        return (thresh_pos_mat, -thresh_neg_mat), (thresh_pos_mat, thresh_neg_mat), mask_pos, mask_neg, permutation_t_matrices
    end
end

function _validate_and_mask_analytic(
    t_matrix::Array{Float64,N},
    p_matrix::Array{Float64,N},
    df::Float64,
    alpha_f::Float64,
    tail::Symbol,
    correction_method::Symbol,
) where {N}
    (!isfinite(df) || df <= 0) && @minimal_error "Invalid degrees of freedom: df=$df!"

    thresh = if correction_method == :bonferroni
        n_comp = count(!isnan, p_matrix)
        n_sig_raw = count(p -> !isnan(p) && p <= alpha_f, p_matrix)
        bonf_alpha = n_comp > 0 ? alpha_f / n_comp : 0.0
        n_sig_bonf = count(p -> !isnan(p) && p <= bonf_alpha, p_matrix)
        @info "Bonferroni correction: n_comparisons=$n_comp, bonferroni_alpha=$bonf_alpha, raw_sig=$n_sig_raw, bonferroni_sig=$n_sig_bonf"
        bonf_alpha
    else
        alpha_f
    end

    mask_positive = BitArray{N}(undef, size(t_matrix))
    mask_negative = BitArray{N}(undef, size(t_matrix))
    if tail == :both
        @inbounds for i in eachindex(t_matrix)
            p, t = p_matrix[i], t_matrix[i]
            sig = !isnan(p) && p <= thresh && !isnan(t)
            mask_positive[i] = sig && t > 0.0
            mask_negative[i] = sig && t < 0.0
        end
    elseif tail == :right
        fill!(mask_negative, false)
        @inbounds for i in eachindex(t_matrix)
            p = p_matrix[i]
            mask_positive[i] = !isnan(p) && p <= thresh && !isnan(t_matrix[i])
        end
    else
        fill!(mask_positive, false)
        @inbounds for i in eachindex(t_matrix)
            p = p_matrix[i]
            mask_negative[i] = !isnan(p) && p <= thresh && !isnan(t_matrix[i])
        end
    end

    return mask_positive, mask_negative, _critical_t_scalar(df, alpha_f, tail)
end

function _validate_analytic_params(alpha::Real, tail::Symbol, correction_method::Symbol)
    correction_method == :no && (correction_method = :none)
    correction_method ∈ (:none, :bonferroni) || @minimal_error "correction_method must be :none or :bonferroni. Got :$correction_method"
    tail ∈ (:both, :left, :right) || @minimal_error "tail must be :both, :left, or :right. Got :$tail"
    alpha_f = Float64(alpha)
    (0.0 < alpha_f < 1.0) || @minimal_error "alpha must be between 0 and 1, got $alpha_f"
    return alpha_f, correction_method
end
