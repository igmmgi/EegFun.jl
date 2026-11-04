# Public statistical testing entry points: `permutation_test` and `analytic_test` (ERP and TF).

const _VALID_THRESHOLD_METHODS = (:parametric, :nonparametric_common, :nonparametric_individual)

function _resolve_erp_thresholds(
    prepared::StatisticalData,
    t_matrix::Matrix{Float64},
    df::Float64,
    plan::AbstractMatrix,
    threshold_method::Symbol,
    cluster_alpha_f::Float64,
    tail::Symbol,
    show_progress::Bool,
    permutation_t_matrices::Union{Nothing,Array{Float64,3}},
)
    if threshold_method == :parametric
        @info "Computing parametric critical t-values..."
        crit = _compute_critical_t_values(df, size(t_matrix), cluster_alpha_f, tail)
        @info "Thresholding observed data (parametric)..."
        mask_pos, mask_neg = _threshold_t_matrix_parametric(t_matrix, crit, tail)
        return crit, crit, mask_pos, mask_neg, permutation_t_matrices
    end

    if isnothing(permutation_t_matrices)
        @info "Collecting permutation t-matrices for $(threshold_method) thresholding..."
        permutation_t_matrices = _collect_permutation_t_matrices(prepared, plan, show_progress)
    end

    if threshold_method == :nonparametric_common
        @info "Computing non-parametric common threshold..."
        thresh_pos, thresh_neg = _compute_nonparametric_threshold_common(permutation_t_matrices, cluster_alpha_f, tail)
        @info "Thresholding observed data (non-parametric common)..."
        mask_pos, mask_neg = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos, thresh_neg, tail)
        return (thresh_pos, -thresh_neg), (thresh_pos, thresh_neg), mask_pos, mask_neg, permutation_t_matrices
    else
        @info "Computing non-parametric individual thresholds..."
        thresh_pos_mat, thresh_neg_mat = _compute_nonparametric_threshold_individual(permutation_t_matrices, cluster_alpha_f, tail)
        @info "Thresholding observed data (non-parametric individual)..."
        mask_pos, mask_neg = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos_mat, thresh_neg_mat, tail)
        return (thresh_pos_mat, -thresh_neg_mat), (thresh_pos_mat, thresh_neg_mat), mask_pos, mask_neg, permutation_t_matrices
    end
end

function _resolve_tf_thresholds(
    prepared::TFStatisticalData,
    t_matrix::Array{Float64,3},
    df::Float64,
    plan::AbstractMatrix,
    threshold_method::Symbol,
    cluster_alpha_f::Float64,
    tail::Symbol,
    show_progress::Bool,
)
    if threshold_method == :parametric
        @info "Computing parametric critical t-values..."
        crit = _compute_critical_t_values_tf(df, size(t_matrix), cluster_alpha_f, tail)
        @info "Thresholding observed data (parametric)..."
        mask_pos, mask_neg = _threshold_t_matrix_parametric_tf(t_matrix, crit, tail)
        return crit, mask_pos, mask_neg, nothing
    end

    total_floats = length(t_matrix) * size(plan, 2)
    permutation_t_matrices = nothing

    if threshold_method == :nonparametric_common
        if total_floats <= _TF_PERM_CACHE_MAX_FLOATS
            @info "Collecting TF permutation t-matrices in parallel..."
            permutation_t_matrices = _collect_permutation_t_matrices_tf(
                prepared.analysis.data[1],
                prepared.analysis.data[2],
                prepared.analysis.design,
                plan,
                show_progress,
            )
            @info "Computing non-parametric common threshold..."
            thresh_pos, thresh_neg = _compute_nonparametric_threshold_common_tf(permutation_t_matrices, cluster_alpha_f, tail)
        else
            @info "Streaming non-parametric common threshold across chunks (avoiding large 4D allocation)..."
            thresh_pos, thresh_neg = _stream_nonparametric_threshold_common_tf(prepared, plan, cluster_alpha_f, tail, show_progress)
        end
        @info "Thresholding observed data (non-parametric common)..."
        mask_pos, mask_neg = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos, thresh_neg, tail)
        return (thresh_pos, thresh_neg), mask_pos, mask_neg, permutation_t_matrices
    else
        if total_floats <= _TF_PERM_CACHE_MAX_FLOATS
            @info "Collecting TF permutation t-matrices in parallel..."
            permutation_t_matrices = _collect_permutation_t_matrices_tf(
                prepared.analysis.data[1],
                prepared.analysis.data[2],
                prepared.analysis.design,
                plan,
                show_progress,
            )
            @info "Computing non-parametric individual thresholds..."
            thresh_pos_mat, thresh_neg_mat = _compute_nonparametric_threshold_individual_tf(permutation_t_matrices, cluster_alpha_f, tail)
        else
            @info "Streaming non-parametric individual thresholds per electrode (avoiding large 4D allocation)..."
            thresh_pos_mat, thresh_neg_mat =
                _stream_nonparametric_threshold_individual_tf(prepared, plan, cluster_alpha_f, tail, show_progress)
        end
        @info "Thresholding observed data (non-parametric individual)..."
        mask_pos, mask_neg = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos_mat, thresh_neg_mat, tail)
        return (thresh_pos_mat, thresh_neg_mat), mask_pos, mask_neg, permutation_t_matrices
    end
end

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
    alpha_f = Float64(alpha)
    cluster_alpha_f = Float64(cluster_alpha)

    _validate_permutation_inputs(prepared, n_permutations, alpha_f, cluster_type, tail; cluster_alpha = cluster_alpha_f)
    threshold_method ∈ _VALID_THRESHOLD_METHODS ||
        @minimal_error("threshold_method must be :parametric, :nonparametric_common, or :nonparametric_individual. Got :$threshold_method")

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

    @info "Computing t-statistics..."
    t_matrix, df, _, se_matrix = _compute_t_matrix(prepared; compute_p_values = false)

    critical_t_values, threshold_for_permutations, mask_positive, mask_negative, permutation_t_matrices = _resolve_erp_thresholds(
        prepared,
        t_matrix,
        df,
        plan,
        threshold_method,
        cluster_alpha_f,
        tail,
        show_progress,
        permutation_t_matrices,
    )

    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points
    layout = prepared.data[1].layout
    if cluster_type ∈ (:spatial, :spatiotemporal) && isnothing(layout.neighbors)
        layout = copy(layout)
        get_neighbors_xy!(layout, 0.35)
    end

    @info "Building connectivity matrix..."
    spatial_connectivity = _build_connectivity_matrix(electrodes, layout, cluster_type)

    if min_num_neighbors > 0
        @info "Pre-filtering masks (min_num_neighbors=$min_num_neighbors)..."
        _prefilter_mask_by_neighbors!(mask_positive, spatial_connectivity, min_num_neighbors)
        _prefilter_mask_by_neighbors!(mask_negative, spatial_connectivity, min_num_neighbors)
    end

    @info "Finding observed clusters..."
    positive_clusters, negative_clusters =
        _find_clusters(mask_positive, mask_negative, electrodes, time_points, spatial_connectivity, cluster_type)

    cluster_stats_positive = Float64[]
    cluster_stats_negative = Float64[]
    if !isempty(positive_clusters)
        positive_clusters, cluster_stats_positive = _compute_cluster_statistics(positive_clusters, t_matrix)
    end
    if !isempty(negative_clusters)
        negative_clusters, cluster_stats_negative = _compute_cluster_statistics(negative_clusters, t_matrix)
    end

    @info "Running $n_permutations permutations for cluster-level inference..."
    permutation_max_positive, permutation_max_negative = _run_permutations(
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

    @info "Computing p-values..."
    if !isempty(positive_clusters)
        positive_clusters =
            _compute_cluster_pvalues(positive_clusters, cluster_stats_positive, permutation_max_positive, n_permutations, alpha_f, tail)
        sort!(positive_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end
    if !isempty(negative_clusters)
        negative_clusters =
            _compute_cluster_pvalues(negative_clusters, cluster_stats_negative, permutation_max_negative, n_permutations, alpha_f, tail)
        sort!(negative_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end

    significant_mask_positive = falses(size(mask_positive))
    significant_mask_negative = falses(size(mask_negative))
    for cluster in positive_clusters
        cluster.is_significant || continue
        @inbounds for (e_idx, t_idx) in cluster.members
            significant_mask_positive[e_idx, t_idx] = true
        end
    end
    for cluster in negative_clusters
        cluster.is_significant || continue
        @inbounds for (e_idx, t_idx) in cluster.members
            significant_mask_negative[e_idx, t_idx] = true
        end
    end

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
        Masks(significant_mask_positive, significant_mask_negative),
        Clusters(positive_clusters, negative_clusters),
        PermutationDistribution(permutation_max_positive, permutation_max_negative),
        electrodes,
        time_points,
        critical_t_values,
        se_matrix,
        prepared.se_cond1,
        prepared.se_cond2,
        prepared.se_diff,
    )

    @info "Permutation test complete. Found $(length(positive_clusters)) positive and $(length(negative_clusters)) negative clusters."
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
    alpha_f = Float64(alpha)
    cluster_alpha_f = Float64(cluster_alpha)

    n_permutations < 1 && @minimal_error "n_permutations must be >= 1, got $n_permutations"
    (alpha_f <= 0.0 || alpha_f >= 1.0) && @minimal_error "alpha must be between 0 and 1, got $alpha_f"
    (cluster_alpha_f <= 0.0 || cluster_alpha_f >= 1.0) && @minimal_error "cluster_alpha must be between 0 and 1, got $cluster_alpha_f"
    cluster_type ∉ (:temporal, :spatial, :spectral, :spatiotemporal, :spectrotemporal, :full) &&
        @minimal_error "cluster_type must be :temporal, :spatial, :spectral, :spatiotemporal, :spectrotemporal, or :full, got :$cluster_type"
    tail ∉ (:both, :left, :right) && @minimal_error "tail must be :both, :left, or :right, got :$tail"
    threshold_method ∈ _VALID_THRESHOLD_METHODS ||
        @minimal_error "threshold_method must be :parametric, :nonparametric_common, or :nonparametric_individual, got :$threshold_method"

    plan = _draw_permutation_plan(rng, prepared, n_permutations)

    @info "Computing TF t-statistics..."
    t_matrix, df, _ = _compute_t_matrix_tf(prepared; tail = tail, compute_p_values = false)

    electrodes = channel_labels(prepared.data[1])
    frequencies = prepared.analysis.frequencies
    time_points = prepared.analysis.time_points

    threshold_for_permutations, mask_positive, mask_negative, permutation_t_matrices =
        _resolve_tf_thresholds(prepared, t_matrix, df, plan, threshold_method, cluster_alpha_f, tail, show_progress)

    layout = prepared.data[1].layout
    if cluster_type ∈ (:spatial, :spatiotemporal, :full) && isnothing(layout.neighbors)
        @minimal_warning "Layout.neighbors is not set. Computing with default distance criterion (0.35)."
        layout = copy(layout)
        get_neighbors_xy!(layout, 0.35)
    end

    @info "Building connectivity matrix..."
    spatial_connectivity = _build_connectivity_matrix(electrodes, layout, cluster_type)

    if min_num_neighbors > 0
        @info "Pre-filtering masks (min_num_neighbors=$min_num_neighbors)..."
        _prefilter_mask_by_neighbors_tf!(mask_positive, spatial_connectivity, min_num_neighbors)
        _prefilter_mask_by_neighbors_tf!(mask_negative, spatial_connectivity, min_num_neighbors)
    end

    @info "Finding observed clusters..."
    positive_clusters, negative_clusters =
        _find_clusters_tf(mask_positive, mask_negative, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)

    cluster_stats_positive = Float64[]
    cluster_stats_negative = Float64[]
    if !isempty(positive_clusters)
        positive_clusters, cluster_stats_positive = _compute_cluster_statistics_tf(positive_clusters, t_matrix, electrodes)
    end
    if !isempty(negative_clusters)
        negative_clusters, cluster_stats_negative = _compute_cluster_statistics_tf(negative_clusters, t_matrix, electrodes)
    end

    @info "Running $n_permutations permutations for cluster-level inference..."
    permutation_max_positive, permutation_max_negative = _run_permutations_tf(
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

    @info "Computing p-values..."
    if !isempty(positive_clusters)
        positive_clusters =
            _compute_tf_cluster_pvalues(positive_clusters, cluster_stats_positive, permutation_max_positive, n_permutations, alpha_f, tail)
        sort!(positive_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end
    if !isempty(negative_clusters)
        negative_clusters =
            _compute_tf_cluster_pvalues(negative_clusters, cluster_stats_negative, permutation_max_negative, n_permutations, alpha_f, tail)
        sort!(negative_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end

    significant_mask_positive = falses(size(mask_positive))
    significant_mask_negative = falses(size(mask_negative))
    for cluster in positive_clusters
        cluster.is_significant || continue
        @inbounds for pixel in cluster.pixels
            significant_mask_positive[pixel] = true
        end
    end
    for cluster in negative_clusters
        cluster.is_significant || continue
        @inbounds for pixel in cluster.pixels
            significant_mask_negative[pixel] = true
        end
    end

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
        TFMasks(significant_mask_positive, significant_mask_negative),
        TFClusters(positive_clusters, negative_clusters),
        PermutationDistribution(permutation_max_positive, permutation_max_negative),
        electrodes,
        frequencies,
        time_points,
    )

    @info "TF Permutation test complete. Found $(length(positive_clusters)) positive and $(length(negative_clusters)) negative clusters."
    return result
end

# === ANALYTIC TEST (ERP & TF) ===

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
        n_comp = 0
        n_sig_raw = 0
        @inbounds for i in eachindex(p_matrix)
            p = p_matrix[i]
            if !isnan(p)
                n_comp += 1
                p <= alpha_f && (n_sig_raw += 1)
            end
        end
        bonf_alpha = n_comp > 0 ? alpha_f / n_comp : 0.0
        n_sig_bonf = 0
        @inbounds for i in eachindex(p_matrix)
            p = p_matrix[i]
            (!isnan(p) && p <= bonf_alpha) && (n_sig_bonf += 1)
        end
        @info "Bonferroni correction: n_comparisons=$n_comp, bonferroni_alpha=$bonf_alpha, raw_sig=$n_sig_raw, bonferroni_sig=$n_sig_bonf"
        bonf_alpha
    else
        alpha_f
    end

    mask_positive = BitArray{N}(undef, size(t_matrix))
    mask_negative = BitArray{N}(undef, size(t_matrix))
    if tail == :both
        @inbounds for i in eachindex(t_matrix)
            p = p_matrix[i]
            t = t_matrix[i]
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

const _compute_tf_cluster_pvalues = _compute_cluster_pvalues
const _create_significance_mask_tf = _create_significance_mask
