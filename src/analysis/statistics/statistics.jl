# ====================================================================================
# MAIN TEST FUNCTIONS
# ====================================================================================
# This module contains the main statistical test functions that users call:
# cluster_permutation_test() and analytic_ttest().

"""
    permutation_test(prepared::StatisticalData; kwargs...)

Perform cluster-based permutation test on prepared ERP data.

# Arguments
- `prepared::StatisticalData`: Prepared data from `prepare_stats`
- `n_permutations::Int`: Number of permutations (default: 1000)
- `threshold::Float64`: P-value threshold (default: 0.05)
- `threshold_method::Symbol`: Threshold method - `:parametric` (default), `:nonparametric_individual`, or `:nonparametric_common`
- `cluster_type::Symbol`: Type of clustering - `:spatial`, `:temporal`, or `:spatiotemporal` (default)
- `min_num_neighbors::Int`: Minimum number of neighboring significant channels required (FieldTrip's minNumChannels, default: 0). Points with fewer neighbors are removed before clustering.
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`
- `random_seed::Union{Int, Nothing}`: Random seed for reproducibility (default: nothing)
- `show_progress::Bool`: Whether to show progress bar (default: true)

# Returns
- `PermutationResult`: Complete results structure

# Notes
Cluster statistics are computed using the **maxsum** method (sum of t-values within cluster), which is the most sensitive and standard approach in EEG research.
"""
function permutation_test(
    prepared::StatisticalData;
    n_permutations::Int = 1000,
    threshold::Float64 = 0.05,
    threshold_method::Symbol = :parametric,
    cluster_type::Symbol = :spatiotemporal,
    min_num_neighbors::Int = 0,
    tail::Symbol = :both,
    random_seed::Union{Int,Nothing} = nothing,
    show_progress::Bool = true,
)
    # Validate inputs
    _validate_permutation_inputs(prepared, n_permutations, threshold, cluster_type, tail)

    # Validate threshold method
    if !(threshold_method in [:parametric, :nonparametric_common, :nonparametric_individual])
        error(
            "threshold_method must be :parametric, :nonparametric_common, or :nonparametric_individual. " *
            "Got :$threshold_method",
        )
    end

    # Compute observed t-matrix and df
    @info "Computing t-statistics..."
    t_matrix, df, _ = _compute_t_matrix(prepared)

    # Handle different thresholding methods
    if threshold_method == :parametric
        # Parametric thresholding: compute critical t-values from t-distribution
        @info "Computing parametric critical t-values..."
        critical_t_values = _compute_critical_t_values(df, size(t_matrix), threshold, tail)

        # Threshold observed data
        @info "Thresholding observed data (parametric)..."
        mask_positive, mask_negative = _threshold_t_matrix_parametric(t_matrix, critical_t_values, tail)

        # Store for later use in permutations
        threshold_for_permutations = critical_t_values
        permutation_t_matrices = nothing

    elseif threshold_method == :nonparametric_common
        # Non-parametric common: run all permutations first to get threshold
        @info "Collecting permutation t-matrices for non-parametric common thresholding..."
        permutation_t_matrices = _collect_permutation_t_matrices(prepared, n_permutations, random_seed, show_progress)

        # Compute common threshold from permutation distribution
        @info "Computing non-parametric common threshold..."
        thresh_pos, thresh_neg = _compute_nonparametric_threshold_common(permutation_t_matrices, threshold, tail)
        critical_t_values = (thresh_pos, thresh_neg)

        # Threshold observed data
        @info "Thresholding observed data (non-parametric common)..."
        mask_positive, mask_negative = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos, thresh_neg, tail)

        # Store for later use in permutations
        threshold_for_permutations = critical_t_values

    elseif threshold_method == :nonparametric_individual
        # Non-parametric individual: run all permutations first to get thresholds
        @info "Collecting permutation t-matrices for non-parametric individual thresholding..."
        permutation_t_matrices = _collect_permutation_t_matrices(prepared, n_permutations, random_seed, show_progress)

        # Compute individual thresholds from permutation distribution
        @info "Computing non-parametric individual thresholds..."
        thresh_pos_mat, thresh_neg_mat =
            _compute_nonparametric_threshold_individual(permutation_t_matrices, threshold, tail)
        critical_t_values = (thresh_pos_mat, thresh_neg_mat)

        # Threshold observed data
        @info "Thresholding observed data (non-parametric individual)..."
        mask_positive, mask_negative = _threshold_t_matrix_nonparametric(t_matrix, thresh_pos_mat, thresh_neg_mat, tail)

        # Store for later use in permutations
        threshold_for_permutations = critical_t_values
    end

    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points
    layout = prepared.data[1].layout

    # Build connectivity matrix
    @info "Building connectivity matrix..."
    spatial_connectivity, _, _ = _build_connectivity_matrix(electrodes, layout, cluster_type)

    # Pre-filter masks to remove isolated points (FieldTrip's minNumChannels approach)
    if min_num_neighbors > 0
        @info "Pre-filtering masks (min_num_neighbors=$min_num_neighbors)..."
        mask_positive = _prefilter_mask_by_neighbors(mask_positive, spatial_connectivity, min_num_neighbors)
        mask_negative = _prefilter_mask_by_neighbors(mask_negative, spatial_connectivity, min_num_neighbors)
    end

    # Find observed clusters
    @info "Finding observed clusters..."
    positive_clusters, negative_clusters =
        _find_clusters(mask_positive, mask_negative, electrodes, time_points, spatial_connectivity, cluster_type)

    # Compute observed cluster statistics (maxsum)
    cluster_stats_positive = Float64[]
    cluster_stats_negative = Float64[]

    if !isempty(positive_clusters)
        @info "Computing cluster statistics for $(length(positive_clusters)) positive clusters..."
        positive_clusters, cluster_stats_positive = _compute_cluster_statistics(positive_clusters, t_matrix, electrodes)
    end

    if !isempty(negative_clusters)
        @info "Computing cluster statistics for $(length(negative_clusters)) negative clusters..."
        negative_clusters, cluster_stats_negative = _compute_cluster_statistics(negative_clusters, t_matrix, electrodes)
    end

    # Run permutations for cluster-level inference
    if threshold_method == :parametric
        @info "Running $n_permutations permutations for cluster-level inference..."
    else
        @info "Running $n_permutations permutations for cluster-level inference (reusing stored t-matrices)..."
    end
    permutation_max_positive, permutation_max_negative = _run_permutations(
        prepared,
        n_permutations,
        threshold,
        threshold_for_permutations,
        spatial_connectivity,
        cluster_type,
        tail,
        min_num_neighbors,
        random_seed,
        show_progress;
        permutation_t_matrices = permutation_t_matrices,
    )

    # Compute p-values
    @info "Computing p-values..."
    if !isempty(positive_clusters)
        positive_clusters = _compute_cluster_pvalues(
            positive_clusters,
            cluster_stats_positive,
            permutation_max_positive,
            n_permutations,
            threshold,
        )
        sort!(positive_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end

    if !isempty(negative_clusters)
        negative_clusters = _compute_cluster_pvalues(
            negative_clusters,
            cluster_stats_negative,
            permutation_max_negative,
            n_permutations,
            threshold,
        )
        sort!(negative_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end

    # Create significance masks
    significant_mask_positive = falses(size(mask_positive))
    significant_mask_negative = falses(size(mask_negative))

    for cluster in positive_clusters
        if cluster.is_significant
            for electrode in cluster.electrodes
                e_idx = findfirst(==(electrode), electrodes)
                if e_idx !== nothing
                    for t_idx in cluster.time_indices
                        if 1 <= t_idx <= size(significant_mask_positive, 2)
                            significant_mask_positive[e_idx, t_idx] = true
                        end
                    end
                end
            end
        end
    end

    for cluster in negative_clusters
        if cluster.is_significant
            for electrode in cluster.electrodes
                e_idx = findfirst(==(electrode), electrodes)
                if e_idx !== nothing
                    for t_idx in cluster.time_indices
                        if 1 <= t_idx <= size(significant_mask_negative, 2)
                            significant_mask_negative[e_idx, t_idx] = true
                        end
                    end
                end
            end
        end
    end

    # Assemble nested structs
    cluster_info = ClusterInfo(threshold_method, cluster_type, n_permutations, random_seed)

    test_info = TestInfo(prepared.analysis.design, df, threshold, :both, :cluster_permutation, cluster_info)

    stat_matrix = StatMatrix(t_matrix, nothing)
    masks = Masks(significant_mask_positive, significant_mask_negative)
    clusters = Clusters(positive_clusters, negative_clusters)
    permutation_dist = PermutationDistribution(permutation_max_positive, permutation_max_negative)

    result = PermutationResult(
        test_info,
        prepared.data,
        stat_matrix,
        masks,
        clusters,
        permutation_dist,
        electrodes,
        time_points,
        critical_t_values,
    )

    @info "Permutation test complete. Found $(length(clusters.positive)) positive and $(length(clusters.negative)) negative clusters."

    return result
end
# ===================
# ANALYTIC TEST
# ===================

"""
    analytic_test(prepared::StatisticalData;
                  alpha::Float64 = 0.05,
                  tail::Symbol = :both,
                  correction_method::Symbol = :no)

Perform analytic (parametric) t-test without permutation (FieldTrip's 'analytic' method).

# Arguments
- `prepared::StatisticalData`: Prepared data from `prepare_stats`
- `alpha::Float64`: Significance threshold (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`
- `correction_method::Symbol`: Multiple comparison correction - `:no` (default) or `:bonferroni`

# Returns
- `AnalyticResult`: Results structure with t-statistics, p-values, and significant masks

# Examples
```julia
result = analytic_test(prepared, alpha=0.05, correction_method=:no)
result = analytic_test(prepared, alpha=0.05, correction_method=:bonferroni)
```
"""
function analytic_test(
    prepared::StatisticalData;
    alpha::Float64 = 0.05,
    tail::Symbol = :both,
    correction_method::Symbol = :no,
)

    correction_method ∉ (:no, :bonferroni) &&
        @minimal_error "correction_method must be :no or :bonferroni. Got :$correction_method"
    tail ∉ (:both, :left, :right) && @minimal_error "tail must be :both, :left, or :right. Got :$tail"

    # Compute t-statistics, degrees of freedom, and p-values in one pass
    t_matrix, df, p_matrix = _compute_t_matrix(prepared, tail = tail)

    # Create significance mask (with optional correction)
    corrected_mask = _create_significance_mask(p_matrix, alpha, correction_method)

    if correction_method == :bonferroni
        n_comparisons = count(!isnan, p_matrix)
        bonferroni_alpha = n_comparisons > 0 ? alpha / n_comparisons : 0.0
        n_sig_uncorrected = count(p -> !isnan(p) && p <= alpha, p_matrix)
        n_sig_bonferroni = count(corrected_mask)
        @info "Bonferroni correction: n_comparisons=$n_comparisons, bonferroni_alpha=$bonferroni_alpha, uncorrected_sig=$n_sig_uncorrected, bonferroni_sig=$n_sig_bonferroni"
    end

    # Create positive and negative masks based on t-values (vectorized)
    if tail == :both
        mask_positive = corrected_mask .& .!isnan.(t_matrix) .& (t_matrix .> 0)
        mask_negative = corrected_mask .& .!isnan.(t_matrix) .& (t_matrix .< 0)
    elseif tail == :right
        mask_positive = corrected_mask .& .!isnan.(t_matrix)
        mask_negative = falses(size(t_matrix))
    elseif tail == :left
        mask_positive = falses(size(t_matrix))
        mask_negative = corrected_mask .& .!isnan.(t_matrix)
    end

    test_info = TestInfo(
        prepared.analysis.design,
        df,
        alpha,
        tail,
        correction_method,
        nothing,  # cluster_info is nothing for analytic tests
    )

    stat_matrix = StatMatrix(t_matrix, p_matrix)
    masks = Masks(mask_positive, mask_negative)

    # Validate df before computing critical t-value
    if isnan(df) || isinf(df) || df <= 0
        @minimal_error "Invalid degrees of freedom: df=$df!"
    end

    # Compute critical t-value
    dist = TDist(df)
    if tail == :both
        critical_t = quantile(dist, 1.0 - alpha / 2.0)
    elseif tail == :right
        critical_t = quantile(dist, 1.0 - alpha)
    else  # :left
        critical_t = quantile(dist, alpha)
    end

    result = AnalyticResult(
        test_info,
        prepared.data,
        stat_matrix,
        masks,
        channel_labels(prepared.data[1]),
        prepared.analysis.time_points,
        critical_t,
    )

    @info "Analytic test complete: $(count(mask_positive)) +ve, $(count(mask_negative)) -ve significant points."

    return result
end
