# This file contains permutation logic for Monte Carlo statistical testing,
# including label shuffling, t-matrix collection, and permutation loop execution.
"""
    _generate_swap_mask(n_participants)

Generate a random swap mask for paired design.

# Arguments
- `n_participants::Int`: Number of participants

# Returns
- `BitVector`: Swap mask (true = swap, false = don't swap)
"""
function _generate_swap_mask(n_participants::Int)
    return BitVector([rand(Bool) for _ = 1:n_participants])
end

"""
    _shuffle_labels!(shuffled_A, shuffled_B, data1, data2, design)

In-place version: Shuffle condition/group labels for permutation test.

# Arguments
- `shuffled_A::Array{Float64, 3}`: Output buffer for shuffled condition A
- `shuffled_B::Array{Float64, 3}`: Output buffer for shuffled condition B
- `data1::Array{Float64, 3}`: Original data for condition 1
- `data2::Array{Float64, 3}`: Original data for condition 2
- `design::Symbol`: Design type (`:paired` or `:independent`)

# Returns
- `(shuffled_A, shuffled_B)`: Shuffled data arrays

# Notes
For paired design: randomly swap A/B labels within each participant.
For independent design: randomly shuffle participants between groups.
"""
function _shuffle_labels!(
    shuffled_A::Array{Float64,3},
    shuffled_B::Array{Float64,3},
    data1::Array{Float64,3},
    data2::Array{Float64,3},
    design::Symbol,
)
    if design == :paired
        # Paired: copy data first (copyto! is optimized), then swap slices in-place
        copyto!(shuffled_A, data1)
        copyto!(shuffled_B, data2)

        n_participants, n_electrodes, n_time = size(data1)
        swaps = rand(Bool, n_participants)
        
        for t_idx in 1:n_time
            for e_idx in 1:n_electrodes
                @inbounds @simd for p_idx = 1:n_participants
                    if swaps[p_idx]
                        temp = shuffled_A[p_idx, e_idx, t_idx]
                        shuffled_A[p_idx, e_idx, t_idx] = shuffled_B[p_idx, e_idx, t_idx]
                        shuffled_B[p_idx, e_idx, t_idx] = temp
                    end
                end
            end
        end

    elseif design == :independent
        n_A = size(data1, 1)
        n_B = size(data2, 1)
        n_total = n_A + n_B
        n_electrodes = size(data1, 2)
        n_time = size(data1, 3)

        # Shuffle indices
        shuffled_indices = collect(1:n_total)
        shuffle!(shuffled_indices)

        @inbounds for t_idx in 1:n_time
            for e_idx in 1:n_electrodes
                @simd for i = 1:n_A
                    src_idx = shuffled_indices[i]
                    if src_idx <= n_A
                        shuffled_A[i, e_idx, t_idx] = data1[src_idx, e_idx, t_idx]
                    else
                        shuffled_A[i, e_idx, t_idx] = data2[src_idx - n_A, e_idx, t_idx]
                    end
                end
                
                @simd for i = 1:n_B
                    src_idx = shuffled_indices[n_A + i]
                    if src_idx <= n_A
                        shuffled_B[i, e_idx, t_idx] = data1[src_idx, e_idx, t_idx]
                    else
                        shuffled_B[i, e_idx, t_idx] = data2[src_idx - n_A, e_idx, t_idx]
                    end
                end
            end
        end
    end
    return shuffled_A, shuffled_B
end

"""
    _shuffle_labels(prepared)

Shuffle condition/group labels for permutation test.
"""
function _shuffle_labels(prepared::StatisticalData)
    shuffled_A = similar(prepared.analysis.data[1])
    shuffled_B = similar(prepared.analysis.data[2])
    return _shuffle_labels!(shuffled_A, shuffled_B, prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design)
end

"""
    _collect_permutation_t_matrices(prepared, n_permutations, show_progress)

Run permutations and collect t-matrices for non-parametric thresholding.

# Arguments
- `prepared::StatisticalData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]

# Notes
For reproducible results, call `Random.seed!(xxx)` in your Julia session before running this function.

# Examples
```julia
perm_t_matrices = _collect_permutation_t_matrices(prepared, 1000, nothing, true)
```
"""
function _collect_permutation_t_matrices(prepared::StatisticalData, n_permutations::Int, show_progress::Bool = true)
    # Get dimensions directly from data
    n_electrodes = size(prepared.analysis.data[1], 2)
    n_time = size(prepared.analysis.data[1], 3)

    # Pre-allocate array for all t-matrices
    permutation_t_matrices = Array{Float64,3}(undef, n_electrodes, n_time, n_permutations)

    # Progress bar
    if show_progress
        progress = Progress(n_permutations, desc = "Collecting permutation t-matrices: ", showspeed = true)
    end

    # Pre-allocate buffers per thread to reuse memory safely across concurrent threads
    n_threads = Threads.maxthreadid()
    shuffled_A_buffers = [similar(prepared.analysis.data[1]) for _ in 1:n_threads]
    shuffled_B_buffers = [similar(prepared.analysis.data[2]) for _ in 1:n_threads]
    
    mean1_buffers = prepared.analysis.design == :paired ? [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads] : nothing
    mean2_buffers = prepared.analysis.design == :paired ? [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads] : nothing
    mean_diff_buffers = prepared.analysis.design == :paired ? [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads] : nothing
    std_diff_buffers = prepared.analysis.design == :paired ? [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads] : nothing
    t_matrix_buffers = [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads]

    Threads.@threads for perm_idx = 1:n_permutations
        tid = Threads.threadid()
        shuffled_A_buffer = shuffled_A_buffers[tid]
        shuffled_B_buffer = shuffled_B_buffers[tid]
        t_matrix_buffer = t_matrix_buffers[tid]

        # Shuffle labels using pre-allocated thread-local buffers
        _shuffle_labels!(
            shuffled_A_buffer,
            shuffled_B_buffer,
            prepared.analysis.data[1],
            prepared.analysis.data[2],
            prepared.analysis.design,
        )

        mean1_buffer = prepared.analysis.design == :paired ? mean1_buffers[tid] : nothing
        mean2_buffer = prepared.analysis.design == :paired ? mean2_buffers[tid] : nothing
        mean_diff_buffer = prepared.analysis.design == :paired ? mean_diff_buffers[tid] : nothing
        std_diff_buffer = prepared.analysis.design == :paired ? std_diff_buffers[tid] : nothing

        # Compute t-matrix directly from arrays using thread-local memory
        t_matrix_perm, _, _, _ = _compute_t_matrix(
            shuffled_A_buffer, 
            shuffled_B_buffer, 
            prepared.analysis.design,
            mean1_buffer = mean1_buffer,
            mean2_buffer = mean2_buffer,
            mean_diff_buffer = mean_diff_buffer,
            std_diff_buffer = std_diff_buffer,
            t_matrix_buffer = t_matrix_buffer
        )
        permutation_t_matrices[:, :, perm_idx] = t_matrix_perm

        if show_progress
            next!(progress)
        end
    end

    return permutation_t_matrices
end


"""
    _run_permutations(prepared, n_permutations, critical_t_values, spatial_connectivity, cluster_type, tail, min_num_neighbors, show_progress; permutation_t_matrices)

Run permutation loop to generate distribution of maximum cluster statistics (maxsum).

# Arguments
- `prepared::StatisticalData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `threshold::Float64`: P-value threshold
- `critical_t_values`: Critical t-values - can be:
  - `Array{Float64, 2}` for parametric (common critical t-values)
  - `Tuple{Float64, Float64}` for non-parametric common (positive, negative thresholds)
  - `Tuple{Array{Float64, 2}, Array{Float64, 2}}` for non-parametric individual (positive, negative threshold matrices)
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Connectivity matrix
- `cluster_type::Symbol`: Type of clustering
- `tail::Symbol`: Test tail
- `min_num_neighbors::Int`: Minimum number of neighbors for pre-filtering
- `show_progress::Bool`: Show progress bar (default: true)
- `permutation_t_matrices::Union{Nothing, Array{Float64, 3}}`: Pre-computed t-matrices from permutations (optional, for non-parametric)

# Returns
- `permutation_max_positive::Vector{Float64}`: Max cluster stats from permutations (positive)
- `permutation_max_negative::Vector{Float64}`: Max cluster stats from permutations (negative)

# Notes
For reproducible results, call `Random.seed!(xxx)` in your Julia session before running this function.

# Examples
```julia
perm_pos, perm_neg = _run_permutations(prepared, 1000, 0.05, critical_t, conn, 
                                     :spatiotemporal, :both, 0, 0, true)
```
"""
function _run_permutations(
    prepared::StatisticalData,
    n_permutations::Int,
    critical_t_values,
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    tail::Symbol,
    min_num_neighbors::Int,
    show_progress::Bool = true;
    permutation_t_matrices::Union{Nothing,Array{Float64,3}} = nothing,
)
    # Determine thresholding type from critical_t_values type
    is_parametric = isa(critical_t_values, Array{Float64,2})
    is_nonparametric_common = isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && isa(critical_t_values[1], Float64)
    is_nonparametric_individual =
        isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && isa(critical_t_values[1], Array{Float64,2})

    permutation_max_positive = zeros(Float64, n_permutations)
    permutation_max_negative = zeros(Float64, n_permutations)

    n_electrodes = size(prepared.analysis.data[1], 2)
    n_time = size(prepared.analysis.data[1], 3)
    n_threads = Threads.maxthreadid()

    # Pre-allocate buffers per thread to reuse memory across concurrent permutations
    shuffled_A_buffers = [similar(prepared.analysis.data[1]) for _ in 1:n_threads]
    shuffled_B_buffers = [similar(prepared.analysis.data[2]) for _ in 1:n_threads]
    mask_pos_buffers = [BitArray{2}(undef, n_electrodes, n_time) for _ in 1:n_threads]
    mask_neg_buffers = [BitArray{2}(undef, n_electrodes, n_time) for _ in 1:n_threads]

    if prepared.analysis.design == :paired
        mean1_buffers = [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads]
        mean2_buffers = [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads]
        mean_diff_buffers = [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads]
        std_diff_buffers = [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads]
    else
        mean1_buffers = mean2_buffers = mean_diff_buffers = std_diff_buffers = nothing
    end
    t_matrix_buffers = [Array{Float64,2}(undef, n_electrodes, n_time) for _ in 1:n_threads]

    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points

    # Pre-allocate electrode lookup (reused across all permutations)
    electrode_to_idx = Dict(e => i for (i, e) in enumerate(electrodes))

    # Progress bar
    if show_progress
        progress = Progress(n_permutations, desc = "Permutations: ", showspeed = true)
    end

    Threads.@threads for perm_idx = 1:n_permutations
        tid = Threads.threadid()
        
        # Select thread-local buffers
        shuffled_A_buffer = shuffled_A_buffers[tid]
        shuffled_B_buffer = shuffled_B_buffers[tid]
        mask_pos_buffer = mask_pos_buffers[tid]
        mask_neg_buffer = mask_neg_buffers[tid]
        mean1_buffer = prepared.analysis.design == :paired ? mean1_buffers[tid] : nothing
        mean2_buffer = prepared.analysis.design == :paired ? mean2_buffers[tid] : nothing
        mean_diff_buffer = prepared.analysis.design == :paired ? mean_diff_buffers[tid] : nothing
        std_diff_buffer = prepared.analysis.design == :paired ? std_diff_buffers[tid] : nothing
        t_matrix_buffer = t_matrix_buffers[tid]

        # Get t-matrix: either from pre-computed or compute new
        if !isnothing(permutation_t_matrices)
            t_matrix_perm = permutation_t_matrices[:, :, perm_idx]
        else
            # Shuffle and compute t-matrix
            _shuffle_labels!(
                shuffled_A_buffer,
                shuffled_B_buffer,
                prepared.analysis.data[1],
                prepared.analysis.data[2],
                prepared.analysis.design,
            )
            t_matrix_perm, _, _, _ = _compute_t_matrix(
                shuffled_A_buffer,
                shuffled_B_buffer,
                prepared.analysis.design,
                mean1_buffer = mean1_buffer,
                mean2_buffer = mean2_buffer,
                mean_diff_buffer = mean_diff_buffer,
                std_diff_buffer = std_diff_buffer,
                t_matrix_buffer = t_matrix_buffer,
                compute_p_values = false
            )
        end

        # Threshold (fill pre-allocated masks)
        if is_parametric
            _threshold_t_matrix_parametric!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, critical_t_values, tail)
        elseif is_nonparametric_common
            thresh_pos, thresh_neg = critical_t_values
            _threshold_t_matrix_nonparametric!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, thresh_pos, thresh_neg, tail)
        elseif is_nonparametric_individual
            thresh_pos_mat, thresh_neg_mat = critical_t_values
            _threshold_t_matrix_nonparametric!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, thresh_pos_mat, thresh_neg_mat, tail)
        end

        # Pre-filter masks (modify in place)
        if min_num_neighbors > 0
            _prefilter_mask_by_neighbors!(mask_pos_buffer, spatial_connectivity, min_num_neighbors)
            _prefilter_mask_by_neighbors!(mask_neg_buffer, spatial_connectivity, min_num_neighbors)
        end

        # Find clusters
        pos_clusters_perm, neg_clusters_perm =
            _find_clusters(mask_pos_buffer, mask_neg_buffer, electrodes, time_points, spatial_connectivity, cluster_type)

        # Compute statistics without creating new Cluster objects (use pre-allocated lookup)
        if !isempty(pos_clusters_perm)
            pos_stats_perm = _compute_cluster_statistics(pos_clusters_perm, t_matrix_perm, electrode_to_idx, return_clusters = false)
            max_pos = isempty(pos_stats_perm) ? 0.0 : maximum(pos_stats_perm)
        else
            max_pos = 0.0
        end

        if !isempty(neg_clusters_perm)
            neg_stats_perm = _compute_cluster_statistics(neg_clusters_perm, t_matrix_perm, electrode_to_idx, return_clusters = false)
            # For negative clusters, the statistic is negative (sum of negative t-values)
            # The most extreme negative value is the MINIMUM (most negative), not maximum
            max_neg = isempty(neg_stats_perm) ? 0.0 : minimum(neg_stats_perm)
        else
            max_neg = 0.0
        end

        permutation_max_positive[perm_idx] = max_pos
        permutation_max_negative[perm_idx] = max_neg

        if show_progress
            next!(progress)
        end
    end

    return permutation_max_positive, permutation_max_negative
end


# ===================
# TF LABEL SHUFFLING (4D)
# ===================

"""
    _shuffle_labels_tf!(shuffled_A, shuffled_B, data1, data2, design)

In-place label shuffling for TF permutation tests.

For paired design: randomly swap conditions for each participant.
For independent design: randomly shuffle participants between groups.
"""
function _shuffle_labels_tf!(
    shuffled_A::Array{Float64,4},
    shuffled_B::Array{Float64,4},
    data1::Array{Float64,4},
    data2::Array{Float64,4},
    design::Symbol,
)
    if design == :paired
        copyto!(shuffled_A, data1)
        copyto!(shuffled_B, data2)

        n_participants, n_electrodes, n_freqs, n_time = size(data1)
        swaps = rand(Bool, n_participants)
        
        for t_idx in 1:n_time
            for f_idx in 1:n_freqs
                for e_idx in 1:n_electrodes
                    @inbounds @simd for p_idx = 1:n_participants
                        if swaps[p_idx]
                            temp = shuffled_A[p_idx, e_idx, f_idx, t_idx]
                            shuffled_A[p_idx, e_idx, f_idx, t_idx] = shuffled_B[p_idx, e_idx, f_idx, t_idx]
                            shuffled_B[p_idx, e_idx, f_idx, t_idx] = temp
                        end
                    end
                end
            end
        end

    elseif design == :independent
        n_A = size(data1, 1)
        n_B = size(data2, 1)
        n_total = n_A + n_B

        shuffled_indices = randperm(n_total)

        n_electrodes = size(data1, 2)
        n_freqs = size(data1, 3)
        n_time = size(data1, 4)

        @inbounds for t_idx in 1:n_time
            for f_idx in 1:n_freqs
                for e_idx in 1:n_electrodes
                    @simd for i = 1:n_A
                        src_idx = shuffled_indices[i]
                        if src_idx <= n_A
                            shuffled_A[i, e_idx, f_idx, t_idx] = data1[src_idx, e_idx, f_idx, t_idx]
                        else
                            shuffled_A[i, e_idx, f_idx, t_idx] = data2[src_idx - n_A, e_idx, f_idx, t_idx]
                        end
                    end
                    @simd for i = 1:n_B
                        src_idx = shuffled_indices[n_A+i]
                        if src_idx <= n_A
                            shuffled_B[i, e_idx, f_idx, t_idx] = data1[src_idx, e_idx, f_idx, t_idx]
                        else
                            shuffled_B[i, e_idx, f_idx, t_idx] = data2[src_idx - n_A, e_idx, f_idx, t_idx]
                        end
                    end
                end
            end
        end
    end
    return shuffled_A, shuffled_B
end


# ===================
# TF PERMUTATION LOOP
# ===================

"""
    _run_permutations_tf(prepared, n_permutations, critical_t_values, spatial_connectivity,
                          cluster_type, tail, min_num_neighbors, electrodes, frequencies,
                          time_points, show_progress)

Run Monte Carlo permutations for TF cluster-based permutation test.

# Returns
- `permutation_max_positive::Vector{Float64}`: Max positive cluster stats per permutation
- `permutation_max_negative::Vector{Float64}`: Max negative cluster stats per permutation
"""
function _run_permutations_tf(
    prepared::TFStatisticalData,
    n_permutations::Int,
    critical_t_values,  # Can be Array{Float64,3} or Tuple{Float64,Float64}
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    tail::Symbol,
    min_num_neighbors::Int,
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    show_progress::Bool = true,
)
    is_parametric = isa(critical_t_values, Array{Float64,3})
    is_nonparametric_common = isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && isa(critical_t_values[1], Float64)

    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    permutation_max_positive = zeros(Float64, n_permutations)
    permutation_max_negative = zeros(Float64, n_permutations)

    n_threads = Threads.maxthreadid()

    # Pre-allocate thread-local buffers
    mask_pos_buffers = [BitArray{3}(undef, n_electrodes, n_freqs, n_time) for _ in 1:n_threads]
    mask_neg_buffers = [BitArray{3}(undef, n_electrodes, n_freqs, n_time) for _ in 1:n_threads]
    shuffled_A_buffers = [similar(prepared.analysis.data[1]) for _ in 1:n_threads]
    shuffled_B_buffers = [similar(prepared.analysis.data[2]) for _ in 1:n_threads]
    t_matrix_buffers = [Array{Float64,3}(undef, n_electrodes, n_freqs, n_time) for _ in 1:n_threads]

    electrode_to_idx = Dict(e => i for (i, e) in enumerate(electrodes))

    if show_progress
        progress = Progress(n_permutations, desc = "TF Permutations: ", showspeed = true)
    end

    Threads.@threads for perm_idx = 1:n_permutations
        tid = Threads.threadid()
        mask_pos_buffer = mask_pos_buffers[tid]
        mask_neg_buffer = mask_neg_buffers[tid]
        shuffled_A = shuffled_A_buffers[tid]
        shuffled_B = shuffled_B_buffers[tid]
        t_matrix_buffer = t_matrix_buffers[tid]

        # Shuffle labels
        _shuffle_labels_tf!(shuffled_A, shuffled_B, prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design)

        # Compute t-matrix for shuffled data
        t_matrix_perm, _, _ = _compute_t_matrix_tf(shuffled_A, shuffled_B, prepared.analysis.design, tail = tail, compute_p_values = false, t_matrix_buffer = t_matrix_buffer)

        # Threshold
        if is_parametric
            _threshold_t_matrix_parametric_tf!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, critical_t_values, tail)
        elseif is_nonparametric_common
            thresh_pos, thresh_neg = critical_t_values
            _threshold_t_matrix_nonparametric_tf!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, thresh_pos, thresh_neg, tail)
        end

        # Pre-filter
        if min_num_neighbors > 0
            _prefilter_mask_by_neighbors_tf!(mask_pos_buffer, spatial_connectivity, min_num_neighbors)
            _prefilter_mask_by_neighbors_tf!(mask_neg_buffer, spatial_connectivity, min_num_neighbors)
        end

        # Find clusters
        pos_clusters_perm, neg_clusters_perm =
            _find_clusters_tf(mask_pos_buffer, mask_neg_buffer, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)

        # Get max cluster stat (positive)
        if !isempty(pos_clusters_perm)
            pos_stats = _compute_cluster_statistics_tf(pos_clusters_perm, t_matrix_perm, false)
            max_pos = isempty(pos_stats) ? 0.0 : maximum(pos_stats)
        else
            max_pos = 0.0
        end

        # Get max cluster stat (negative)
        if !isempty(neg_clusters_perm)
            neg_stats = _compute_cluster_statistics_tf(neg_clusters_perm, t_matrix_perm, false)
            max_neg = isempty(neg_stats) ? 0.0 : minimum(neg_stats)
        else
            max_neg = 0.0
        end

        permutation_max_positive[perm_idx] = max_pos
        permutation_max_negative[perm_idx] = max_neg

        if show_progress
            next!(progress)
        end
    end

    return permutation_max_positive, permutation_max_negative
end
