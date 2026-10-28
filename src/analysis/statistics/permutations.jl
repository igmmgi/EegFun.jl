# This file contains permutation logic for Monte Carlo statistical testing,
# including label shuffling, t-matrix collection, and permutation loop execution.
#
# Randomness: all permutations are drawn up front on the calling task via
# `_draw_permutation_plan(rng, ...)`. The (multi-threaded) permutation loops then only
# *consume* columns of that plan, so results depend solely on `rng` and not on the number
# of threads or on task scheduling.

"""
    _draw_permutation_plan(rng, design, n_A, n_B, n_permutations)

Pre-draw all label permutations for a two-condition permutation test.

# Returns
- `:paired` → `Matrix{Bool}` of size `n_A × n_permutations`; `true` = swap the two conditions
  for that participant in that permutation (`n_A == n_B` = number of participants).
- `:independent` → `Matrix{Int}` of size `(n_A + n_B) × n_permutations`; each column is a
  permutation of the pooled participants (first `n_A` entries form group A).
"""
function _draw_permutation_plan(rng::AbstractRNG, design::Symbol, n_A::Int, n_B::Int, n_permutations::Int)
    if design == :paired
        return rand(rng, Bool, n_A, n_permutations)
    elseif design == :independent
        n_total = n_A + n_B
        plan = Matrix{Int}(undef, n_total, n_permutations)
        buf = Vector{Int}(undef, n_total)
        for p = 1:n_permutations
            randperm!(rng, buf)
            plan[:, p] .= buf
        end
        return plan
    end
    error("Unknown design: $design (expected :paired or :independent)")
end

_draw_permutation_plan(rng::AbstractRNG, prepared, n_permutations::Int) = _draw_permutation_plan(
    rng,
    prepared.analysis.design,
    size(prepared.analysis.data[1], 1),
    size(prepared.analysis.data[2], 1),
    n_permutations,
)

# Draw a single permutation (one plan column) — used by the non-plan convenience methods.
_draw_permutation(rng::AbstractRNG, design::Symbol, n_A::Int, n_B::Int) = vec(_draw_permutation_plan(rng, design, n_A, n_B, 1))

"""
    _shuffle_labels!(shuffled_A, shuffled_B, data1, data2, design, permutation)
    _shuffle_labels!(shuffled_A, shuffled_B, data1, data2, design; rng = Random.default_rng())

In-place version: Shuffle condition/group labels for permutation test.

# Arguments
- `shuffled_A::Array{Float64, 3}`: Output buffer for shuffled condition A
- `shuffled_B::Array{Float64, 3}`: Output buffer for shuffled condition B
- `data1::Array{Float64, 3}`: Original data for condition 1
- `data2::Array{Float64, 3}`: Original data for condition 2
- `design::Symbol`: Design type (`:paired` or `:independent`)
- `permutation::AbstractVector`: One column of a plan from [`_draw_permutation_plan`](@ref).
  If omitted, a fresh permutation is drawn from `rng`.

# Returns
- `(shuffled_A, shuffled_B)`: Shuffled data arrays

# Notes
For paired design: swap A/B labels within each participant where `permutation[p]` is `true`.
For independent design: reassign pooled participants to groups according to `permutation`.
"""
function _shuffle_labels!(
    shuffled_A::Array{Float64,N},
    shuffled_B::Array{Float64,N},
    data1::Array{Float64,N},
    data2::Array{Float64,N},
    design::Symbol;
    rng::AbstractRNG = Random.default_rng(),
) where {N}
    permutation = _draw_permutation(rng, design, size(data1, 1), size(data2, 1))
    return _shuffle_labels!(shuffled_A, shuffled_B, data1, data2, design, permutation)
end

function _shuffle_labels!(
    shuffled_A::Array{Float64,N},
    shuffled_B::Array{Float64,N},
    data1::Array{Float64,N},
    data2::Array{Float64,N},
    design::Symbol,
    permutation::AbstractVector,
) where {N}
    n_features = prod(size(data1)[2:end])
    n_A = size(data1, 1)
    n_B = size(data2, 1)

    if design == :paired
        copyto!(shuffled_A, data1)
        copyto!(shuffled_B, data2)

        A_mat = reshape(shuffled_A, n_A, n_features)
        B_mat = reshape(shuffled_B, n_B, n_features)
        swaps = permutation

        for feat = 1:n_features
            @inbounds @simd for p_idx = 1:n_A
                if swaps[p_idx]
                    temp = A_mat[p_idx, feat]
                    A_mat[p_idx, feat] = B_mat[p_idx, feat]
                    B_mat[p_idx, feat] = temp
                end
            end
        end

    elseif design == :independent
        A_mat = reshape(shuffled_A, n_A, n_features)
        B_mat = reshape(shuffled_B, n_B, n_features)
        d1_mat = reshape(data1, n_A, n_features)
        d2_mat = reshape(data2, n_B, n_features)
        shuffled_indices = permutation

        @inbounds for feat = 1:n_features
            @simd for i = 1:n_A
                src_idx = shuffled_indices[i]
                A_mat[i, feat] = src_idx <= n_A ? d1_mat[src_idx, feat] : d2_mat[src_idx-n_A, feat]
            end
            @simd for i = 1:n_B
                src_idx = shuffled_indices[n_A+i]
                B_mat[i, feat] = src_idx <= n_A ? d1_mat[src_idx, feat] : d2_mat[src_idx-n_A, feat]
            end
        end
    end
    return shuffled_A, shuffled_B
end

"""
    _shuffle_labels(prepared; rng = Random.default_rng())

Shuffle condition/group labels for permutation test.
"""
function _shuffle_labels(prepared::StatisticalData; rng::AbstractRNG = Random.default_rng())
    shuffled_A = similar(prepared.analysis.data[1])
    shuffled_B = similar(prepared.analysis.data[2])
    return _shuffle_labels!(
        shuffled_A,
        shuffled_B,
        prepared.analysis.data[1],
        prepared.analysis.data[2],
        prepared.analysis.design;
        rng = rng,
    )
end

# Start an (optional) progress bar fed by an atomic counter; returns the counter or `nothing`.
function _start_permutation_progress(show_progress::Bool, n_permutations::Int, desc::String)
    show_progress || return nothing
    progress = Progress(n_permutations, desc = desc, showspeed = true)
    completed = Threads.Atomic{Int}(0)
    track_progress_async(progress, completed, n_permutations)
    return completed
end

_tick!(::Nothing) = nothing
_tick!(completed::Threads.Atomic{Int}) = (Threads.atomic_add!(completed, 1); nothing)

"""
    _collect_permutation_t_matrices(prepared, plan, show_progress)

Run permutations and collect t-matrices for non-parametric thresholding.

# Arguments
- `prepared::StatisticalData`: Prepared data
- `plan::AbstractMatrix`: Pre-drawn permutations from [`_draw_permutation_plan`](@ref) (one column per permutation)
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]
"""
function _collect_permutation_t_matrices(prepared::StatisticalData, plan::AbstractMatrix, show_progress::Bool = true)
    n_permutations = size(plan, 2)
    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design
    paired = design == :paired
    n_electrodes = size(data1, 2)
    n_time = size(data1, 3)

    permutation_t_matrices = Array{Float64,3}(undef, n_electrodes, n_time, n_permutations)
    completed = _start_permutation_progress(show_progress, n_permutations, "Collecting permutation t-matrices: ")

    _foreach_chunk(n_permutations) do perm_range
        # Scratch buffers owned exclusively by this chunk/task
        shuffled_A = similar(data1)
        shuffled_B = similar(data2)
        t_matrix_buffer = Matrix{Float64}(undef, n_electrodes, n_time)
        mean1_buffer = paired ? similar(t_matrix_buffer) : nothing
        mean2_buffer = paired ? similar(t_matrix_buffer) : nothing
        mean_diff_buffer = paired ? similar(t_matrix_buffer) : nothing
        std_diff_buffer = paired ? similar(t_matrix_buffer) : nothing

        for perm_idx in perm_range
            _shuffle_labels!(shuffled_A, shuffled_B, data1, data2, design, view(plan, :, perm_idx))
            t_matrix_perm, _, _, _ = _compute_t_matrix(
                shuffled_A,
                shuffled_B,
                design,
                mean1_buffer = mean1_buffer,
                mean2_buffer = mean2_buffer,
                mean_diff_buffer = mean_diff_buffer,
                std_diff_buffer = std_diff_buffer,
                t_matrix_buffer = t_matrix_buffer,
                compute_p_values = false,
                compute_se = false,
            )
            permutation_t_matrices[:, :, perm_idx] = t_matrix_perm
            _tick!(completed)
        end
    end

    return permutation_t_matrices
end

"""
    _collect_permutation_t_matrices_gpu(prepared, plan)

GPU-accelerated collection of permutation t-matrices.
Computes permutations using massive batched matrix multiplications instead of CPU loops.
Uses the same pre-drawn `plan` as the CPU path, so CPU and GPU see identical permutations.
"""
function _collect_permutation_t_matrices_gpu(prepared::StatisticalData, plan::AbstractMatrix)
    n_permutations = size(plan, 2)
    n_electrodes = size(prepared.analysis.data[1], 2)
    n_time = size(prepared.analysis.data[1], 3)
    n_features = n_electrodes * n_time

    if prepared.analysis.design == :paired
        n_participants = size(prepared.analysis.data[1], 1)

        # Calculate D = data1 - data2
        D_cpu = reshape(prepared.analysis.data[1] .- prepared.analysis.data[2], n_participants, n_features)

        D_gpu = gpu_array(Float32.(D_cpu)) # n_participants x n_features

        # Sum of squares is invariant
        sum_sq_D = sum(D_gpu .^ 2, dims = 1) # 1 x n_features

        # Sign matrix S (+1 or -1): swapping conditions for a participant flips the sign of D
        S_cpu = Float32[plan[p_idx, perm_idx] ? -1.0f0 : 1.0f0 for perm_idx = 1:n_permutations, p_idx = 1:n_participants]
        S_gpu = gpu_array(S_cpu)

        # Matrix multiplication computes permutation sums simultaneously!
        sum_D_perm = S_gpu * D_gpu # n_permutations x n_features

        mean_D_perm = sum_D_perm ./ Float32(n_participants)

        # var = (sum(x^2)/n - mean^2) * n/(n-1)
        var_D_perm = @. (sum_sq_D / Float32(n_participants) - mean_D_perm^2) * (Float32(n_participants) / Float32(n_participants - 1))

        # Avoid negative variance due to floating point precision
        var_D_perm = max.(var_D_perm, 0.0f0)

        std_D_perm = sqrt.(var_D_perm)
        se_D_perm = std_D_perm ./ sqrt(Float32(n_participants))

        # Avoid division by zero
        t_perm_gpu = @. ifelse(se_D_perm == 0.0f0, ifelse(mean_D_perm == 0.0f0, 0.0f0, Inf32), mean_D_perm / se_D_perm)

        # Download and reshape back to [n_electrodes, n_time, n_permutations]
        t_perm_cpu = Float64.(Array(t_perm_gpu)) # n_permutations x n_features
        return copy(reshape(t_perm_cpu', n_electrodes, n_time, n_permutations))

    elseif prepared.analysis.design == :independent
        n_A = size(prepared.analysis.data[1], 1)
        n_B = size(prepared.analysis.data[2], 1)
        n_total = n_A + n_B

        data1_flat = reshape(prepared.analysis.data[1], n_A, n_features)
        data2_flat = reshape(prepared.analysis.data[2], n_B, n_features)
        X_cpu = vcat(data1_flat, data2_flat)

        X_gpu = gpu_array(Float32.(X_cpu)) # n_total x n_features
        X2_gpu = X_gpu .^ 2

        sum_tot = sum(X_gpu, dims = 1) # 1 x n_features
        sum_sq_tot = sum(X2_gpu, dims = 1) # 1 x n_features

        # Generate binary indicator matrix M on CPU
        M_cpu = zeros(Float32, n_permutations, n_total)
        for p = 1:n_permutations
            M_cpu[p, view(plan, 1:n_A, p)] .= 1.0f0
        end
        M_gpu = gpu_array(M_cpu)

        # Matrix multiplications to compute partitioned sums
        sum_A = M_gpu * X_gpu # n_permutations x n_features
        sum_sq_A = M_gpu * X2_gpu # n_permutations x n_features

        sum_B = sum_tot .- sum_A
        sum_sq_B = sum_sq_tot .- sum_sq_A

        mean_A = sum_A ./ Float32(n_A)
        mean_B = sum_B ./ Float32(n_B)

        var_A = @. (sum_sq_A / Float32(n_A) - mean_A^2) * (Float32(n_A) / Float32(n_A - 1))
        var_B = @. (sum_sq_B / Float32(n_B) - mean_B^2) * (Float32(n_B) / Float32(n_B - 1))

        var_A = max.(var_A, 0.0f0)
        var_B = max.(var_B, 0.0f0)

        df = Float32(n_total - 2)
        pooled_var = @. ((Float32(n_A - 1) * var_A) + (Float32(n_B - 1) * var_B)) / df
        se = sqrt.(pooled_var .* (1.0f0 / Float32(n_A) + 1.0f0 / Float32(n_B)))

        diff_mean = mean_A .- mean_B
        t_perm_gpu = @. ifelse(se == 0.0f0, ifelse(diff_mean == 0.0f0, 0.0f0, Inf32), diff_mean / se)

        t_perm_cpu = Float64.(Array(t_perm_gpu))
        return copy(reshape(t_perm_cpu', n_electrodes, n_time, n_permutations))
    end
end


"""
    _run_permutations(prepared, plan, critical_t_values, spatial_connectivity, cluster_type, tail, min_num_neighbors, show_progress; permutation_t_matrices)

Run permutation loop to generate distribution of maximum cluster statistics (maxsum).

# Arguments
- `prepared::StatisticalData`: Prepared data
- `plan::AbstractMatrix`: Pre-drawn permutations from [`_draw_permutation_plan`](@ref) (one column per permutation)
- `critical_t_values`: Critical t-values - can be:
  - `Array{Float64, 2}` for parametric (common critical t-values)
  - `Tuple{Float64, Float64}` for non-parametric common (positive, negative thresholds)
  - `Tuple{Array{Float64, 2}, Array{Float64, 2}}` for non-parametric individual (positive, negative threshold matrices)
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Connectivity matrix
- `cluster_type::Symbol`: Type of clustering
- `tail::Symbol`: Test tail
- `min_num_neighbors::Int`: Minimum number of neighbors for pre-filtering
- `show_progress::Bool`: Show progress bar (default: true)
- `permutation_t_matrices::Union{Nothing, Array{Float64, 3}}`: Pre-computed t-matrices from permutations
  (optional; must have been computed from the same `plan`)

# Returns
- `permutation_max_positive::Vector{Float64}`: Max cluster stats from permutations (positive)
- `permutation_max_negative::Vector{Float64}`: Max cluster stats from permutations (negative)
"""
function _run_permutations(
    prepared::StatisticalData,
    plan::AbstractMatrix,
    critical_t_values,
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    tail::Symbol,
    min_num_neighbors::Int,
    show_progress::Bool = true;
    permutation_t_matrices::Union{Nothing,Array{Float64,3}} = nothing,
)
    n_permutations = size(plan, 2)

    # Determine thresholding type from critical_t_values type
    is_parametric = isa(critical_t_values, Array{Float64,2})
    is_nonparametric_common = isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && isa(critical_t_values[1], Float64)
    is_nonparametric_individual =
        isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && isa(critical_t_values[1], Array{Float64,2})

    permutation_max_positive = zeros(Float64, n_permutations)
    permutation_max_negative = zeros(Float64, n_permutations)

    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design
    paired = design == :paired
    n_electrodes = size(data1, 2)
    n_time = size(data1, 3)

    # Ensure spatial connectivity is symmetric once outside permutation loop
    spatial_conn_sym =
        issymmetric(spatial_connectivity) ? spatial_connectivity : (spatial_connectivity .| copy(transpose(spatial_connectivity)))

    completed = _start_permutation_progress(show_progress, n_permutations, "Permutations: ")

    _foreach_chunk(n_permutations) do perm_range
        # Scratch buffers owned exclusively by this chunk/task
        shuffled_A = similar(data1)
        shuffled_B = similar(data2)
        mask_pos_buffer = BitArray{2}(undef, n_electrodes, n_time)
        mask_neg_buffer = BitArray{2}(undef, n_electrodes, n_time)
        t_matrix_buffer = Matrix{Float64}(undef, n_electrodes, n_time)
        mean1_buffer = paired ? similar(t_matrix_buffer) : nothing
        mean2_buffer = paired ? similar(t_matrix_buffer) : nothing
        mean_diff_buffer = paired ? similar(t_matrix_buffer) : nothing
        std_diff_buffer = paired ? similar(t_matrix_buffer) : nothing
        cluster_labels = zeros(Int, n_electrodes, n_time)
        cluster_queue = Vector{Tuple{Int,Int}}(undef, n_electrodes * n_time)

        for perm_idx in perm_range
            # Get t-matrix: either from pre-computed or compute new
            if !isnothing(permutation_t_matrices)
                copyto!(t_matrix_buffer, view(permutation_t_matrices, :, :, perm_idx))
                t_matrix_perm = t_matrix_buffer
            else
                _shuffle_labels!(shuffled_A, shuffled_B, data1, data2, design, view(plan, :, perm_idx))
                t_matrix_perm, _, _, _ = _compute_t_matrix(
                    shuffled_A,
                    shuffled_B,
                    design,
                    mean1_buffer = mean1_buffer,
                    mean2_buffer = mean2_buffer,
                    mean_diff_buffer = mean_diff_buffer,
                    std_diff_buffer = std_diff_buffer,
                    t_matrix_buffer = t_matrix_buffer,
                    compute_p_values = false,
                    compute_se = false,
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
                _prefilter_mask_by_neighbors!(mask_pos_buffer, spatial_conn_sym, min_num_neighbors)
                _prefilter_mask_by_neighbors!(mask_neg_buffer, spatial_conn_sym, min_num_neighbors)
            end

            # Find extreme cluster statistics directly (zero allocation)
            permutation_max_positive[perm_idx] = _max_cluster_stat(
                mask_pos_buffer,
                t_matrix_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = false,
            )
            permutation_max_negative[perm_idx] = _max_cluster_stat(
                mask_neg_buffer,
                t_matrix_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = true,
            )
            _tick!(completed)
        end
    end

    return permutation_max_positive, permutation_max_negative
end


# Forward TF label shuffling to unified N-D _shuffle_labels!
const _shuffle_labels_tf! = _shuffle_labels!


# === TF PERMUTATION LOOP ===

"""
    _run_permutations_tf(prepared, plan, critical_t_values, spatial_connectivity,
                          cluster_type, tail, min_num_neighbors, electrodes, frequencies,
                          time_points, show_progress)

Run Monte Carlo permutations for TF cluster-based permutation test.

`plan` holds the pre-drawn permutations from [`_draw_permutation_plan`](@ref) (one column per permutation).

# Returns
- `permutation_max_positive::Vector{Float64}`: Max positive cluster stats per permutation
- `permutation_max_negative::Vector{Float64}`: Max negative cluster stats per permutation
"""
function _run_permutations_tf(
    prepared::TFStatisticalData,
    plan::AbstractMatrix,
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
    n_permutations = size(plan, 2)
    is_parametric = isa(critical_t_values, Array{Float64,3})
    is_nonparametric_common = isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && isa(critical_t_values[1], Float64)

    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    permutation_max_positive = zeros(Float64, n_permutations)
    permutation_max_negative = zeros(Float64, n_permutations)

    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design

    # Ensure spatial connectivity is symmetric once outside permutation loop
    spatial_conn_sym =
        issymmetric(spatial_connectivity) ? spatial_connectivity : (spatial_connectivity .| copy(transpose(spatial_connectivity)))

    completed = _start_permutation_progress(show_progress, n_permutations, "TF Permutations: ")

    _foreach_chunk(n_permutations) do perm_range
        # Scratch buffers owned exclusively by this chunk/task
        mask_pos_buffer = BitArray{3}(undef, n_electrodes, n_freqs, n_time)
        mask_neg_buffer = BitArray{3}(undef, n_electrodes, n_freqs, n_time)
        shuffled_A = similar(data1)
        shuffled_B = similar(data2)
        t_matrix_buffer = Array{Float64,3}(undef, n_electrodes, n_freqs, n_time)
        cluster_labels = zeros(Int, n_electrodes, n_freqs, n_time)
        cluster_queue = Vector{Tuple{Int,Int,Int}}(undef, n_electrodes * n_freqs * n_time)

        for perm_idx in perm_range
            _shuffle_labels_tf!(shuffled_A, shuffled_B, data1, data2, design, view(plan, :, perm_idx))

            t_matrix_perm, _, _ = _compute_t_matrix_tf(
                shuffled_A,
                shuffled_B,
                design,
                tail = tail,
                compute_p_values = false,
                t_matrix_buffer = t_matrix_buffer,
            )

            # Threshold
            if is_parametric
                _threshold_t_matrix_parametric_tf!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, critical_t_values, tail)
            elseif is_nonparametric_common
                thresh_pos, thresh_neg = critical_t_values
                _threshold_t_matrix_nonparametric_tf!(mask_pos_buffer, mask_neg_buffer, t_matrix_perm, thresh_pos, thresh_neg, tail)
            end

            # Pre-filter
            if min_num_neighbors > 0
                _prefilter_mask_by_neighbors_tf!(mask_pos_buffer, spatial_conn_sym, min_num_neighbors)
                _prefilter_mask_by_neighbors_tf!(mask_neg_buffer, spatial_conn_sym, min_num_neighbors)
            end

            # Find extreme cluster statistics directly (zero allocation)
            permutation_max_positive[perm_idx] = _max_cluster_stat_tf(
                mask_pos_buffer,
                t_matrix_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = false,
            )
            permutation_max_negative[perm_idx] = _max_cluster_stat_tf(
                mask_neg_buffer,
                t_matrix_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = true,
            )
            _tick!(completed)
        end
    end

    return permutation_max_positive, permutation_max_negative
end

"""
    generate_permutation_matrix(df, block_col; sync_col=nothing, n_perms=1000, type=:within, rng=Random.GLOBAL_RNG)

Generate a mathematically valid `N × n_perms` permutation matrix of integer row indices for Monte Carlo / permutation tests.

# Arguments
- `df`: A `DataFrame` containing the experimental design / trial metadata.
- `block_col`: Column identifier (Symbol or String) for the grouping block (e.g. `:Subject` or `:Item`).

# Keyword Arguments
- `sync_col=nothing`: Column identifier to synchronize across for crossed designs (e.g. `sync_col=:Subject, block_col=:Item`).
- `n_perms=1000`: Number of permutations to generate.
- `type=:within`: Permutation scheme:
  - `:within`: Independent shuffling within blocks (e.g. within Subjects). Standard for uncrossed repeated measures.
  - `:synchronized`: Shuffles `block_col` (e.g. Items) and perfectly synchronizes that shuffle across `sync_col` (e.g. Subjects).
    Requires a balanced crossed design where every subject sees the items.
  - `:global`: Unrestricted random permutation across all rows.
- `rng=Random.GLOBAL_RNG`: Random number generator for reproducibility.

# Returns
- `Matrix{Int}`: An `N × n_perms` matrix where each column is a valid permutation of `1:N`.
"""
function _permute_within!(perm_matrix::Matrix{Int}, block_idxs::Vector{Vector{Int}}, n_perms::Int, rng)
    scratch_bufs = [copy(b) for b in block_idxs]
    for p = 1:n_perms
        for (k, b_idx) in enumerate(block_idxs)
            buf = scratch_bufs[k]
            m = length(buf)
            for i = m:-1:2
                j = rand(rng, 1:i)
                buf[i], buf[j] = buf[j], buf[i]
            end
            @inbounds for i = 1:m
                perm_matrix[b_idx[i], p] = buf[i]
            end
        end
    end
end

function _permute_synchronized!(perm_matrix::Matrix{Int}, grid::Matrix{Vector{Int}}, cell_len::Int, n_perms::Int, rng)
    n_syncs, n_blocks = size(grid)
    shuffled_idx = collect(1:n_blocks)
    for p = 1:n_perms
        for i = n_blocks:-1:2
            j = rand(rng, 1:i)
            shuffled_idx[i], shuffled_idx[j] = shuffled_idx[j], shuffled_idx[i]
        end
        for s = 1:n_syncs
            for b = 1:n_blocks
                target_b = shuffled_idx[b]
                orig_rows = grid[s, b]
                target_rows = grid[s, target_b]
                @inbounds for r = 1:cell_len
                    perm_matrix[orig_rows[r], p] = target_rows[r]
                end
            end
        end
    end
end

function _permute_global!(perm_matrix::Matrix{Int}, n_perms::Int, rng)
    N = size(perm_matrix, 1)
    for p = 1:n_perms
        col = @view perm_matrix[:, p]
        for i = 1:N
            col[i] = i
        end
        for i = N:-1:2
            j = rand(rng, 1:i)
            col[i], col[j] = col[j], col[i]
        end
    end
end

function generate_permutation_matrix(df, block_col; sync_col = nothing, n_perms = 1000, type = :within, rng = Random.GLOBAL_RNG)
    N = nrow(df)
    n_perms <= 0 && return zeros(Int, N, 0)

    perm_matrix = zeros(Int, N, n_perms)
    b_col = Symbol(block_col)
    hasproperty(df, b_col) || @minimal_error("Column '$b_col' not found in DataFrame. Available columns: $(names(df))")
    blocks = df[!, b_col]
    unique_blocks = unique(blocks)
    n_blocks = length(unique_blocks)

    if type == :within
        block_idxs = [findall(==(b), blocks) for b in unique_blocks]
        _permute_within!(perm_matrix, block_idxs, n_perms, rng)
    elseif type == :synchronized
        isnothing(sync_col) && @minimal_error("Synchronized shuffling requires `sync_col` (e.g. sync_col=:Subject, block_col=:Item).")
        s_col = Symbol(sync_col)
        hasproperty(df, s_col) || @minimal_error("Column '$s_col' not found in DataFrame. Available columns: $(names(df))")

        syncs = df[!, s_col]
        unique_syncs = unique(syncs)
        n_syncs = length(unique_syncs)

        sync_map = Dict(s => i for (i, s) in enumerate(unique_syncs))
        block_map = Dict(b => i for (i, b) in enumerate(unique_blocks))

        grid = [Int[] for _ = 1:n_syncs, _ = 1:n_blocks]
        for i = 1:N
            s_idx = sync_map[syncs[i]]
            b_idx = block_map[blocks[i]]
            push!(grid[s_idx, b_idx], i)
        end

        cell_len = length(grid[1, 1])
        for s = 1:n_syncs, b = 1:n_blocks
            if length(grid[s, b]) != cell_len
                @minimal_error(
                    "Synchronized shuffling requires a balanced design where every sync group sees all items with equal trial count. Row maps unbalanced cell count for Sync: $(unique_syncs[s]), Block: $(unique_blocks[b]). For unbalanced designs, please supply a custom `perm_matrix`."
                )
            end
        end

        _permute_synchronized!(perm_matrix, grid, cell_len, n_perms, rng)
    elseif type == :global
        _permute_global!(perm_matrix, n_perms, rng)
    else
        @minimal_error("Unknown permutation type: $type. Supported types: :within, :synchronized, :global")
    end

    return perm_matrix
end
