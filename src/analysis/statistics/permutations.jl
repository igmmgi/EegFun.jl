# Monte Carlo permutation planning, permutation t-statistic kernels, and block permutation matrices.
#
# All permutations are drawn up front on the calling task via `_draw_permutation_plan(rng, ...)`,
# so threaded execution is deterministic for a given `rng` regardless of thread count.

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

_draw_permutation(rng::AbstractRNG, design::Symbol, n_A::Int, n_B::Int) = vec(_draw_permutation_plan(rng, design, n_A, n_B, 1))

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
        @inbounds for feat = 1:n_features
            @simd for p_idx = 1:n_A
                if permutation[p_idx]
                    A_mat[p_idx, feat], B_mat[p_idx, feat] = B_mat[p_idx, feat], A_mat[p_idx, feat]
                end
            end
        end
    elseif design == :independent
        A_mat = reshape(shuffled_A, n_A, n_features)
        B_mat = reshape(shuffled_B, n_B, n_features)
        d1_mat = reshape(data1, n_A, n_features)
        d2_mat = reshape(data2, n_B, n_features)
        @inbounds for feat = 1:n_features
            @simd for i = 1:n_A
                src_idx = permutation[i]
                A_mat[i, feat] = src_idx <= n_A ? d1_mat[src_idx, feat] : d2_mat[src_idx-n_A, feat]
            end
            @simd for i = 1:n_B
                src_idx = permutation[n_A+i]
                B_mat[i, feat] = src_idx <= n_A ? d1_mat[src_idx, feat] : d2_mat[src_idx-n_A, feat]
            end
        end
    end
    return shuffled_A, shuffled_B
end

function _start_permutation_progress(show_progress::Bool, n_permutations::Int, desc::String)
    show_progress || return nothing
    progress = Progress(n_permutations, desc = desc, showspeed = true)
    completed = Threads.Atomic{Int}(0)
    track_progress_async(progress, completed, n_permutations)
    return completed
end

_tick!(::Nothing) = nothing
_tick!(completed::Threads.Atomic{Int}) = (Threads.atomic_add!(completed, 1); nothing)

# === CPU PERMUTATION T-STATISTIC KERNELS ===

function _perm_ttest_paired!(t_flat::AbstractVector{Float64}, D::Matrix{Float64}, signs::Vector{Float64}, swaps::AbstractVector)
    n_p, n_voxels = size(D)
    @inbounds @simd for p = 1:n_p
        signs[p] = ifelse(Bool(swaps[p]), -1.0, 1.0)
    end
    inv_n = 1.0 / n_p
    inv_var_se = 1.0 / ((n_p - 1) * n_p)

    @inbounds for v = 1:n_voxels
        sum_d = 0.0
        @simd for p = 1:n_p
            sum_d += signs[p] * D[p, v]
        end
        mean_d = sum_d * inv_n

        ss_d = 0.0
        @simd for p = 1:n_p
            delta = signs[p] * D[p, v] - mean_d
            ss_d += delta * delta
        end

        se = sqrt(max(ss_d * inv_var_se, 0.0))
        t_flat[v] = se == 0.0 ? (mean_d == 0.0 ? NaN : copysign(Inf, mean_d)) : (mean_d / se)
    end
    return nothing
end

function _perm_ttest_indep!(t_flat::AbstractVector{Float64}, X::Matrix{Float64}, n_A::Int, n_B::Int, perm::AbstractVector)
    n_voxels = size(X, 2)
    inv_A = 1.0 / n_A
    inv_B = 1.0 / n_B
    se_scale = (inv_A + inv_B) / (n_A + n_B - 2)

    @inbounds for v = 1:n_voxels
        sum1 = 0.0
        for i = 1:n_A
            sum1 += X[perm[i], v]
        end
        sum2 = 0.0
        for i = 1:n_B
            sum2 += X[perm[n_A+i], v]
        end
        mean1 = sum1 * inv_A
        mean2 = sum2 * inv_B

        ss1 = 0.0
        for i = 1:n_A
            d1 = X[perm[i], v] - mean1
            ss1 += d1 * d1
        end
        ss2 = 0.0
        for i = 1:n_B
            d2 = X[perm[n_A+i], v] - mean2
            ss2 += d2 * d2
        end

        se = sqrt(max((ss1 + ss2) * se_scale, 0.0))
        diff_mean = mean1 - mean2
        t_flat[v] = se == 0.0 ? (diff_mean == 0.0 ? NaN : copysign(Inf, diff_mean)) : (diff_mean / se)
    end
    return nothing
end

function _build_perm_source(data1::Array{Float64,M}, data2::Array{Float64,M}, design::Symbol) where {M}
    n_A = size(data1, 1)
    n_B = size(data2, 1)
    n_voxels = div(length(data1), n_A)
    if design == :paired
        D = Matrix{Float64}(undef, n_A, n_voxels)
        d1 = reshape(data1, n_A, n_voxels)
        d2 = reshape(data2, n_A, n_voxels)
        @inbounds @simd for i in eachindex(D)
            D[i] = d1[i] - d2[i]
        end
        return D
    else
        X = Matrix{Float64}(undef, n_A + n_B, n_voxels)
        d1 = reshape(data1, n_A, n_voxels)
        d2 = reshape(data2, n_B, n_voxels)
        @inbounds for v = 1:n_voxels
            @simd for p = 1:n_A
                X[p, v] = d1[p, v]
            end
            @simd for p = 1:n_B
                X[n_A+p, v] = d2[p, v]
            end
        end
        return X
    end
end

function _collect_perm_t_generic(
    data1::Array{Float64,M},
    data2::Array{Float64,M},
    design::Symbol,
    plan::AbstractMatrix,
    show_progress::Bool,
    desc::String,
) where {M}
    n_permutations = size(plan, 2)
    n_A = size(data1, 1)
    n_B = size(data2, 1)
    grid_dims = ntuple(d -> size(data1, d + 1), Val(M - 1))
    n_voxels = prod(grid_dims)

    out = Array{Float64,M}(undef, grid_dims..., n_permutations)
    out_flat = reshape(out, n_voxels, n_permutations)
    src = _build_perm_source(data1, data2, design)
    paired = design == :paired
    completed = _start_permutation_progress(show_progress, n_permutations, desc)

    _foreach_chunk(n_permutations) do perm_range
        signs = paired ? Vector{Float64}(undef, n_A) : Float64[]
        for perm_idx in perm_range
            col_out = view(out_flat, :, perm_idx)
            col_plan = view(plan, :, perm_idx)
            if paired
                _perm_ttest_paired!(col_out, src, signs, col_plan)
            else
                _perm_ttest_indep!(col_out, src, n_A, n_B, col_plan)
            end
            _tick!(completed)
        end
    end

    return out
end

_collect_permutation_t_matrices(prepared::StatisticalData, plan::AbstractMatrix, show_progress::Bool = true) = _collect_perm_t_generic(
    prepared.analysis.data[1],
    prepared.analysis.data[2],
    prepared.analysis.design,
    plan,
    show_progress,
    "Collecting permutation t-matrices: ",
)

_collect_permutation_t_matrices_tf(
    data1::Array{Float64,4},
    data2::Array{Float64,4},
    design::Symbol,
    plan::AbstractMatrix,
    show_progress::Bool = true,
) = _collect_perm_t_generic(data1, data2, design, plan, show_progress, "Collecting TF t-matrices: ")

# === GPU PERMUTATION KERNELS ===

@kernel function _paired_perm_ttest_kernel!(t_out, @Const(D), @Const(plan), n_participants::Int)
    feat_idx, perm_idx = @index(Global, NTuple)

    sum_d = 0.0f0
    @inbounds for p = 1:n_participants
        val = D[p, feat_idx]
        sum_d += plan[p, perm_idx] ? -val : val
    end
    mean_d = sum_d / Float32(n_participants)

    ss_d = 0.0f0
    @inbounds for p = 1:n_participants
        val = D[p, feat_idx]
        delta = (plan[p, perm_idx] ? -val : val) - mean_d
        ss_d += delta * delta
    end

    var_d = max(ss_d / Float32(n_participants - 1), 0.0f0)
    se_d = sqrt(var_d / Float32(n_participants))
    @inbounds t_out[feat_idx, perm_idx] = se_d == 0.0f0 ? (mean_d == 0.0f0 ? NaN32 : copysign(Inf32, mean_d)) : (mean_d / se_d)
end

@kernel function _indep_perm_ttest_kernel!(t_out, @Const(X), @Const(plan), n_A::Int, n_B::Int)
    feat_idx, perm_idx = @index(Global, NTuple)

    sum1 = 0.0f0
    @inbounds for i = 1:n_A
        sum1 += X[plan[i, perm_idx], feat_idx]
    end
    sum2 = 0.0f0
    @inbounds for i = 1:n_B
        sum2 += X[plan[n_A+i, perm_idx], feat_idx]
    end

    mean1 = sum1 / Float32(n_A)
    mean2 = sum2 / Float32(n_B)

    ss1 = 0.0f0
    @inbounds for i = 1:n_A
        d1 = X[plan[i, perm_idx], feat_idx] - mean1
        ss1 += d1 * d1
    end
    ss2 = 0.0f0
    @inbounds for i = 1:n_B
        d2 = X[plan[n_A+i, perm_idx], feat_idx] - mean2
        ss2 += d2 * d2
    end

    df = Float32(n_A + n_B - 2)
    pooled_var = max((ss1 + ss2) / df, 0.0f0)
    se = sqrt(pooled_var * (1.0f0 / Float32(n_A) + 1.0f0 / Float32(n_B)))
    diff_mean = mean1 - mean2
    @inbounds t_out[feat_idx, perm_idx] = se == 0.0f0 ? (diff_mean == 0.0f0 ? NaN32 : copysign(Inf32, diff_mean)) : (diff_mean / se)
end

function _collect_permutation_t_matrices_gpu(prepared::StatisticalData, plan::AbstractMatrix)
    n_permutations = size(plan, 2)
    n_electrodes = size(prepared.analysis.data[1], 2)
    n_time = size(prepared.analysis.data[1], 3)
    n_features = n_electrodes * n_time
    backend = gpu_backend()

    if prepared.analysis.design == :paired
        n_participants = size(prepared.analysis.data[1], 1)
        D_cpu = _build_perm_source(prepared.analysis.data[1], prepared.analysis.data[2], :paired)
        D_gpu = gpu_array(Float32.(D_cpu))
        plan_gpu = gpu_array(Matrix{Bool}(plan))
        t_perm_gpu = KernelAbstractions.allocate(backend, Float32, n_features, n_permutations)

        _paired_perm_ttest_kernel!(backend)(t_perm_gpu, D_gpu, plan_gpu, n_participants; ndrange = (n_features, n_permutations))
        KernelAbstractions.synchronize(backend)
        return reshape(Float64.(Array(t_perm_gpu)), n_electrodes, n_time, n_permutations)

    elseif prepared.analysis.design == :independent
        n_A = size(prepared.analysis.data[1], 1)
        n_B = size(prepared.analysis.data[2], 1)
        X_cpu = _build_perm_source(prepared.analysis.data[1], prepared.analysis.data[2], :independent)
        X_centered = X_cpu .- mean(X_cpu, dims = 1)

        X_gpu = gpu_array(Float32.(X_centered))
        plan_gpu = gpu_array(Matrix{Int32}(plan))
        t_perm_gpu = KernelAbstractions.allocate(backend, Float32, n_features, n_permutations)

        _indep_perm_ttest_kernel!(backend)(t_perm_gpu, X_gpu, plan_gpu, n_A, n_B; ndrange = (n_features, n_permutations))
        KernelAbstractions.synchronize(backend)
        return reshape(Float64.(Array(t_perm_gpu)), n_electrodes, n_time, n_permutations)
    end
end

# === CLUSTER PERMUTATION LOOPS (ERP & TF) ===

@inline _apply_perm_threshold!(
    mask_pos::BitArray{N},
    mask_neg::BitArray{N},
    t_mat::AbstractArray{Float64,N},
    crit::Array{Float64,N},
    tail::Symbol,
) where {N} = _threshold_t_matrix_parametric!(mask_pos, mask_neg, t_mat, crit, tail)

@inline _apply_perm_threshold!(
    mask_pos::BitArray{N},
    mask_neg::BitArray{N},
    t_mat::AbstractArray{Float64,N},
    crit::Tuple{Float64,Float64},
    tail::Symbol,
) where {N} = _threshold_t_matrix_nonparametric!(mask_pos, mask_neg, t_mat, crit[1], crit[2], tail)

@inline _apply_perm_threshold!(
    mask_pos::BitArray{N},
    mask_neg::BitArray{N},
    t_mat::AbstractArray{Float64,N},
    crit::Tuple{Array{Float64,N},Array{Float64,N}},
    tail::Symbol,
) where {N} = _threshold_t_matrix_nonparametric!(mask_pos, mask_neg, t_mat, crit[1], crit[2], tail)

function _run_permutations(
    prepared::StatisticalData,
    plan::AbstractMatrix,
    critical_t_values::C,
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    tail::Symbol,
    min_num_neighbors::Int,
    show_progress::Bool = true;
    permutation_t_matrices::Union{Nothing,Array{Float64,3}} = nothing,
) where {C}
    n_permutations = size(plan, 2)
    permutation_max_positive = zeros(Float64, n_permutations)
    permutation_max_negative = zeros(Float64, n_permutations)

    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design
    paired = design == :paired
    n_A, n_electrodes, n_time = size(data1)
    n_B = size(data2, 1)

    src = isnothing(permutation_t_matrices) ? _build_perm_source(data1, data2, design) : Matrix{Float64}(undef, 0, 0)
    spatial_conn_sym = _symmetrize_connectivity(spatial_connectivity)
    completed = _start_permutation_progress(show_progress, n_permutations, "Permutations: ")

    _foreach_chunk(n_permutations) do perm_range
        mask_pos = BitMatrix(undef, n_electrodes, n_time)
        mask_neg = BitMatrix(undef, n_electrodes, n_time)
        t_buf = isnothing(permutation_t_matrices) ? Matrix{Float64}(undef, n_electrodes, n_time) : Matrix{Float64}(undef, 0, 0)
        signs = (isnothing(permutation_t_matrices) && paired) ? Vector{Float64}(undef, n_A) : Float64[]
        cluster_labels = zeros(Int, n_electrodes, n_time)
        cluster_queue = Vector{Tuple{Int,Int}}(undef, n_electrodes * n_time)

        for perm_idx in perm_range
            if !isnothing(permutation_t_matrices)
                t_perm = view(permutation_t_matrices, :, :, perm_idx)
            else
                col_plan = view(plan, :, perm_idx)
                if paired
                    _perm_ttest_paired!(vec(t_buf), src, signs, col_plan)
                else
                    _perm_ttest_indep!(vec(t_buf), src, n_A, n_B, col_plan)
                end
                t_perm = t_buf
            end

            _apply_perm_threshold!(mask_pos, mask_neg, t_perm, critical_t_values, tail)

            if min_num_neighbors > 0 && cluster_type ∈ (:spatial, :spatiotemporal, :full)
                _prefilter_mask_by_neighbors!(mask_pos, spatial_conn_sym, min_num_neighbors)
                _prefilter_mask_by_neighbors!(mask_neg, spatial_conn_sym, min_num_neighbors)
            end

            permutation_max_positive[perm_idx] = _max_cluster_stat(
                mask_pos,
                t_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = false,
                labels_clean = true,
            )
            permutation_max_negative[perm_idx] = _max_cluster_stat(
                mask_neg,
                t_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = true,
                labels_clean = true,
            )
            _tick!(completed)
        end
    end

    return permutation_max_positive, permutation_max_negative
end

const _TF_PERM_CACHE_MAX_FLOATS = 16_000_000 # 128 MB cap for 4D TF permutation t-matrix cache

function _stream_nonparametric_threshold_common_tf(
    prepared::TFStatisticalData,
    plan::AbstractMatrix,
    alpha::Float64,
    tail::Symbol,
    show_progress::Bool = true,
)
    n_permutations = size(plan, 2)
    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design
    paired = design == :paired
    n_A, n_electrodes, n_freqs, n_time = size(data1)
    n_B = size(data2, 1)
    n_voxels = n_electrodes * n_freqs * n_time

    stride = max(1, cld(n_voxels * n_permutations, _TF_PERM_CACHE_MAX_FLOATS))
    samples_per_perm = cld(n_voxels, stride)
    pooled = fill(NaN, samples_per_perm, n_permutations)
    src = _build_perm_source(data1, data2, design)
    completed = _start_permutation_progress(show_progress, n_permutations, "Streaming TF thresholds: ")

    _foreach_chunk(n_permutations) do perm_range
        t_buf = Vector{Float64}(undef, n_voxels)
        signs = paired ? Vector{Float64}(undef, n_A) : Float64[]
        for perm_idx in perm_range
            col_plan = view(plan, :, perm_idx)
            if paired
                _perm_ttest_paired!(t_buf, src, signs, col_plan)
            else
                _perm_ttest_indep!(t_buf, src, n_A, n_B, col_plan)
            end
            offset = mod(perm_idx - 1, stride) + 1
            out_idx = 0
            @inbounds for v = offset:stride:n_voxels
                out_idx += 1
                pooled[out_idx, perm_idx] = t_buf[v]
            end
            _tick!(completed)
        end
    end

    return _compute_nonparametric_threshold_common_tf(pooled, alpha, tail)
end

function _stream_nonparametric_threshold_individual_tf(
    prepared::TFStatisticalData,
    plan::AbstractMatrix,
    alpha::Float64,
    tail::Symbol,
    show_progress::Bool = true,
)
    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design
    _, n_electrodes, n_freqs, n_time = size(data1)

    thresholds_pos = Array{Float64,3}(undef, n_electrodes, n_freqs, n_time)
    thresholds_neg = Array{Float64,3}(undef, n_electrodes, n_freqs, n_time)

    for e = 1:n_electrodes
        perm_e = _collect_permutation_t_matrices_tf(data1[:, e:e, :, :], data2[:, e:e, :, :], design, plan, show_progress && e == 1)
        tp_e, tn_e = _compute_nonparametric_threshold_individual_tf(perm_e, alpha, tail)
        thresholds_pos[e:e, :, :] .= tp_e
        thresholds_neg[e:e, :, :] .= tn_e
    end

    return thresholds_pos, thresholds_neg
end

function _run_permutations_tf(
    prepared::TFStatisticalData,
    plan::AbstractMatrix,
    critical_t_values::C,
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    tail::Symbol,
    min_num_neighbors::Int,
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    show_progress::Bool = true;
    permutation_t_matrices::Union{Nothing,Array{Float64,4}} = nothing,
) where {C}
    n_permutations = size(plan, 2)
    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    permutation_max_positive = zeros(Float64, n_permutations)
    permutation_max_negative = zeros(Float64, n_permutations)

    data1, data2 = prepared.analysis.data[1], prepared.analysis.data[2]
    design = prepared.analysis.design
    paired = design == :paired
    n_A = size(data1, 1)
    n_B = size(data2, 1)

    src = isnothing(permutation_t_matrices) ? _build_perm_source(data1, data2, design) : Matrix{Float64}(undef, 0, 0)
    spatial_conn_sym = _symmetrize_connectivity(spatial_connectivity)
    completed = _start_permutation_progress(show_progress, n_permutations, "TF Permutations: ")

    _foreach_chunk(n_permutations) do perm_range
        mask_pos = BitArray{3}(undef, n_electrodes, n_freqs, n_time)
        mask_neg = BitArray{3}(undef, n_electrodes, n_freqs, n_time)
        t_buf =
            isnothing(permutation_t_matrices) ? Array{Float64,3}(undef, n_electrodes, n_freqs, n_time) : Array{Float64,3}(undef, 0, 0, 0)
        signs = (isnothing(permutation_t_matrices) && paired) ? Vector{Float64}(undef, n_A) : Float64[]
        cluster_labels = zeros(Int, n_electrodes, n_freqs, n_time)
        cluster_queue = Vector{Tuple{Int,Int,Int}}(undef, n_electrodes * n_freqs * n_time)

        for perm_idx in perm_range
            if !isnothing(permutation_t_matrices)
                t_perm = view(permutation_t_matrices, :, :, :, perm_idx)
            else
                col_plan = view(plan, :, perm_idx)
                if paired
                    _perm_ttest_paired!(vec(t_buf), src, signs, col_plan)
                else
                    _perm_ttest_indep!(vec(t_buf), src, n_A, n_B, col_plan)
                end
                t_perm = t_buf
            end

            _apply_perm_threshold!(mask_pos, mask_neg, t_perm, critical_t_values, tail)

            if min_num_neighbors > 0 && cluster_type ∈ (:spatial, :spatiotemporal, :full)
                _prefilter_mask_by_neighbors_tf!(mask_pos, spatial_conn_sym, min_num_neighbors)
                _prefilter_mask_by_neighbors_tf!(mask_neg, spatial_conn_sym, min_num_neighbors)
            end

            permutation_max_positive[perm_idx] = _max_cluster_stat_tf(
                mask_pos,
                t_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = false,
                labels_clean = true,
            )
            permutation_max_negative[perm_idx] = _max_cluster_stat_tf(
                mask_neg,
                t_perm,
                spatial_conn_sym,
                cluster_type,
                cluster_labels,
                cluster_queue;
                is_negative = true,
                labels_clean = true,
            )
            _tick!(completed)
        end
    end

    return permutation_max_positive, permutation_max_negative
end

# === DESIGN-AWARE BLOCK PERMUTATION MATRIX GENERATOR ===

function _permute_within!(perm_matrix::Matrix{Int}, block_idxs::Vector{Vector{Int}}, n_perms::Int, rng::AbstractRNG)
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

function _permute_synchronized!(perm_matrix::Matrix{Int}, grid::Matrix{Vector{Int}}, cell_len::Int, n_perms::Int, rng::AbstractRNG)
    n_syncs, n_blocks = size(grid)
    shuffled_idx = collect(1:n_blocks)
    for p = 1:n_perms
        for i = n_blocks:-1:2
            j = rand(rng, 1:i)
            shuffled_idx[i], shuffled_idx[j] = shuffled_idx[j], shuffled_idx[i]
        end
        for s = 1:n_syncs, b = 1:n_blocks
            orig_rows = grid[s, b]
            target_rows = grid[s, shuffled_idx[b]]
            @inbounds for r = 1:cell_len
                perm_matrix[orig_rows[r], p] = target_rows[r]
            end
        end
    end
end

function _permute_global!(perm_matrix::Matrix{Int}, n_perms::Int, rng::AbstractRNG)
    N = size(perm_matrix, 1)
    for p = 1:n_perms
        col = @view perm_matrix[:, p]
        @inbounds for i = 1:N
            col[i] = i
        end
        for i = N:-1:2
            j = rand(rng, 1:i)
            col[i], col[j] = col[j], col[i]
        end
    end
end

"""
    generate_permutation_matrix(df, block_col; sync_col = nothing, n_perms = 1000, type = :within, rng = Random.default_rng())

Generate an `N × n_perms` permutation matrix of integer row indices for Monte Carlo / permutation tests.

# Arguments
- `df`: A `DataFrame` containing the experimental design / trial metadata.
- `block_col`: Column identifier (`Symbol` or `String`) for the grouping block (e.g. `:Subject` or `:Item`).

# Keyword Arguments
- `sync_col = nothing`: Column identifier to synchronize across for balanced crossed designs.
- `n_perms = 1000`: Number of permutations to generate.
- `type = :within`: `:within` (shuffle within blocks), `:synchronized` (synchronized shuffle across `sync_col`), or `:global`.
- `rng = Random.default_rng()`: Random number generator for reproducibility.

# Returns
- `Matrix{Int}`: An `N × n_perms` matrix where each column is a valid permutation of `1:N`.
"""
function generate_permutation_matrix(
    df,
    block_col;
    sync_col = nothing,
    n_perms = 1000,
    type = :within,
    rng::AbstractRNG = Random.default_rng(),
)
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
            push!(grid[sync_map[syncs[i]], block_map[blocks[i]]], i)
        end

        cell_len = length(grid[1, 1])
        for s = 1:n_syncs, b = 1:n_blocks
            length(grid[s, b]) == cell_len || @minimal_error(
                "Synchronized shuffling requires a balanced design where every sync group sees all items with equal trial count. Row maps unbalanced cell count for Sync: $(unique_syncs[s]), Block: $(unique_blocks[b]). For unbalanced designs, please supply a custom `perm_matrix`."
            )
        end

        _permute_synchronized!(perm_matrix, grid, cell_len, n_perms, rng)
    elseif type == :global
        _permute_global!(perm_matrix, n_perms, rng)
    else
        @minimal_error("Unknown permutation type: $type. Supported types: :within, :synchronized, :global")
    end

    return perm_matrix
end
