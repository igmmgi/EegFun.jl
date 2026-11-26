# Connected-component clustering, neighbor pre-filtering, cluster statistics, and TFCE.

######################################
######### Internal Functions #########
######################################

_symmetrize_connectivity(conn::SparseMatrixCSC{Bool}) = issymmetric(conn) ? conn : (conn .| copy(transpose(conn)))

function _build_connectivity_matrix(electrodes::Vector{Symbol}, layout::Layout, cluster_type::Symbol)
    n_electrodes = length(electrodes)

    if cluster_type in (:spatial, :spatiotemporal, :full)
        if isnothing(layout.neighbors)
            @minimal_warning "Layout.neighbors is not set. Computing with default distance criterion (0.35)."
            layout = copy(layout)
            get_neighbors_xy!(layout, 0.35)
        end

        I = Int[]
        J = Int[]
        elec_to_idx = Dict{Symbol,Int}(e => i for (i, e) in enumerate(electrodes))

        for (e_idx, electrode) in enumerate(electrodes)
            if haskey(layout.neighbors, electrode)
                for neighbor in layout.neighbors[electrode].channels
                    n_idx = get(elec_to_idx, neighbor, 0)
                    if n_idx > 0
                        push!(I, e_idx, n_idx)
                        push!(J, n_idx, e_idx)
                    end
                end
            end
        end
        return sparse(I, J, true, n_electrodes, n_electrodes)
    elseif cluster_type in (:temporal, :spectral, :spectrotemporal)
        return sparse(Int[], Int[], Bool[], n_electrodes, n_electrodes)
    else
        @minimal_error(
            "cluster_type must be :spatial, :temporal, :spectral, :spatiotemporal, :spectrotemporal, or :full, got :$cluster_type"
        )
    end
end

# === SPATIAL NEIGHBOR PRE-FILTERING (2D & 3D) ===
function _prefilter_mask_by_neighbors!(mask::BitArray{N}, spatial_connectivity::SparseMatrixCSC{Bool}, min_num_neighbors::Int) where {N}
    min_num_neighbors <= 0 && return nothing

    n_electrodes = size(mask, 1)
    n_cols = div(length(mask), n_electrodes)
    spatial_conn_sym = _symmetrize_connectivity(spatial_connectivity)
    rv = rowvals(spatial_conn_sym)

    @inbounds for col = 0:(n_cols-1)
        offset = col * n_electrodes
        changed = true
        while changed
            changed = false
            for e_idx = 1:n_electrodes
                mask[offset+e_idx] || continue
                neighbor_count = 0
                for ptr in nzrange(spatial_conn_sym, e_idx)
                    n_e_idx = rv[ptr]
                    if e_idx != n_e_idx && mask[offset+n_e_idx]
                        neighbor_count += 1
                    end
                end
                if neighbor_count < min_num_neighbors
                    mask[offset+e_idx] = false
                    changed = true
                end
            end
        end
    end
    return nothing
end

const _prefilter_mask_by_neighbors_tf! = _prefilter_mask_by_neighbors!

function _prefilter_mask_by_neighbors(mask::BitArray{N}, spatial_connectivity::SparseMatrixCSC{Bool}, min_num_neighbors::Int) where {N}
    min_num_neighbors <= 0 && return mask
    filtered_mask = copy(mask)
    _prefilter_mask_by_neighbors!(filtered_mask, spatial_connectivity, min_num_neighbors)
    return filtered_mask
end

# === OBSERVED CLUSTER FINDING & N-D BFS ===
@inline _step_dim(ci::CartesianIndex{N}, d::Int, val::Int) where {N} = CartesianIndex(Base.setindex(ci.I, val, d))

@inline _cluster_axis_flags(::Val{2}, cluster_type::Symbol) =
    (cluster_type in (:spatial, :spatiotemporal, :full), (cluster_type in (:temporal, :spatiotemporal, :full),))

@inline _cluster_axis_flags(::Val{3}, cluster_type::Symbol) = (
    cluster_type in (:spatial, :spatiotemporal, :full),
    (cluster_type in (:spectral, :spectrotemporal, :full), cluster_type in (:temporal, :spatiotemporal, :spectrotemporal, :full)),
)

function _bfs_components(mask::BitArray{N}, spatial_connectivity::SparseMatrixCSC{Bool}, cluster_type::Symbol) where {N}
    components = Vector{CartesianIndex{N}}[]
    any(mask) || return components

    dims = size(mask)
    cluster_labels = zeros(Int, dims)
    rv = rowvals(spatial_connectivity)
    has_spatial, axis_flags = _cluster_axis_flags(Val(N), cluster_type)
    cid = 0

    @inbounds for start_ci in CartesianIndices(mask)
        (!mask[start_ci] || cluster_labels[start_ci] > 0) && continue
        cid += 1
        pixels = CartesianIndex{N}[start_ci]
        cluster_labels[start_ci] = cid
        q_head = 1

        while q_head <= length(pixels)
            ci = pixels[q_head]
            q_head += 1

            if has_spatial
                tail_idx = Base.tail(ci.I)
                for ptr in nzrange(spatial_connectivity, ci[1])
                    nci = CartesianIndex(rv[ptr], tail_idx...)
                    if mask[nci] && cluster_labels[nci] == 0
                        cluster_labels[nci] = cid
                        push!(pixels, nci)
                    end
                end
            end

            for d = 2:N
                if axis_flags[d-1]
                    coord = ci[d]
                    if coord > 1
                        nci = _step_dim(ci, d, coord - 1)
                        if mask[nci] && cluster_labels[nci] == 0
                            cluster_labels[nci] = cid
                            push!(pixels, nci)
                        end
                    end
                    if coord < dims[d]
                        nci = _step_dim(ci, d, coord + 1)
                        if mask[nci] && cluster_labels[nci] == 0
                            cluster_labels[nci] = cid
                            push!(pixels, nci)
                        end
                    end
                end
            end
        end
        push!(components, pixels)
    end
    return components
end

function _seen_indices!(seen::BitVector, pixels::Vector{CartesianIndex{N}}, d::Int) where {N}
    fill!(seen, false)
    @inbounds for ci in pixels
        seen[ci[d]] = true
    end
    return findall(seen)
end

function _find_clusters_connected_components(
    mask::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol;
    polarity::Symbol = :positive,
)
    comps = _bfs_components(mask, spatial_connectivity, cluster_type)
    isempty(comps) && return Cluster[]
    elec_seen = falses(size(mask, 1))
    time_seen = falses(size(mask, 2))
    clusters = Vector{Cluster}(undef, length(comps))
    for (id, pixels) in enumerate(comps)
        t_idx = _seen_indices!(time_seen, pixels, 2)
        e_vec = electrodes[_seen_indices!(elec_seen, pixels, 1)]
        members = Tuple{Int,Int}[ci.I for ci in pixels]
        clusters[id] = Cluster(id, e_vec, t_idx, (time_points[first(t_idx)], time_points[last(t_idx)]), 0.0, 1.0, false, polarity, members)
    end
    return clusters
end

function _find_clusters_connected_components_tf(
    mask::BitArray{3},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol;
    polarity::Symbol = :positive,
)
    comps = _bfs_components(mask, spatial_connectivity, cluster_type)
    isempty(comps) && return TFCluster[]
    elec_seen = falses(size(mask, 1))
    freq_seen = falses(size(mask, 2))
    time_seen = falses(size(mask, 3))
    clusters = Vector{TFCluster}(undef, length(comps))
    for (id, pixels) in enumerate(comps)
        f_idx = _seen_indices!(freq_seen, pixels, 2)
        t_idx = _seen_indices!(time_seen, pixels, 3)
        e_vec = electrodes[_seen_indices!(elec_seen, pixels, 1)]
        f_rng = (frequencies[first(f_idx)], frequencies[last(f_idx)])
        t_rng = (time_points[first(t_idx)], time_points[last(t_idx)])
        clusters[id] = TFCluster(id, e_vec, f_idx, t_idx, f_rng, t_rng, 0.0, 1.0, false, polarity, pixels)
    end
    return clusters
end

function _find_clusters(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
)
    pos =
        isempty(mask_positive) ? Cluster[] :
        _find_clusters_connected_components(
            mask_positive,
            electrodes,
            time_points,
            spatial_connectivity,
            cluster_type;
            polarity = :positive,
        )
    neg =
        isempty(mask_negative) ? Cluster[] :
        _find_clusters_connected_components(
            mask_negative,
            electrodes,
            time_points,
            spatial_connectivity,
            cluster_type;
            polarity = :negative,
        )
    return pos, neg
end

function _find_clusters_tf(
    mask_positive::BitArray{3},
    mask_negative::BitArray{3},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
)
    pos =
        isempty(mask_positive) ? TFCluster[] :
        _find_clusters_connected_components_tf(
            mask_positive,
            electrodes,
            frequencies,
            time_points,
            spatial_connectivity,
            cluster_type;
            polarity = :positive,
        )
    neg =
        isempty(mask_negative) ? TFCluster[] :
        _find_clusters_connected_components_tf(
            mask_negative,
            electrodes,
            frequencies,
            time_points,
            spatial_connectivity,
            cluster_type;
            polarity = :negative,
        )
    return pos, neg
end

# === CLUSTER STATISTICS (MAXSUM) ===
@inline _cluster_points(c::Cluster) = c.members
@inline _cluster_points(c::TFCluster) = c.pixels

@inline _to_ci(idx::Tuple) = CartesianIndex(idx)
@inline _to_ci(ci::CartesianIndex) = ci

@inline _rebuild_cluster(c::Cluster, stat::Float64, p_val::Float64 = c.p_value, is_sig::Bool = c.is_significant) =
    Cluster(c.id, c.electrodes, c.time_indices, c.time_range, stat, p_val, is_sig, c.polarity, c.members)

@inline _rebuild_cluster(c::TFCluster, stat::Float64, p_val::Float64 = c.p_value, is_sig::Bool = c.is_significant) =
    TFCluster(c.id, c.electrodes, c.freq_indices, c.time_indices, c.freq_range, c.time_range, stat, p_val, is_sig, c.polarity, c.pixels)

function _compute_cluster_statistics(
    clusters::Vector{C},
    t_matrix::AbstractArray{<:AbstractFloat};
    return_clusters::Bool = true,
) where {C<:Union{Cluster,TFCluster}}
    isempty(clusters) && return (return_clusters ? (C[], Float64[]) : Float64[])

    n = length(clusters)
    cluster_stats = Vector{Float64}(undef, n)
    updated_clusters = return_clusters ? Vector{C}(undef, n) : nothing

    @inbounds for (i, c) in enumerate(clusters)
        stat = 0.0
        for pt in _cluster_points(c)
            t_val = t_matrix[_to_ci(pt)]
            isfinite(t_val) && (stat += t_val)
        end
        cluster_stats[i] = stat
        if return_clusters
            updated_clusters[i] = _rebuild_cluster(c, stat)
        end
    end

    return return_clusters ? (updated_clusters, cluster_stats) : cluster_stats
end

_compute_cluster_statistics_tf(clusters::Vector{TFCluster}, t_matrix::AbstractArray{<:AbstractFloat,3}, return_clusters::Bool = true) =
    _compute_cluster_statistics(clusters, t_matrix; return_clusters = return_clusters)

_compute_cluster_statistics_tf(clusters::Vector{TFCluster}, t_matrix::AbstractArray{<:AbstractFloat,3}, ::Vector{Symbol}) =
    _compute_cluster_statistics(clusters, t_matrix; return_clusters = true)

# === PERMUTATION CLUSTER STATISTIC (N-D) ===
@inline function _enqueue_nd!(
    ci::CartesianIndex{N},
    cid::Int,
    t_matrix::AbstractArray{<:AbstractFloat,N},
    labels_buf::Array{Int,N},
    queue_buf::AbstractVector{CartesianIndex{N}},
    q_tail::Int,
    cluster_stat::Float64,
) where {N}
    @inbounds begin
        labels_buf[ci] = cid
        val = Float64(t_matrix[ci])
        isfinite(val) && (cluster_stat += val)
        q_tail += 1
        queue_buf[q_tail] = ci
    end
    return q_tail, cluster_stat
end

function _max_cluster_stat(
    mask::BitArray{N},
    t_matrix::AbstractArray{<:AbstractFloat,N},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    labels_buf::Array{Int,N},
    queue_buf::AbstractVector{CartesianIndex{N}};
    is_negative::Bool = false,
    labels_clean::Bool = false,
)::Float64 where {N}
    any(mask) || return 0.0
    labels_clean || fill!(labels_buf, 0)

    dims = size(mask)
    current_cluster_id = 0
    extreme_stat = 0.0
    total_visited = 0

    rv = rowvals(spatial_connectivity)
    has_spatial, axis_flags = _cluster_axis_flags(Val(N), cluster_type)

    @inbounds for start_ci in CartesianIndices(mask)
        (!mask[start_ci] || labels_buf[start_ci] > 0) && continue

        current_cluster_id += 1
        q_head = total_visited + 1
        q_tail, cluster_stat = _enqueue_nd!(start_ci, current_cluster_id, t_matrix, labels_buf, queue_buf, total_visited, 0.0)

        while q_head <= q_tail
            ci = queue_buf[q_head]
            q_head += 1

            if has_spatial
                tail_idx = Base.tail(ci.I)
                for ptr in nzrange(spatial_connectivity, ci[1])
                    nci = CartesianIndex(rv[ptr], tail_idx...)
                    if mask[nci] && labels_buf[nci] == 0
                        q_tail, cluster_stat = _enqueue_nd!(nci, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                    end
                end
            end

            for d = 2:N
                if axis_flags[d-1]
                    coord = ci[d]
                    if coord > 1
                        nci = _step_dim(ci, d, coord - 1)
                        if mask[nci] && labels_buf[nci] == 0
                            q_tail, cluster_stat =
                                _enqueue_nd!(nci, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                        end
                    end
                    if coord < dims[d]
                        nci = _step_dim(ci, d, coord + 1)
                        if mask[nci] && labels_buf[nci] == 0
                            q_tail, cluster_stat =
                                _enqueue_nd!(nci, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                        end
                    end
                end
            end
        end

        total_visited = q_tail
        if is_negative
            cluster_stat < extreme_stat && (extreme_stat = cluster_stat)
        else
            cluster_stat > extreme_stat && (extreme_stat = cluster_stat)
        end
    end

    if labels_clean
        @inbounds for k = 1:total_visited
            labels_buf[queue_buf[k]] = 0
        end
    end

    return extreme_stat
end

_max_cluster_stat(
    mask::BitArray{N},
    t_matrix::AbstractArray{<:AbstractFloat,N},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    labels_buf::Array{Int,N},
    queue_buf::Vector{NTuple{N,Int}};
    kwargs...,
) where {N} =
    _max_cluster_stat(mask, t_matrix, spatial_connectivity, cluster_type, labels_buf, reinterpret(CartesianIndex{N}, queue_buf); kwargs...)

const _max_cluster_stat_tf = _max_cluster_stat

# === THRESHOLD-FREE CLUSTER ENHANCEMENT (TFCE) ===
function _accumulate_tfce_polarity!(
    tfce_map::Matrix{Float64},
    t_map::AbstractMatrix{<:AbstractFloat},
    threshold_max::Float64,
    dh::Float64,
    E::Float64,
    H::Float64,
    sign_mult::Float64,
    spatial_connectivity::SparseMatrixCSC{Bool},
    has_spatial::Bool,
    has_temporal::Bool,
    labels_buf::Matrix{Int},
    queue_buf::Vector{Tuple{Int,Int}},
)
    n_electrodes, n_time = size(t_map)
    rv = rowvals(spatial_connectivity)
    h = dh
    while h <= threshold_max
        fill!(labels_buf, 0)
        h_factor = (h^H) * dh * sign_mult
        cid = 0
        @inbounds for t_idx = 1:n_time, e_idx = 1:n_electrodes
            val = Float64(t_map[e_idx, t_idx]) * sign_mult
            (!(val >= h) || labels_buf[e_idx, t_idx] > 0) && continue

            cid += 1
            labels_buf[e_idx, t_idx] = cid
            queue_buf[1] = (e_idx, t_idx)
            q_head = 1
            q_tail = 1

            while q_head <= q_tail
                curr_e, curr_t = queue_buf[q_head]
                q_head += 1

                if has_spatial
                    for ptr in nzrange(spatial_connectivity, curr_e)
                        ne = rv[ptr]
                        if labels_buf[ne, curr_t] == 0 && Float64(t_map[ne, curr_t]) * sign_mult >= h
                            labels_buf[ne, curr_t] = cid
                            q_tail += 1
                            queue_buf[q_tail] = (ne, curr_t)
                        end
                    end
                end

                if has_temporal
                    if curr_t > 1 && labels_buf[curr_e, curr_t-1] == 0 && Float64(t_map[curr_e, curr_t-1]) * sign_mult >= h
                        labels_buf[curr_e, curr_t-1] = cid
                        q_tail += 1
                        queue_buf[q_tail] = (curr_e, curr_t - 1)
                    end
                    if curr_t < n_time && labels_buf[curr_e, curr_t+1] == 0 && Float64(t_map[curr_e, curr_t+1]) * sign_mult >= h
                        labels_buf[curr_e, curr_t+1] = cid
                        q_tail += 1
                        queue_buf[q_tail] = (curr_e, curr_t + 1)
                    end
                end
            end

            contrib = (Float64(q_tail)^E) * h_factor
            for k = 1:q_tail
                ce, ct = queue_buf[k]
                tfce_map[ce, ct] += contrib
            end
        end
        h += dh
    end
    return tfce_map
end

function _compute_tfce!(
    tfce_map::Matrix{Float64},
    t_map::AbstractMatrix{<:AbstractFloat},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    labels_buf::Matrix{Int},
    queue_buf::Vector{Tuple{Int,Int}};
    E::Real = 0.5,
    H::Real = 2.0,
    dh::Real = 0.1,
)
    dh_f = Float64(dh)
    dh_f > 0.0 || @minimal_error("tfce_dh must be positive, got $dh")
    E_f = Float64(E)
    H_f = Float64(H)
    fill!(tfce_map, 0.0)

    max_pos = 0.0
    max_neg = 0.0
    @inbounds for x in t_map
        v = Float64(x)
        if isfinite(v)
            v > max_pos && (max_pos = v)
            -v > max_neg && (max_neg = -v)
        end
    end

    has_spatial = cluster_type in (:spatial, :spatiotemporal, :full)
    has_temporal = cluster_type in (:temporal, :spatiotemporal, :full)

    if max_pos >= dh_f
        _accumulate_tfce_polarity!(
            tfce_map,
            t_map,
            max_pos,
            dh_f,
            E_f,
            H_f,
            1.0,
            spatial_connectivity,
            has_spatial,
            has_temporal,
            labels_buf,
            queue_buf,
        )
    end
    if max_neg >= dh_f
        _accumulate_tfce_polarity!(
            tfce_map,
            t_map,
            max_neg,
            dh_f,
            E_f,
            H_f,
            -1.0,
            spatial_connectivity,
            has_spatial,
            has_temporal,
            labels_buf,
            queue_buf,
        )
    end
    return tfce_map
end

"""
    _compute_tfce(t_map, electrodes, time_points, spatial_connectivity, cluster_type; E = 0.5, H = 2.0, dh = 0.1)

Compute the Threshold-Free Cluster Enhancement (TFCE) map for a 2D statistic matrix (`electrodes × time`).
Non-finite entries (`NaN`, `Inf`) are ignored and receive a TFCE score of `0.0`.
"""
function _compute_tfce(
    t_map::AbstractMatrix{<:AbstractFloat},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol;
    E::Real = 0.5,
    H::Real = 2.0,
    dh::Real = 0.1,
)
    n_electrodes, n_time = size(t_map)
    tfce_map = zeros(Float64, n_electrodes, n_time)
    labels_buf = Matrix{Int}(undef, n_electrodes, n_time)
    queue_buf = Vector{Tuple{Int,Int}}(undef, n_electrodes * n_time)
    return _compute_tfce!(tfce_map, t_map, spatial_connectivity, cluster_type, labels_buf, queue_buf; E = E, H = H, dh = dh)
end
