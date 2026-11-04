# Connected-component clustering, neighbor pre-filtering, cluster statistics, and TFCE.

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

# === OBSERVED CLUSTER FINDING (2D & 3D BFS) ===

function _find_clusters_connected_components(
    mask::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol;
    polarity::Symbol = :positive,
)
    n_electrodes, n_time = size(mask)
    clusters = Cluster[]
    any(mask) || return clusters

    cluster_labels = zeros(Int, n_electrodes, n_time)
    elec_seen = falses(n_electrodes)
    time_seen = falses(n_time)
    rv = rowvals(spatial_connectivity)
    use_spatial = cluster_type in (:spatial, :spatiotemporal)
    use_temporal = cluster_type in (:temporal, :spatiotemporal)
    current_cluster_id = 0

    @inbounds for e_idx = 1:n_electrodes, t_idx = 1:n_time
        (!mask[e_idx, t_idx] || cluster_labels[e_idx, t_idx] > 0) && continue

        current_cluster_id += 1
        fill!(elec_seen, false)
        fill!(time_seen, false)

        members = Tuple{Int,Int}[(e_idx, t_idx)]
        cluster_labels[e_idx, t_idx] = current_cluster_id
        elec_seen[e_idx] = true
        time_seen[t_idx] = true
        q_head = 1

        while q_head <= length(members)
            ce, ct = members[q_head]
            q_head += 1

            if use_spatial
                for ptr in nzrange(spatial_connectivity, ce)
                    ne = rv[ptr]
                    if mask[ne, ct] && cluster_labels[ne, ct] == 0
                        cluster_labels[ne, ct] = current_cluster_id
                        elec_seen[ne] = true
                        push!(members, (ne, ct))
                    end
                end
            end

            if use_temporal
                if ct > 1 && mask[ce, ct-1] && cluster_labels[ce, ct-1] == 0
                    cluster_labels[ce, ct-1] = current_cluster_id
                    time_seen[ct-1] = true
                    push!(members, (ce, ct-1))
                end
                if ct < n_time && mask[ce, ct+1] && cluster_labels[ce, ct+1] == 0
                    cluster_labels[ce, ct+1] = current_cluster_id
                    time_seen[ct+1] = true
                    push!(members, (ce, ct+1))
                end
            end
        end

        time_indices_vec = findall(time_seen)
        time_range = (time_points[first(time_indices_vec)], time_points[last(time_indices_vec)])
        electrodes_vec = electrodes[findall(elec_seen)]

        push!(clusters, Cluster(current_cluster_id, electrodes_vec, time_indices_vec, time_range, 0.0, 1.0, false, polarity, members))
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

function _find_clusters_connected_components_tf(
    mask::BitArray{3},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol;
    polarity::Symbol = :positive,
)
    n_electrodes, n_freqs, n_time = size(mask)
    clusters = TFCluster[]
    any(mask) || return clusters

    cluster_labels = zeros(Int, n_electrodes, n_freqs, n_time)
    elec_seen = falses(n_electrodes)
    freq_seen = falses(n_freqs)
    time_seen = falses(n_time)
    rv = rowvals(spatial_connectivity)
    use_spatial = cluster_type in (:spatial, :spatiotemporal, :full)
    use_spectral = cluster_type in (:spectral, :spectrotemporal, :full)
    use_temporal = cluster_type in (:temporal, :spatiotemporal, :spectrotemporal, :full)
    current_cluster_id = 0

    @inbounds for e_idx = 1:n_electrodes, f_idx = 1:n_freqs, t_idx = 1:n_time
        (!mask[e_idx, f_idx, t_idx] || cluster_labels[e_idx, f_idx, t_idx] > 0) && continue

        current_cluster_id += 1
        fill!(elec_seen, false)
        fill!(freq_seen, false)
        fill!(time_seen, false)

        pixels = CartesianIndex{3}[CartesianIndex(e_idx, f_idx, t_idx)]
        cluster_labels[e_idx, f_idx, t_idx] = current_cluster_id
        elec_seen[e_idx] = true
        freq_seen[f_idx] = true
        time_seen[t_idx] = true
        q_head = 1

        while q_head <= length(pixels)
            ci = pixels[q_head]
            ce, cf, ct = ci[1], ci[2], ci[3]
            q_head += 1

            if use_spatial
                for ptr in nzrange(spatial_connectivity, ce)
                    ne = rv[ptr]
                    if mask[ne, cf, ct] && cluster_labels[ne, cf, ct] == 0
                        cluster_labels[ne, cf, ct] = current_cluster_id
                        elec_seen[ne] = true
                        push!(pixels, CartesianIndex(ne, cf, ct))
                    end
                end
            end

            if use_spectral
                if cf > 1 && mask[ce, cf-1, ct] && cluster_labels[ce, cf-1, ct] == 0
                    cluster_labels[ce, cf-1, ct] = current_cluster_id
                    freq_seen[cf-1] = true
                    push!(pixels, CartesianIndex(ce, cf-1, ct))
                end
                if cf < n_freqs && mask[ce, cf+1, ct] && cluster_labels[ce, cf+1, ct] == 0
                    cluster_labels[ce, cf+1, ct] = current_cluster_id
                    freq_seen[cf+1] = true
                    push!(pixels, CartesianIndex(ce, cf+1, ct))
                end
            end

            if use_temporal
                if ct > 1 && mask[ce, cf, ct-1] && cluster_labels[ce, cf, ct-1] == 0
                    cluster_labels[ce, cf, ct-1] = current_cluster_id
                    time_seen[ct-1] = true
                    push!(pixels, CartesianIndex(ce, cf, ct-1))
                end
                if ct < n_time && mask[ce, cf, ct+1] && cluster_labels[ce, cf, ct+1] == 0
                    cluster_labels[ce, cf, ct+1] = current_cluster_id
                    time_seen[ct+1] = true
                    push!(pixels, CartesianIndex(ce, cf, ct+1))
                end
            end
        end

        freq_indices_vec = findall(freq_seen)
        time_indices_vec = findall(time_seen)
        freq_range = (frequencies[first(freq_indices_vec)], frequencies[last(freq_indices_vec)])
        time_range = (time_points[first(time_indices_vec)], time_points[last(time_indices_vec)])
        electrodes_vec = electrodes[findall(elec_seen)]

        push!(
            clusters,
            TFCluster(
                current_cluster_id,
                electrodes_vec,
                freq_indices_vec,
                time_indices_vec,
                freq_range,
                time_range,
                0.0,
                1.0,
                false,
                polarity,
                pixels,
            ),
        )
    end

    return clusters
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

function _compute_cluster_statistics(clusters::Vector{Cluster}, t_matrix::AbstractMatrix{<:AbstractFloat}; return_clusters::Bool = true)
    isempty(clusters) && return (return_clusters ? (Cluster[], Float64[]) : Float64[])

    n = length(clusters)
    cluster_stats = Vector{Float64}(undef, n)
    updated_clusters = return_clusters ? Vector{Cluster}(undef, n) : nothing

    @inbounds for (i, c) in enumerate(clusters)
        stat = 0.0
        for (e_idx, t_idx) in c.members
            t_val = t_matrix[e_idx, t_idx]
            isfinite(t_val) && (stat += t_val)
        end
        cluster_stats[i] = stat
        if return_clusters
            updated_clusters[i] =
                Cluster(c.id, c.electrodes, c.time_indices, c.time_range, stat, c.p_value, c.is_significant, c.polarity, c.members)
        end
    end

    return return_clusters ? (updated_clusters, cluster_stats) : cluster_stats
end

function _compute_cluster_statistics_tf(
    clusters::Vector{TFCluster},
    t_matrix::AbstractArray{<:AbstractFloat,3},
    return_clusters::Bool = true,
)
    isempty(clusters) && return (return_clusters ? (TFCluster[], Float64[]) : Float64[])

    n = length(clusters)
    cluster_stats = Vector{Float64}(undef, n)
    updated_clusters = return_clusters ? Vector{TFCluster}(undef, n) : nothing

    @inbounds for (i, c) in enumerate(clusters)
        stat = 0.0
        for pixel in c.pixels
            t_val = t_matrix[pixel]
            isfinite(t_val) && (stat += t_val)
        end
        cluster_stats[i] = stat
        if return_clusters
            updated_clusters[i] = TFCluster(
                c.id,
                c.electrodes,
                c.freq_indices,
                c.time_indices,
                c.freq_range,
                c.time_range,
                stat,
                c.p_value,
                c.is_significant,
                c.polarity,
                c.pixels,
            )
        end
    end

    return return_clusters ? (updated_clusters, cluster_stats) : cluster_stats
end

_compute_cluster_statistics_tf(clusters::Vector{TFCluster}, t_matrix::AbstractArray{<:AbstractFloat,3}, ::Vector{Symbol}) =
    _compute_cluster_statistics_tf(clusters, t_matrix, true)

# === ZERO-ALLOCATION PERMUTATION EXTREME CLUSTER STATISTIC ===

@inline function _enqueue_2d!(
    ne::Int,
    nt::Int,
    cid::Int,
    t_matrix::AbstractMatrix{<:AbstractFloat},
    labels_buf::Matrix{Int},
    queue_buf::Vector{Tuple{Int,Int}},
    q_tail::Int,
    cluster_stat::Float64,
)
    @inbounds begin
        labels_buf[ne, nt] = cid
        val = Float64(t_matrix[ne, nt])
        if isfinite(val)
            cluster_stat += val
        end
        q_tail += 1
        queue_buf[q_tail] = (ne, nt)
    end
    return q_tail, cluster_stat
end

function _max_cluster_stat(
    mask::BitArray{2},
    t_matrix::AbstractMatrix{<:AbstractFloat},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    labels_buf::Matrix{Int},
    queue_buf::Vector{Tuple{Int,Int}};
    is_negative::Bool = false,
    labels_clean::Bool = false,
)::Float64
    any(mask) || return 0.0
    labels_clean || fill!(labels_buf, 0)

    n_electrodes, n_time = size(mask)
    current_cluster_id = 0
    extreme_stat = 0.0
    total_visited = 0

    rv = rowvals(spatial_connectivity)
    has_spatial = cluster_type in (:spatial, :spatiotemporal, :full)
    has_temporal = cluster_type in (:temporal, :spatiotemporal, :full)

    @inbounds for t_idx = 1:n_time, e_idx = 1:n_electrodes
        (!mask[e_idx, t_idx] || labels_buf[e_idx, t_idx] > 0) && continue

        current_cluster_id += 1
        q_head = total_visited + 1
        q_tail, cluster_stat = _enqueue_2d!(e_idx, t_idx, current_cluster_id, t_matrix, labels_buf, queue_buf, total_visited, 0.0)

        while q_head <= q_tail
            curr_e, curr_t = queue_buf[q_head]
            q_head += 1

            if has_spatial
                for ptr in nzrange(spatial_connectivity, curr_e)
                    ne = rv[ptr]
                    if mask[ne, curr_t] && labels_buf[ne, curr_t] == 0
                        q_tail, cluster_stat =
                            _enqueue_2d!(ne, curr_t, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                    end
                end
            end

            if has_temporal
                if curr_t > 1 && mask[curr_e, curr_t-1] && labels_buf[curr_e, curr_t-1] == 0
                    q_tail, cluster_stat =
                        _enqueue_2d!(curr_e, curr_t - 1, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                end
                if curr_t < n_time && mask[curr_e, curr_t+1] && labels_buf[curr_e, curr_t+1] == 0
                    q_tail, cluster_stat =
                        _enqueue_2d!(curr_e, curr_t + 1, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
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
            e, t = queue_buf[k]
            labels_buf[e, t] = 0
        end
    end

    return extreme_stat
end

@inline function _enqueue_3d!(
    ne::Int,
    nf::Int,
    nt::Int,
    cid::Int,
    t_matrix::AbstractArray{<:AbstractFloat,3},
    labels_buf::Array{Int,3},
    queue_buf::Vector{Tuple{Int,Int,Int}},
    q_tail::Int,
    cluster_stat::Float64,
)
    @inbounds begin
        labels_buf[ne, nf, nt] = cid
        val = Float64(t_matrix[ne, nf, nt])
        if isfinite(val)
            cluster_stat += val
        end
        q_tail += 1
        queue_buf[q_tail] = (ne, nf, nt)
    end
    return q_tail, cluster_stat
end

function _max_cluster_stat_tf(
    mask::BitArray{3},
    t_matrix::AbstractArray{<:AbstractFloat,3},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    labels_buf::Array{Int,3},
    queue_buf::Vector{Tuple{Int,Int,Int}};
    is_negative::Bool = false,
    labels_clean::Bool = false,
)::Float64
    any(mask) || return 0.0
    labels_clean || fill!(labels_buf, 0)

    n_electrodes, n_freqs, n_time = size(mask)
    current_cluster_id = 0
    extreme_stat = 0.0
    total_visited = 0

    rv = rowvals(spatial_connectivity)
    use_spatial = cluster_type in (:spatial, :spatiotemporal, :full)
    use_spectral = cluster_type in (:spectral, :spectrotemporal, :full)
    use_temporal = cluster_type in (:temporal, :spatiotemporal, :spectrotemporal, :full)

    @inbounds for t_idx = 1:n_time, f_idx = 1:n_freqs, e_idx = 1:n_electrodes
        (!mask[e_idx, f_idx, t_idx] || labels_buf[e_idx, f_idx, t_idx] > 0) && continue

        current_cluster_id += 1
        q_head = total_visited + 1
        q_tail, cluster_stat = _enqueue_3d!(e_idx, f_idx, t_idx, current_cluster_id, t_matrix, labels_buf, queue_buf, total_visited, 0.0)

        while q_head <= q_tail
            curr_e, curr_f, curr_t = queue_buf[q_head]
            q_head += 1

            if use_spatial
                for ptr in nzrange(spatial_connectivity, curr_e)
                    ne = rv[ptr]
                    if mask[ne, curr_f, curr_t] && labels_buf[ne, curr_f, curr_t] == 0
                        q_tail, cluster_stat =
                            _enqueue_3d!(ne, curr_f, curr_t, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                    end
                end
            end

            if use_spectral
                if curr_f > 1 && mask[curr_e, curr_f-1, curr_t] && labels_buf[curr_e, curr_f-1, curr_t] == 0
                    q_tail, cluster_stat =
                        _enqueue_3d!(curr_e, curr_f - 1, curr_t, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                end
                if curr_f < n_freqs && mask[curr_e, curr_f+1, curr_t] && labels_buf[curr_e, curr_f+1, curr_t] == 0
                    q_tail, cluster_stat =
                        _enqueue_3d!(curr_e, curr_f + 1, curr_t, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                end
            end

            if use_temporal
                if curr_t > 1 && mask[curr_e, curr_f, curr_t-1] && labels_buf[curr_e, curr_f, curr_t-1] == 0
                    q_tail, cluster_stat =
                        _enqueue_3d!(curr_e, curr_f, curr_t - 1, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
                end
                if curr_t < n_time && mask[curr_e, curr_f, curr_t+1] && labels_buf[curr_e, curr_f, curr_t+1] == 0
                    q_tail, cluster_stat =
                        _enqueue_3d!(curr_e, curr_f, curr_t + 1, current_cluster_id, t_matrix, labels_buf, queue_buf, q_tail, cluster_stat)
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
            e, f, t = queue_buf[k]
            labels_buf[e, f, t] = 0
        end
    end

    return extreme_stat
end
