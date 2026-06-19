# This file contains clustering functions for finding connected components
"""
    _build_connectivity_matrix(electrodes, layout, cluster_type)

Build connectivity matrix for clustering.

# Arguments
- `electrodes::Vector{Symbol}`: Electrode labels
- `layout::Layout`: Layout with neighbours information
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`

# Returns
- `connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix [electrodes × electrodes]

# Notes
For all cluster types, returns the spatial connectivity matrix.
Temporal connectivity is handled directly in the clustering algorithm.

# Examples
```julia
conn = _build_connectivity_matrix(electrodes, layout, :spatiotemporal)
```
"""
function _build_connectivity_matrix(electrodes::Vector{Symbol}, layout::Layout, cluster_type::Symbol)
    n_electrodes = length(electrodes)

    # Types that use spatial connectivity
    uses_spatial = cluster_type in (:spatial, :spatiotemporal, :full)

    # Types that don't use spatial connectivity (temporal/spectral connections handled in BFS)
    no_spatial = cluster_type in (:temporal, :spectral, :spectrotemporal)

    if uses_spatial
        if isnothing(layout.neighbours)
            @minimal_warning "Layout.neighbours is not set. Computing with default distance criterion (0.25)."
            get_neighbours_xy!(layout, 0.25)
        end

        # Build adjacency matrix
        I = Int[]
        J = Int[]

        for (e_idx, electrode) in enumerate(electrodes)
            # Get neighbours from layout
            if haskey(layout.neighbours, electrode)
                neighbours = layout.neighbours[electrode]
                for neighbour in neighbours.channels
                    # Find index of neighbour in electrodes list
                    n_idx = findfirst(==(neighbour), electrodes)
                    if !isnothing(n_idx)
                        push!(I, e_idx)
                        push!(J, n_idx)
                        # Also add reverse
                        push!(I, n_idx)
                        push!(J, e_idx)
                    end
                end
            end
        end

        # Create sparse matrix and return
        return sparse(I, J, true, n_electrodes, n_electrodes)

    elseif no_spatial
        # No spatial connectivity needed - return empty sparse matrix
        return sparse(Int[], Int[], Bool[], n_electrodes, n_electrodes)

    else
        error("cluster_type must be :spatial, :temporal, :spectral, :spatiotemporal, :spectrotemporal, or :full, got :$cluster_type")
    end
end

# ===================
# PRE-FILTERING
# ===================

"""
    _prefilter_mask_by_neighbors!(mask, spatial_connectivity, min_num_neighbors)

In-place version: Pre-filter mask to remove isolated points 

For each significant point, count how many *neighbouring* significant channels it has
(not counting itself). If a point has fewer than `min_num_neighbors` significant
neighbours, remove it. This is a single-pass filter

# Arguments
- `mask::BitArray{2}`: Significant points mask [electrodes × time]
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix [electrodes × electrodes]
- `min_num_neighbors::Int`: Minimum number of neighboring significant channels required
"""
function _prefilter_mask_by_neighbors!(mask::BitArray{2}, spatial_connectivity::SparseMatrixCSC{Bool}, min_num_neighbors::Int)
    if min_num_neighbors <= 0
        return  # No filtering if min_num_neighbors is 0 or negative
    end

    n_electrodes, n_time = size(mask)

    # Make spatial connectivity symmetric
    spatial_conn_sym = spatial_connectivity .| spatial_connectivity'

    # removing a channel may cause its neighbors to now have too few neighbors
    nremoved = 1
    while nremoved > 0
        nremoved = 0
        for t_idx = 1:n_time
            for e_idx = 1:n_electrodes
                if !mask[e_idx, t_idx]
                    continue
                end

                # Count significant spatial neighbours (NOT counting self)
                neighbor_count = 0
                rv = rowvals(spatial_conn_sym)
                for ptr in nzrange(spatial_conn_sym, e_idx)
                    n_e_idx = rv[ptr]
                    if e_idx != n_e_idx && mask[n_e_idx, t_idx]
                        neighbor_count += 1
                    end
                end

                if neighbor_count < min_num_neighbors
                    mask[e_idx, t_idx] = false
                    nremoved += 1
                end
            end
        end
    end
end

"""
    _prefilter_mask_by_neighbors(mask, spatial_connectivity, min_num_neighbors)

Pre-filter mask to remove isolated points 

Creates a copy of the mask and calls the in-place version.

# Arguments
- `mask::BitArray{2}`: Significant points mask [electrodes × time]
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix [electrodes × electrodes]
- `min_num_neighbors::Int`: Minimum number of neighboring significant channels required

# Returns
- `filtered_mask::BitArray{2}`: Filtered mask with isolated points removed

# Examples
```julia
filtered_mask = _prefilter_mask_by_neighbors(mask, spatial_connectivity, 3)
```
"""
function _prefilter_mask_by_neighbors(mask::BitArray{2}, spatial_connectivity::SparseMatrixCSC{Bool}, min_num_neighbors::Int)
    if min_num_neighbors <= 0
        return mask  # No filtering if min_num_neighbors is 0 or negative
    end

    filtered_mask = copy(mask)
    _prefilter_mask_by_neighbors!(filtered_mask, spatial_connectivity, min_num_neighbors)
    return filtered_mask
end

# ===================
# CLUSTER FINDING (BFS)
# ===================

"""
    _set_cluster_polarity(clusters, polarity)

Create new Cluster objects with the specified polarity, copying all other fields.

# Arguments
- `clusters::Vector{Cluster}`: Clusters to update
- `polarity::Symbol`: `:positive` or `:negative`

# Returns
- `Vector{Cluster}`: New clusters with updated polarity
"""
function _set_cluster_polarity(clusters::Vector{Cluster}, polarity::Symbol)
    polarity in (:positive, :negative) || @minimal_error "polarity must be :positive or :negative"
    return [
        Cluster(c.id, c.electrodes, c.time_indices, c.time_range, c.cluster_stat, c.p_value, c.is_significant, polarity, c.members) for
        c in clusters
    ]
end

"""
    _find_clusters_connected_components(mask, electrodes, time_points, spatial_connectivity, cluster_type)

Find connected clusters in thresholded data using BFS.

# Arguments
- `mask::BitArray{2}`: Significant points mask [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`

# Returns
- `clusters::Vector{Cluster}`: Found clusters (without statistics)

# Examples
```julia
clusters = _find_clusters_connected_components(mask, electrodes, time_points, conn, :spatiotemporal)
```
"""
function _find_clusters_connected_components(
    mask::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
)
    n_electrodes, n_time = size(mask)
    clusters = Cluster[]

    if !any(mask)
        # No significant points
        return clusters
    end

    # Label matrix to track which cluster each point belongs to
    cluster_labels = zeros(Int, n_electrodes, n_time)
    current_cluster_id = 0


    # BFS to find connected components
    for e_idx = 1:n_electrodes
        for t_idx = 1:n_time
            # Skip if not significant or already labeled
            if !mask[e_idx, t_idx] || cluster_labels[e_idx, t_idx] > 0
                continue
            end

            # Start new cluster
            current_cluster_id += 1
            cluster_electrodes = Set{Symbol}()
            cluster_time_indices = Set{Int}()
            cluster_members = Tuple{Int,Int}[]  # exact (electrode_idx, time_idx) pairs

            # BFS queue
            queue = [(e_idx, t_idx)]
            queue_head = 1
            cluster_labels[e_idx, t_idx] = current_cluster_id
            push!(cluster_electrodes, electrodes[e_idx])
            push!(cluster_time_indices, t_idx)
            push!(cluster_members, (e_idx, t_idx))

            # BFS traversal
            while queue_head <= length(queue)
                current_e, current_t = queue[queue_head]
                queue_head += 1

                # Spatial connectivity
                if cluster_type in (:spatial, :spatiotemporal)
                    # EegFun's spatial connectivity is symmetric, so we can iterate rows
                    rv = rowvals(spatial_connectivity)
                    for ptr in nzrange(spatial_connectivity, current_e)
                        n_e = rv[ptr]
                        if mask[n_e, current_t] && cluster_labels[n_e, current_t] == 0
                            cluster_labels[n_e, current_t] = current_cluster_id
                            push!(queue, (n_e, current_t))
                            push!(cluster_electrodes, electrodes[n_e])
                            push!(cluster_time_indices, current_t)
                            push!(cluster_members, (n_e, current_t))
                        end
                    end
                end

                # Temporal connectivity
                if cluster_type in (:temporal, :spatiotemporal)
                    if current_t > 1 && mask[current_e, current_t-1] && cluster_labels[current_e, current_t-1] == 0
                        cluster_labels[current_e, current_t-1] = current_cluster_id
                        push!(queue, (current_e, current_t-1))
                        push!(cluster_electrodes, electrodes[current_e])
                        push!(cluster_time_indices, current_t-1)
                        push!(cluster_members, (current_e, current_t-1))
                    end
                    if current_t < n_time && mask[current_e, current_t+1] && cluster_labels[current_e, current_t+1] == 0
                        cluster_labels[current_e, current_t+1] = current_cluster_id
                        push!(queue, (current_e, current_t+1))
                        push!(cluster_electrodes, electrodes[current_e])
                        push!(cluster_time_indices, current_t+1)
                        push!(cluster_members, (current_e, current_t+1))
                    end
                end
            end

            # Convert Set to sorted Vector efficiently
            time_indices_vec = sort([t for t in cluster_time_indices])
            time_range = (time_points[time_indices_vec[1]], time_points[time_indices_vec[end]])

            # Convert Set to Vector efficiently
            electrodes_vec = [e for e in cluster_electrodes]

            cluster = Cluster(
                current_cluster_id,
                electrodes_vec,
                time_indices_vec,
                time_range,
                0.0,  # cluster_stat - will be computed later
                1.0,  # p_value - will be computed later
                false,  # is_significant - will be set by caller (temporary)
                :positive,  # polarity - will be set by caller (temporary)
                cluster_members,
            )
            push!(clusters, cluster)
        end
    end

    return clusters
end

"""
    find_clusters(mask_positive, mask_negative, electrodes, time_points, spatial_connectivity, cluster_type)

Find positive and negative clusters.

# Arguments
- `mask_positive::BitArray{2}`: Positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Negative significant points [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`

# Returns
- `positive_clusters::Vector{Cluster}`: Positive clusters
- `negative_clusters::Vector{Cluster}`: Negative clusters

# Examples
```julia
pos_clusters, neg_clusters = _find_clusters(mask_pos, mask_neg, electrodes, time_points, conn, :spatiotemporal)
```
"""
function _find_clusters(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
)
    positive_clusters = Cluster[]
    negative_clusters = Cluster[]

    if !isempty(mask_positive)
        positive_clusters_raw =
            _find_clusters_connected_components(mask_positive, electrodes, time_points, spatial_connectivity, cluster_type)
        positive_clusters = _set_cluster_polarity(positive_clusters_raw, :positive)
    end

    if !isempty(mask_negative)
        negative_clusters_raw =
            _find_clusters_connected_components(mask_negative, electrodes, time_points, spatial_connectivity, cluster_type)
        negative_clusters = _set_cluster_polarity(negative_clusters_raw, :negative)
    end
    return positive_clusters, negative_clusters
end

# ===================
# CLUSTER STATISTICS
# ===================

"""
    _compute_cluster_statistics(clusters, t_matrix, electrode_to_idx; return_clusters)

Compute cluster-level statistics (maxsum).

Consolidated version that can either return full Cluster objects or just statistics.

# Arguments
- `clusters::Vector{Cluster}`: Clusters to compute statistics for
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `electrode_to_idx::Dict{Symbol,Int}`: Electrode label to index mapping
- `return_clusters::Bool`: If true, return (clusters, stats); if false, return just stats

# Returns
- If `return_clusters=true`: `(updated_clusters::Vector{Cluster}, cluster_stats::Vector{Float64})`
- If `return_clusters=false`: `cluster_stats::Vector{Float64}`
"""
function _compute_cluster_statistics(
    clusters::Vector{Cluster},
    t_matrix::Array{Float64,2},
    electrode_to_idx::Dict{Symbol,Int};
    return_clusters::Bool = true,
)
    if isempty(clusters)
        return return_clusters ? (Cluster[], Float64[]) : Float64[]
    end

    cluster_stats = Float64[]
    sizehint!(cluster_stats, length(clusters))

    updated_clusters = return_clusters ? Cluster[] : nothing

    for cluster in clusters
        cluster_stat = 0.0

        # Sum of t-values only at exact member points 
        for (e_idx, t_idx) in cluster.members
            t_val = t_matrix[e_idx, t_idx]
            if !isnan(t_val) && !isinf(t_val)
                cluster_stat += t_val
            end
        end

        push!(cluster_stats, cluster_stat)

        if return_clusters
            # Update cluster with statistic
            updated_cluster = Cluster(
                cluster.id,
                cluster.electrodes,
                cluster.time_indices,
                cluster.time_range,
                cluster_stat,
                cluster.p_value,
                cluster.is_significant,
                cluster.polarity,
                cluster.members,
            )
            push!(updated_clusters, updated_cluster)
        end
    end

    return return_clusters ? (updated_clusters, cluster_stats) : cluster_stats
end

"""
    _compute_cluster_statistics(clusters, t_matrix, electrodes)

Compute cluster-level statistics (maxsum).
"""
function _compute_cluster_statistics(clusters::Vector{Cluster}, t_matrix::Array{Float64,2}, electrodes::Vector{Symbol})
    # Create electrode index lookup
    electrode_to_idx = Dict(e => i for (i, e) in enumerate(electrodes))
    return _compute_cluster_statistics(clusters, t_matrix, electrode_to_idx, return_clusters = true)
end


# ===================
# TF PRE-FILTERING (3D)
# ===================

"""
    _prefilter_mask_by_neighbors_tf!(mask, spatial_connectivity, min_num_neighbors)

In-place pre-filtering of 3D mask to remove spatially isolated points.
"""
function _prefilter_mask_by_neighbors_tf!(mask::BitArray{3}, spatial_connectivity::SparseMatrixCSC{Bool}, min_num_neighbors::Int)
    if min_num_neighbors <= 0
        return
    end

    n_electrodes, n_freqs, n_time = size(mask)
    spatial_conn_sym = spatial_connectivity .| spatial_connectivity'

    nremoved = 1
    while nremoved > 0
        nremoved = 0
        for t_idx = 1:n_time
            for f_idx = 1:n_freqs
                for e_idx = 1:n_electrodes
                    if !mask[e_idx, f_idx, t_idx]
                        continue
                    end
                    # Count significant spatial neighbours (NOT counting self)
                    neighbor_count = 0
                    rv = rowvals(spatial_conn_sym)
                    for ptr in nzrange(spatial_conn_sym, e_idx)
                        n_e_idx = rv[ptr]
                        if e_idx != n_e_idx && mask[n_e_idx, f_idx, t_idx]
                            neighbor_count += 1
                        end
                    end
                    if neighbor_count < min_num_neighbors
                        mask[e_idx, f_idx, t_idx] = false
                        nremoved += 1
                    end
                end
            end
        end
    end
end


# ===================
# TF CLUSTER FINDING (3D BFS)
# ===================

"""
    _find_clusters_connected_components_tf(mask, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)

Find connected clusters in 3D thresholded TF data using BFS.

Connectivity is determined by `cluster_type`:
- `:temporal`: adjacent time points only
- `:spatial`: spatially adjacent electrodes only (same freq and time)
- `:spectral`: adjacent frequency bins only (same electrode and time)
- `:spatiotemporal`: spatial + temporal neighbours
- `:spectrotemporal`: spectral + temporal neighbours
- `:full`: spatial + spectral + temporal neighbours (all 3 dimensions)
"""
function _find_clusters_connected_components_tf(
    mask::BitArray{3},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
)
    n_electrodes, n_freqs, n_time = size(mask)
    clusters = TFCluster[]

    if !any(mask)
        return clusters
    end

    cluster_labels = zeros(Int, n_electrodes, n_freqs, n_time)
    current_cluster_id = 0

    # Determine which connectivity dimensions are active
    use_spatial = cluster_type in (:spatial, :spatiotemporal, :full)
    use_spectral = cluster_type in (:spectral, :spectrotemporal, :full)
    use_temporal = cluster_type in (:temporal, :spatiotemporal, :spectrotemporal, :full)



    # BFS to find connected components
    for e_idx = 1:n_electrodes
        for f_idx = 1:n_freqs
            for t_idx = 1:n_time
                if !mask[e_idx, f_idx, t_idx] || cluster_labels[e_idx, f_idx, t_idx] > 0
                    continue
                end

                # Start new cluster
                current_cluster_id += 1
                cluster_electrodes = Set{Symbol}()
                cluster_freq_indices = Set{Int}()
                cluster_time_indices = Set{Int}()
                cluster_pixels = CartesianIndex{3}[]
                cluster_members = Tuple{Int,Int,Int}[]

                # BFS queue
                queue = [(e_idx, f_idx, t_idx)]
                queue_head = 1
                cluster_labels[e_idx, f_idx, t_idx] = current_cluster_id
                push!(cluster_electrodes, electrodes[e_idx])
                push!(cluster_freq_indices, f_idx)
                push!(cluster_time_indices, t_idx)
                push!(cluster_pixels, CartesianIndex(e_idx, f_idx, t_idx))
                push!(cluster_members, (e_idx, f_idx, t_idx))

                # BFS traversal
                while queue_head <= length(queue)
                    current_e, current_f, current_t = queue[queue_head]
                    queue_head += 1

                    # Spatial connectivity (same freq, same time)
                    if use_spatial
                        rv = rowvals(spatial_connectivity)
                        for ptr in nzrange(spatial_connectivity, current_e)
                            n_e = rv[ptr]
                            if mask[n_e, current_f, current_t] && cluster_labels[n_e, current_f, current_t] == 0
                                cluster_labels[n_e, current_f, current_t] = current_cluster_id
                                push!(queue, (n_e, current_f, current_t))
                                push!(cluster_electrodes, electrodes[n_e])
                                push!(cluster_freq_indices, current_f)
                                push!(cluster_time_indices, current_t)
                                push!(cluster_pixels, CartesianIndex(n_e, current_f, current_t))
                                push!(cluster_members, (n_e, current_f, current_t))
                            end
                        end
                    end

                    # Spectral connectivity (same electrode, same time, adjacent freq)
                    if use_spectral
                        if current_f > 1 &&
                           mask[current_e, current_f-1, current_t] &&
                           cluster_labels[current_e, current_f-1, current_t] == 0
                            cluster_labels[current_e, current_f-1, current_t] = current_cluster_id
                            push!(queue, (current_e, current_f-1, current_t))
                            push!(cluster_electrodes, electrodes[current_e])
                            push!(cluster_freq_indices, current_f-1)
                            push!(cluster_time_indices, current_t)
                            push!(cluster_pixels, CartesianIndex(current_e, current_f-1, current_t))
                            push!(cluster_members, (current_e, current_f-1, current_t))
                        end
                        if current_f < n_freqs &&
                           mask[current_e, current_f+1, current_t] &&
                           cluster_labels[current_e, current_f+1, current_t] == 0
                            cluster_labels[current_e, current_f+1, current_t] = current_cluster_id
                            push!(queue, (current_e, current_f+1, current_t))
                            push!(cluster_electrodes, electrodes[current_e])
                            push!(cluster_freq_indices, current_f+1)
                            push!(cluster_time_indices, current_t)
                            push!(cluster_pixels, CartesianIndex(current_e, current_f+1, current_t))
                            push!(cluster_members, (current_e, current_f+1, current_t))
                        end
                    end

                    # Temporal connectivity (same electrode, same freq, adjacent time)
                    if use_temporal
                        if current_t > 1 &&
                           mask[current_e, current_f, current_t-1] &&
                           cluster_labels[current_e, current_f, current_t-1] == 0
                            cluster_labels[current_e, current_f, current_t-1] = current_cluster_id
                            push!(queue, (current_e, current_f, current_t-1))
                            push!(cluster_electrodes, electrodes[current_e])
                            push!(cluster_freq_indices, current_f)
                            push!(cluster_time_indices, current_t-1)
                            push!(cluster_pixels, CartesianIndex(current_e, current_f, current_t-1))
                            push!(cluster_members, (current_e, current_f, current_t-1))
                        end
                        if current_t < n_time &&
                           mask[current_e, current_f, current_t+1] &&
                           cluster_labels[current_e, current_f, current_t+1] == 0
                            cluster_labels[current_e, current_f, current_t+1] = current_cluster_id
                            push!(queue, (current_e, current_f, current_t+1))
                            push!(cluster_electrodes, electrodes[current_e])
                            push!(cluster_freq_indices, current_f)
                            push!(cluster_time_indices, current_t+1)
                            push!(cluster_pixels, CartesianIndex(current_e, current_f, current_t+1))
                            push!(cluster_members, (current_e, current_f, current_t+1))
                        end
                    end
                end

                freq_indices_vec = sort(collect(cluster_freq_indices))
                time_indices_vec = sort(collect(cluster_time_indices))
                freq_range = (frequencies[freq_indices_vec[1]], frequencies[freq_indices_vec[end]])
                time_range = (time_points[time_indices_vec[1]], time_points[time_indices_vec[end]])
                electrodes_vec = collect(cluster_electrodes)

                cluster = TFCluster(
                    current_cluster_id,
                    electrodes_vec,
                    freq_indices_vec,
                    time_indices_vec,
                    freq_range,
                    time_range,
                    0.0,  # cluster_stat - computed later
                    1.0,  # p_value - computed later
                    false,  # is_significant
                    :positive,  # polarity - set by caller
                    cluster_pixels,
                )
                push!(clusters, cluster)
            end
        end
    end

    return clusters
end


"""
    _set_cluster_polarity_tf(clusters, polarity)

Create new TFCluster objects with the specified polarity.
"""
function _set_cluster_polarity_tf(clusters::Vector{TFCluster}, polarity::Symbol)
    polarity in (:positive, :negative) || @minimal_error "polarity must be :positive or :negative"
    return [
        TFCluster(
            c.id,
            c.electrodes,
            c.freq_indices,
            c.time_indices,
            c.freq_range,
            c.time_range,
            c.cluster_stat,
            c.p_value,
            c.is_significant,
            polarity,
            c.pixels,
        ) for c in clusters
    ]
end


"""
    _find_clusters_tf(mask_positive, mask_negative, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)

Find positive and negative TF clusters.
"""
function _find_clusters_tf(
    mask_positive::BitArray{3},
    mask_negative::BitArray{3},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
)
    positive_clusters = TFCluster[]
    negative_clusters = TFCluster[]

    if !isempty(mask_positive)
        raw =
            _find_clusters_connected_components_tf(mask_positive, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)
        positive_clusters = _set_cluster_polarity_tf(raw, :positive)
    end

    if !isempty(mask_negative)
        raw =
            _find_clusters_connected_components_tf(mask_negative, electrodes, frequencies, time_points, spatial_connectivity, cluster_type)
        negative_clusters = _set_cluster_polarity_tf(raw, :negative)
    end

    return positive_clusters, negative_clusters
end


# ===================
# TF CLUSTER STATISTICS (3D)
# ===================

"""
    _compute_cluster_statistics_tf(clusters, t_matrix, electrode_to_idx; return_clusters)

Compute cluster-level statistics (maxsum) for TF clusters.
Sums t-values at each cluster's exact pixel locations.
"""
function _compute_cluster_statistics_tf(clusters::Vector{TFCluster}, t_matrix::Array{Float64,3}, return_clusters::Bool = true)
    if isempty(clusters)
        return return_clusters ? (TFCluster[], Float64[]) : Float64[]
    end

    cluster_stats = Float64[]
    sizehint!(cluster_stats, length(clusters))
    updated_clusters = return_clusters ? TFCluster[] : nothing

    for cluster in clusters
        cluster_stat = 0.0
        for pixel in cluster.pixels
            t_val = t_matrix[pixel]
            if !isnan(t_val) && !isinf(t_val)
                cluster_stat += t_val
            end
        end

        push!(cluster_stats, cluster_stat)

        if return_clusters
            updated_cluster = TFCluster(
                cluster.id,
                cluster.electrodes,
                cluster.freq_indices,
                cluster.time_indices,
                cluster.freq_range,
                cluster.time_range,
                cluster_stat,
                cluster.p_value,
                cluster.is_significant,
                cluster.polarity,
                cluster.pixels,
            )
            push!(updated_clusters, updated_cluster)
        end
    end

    return return_clusters ? (updated_clusters, cluster_stats) : cluster_stats
end

function _compute_cluster_statistics_tf(clusters::Vector{TFCluster}, t_matrix::Array{Float64,3}, electrodes::Vector{Symbol})
    return _compute_cluster_statistics_tf(clusters, t_matrix, true)
end
