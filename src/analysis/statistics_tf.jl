"""
Time-Frequency (3D) Statistics Functions.
Adapts the permutation framework for [Electrodes × Frequencies × Time] data.
"""

using SparseArrays
using Statistics
using Distributions
using Random
using ProgressMeter: Progress, next!

"""
    compute_t_map_3d(data1::Array{Float64, 4}, data2::Array{Float64, 4}; tail::Symbol=:both)

Compute 3D T-statistic map [Electrodes x Frequencies x Time] from 4D data [Participants x Elec x Freq x Time].
Only supports paired t-test for now (data1 - data2).
"""
function compute_t_map_3d(data1::Array{Float64, 4}, data2::Array{Float64, 4}; tail::Symbol=:both)
    n_part, n_elec, n_freq, n_time = size(data1)
    t_map = Array{Float64}(undef, n_elec, n_freq, n_time)
    
    # Paired T-Test
    diffs = data1 .- data2
    
    # Mean and Std over participants (dim 1)
    mean_diff = dropdims(mean(diffs, dims=1), dims=1)
    std_diff = dropdims(std(diffs, dims=1), dims=1)
    
    sqrt_n = sqrt(n_part)
    
    @inbounds for i in eachindex(t_map)
        if std_diff[i] > 1e-12
            t_map[i] = mean_diff[i] / (std_diff[i] / sqrt_n)
        else
            t_map[i] = 0.0
        end
    end
    
    return t_map
end

"""
    cluster_permutation_test_tf(data1::Array{Float64, 4}, 
                                data2::Array{Float64, 4},
                                spatial_conn::SparseMatrixCSC{Bool},
                                electrodes::Vector{Symbol},
                                freqs::Vector{Float64},
                                times::Vector{Float64};
                                n_permutations::Int=1000,
                                cluster_alpha::Float64=0.05,
                                test_alpha::Float64=0.05,
                                tail::Symbol=:both)

Perform a Cluster-Based Permutation Test on Time-Frequency Data.
"""
function cluster_permutation_test_tf(data1::Array{Float64, 4}, 
                                     data2::Array{Float64, 4},
                                     spatial_conn::SparseMatrixCSC{Bool},
                                     electrodes::Vector{Symbol},
                                     freqs::Vector{Float64},
                                     times::Vector{Float64};
                                     n_permutations::Int=1000,
                                     cluster_alpha::Float64=0.05,
                                     test_alpha::Float64=0.05,
                                     tail::Symbol=:both)
                                     
    n_part = size(data1, 1)
    df = n_part - 1
    
    # 1. Compute Observed T-Map
    t_obs = compute_t_map_3d(data1, data2; tail=tail)
    
    # 2. Determine Critical T-Value (Parametric) for Cluster Definition
    tdist = TDist(df)
    if tail == :both
        crit_t = quantile(tdist, 1 - cluster_alpha/2)
    else
        crit_t = quantile(tdist, 1 - cluster_alpha)
    end
    
    # 3. Find Observed Clusters
    mask_obs = abs.(t_obs) .> crit_t
    clusters_obs = find_clusters_tf(mask_obs, electrodes, freqs, times, spatial_conn)
    updated_clusters, stats_obs = compute_cluster_statistics_tf(clusters_obs, t_obs; statistic_type=:sum)
    
    # 4. Permutation Loop
    null_distribution = Float64[]
    sizehint!(null_distribution, n_permutations)
    
    # Pre-allocate buffer for permutation
    # For paired test, we flip signs of difference
    diff_data = data1 .- data2
    # Pre-compute sum of squares (constant for sign flips)
    sum_sq = dropdims(sum(diff_data.^2, dims=1), dims=1)
    
    # Progress bar
    p = Progress(n_permutations, 1, "Permuting...")
    
    # Reusable buffers for t-calculation
    sum_perm = zeros(size(t_obs))
    
    for i in 1:n_permutations
        # Random sign flip
        signs = rand([-1.0, 1.0], n_part)
        
        # Compute sum of sign-flipped differences
        fill!(sum_perm, 0.0)
        @inbounds for p_idx in 1:n_part
            s = signs[p_idx]
            for idx in CartesianIndices(sum_perm)
                sum_perm[idx] += diff_data[p_idx, idx] * s
            end
        end
        
        # Fast T-calculation: t = mean / (std / sqrt(n))
        # std = sqrt( (sum_sq - (sum^2)/n) / (n-1) )
        t_perm = similar(sum_perm)
        @inbounds for idx in CartesianIndices(sum_perm)
            s_val = sum_perm[idx]
            m_val = s_val / n_part
            # variance = (sum_sq - (sum^2)/n) / (n-1)
            var_val = (sum_sq[idx] - (s_val^2)/n_part) / (n_part - 1)
            
            if var_val > 1e-12
                t_perm[idx] = m_val / (sqrt(var_val) / sqrt(n_part))
            else
                t_perm[idx] = 0.0
            end
        end
        
        # Find Max Cluster Stat
        mask_perm = abs.(t_perm) .> crit_t
        clusters_perm = find_clusters_tf(mask_perm, electrodes, freqs, times, spatial_conn)
        
        if isempty(clusters_perm)
            push!(null_distribution, 0.0)
        else
            # Get max absolute statistic
            _, s = compute_cluster_statistics_tf(clusters_perm, t_perm; statistic_type=:sum)
            push!(null_distribution, maximum(abs.(s)))
        end
        next!(p)
    end
    
    # 5. Assign P-Values
    final_clusters = TFCluster[]
    
    for c_idx in eachindex(updated_clusters)
        c = updated_clusters[c_idx]
        obs_stat = abs(c.cluster_stat)
        
        # P-value = proportion of null >= observed
        p_val = mean(null_distribution .>= obs_stat)
        
        is_sig = p_val < test_alpha
        
        new_c = TFCluster(c.id, c.electrodes, c.freq_indices, c.time_indices,
                          c.freq_range, c.time_range, c.cluster_stat, p_val, is_sig, c.polarity, c.pixels)
        push!(final_clusters, new_c)
    end
    
    return final_clusters, null_distribution, t_obs
end

# Re-include the primitives (find_clusters_tf etc.)
# I need to duplicate them or assume this file REPLACES the previous one.
# Since I used write_to_file with Overwrite=true, I must include EVERYTHING in this file.

"""
    find_clusters_tf(...)
"""
function find_clusters_tf(mask::BitArray{3},
                          electrodes::Vector{Symbol},
                          freqs::Vector{Float64},
                          times::Vector{Float64},
                          spatial_conn::SparseMatrixCSC{Bool})
                          
    n_elec, n_freq, n_time = size(mask)
    clusters = TFCluster[]
    
    if !any(mask)
        return clusters
    end
    
    visited = falses(n_elec, n_freq, n_time)
    current_id = 0
    q = CartesianIndex{3}[]
    sizehint!(q, 1000)
    
    conn_rows = rowvals(spatial_conn)
    conn_nzrange = [nzrange(spatial_conn, col) for col in 1:n_elec]

    for e in 1:n_elec, f in 1:n_freq, t in 1:n_time
        if mask[e, f, t] && !visited[e, f, t]
            current_id += 1
            empty!(q)
            push!(q, CartesianIndex(e, f, t))
            visited[e, f, t] = true
            
            c_elec_set = Set{Symbol}()
            c_freq_set = Set{Int}()
            c_time_set = Set{Int}()
            c_pixels = CartesianIndex{3}[]
            
            push!(c_elec_set, electrodes[e])
            push!(c_freq_set, f)
            push!(c_time_set, t)
            push!(c_pixels, CartesianIndex(e, f, t))
            
            head = 1
            while head <= length(q)
                curr = q[head]
                head += 1
                ce, cf, ct = curr[1], curr[2], curr[3]
                
                # Spatial
                nzr = conn_nzrange[ce]
                for i in nzr
                     ne = conn_rows[i]
                     if mask[ne, cf, ct] && !visited[ne, cf, ct]
                         visited[ne, cf, ct] = true
                         cart = CartesianIndex(ne, cf, ct)
                         push!(q, cart)
                         push!(c_elec_set, electrodes[ne])
                         push!(c_freq_set, cf)
                         push!(c_time_set, ct)
                         push!(c_pixels, cart)
                     end
                end

                # Spectral
                for df in (-1, 1)
                    nf = cf + df
                    if nf >= 1 && nf <= n_freq
                        if mask[ce, nf, ct] && !visited[ce, nf, ct]
                            visited[ce, nf, ct] = true
                            cart = CartesianIndex(ce, nf, ct)
                            push!(q, cart)
                            push!(c_elec_set, electrodes[ce])
                            push!(c_freq_set, nf)
                            push!(c_time_set, ct)
                            push!(c_pixels, cart)
                        end
                    end
                end
                
                # Temporal
                for dt in (-1, 1)
                    nt = ct + dt
                    if nt >= 1 && nt <= n_time
                        if mask[ce, cf, nt] && !visited[ce, cf, nt]
                            visited[ce, cf, nt] = true
                            cart = CartesianIndex(ce, cf, nt)
                            push!(q, cart)
                            push!(c_elec_set, electrodes[ce])
                            push!(c_freq_set, cf)
                            push!(c_time_set, nt)
                            push!(c_pixels, cart)
                        end
                    end
                end
            end
            
            c_elec_vec = collect(c_elec_set) 
            c_freq_vec = sort(collect(c_freq_set))
            c_time_vec = sort(collect(c_time_set))
            
            f_range = (freqs[c_freq_vec[1]], freqs[c_freq_vec[end]])
            t_range = (times[c_time_vec[1]], times[c_time_vec[end]])
            
            push!(clusters, TFCluster(current_id, c_elec_vec, c_freq_vec, c_time_vec, f_range, t_range, 
                                      0.0, 1.0, false, :positive, c_pixels))
        end
    end
    
    return clusters
end

"""
    compute_cluster_statistics_tf(...)
"""
function compute_cluster_statistics_tf(clusters::Vector{TFCluster},
                                       t_matrix::Array{Float64, 3};
                                       statistic_type::Symbol=:sum)
    updated_clusters = TFCluster[]
    stats = Float64[]
    
    for c in clusters
        cluster_stat = 0.0
        if statistic_type == :sum
            for pix in c.pixels
                val = t_matrix[pix]
                if !isnan(val) && !isinf(val)
                    cluster_stat += val
                end
            end
        elseif statistic_type == :max
            cluster_stat = -Inf
            for pix in c.pixels
                val = t_matrix[pix]
                if !isnan(val) && !isinf(val)
                    if val > cluster_stat cluster_stat = val end
                end
            end
            if cluster_stat == -Inf cluster_stat = 0.0 end
        elseif statistic_type == :size
             cluster_stat = Float64(length(c.pixels))
        end
        
        updated_cluster = TFCluster(c.id, c.electrodes, c.freq_indices, c.time_indices,
                                    c.freq_range, c.time_range, cluster_stat, c.p_value, c.is_significant, c.polarity, c.pixels)
        push!(updated_clusters, updated_cluster)
        push!(stats, cluster_stat)
    end
    return updated_clusters, stats
end
