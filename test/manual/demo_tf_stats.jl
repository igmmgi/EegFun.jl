using eegfun
using SparseArrays
using Random

println("Demonstrating 3D Time-Frequency Statistics...")

# 1. Setup Mock Data
n_elec = 32
n_freq = 40
n_time = 50
electrodes = [Symbol("E$i") for i in 1:n_elec]
freqs = collect(1.0:1.0:40.0)
times = collect(range(-0.5, 1.5, length=n_time))

# 2. Create Spatial Connectivity (Nearest neighbor chain for demo)
# In reality, this comes from Layout
I = Int[]
J = Int[]
for i in 1:n_elec-1
    push!(I, i); push!(J, i+1)
    push!(I, i+1); push!(J, i)
end
spatial_conn = sparse(I, J, true, n_elec, n_elec)

# 3. Create Synthetic Significant Mask with a known cluster
mask = falses(n_elec, n_freq, n_time)
t_matrix = zeros(n_elec, n_freq, n_time)

# Make a "blob" significant
# Electrodes 1-3, Freqs 10-15Hz, Time 0.0-0.5s
f_idxs = 10:15
t_idxs = 25:35 # approx 0.0 to 0.5s
e_idxs = 1:3

for e in e_idxs, f in f_idxs, t in t_idxs
    mask[e, f, t] = true
    t_matrix[e, f, t] = 2.5 # t-value
end

println("Created significant cluster with $(sum(mask)) pixels.")

# 4. Find Clusters
println("Finding clusters...")
clusters = eegfun.find_clusters_tf(mask, electrodes, freqs, times, spatial_conn)

println("Found $(length(clusters)) clusters.")
if !isempty(clusters)
    c = clusters[1]
    println("Cluster 1:")
    println("  Electrodes: $(c.electrodes)")
    println("  Freq Range: $(c.freq_range) Hz")
    println("  Time Range: $(c.time_range) s")
    println("  Pixels: $(length(c.pixels))")
end

# 5. Compute Stats
println("Computing stats...")
updated_clusters, stats = eegfun.compute_cluster_statistics_tf(clusters, t_matrix; statistic_type=:sum)

if !isempty(updated_clusters)
    println("Cluster 1 Stat: $(updated_clusters[1].cluster_stat) (Expected: $(sum(mask)*2.5))")
end
