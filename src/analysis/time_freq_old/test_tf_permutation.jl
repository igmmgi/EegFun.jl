using eegfun
using SparseArrays
using Random
using Statistics

println("Testing High-Level TF Cluster Permutation Test...")

# 1. Setup Parameters
n_part = 10
n_elec = 5
n_freq = 10
n_time = 20

electrodes = [Symbol("E$i") for i in 1:n_elec]
freqs = collect(1.0:1.0:n_freq)
times = collect(range(-0.5, 0.5, length=n_time))

# 2. Generate Synthetic 4D Data [P, E, F, T]
# Condition 1: Baseline noise
data1 = randn(n_part, n_elec, n_freq, n_time)

# Condition 2: Noise + Signal in a specific cluster
data2 = randn(n_part, n_elec, n_freq, n_time)
# Add signal to data2 for participants 1:n_part, electrodes 1:2, freqs 3:5, time 10:15
# This simulates a consistent effect across participants
signal_magnitude = 2.0
data2[:, 1:2, 3:5, 10:15] .+= signal_magnitude

# 3. Spatial Connectivity (Chain)
I, J = Int[], Int[]
for i in 1:n_elec-1
    push!(I, i); push!(J, i+1)
    push!(I, i+1); push!(J, i)
end
spatial_conn = sparse(I, J, true, n_elec, n_elec)

# 4. Run Permutation Test
println("Running permutation test (100 perms)...")
# Note: Using 100 for speed in verification
clusters, null_dist, t_obs = eegfun.cluster_permutation_test_tf(
    data1, data2, spatial_conn, electrodes, freqs, times;
    n_permutations=100,
    cluster_alpha=0.05,
    test_alpha=0.05,
    tail=:both
)

println("\nTest Results:")
println("Observed Max T: $(maximum(abs.(t_obs)))")
println("Found $(length(clusters)) clusters.")

sig_clusters = filter(c -> c.is_significant, clusters)
println("Found $(length(sig_clusters)) significant clusters.")

if !isempty(sig_clusters)
    c = sig_clusters[1]
    println("Significant Cluster 1:")
    println("  P-value: $(c.p_value)")
    println("  Freq Range: $(c.freq_range)")
    println("  Time Range: $(c.time_range)")
    println("  Size (Pixels): $(length(c.pixels))")
end
