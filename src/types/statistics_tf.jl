"""
Time-Frequency Statistical Types for EEG/ERP data analysis.

Separated from main statistics types to handle 3D data (Electrodes × Frequencies × Time).
"""

"""
    TFCluster

Stores information about a single cluster found in Time-Frequency data.

# Fields
- `id::Int`: Unique cluster identifier
- `electrodes::Vector{Symbol}`: Electrode labels
- `freq_indices::Vector{Int}`: Frequency indices
- `time_indices::Vector{Int}`: Time indices
- `freq_range::Tuple{Float64, Float64}`: Frequency range in Hz
- `time_range::Tuple{Float64, Float64}`: Time range in seconds
- `cluster_stat::Float64`: Cluster-level statistic
- `p_value::Float64`: P-value
- `is_significant::Bool`: Whether cluster is significant
- `polarity::Symbol`: `:positive` or `:negative`
"""
struct TFCluster
    id::Int
    electrodes::Vector{Symbol}
    freq_indices::Vector{Int}
    time_indices::Vector{Int}
    freq_range::Tuple{Float64,Float64}
    time_range::Tuple{Float64,Float64}
    cluster_stat::Float64
    p_value::Float64
    is_significant::Bool
    polarity::Symbol
    pixels::Vector{CartesianIndex{3}} # Exact points (e, f, t)
end

"""
    TFClusters

Stores positive and negative TF clusters.
"""
struct TFClusters
    positive::Vector{TFCluster}
    negative::Vector{TFCluster}
end

"""
    TFClusterPermutationResult

Stores complete results from a TF cluster-based permutation test.

# Fields
- `test_info::TestInfo`: Test configuration
- `stat_matrix::Array{Float64, 3}`: T-statistics [electrodes × freqs × time]
- `masks_positive::BitArray{3}`: Significant mask (positive)
- `masks_negative::BitArray{3}`: Significant mask (negative)
- `clusters::TFClusters`: Found clusters
- `permutation_distribution::PermutationDistribution`: Null distribution
- `electrodes::Vector{Symbol}`: Channels
- `frequencies::Vector{Float64}`: Frequencies
- `time_points::Vector{Float64}`: Time points
"""
struct TFClusterPermutationResult
    test_info::TestInfo
    stat_matrix::Array{Float64,3}
    masks_positive::BitArray{3}
    masks_negative::BitArray{3}
    clusters::TFClusters
    permutation_distribution::PermutationDistribution
    electrodes::Vector{Symbol}
    frequencies::Vector{Float64}
    time_points::Vector{Float64}
end
