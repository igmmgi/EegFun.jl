"""
Time-Frequency Statistical Types for EEG/ERP data analysis.

Separated from main statistics types to handle 3D data (Electrodes × Frequencies × Time).
"""

"""
    TFStatsResult

Abstract supertype for TF statistical results, used for dispatch in plotting functions.
"""
abstract type TFStatsResult <: StatsResult end

# ==============
# DATA PREPARATION TYPES
# ==============

"""
    TFAnalysisData

Stores core analysis data for TF statistical tests.

# Fields
- `design::Symbol`: Design type - `:paired` or `:independent`
- `data::Vector{Array{Float64, 4}}`: Data for conditions 1 and 2, each [participants × electrodes × frequencies × time]
- `frequencies::Vector{Float64}`: Frequency values in Hz
- `time_points::Vector{Float64}`: Time points in seconds for the analysis interval
"""
struct TFAnalysisData
    design::Symbol
    data::Vector{Array{Float64,4}}  # [condition1, condition2] - each [participants × electrodes × freqs × time]
    frequencies::Vector{Float64}
    time_points::Vector{Float64}
end

"""
    TFStatisticalData

Stores prepared TF data for statistical tests (both permutation and analytic tests).

# Fields
- `data::Vector{TimeFreqData}`: Grand average TF data for conditions 1 and 2 (for visualization/storage)
- `analysis::TFAnalysisData`: Core analysis data (design, 4D data arrays, frequencies, time points)
"""
struct TFStatisticalData
    data::Vector{TimeFreqData}
    analysis::TFAnalysisData
end

function Base.show(io::IO, data::TFStatisticalData)
    time_range_str(times) = isempty(times) ? "N/A" : "$(first(times)) to $(last(times)) s"
    freq_range_str(freqs) = isempty(freqs) ? "N/A" : "$(first(freqs)) to $(last(freqs)) Hz"

    n_participants1 = size(data.analysis.data[1], 1)
    n_participants2 = size(data.analysis.data[2], 1)
    n_electrodes = size(data.analysis.data[1], 2)
    n_freqs = size(data.analysis.data[1], 3)
    n_time = size(data.analysis.data[1], 4)

    println(io, "TFStatisticalData")
    println(io, "├─ Design: $(data.analysis.design)")
    println(io, "├─ Condition 1 ($(data.data[1].condition_name)): $n_participants1 participants")
    println(io, "├─ Condition 2 ($(data.data[2].condition_name)): $n_participants2 participants")
    println(io, "├─ Frequencies: $(freq_range_str(data.analysis.frequencies)) ($(n_freqs) bins)")
    println(io, "├─ Time points: $(time_range_str(data.analysis.time_points))")
    println(
        io,
        "└─ Analysis dimensions: $(n_participants1) × $(n_electrodes) × $(n_freqs) × $(n_time) (participants × electrodes × freqs × time)",
    )
end

# ==============
# STAT MATRIX AND MASK TYPES
# ==============

"""
    TFStatMatrix

Stores t-statistics and p-values matrices for TF data.

# Fields
- `t::Array{Float64, 3}`: T-statistics [electrodes × frequencies × time]
- `p::Union{Array{Float64, 3}, Nothing}`: P-values [electrodes × frequencies × time], or `nothing` for cluster permutation
"""
struct TFStatMatrix
    t::Array{Float64,3}
    p::Union{Array{Float64,3},Nothing}
end

"""
    TFMasks

Stores significance masks for positive and negative effects in TF data.

# Fields
- `positive::BitArray{3}`: Positive significant points [electrodes × frequencies × time]
- `negative::BitArray{3}`: Negative significant points [electrodes × frequencies × time]
"""
struct TFMasks
    positive::BitArray{3}
    negative::BitArray{3}
end

# ==============
# CLUSTER TYPES
# ==============

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
- `pixels::Vector{CartesianIndex{3}}`: Exact points (electrode, frequency, time)
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

# ==============
# RESULT TYPES
# ==============

"""
    TFClusterPermutationResult

Stores complete results from a TF cluster-based permutation test.

# Fields
- `test_info::TestInfo`: Test configuration
- `data::Vector{TimeFreqData}`: Grand average TF data for conditions 1 and 2
- `stat_matrix::TFStatMatrix`: T-statistics and optionally p-values [electrodes × freqs × time]
- `masks::TFMasks`: Significance masks for positive and negative effects
- `clusters::TFClusters`: Found clusters
- `permutation_distribution::PermutationDistribution`: Null distribution
- `electrodes::Vector{Symbol}`: Channels
- `frequencies::Vector{Float64}`: Frequencies
- `time_points::Vector{Float64}`: Time points
"""
struct TFClusterPermutationResult <: TFStatsResult
    test_info::TestInfo
    data::Vector{TimeFreqData}
    stat_matrix::TFStatMatrix
    masks::TFMasks
    clusters::TFClusters
    permutation_distribution::PermutationDistribution
    electrodes::Vector{Symbol}
    frequencies::Vector{Float64}
    time_points::Vector{Float64}
end

function Base.show(io::IO, result::TFClusterPermutationResult)
    n_electrodes = length(result.electrodes)
    n_freqs = length(result.frequencies)
    n_time = length(result.time_points)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"
    freq_range = isempty(result.frequencies) ? "N/A" : "$(first(result.frequencies)) to $(last(result.frequencies)) Hz"

    n_pos_clusters = length(result.clusters.positive)
    n_neg_clusters = length(result.clusters.negative)
    n_sig_pos = count(c -> c.is_significant, result.clusters.positive)
    n_sig_neg = count(c -> c.is_significant, result.clusters.negative)

    n_sig_pos_points = count(result.masks.positive)
    n_sig_neg_points = count(result.masks.negative)

    test_info = result.test_info
    cluster_info = test_info.cluster_info

    println(io, "TFClusterPermutationResult")
    println(io, "├─ Design: $(test_info.type)")
    println(io, "├─ Degrees of freedom: $(round(Int, test_info.df))")
    println(io, "├─ Permutations: $(cluster_info.n_permutations)")
    println(io, "├─ Threshold: $(test_info.alpha) ($(cluster_info.threshold_method))")
    println(io, "├─ Cluster type: $(cluster_info.cluster_type)")
    println(io, "├─ Data dimensions: $n_electrodes electrodes × $n_freqs freqs × $n_time time points")
    println(io, "│  └─ $freq_range, $time_range")
    println(io, "├─ Significant points: $n_sig_pos_points positive, $n_sig_neg_points negative")
    println(io, "├─ Clusters found: $n_pos_clusters positive, $n_neg_clusters negative")

    if n_sig_pos > 0 || n_sig_neg > 0
        print(io, "├─ Significant clusters: ")
        if n_sig_pos > 0
            print(io, "$n_sig_pos positive")
            n_sig_neg > 0 && print(io, ", $n_sig_neg negative")
        elseif n_sig_neg > 0
            print(io, "$n_sig_neg negative")
        end
        println(io)

        println(io, "└─ Significant cluster details:")
        if n_sig_pos > 0
            println(io, "   Positive clusters:")
            for cluster in filter(c -> c.is_significant, result.clusters.positive)
                p_str = cluster.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", cluster.p_value)
                stat_str = Printf.@sprintf("%.2f", cluster.cluster_stat)
                time_str = "$(cluster.time_range[1])-$(cluster.time_range[2]) s"
                freq_str = "$(cluster.freq_range[1])-$(cluster.freq_range[2]) Hz"
                n_elec = length(cluster.electrodes)
                println(io, "     [✓] Cluster $(cluster.id): stat=$stat_str, p=$p_str, $n_elec electrodes, $freq_str, $time_str")
            end
        end
        if n_sig_neg > 0
            println(io, "   Negative clusters:")
            for cluster in filter(c -> c.is_significant, result.clusters.negative)
                p_str = cluster.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", cluster.p_value)
                stat_str = Printf.@sprintf("%.2f", cluster.cluster_stat)
                time_str = "$(cluster.time_range[1])-$(cluster.time_range[2]) s"
                freq_str = "$(cluster.freq_range[1])-$(cluster.freq_range[2]) Hz"
                n_elec = length(cluster.electrodes)
                println(io, "     [✓] Cluster $(cluster.id): stat=$stat_str, p=$p_str, $n_elec electrodes, $freq_str, $time_str")
            end
        end
    else
        println(io, "└─ Significant clusters: 0")
    end
end

"""
    TFAnalyticResult

Stores results from an analytic (parametric) t-test on TF data.

# Fields
- `test_info::TestInfo`: Test configuration
- `data::Vector{TimeFreqData}`: Grand average TF data for conditions 1 and 2
- `stat_matrix::TFStatMatrix`: T-statistics and p-values [electrodes × freqs × time]
- `masks::TFMasks`: Significance masks for positive and negative effects
- `electrodes::Vector{Symbol}`: Electrode labels
- `frequencies::Vector{Float64}`: Frequency values
- `time_points::Vector{Float64}`: Time points in seconds
- `critical_t::Float64`: Critical t-value for significance
"""
struct TFAnalyticResult <: TFStatsResult
    test_info::TestInfo
    data::Vector{TimeFreqData}
    stat_matrix::TFStatMatrix
    masks::TFMasks
    electrodes::Vector{Symbol}
    frequencies::Vector{Float64}
    time_points::Vector{Float64}
    critical_t::Float64
end

function Base.show(io::IO, result::TFAnalyticResult)
    n_electrodes = length(result.electrodes)
    n_freqs = length(result.frequencies)
    n_time = length(result.time_points)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"
    freq_range = isempty(result.frequencies) ? "N/A" : "$(first(result.frequencies)) to $(last(result.frequencies)) Hz"

    n_sig_pos = count(result.masks.positive)
    n_sig_neg = count(result.masks.negative)

    println(io, "TFAnalyticResult")
    println(io, "├─ Test info")
    println(io, "│  ├─ Type: $(result.test_info.type)")
    println(io, "│  ├─ DF: $(result.test_info.df)")
    println(io, "│  ├─ Alpha: $(result.test_info.alpha)")
    println(io, "│  ├─ Tail: $(result.test_info.tail)")
    println(io, "│  └─ Correction method: $(result.test_info.correction_method)")
    println(io, "├─ Data dimensions: $n_electrodes electrodes × $n_freqs freqs × $n_time time points")
    println(io, "│  └─ $freq_range, $time_range")
    println(io, "└─ Significant points: $n_sig_pos positive, $n_sig_neg negative")
end
