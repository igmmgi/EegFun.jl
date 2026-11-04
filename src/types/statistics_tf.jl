"""
    TFStatsResult <: StatsResult

Abstract supertype for 3D time-frequency statistical test results (`[electrodes × frequencies × time]`).
"""
abstract type TFStatsResult <: StatsResult end

struct TFAnalysisData
    design::Symbol
    data::Vector{Array{Float64,4}}  # [condition1, condition2], each [participants × electrodes × freqs × time]
    frequencies::Vector{Float64}
    time_points::Vector{Float64}
end

"""
    TFStatisticalData

Prepared time-frequency data for permutation or analytic statistical testing (`prepare_stats`).

# Fields
- `data::Vector{TimeFreqData}`: Grand-average TF data for conditions 1 and 2
- `analysis::TFAnalysisData`: Design, 4D participant arrays `[participants × electrodes × frequencies × time]`, frequencies, and time points
"""
struct TFStatisticalData
    data::Vector{TimeFreqData}
    analysis::TFAnalysisData
end

function Base.show(io::IO, data::TFStatisticalData)
    n_p1, n_elec, n_freqs, n_time = size(data.analysis.data[1])
    n_p2 = size(data.analysis.data[2], 1)

    println(io, "TFStatisticalData")
    println(io, "├─ Design: $(data.analysis.design)")
    println(io, "├─ Condition 1 ($(data.data[1].condition_name)): $n_p1 participants")
    println(io, "├─ Condition 2 ($(data.data[2].condition_name)): $n_p2 participants")
    println(io, "├─ Frequencies: $(_range_str(data.analysis.frequencies, "Hz")) ($n_freqs bins)")
    println(io, "├─ Time points: $(_range_str(data.analysis.time_points, "s"))")
    println(io, "└─ Analysis dimensions: $n_p1 × $n_elec × $n_freqs × $n_time (participants × electrodes × freqs × time)")
end

struct TFStatMatrix
    t::Array{Float64,3}
    p::Union{Array{Float64,3},Nothing}
end

struct TFMasks
    positive::BitArray{3}
    negative::BitArray{3}
end

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
    pixels::Vector{CartesianIndex{3}}
end

struct TFClusters
    positive::Vector{TFCluster}
    negative::Vector{TFCluster}
end

"""
    TFClusterPermutationResult <: TFStatsResult

Cluster-based permutation test results on 3D time-frequency data (`permutation_test`).

# Fields
- `test_info::TestInfo`: Test design, degrees of freedom, alpha, tail, and `ClusterInfo`
- `data::Vector{TimeFreqData}`: Grand-average TF data for conditions 1 and 2
- `stat_matrix::TFStatMatrix`: T-statistic matrix `[electrodes × frequencies × time]`
- `masks::TFMasks`: Significant positive and negative masks `[electrodes × frequencies × time]`
- `clusters::TFClusters`: Detected positive and negative 3D TF clusters
- `permutation_distribution::PermutationDistribution`: Null distribution of maximum cluster statistics
- `electrodes::Vector{Symbol}`: Electrode labels
- `frequencies::Vector{Float64}`: Frequency bins in Hz
- `time_points::Vector{Float64}`: Time points in seconds
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

function _print_cluster_line(io::IO, c::TFCluster)
    p_str = c.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", c.p_value)
    stat_str = Printf.@sprintf("%.2f", c.cluster_stat)
    println(
        io,
        "     [✓] Cluster $(c.id): stat=$stat_str, p=$p_str, $(length(c.electrodes)) electrodes, $(c.freq_range[1])-$(c.freq_range[2]) Hz, $(c.time_range[1])-$(c.time_range[2]) s",
    )
end

function Base.show(io::IO, result::TFClusterPermutationResult)
    n_sig_pos = count(c -> c.is_significant, result.clusters.positive)
    n_sig_neg = count(c -> c.is_significant, result.clusters.negative)
    ti = result.test_info
    ci = ti.cluster_info

    println(io, "TFClusterPermutationResult")
    println(io, "├─ Design: $(ti.type)")
    println(io, "├─ Degrees of freedom: $(round(Int, ti.df))")
    println(io, "├─ Permutations: $(ci.n_permutations)")
    println(io, "├─ Alpha: $(ti.alpha), Cluster alpha: $(ci.cluster_alpha) ($(ci.threshold_method))")
    println(io, "├─ Cluster type: $(ci.cluster_type)")
    println(
        io,
        "├─ Data dimensions: $(length(result.electrodes)) electrodes × $(length(result.frequencies)) freqs × $(length(result.time_points)) time points",
    )
    println(io, "│  └─ $(_range_str(result.frequencies, "Hz")), $(_range_str(result.time_points, "s"))")
    println(io, "├─ Significant points: $(count(result.masks.positive)) positive, $(count(result.masks.negative)) negative")
    println(io, "├─ Clusters found: $(length(result.clusters.positive)) positive, $(length(result.clusters.negative)) negative")
    _print_sig_cluster_counts(io, n_sig_pos, n_sig_neg, "└─")

    if n_sig_pos > 0 || n_sig_neg > 0
        println(io, "└─ Significant cluster details:")
        if n_sig_pos > 0
            println(io, "   Positive clusters:")
            for c in result.clusters.positive
                c.is_significant && _print_cluster_line(io, c)
            end
        end
        if n_sig_neg > 0
            println(io, "   Negative clusters:")
            for c in result.clusters.negative
                c.is_significant && _print_cluster_line(io, c)
            end
        end
    end
end

"""
    TFAnalyticResult <: TFStatsResult

Parametric point-wise t-test results on 3D time-frequency data (`analytic_test`).

# Fields
- `test_info::TestInfo`: Test design, degrees of freedom, alpha, tail, and multiple-comparison correction
- `data::Vector{TimeFreqData}`: Grand-average TF data for conditions 1 and 2
- `stat_matrix::TFStatMatrix`: T-statistic and p-value matrices `[electrodes × frequencies × time]`
- `masks::TFMasks`: Significant positive and negative masks `[electrodes × frequencies × time]`
- `electrodes::Vector{Symbol}`: Electrode labels
- `frequencies::Vector{Float64}`: Frequency bins in Hz
- `time_points::Vector{Float64}`: Time points in seconds
- `critical_t::Float64`: Critical t-threshold
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
    println(io, "TFAnalyticResult")
    _print_test_info_tree(io, result.test_info)
    println(
        io,
        "├─ Data dimensions: $(length(result.electrodes)) electrodes × $(length(result.frequencies)) freqs × $(length(result.time_points)) time points",
    )
    println(io, "│  └─ $(_range_str(result.frequencies, "Hz")), $(_range_str(result.time_points, "s"))")
    println(io, "└─ Significant points: $(count(result.masks.positive)) positive, $(count(result.masks.negative)) negative")
end
