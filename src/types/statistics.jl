struct TTestResult
    df::Float64
    t::Float64
    p::Float64
end

function Base.show(io::IO, result::TTestResult)
    t_str = isnan(result.t) ? "NaN" : (isinf(result.t) ? (result.t > 0 ? "Inf" : "-Inf") : Printf.@sprintf("%.4f", result.t))
    p_str = isnan(result.p) ? "NaN" : Printf.@sprintf("%.4f", result.p)
    df_str = isnan(result.df) ? "NaN" : string(round(Int, result.df))
    print(io, "t($df_str) = $t_str, p = $p_str")
end

struct AnalysisData
    design::Symbol
    data::Vector{Array{Float64,3}}  # [condition1, condition2], each [participants × electrodes × time]
    time_points::Vector{Float64}
end

"""
    StatisticalData

Prepared ERP data for permutation or analytic statistical testing (`prepare_stats`).

# Fields
- `data::Vector{ErpData}`: Grand-average ERPs for conditions 1 and 2
- `analysis::AnalysisData`: Design, 3D participant arrays `[participants × electrodes × time]`, and analysis time points
- `se_cond1::Array{Float64,2}`: Condition 1 standard error `[electrodes × time]`
- `se_cond2::Array{Float64,2}`: Condition 2 standard error `[electrodes × time]`
- `se_diff::Array{Float64,2}`: Mean-difference standard error `[electrodes × time]`
"""
struct StatisticalData
    data::Vector{ErpData}
    analysis::AnalysisData
    se_cond1::Array{Float64,2}
    se_cond2::Array{Float64,2}
    se_diff::Array{Float64,2}
end

_range_str(v::AbstractVector, unit::AbstractString) = isempty(v) ? "N/A" : "$(first(v)) to $(last(v)) $unit"

function Base.show(io::IO, data::StatisticalData)
    n_p1 = size(data.analysis.data[1], 1)
    n_p2 = size(data.analysis.data[2], 1)
    d1, d2 = data.data[1], data.data[2]

    println(io, "StatisticalData")
    println(io, "├─ Design: $(data.analysis.design)")
    println(io, "├─ Condition 1 ($(d1.condition_name)): $n_p1 participants")
    println(io, "│  └─ $(n_channels(d1)) channels, $(_range_str(d1.data[!, :time], "s")), $(d1.sample_rate) Hz")
    println(io, "├─ Condition 2 ($(d2.condition_name)): $n_p2 participants")
    println(io, "│  └─ $(n_channels(d2)) channels, $(_range_str(d2.data[!, :time], "s")), $(d2.sample_rate) Hz")
    println(io, "├─ Analysis time points: $(_range_str(data.analysis.time_points, "s"))")
    println(io, "└─ Analysis dimensions: $(size(data.analysis.data[1])) (participants × electrodes × time)")
end

struct Cluster
    id::Int
    electrodes::Vector{Symbol}
    time_indices::Vector{Int}
    time_range::Tuple{Float64,Float64}
    cluster_stat::Float64
    p_value::Float64
    is_significant::Bool
    polarity::Symbol
    members::Vector{Tuple{Int,Int}}
end

struct ClusterInfo
    threshold_method::Symbol
    cluster_type::Symbol
    n_permutations::Int
    cluster_alpha::Float64
end

ClusterInfo(threshold_method::Symbol, cluster_type::Symbol, n_permutations::Int) =
    ClusterInfo(threshold_method, cluster_type, n_permutations, 0.05)

struct TestInfo
    type::Symbol
    df::Float64
    alpha::Float64
    tail::Symbol
    correction_method::Symbol
    cluster_info::Union{ClusterInfo,Nothing}
end

struct StatMatrix
    t::Array{Float64,2}
    p::Union{Array{Float64,2},Nothing}
end

struct Masks
    positive::BitArray{2}
    negative::BitArray{2}
end

struct Clusters
    positive::Vector{Cluster}
    negative::Vector{Cluster}
end

struct PermutationDistribution
    positive::Vector{Float64}
    negative::Vector{Float64}
end

"""
    StatsResult <: EegFunData

Abstract supertype for ERP and time-frequency statistical test results.
"""
abstract type StatsResult <: EegFunData end

"""
    PermutationResult <: StatsResult

Cluster-based permutation test results on ERP data (`permutation_test`).

# Fields
- `test_info::TestInfo`: Test design, degrees of freedom, alpha, tail, and `ClusterInfo`
- `data::Vector{ErpData}`: Grand-average ERPs for conditions 1 and 2
- `stat_matrix::StatMatrix`: T-statistic matrix `[electrodes × time]`
- `masks::Masks`: Significant positive and negative masks `[electrodes × time]`
- `clusters::Clusters`: Detected positive and negative clusters
- `permutation_distribution::PermutationDistribution`: Null distribution of maximum cluster statistics
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Analysis time points in seconds
- `critical_t::Union{Array{Float64,2},Tuple{Float64,Float64},Tuple{Array{Float64,2},Array{Float64,2}}}`: Cluster-forming critical t-values
- `se::Array{Float64,2}`: Mean-difference standard error over the analysis interval
- `se_cond1::Array{Float64,2}`: Condition 1 standard error over the display interval
- `se_cond2::Array{Float64,2}`: Condition 2 standard error over the display interval
- `se_diff::Array{Float64,2}`: Mean-difference standard error over the display interval
"""
struct PermutationResult <: StatsResult
    test_info::TestInfo
    data::Vector{ErpData}
    stat_matrix::StatMatrix
    masks::Masks
    clusters::Clusters
    permutation_distribution::PermutationDistribution
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
    critical_t::Union{Array{Float64,2},Tuple{Float64,Float64},Tuple{Array{Float64,2},Array{Float64,2}}}
    se::Array{Float64,2}
    se_cond1::Array{Float64,2}
    se_cond2::Array{Float64,2}
    se_diff::Array{Float64,2}
end

function _print_cluster_line(io::IO, c::Cluster)
    p_str = c.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", c.p_value)
    stat_str = Printf.@sprintf("%.2f", c.cluster_stat)
    println(
        io,
        "     [✓] Cluster $(c.id): stat=$stat_str, p=$p_str, $(length(c.electrodes)) electrodes, $(c.time_range[1])-$(c.time_range[2]) s",
    )
end

function _print_sig_cluster_counts(io::IO, n_sig_pos::Int, n_sig_neg::Int, zero_prefix::AbstractString)
    if n_sig_pos > 0 || n_sig_neg > 0
        print(io, "├─ Significant clusters: ")
        if n_sig_pos > 0
            print(io, "$n_sig_pos positive")
            n_sig_neg > 0 && print(io, ", $n_sig_neg negative")
        else
            print(io, "$n_sig_neg negative")
        end
        println(io)
    else
        println(io, "$zero_prefix Significant clusters: 0")
    end
end

function Base.show(io::IO, result::PermutationResult)
    n_sig_pos = count(c -> c.is_significant, result.clusters.positive)
    n_sig_neg = count(c -> c.is_significant, result.clusters.negative)
    ti = result.test_info
    ci = ti.cluster_info

    println(io, "PermutationResult")
    println(io, "├─ Design: $(ti.type)")
    println(io, "├─ Degrees of freedom: $(round(Int, ti.df))")
    println(io, "├─ Permutations: $(ci.n_permutations)")
    println(io, "├─ Alpha: $(ti.alpha), Cluster alpha: $(ci.cluster_alpha) ($(ci.threshold_method))")
    println(io, "├─ Cluster type: $(ci.cluster_type)")
    println(
        io,
        "├─ Data dimensions: $(length(result.electrodes)) electrodes × $(length(result.time_points)) time points ($(_range_str(result.time_points, "s")))",
    )
    println(io, "├─ Significant points: $(count(result.masks.positive)) positive, $(count(result.masks.negative)) negative")
    println(io, "├─ Clusters found: $(length(result.clusters.positive)) positive, $(length(result.clusters.negative)) negative")
    _print_sig_cluster_counts(io, n_sig_pos, n_sig_neg, "├─")

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
    AnalyticResult <: StatsResult

Parametric point-wise t-test results on ERP data (`analytic_test`).

# Fields
- `test_info::TestInfo`: Test design, degrees of freedom, alpha, tail, and multiple-comparison correction
- `data::Vector{ErpData}`: Grand-average ERPs for conditions 1 and 2
- `stat_matrix::StatMatrix`: T-statistic and p-value matrices `[electrodes × time]`
- `masks::Masks`: Significant positive and negative masks `[electrodes × time]`
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Analysis time points in seconds
- `critical_t::Float64`: Critical t-threshold
- `se::Array{Float64,2}`: Mean-difference standard error over the analysis interval
- `se_cond1::Array{Float64,2}`: Condition 1 standard error over the display interval
- `se_cond2::Array{Float64,2}`: Condition 2 standard error over the display interval
- `se_diff::Array{Float64,2}`: Mean-difference standard error over the display interval
"""
struct AnalyticResult <: StatsResult
    test_info::TestInfo
    data::Vector{ErpData}
    stat_matrix::StatMatrix
    masks::Masks
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
    critical_t::Float64
    se::Array{Float64,2}
    se_cond1::Array{Float64,2}
    se_cond2::Array{Float64,2}
    se_diff::Array{Float64,2}
end

function _print_test_info_tree(io::IO, ti::TestInfo)
    println(io, "├─ Test info")
    println(io, "│  ├─ Type: $(ti.type)")
    println(io, "│  ├─ DF: $(ti.df)")
    println(io, "│  ├─ Alpha: $(ti.alpha)")
    println(io, "│  ├─ Tail: $(ti.tail)")
    println(io, "│  └─ Correction method: $(ti.correction_method)")
end

function Base.show(io::IO, result::AnalyticResult)
    println(io, "AnalyticResult")
    _print_test_info_tree(io, result.test_info)
    println(
        io,
        "├─ Data dimensions: $(length(result.electrodes)) electrodes × $(length(result.time_points)) time points ($(_range_str(result.time_points, "s")))",
    )
    println(io, "└─ Significant points: $(count(result.masks.positive)) positive, $(count(result.masks.negative)) negative")
end
