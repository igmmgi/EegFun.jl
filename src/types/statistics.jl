# Statistical test types for EEG/ERP data analysis.
#
# Defines all types used for statistical testing, including
# t-test results, cluster permutation test results, and analytic t-test results.

# === BASIC T-TEST RESULT ===

"""
    TTestResult(df::Float64, t::Float64, p::Float64)

Struct containing the results of a basic t-test.
"""
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

# === DATA PREPARATION TYPES ===

"""
    AnalysisData

Stores core analysis data for statistical tests.

# Fields
- `design::Symbol`: Design type - `:paired` or `:independent`
- `data::Vector{Array{Float64, 3}}`: Data for conditions 1 and 2 [condition1, condition2], each [participants × electrodes × time]
- `time_points::Vector{Float64}`: Time points in seconds for the analysis interval
"""
struct AnalysisData
    design::Symbol                   # :paired vs. :independent
    data::Vector{Array{Float64,3}}  # [condition1, condition2] - each [participants × electrodes × time]
    time_points::Vector{Float64}     # analysis interval
end


"""
    StatisticalData

Stores prepared data for statistical tests (both permutation and analytic tests).
Created by `prepare_stats`.

# Fields
- `data::Vector{ErpData}`: Grand average ERPs for conditions 1 and 2 (for visualization/storage)
- `analysis::AnalysisData`: Core analysis data (design, data arrays, time points)
- `se_cond1::Array{Float64, 2}`: Standard error for condition 1 over the full display interval [electrodes × time]
- `se_cond2::Array{Float64, 2}`: Standard error for condition 2 over the full display interval [electrodes × time]
- `se_diff::Array{Float64, 2}`: Standard error of the mean difference over the full display interval [electrodes × time]
"""
struct StatisticalData
    data::Vector{ErpData}
    analysis::AnalysisData
    se_cond1::Array{Float64,2}
    se_cond2::Array{Float64,2}
    se_diff::Array{Float64,2}
end

function Base.show(io::IO, data::StatisticalData)
    time_range_str(times) = isempty(times) ? "N/A" : "$(first(times)) to $(last(times)) s"

    # Get dimensions
    n_participants1 = size(data.analysis.data[1], 1)
    n_participants2 = size(data.analysis.data[2], 1)

    # Get summary info from grand averages
    n_electrodes1 = n_channels(data.data[1])
    n_electrodes2 = n_channels(data.data[2])
    time_range1_str = time_range_str(data.data[1].data[!, :time])
    time_range2_str = time_range_str(data.data[2].data[!, :time])
    analysis_time_range = time_range_str(data.analysis.time_points)

    println(io, "StatisticalData")
    println(io, "├─ Design: $(data.analysis.design)")
    println(io, "├─ Condition 1 ($(data.data[1].condition_name)): $n_participants1 participants")
    println(io, "│  └─ $(n_electrodes1) channels, $time_range1_str, $(data.data[1].sample_rate) Hz")
    println(io, "├─ Condition 2 ($(data.data[2].condition_name)): $n_participants2 participants")
    println(io, "│  └─ $(n_electrodes2) channels, $time_range2_str, $(data.data[2].sample_rate) Hz")
    println(io, "├─ Analysis time points: $analysis_time_range")
    println(io, "└─ Analysis dimensions: $(size(data.analysis.data[1])) (participants × electrodes × time)")
end

# === CLUSTER TYPES ===

"""
    Cluster

Stores information about a single cluster found in the thresholded data.

# Fields
- `id::Int`: Unique cluster identifier
- `electrodes::Vector{Symbol}`: Electrode labels in this cluster
- `time_indices::Vector{Int}`: Time point indices in this cluster
- `time_range::Tuple{Float64, Float64}`: Time range in seconds (start, end)
- `cluster_stat::Float64`: Cluster-level statistic (e.g., sum of t-values)
- `p_value::Float64`: P-value from permutation test
- `is_significant::Bool`: Whether cluster is significant (p < alpha)
- `polarity::Symbol`: `:positive` or `:negative`
- `members::Vector{Tuple{Int,Int}}`: Exact (electrode_index, time_index) pairs in this cluster
"""
struct Cluster
    id::Int
    electrodes::Vector{Symbol}
    time_indices::Vector{Int}
    time_range::Tuple{Float64,Float64}
    cluster_stat::Float64
    p_value::Float64
    is_significant::Bool
    polarity::Symbol
    members::Vector{Tuple{Int,Int}}  # (electrode_idx, time_idx) pairs
end

"""
    ClusterInfo

Stores cluster-specific parameters for cluster permutation tests.

# Fields
- `threshold_method::Symbol`: `:parametric`, `:nonparametric_individual`, or `:nonparametric_common`
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`
- `n_permutations::Int`: Number of permutations performed
"""
struct ClusterInfo
    threshold_method::Symbol
    cluster_type::Symbol
    n_permutations::Int
end

# === TEST CONFIGURATION TYPES ===

"""
    TestInfo

Stores test configuration and parameters.

# Fields
- `type::Symbol`: Design/test type - `:paired` or `:independent`
- `df::Float64`: Degrees of freedom
- `alpha::Float64`: Significance threshold
- `tail::Symbol`: Test tail - `:both`, `:left`, or `:right`
- `correction_method::Symbol`: Multiple comparison correction - `:none`, `:bonferroni`, or `:cluster_permutation`
- `cluster_info::Union{ClusterInfo, Nothing}`: Cluster-specific info (nothing for analytic tests)
"""
struct TestInfo
    type::Symbol
    df::Float64
    alpha::Float64
    tail::Symbol
    correction_method::Symbol
    cluster_info::Union{ClusterInfo,Nothing}
end

"""
    StatMatrix

Stores t-statistics and p-values matrices.

# Fields
- `t::Array{Float64, 2}`: T-statistics [electrodes × time]
- `p::Union{Array{Float64, 2}, Nothing}`: P-values [electrodes × time], or `nothing` for cluster permutation
"""
struct StatMatrix
    t::Array{Float64,2}
    p::Union{Array{Float64,2},Nothing}
end

"""
    Masks

Stores significance masks for positive and negative effects.

# Fields
- `positive::BitArray{2}`: Positive significant points [electrodes × time]
- `negative::BitArray{2}`: Negative significant points [electrodes × time]
"""
struct Masks
    positive::BitArray{2}
    negative::BitArray{2}
end

"""
    Clusters

Stores positive and negative clusters from cluster permutation tests.

# Fields
- `positive::Vector{Cluster}`: Positive clusters
- `negative::Vector{Cluster}`: Negative clusters
"""
struct Clusters
    positive::Vector{Cluster}
    negative::Vector{Cluster}
end

"""
    PermutationDistribution

Stores the null distribution of maximum cluster statistics from permutations.

# Fields
- `positive::Vector{Float64}`: Maximum positive cluster stats from each permutation
- `negative::Vector{Float64}`: Maximum negative cluster stats from each permutation
"""
struct PermutationDistribution
    positive::Vector{Float64}
    negative::Vector{Float64}
end

# === RESULT TYPES ===

"""
    StatsResult

Abstract type for statistical test results. All statistical test results share common fields:
- `test_info::TestInfo`: Test configuration and parameters
- `data::Vector{ErpData}`: Grand average ERPs for conditions 1 and 2
- `stat_matrix::StatMatrix`: T-statistics and optionally p-values
- `masks::Masks`: Significance masks for positive and negative effects
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `critical_t`: Critical t-values (type varies by test method)
- `se::Array{Float64, 2}`: Standard error of the mean difference [electrodes × time]
"""
abstract type StatsResult <: EegFunData end

"""
    PermutationResult

Stores complete results from a cluster-based permutation test.

# Fields
- `test_info::TestInfo`: Test configuration and parameters (includes ClusterInfo)
- `data::Vector{ErpData}`: Grand average ERPs for conditions 1 and 2 (for visualization)
- `stat_matrix::StatMatrix`: T-statistics matrix. Note: p is `nothing` for cluster permutation
- `masks::Masks`: Significance masks for positive and negative effects
- `clusters::Clusters`: Positive and negative clusters
- `permutation_distribution::PermutationDistribution`: Null distribution of max cluster stats
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `critical_t::Union{Array{Float64, 2}, Tuple{Float64, Float64}, Tuple{Array{Float64, 2}, Array{Float64, 2}}}`: Critical t-values used
- `se::Array{Float64, 2}`: Standard error of the mean difference [electrodes × time] (analysis interval)
- `se_cond1::Array{Float64, 2}`: Standard error for condition 1 [electrodes × time] (full display interval)
- `se_cond2::Array{Float64, 2}`: Standard error for condition 2 [electrodes × time] (full display interval)
- `se_diff::Array{Float64, 2}`: Standard error of the mean difference [electrodes × time] (full display interval)
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

function Base.show(io::IO, result::PermutationResult)
    n_electrodes = length(result.electrodes)
    n_time_points = length(result.time_points)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"

    n_pos_clusters = length(result.clusters.positive)
    n_neg_clusters = length(result.clusters.negative)
    n_sig_pos = count(c -> c.is_significant, result.clusters.positive)
    n_sig_neg = count(c -> c.is_significant, result.clusters.negative)

    # Count significant points
    n_sig_pos_points = count(result.masks.positive)
    n_sig_neg_points = count(result.masks.negative)

    # Extract nested fields for readability
    test_info = result.test_info
    cluster_info = test_info.cluster_info

    println(io, "PermutationResult")
    println(io, "├─ Design: $(test_info.type)")
    println(io, "├─ Degrees of freedom: $(round(Int, test_info.df))")
    println(io, "├─ Permutations: $(cluster_info.n_permutations)")
    println(io, "├─ Threshold: $(test_info.alpha) ($(cluster_info.threshold_method))")
    println(io, "├─ Cluster type: $(cluster_info.cluster_type)")
    println(io, "├─ Data dimensions: $n_electrodes electrodes × $n_time_points time points ($time_range)")
    println(io, "├─ Significant points: $n_sig_pos_points positive, $n_sig_neg_points negative")
    println(io, "├─ Clusters found: $n_pos_clusters positive, $n_neg_clusters negative")

    if n_sig_pos > 0 || n_sig_neg > 0
        print(io, "├─ Significant clusters: ")
        if n_sig_pos > 0
            print(io, "$n_sig_pos positive")
            if n_sig_neg > 0
                print(io, ", $n_sig_neg negative")
            end
        elseif n_sig_neg > 0
            print(io, "$n_sig_neg negative")
        end
        println(io)
    else
        println(io, "├─ Significant clusters: 0")
    end

    # Show only significant cluster details
    if n_sig_pos > 0 || n_sig_neg > 0
        println(io, "└─ Significant cluster details:")

        if n_sig_pos > 0
            println(io, "   Positive clusters:")
            sig_clusters = [c for c in result.clusters.positive if c.is_significant]

            for cluster in sig_clusters
                sig_marker = "✓"
                p_str = cluster.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", cluster.p_value)
                stat_str = Printf.@sprintf("%.2f", cluster.cluster_stat)
                time_str = "$(cluster.time_range[1])-$(cluster.time_range[2]) s"
                n_elec = length(cluster.electrodes)
                println(io, "     [$sig_marker] Cluster $(cluster.id): stat=$stat_str, p=$p_str, $n_elec electrodes, $time_str")
            end
        end

        if n_sig_neg > 0
            println(io, "   Negative clusters:")
            sig_clusters = [c for c in result.clusters.negative if c.is_significant]

            for cluster in sig_clusters
                sig_marker = "✓"
                p_str = cluster.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", cluster.p_value)
                stat_str = Printf.@sprintf("%.2f", cluster.cluster_stat)
                time_str = "$(cluster.time_range[1])-$(cluster.time_range[2]) s"
                n_elec = length(cluster.electrodes)
                println(io, "     [$sig_marker] Cluster $(cluster.id): stat=$stat_str, p=$p_str, $n_elec electrodes, $time_str")
            end
        end
    end
end

"""
    AnalyticResult

Stores results from an analytic (parametric) t-test without permutation.

# Fields
- `test_info::TestInfo`: Test configuration (type, df, alpha, tail, correction_method, cluster_info=nothing)
- `data::Vector{ErpData}`: Grand average ERPs for conditions 1 and 2 (for visualization/storage)
- `stat_matrix::StatMatrix`: Statistical results containing `t` (t-statistics) and `p` (p-values) [electrodes × time]
- `masks::Masks`: Significance masks containing `positive` and `negative` significant points [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `critical_t::Float64`: Critical t-value for significance (uniform across all points)
- `se::Array{Float64, 2}`: Standard error of the mean difference [electrodes × time] (analysis interval)
- `se_cond1::Array{Float64, 2}`: Standard error for condition 1 [electrodes × time] (full display interval)
- `se_cond2::Array{Float64, 2}`: Standard error for condition 2 [electrodes × time] (full display interval)
- `se_diff::Array{Float64, 2}`: Standard error of the mean difference [electrodes × time] (full display interval)
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

function Base.show(io::IO, result::AnalyticResult)
    n_electrodes = length(result.electrodes)
    n_time_points = length(result.time_points)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"

    n_sig_pos = count(result.masks.positive)
    n_sig_neg = count(result.masks.negative)

    println(io, "AnalyticResult")
    println(io, "├─ Test info")
    println(io, "│  ├─ Type: $(result.test_info.type)")
    println(io, "│  ├─ DF: $(result.test_info.df)")
    println(io, "│  ├─ Alpha: $(result.test_info.alpha)")
    println(io, "│  ├─ Tail: $(result.test_info.tail)")
    println(io, "│  └─ Correction method: $(result.test_info.correction_method)")
    println(io, "├─ Data dimensions: $n_electrodes electrodes × $n_time_points time points ($time_range)")
    println(io, "└─ Significant points: $n_sig_pos positive, $n_sig_neg negative")
end

"""
    LmmStatsResult

Result of a mass-univariate linear mixed model fit across a channels × time grid
(see [`fit_mass_lmm`](@ref)).

# Fields
- `coefficients::Vector{String}`: Fixed-effect names (e.g. "(Intercept)", "condition: B")
- `channels::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points
- `beta::Array{Float64,3}`: Fixed-effect estimates [channels × time × coefficients]
- `se::Array{Float64,3}`: Standard errors [channels × time × coefficients]
- `z_values::Array{Float64,3}`: Wald z-statistics (β̂ / SE) [channels × time × coefficients]
- `p_parametric::Array{Float64,3}`: Uncorrected two-sided Wald p-values (Normal reference)
- `p_corrected::Union{Array{Float64,3},Nothing}`: Permutation FWER-corrected p-values (max-statistic,
  cluster-mass or TFCE, see `inference.correction`); `NaN` for untested coefficients; `nothing` without permutations
- `null_max::Matrix{Float64}`: Null distribution of the maximum statistic [permutations × coefficients]
  (`NaN` columns for untested coefficients)
- `tested::Vector{Int}`: Indices of the permutation-tested coefficients
- `singular_fits::Matrix{Bool}`: Whether the fit at each [channel, time] point is singular
- `inference::NamedTuple`: Settings used for inference (`correction`, `n_perms`, `exhaustive`, `flip_by`,
  `n_blocks`, `cluster_threshold`, `tfce`, `spatial_connectivity`, `REML`, `n_obs`, `formula`)
- `epochs::Union{EpochData,Nothing}`: The input epochs (EpochData method only; used for plotting)
"""
struct LmmStatsResult <: EegFunData
    coefficients::Vector{String}
    channels::Vector{Symbol}
    time_points::Vector{Float64}
    beta::Array{Float64, 3}
    se::Array{Float64, 3}
    z_values::Array{Float64, 3}
    p_parametric::Array{Float64, 3}
    p_corrected::Union{Array{Float64, 3}, Nothing}
    null_max::Matrix{Float64}
    tested::Vector{Int}
    singular_fits::Matrix{Bool}
    inference::NamedTuple
    epochs::Union{EpochData, Nothing}
end

function Base.show(io::IO, result::LmmStatsResult)
    n_chans = length(result.channels)
    n_times = length(result.time_points)
    n_coefs = length(result.coefficients)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"
    inf = result.inference
    n_singular = count(result.singular_fits)
    total_fits = n_chans * n_times

    println(io, "LmmStatsResult")
    println(io, "├─ Dimensions: $n_chans channels × $n_times time points ($time_range) × $n_coefs coefficients")
    println(io, "├─ Predictors: ", join(result.coefficients, ", "))
    sing_pct = total_fits > 0 ? round(100 * n_singular / total_fits, digits = 1) : 0.0
    println(io, "├─ Singular fits: $n_singular / $total_fits ($sing_pct%)")
    if isnothing(result.p_corrected)
        println(io, "└─ Inference: Wald z only (no permutation correction; p_parametric is uncorrected)")
        return
    end
    method = inf.correction === :maxt ? "max-|z|" : inf.correction === :cluster ? "cluster mass (threshold $(inf.cluster_threshold))" : "TFCE"
    perms = inf.exhaustive ? "all $(inf.n_perms) sign-flip patterns" : "$(inf.n_perms) random sign flips"
    println(io, "├─ Inference: FWER $method, whole-block sign flipping by :$(inf.flip_by) ($(inf.n_blocks) blocks, $perms)")
    sig_strs = ["$(result.coefficients[i]): $(count(<(0.05), @view(result.p_corrected[:, :, i])))" for i in result.tested]
    print(io, "└─ Significant points (corrected p < 0.05): ", join(sig_strs, ", "))
end

Base.show(io::IO, ::MIME"text/plain", result::LmmStatsResult) = show(io, result)

