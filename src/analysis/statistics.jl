"""
Statistical test functions for EEG/ERP data analysis.

This module provides t-test functions for use in cluster-based permutation tests
and other statistical analyses. Functions are provided to compute t-values and/or
p-values for paired and independent sample t-tests.
"""

# this seems a bit lighter than a full struct
const TTestResult = NamedTuple{(:df, :t, :p), Tuple{Float64, Float64, Float64}}

# Format t-value and p-value for display
function Base.show(io::IO, result::TTestResult)
    t_str = isnan(result.t) ? "NaN" : (isinf(result.t) ? (result.t > 0 ? "Inf" : "-Inf") : Printf.@sprintf("%.4f", result.t))
    p_str = isnan(result.p) ? "NaN" : Printf.@sprintf("%.4f", result.p)
    df_str = isnan(result.df) ? "NaN" : string(round(Int, result.df))
    print(io, "t($df_str) = $t_str, p = $p_str")
end


# ==============
# PAIRED T-TEST
# ==============
"""
    paired_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

Compute paired t-test with degrees of freedom, t-statistic, and p-value.

# Arguments
- `x::AbstractVector`: First condition data (must have same length as y)
- `y::AbstractVector`: Second condition data (must have same length as x)
- `tail::Symbol`: Type of test - `:both` (two-tailed, default), `:left` (one-tailed, A < B), or `:right` (one-tailed, A > B)

# Returns
- `TTestResult`: Struct containing `df` (degrees of freedom), `t` (t-statistic), and `p` (p-value).
  Returns `NaN` for p-value if t-value is `NaN` or `Inf`.
"""
function paired_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

    # validate equal lengths
    length(x) == length(y) || error("Paired t-test requires equal sample sizes")
    
    n = length(x)
    n < 2 && return (df = NaN, t = NaN, p = NaN)  # Need at least 2 observations
    
    df = n - 1  # degrees of freedom for paired t-test
    
    # Compute mean and std of differences 
    diff = x .- y
    mean_diff = mean(diff)
    std_diff = std(diff, corrected = true)  # Sample standard deviation (n-1)
    
    # Compute t-value
    if std_diff == 0.0
        if mean_diff == 0.0
            t = NaN
        else
            t = Inf * sign(mean_diff)
        end
    else
        t = mean_diff / (std_diff / sqrt(n))
    end
    
    # Handle edge cases
    isnan(t) || isinf(t) && return (df = df, t = t, p = NaN)
    
    # Compute p-value
    dist = TDist(df)
    if tail == :both # Two-tailed test
        p = 2 * (1 - cdf(dist, abs(t)))
    elseif tail == :left # One-tailed: H0: A >= B, H1: A < B
        p = cdf(dist, t)
    elseif tail == :right # One-tailed: H0: A <= B, H1: A > B
        p = 1 - cdf(dist, t)
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return (df = df, t = t, p = p)
end

# ==================
# INDEPENDENT T-TEST
# ==================
"""
    independent_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

Compute independent sample t-test with degrees of freedom, t-statistic, and p-value.
Assumes equal variances (standard independent t-test).

# Arguments
- `x::AbstractVector`: First group data (must have at least 2 observations)
- `y::AbstractVector`: Second group data (must have at least 2 observations)
- `tail::Symbol`: Type of test - `:both` (two-tailed, default), `:left` (one-tailed, A < B), or `:right` (one-tailed, A > B)

# Returns
- `TTestResult`: Struct containing `df` (degrees of freedom), `t` (t-statistic), and `p` (p-value).
  Returns `NaN` for p-value if t-value is `NaN` or `Inf`.
"""
function independent_ttest(x::AbstractVector, y::AbstractVector; tail::Symbol = :both)

    # Validate input lengths
    n_A, n_B = length(x), length(y)
    n_A < 2 || n_B < 2 && error("Independent t-test requires at least 2 observations per group")
    
    df = n_A + n_B - 2  # degrees of freedom for independent t-test
    
    # Compute means and variances (Statistics.jl is already optimized with SIMD)
    mean_x, mean_y = mean(x), mean(y)
    
    # Pooled variance (assuming equal variances)
    var_x = var(x, corrected = true)
    var_y = var(y, corrected = true)
    pooled_var = ((n_A - 1) * var_x + (n_B - 1) * var_y) / df
    
    # Compute t-value
    if pooled_var == 0.0
        if mean_x == mean_y
            t = NaN  # Both groups identical
        else
            t = Inf * sign(mean_x - mean_y)
        end
    else
        t = (mean_x - mean_y) / sqrt(pooled_var * (1/n_A + 1/n_B))
    end
    
    # Handle edge cases
    isnan(t) || isinf(t) && return (df = df, t = t, p = NaN)
    
    # Compute p-value
    dist = TDist(df)
    if tail == :both # Two-tailed test
        p = 2 * (1 - cdf(dist, abs(t)))
    elseif tail == :left # One-tailed: H0: A >= B, H1: A < B
        p = cdf(dist, t)
    elseif tail == :right # One-tailed: H0: A <= B, H1: A > B
        p = 1 - cdf(dist, t)
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return (df = df, t = t, p = p)
end

"""
    AnalysisData

Core data structure for statistical test computations.

# Fields
- `design::Symbol`: Design type - `:paired` or `:independent`
- `data::Vector{Array{Float64, 3}}`: Data for conditions 1 and 2 [condition1, condition2], each [participants × electrodes × time]
- `time_points::Vector{Float64}`: Time points in seconds for the analysis window
"""
struct AnalysisData
    design::Symbol                   # :paired vs. :independent
    data::Vector{Array{Float64, 3}}  # [condition1, condition2] - each [participants × electrodes × time]
    time_points::Vector{Float64}     # analysis window
end


"""
    StatisticalTestData

Stores prepared data for statistical tests (both permutation and analytic tests).

# Fields
- `data::Vector{ErpData}`: Grand average ERPs for conditions 1 and 2 (for visualization/storage)
- `analysis::AnalysisData`: Core analysis data (design, data arrays, time points)
"""
struct StatisticalTestData
    data::Vector{ErpData}
    analysis::AnalysisData
end

function Base.show(io::IO, data::StatisticalTestData)
    time_range_str(times) = isempty(times) ? "N/A" : "$(first(times)) to $(last(times)) s"
    
    # Get dimensions
    n_participants1 = size(data.analysis.data[1], 1)
    n_participants2 = size(data.analysis.data[2], 1)
    
    # Get summary info from grand averages
    n_electrodes1 = length(channel_labels(data.data[1]))
    n_electrodes2 = length(channel_labels(data.data[2]))
    time_range1_str = time_range_str(data.data[1].data[!, :time])
    time_range2_str = time_range_str(data.data[2].data[!, :time])
    analysis_time_range = time_range_str(data.analysis.time_points)

    println(io, "StatisticalTestData")
    println(io, "├─ Design: $(data.analysis.design)")
    println(io, "├─ Condition 1 ($(data.data[1].condition_name)): $n_participants1 participants")
    println(io, "│  └─ $(n_electrodes1) channels, $time_range1_str, $(data.data[1].sample_rate) Hz")
    println(io, "├─ Condition 2 ($(data.data[2].condition_name)): $n_participants2 participants")
    println(io, "│  └─ $(n_electrodes2) channels, $time_range2_str, $(data.data[2].sample_rate) Hz")
    println(io, "├─ Analysis time points: $analysis_time_range")
    println(io, "└─ Analysis dimensions: $(size(data.analysis.data[1])) (participants × electrodes × time)")
end


"""
    prepare_statistical_test_data(erps::Vector{ErpData};
                             design::Symbol = :paired,
                             condition_selection::Function = conditions([1, 2]),
                             channel_selection::Function = channels(),
                             sample_selection::Function = samples(),
                             baseline_window::Function = samples(),
                             analysis_window::Function = samples())

    prepare_statistical_test_data(file_pattern::String, design::Symbol;
                             input_dir::String = pwd(),
                             participant_selection::Function = participants(),
                             condition_selection::Function = conditions([1, 2]),
                             channel_selection::Function = channels(),
                             sample_selection::Function = samples(),
                             baseline_window::Function = samples(),
                             analysis_window::Function = samples())

Prepare ErpData for statistical tests (permutation and analytic tests).

Organizes ErpData into participant × electrode × time arrays for statistical analysis.
Validates the design and ensures data consistency across conditions.

# Arguments (Direct data version)
- `erps::Vector{ErpData}`: ERPs containing data for multiple conditions/participants
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `condition_selection::Function`: Predicate to select exactly 2 conditions for comparison (default: `conditions([1, 2])`)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()` - all channels)
- `sample_selection::Function`: Predicate to select time points (default: `samples()` - all samples). Use `samples((start, end))` for time windows.
- `baseline_window::Function`: Baseline window sample selection predicate (default: `samples()` - all samples, baseline skipped). Use `samples((start, end))` for baseline window.
- `analysis_window::Function`: Analysis window sample selection predicate (default: `samples()` - all samples). Use `samples((start, end))` for analysis time windows.

# Arguments (File-based version)
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned")
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `input_dir::String`: Directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Predicate to filter participants (default: `participants()` - all participants)
- `condition_selection::Function`: Predicate to select exactly 2 conditions for comparison (default: `conditions([1, 2])`)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()` - all channels)
- `sample_selection::Function`: Predicate to select time points (default: `samples()` - all samples). Use `samples((start, end))` for time windows.
- `baseline_window::Function`: Baseline window sample selection predicate (default: `samples()` - all samples, baseline skipped). Use `samples((start, end))` for baseline window.
- `analysis_window::Function`: Analysis window sample selection predicate (default: `samples()` - all samples). Use `samples((start, end))` for analysis time windows.

# Returns
- `StatisticalTestData`: Prepared data structure ready for statistical testing
"""

# Direct data version
function prepare_statistical_test_data(
    erps::Vector{ErpData};
    design::Symbol = :paired,
    condition_selection::Function = conditions([1, 2]),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_window::Function = samples(),
    analysis_window::Function = samples(),
)
    # Group all ERPs by condition first
    erps_by_condition = group_by_condition(erps)
    
    # Apply condition selection to the sorted condition numbers
    all_cond_nums = collect(keys(erps_by_condition))  # Already sorted by group_by_condition
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]
    
    # Validate exactly 2 conditions
    length(selected_cond_nums) == 2 || @minimal_error_throw "Statistical tests require exactly 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select exactly 2 conditions."
    
    condition1 = erps_by_condition[selected_cond_nums[1]]
    condition2 = erps_by_condition[selected_cond_nums[2]]
    
    # Validate design
    design in (:paired, :independent) || @minimal_error "design must be :paired or :independent, got :$design"
    
    # Extract participant IDs from filenames (using utility from batch.jl)
    participants1 = [_extract_participant_id(basename(data.file)) for data in condition1]
    participants2 = [_extract_participant_id(basename(data.file)) for data in condition2]
    
    # Validate design
    if design == :paired # Paired design: same participants in both conditions, in the same order
        if participants1 != participants2
            @minimal_error "Paired design requires same participants in both conditions"
        end
    elseif design == :independent # Independent design: different participants (or allow overlap)
        if length(participants1) < 2 || length(participants2) < 2
            @minimal_error "Independent design requires at least 2 participants per group"
        end
    end
    
    # Validate all ERPs have same structure within each condition
    have_same_structure(condition1) || @minimal_error("Condition 1: ERPs have inconsistent structure")
    have_same_structure(condition2) || @minimal_error("Condition 2: ERPs have inconsistent structure")
    have_same_structure(condition1[1], condition2[1]) || @minimal_error("Condition 1 vs. 2: ERPs have inconsistent structure")

    condition1 = subset(
        condition1;
        channel_selection = channel_selection,  
        sample_selection = sample_selection,
    )
    isempty(condition1) && @minimal_error_throw "No data matched the selection criteria!"

    condition2 = subset(
        condition2;
        channel_selection = channel_selection,  
        sample_selection = sample_selection,
    )
    isempty(condition2) && @minimal_error_throw "No data matched the selection criteria!"
   
    # baseline 
    baseline!.(condition1, Ref(baseline_window))
    baseline!.(condition2, Ref(baseline_window))

    # create grand averages for ease of use in plotting results
    condition1_avg = grand_average(condition1, selected_cond_nums[1])
    condition2_avg = grand_average(condition2, selected_cond_nums[2])

    # create second subset with analysis_window for statistical tests
    condition1 = subset(
        condition1;
        channel_selection = channel_selection,
        sample_selection = analysis_window,
    )
    isempty(condition1) && @minimal_error_throw "No data matched the analysis window criteria!"

    condition2 = subset(
        condition2;
        channel_selection = channel_selection,
        sample_selection = analysis_window,
    )
    isempty(condition2) && @minimal_error_throw "No data matched the analysis window criteria!"

    # Get dimensions and metadata from analysis subset
    electrodes = channel_labels(condition1[1])
    n_electrodes = length(electrodes)
    time_points = condition1[1].data[!, :time]
    n_time = length(time_points)
    
    # Extract data arrays: [participants × electrodes × time]
    condition1 = cat([reshape(Matrix(erp.data[!, electrodes])', 1, n_electrodes, n_time) for erp in condition1]..., dims=1)
    condition2 = cat([reshape(Matrix(erp.data[!, electrodes])', 1, n_electrodes, n_time) for erp in condition2]..., dims=1)
    
    return StatisticalTestData(
        [condition1_avg, condition2_avg],
        AnalysisData(design, [condition1, condition2], time_points),
    )

end


# File-based version (convenience wrapper)
function prepare_statistical_test_data(
    file_pattern::String,
    design::Symbol;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions([1, 2]),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_window::Function = samples(),
    analysis_window::Function = samples(),
)
    # Load all ERPs
    all_erps = load_all_data(ErpData, file_pattern; input_dir = input_dir, participant_selection = participant_selection)
    isempty(all_erps) && @minimal_error_throw "No valid ERP data found matching pattern '$file_pattern' in $input_dir"
    
    # Call the main preparation function (handles grouping and condition selection)
    return prepare_statistical_test_data(
        all_erps;
        design = design,
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        baseline_window = baseline_window,
        analysis_window = analysis_window,
    )
end

# ======================================
# CLUSTER PERMUTATION TESTS
# ======================================

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
"""
struct Cluster
    id::Int
    electrodes::Vector{Symbol}
    time_indices::Vector{Int}
    time_range::Tuple{Float64, Float64}
    cluster_stat::Float64
    p_value::Float64
    is_significant::Bool
    polarity::Symbol
end

"""
    ClusterPermutationResult

Stores complete results from a cluster-based permutation test.

# Fields
- `design::Symbol`: Design type - `:paired` or `:independent`
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df::Float64`: Degrees of freedom (constant across all electrodes/time points)
- `significant_mask_positive::BitArray{2}`: Positive significant points [electrodes × time]
- `significant_mask_negative::BitArray{2}`: Negative significant points [electrodes × time]
- `positive_clusters::Vector{Cluster}`: Positive clusters (each Cluster contains its cluster_stat)
- `negative_clusters::Vector{Cluster}`: Negative clusters (each Cluster contains its cluster_stat)
- `n_permutations::Int`: Number of permutations performed
- `permutation_max_positive::Vector{Float64}`: Maximum cluster stats from permutations (positive)
- `permutation_max_negative::Vector{Float64}`: Maximum cluster stats from permutations (negative)
- `random_seed::Union{Int, Nothing}`: Random seed used (if any)
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `threshold::Float64`: Threshold used (alpha level)
- `threshold_method::Symbol`: `:parametric`, `:nonparametric_individual`, or `:nonparametric_common`
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`
- `cluster_statistic::Symbol`: `:sum`, `:max`, `:size`, or `:wcm`
- `critical_t_values::Array{Float64, 2}`: Critical t-values used [electrodes × time]
"""
struct ClusterPermutationResult
    design::Symbol
    t_matrix::Array{Float64, 2}
    df::Float64
    significant_mask_positive::BitArray{2}
    significant_mask_negative::BitArray{2}
    positive_clusters::Vector{Cluster}
    negative_clusters::Vector{Cluster}
    n_permutations::Int
    permutation_max_positive::Vector{Float64}
    permutation_max_negative::Vector{Float64}
    random_seed::Union{Int, Nothing}
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
    threshold::Float64
    threshold_method::Symbol
    cluster_type::Symbol
    cluster_statistic::Symbol
    critical_t_values::Union{Array{Float64, 2}, Tuple{Float64, Float64}, Tuple{Array{Float64, 2}, Array{Float64, 2}}}
end

function Base.show(io::IO, result::ClusterPermutationResult)
    n_electrodes = length(result.electrodes)
    n_time_points = length(result.time_points)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"
    
    n_pos_clusters = length(result.positive_clusters)
    n_neg_clusters = length(result.negative_clusters)
    n_sig_pos = count(c -> c.is_significant, result.positive_clusters)
    n_sig_neg = count(c -> c.is_significant, result.negative_clusters)
    
    # Count significant points
    n_sig_pos_points = count(result.masks.positive)
    n_sig_neg_points = count(result.masks.negative)
    
    println(io, "ClusterPermutationResult")
    println(io, "├─ Design: $(result.design)")
    println(io, "├─ Permutations: $(result.n_permutations)")
    println(io, "├─ Threshold: $(result.threshold) ($(result.threshold_method))")
    println(io, "├─ Cluster type: $(result.cluster_type)")
    println(io, "├─ Cluster statistic: $(result.cluster_statistic)")
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
            sig_clusters = [c for c in result.positive_clusters if c.is_significant]
            
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
            sig_clusters = [c for c in result.negative_clusters if c.is_significant]
            
            for cluster in sig_clusters
                sig_marker = "✓"
                p_str = cluster.p_value < 0.001 ? "<0.001" : Printf.@sprintf("%.3f", cluster.p_value)
                stat_str = Printf.@sprintf("%.2f", cluster.cluster_stat)
                time_str = "$(cluster.time_range[1])-$(cluster.time_range[2]) s"
                n_elec = length(cluster.electrodes)
                println(io, "     [$sig_marker] Cluster $(cluster.id): stat=$stat_str, p=$p_str, $n_elec electrodes, $time_str")
            end
        end
        
        if result.random_seed !== nothing
            println(io, "   Random seed: $(result.random_seed)")
        end
    else
        if result.random_seed !== nothing
            println(io, "└─ Random seed: $(result.random_seed)")
        end
    end
end






# ===================
# PHASE 0: VALIDATION
# ===================

"""
    validate_permutation_inputs(prepared::StatisticalTestData, 
                                n_permutations::Int,
                                threshold::Float64,
                                cluster_type::Symbol,
                                cluster_statistic::Symbol,
                                tail::Symbol)

Validate inputs for cluster permutation test.

# Arguments
- `prepared::StatisticalTestData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `threshold::Float64`: P-value threshold
- `cluster_type::Symbol`: Type of clustering
- `cluster_statistic::Symbol`: Cluster statistic type
- `tail::Symbol`: Test tail

# Throws
- `ArgumentError`: If any validation fails
"""
function validate_permutation_inputs(
    prepared::StatisticalTestData,
    n_permutations::Int,
    threshold::Float64,
    cluster_type::Symbol,
    cluster_statistic::Symbol,
    tail::Symbol
)
    # Validate design
    if prepared.analysis.design ∉ (:paired, :independent)
        error("Design must be :paired or :independent, got :$(prepared.analysis.design)")
    end
    
    # Validate minimum participants
    n_participants1 = size(prepared.analysis.data[1], 1)
    n_participants2 = size(prepared.analysis.data[2], 1)
    if prepared.analysis.design == :paired
        if n_participants1 < 2
            error("Paired design requires at least 2 participants, got $n_participants1")
        end
    else
        if n_participants1 < 2 || n_participants2 < 2
            error("Independent design requires at least 2 participants per group, " *
                  "got $n_participants1 and $n_participants2")
        end
    end
    
    # Validate n_permutations
    if n_permutations < 1
        error("n_permutations must be >= 1, got $n_permutations")
    end
    
    # Validate threshold
    if threshold <= 0.0 || threshold >= 1.0
        error("threshold must be between 0 and 1, got $threshold")
    end
    
    # Validate cluster_type
    if cluster_type ∉ (:spatial, :temporal, :spatiotemporal)
        error("cluster_type must be :spatial, :temporal, or :spatiotemporal, got :$cluster_type")
    end
    
    # Validate cluster_statistic
    if cluster_statistic ∉ (:sum, :max, :size, :wcm)
        error("cluster_statistic must be :sum, :max, :size, or :wcm, got :$cluster_statistic")
    end
    
    # Validate tail
    if tail ∉ (:both, :left, :right)
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points
    layout = prepared.data[1].layout
    
    # Validate Layout.neighbours for spatial clustering
    if cluster_type ∈ (:spatial, :spatiotemporal)
        if isnothing(layout.neighbours)
            @minimal_warning "Layout.neighbours is not set. Computing with default distance criterion (40.0)."
            # Compute neighbours if missing (using default criterion)
            get_layout_neighbours_xy!(layout, 40.0)
        end
    end
    
    return nothing
end

# ===================
# PHASE 1: CORE STATISTICS
# ===================

"""
    compute_t_matrix(prepared::StatisticalTestData; tail::Symbol = :both)

Compute t-statistics and p-values for all electrode × time points.

# Arguments
- `prepared::StatisticalTestData`: Prepared data for statistical test
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`

# Returns
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df::Float64`: Degrees of freedom (constant across all electrodes/time points)
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]

# Examples
```julia
t_matrix, df, p_matrix = compute_t_matrix(prepared, tail=:both)
```
"""
# Optimized version that accepts raw arrays (for permutation loop - avoids StatisticalTestData creation)
# Vectorized computation to avoid function call overhead and allocations
# Can accept pre-allocated buffers to avoid allocations
function compute_t_matrix(
    data1::Array{Float64, 3},
    data2::Array{Float64, 3},
    design::Symbol;
    tail::Symbol = :both,
    diff_buffer::Union{Nothing, Array{Float64, 3}} = nothing,
    mean1_buffer::Union{Nothing, Array{Float64, 2}} = nothing,
    mean2_buffer::Union{Nothing, Array{Float64, 2}} = nothing,
    mean_diff_buffer::Union{Nothing, Array{Float64, 2}} = nothing,
    std_diff_buffer::Union{Nothing, Array{Float64, 2}} = nothing
)
    n_participants, n_electrodes, n_time = size(data1)
    t_matrix = Array{Float64, 2}(undef, n_electrodes, n_time)
    p_matrix = Array{Float64, 2}(undef, n_electrodes, n_time)
    
    if design == :paired
        # Vectorized paired t-test computation without creating intermediate diff array
        # Compute mean and std of differences directly using Statistics functions
        # This avoids allocating the full diff array [participants × electrodes × time]
        
        # Compute mean of differences: mean(data1 - data2) = mean(data1) - mean(data2)
        # Use pre-allocated buffers if provided, otherwise allocate
        if mean1_buffer !== nothing
            mean1 = mean1_buffer
            mean2 = mean2_buffer
            mean_diff = mean_diff_buffer
        else
            mean1 = Array{Float64, 2}(undef, n_electrodes, n_time)
            mean2 = Array{Float64, 2}(undef, n_electrodes, n_time)
            mean_diff = Array{Float64, 2}(undef, n_electrodes, n_time)
        end
        
        # Compute means manually to avoid allocations from mean(..., dims=1)
        for e_idx in 1:n_electrodes
            for t_idx in 1:n_time
                sum1 = 0.0
                sum2 = 0.0
                for p_idx in 1:n_participants
                    sum1 += data1[p_idx, e_idx, t_idx]
                    sum2 += data2[p_idx, e_idx, t_idx]
                end
                mean1[e_idx, t_idx] = sum1 / n_participants
                mean2[e_idx, t_idx] = sum2 / n_participants
                mean_diff[e_idx, t_idx] = mean1[e_idx, t_idx] - mean2[e_idx, t_idx]
            end
        end
        
        # Compute std of differences: need to compute variance of (data1 - data2)
        # var(data1 - data2) = var(data1) + var(data2) - 2*cov(data1, data2)
        # For paired data, we can compute this more efficiently
        # But std(diff) requires the diff array... so we need to compute it differently
        
        # Alternative: compute variance directly without full diff array
        # var(diff) = mean((diff - mean_diff)^2) * n/(n-1)
        # We can compute this by iterating or using a more efficient approach
        
        # Compute std of differences using pre-allocated buffer if provided
        if diff_buffer !== nothing
            diff = diff_buffer
        else
            diff = Array{Float64, 3}(undef, n_participants, n_electrodes, n_time)
        end
        # Compute diff array (subtraction is fast)
        diff .= data1 .- data2
        
        # Compute std manually to avoid allocations from std(..., dims=1)
        if std_diff_buffer !== nothing
            std_diff = std_diff_buffer
        else
            std_diff = Array{Float64, 2}(undef, n_electrodes, n_time)
        end
        
        # Compute std of differences manually (avoids allocations from Statistics.jl)
        for e_idx in 1:n_electrodes
            for t_idx in 1:n_time
                mean_d = mean_diff[e_idx, t_idx]
                sum_sq_diff = 0.0
                for p_idx in 1:n_participants
                    diff_val = diff[p_idx, e_idx, t_idx]
                    sum_sq_diff += (diff_val - mean_d)^2
                end
                std_diff[e_idx, t_idx] = sqrt(sum_sq_diff / (n_participants - 1))
            end
        end
        
        # Compute t-statistics: t = mean_diff / (std_diff / sqrt(n))
        # Handle division by zero
        zero_std_mask = std_diff .== 0.0
        zero_mean_mask = mean_diff .== 0.0
        
        t_matrix = mean_diff ./ (std_diff ./ sqrt(n_participants))
        # Where std is zero: NaN if mean is also zero, Inf otherwise
        t_matrix[zero_std_mask .& zero_mean_mask] .= NaN
        t_matrix[zero_std_mask .& .!zero_mean_mask] .= Inf
        
        # Degrees of freedom (same for all points in paired design)
        df = Float64(n_participants - 1)
        
        # Compute p-values from t-values and df
        # Note: For permutation loop, we could skip this, but it's needed for observed data
        # The overhead is minimal compared to the t-statistic computation
        dist = TDist(df)
        
        if tail == :both
            # Two-tailed: p = 2 * (1 - cdf(|t|))
            for i in eachindex(t_matrix)
                t_val = t_matrix[i]
                if isnan(t_val) || isinf(t_val)
                    p_matrix[i] = NaN
                else
                    p_matrix[i] = 2 * (1 - cdf(dist, abs(t_val)))
                end
            end
        elseif tail == :right
            # One-tailed right: p = 1 - cdf(t)
            for i in eachindex(t_matrix)
                t_val = t_matrix[i]
                if isnan(t_val) || isinf(t_val)
                    p_matrix[i] = NaN
                else
                    p_matrix[i] = 1 - cdf(dist, t_val)
                end
            end
        elseif tail == :left
            # One-tailed left: p = cdf(t)
            for i in eachindex(t_matrix)
                t_val = t_matrix[i]
                if isnan(t_val) || isinf(t_val)
                    p_matrix[i] = NaN
                else
                    p_matrix[i] = cdf(dist, t_val)
                end
            end
        end
        
    else
        # Independent design: need to loop (but df is constant across all points)
        # Compute df once from first electrode/time point
        data_A_first = view(data1, :, 1, 1)
        data_B_first = view(data2, :, 1, 1)
        result_first = independent_ttest(data_A_first, data_B_first, tail=tail)
        df = result_first.df
        
        for e_idx in 1:n_electrodes
            for t_idx in 1:n_time
                data_A = view(data1, :, e_idx, t_idx)
                data_B = view(data2, :, e_idx, t_idx)
                result = independent_ttest(data_A, data_B, tail=tail)
                t_matrix[e_idx, t_idx] = result.t
                p_matrix[e_idx, t_idx] = result.p
            end
        end
    end
    
    return t_matrix, df, p_matrix
end

# Original version for backward compatibility
function compute_t_matrix(prepared::StatisticalTestData; tail::Symbol = :both)
    return compute_t_matrix(prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, tail=tail)
end

"""
    compute_critical_t_values(df::Float64, 
                              matrix_size::Tuple{Int, Int},
                              alpha::Float64 = 0.05, 
                              tail::Symbol = :both)

Compute critical t-values for parametric thresholding.

# Arguments
- `df::Float64`: Degrees of freedom (constant across all electrodes/time points)
- `matrix_size::Tuple{Int, Int}`: Size of the t_matrix (electrodes × time)
- `alpha::Float64`: Significance level (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time] (all values are the same)

# Examples
```julia
critical_t = compute_critical_t_values(df, size(t_matrix), 0.05, :both)
```
"""
function compute_critical_t_values(
    df::Float64,
    matrix_size::Tuple{Int, Int},
    alpha::Float64 = 0.05,
    tail::Symbol = :both
)
    critical_t_values = Array{Float64, 2}(undef, matrix_size)
    
    if isnan(df) || isinf(df) || df <= 0
        fill!(critical_t_values, NaN)
        return critical_t_values
    end
    
    dist = TDist(df)
    if tail == :both
        alpha_per_tail = alpha / 2.0
        crit_t = quantile(dist, 1.0 - alpha_per_tail)
        fill!(critical_t_values, crit_t)
    elseif tail == :right
        crit_t = quantile(dist, 1.0 - alpha)
        fill!(critical_t_values, crit_t)
    elseif tail == :left
        crit_t = quantile(dist, alpha)
        fill!(critical_t_values, crit_t)
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return critical_t_values
end

# ===================
# PHASE 2: THRESHOLDING
# ===================

"""
    threshold_t_matrix_parametric(t_matrix::Array{Float64, 2},
                                   critical_t_values::Array{Float64, 2},
                                   tail::Symbol = :both)

Threshold t-matrix using parametric critical values.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time]
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `mask_positive::BitArray{2}`: Mask for positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Mask for negative significant points [electrodes × time]

# Examples
```julia
mask_pos, mask_neg = threshold_t_matrix_parametric(t_matrix, critical_t, :both)
```
"""
function threshold_t_matrix_parametric(
    t_matrix::Array{Float64, 2},
    critical_t_values::Array{Float64, 2},
    tail::Symbol = :both
)
    n_electrodes, n_time = size(t_matrix)
    
    if tail == :both
        # Two-tailed: significant if |t| > critical_t
        mask_positive = BitArray{2}(undef, n_electrodes, n_time)
        mask_negative = BitArray{2}(undef, n_electrodes, n_time)
        
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            
            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = t_val > crit_t
                mask_negative[i] = t_val < -crit_t
            end
        end
        
        return mask_positive, mask_negative
        
    elseif tail == :right
        # One-tailed right: significant if t > critical_t
        mask_positive = BitArray{2}(undef, n_electrodes, n_time)
        mask_negative = falses(n_electrodes, n_time)
        
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            
            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
            else
                mask_positive[i] = t_val > crit_t
            end
        end
        
        return mask_positive, mask_negative
        
    elseif tail == :left
        # One-tailed left: significant if t < critical_t
        mask_positive = falses(n_electrodes, n_time)
        mask_negative = BitArray{2}(undef, n_electrodes, n_time)
        
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            
            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_negative[i] = false
            else
                mask_negative[i] = t_val < crit_t
            end
        end
        
        return mask_positive, mask_negative
        
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

# In-place version for performance (fills pre-allocated masks)
function threshold_t_matrix_parametric!(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    t_matrix::Array{Float64, 2},
    critical_t_values::Array{Float64, 2},
    tail::Symbol = :both
)
    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            
            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                mask_positive[i] = t_val > crit_t
                mask_negative[i] = t_val < -crit_t
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            
            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_positive[i] = false
            else
                mask_positive[i] = t_val > crit_t
            end
        end
    elseif tail == :left
        fill!(mask_positive, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            crit_t = critical_t_values[i]
            
            if isnan(t_val) || isnan(crit_t) || isinf(t_val) || isinf(crit_t)
                mask_negative[i] = false
            else
                mask_negative[i] = t_val < crit_t
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

"""
    compute_nonparametric_threshold_common(permutation_t_matrices::Array{Float64, 3},
                                          alpha::Float64 = 0.05,
                                          tail::Symbol = :both)

Compute a single non-parametric threshold from pooled permutation distribution (common threshold).

# Arguments
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]
- `alpha::Float64`: Significance level (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `threshold_positive::Float64`: Threshold for positive t-values
- `threshold_negative::Float64`: Threshold for negative t-values (absolute value)

# Examples
```julia
thresh_pos, thresh_neg = compute_nonparametric_threshold_common(perm_t_matrices, 0.05, :both)
```
"""
function compute_nonparametric_threshold_common(
    permutation_t_matrices::Array{Float64, 3},
    alpha::Float64 = 0.05,
    tail::Symbol = :both
)
    n_electrodes, n_time, n_permutations = size(permutation_t_matrices)
    
    if tail == :both
        # Two-tailed: collect all absolute t-values
        all_t_values = Float64[]
        sizehint!(all_t_values, n_electrodes * n_time * n_permutations)
        
        for perm_idx in 1:n_permutations
            for i in 1:n_electrodes
                for j in 1:n_time
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val)
                        push!(all_t_values, abs(t_val))
                    end
                end
            end
        end
        
        if isempty(all_t_values)
            error("No valid t-values found in permutation distribution")
        end
        
        # Compute (1 - alpha/2) percentile for two-tailed
        percentile_level = 1.0 - (alpha / 2.0)
        threshold = quantile(all_t_values, percentile_level)
        
        return threshold, threshold
        
    elseif tail == :right
        # One-tailed right: collect all positive t-values
        all_t_values = Float64[]
        sizehint!(all_t_values, n_electrodes * n_time * n_permutations)
        
        for perm_idx in 1:n_permutations
            for i in 1:n_electrodes
                for j in 1:n_time
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val > 0
                        push!(all_t_values, t_val)
                    end
                end
            end
        end
        
        if isempty(all_t_values)
            error("No valid positive t-values found in permutation distribution")
        end
        
        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)
        
        return threshold, NaN
        
    elseif tail == :left
        # One-tailed left: collect all negative t-values (absolute values)
        all_t_values = Float64[]
        sizehint!(all_t_values, n_electrodes * n_time * n_permutations)
        
        for perm_idx in 1:n_permutations
            for i in 1:n_electrodes
                for j in 1:n_time
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val < 0
                        push!(all_t_values, abs(t_val))
                    end
                end
            end
        end
        
        if isempty(all_t_values)
            error("No valid negative t-values found in permutation distribution")
        end
        
        percentile_level = 1.0 - alpha
        threshold = quantile(all_t_values, percentile_level)
        
        return NaN, threshold
        
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

"""
    compute_nonparametric_threshold_individual(permutation_t_matrices::Array{Float64, 3},
                                               alpha::Float64 = 0.05,
                                               tail::Symbol = :both)

Compute point-specific non-parametric thresholds from permutation distribution (individual thresholds).

# Arguments
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]
- `alpha::Float64`: Significance level (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `thresholds_positive::Array{Float64, 2}`: Point-specific thresholds for positive t-values [electrodes × time]
- `thresholds_negative::Array{Float64, 2}`: Point-specific thresholds for negative t-values [electrodes × time]

# Examples
```julia
thresh_pos, thresh_neg = compute_nonparametric_threshold_individual(perm_t_matrices, 0.05, :both)
```
"""
function compute_nonparametric_threshold_individual(
    permutation_t_matrices::Array{Float64, 3},
    alpha::Float64 = 0.05,
    tail::Symbol = :both
)
    n_electrodes, n_time, n_permutations = size(permutation_t_matrices)
    
    thresholds_positive = Array{Float64, 2}(undef, n_electrodes, n_time)
    thresholds_negative = Array{Float64, 2}(undef, n_electrodes, n_time)
    
    if tail == :both
        # Two-tailed: for each point, compute (1 - alpha/2) percentile
        percentile_level = 1.0 - (alpha / 2.0)
        
        for i in 1:n_electrodes
            for j in 1:n_time
                # Collect t-values at this point across all permutations
                t_values = Float64[]
                sizehint!(t_values, n_permutations)
                
                for perm_idx in 1:n_permutations
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val)
                        push!(t_values, abs(t_val))
                    end
                end
                
                if isempty(t_values)
                    thresholds_positive[i, j] = NaN
                    thresholds_negative[i, j] = NaN
                else
                    threshold = quantile(t_values, percentile_level)
                    thresholds_positive[i, j] = threshold
                    thresholds_negative[i, j] = threshold
                end
            end
        end
        
    elseif tail == :right
        # One-tailed right: for each point, compute (1 - alpha) percentile of positive values
        percentile_level = 1.0 - alpha
        
        for i in 1:n_electrodes
            for j in 1:n_time
                # Collect positive t-values at this point
                t_values = Float64[]
                sizehint!(t_values, n_permutations)
                
                for perm_idx in 1:n_permutations
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val > 0
                        push!(t_values, t_val)
                    end
                end
                
                if isempty(t_values)
                    thresholds_positive[i, j] = NaN
                else
                    thresholds_positive[i, j] = quantile(t_values, percentile_level)
                end
                thresholds_negative[i, j] = NaN
            end
        end
        
    elseif tail == :left
        # One-tailed left: for each point, compute (1 - alpha) percentile of negative values (absolute)
        percentile_level = 1.0 - alpha
        
        for i in 1:n_electrodes
            for j in 1:n_time
                # Collect negative t-values at this point (absolute values)
                t_values = Float64[]
                sizehint!(t_values, n_permutations)
                
                for perm_idx in 1:n_permutations
                    t_val = permutation_t_matrices[i, j, perm_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val < 0
                        push!(t_values, abs(t_val))
                    end
                end
                
                thresholds_positive[i, j] = NaN
                if isempty(t_values)
                    thresholds_negative[i, j] = NaN
                else
                    thresholds_negative[i, j] = quantile(t_values, percentile_level)
                end
            end
        end
        
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return thresholds_positive, thresholds_negative
end

"""
    threshold_t_matrix_nonparametric(t_matrix::Array{Float64, 2},
                                     thresholds_positive::Union{Float64, Array{Float64, 2}},
                                     thresholds_negative::Union{Float64, Array{Float64, 2}},
                                     tail::Symbol = :both)

Threshold t-matrix using non-parametric thresholds.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `thresholds_positive::Union{Float64, Array{Float64, 2}}`: Threshold(s) for positive t-values (scalar for common, matrix for individual)
- `thresholds_negative::Union{Float64, Array{Float64, 2}}`: Threshold(s) for negative t-values (scalar for common, matrix for individual)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `mask_positive::BitArray{2}`: Mask for positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Mask for negative significant points [electrodes × time]

# Examples
```julia
# Common threshold (scalar)
mask_pos, mask_neg = threshold_t_matrix_nonparametric(t_matrix, thresh_pos, thresh_neg, :both)

# Individual thresholds (matrices)
mask_pos, mask_neg = threshold_t_matrix_nonparametric(t_matrix, thresh_pos_mat, thresh_neg_mat, :both)
```
"""
function threshold_t_matrix_nonparametric(
    t_matrix::Array{Float64, 2},
    thresholds_positive::Union{Float64, Array{Float64, 2}},
    thresholds_negative::Union{Float64, Array{Float64, 2}},
    tail::Symbol = :both
)
    n_electrodes, n_time = size(t_matrix)
    mask_positive = BitArray{2}(undef, n_electrodes, n_time)
    mask_negative = BitArray{2}(undef, n_electrodes, n_time)
    
    # Determine if thresholds are scalar (common) or matrix (individual)
    is_common = isa(thresholds_positive, Float64)
    
    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            
            if isnan(t_val) || isinf(t_val)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                # Get threshold for this point
                if is_common
                    thresh_pos = thresholds_positive::Float64
                    thresh_neg = thresholds_negative::Float64
                else
                    thresh_pos = thresholds_positive[i]
                    thresh_neg = thresholds_negative[i]
                end
                
                if isnan(thresh_pos) || isnan(thresh_neg)
                    mask_positive[i] = false
                    mask_negative[i] = false
                else
                    mask_positive[i] = t_val > thresh_pos
                    mask_negative[i] = t_val < -thresh_neg
                end
            end
        end
        
    elseif tail == :right
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            
            if isnan(t_val) || isinf(t_val)
                mask_positive[i] = false
            else
                if is_common
                    thresh_pos = thresholds_positive::Float64
                else
                    thresh_pos = thresholds_positive[i]
                end
                
                if isnan(thresh_pos)
                    mask_positive[i] = false
                else
                    mask_positive[i] = t_val > thresh_pos
                end
            end
            mask_negative[i] = false
        end
        
    elseif tail == :left
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            
            if isnan(t_val) || isinf(t_val)
                mask_negative[i] = false
            else
                if is_common
                    thresh_neg = thresholds_negative::Float64
                else
                    thresh_neg = thresholds_negative[i]
                end
                
                if isnan(thresh_neg)
                    mask_negative[i] = false
                else
                    mask_negative[i] = t_val < -thresh_neg
                end
            end
            mask_positive[i] = false
        end
        
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return mask_positive, mask_negative
end

# In-place version for performance (fills pre-allocated masks)
function threshold_t_matrix_nonparametric!(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    t_matrix::Array{Float64, 2},
    thresholds_positive::Union{Float64, Array{Float64, 2}},
    thresholds_negative::Union{Float64, Array{Float64, 2}},
    tail::Symbol = :both
)
    # Determine if thresholds are scalar (common) or matrix (individual)
    is_common = isa(thresholds_positive, Float64)
    
    if tail == :both
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            
            if isnan(t_val) || isinf(t_val)
                mask_positive[i] = false
                mask_negative[i] = false
            else
                # Get threshold for this point
                if is_common
                    thresh_pos = thresholds_positive::Float64
                    thresh_neg = thresholds_negative::Float64
                else
                    thresh_pos = thresholds_positive[i]
                    thresh_neg = thresholds_negative[i]
                end
                
                if isnan(thresh_pos) || isnan(thresh_neg)
                    mask_positive[i] = false
                    mask_negative[i] = false
                else
                    mask_positive[i] = t_val > thresh_pos
                    mask_negative[i] = t_val < -thresh_neg
                end
            end
        end
    elseif tail == :right
        fill!(mask_negative, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            
            if isnan(t_val) || isinf(t_val)
                mask_positive[i] = false
            else
                thresh_pos = is_common ? thresholds_positive::Float64 : thresholds_positive[i]
                if isnan(thresh_pos)
                    mask_positive[i] = false
                else
                    mask_positive[i] = t_val > thresh_pos
                end
            end
        end
    elseif tail == :left
        fill!(mask_positive, false)
        for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            
            if isnan(t_val) || isinf(t_val)
                mask_negative[i] = false
            else
                thresh_neg = is_common ? thresholds_negative::Float64 : thresholds_negative[i]
                if isnan(thresh_neg)
                    mask_negative[i] = false
                else
                    mask_negative[i] = t_val < -thresh_neg
                end
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
end

# ===================
# PHASE 3: CONNECTIVITY MATRIX
# ===================

"""
    build_connectivity_matrix(electrodes::Vector{Symbol},
                              layout::Layout,
                              cluster_type::Symbol)

Build connectivity matrix for clustering.

# Arguments
- `electrodes::Vector{Symbol}`: Electrode labels
- `layout::Layout`: Layout with neighbours information
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`

# Returns
- `connectivity::SparseMatrixCSC{Bool}`: Connectivity matrix
- `n_electrodes::Int`: Number of electrodes
- `n_time::Int`: Number of time points (1 for spatial, actual for temporal/spatiotemporal)

# Notes
For spatiotemporal clustering, the connectivity matrix is for [electrodes × time] space.
For spatial only, it's just [electrodes × electrodes].
For temporal only, it's just consecutive time points.

# Examples
```julia
conn, n_elec, n_time = build_connectivity_matrix(electrodes, layout, :spatiotemporal)
```
"""
function build_connectivity_matrix(
    electrodes::Vector{Symbol},
    layout::Layout,
    cluster_type::Symbol
)
    n_electrodes = length(electrodes)
    
    if cluster_type == :spatial
        # Spatial only: [electrodes × electrodes] matrix
        # Build from Layout.neighbours
        if isnothing(layout.neighbours)
            error("Layout.neighbours is required for spatial clustering. " *
                  "Compute neighbours first using get_layout_neighbours_xy! or get_layout_neighbours_xyz!.")
        end
        
        # Build adjacency matrix
        I = Int[]
        J = Int[]
        
        for (e_idx, electrode) in enumerate(electrodes)
            # Note: No self-connections - FieldTrip doesn't include them in clustering
            # Self-connections are only used for pre-filtering (minNumChannels)
            
            # Get neighbours from layout
            if haskey(layout.neighbours, electrode)
                neighbours = layout.neighbours[electrode]
                for neighbour in neighbours.channels
                    # Find index of neighbour in electrodes list
                    n_idx = findfirst(==(neighbour), electrodes)
                    if n_idx !== nothing
                        push!(I, e_idx)
                        push!(J, n_idx)
                        # Also add reverse (undirected graph)
                        push!(I, n_idx)
                        push!(J, e_idx)
                    end
                end
            end
        end
        
        # Create sparse matrix
        connectivity = sparse(I, J, true, n_electrodes, n_electrodes)
        return connectivity, n_electrodes, 1
        
    elseif cluster_type == :temporal
        # Temporal only: consecutive time points
        # This will be handled differently - we'll build it per time dimension
        # For now, return identity (will be handled in cluster finding)
        n_time = 1  # Placeholder, actual time dimension handled separately
        connectivity = sparse(Bool[1], Bool[1], Bool[true], 1, 1)
        return connectivity, n_electrodes, n_time
        
    elseif cluster_type == :spatiotemporal
        # Spatiotemporal: combine spatial and temporal
        # Matrix size: [electrodes × time] × [electrodes × time]
        # This is complex - we'll build it as needed in cluster finding
        # For now, return spatial connectivity and note that temporal is handled separately
        if isnothing(layout.neighbours)
            error("Layout.neighbours is required for spatiotemporal clustering. " *
                  "Compute neighbours first using get_layout_neighbours_xy! or get_layout_neighbours_xyz!.")
        end
        
        # Build spatial adjacency (same as spatial case)
        I = Int[]
        J = Int[]
        
        for (e_idx, electrode) in enumerate(electrodes)
            # Note: No self-connections - FieldTrip doesn't include them in clustering
            # Self-connections are only used for pre-filtering (minNumChannels)
            
            if haskey(layout.neighbours, electrode)
                neighbours = layout.neighbours[electrode]
                for neighbour in neighbours.channels
                    n_idx = findfirst(==(neighbour), electrodes)
                    if n_idx !== nothing
                        push!(I, e_idx)
                        push!(J, n_idx)
                        push!(I, n_idx)
                        push!(J, e_idx)
                    end
                end
            end
        end
        
        spatial_connectivity = sparse(I, J, true, n_electrodes, n_electrodes)
        return spatial_connectivity, n_electrodes, 1  # Time dimension handled separately
        
    else
        error("cluster_type must be :spatial, :temporal, or :spatiotemporal, got :$cluster_type")
    end
end

# ===================
# PHASE 3.5: PRE-FILTERING (FieldTrip minNumChannels)
# ===================

"""
    prefilter_mask_by_neighbors(mask::BitArray{2},
                                spatial_connectivity::SparseMatrixCSC{Bool},
                                min_num_neighbors::Int)

Pre-filter mask to remove isolated points (FieldTrip's minNumChannels approach).

For each significant point, count how many neighboring significant channels it has.
If a point has fewer than `min_num_neighbors` neighbors, remove it.
This is done iteratively until no more points are removed.

# Arguments
- `mask::BitArray{2}`: Significant points mask [electrodes × time]
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix [electrodes × electrodes]
- `min_num_neighbors::Int`: Minimum number of neighboring significant channels required

# Returns
- `filtered_mask::BitArray{2}`: Filtered mask with isolated points removed

# Examples
```julia
filtered_mask = prefilter_mask_by_neighbors(mask, spatial_connectivity, 3)
```
"""
function prefilter_mask_by_neighbors(
    mask::BitArray{2},
    spatial_connectivity::SparseMatrixCSC{Bool},
    min_num_neighbors::Int
)
    if min_num_neighbors <= 0
        return mask  # No filtering if min_num_neighbors is 0 or negative
    end
    
    n_electrodes, n_time = size(mask)
    filtered_mask = copy(mask)
    
    # Make spatial connectivity symmetric (as FieldTrip does)
    # spatial_connectivity should already be symmetric, but ensure it
    spatial_conn_sym = spatial_connectivity .| spatial_connectivity'
    
    # Iterative removal until no more points are removed
    # FieldTrip approach: for each (time, frequency) element, count how many
    # neighboring significant channels it has. If fewer than min_num_neighbors, remove it.
    n_removed = 1
    while n_removed > 0
        n_removed = 0
        
        # For each time point
        for t_idx in 1:n_time
            # Count neighbors for each electrode at this time point
            for e_idx in 1:n_electrodes
                if !filtered_mask[e_idx, t_idx]
                    continue  # Skip non-significant points
                end
                
                # Count how many neighboring significant channels this point has
                # FieldTrip's minNumChannels includes the channel itself in the count
                # Note: spatial_conn_sym does NOT include self-connections (removed for clustering)
                # So we need to explicitly count self
                neighbor_count = filtered_mask[e_idx, t_idx] ? 1 : 0  # Count self if significant
                for n_e_idx in 1:n_electrodes
                    if e_idx != n_e_idx && spatial_conn_sym[e_idx, n_e_idx] && filtered_mask[n_e_idx, t_idx]
                        neighbor_count += 1
                    end
                end
                
                # Remove point if total neighbor count (including self) is less than min_num_neighbors
                # This matches FieldTrip: (onoff.*nsigneighb) < minnbchan
                if neighbor_count < min_num_neighbors
                    filtered_mask[e_idx, t_idx] = false
                    n_removed += 1
                end
            end
        end
    end
    
    return filtered_mask
end

# In-place version for performance (modifies mask directly)
function prefilter_mask_by_neighbors!(
    mask::BitArray{2},
    spatial_connectivity::SparseMatrixCSC{Bool},
    min_num_neighbors::Int
)
    if min_num_neighbors <= 0
        return  # No filtering if min_num_neighbors is 0 or negative
    end
    
    n_electrodes, n_time = size(mask)
    
    # Make spatial connectivity symmetric (as FieldTrip does)
    spatial_conn_sym = spatial_connectivity .| spatial_connectivity'
    
    # Iterative removal until no more points are removed
    n_removed = 1
    while n_removed > 0
        n_removed = 0
        
        # For each time point
        for t_idx in 1:n_time
            # Count neighbors for each electrode at this time point
            for e_idx in 1:n_electrodes
                if !mask[e_idx, t_idx]
                    continue  # Skip non-significant points
                end
                
                # Count how many neighboring significant channels this point has
                # FieldTrip's minNumChannels includes the channel itself in the count
                # Note: spatial_conn_sym does NOT include self-connections (removed for clustering)
                # So we need to explicitly count self
                neighbor_count = mask[e_idx, t_idx] ? 1 : 0  # Count self if significant
                for n_e_idx in 1:n_electrodes
                    if e_idx != n_e_idx && spatial_conn_sym[e_idx, n_e_idx] && mask[n_e_idx, t_idx]
                        neighbor_count += 1
                    end
                end
                
                # Remove point if total neighbor count is less than min_num_neighbors
                if neighbor_count < min_num_neighbors
                    mask[e_idx, t_idx] = false
                    n_removed += 1
                end
            end
        end
    end
end

# ===================
# PHASE 4: CLUSTER FINDING
# ===================

"""
    find_clusters_connected_components(mask::BitArray{2},
                                       electrodes::Vector{Symbol},
                                       time_points::Vector{Float64},
                                       spatial_connectivity::SparseMatrixCSC{Bool},
                                       cluster_type::Symbol)

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
clusters = find_clusters_connected_components(mask, electrodes, time_points, conn, :spatiotemporal)
```
"""
function find_clusters_connected_components(
    mask::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol
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
    
    # Helper function to get neighbors of a point
    function get_neighbors(e_idx::Int, t_idx::Int)
        neighbors = Tuple{Int, Int}[]
        
        if cluster_type == :spatial
            # Spatial: neighbors are spatially adjacent electrodes at same time
            for n_e_idx in 1:n_electrodes
                if spatial_connectivity[e_idx, n_e_idx] && mask[n_e_idx, t_idx]
                    push!(neighbors, (n_e_idx, t_idx))
                end
            end
            
        elseif cluster_type == :temporal
            # Temporal: neighbors are same electrode at adjacent time points
            if t_idx > 1 && mask[e_idx, t_idx - 1]
                push!(neighbors, (e_idx, t_idx - 1))
            end
            if t_idx < n_time && mask[e_idx, t_idx + 1]
                push!(neighbors, (e_idx, t_idx + 1))
            end
            
        elseif cluster_type == :spatiotemporal
            # Spatiotemporal: both spatial and temporal neighbors
            # Spatial neighbors at same time
            for n_e_idx in 1:n_electrodes
                if spatial_connectivity[e_idx, n_e_idx] && mask[n_e_idx, t_idx]
                    push!(neighbors, (n_e_idx, t_idx))
                end
            end
            # Temporal neighbors at same electrode
            if t_idx > 1 && mask[e_idx, t_idx - 1]
                push!(neighbors, (e_idx, t_idx - 1))
            end
            if t_idx < n_time && mask[e_idx, t_idx + 1]
                push!(neighbors, (e_idx, t_idx + 1))
            end
        end
        
        return neighbors
    end
    
    # BFS to find connected components
    for e_idx in 1:n_electrodes
        for t_idx in 1:n_time
            # Skip if not significant or already labeled
            if !mask[e_idx, t_idx] || cluster_labels[e_idx, t_idx] > 0
                continue
            end
            
            # Start new cluster
            current_cluster_id += 1
            cluster_electrodes = Set{Symbol}()
            cluster_time_indices = Set{Int}()
            
            # BFS queue
            queue = [(e_idx, t_idx)]
            cluster_labels[e_idx, t_idx] = current_cluster_id
            push!(cluster_electrodes, electrodes[e_idx])
            push!(cluster_time_indices, t_idx)
            
            # BFS traversal
            while !isempty(queue)
                current_e, current_t = popfirst!(queue)
                
                # Get neighbors
                neighbors = get_neighbors(current_e, current_t)
                
                for (n_e, n_t) in neighbors
                    if cluster_labels[n_e, n_t] == 0
                        # Not yet visited
                        cluster_labels[n_e, n_t] = current_cluster_id
                        push!(queue, (n_e, n_t))
                        push!(cluster_electrodes, electrodes[n_e])
                        push!(cluster_time_indices, n_t)
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
                :positive  # polarity - will be set by caller (temporary)
            )
            push!(clusters, cluster)
        end
    end
    
    return clusters
end

"""
    find_clusters(mask_positive::BitArray{2},
                  mask_negative::BitArray{2},
                  electrodes::Vector{Symbol},
                  time_points::Vector{Float64},
                  spatial_connectivity::SparseMatrixCSC{Bool},
                  cluster_type::Symbol)

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
pos_clusters, neg_clusters = find_clusters(mask_pos, mask_neg, electrodes, time_points, conn, :spatiotemporal)
```
"""
function find_clusters(
    mask_positive::BitArray{2},
    mask_negative::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol
)
    # Find positive clusters
    positive_clusters_raw = find_clusters_connected_components(
        mask_positive, electrodes, time_points, spatial_connectivity, cluster_type
    )
    # Set polarity to positive (reuse existing Cluster objects, just update polarity)
    positive_clusters = Vector{Cluster}(undef, length(positive_clusters_raw))
    for (i, c) in enumerate(positive_clusters_raw)
        positive_clusters[i] = Cluster(
            c.id, c.electrodes, c.time_indices, c.time_range,
            c.cluster_stat, c.p_value, c.is_significant, :positive
        )
    end
    
    # Find negative clusters
    negative_clusters_raw = find_clusters_connected_components(
        mask_negative, electrodes, time_points, spatial_connectivity, cluster_type
    )
    # Set polarity to negative (reuse existing Cluster objects, just update polarity)
    negative_clusters = Vector{Cluster}(undef, length(negative_clusters_raw))
    for (i, c) in enumerate(negative_clusters_raw)
        negative_clusters[i] = Cluster(
            c.id, c.electrodes, c.time_indices, c.time_range,
            c.cluster_stat, c.p_value, c.is_significant, :negative
        )
    end
    
    return positive_clusters, negative_clusters
end

# ===================
# PHASE 5: CLUSTER STATISTICS
# ===================

"""
    compute_cluster_statistics(clusters::Vector{Cluster},
                                t_matrix::Array{Float64, 2},
                                electrodes::Vector{Symbol},
                                statistic_type::Symbol = :sum)

Compute cluster-level statistics.

# Arguments
- `clusters::Vector{Cluster}`: Clusters to compute statistics for
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels (must match t_matrix row order)
- `statistic_type::Symbol`: Type of statistic - `:sum`, `:max`, `:size`, or `:wcm` (default: `:sum`)

# Returns
- `updated_clusters::Vector{Cluster}`: Clusters with computed statistics
- `cluster_stats::Vector{Float64}`: Vector of cluster statistics

# Examples
```julia
clusters, stats = compute_cluster_statistics(clusters, t_matrix, electrodes, :sum)
```
"""
# Helper function with electrode mapping
function _compute_cluster_statistics(
    clusters::Vector{Cluster},
    t_matrix::Array{Float64, 2},
    electrodes::Vector{Symbol},
    statistic_type::Symbol = :sum
)
    if isempty(clusters)
        return Cluster[], Float64[]
    end
    
    updated_clusters = Cluster[]
    cluster_stats = Float64[]
    
    # Create electrode index lookup
    electrode_to_idx = Dict(e => i for (i, e) in enumerate(electrodes))
    
    for cluster in clusters
        cluster_stat = 0.0
        
        if statistic_type == :sum
            # Sum of t-values within cluster
            for electrode in cluster.electrodes
                e_idx = electrode_to_idx[electrode]
                for t_idx in cluster.time_indices
                    if !isnan(t_matrix[e_idx, t_idx]) && !isinf(t_matrix[e_idx, t_idx])
                        cluster_stat += t_matrix[e_idx, t_idx]
                    end
                end
            end
            
        elseif statistic_type == :max
            # Maximum t-value in cluster
            cluster_stat = -Inf
            for electrode in cluster.electrodes
                e_idx = electrode_to_idx[electrode]
                for t_idx in cluster.time_indices
                    t_val = t_matrix[e_idx, t_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val > cluster_stat
                        cluster_stat = t_val
                    end
                end
            end
            if cluster_stat == -Inf
                cluster_stat = 0.0
            end
            
        elseif statistic_type == :size
            # Number of points in cluster
            cluster_stat = Float64(length(cluster.electrodes) * length(cluster.time_indices))
            
        elseif statistic_type == :wcm
            # Weighted cluster mass (simplified version)
            # WCM = sum(|t|^weight), default weight = 1.0
            weight = 1.0
            for electrode in cluster.electrodes
                e_idx = electrode_to_idx[electrode]
                for t_idx in cluster.time_indices
                    t_val = t_matrix[e_idx, t_idx]
                    if !isnan(t_val) && !isinf(t_val)
                        cluster_stat += abs(t_val)^weight
                    end
                end
            end
        else
            error("statistic_type must be :sum, :max, :size, or :wcm, got :$statistic_type")
        end
        
        # Update cluster with statistic
        updated_cluster = Cluster(
            cluster.id, cluster.electrodes, cluster.time_indices, cluster.time_range,
            cluster_stat, cluster.p_value, cluster.is_significant, cluster.polarity
        )
        push!(updated_clusters, updated_cluster)
        push!(cluster_stats, cluster_stat)
    end
    
    return updated_clusters, cluster_stats
end

# Optimized version that only returns statistics (for permutation loop)
# Avoids creating new Cluster objects
function _compute_cluster_statistics_only(
    clusters::Vector{Cluster},
    t_matrix::Array{Float64, 2},
    electrode_to_idx::Dict{Symbol, Int},
    statistic_type::Symbol = :sum
)
    if isempty(clusters)
        return Float64[]
    end
    
    cluster_stats = Float64[]
    sizehint!(cluster_stats, length(clusters))
    
    for cluster in clusters
        cluster_stat = 0.0
        
        if statistic_type == :sum
            # Sum of t-values within cluster
            for electrode in cluster.electrodes
                e_idx = electrode_to_idx[electrode]
                for t_idx in cluster.time_indices
                    t_val = t_matrix[e_idx, t_idx]
                    if !isnan(t_val) && !isinf(t_val)
                        cluster_stat += t_val
                    end
                end
            end
            
        elseif statistic_type == :max
            # Maximum t-value in cluster
            cluster_stat = -Inf
            for electrode in cluster.electrodes
                e_idx = electrode_to_idx[electrode]
                for t_idx in cluster.time_indices
                    t_val = t_matrix[e_idx, t_idx]
                    if !isnan(t_val) && !isinf(t_val) && t_val > cluster_stat
                        cluster_stat = t_val
                    end
                end
            end
            if cluster_stat == -Inf
                cluster_stat = 0.0
            end
            
        elseif statistic_type == :size
            # Number of points in cluster
            cluster_stat = Float64(length(cluster.electrodes) * length(cluster.time_indices))
            
        elseif statistic_type == :wcm
            # Weighted cluster mass
            weight = 1.0
            for electrode in cluster.electrodes
                e_idx = electrode_to_idx[electrode]
                for t_idx in cluster.time_indices
                    t_val = t_matrix[e_idx, t_idx]
                    if !isnan(t_val) && !isinf(t_val)
                        cluster_stat += abs(t_val)^weight
                    end
                end
            end
        else
            error("statistic_type must be :sum, :max, :size, or :wcm, got :$statistic_type")
        end
        
        push!(cluster_stats, cluster_stat)
    end
    
    return cluster_stats
end

# Public wrapper
function compute_cluster_statistics(
    clusters::Vector{Cluster},
    t_matrix::Array{Float64, 2},
    electrodes::Vector{Symbol},
    statistic_type::Symbol = :sum
)
    return _compute_cluster_statistics(clusters, t_matrix, electrodes, statistic_type)
end

# ===================
# PHASE 6: PERMUTATION LOOP
# ===================

"""
    shuffle_labels(prepared::StatisticalTestData, rng::AbstractRNG = Random.GLOBAL_RNG)

Shuffle condition/group labels for permutation test.

# Arguments
- `prepared::StatisticalTestData`: Prepared data
- `rng::AbstractRNG`: Random number generator (default: global RNG)

# Returns
- `shuffled_data_A::Array{Float64, 3}`: Shuffled data for condition A
- `shuffled_data_B::Array{Float64, 3}`: Shuffled data for condition B

# Notes
For paired design: randomly swap A/B labels within each participant.
For independent design: randomly shuffle participants between groups.
"""
# Optimized version that accepts raw arrays and reuses buffers
# Generate swap mask for paired design (avoids copying arrays)
function generate_swap_mask(n_participants::Int, rng::AbstractRNG = Random.GLOBAL_RNG)
    return BitVector([rand(rng, Bool) for _ in 1:n_participants])
end

function shuffle_labels!(
    shuffled_A::Array{Float64, 3},
    shuffled_B::Array{Float64, 3},
    data1::Array{Float64, 3},
    data2::Array{Float64, 3},
    design::Symbol,
    rng::AbstractRNG = Random.GLOBAL_RNG
)
    if design == :paired
        # Paired: copy data first (copyto! is optimized), then swap slices in-place
        copyto!(shuffled_A, data1)
        copyto!(shuffled_B, data2)
        
        n_participants = size(data1, 1)
        for p_idx in 1:n_participants
            if rand(rng, Bool)
                # Swap conditions for this participant in-place
                shuffled_A[p_idx, :, :], shuffled_B[p_idx, :, :] = 
                    shuffled_B[p_idx, :, :], shuffled_A[p_idx, :, :]
            end
        end
        
    else
        # Independent: randomly shuffle participants between groups
        n_A = size(data1, 1)
        n_B = size(data2, 1)
        n_total = n_A + n_B
        n_electrodes = size(data1, 2)
        n_time = size(data1, 3)
        
        # Use shuffled_A as temporary storage for combined data (needs to be large enough)
        # Ensure shuffled_A is large enough (it should be, but check)
        if size(shuffled_A, 1) < n_total
            error("shuffled_A buffer too small for independent design")
        end
        
        # Copy data1 and data2 into shuffled_A (temporary combined storage)
        copyto!(shuffled_A, 1, data1, 1, n_A * n_electrodes * n_time)
        copyto!(shuffled_A, n_A * n_electrodes * n_time + 1, data2, 1, n_B * n_electrodes * n_time)
        
        # Shuffle indices
        shuffled_indices = collect(1:n_total)
        shuffle!(rng, shuffled_indices)
        
        # Copy shuffled data to output arrays
        # Use shuffled_B as temporary storage to avoid overwriting while reading
        for i in 1:n_total
            src_idx = shuffled_indices[i]
            for e_idx in 1:n_electrodes
                for t_idx in 1:n_time
                    shuffled_B[i, e_idx, t_idx] = shuffled_A[src_idx, e_idx, t_idx]
                end
            end
        end
        
        # Split into groups
        for i in 1:n_A
            for e_idx in 1:n_electrodes
                for t_idx in 1:n_time
                    shuffled_A[i, e_idx, t_idx] = shuffled_B[i, e_idx, t_idx]
                end
            end
        end
        for i in 1:n_B
            for e_idx in 1:n_electrodes
                for t_idx in 1:n_time
                    shuffled_B[i, e_idx, t_idx] = shuffled_B[n_A + i, e_idx, t_idx]
                end
            end
        end
    end
    return shuffled_A, shuffled_B
end

# Original version for backward compatibility
function shuffle_labels(prepared::StatisticalTestData, rng::AbstractRNG = Random.GLOBAL_RNG)
    shuffled_A = similar(prepared.analysis.data[1])
    shuffled_B = similar(prepared.analysis.data[2])
    return shuffle_labels!(shuffled_A, shuffled_B, prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, rng)
end

"""
    collect_permutation_t_matrices(prepared::StatisticalTestData,
                                   n_permutations::Int,
                                   random_seed::Union{Int, Nothing} = nothing,
                                   show_progress::Bool = true)

Run permutations and collect t-matrices for non-parametric thresholding.

# Arguments
- `prepared::StatisticalTestData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `random_seed::Union{Int, Nothing}`: Random seed (default: nothing)
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `permutation_t_matrices::Array{Float64, 3}`: T-statistics from all permutations [electrodes × time × permutations]

# Examples
```julia
perm_t_matrices = collect_permutation_t_matrices(prepared, 1000)
```
"""
function collect_permutation_t_matrices(
    prepared::StatisticalTestData,
    n_permutations::Int,
    random_seed::Union{Int, Nothing} = nothing,
    show_progress::Bool = true
)
    # Set random seed if provided
    if random_seed !== nothing
        rng = MersenneTwister(random_seed)
    else
        rng = Random.GLOBAL_RNG
    end
    
    # Get dimensions from first permutation (using optimized functions)
    shuffled_A_buffer = similar(prepared.analysis.data[1])
    shuffled_B_buffer = similar(prepared.analysis.data[2])
    shuffle_labels!(shuffled_A_buffer, shuffled_B_buffer,
                   prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, rng)
    t_matrix_sample, _, _ = compute_t_matrix(
        shuffled_A_buffer, shuffled_B_buffer, prepared.analysis.design
    )
    n_electrodes, n_time = size(t_matrix_sample)
    
    # Pre-allocate array for all t-matrices
    permutation_t_matrices = Array{Float64, 3}(undef, n_electrodes, n_time, n_permutations)
    
    # Progress bar
    if show_progress
        progress = Progress(n_permutations, desc="Collecting permutation t-matrices: ", showspeed=true)
    end
    
    # Pre-allocate buffers for shuffling (reuse across permutations)
    shuffled_A_buffer = similar(prepared.analysis.data[1])
    shuffled_B_buffer = similar(prepared.analysis.data[2])
    
    for perm_idx in 1:n_permutations
        # Shuffle labels using pre-allocated buffers
        shuffle_labels!(shuffled_A_buffer, shuffled_B_buffer,
                       prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, rng)
        
        # Compute t-matrix directly from arrays (no StatisticalTestData needed)
        t_matrix_perm, _, _ = compute_t_matrix(
            shuffled_A_buffer, shuffled_B_buffer, prepared.analysis.design
        )
        permutation_t_matrices[:, :, perm_idx] = t_matrix_perm
        
        if show_progress
            next!(progress)
        end
    end
    
    return permutation_t_matrices
end

"""
    run_permutations(prepared::StatisticalTestData,
                     n_permutations::Int,
                     threshold::Float64,
                     critical_t_values::Union{Array{Float64, 2}, Tuple{Float64, Float64}, Tuple{Array{Float64, 2}, Array{Float64, 2}}},
                     spatial_connectivity::SparseMatrixCSC{Bool},
                     cluster_type::Symbol,
                     cluster_statistic::Symbol,
                     tail::Symbol,
                     min_num_neighbors::Int,
                     random_seed::Union{Int, Nothing} = nothing,
                     show_progress::Bool = true;
                     permutation_t_matrices::Union{Nothing, Array{Float64, 3}} = nothing)

Run permutation loop to generate distribution of maximum cluster statistics.

# Arguments
- `prepared::StatisticalTestData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `threshold::Float64`: P-value threshold
- `critical_t_values`: Critical t-values - can be:
  - `Array{Float64, 2}` for parametric (common critical t-values)
  - `Tuple{Float64, Float64}` for non-parametric common (positive, negative thresholds)
  - `Tuple{Array{Float64, 2}, Array{Float64, 2}}` for non-parametric individual (positive, negative threshold matrices)
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Connectivity matrix
- `cluster_type::Symbol`: Type of clustering
- `cluster_statistic::Symbol`: Cluster statistic type
- `tail::Symbol`: Test tail
- `min_num_neighbors::Int`: Minimum number of neighbors for pre-filtering
- `random_seed::Union{Int, Nothing}`: Random seed (default: nothing)
- `show_progress::Bool`: Show progress bar (default: true)
- `permutation_t_matrices::Union{Nothing, Array{Float64, 3}}`: Pre-computed t-matrices from permutations (optional, for non-parametric)

# Returns
- `permutation_max_positive::Vector{Float64}`: Max cluster stats from permutations (positive)
- `permutation_max_negative::Vector{Float64}`: Max cluster stats from permutations (negative)

# Examples
```julia
perm_pos, perm_neg = run_permutations(prepared, 1000, 0.05, critical_t, conn, 
                                     :spatiotemporal, :sum, :both, 0, 0)
```
"""
function run_permutations(
    prepared::StatisticalTestData,
    n_permutations::Int,
    threshold::Float64,
    critical_t_values,
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    cluster_statistic::Symbol,
    tail::Symbol,
    min_num_neighbors::Int,
    random_seed::Union{Int, Nothing} = nothing,
    show_progress::Bool = true;
    permutation_t_matrices::Union{Nothing, Array{Float64, 3}} = nothing
)
    # Determine thresholding type from critical_t_values type
    is_parametric = isa(critical_t_values, Array{Float64, 2})
    is_nonparametric_common = isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && 
                               isa(critical_t_values[1], Float64)
    is_nonparametric_individual = isa(critical_t_values, Tuple) && length(critical_t_values) == 2 && 
                                  isa(critical_t_values[1], Array{Float64, 2})
    
    # Set random seed if provided (only if not using pre-computed t-matrices)
    if random_seed !== nothing && permutation_t_matrices === nothing
        rng = MersenneTwister(random_seed)
    else
        rng = Random.GLOBAL_RNG
    end
    
    permutation_max_positive = Float64[]
    permutation_max_negative = Float64[]
    sizehint!(permutation_max_positive, n_permutations)
    sizehint!(permutation_max_negative, n_permutations)
    
    # Pre-allocate buffers (reuse across all permutations)
    shuffled_A_buffer = similar(prepared.analysis.data[1])
    shuffled_B_buffer = similar(prepared.analysis.data[2])
    mask_pos_buffer = BitArray{2}(undef, size(prepared.analysis.data[1], 2), size(prepared.analysis.data[1], 3))
    mask_neg_buffer = BitArray{2}(undef, size(prepared.analysis.data[1], 2), size(prepared.analysis.data[1], 3))
    
    # Pre-allocate all buffers for paired design (reuse across permutations)
    if prepared.analysis.design == :paired
        diff_buffer = similar(prepared.analysis.data[1])
        mean1_buffer = Array{Float64, 2}(undef, size(prepared.analysis.data[1], 2), size(prepared.analysis.data[1], 3))
        mean2_buffer = Array{Float64, 2}(undef, size(prepared.analysis.data[1], 2), size(prepared.analysis.data[1], 3))
        mean_diff_buffer = Array{Float64, 2}(undef, size(prepared.analysis.data[1], 2), size(prepared.analysis.data[1], 3))
        std_diff_buffer = Array{Float64, 2}(undef, size(prepared.analysis.data[1], 2), size(prepared.analysis.data[1], 3))
    else
        diff_buffer = mean1_buffer = mean2_buffer = mean_diff_buffer = std_diff_buffer = nothing
    end
    
    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points
    
    # Pre-allocate electrode lookup (reused across all permutations)
    electrode_to_idx = Dict(e => i for (i, e) in enumerate(electrodes))
    
    # Progress bar
    if show_progress
        progress = Progress(n_permutations, desc="Permutations: ", showspeed=true)
    end
    
    for perm_idx in 1:n_permutations
        # Get t-matrix: either from pre-computed or compute new
        if permutation_t_matrices !== nothing
            t_matrix_perm = permutation_t_matrices[:, :, perm_idx]
        else
            # Shuffle and compute t-matrix
            shuffle_labels!(shuffled_A_buffer, shuffled_B_buffer,
                           prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, rng)
            t_matrix_perm, _, _ = compute_t_matrix(
                shuffled_A_buffer, shuffled_B_buffer, prepared.analysis.design,
                diff_buffer=diff_buffer,
                mean1_buffer=mean1_buffer,
                mean2_buffer=mean2_buffer,
                mean_diff_buffer=mean_diff_buffer,
                std_diff_buffer=std_diff_buffer
            )
        end
        
        # Threshold (fill pre-allocated masks)
        if is_parametric
            threshold_t_matrix_parametric!(
                mask_pos_buffer, mask_neg_buffer,
                t_matrix_perm, critical_t_values, tail
            )
        elseif is_nonparametric_common
            thresh_pos, thresh_neg = critical_t_values
            threshold_t_matrix_nonparametric!(
                mask_pos_buffer, mask_neg_buffer,
                t_matrix_perm, thresh_pos, thresh_neg, tail
            )
        elseif is_nonparametric_individual
            thresh_pos_mat, thresh_neg_mat = critical_t_values
            threshold_t_matrix_nonparametric!(
                mask_pos_buffer, mask_neg_buffer,
                t_matrix_perm, thresh_pos_mat, thresh_neg_mat, tail
            )
        end
        
        # Pre-filter masks (modify in place)
        if min_num_neighbors > 0
            prefilter_mask_by_neighbors!(mask_pos_buffer, spatial_connectivity, min_num_neighbors)
            prefilter_mask_by_neighbors!(mask_neg_buffer, spatial_connectivity, min_num_neighbors)
        end
        
        # Find clusters
        pos_clusters_perm, neg_clusters_perm = find_clusters(
            mask_pos_buffer, mask_neg_buffer, electrodes, time_points,
            spatial_connectivity, cluster_type
        )
        
        # Compute statistics without creating new Cluster objects (use pre-allocated lookup)
        if !isempty(pos_clusters_perm)
            pos_stats_perm = _compute_cluster_statistics_only(
                pos_clusters_perm, t_matrix_perm, electrode_to_idx, cluster_statistic
            )
            max_pos = isempty(pos_stats_perm) ? 0.0 : maximum(pos_stats_perm)
        else
            max_pos = 0.0
        end
        
        if !isempty(neg_clusters_perm)
            neg_stats_perm = _compute_cluster_statistics_only(
                neg_clusters_perm, t_matrix_perm, electrode_to_idx, cluster_statistic
            )
            # For negative clusters, the statistic is negative (sum of negative t-values)
            # The most extreme negative value is the MINIMUM (most negative), not maximum
            max_neg = isempty(neg_stats_perm) ? 0.0 : minimum(neg_stats_perm)
        else
            max_neg = 0.0
        end
        
        push!(permutation_max_positive, max_pos)
        push!(permutation_max_negative, max_neg)
        
        if show_progress
            next!(progress)
        end
    end
    
    return permutation_max_positive, permutation_max_negative
end

# ===================
# PHASE 7: INFERENCE
# ===================

"""
    compute_cluster_pvalues(clusters::Vector{Cluster},
                            cluster_stats::Vector{Float64},
                            permutation_max::Vector{Float64},
                            n_permutations::Int,
                            alpha::Float64 = 0.05)

Compute p-values for clusters by comparing to permutation distribution.

# Arguments
- `clusters::Vector{Cluster}`: Observed clusters
- `cluster_stats::Vector{Float64}`: Observed cluster statistics
- `permutation_max::Vector{Float64}`: Maximum cluster stats from permutations
- `n_permutations::Int`: Number of permutations
- `alpha::Float64`: Significance level (default: 0.05)

# Returns
- `updated_clusters::Vector{Cluster}`: Clusters with p-values and significance flags

# Examples
```julia
clusters = compute_cluster_pvalues(clusters, stats, perm_max, 1000, 0.05)
```
"""
function compute_cluster_pvalues(
    clusters::Vector{Cluster},
    cluster_stats::Vector{Float64},
    permutation_max::Vector{Float64},
    n_permutations::Int,
    alpha::Float64 = 0.05
)
    if isempty(clusters)
        return Cluster[]
    end
    
    updated_clusters = Cluster[]
    
    for (i, cluster) in enumerate(clusters)
        cluster_stat = cluster_stats[i]
        
        # For positive clusters: count permutations with max >= observed
        # For negative clusters: count permutations with max <= observed (more extreme negative)
        # FieldTrip approach: compare observed to permutation distribution
        if cluster.polarity == :positive
            # Positive clusters: larger is more extreme
            count_exceed = sum(permutation_max .>= cluster_stat) + 1
        else
            # Negative clusters: more negative (smaller) is more extreme
            # permutation_max contains minimum (most negative) values from each permutation
            count_exceed = sum(permutation_max .<= cluster_stat) + 1
        end
        p_value = count_exceed / (n_permutations + 1)
        
        is_significant = p_value < alpha
        
        updated_cluster = Cluster(
            cluster.id, cluster.electrodes, cluster.time_indices, cluster.time_range,
            cluster_stat, p_value, is_significant, cluster.polarity
        )
        push!(updated_clusters, updated_cluster)
    end
    
    return updated_clusters
end

# ===================
# PHASE 8: MAIN FUNCTION
# ===================

"""
    cluster_permutation_test(prepared::StatisticalTestData;
                             n_permutations::Int = 1000,
                             threshold::Float64 = 0.05,
                             threshold_method::Symbol = :parametric,
                             cluster_type::Symbol = :spatiotemporal,
                             cluster_statistic::Symbol = :sum,
                             min_num_neighbors::Int = 0,
                             tail::Symbol = :both,
                             random_seed::Union{Int, Nothing} = nothing,
                             show_progress::Bool = true)

Perform cluster-based permutation test on prepared ERP data.

# Arguments
- `prepared::StatisticalTestData`: Prepared data from `prepare_statistical_test_data`
- `n_permutations::Int`: Number of permutations (default: 1000)
- `threshold::Float64`: P-value threshold (default: 0.05)
- `threshold_method::Symbol`: Threshold method - `:parametric` (default), `:nonparametric_individual`, or `:nonparametric_common`
- `cluster_type::Symbol`: Type of clustering - `:spatial`, `:temporal`, or `:spatiotemporal` (default)
- `cluster_statistic::Symbol`: Cluster statistic - `:sum` (default), `:max`, `:size`, or `:wcm`
- `min_num_neighbors::Int`: Minimum number of neighboring significant channels required (FieldTrip's minNumChannels, default: 0). Points with fewer neighbors are removed before clustering.
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`
- `random_seed::Union{Int, Nothing}`: Random seed for reproducibility (default: nothing)
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `ClusterPermutationResult`: Complete results structure

# Examples
```julia
# Prepare data
prepared = prepare_statistical_test_data("erps_good", 1, 2, design=:paired, input_dir="data/")

# Run permutation test
result = cluster_permutation_test(prepared, n_permutations=1000, threshold=0.05)

# Access results
println("Found ", length(result.positive_clusters), " positive clusters")
println("Found ", length(result.negative_clusters), " negative clusters")
```
"""
function cluster_permutation_test(
    prepared::StatisticalTestData;
    n_permutations::Int = 1000,
    threshold::Float64 = 0.05,
    threshold_method::Symbol = :parametric,
    cluster_type::Symbol = :spatiotemporal,
    cluster_statistic::Symbol = :sum,
    min_num_neighbors::Int = 0,
    tail::Symbol = :both,
    random_seed::Union{Int, Nothing} = nothing,
    show_progress::Bool = true
)
    # Validate inputs
    validate_permutation_inputs(
        prepared, n_permutations, threshold, cluster_type, cluster_statistic, tail
    )
    
    # Validate threshold method
    if !(threshold_method in [:parametric, :nonparametric_common, :nonparametric_individual])
        error("threshold_method must be :parametric, :nonparametric_common, or :nonparametric_individual. " *
              "Got :$threshold_method")
    end
    
    # Compute observed t-matrix and df
    @info "Computing t-statistics..."
    t_matrix, df, _ = compute_t_matrix(prepared)
    
    # Handle different thresholding methods
    if threshold_method == :parametric
        # Parametric thresholding: compute critical t-values from t-distribution
        @info "Computing parametric critical t-values..."
        critical_t_values = compute_critical_t_values(df, size(t_matrix), threshold, tail)
        
        # Threshold observed data
        @info "Thresholding observed data (parametric)..."
        mask_positive, mask_negative = threshold_t_matrix_parametric(
            t_matrix, df, critical_t_values, tail
        )
        
        # Store for later use in permutations
        threshold_for_permutations = critical_t_values
        permutation_t_matrices = nothing
        
    elseif threshold_method == :nonparametric_common
        # Non-parametric common: run all permutations first to get threshold
        @info "Collecting permutation t-matrices for non-parametric common thresholding..."
        permutation_t_matrices = collect_permutation_t_matrices(
            prepared, n_permutations, random_seed, show_progress
        )
        
        # Compute common threshold from permutation distribution
        @info "Computing non-parametric common threshold..."
        thresh_pos, thresh_neg = compute_nonparametric_threshold_common(
            permutation_t_matrices, threshold, tail
        )
        critical_t_values = (thresh_pos, thresh_neg)
        
        # Threshold observed data
        @info "Thresholding observed data (non-parametric common)..."
        mask_positive, mask_negative = threshold_t_matrix_nonparametric(
            t_matrix, thresh_pos, thresh_neg, tail
        )
        
        # Store for later use in permutations
        threshold_for_permutations = critical_t_values
        
    elseif threshold_method == :nonparametric_individual
        # Non-parametric individual: run all permutations first to get thresholds
        @info "Collecting permutation t-matrices for non-parametric individual thresholding..."
        permutation_t_matrices = collect_permutation_t_matrices(
            prepared, n_permutations, random_seed, show_progress
        )
        
        # Compute individual thresholds from permutation distribution
        @info "Computing non-parametric individual thresholds..."
        thresh_pos_mat, thresh_neg_mat = compute_nonparametric_threshold_individual(
            permutation_t_matrices, threshold, tail
        )
        critical_t_values = (thresh_pos_mat, thresh_neg_mat)
        
        # Threshold observed data
        @info "Thresholding observed data (non-parametric individual)..."
        mask_positive, mask_negative = threshold_t_matrix_nonparametric(
            t_matrix, thresh_pos_mat, thresh_neg_mat, tail
        )
        
        # Store for later use in permutations
        threshold_for_permutations = critical_t_values
    end
    
    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    time_points = prepared.analysis.time_points
    layout = prepared.data[1].layout
    
    # Build connectivity matrix
    @info "Building connectivity matrix..."
    spatial_connectivity, _, _ = build_connectivity_matrix(
        electrodes, layout, cluster_type
    )
    
    # Pre-filter masks to remove isolated points (FieldTrip's minNumChannels approach)
    if min_num_neighbors > 0
        @info "Pre-filtering masks (min_num_neighbors=$min_num_neighbors)..."
        mask_positive = prefilter_mask_by_neighbors(mask_positive, spatial_connectivity, min_num_neighbors)
        mask_negative = prefilter_mask_by_neighbors(mask_negative, spatial_connectivity, min_num_neighbors)
    end
    
    # Find observed clusters
    @info "Finding observed clusters..."
    positive_clusters, negative_clusters = find_clusters(
        mask_positive, mask_negative, electrodes, time_points,
        spatial_connectivity, cluster_type
    )
    
    # Compute observed cluster statistics
    cluster_stats_positive = Float64[]
    cluster_stats_negative = Float64[]
    
    if !isempty(positive_clusters)
        @info "Computing cluster statistics for $(length(positive_clusters)) positive clusters..."
        positive_clusters, cluster_stats_positive = compute_cluster_statistics(
            positive_clusters, t_matrix, electrodes, cluster_statistic
        )
    end
    
    if !isempty(negative_clusters)
        @info "Computing cluster statistics for $(length(negative_clusters)) negative clusters..."
        negative_clusters, cluster_stats_negative = compute_cluster_statistics(
            negative_clusters, t_matrix, electrodes, cluster_statistic
        )
    end
    
    # Run permutations for cluster-level inference
    # For non-parametric methods, we reuse the stored t-matrices
    if threshold_method == :parametric
        @info "Running $n_permutations permutations for cluster-level inference..."
    else
        @info "Running $n_permutations permutations for cluster-level inference (reusing stored t-matrices)..."
    end
        permutation_max_positive, permutation_max_negative = run_permutations(
        prepared, n_permutations, threshold, threshold_for_permutations, spatial_connectivity,
        cluster_type, cluster_statistic, tail, min_num_neighbors,
        random_seed, show_progress;
        permutation_t_matrices = permutation_t_matrices
    )
    
    # Compute p-values
    @info "Computing p-values..."
    if !isempty(positive_clusters)
        positive_clusters = compute_cluster_pvalues(
            positive_clusters, cluster_stats_positive, permutation_max_positive,
            n_permutations, threshold
        )
        # Sort: significant first, then by cluster statistic (largest first)
        sort!(positive_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end
    
    if !isempty(negative_clusters)
        negative_clusters = compute_cluster_pvalues(
            negative_clusters, cluster_stats_negative, permutation_max_negative,
            n_permutations, threshold
        )
        # Sort: significant first, then by cluster statistic (largest absolute value first)
        sort!(negative_clusters, by = c -> (!c.is_significant, -abs(c.cluster_stat)))
    end
    
    # Create significance masks that only include points from SIGNIFICANT clusters
    # This is the key difference: cluster permutation tests report significance at the cluster level,
    # not the individual point level. Only points that are part of significant clusters are marked as significant.
    significant_mask_positive = falses(size(mask_positive))
    significant_mask_negative = falses(size(mask_negative))
    
    # Mark points from significant positive clusters
    for cluster in positive_clusters
        if cluster.is_significant
            # Find electrode indices
            for electrode in cluster.electrodes
                e_idx = findfirst(==(electrode), electrodes)
                if e_idx !== nothing
                    for t_idx in cluster.time_indices
                        if 1 <= t_idx <= size(significant_mask_positive, 2)
                            significant_mask_positive[e_idx, t_idx] = true
                        end
                    end
                end
            end
        end
    end
    
    # Mark points from significant negative clusters
    for cluster in negative_clusters
        if cluster.is_significant
            # Find electrode indices
            for electrode in cluster.electrodes
                e_idx = findfirst(==(electrode), electrodes)
                if e_idx !== nothing
                    for t_idx in cluster.time_indices
                        if 1 <= t_idx <= size(significant_mask_negative, 2)
                            significant_mask_negative[e_idx, t_idx] = true
                        end
                    end
                end
            end
        end
    end
    
    # Assemble results
    result = ClusterPermutationResult(
        prepared.analysis.design,
        t_matrix,
        df,
        significant_mask_positive,
        significant_mask_negative,
        positive_clusters,
        negative_clusters,
        n_permutations,
        permutation_max_positive,
        permutation_max_negative,
        random_seed,
        electrodes,
        time_points,
        threshold,
        threshold_method,
        cluster_type,
        cluster_statistic,
        critical_t_values
    )
    
    @info "Permutation test complete. Found $(length(positive_clusters)) positive and $(length(negative_clusters)) negative clusters."
    
    return result
end

# ===================
# ANALYTIC T-TEST (Non-permutation)
# ===================

"""
    AnalyticTTestResult

Stores results from an analytic (parametric) t-test without permutation.

# Fields
- `test_type::Symbol`: `:paired` or `:independent`
- `test_info::TestInfo`: Test configuration (type, df, alpha, tail, correction_method)
- `data::Vector{ErpData}`: Grand average ERPs for conditions 1 and 2 (for visualization/storage)
- `stat_matrix::StatMatrix`: Statistical results containing `t` (t-statistics) and `p` (p-values) [electrodes × time]
- `masks::Masks`: Significance masks containing `positive` and `negative` significant points [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
"""
struct TestInfo
    type::Symbol
    df::Float64
    alpha::Float64
    tail::Symbol
    correction_method::Symbol
end

struct StatMatrix
    t::Array{Float64, 2}
    p::Array{Float64, 2}
end

struct Masks
    positive::BitArray{2}
    negative::BitArray{2}
end

struct AnalyticTTestResult
    test_info::TestInfo
    data::Vector{ErpData}
    stat_matrix::StatMatrix
    masks::Masks
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
end

"""
    compute_p_matrix(t_matrix::Array{Float64, 2}, df::Float64, tail::Symbol = :both)

Compute p-values from t-statistics.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df::Float64`: Degrees of freedom (constant across all electrodes/time points)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]
"""
function compute_p_matrix(t_matrix::Array{Float64, 2}, df::Float64, tail::Symbol = :both)
    p_matrix = Array{Float64, 2}(undef, size(t_matrix))
    
    if isnan(df) || isinf(df) || df <= 0
        fill!(p_matrix, NaN)
        return p_matrix
    end
    
    dist = TDist(df)
    
    if tail == :both
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            if isnan(t_val) || isinf(t_val)
                p_matrix[i] = NaN
            else
                p_matrix[i] = 2 * (1 - cdf(dist, abs(t_val)))
            end
        end
    elseif tail == :left
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            if isnan(t_val) || isinf(t_val)
                p_matrix[i] = NaN
            else
                p_matrix[i] = cdf(dist, t_val)
            end
        end
    elseif tail == :right
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            if isnan(t_val) || isinf(t_val)
                p_matrix[i] = NaN
            else
                p_matrix[i] = 1 - cdf(dist, t_val)
            end
        end
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    return p_matrix
end

"""
    apply_multiple_comparison_correction(p_matrix::Array{Float64, 2}, 
                                        alpha::Float64, 
                                        method::Symbol)

Apply multiple comparison correction to p-values.

# Arguments
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]
- `alpha::Float64`: Significance threshold
- `method::Symbol`: Correction method - `:no` or `:bonferroni`

# Returns
- `corrected_mask::BitArray{2}`: Significant points after correction [electrodes × time]
"""
function apply_multiple_comparison_correction(
    p_matrix::Array{Float64, 2},
    alpha::Float64,
    method::Symbol
)
    n_electrodes, n_time = size(p_matrix)
    corrected_mask = BitArray{2}(undef, n_electrodes, n_time)
    
    if method == :no
        # No correction: threshold at alpha
        for i in eachindex(p_matrix)
            corrected_mask[i] = !isnan(p_matrix[i]) && p_matrix[i] <= alpha
        end
        
    elseif method == :bonferroni
        # Bonferroni: alpha / n_comparisons
        n_comparisons = count(!isnan, p_matrix)
        if n_comparisons == 0
            fill!(corrected_mask, false)
        else
            bonferroni_alpha = alpha / n_comparisons
            @info "Bonferroni correction: n_comparisons=$n_comparisons, bonferroni_alpha=$bonferroni_alpha"
            
            # Count how many would be significant with and without correction
            n_sig_uncorrected = 0
            n_sig_bonferroni = 0
            for i in eachindex(p_matrix)
                if !isnan(p_matrix[i])
                    if p_matrix[i] <= alpha
                        n_sig_uncorrected += 1
                    end
                    if p_matrix[i] <= bonferroni_alpha
                        n_sig_bonferroni += 1
                    end
                    corrected_mask[i] = p_matrix[i] <= bonferroni_alpha
                else
                    corrected_mask[i] = false
                end
            end
            @info "Significant points: uncorrected=$n_sig_uncorrected, bonferroni=$n_sig_bonferroni"
        end
        
    else
        error("Unsupported correction method: $method. Use :no or :bonferroni")
    end
    
    return corrected_mask
end

"""
    analytic_ttest(prepared::StatisticalTestData;
                  alpha::Float64 = 0.05,
                  tail::Symbol = :both,
                  correction_method::Symbol = :no)

Perform analytic (parametric) t-test without permutation (FieldTrip's 'analytic' method).

# Arguments
- `prepared::StatisticalTestData`: Prepared data from `prepare_statistical_test_data`
- `alpha::Float64`: Significance threshold (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`
- `correction_method::Symbol`: Multiple comparison correction - `:no` (default) or `:bonferroni`

# Returns
- `AnalyticTTestResult`: Results structure with t-statistics, p-values, and significant masks

# Examples
```julia
# No correction (equivalent to MATLAB 'no')
result = analytic_ttest(prepared, alpha=0.05, correction_method=:no)

# Bonferroni correction
result = analytic_ttest(prepared, alpha=0.05, correction_method=:bonferroni)
```
"""
function analytic_ttest(
    prepared::StatisticalTestData;
    alpha::Float64 = 0.05,
    tail::Symbol = :both,
    correction_method::Symbol = :no
)
    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    
    # Validate correction method
    if !(correction_method in (:no, :bonferroni))
        error("correction_method must be :no or :bonferroni. Got :$correction_method")
    end
    
    # Compute t-statistics, degrees of freedom, and p-values in one pass
    t_matrix, df, p_matrix = compute_t_matrix(prepared, tail=tail)
    
    # Apply multiple comparison correction
    corrected_mask = apply_multiple_comparison_correction(p_matrix, alpha, correction_method)
    
    if correction_method == :bonferroni
        n_comparisons = count(!isnan, p_matrix)
        bonferroni_alpha = alpha / n_comparisons
        n_sig_uncorrected = count(p -> !isnan(p) && p <= alpha, p_matrix)
        n_sig_bonferroni = count(p -> !isnan(p) && p <= bonferroni_alpha, p_matrix)
        @info "Bonferroni correction: n_comparisons=$n_comparisons, bonferroni_alpha=$bonferroni_alpha, uncorrected_sig=$n_sig_uncorrected, bonferroni_sig=$n_sig_bonferroni"
    end
    
    # Create positive and negative masks based on t-values (vectorized)
    if tail == :both
        mask_positive = corrected_mask .& .!isnan.(t_matrix) .& (t_matrix .> 0)
        mask_negative = corrected_mask .& .!isnan.(t_matrix) .& (t_matrix .< 0)
    elseif tail == :right
        mask_positive = corrected_mask .& .!isnan.(t_matrix)
        mask_negative = falses(size(t_matrix))
    elseif tail == :left
        mask_positive = falses(size(t_matrix))
        mask_negative = corrected_mask .& .!isnan.(t_matrix)
    else
        error("tail must be :both, :left, or :right, got :$tail")
    end
    
    # Extract commonly used fields from grand_average ErpData
    electrodes = channel_labels(prepared.data[1])
    
    test_info = TestInfo(
        prepared.analysis.design,
        df,
        alpha,
        tail,
        correction_method
    )
    
    stat_matrix = StatMatrix(t_matrix, p_matrix)
    masks = Masks(mask_positive, mask_negative)
    
    result = AnalyticTTestResult(
        test_info,
        prepared.data,
        stat_matrix,
        masks,
        electrodes,
        prepared.analysis.time_points
    )
    
    n_sig_pos = count(mask_positive)
    n_sig_neg = count(mask_negative)
    @info "Analytic t-test complete. Found $n_sig_pos positive and $n_sig_neg negative significant points."
    
    return result
end

function Base.show(io::IO, result::AnalyticTTestResult)
    n_electrodes = length(result.electrodes)
    n_time_points = length(result.time_points)
    time_range = isempty(result.time_points) ? "N/A" : "$(first(result.time_points)) to $(last(result.time_points)) s"
    
    n_sig_pos = count(result.masks.positive)
    n_sig_neg = count(result.masks.negative)
    
    println(io, "AnalyticTTestResult")
    println(io, "├─ Test info")
    println(io, "│  ├─ Type: $(result.test_info.type)")
    println(io, "│  ├─ DF: $(result.test_info.df)")
    println(io, "│  ├─ Alpha: $(result.test_info.alpha)")
    println(io, "│  ├─ Tail: $(result.test_info.tail)")
    println(io, "│  └─ Correction method: $(result.test_info.correction_method)")
    println(io, "├─ Data dimensions: $n_electrodes electrodes × $n_time_points time points ($time_range)")
    println(io, "└─ Significant points: $n_sig_pos positive, $n_sig_neg negative")
end
