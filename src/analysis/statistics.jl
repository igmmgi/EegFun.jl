"""
Statistical test functions for EEG/ERP data analysis.

This module provides t-test functions for use in cluster-based permutation tests
and other statistical analyses. Functions are provided to compute t-values and/or
p-values for paired and independent sample t-tests.
"""

# Use NamedTuple for result type as seems lighter than full struct
const TTestResult = NamedTuple{(:df, :t, :p), Tuple{Float64, Float64, Float64}}

# Format t-value and p-value for display
function Base.show(io::IO, result::TTestResult)
    t_str = isnan(result.t) ? "NaN" : (isinf(result.t) ? (result.t > 0 ? "Inf" : "-Inf") : Printf.@sprintf("%.4f", result.t))
    p_str = isnan(result.p) ? "NaN" : Printf.@sprintf("%.4f", result.p)
    print(io, "t($(round(Int, result.df))) = $t_str, p = $p_str")
end

# ==============
# PAIRED T-TEST
# ==============
"""
    paired_ttest(data_A::AbstractVector, data_B::AbstractVector; tail::Symbol = :both)

Compute paired t-test with degrees of freedom, t-statistic, and p-value.

# Arguments
- `data_A::AbstractVector`: First condition data (must have same length as data_B)
- `data_B::AbstractVector`: Second condition data (must have same length as data_A)
- `tail::Symbol`: Type of test - `:both` (two-tailed, default), `:left` (one-tailed, A < B), or `:right` (one-tailed, A > B)

# Returns
- `TTestResult`: Struct containing `df` (degrees of freedom), `t` (t-statistic), and `p` (p-value).
  Returns `NaN` for p-value if t-value is `NaN` or `Inf`.

# Examples
```julia
using eegfun

# Two-tailed test
condition_A = [1.0, 2.0, 3.0, 4.0, 5.0]
condition_B = [1.5, 2.5, 3.5, 4.5, 5.5]
result = paired_ttest(condition_A, condition_B)
# Prints: t(4) = -0.7071, p = 0.5129

# Access individual values
result.df  # degrees of freedom
result.t   # t-statistic
result.p   # p-value

# One-tailed p-value (testing if A < B)
result = paired_ttest(condition_A, condition_B, tail = :left)
```
"""
function paired_ttest(data_A::AbstractVector, data_B::AbstractVector; tail::Symbol = :both)

    # validate equal lengths
    length(data_A) == length(data_B) || error("Paired t-test requires equal sample sizes")
    
    n = length(data_A)
    n < 2 && return (df = NaN, t = NaN, p = NaN)  # Need at least 2 observations
    
    df = n - 1  # degrees of freedom for paired t-test
    
    # Compute mean and std of differences 
    diff = data_A .- data_B
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
    independent_ttest(data_A::AbstractVector, data_B::AbstractVector; tail::Symbol = :both)

Compute independent sample t-test with degrees of freedom, t-statistic, and p-value.
Assumes equal variances (standard independent t-test).

# Arguments
- `data_A::AbstractVector`: First group data (must have at least 2 observations)
- `data_B::AbstractVector`: Second group data (must have at least 2 observations)
- `tail::Symbol`: Type of test - `:both` (two-tailed, default), `:left` (one-tailed, A < B), or `:right` (one-tailed, A > B)

# Returns
- `TTestResult`: Struct containing `df` (degrees of freedom), `t` (t-statistic), and `p` (p-value).
  Returns `NaN` for p-value if t-value is `NaN` or `Inf`.

# Examples
```julia
using eegfun

# Two-tailed test
group_A = [1.0, 2.0, 3.0, 4.0, 5.0]
group_B = [6.0, 7.0, 8.0, 9.0, 10.0]
result = independent_ttest(group_A, group_B)
# Prints: t(8) = -5.0000, p = 0.0011

# Access individual values
result.df  # degrees of freedom
result.t   # t-statistic
result.p   # p-value

# One-tailed p-value (testing if A < B)
result = independent_ttest(group_A, group_B, tail = :left)
```
"""
function independent_ttest(data_A::AbstractVector, data_B::AbstractVector; tail::Symbol = :both)

    # Validate input lengths
    n_A, n_B = length(data_A), length(data_B)
    if n_A < 2 || n_B < 2
        error("Independent t-test requires at least 2 observations per group")
    end
    
    df = n_A + n_B - 2  # degrees of freedom for independent t-test
    
    # Compute means and variances (Statistics.jl is already optimized with SIMD)
    mean_A, mean_B = mean(data_A), mean(data_B)
    
    # Pooled variance (assuming equal variances)
    var_A = var(data_A, corrected = true)
    var_B = var(data_B, corrected = true)
    pooled_var = ((n_A - 1) * var_A + (n_B - 1) * var_B) / df
    
    # Compute t-value
    if pooled_var == 0.0
        if mean_A == mean_B
            t = NaN  # Both groups identical
        else
            t = Inf * sign(mean_A - mean_B)
        end
    else
        t = (mean_A - mean_B) / sqrt(pooled_var * (1/n_A + 1/n_B))
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

# ===================
# DATA PREPARATION
# ===================

"""
    PermutationTestData

Stores prepared data for cluster-based permutation tests.

# Fields
- `design::Symbol`: Design type - `:paired` or `:independent`
- `data_A::Array{Float64, 3}`: Data for condition/group A [participants × electrodes × time]
- `data_B::Array{Float64, 3}`: Data for condition/group B [participants × electrodes × time]
- `participants_A::Vector{Int}`: Participant IDs for condition/group A
- `participants_B::Vector{Int}`: Participant IDs for condition/group B
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `condition_A::String`: Condition/group A name
- `condition_B::String`: Condition/group B name
- `layout::Layout`: Layout object for spatial clustering
"""
struct PermutationTestData
    design::Symbol
    data_A::Array{Float64, 3}  # [participants × electrodes × time]
    data_B::Array{Float64, 3}  # [participants × electrodes × time]
    participants_A::Vector{Int}
    participants_B::Vector{Int}
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
    condition_A::String
    condition_B::String
    layout::Layout
end

function Base.show(io::IO, data::PermutationTestData)
    n_participants_A = length(data.participants_A)
    n_participants_B = length(data.participants_B)
    n_electrodes = length(data.electrodes)
    n_time_points = length(data.time_points)
    time_range = isempty(data.time_points) ? "N/A" : "$(first(data.time_points)) to $(last(data.time_points)) s"
    
    println(io, "PermutationTestData")
    println(io, "├─ Design: $(data.design)")
    println(io, "├─ Condition A ($(data.condition_A)): $n_participants_A participants")
    println(io, "├─ Condition B ($(data.condition_B)): $n_participants_B participants")
    println(io, "├─ Electrodes: $n_electrodes")
    println(io, "├─ Time points: $n_time_points ($time_range)")
    println(io, "├─ Data dimensions: $(size(data.data_A)) (participants × electrodes × time)")
    print(io, "└─ Layout: $(n_electrodes) channels")
end

"""
    prepare_permutation_data(erps_A::Vector{ErpData}, erps_B::Vector{ErpData}; 
                             design::Symbol = :paired,
                             channel_selection::Function = channels(),
                             sample_selection::Function = samples(),
                             baseline_window::Function = samples())

    prepare_permutation_data(file_pattern::String, condition_A, condition_B;
                             design::Symbol = :paired,
                             input_dir::String = pwd(),
                             participant_selection::Function = participants(),
                             channel_selection::Function = channels(),
                             sample_selection::Function = samples(),
                             baseline_window::Function = samples())

Prepare ErpData for cluster-based permutation tests.

Organizes ErpData into participant × electrode × time arrays for statistical analysis.
Validates the design and ensures data consistency across conditions.

# Arguments (Direct data version)
- `erps_A::Vector{ErpData}`: ERPs for condition/group A (one per participant)
- `erps_B::Vector{ErpData}`: ERPs for condition/group B (one per participant)
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `channel_selection::Function`: Function to filter channels (default: all channels)
- `sample_selection::Function`: Function to select time points (default: all samples). Use `samples((start, end))` for time windows.
- `baseline_window::Function`: Baseline window sample selection predicate (default: samples() - all samples, baseline skipped). Use `samples((start, end))` for baseline window.

# Arguments (File-based version)
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned")
- `condition_A`: Condition identifier for group A (Int or String)
- `condition_B`: Condition identifier for group B (Int or String)
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `input_dir::String`: Directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Function to filter participants (default: all)
- `channel_selection::Function`: Function to filter channels (default: all)
- `sample_selection::Function`: Function to select time points (default: all samples). Use `samples((start, end))` for time windows.
- `baseline_window::Function`: Baseline window sample selection predicate (default: samples() - all samples, baseline skipped). Use `samples((start, end))` for baseline window.

# Returns
- `PermutationTestData`: Prepared data structure ready for permutation testing

# Examples
```julia
using eegfun

# Direct data version - work with already-loaded ErpData
erps_condition_A = load_erps("condition_A")  # Vector of ErpData
erps_condition_B = load_erps("condition_B")  # Vector of ErpData
prepared = prepare_permutation_data(erps_condition_A, erps_condition_B, design = :paired)

# File-based version - load from files
prepared = prepare_permutation_data(
    "erps_cleaned",
    condition_A = 1,
    condition_B = 2,
    design = :paired,
    input_dir = "data/erps"
)

# With time window selection and baseline correction
prepared = prepare_permutation_data(
    "erps_cleaned",
    1,
    2;
    design = :paired,
    sample_selection = samples((0.0, 1.0)),  # Analysis window: 0 to 1 second
    baseline_window = samples((-0.2, 0.0))   # Baseline: -0.2 to 0 seconds
)
```
"""
# Direct data version
function prepare_permutation_data(erps_A::Vector{ErpData}, erps_B::Vector{ErpData}; 
                                   design::Symbol = :paired,
                                   channel_selection::Function = channels(),
                                   sample_selection::Function = samples(),
                                   baseline_window::Function = samples())
    # Validate design
    design in (:paired, :independent) || error("design must be :paired or :independent, got :$design")
    
    # Extract participant IDs from filenames (using utility from batch.jl)
    participants_A = [_extract_participant_id(basename(erp.file)) for erp in erps_A]
    participants_B = [_extract_participant_id(basename(erp.file)) for erp in erps_B]
    
    # Validate design
    if design == :paired
        # Paired design: same participants in both conditions
        if sort(participants_A) != sort(participants_B)
            error("Paired design requires same participants in both conditions. " *
                  "Condition A participants: $(sort(participants_A)), " *
                  "Condition B participants: $(sort(participants_B))")
        end
        # Sort to match participants across conditions
        sort_idx_A = sortperm(participants_A)
        sort_idx_B = sortperm(participants_B)
        erps_A = erps_A[sort_idx_A]
        erps_B = erps_B[sort_idx_B]
        participants_A = participants_A[sort_idx_A]
        participants_B = participants_B[sort_idx_B]
    else
        # Independent design: different participants (or allow overlap)
        # Just ensure we have at least 2 participants per group
        if length(participants_A) < 2 || length(participants_B) < 2
            error("Independent design requires at least 2 participants per group. " *
                  "Got $(length(participants_A)) and $(length(participants_B))")
        end
    end
    
    # Get condition identifiers
    condition_A = length(erps_A) > 0 ? erps_A[1].condition_name : "A"
    condition_B = length(erps_B) > 0 ? erps_B[1].condition_name : "B"
    
    # Validate all ERPs have same structure
    first_erp_A = erps_A[1]
    first_erp_B = erps_B[1]
    
    # Check sample rates match
    if first_erp_A.sample_rate != first_erp_B.sample_rate
        error("Sample rates must match across conditions. " *
              "Condition A: $(first_erp_A.sample_rate) Hz, " *
              "Condition B: $(first_erp_B.sample_rate) Hz")
    end
    
    # Apply baseline correction if baseline_window is specified (not default all samples)
    # Make copies of ERPs to avoid modifying originals
    erps_A = [copy(erp) for erp in erps_A]
    erps_B = [copy(erp) for erp in erps_B]
    
    # Use first ERP from copied list for baseline calculation
    first_erp_A_baseline = erps_A[1]
    baseline_mask = baseline_window(first_erp_A_baseline.data)
    baseline_indices = findall(baseline_mask)
    
    # Skip baseline if window matches all samples (default) or no samples
    if length(baseline_indices) < nrow(first_erp_A_baseline.data) && !isempty(baseline_indices)
        # Get channel columns (exclude metadata)
        metadata_cols = meta_labels(first_erp_A_baseline)
        all_channels = setdiff(propertynames(first_erp_A_baseline.data), metadata_cols)
        eeg_channels = all_channels[channels()(all_channels)]
        
        if !isempty(eeg_channels)
            # Use existing baseline infrastructure from baseline.jl
            interval = IntervalIndex(start = first(baseline_indices), stop = last(baseline_indices))
            @info "Applying baseline correction to $(length(eeg_channels)) channels over interval: $(first(baseline_indices)) to $(last(baseline_indices))"
            
            # Apply baseline to all ERPs
            for erp in erps_A
                _apply_baseline!(erp.data, eeg_channels, interval)
            end
            for erp in erps_B
                _apply_baseline!(erp.data, eeg_channels, interval)
            end
        else
            @warn "No EEG channels found for baseline correction"
        end
    elseif isempty(baseline_indices)
        @warn "Baseline window matched no samples. Skipping baseline correction."
    end
    
    # Get all electrodes from first ERP (use copied version after baseline)
    first_erp_A_final = erps_A[1]
    first_erp_B_final = erps_B[1]
    all_electrodes_A = channel_labels(first_erp_A_final)
    all_electrodes_B = channel_labels(first_erp_B_final)
    
    # Check electrodes match
    if all_electrodes_A != all_electrodes_B
        error("Electrodes must match across conditions. " *
              "Condition A: $all_electrodes_A, Condition B: $all_electrodes_B")
    end
    
    # Apply channel selection
    channel_mask = channel_selection(all_electrodes_A)
    electrodes = all_electrodes_A[channel_mask]
    
    if isempty(electrodes)
        error("No electrodes selected after applying channel_selection")
    end
    
    # Apply sample selection to get time points
    # Sample selection predicate operates on DataFrame and returns boolean vector
    sample_mask_A = sample_selection(first_erp_A_final.data)
    sample_mask_B = sample_selection(first_erp_B_final.data)
    
    # Get selected time points
    time_points_A = first_erp_A_final.data[sample_mask_A, :time]
    time_points_B = first_erp_B_final.data[sample_mask_B, :time]
    
    # Check time points match (or align them)
    if time_points_A != time_points_B
        @warn "Time points differ between conditions. Using intersection."
        # Use intersection of time points
        common_times = intersect(time_points_A, time_points_B)
        if isempty(common_times)
            error("No common time points between conditions after sample selection")
        end
        time_points = sort(common_times)
    else
        time_points = time_points_A
    end
    
    if isempty(time_points)
        error("No time points selected after applying sample_selection")
    end
    
    # Get layout (use from first ERP - layout is not modified by baseline)
    layout = first_erp_A_final.layout
    
    # Extract data arrays: [participants × electrodes × time]
    n_participants_A = length(erps_A)
    n_participants_B = length(erps_B)
    n_electrodes = length(electrodes)
    n_time = length(time_points)
    
    data_A = Array{Float64, 3}(undef, n_participants_A, n_electrodes, n_time)
    data_B = Array{Float64, 3}(undef, n_participants_B, n_electrodes, n_time)
    
    # Fill data arrays
    for (p_idx, erp) in enumerate(erps_A)
        # Apply sample selection to this ERP
        sample_mask = sample_selection(erp.data)
        erp_times = erp.data[sample_mask, :time]
        
        # Get time indices for selected time points
        time_indices = [findfirst(==(t), erp_times) for t in time_points]
        
        # Extract electrode data from selected samples
        for (e_idx, electrode) in enumerate(electrodes)
            selected_data = erp.data[sample_mask, electrode]
            for (t_idx, time_idx) in enumerate(time_indices)
                if time_idx !== nothing
                    data_A[p_idx, e_idx, t_idx] = selected_data[time_idx]
                else
                    data_A[p_idx, e_idx, t_idx] = NaN
                end
            end
        end
    end
    
    for (p_idx, erp) in enumerate(erps_B)
        # Apply sample selection to this ERP
        sample_mask = sample_selection(erp.data)
        erp_times = erp.data[sample_mask, :time]
        
        # Get time indices for selected time points
        time_indices = [findfirst(==(t), erp_times) for t in time_points]
        
        # Extract electrode data from selected samples
        for (e_idx, electrode) in enumerate(electrodes)
            selected_data = erp.data[sample_mask, electrode]
            for (t_idx, time_idx) in enumerate(time_indices)
                if time_idx !== nothing
                    data_B[p_idx, e_idx, t_idx] = selected_data[time_idx]
                else
                    data_B[p_idx, e_idx, t_idx] = NaN
                end
            end
        end
    end
    
    return PermutationTestData(
        design,
        data_A,
        data_B,
        participants_A,
        participants_B,
        electrodes,
        time_points,
        condition_A,
        condition_B,
        layout
    )
end

# File-based version (convenience wrapper)
function prepare_permutation_data(
    file_pattern::String,
    condition_A,
    condition_B;
    design::Symbol = :paired,
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_window::Function = samples()
)
    # Find files
    files = _find_batch_files(file_pattern, input_dir; participants = participant_selection)
    if isempty(files)
        error("No JLD2 files found matching pattern '$file_pattern' in $input_dir")
    end
    
    # Load and group ERPs by condition (similar to grand_average.jl)
    all_erps_by_condition = Dict{Int,Vector{ErpData}}()
    
    for file in files
        input_path = joinpath(input_dir, file)
        
        # Load ERP data (using load_data which finds by type)
        erps_data = load_data(input_path)
        
        if isnothing(erps_data)
            @warn "No data variables found in $file. Skipping."
            continue
        end
        
        # Validate that data is Vector{ErpData}
        if !(erps_data isa Vector{<:ErpData})
            @warn "Invalid data type in $file: expected Vector{ErpData}, got $(typeof(erps_data)). Skipping."
            continue
        end
        
        # Group ERPs by condition
        for erp in erps_data
            cond_num = erp.condition
            if !haskey(all_erps_by_condition, cond_num)
                all_erps_by_condition[cond_num] = ErpData[]
            end
            push!(all_erps_by_condition[cond_num], erp)
        end
    end
    
    # Convert condition identifiers to condition numbers if needed
    function get_condition_number(cond_id, erps_by_cond)
        if cond_id isa Int
            return cond_id
        elseif cond_id isa String
            # Find condition number by name
            for (cond_num, erps) in erps_by_cond
                if !isempty(erps) && erps[1].condition_name == cond_id
                    return cond_num
                end
            end
            error("Condition name '$cond_id' not found in data")
        else
            error("condition_A and condition_B must be Int or String, got $(typeof(cond_id))")
        end
    end
    
    cond_num_A = get_condition_number(condition_A, all_erps_by_condition)
    cond_num_B = get_condition_number(condition_B, all_erps_by_condition)
    
    # Get ERPs for each condition
    if !haskey(all_erps_by_condition, cond_num_A)
        error("Condition $cond_num_A not found in data")
    end
    if !haskey(all_erps_by_condition, cond_num_B)
        error("Condition $cond_num_B not found in data")
    end
    
    erps_A = all_erps_by_condition[cond_num_A]
    erps_B = all_erps_by_condition[cond_num_B]
    
    # Call the main preparation function
    return prepare_permutation_data(
        erps_A, 
        erps_B; 
        design = design, 
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        baseline_window = baseline_window
    )
end

# ===================
# CLUSTER PERMUTATION TESTS
# ===================

using SparseArrays

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
- `test_type::Symbol`: `:paired` or `:independent`
- `design::Symbol`: Design type (same as test_type)
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df_matrix::Array{Float64, 2}`: Degrees of freedom [electrodes × time]
- `significant_mask_positive::BitArray{2}`: Positive significant points [electrodes × time]
- `significant_mask_negative::BitArray{2}`: Negative significant points [electrodes × time]
- `positive_clusters::Vector{Cluster}`: Positive clusters
- `negative_clusters::Vector{Cluster}`: Negative clusters
- `cluster_stats_positive::Vector{Float64}`: Cluster statistics for positive clusters
- `cluster_stats_negative::Vector{Float64}`: Cluster statistics for negative clusters
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
- `min_cluster_size::Int`: Minimum cluster size used
- `critical_t_values::Array{Float64, 2}`: Critical t-values used [electrodes × time]
"""
struct ClusterPermutationResult
    test_type::Symbol
    design::Symbol
    t_matrix::Array{Float64, 2}
    df_matrix::Array{Float64, 2}
    significant_mask_positive::BitArray{2}
    significant_mask_negative::BitArray{2}
    positive_clusters::Vector{Cluster}
    negative_clusters::Vector{Cluster}
    cluster_stats_positive::Vector{Float64}
    cluster_stats_negative::Vector{Float64}
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
    min_cluster_size::Int
    critical_t_values::Array{Float64, 2}
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
    n_sig_pos_points = count(result.significant_mask_positive)
    n_sig_neg_points = count(result.significant_mask_negative)
    
    println(io, "ClusterPermutationResult")
    println(io, "├─ Test type: $(result.test_type) ($(result.design))")
    println(io, "├─ Permutations: $(result.n_permutations)")
    println(io, "├─ Threshold: $(result.threshold) ($(result.threshold_method))")
    println(io, "├─ Cluster type: $(result.cluster_type)")
    println(io, "├─ Cluster statistic: $(result.cluster_statistic)")
    println(io, "├─ Min cluster size: $(result.min_cluster_size)")
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
    validate_permutation_inputs(prepared::PermutationTestData, 
                                n_permutations::Int,
                                threshold::Float64,
                                cluster_type::Symbol,
                                cluster_statistic::Symbol,
                                tail::Symbol,
                                min_cluster_size::Int)

Validate inputs for cluster permutation test.

# Arguments
- `prepared::PermutationTestData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `threshold::Float64`: P-value threshold
- `cluster_type::Symbol`: Type of clustering
- `cluster_statistic::Symbol`: Cluster statistic type
- `tail::Symbol`: Test tail
- `min_cluster_size::Int`: Minimum cluster size

# Throws
- `ArgumentError`: If any validation fails
"""
function validate_permutation_inputs(
    prepared::PermutationTestData,
    n_permutations::Int,
    threshold::Float64,
    cluster_type::Symbol,
    cluster_statistic::Symbol,
    tail::Symbol,
    min_cluster_size::Int
)
    # Validate design
    if prepared.design ∉ (:paired, :independent)
        error("Design must be :paired or :independent, got :$(prepared.design)")
    end
    
    # Validate minimum participants
    if prepared.design == :paired
        if length(prepared.participants_A) < 2
            error("Paired design requires at least 2 participants, got $(length(prepared.participants_A))")
        end
    else
        if length(prepared.participants_A) < 2 || length(prepared.participants_B) < 2
            error("Independent design requires at least 2 participants per group, " *
                  "got $(length(prepared.participants_A)) and $(length(prepared.participants_B))")
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
    
    # Validate min_cluster_size
    if min_cluster_size < 0
        error("min_cluster_size must be >= 0, got $min_cluster_size")
    end
    
    # Validate Layout.neighbours for spatial clustering
    if cluster_type ∈ (:spatial, :spatiotemporal)
        if isnothing(prepared.layout.neighbours)
            @warn "Layout.neighbours is not set. Computing with default distance criterion (40.0)."
            # Compute neighbours if missing (using default criterion)
            get_layout_neighbours_xy!(prepared.layout, 40.0)
        end
    end
    
    return nothing
end

# ===================
# PHASE 1: CORE STATISTICS
# ===================

"""
    compute_t_matrix(prepared::PermutationTestData)

Compute t-statistics for all electrode × time points.

# Arguments
- `prepared::PermutationTestData`: Prepared data for permutation test

# Returns
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df_matrix::Array{Float64, 2}`: Degrees of freedom [electrodes × time]

# Examples
```julia
t_matrix, df_matrix = compute_t_matrix(prepared)
```
"""
function compute_t_matrix(prepared::PermutationTestData)
    n_electrodes = length(prepared.electrodes)
    n_time = length(prepared.time_points)
    t_matrix = Array{Float64, 2}(undef, n_electrodes, n_time)
    df_matrix = Array{Float64, 2}(undef, n_electrodes, n_time)
    
    for e_idx in 1:n_electrodes
        for t_idx in 1:n_time
            data_A = prepared.data_A[:, e_idx, t_idx]
            data_B = prepared.data_B[:, e_idx, t_idx]
            
            if prepared.design == :paired
                result = paired_ttest(data_A, data_B)
            else
                result = independent_ttest(data_A, data_B)
            end
            t_matrix[e_idx, t_idx] = result.t
            df_matrix[e_idx, t_idx] = result.df
        end
    end
    return t_matrix, df_matrix
end

"""
    compute_critical_t_values(df_matrix::Array{Float64, 2}, 
                              alpha::Float64 = 0.05, 
                              tail::Symbol = :both)

Compute critical t-values for parametric thresholding.

# Arguments
- `df_matrix::Array{Float64, 2}`: Degrees of freedom [electrodes × time]
- `alpha::Float64`: Significance level (default: 0.05)
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time]

# Examples
```julia
critical_t = compute_critical_t_values(df_matrix, 0.05, :both)
```
"""
function compute_critical_t_values(
    df_matrix::Array{Float64, 2},
    alpha::Float64 = 0.05,
    tail::Symbol = :both
)
    critical_t_values = Array{Float64, 2}(undef, size(df_matrix))
    
    if tail == :both
        # Two-tailed: split alpha
        alpha_per_tail = alpha / 2.0
        for i in eachindex(df_matrix)
            df = df_matrix[i]
            if isnan(df) || isinf(df) || df <= 0
                critical_t_values[i] = NaN
            else
                dist = TDist(df)
                critical_t_values[i] = quantile(dist, 1.0 - alpha_per_tail)
            end
        end
    elseif tail == :right
        # One-tailed right: H1: A > B
        for i in eachindex(df_matrix)
            df = df_matrix[i]
            if isnan(df) || isinf(df) || df <= 0
                critical_t_values[i] = NaN
            else
                dist = TDist(df)
                critical_t_values[i] = quantile(dist, 1.0 - alpha)
            end
        end
    elseif tail == :left
        # One-tailed left: H1: A < B
        for i in eachindex(df_matrix)
            df = df_matrix[i]
            if isnan(df) || isinf(df) || df <= 0
                critical_t_values[i] = NaN
            else
                dist = TDist(df)
                critical_t_values[i] = quantile(dist, alpha)
            end
        end
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
                                   df_matrix::Array{Float64, 2},
                                   critical_t_values::Array{Float64, 2},
                                   tail::Symbol = :both)

Threshold t-matrix using parametric critical values.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df_matrix::Array{Float64, 2}`: Degrees of freedom [electrodes × time]
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time]
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `mask_positive::BitArray{2}`: Mask for positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Mask for negative significant points [electrodes × time]

# Examples
```julia
mask_pos, mask_neg = threshold_t_matrix_parametric(t_matrix, df_matrix, critical_t, :both)
```
"""
function threshold_t_matrix_parametric(
    t_matrix::Array{Float64, 2},
    df_matrix::Array{Float64, 2},
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
            # Self-connection (each point is connected to itself)
            push!(I, e_idx)
            push!(J, e_idx)
            
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
            push!(I, e_idx)
            push!(J, e_idx)
            
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
                # This includes the point itself (if self-connected in connectivity matrix)
                # FieldTrip's logic: count total neighbors including self, then check if < minnbchan
                neighbor_count = 0
                for n_e_idx in 1:n_electrodes
                    if spatial_conn_sym[e_idx, n_e_idx] && filtered_mask[n_e_idx, t_idx]
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

# ===================
# PHASE 4: CLUSTER FINDING
# ===================

"""
    find_clusters_connected_components(mask::BitArray{2},
                                       electrodes::Vector{Symbol},
                                       time_points::Vector{Float64},
                                       spatial_connectivity::SparseMatrixCSC{Bool},
                                       cluster_type::Symbol,
                                       min_cluster_size::Int = 0)

Find connected clusters in thresholded data using BFS.

# Arguments
- `mask::BitArray{2}`: Significant points mask [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`
- `min_cluster_size::Int`: Minimum number of points in cluster (default: 0)

# Returns
- `clusters::Vector{Cluster}`: Found clusters (without statistics)

# Examples
```julia
clusters = find_clusters_connected_components(mask, electrodes, time_points, conn, :spatiotemporal, 0)
```
"""
function find_clusters_connected_components(
    mask::BitArray{2},
    electrodes::Vector{Symbol},
    time_points::Vector{Float64},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    min_cluster_size::Int = 0
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
            
            # Check minimum cluster size
            cluster_size = length(cluster_electrodes) * length(cluster_time_indices)
            if cluster_size >= min_cluster_size
                time_indices_vec = sort(collect(cluster_time_indices))
                time_range = (time_points[time_indices_vec[1]], time_points[time_indices_vec[end]])
                
                cluster = Cluster(
                    current_cluster_id,
                    collect(cluster_electrodes),
                    time_indices_vec,
                    time_range,
                    0.0,  # cluster_stat - will be computed later
                    1.0,  # p_value - will be computed later
                    false,  # is_significant - will be computed later
                    :positive  # polarity - will be set by caller (temporary)
                )
                push!(clusters, cluster)
            end
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
                  cluster_type::Symbol,
                  min_cluster_size::Int = 0)

Find positive and negative clusters.

# Arguments
- `mask_positive::BitArray{2}`: Positive significant points [electrodes × time]
- `mask_negative::BitArray{2}`: Negative significant points [electrodes × time]
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Spatial connectivity matrix
- `cluster_type::Symbol`: `:spatial`, `:temporal`, or `:spatiotemporal`
- `min_cluster_size::Int`: Minimum cluster size (default: 0)

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
    cluster_type::Symbol,
    min_cluster_size::Int = 0
)
    # Find positive clusters
    positive_clusters_raw = find_clusters_connected_components(
        mask_positive, electrodes, time_points, spatial_connectivity, cluster_type, min_cluster_size
    )
    # Set polarity to positive
    positive_clusters = [
        Cluster(
            c.id, c.electrodes, c.time_indices, c.time_range,
            c.cluster_stat, c.p_value, c.is_significant, :positive
        ) for c in positive_clusters_raw
    ]
    
    # Find negative clusters
    negative_clusters_raw = find_clusters_connected_components(
        mask_negative, electrodes, time_points, spatial_connectivity, cluster_type, min_cluster_size
    )
    # Set polarity to negative
    negative_clusters = [
        Cluster(
            c.id, c.electrodes, c.time_indices, c.time_range,
            c.cluster_stat, c.p_value, c.is_significant, :negative
        ) for c in negative_clusters_raw
    ]
    
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
    shuffle_labels(prepared::PermutationTestData, rng::AbstractRNG = Random.GLOBAL_RNG)

Shuffle condition/group labels for permutation test.

# Arguments
- `prepared::PermutationTestData`: Prepared data
- `rng::AbstractRNG`: Random number generator (default: global RNG)

# Returns
- `shuffled_data_A::Array{Float64, 3}`: Shuffled data for condition A
- `shuffled_data_B::Array{Float64, 3}`: Shuffled data for condition B

# Notes
For paired design: randomly swap A/B labels within each participant.
For independent design: randomly shuffle participants between groups.
"""
function shuffle_labels(prepared::PermutationTestData, rng::AbstractRNG = Random.GLOBAL_RNG)
    if prepared.design == :paired
        # Paired: randomly swap A/B within each participant
        shuffled_A = copy(prepared.data_A)
        shuffled_B = copy(prepared.data_B)
        
        n_participants = size(prepared.data_A, 1)
        for p_idx in 1:n_participants
            if rand(rng, Bool)
                # Swap conditions for this participant
                shuffled_A[p_idx, :, :], shuffled_B[p_idx, :, :] = 
                    shuffled_B[p_idx, :, :], shuffled_A[p_idx, :, :]
            end
        end
        return shuffled_A, shuffled_B
        
    else
        # Independent: randomly shuffle participants between groups
        n_A = size(prepared.data_A, 1)
        n_B = size(prepared.data_B, 1)
        n_total = n_A + n_B
        
        # Combine all participants
        all_data = cat(prepared.data_A, prepared.data_B, dims=1)
        all_participants = vcat(prepared.participants_A, prepared.participants_B)
        
        # Shuffle indices
        shuffled_indices = collect(1:n_total)
        shuffle!(rng, shuffled_indices)
        
        # Split back into groups
        shuffled_A = all_data[shuffled_indices[1:n_A], :, :]
        shuffled_B = all_data[shuffled_indices[(n_A+1):end], :, :]
        
        return shuffled_A, shuffled_B
    end
end

"""
    run_permutations(prepared::PermutationTestData,
                     n_permutations::Int,
                     threshold::Float64,
                     critical_t_values::Array{Float64, 2},
                     spatial_connectivity::SparseMatrixCSC{Bool},
                     cluster_type::Symbol,
                     cluster_statistic::Symbol,
                     tail::Symbol,
                     min_cluster_size::Int,
                     random_seed::Union{Int, Nothing} = nothing,
                     show_progress::Bool = true)

Run permutation loop to generate distribution of maximum cluster statistics.

# Arguments
- `prepared::PermutationTestData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `threshold::Float64`: P-value threshold
- `critical_t_values::Array{Float64, 2}`: Critical t-values
- `spatial_connectivity::SparseMatrixCSC{Bool}`: Connectivity matrix
- `cluster_type::Symbol`: Type of clustering
- `cluster_statistic::Symbol`: Cluster statistic type
- `tail::Symbol`: Test tail
- `min_cluster_size::Int`: Minimum cluster size
- `random_seed::Union{Int, Nothing}`: Random seed (default: nothing)
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `permutation_max_positive::Vector{Float64}`: Max cluster stats from permutations (positive)
- `permutation_max_negative::Vector{Float64}`: Max cluster stats from permutations (negative)

# Examples
```julia
perm_pos, perm_neg = run_permutations(prepared, 1000, 0.05, critical_t, conn, 
                                     :spatiotemporal, :sum, :both, 0)
```
"""
function run_permutations(
    prepared::PermutationTestData,
    n_permutations::Int,
    threshold::Float64,
    critical_t_values::Array{Float64, 2},
    spatial_connectivity::SparseMatrixCSC{Bool},
    cluster_type::Symbol,
    cluster_statistic::Symbol,
    tail::Symbol,
    min_cluster_size::Int,
    min_num_neighbors::Int,
    random_seed::Union{Int, Nothing} = nothing,
    show_progress::Bool = true
)
    # Set random seed if provided
    if random_seed !== nothing
        rng = MersenneTwister(random_seed)
    else
        rng = Random.GLOBAL_RNG
    end
    
    permutation_max_positive = Float64[]
    permutation_max_negative = Float64[]
    sizehint!(permutation_max_positive, n_permutations)
    sizehint!(permutation_max_negative, n_permutations)
    
    # Progress bar
    if show_progress
        progress = Progress(n_permutations, desc="Permutations: ", showspeed=true)
    end
    
    for perm_idx in 1:n_permutations
        # Shuffle labels
        shuffled_A, shuffled_B = shuffle_labels(prepared, rng)
        
        # Create temporary PermutationTestData with shuffled data
        shuffled_prepared = PermutationTestData(
            prepared.design,
            shuffled_A,
            shuffled_B,
            prepared.participants_A,
            prepared.participants_B,
            prepared.electrodes,
            prepared.time_points,
            prepared.condition_A,
            prepared.condition_B,
            prepared.layout
        )
        
        # Compute t-matrix for shuffled data
        t_matrix_perm, df_matrix_perm = compute_t_matrix(shuffled_prepared)
        
        # Threshold
        mask_pos_perm, mask_neg_perm = threshold_t_matrix_parametric(
            t_matrix_perm, df_matrix_perm, critical_t_values, tail
        )
        
        # Pre-filter masks (same as observed data)
        if min_num_neighbors > 0
            mask_pos_perm = prefilter_mask_by_neighbors(mask_pos_perm, spatial_connectivity, min_num_neighbors)
            mask_neg_perm = prefilter_mask_by_neighbors(mask_neg_perm, spatial_connectivity, min_num_neighbors)
        end
        
        # Find clusters
        pos_clusters_perm, neg_clusters_perm = find_clusters(
            mask_pos_perm, mask_neg_perm, prepared.electrodes, prepared.time_points,
            spatial_connectivity, cluster_type, min_cluster_size
        )
        
        # Compute cluster statistics
        if !isempty(pos_clusters_perm)
            _, pos_stats_perm = compute_cluster_statistics(
                pos_clusters_perm, t_matrix_perm, prepared.electrodes, cluster_statistic
            )
            max_pos = isempty(pos_stats_perm) ? 0.0 : maximum(pos_stats_perm)
        else
            max_pos = 0.0
        end
        
        if !isempty(neg_clusters_perm)
            _, neg_stats_perm = compute_cluster_statistics(
                neg_clusters_perm, t_matrix_perm, prepared.electrodes, cluster_statistic
            )
            max_neg = isempty(neg_stats_perm) ? 0.0 : maximum(neg_stats_perm)
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
        
        # Count permutations with max >= observed (add 1 for standard correction)
        count_exceed = sum(permutation_max .>= cluster_stat) + 1
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
    cluster_permutation_test(prepared::PermutationTestData;
                             n_permutations::Int = 1000,
                             threshold::Float64 = 0.05,
                             threshold_method::Symbol = :parametric,
                             cluster_type::Symbol = :spatiotemporal,
                             cluster_statistic::Symbol = :sum,
                             min_cluster_size::Int = 0,
                             min_num_neighbors::Int = 0,
                             tail::Symbol = :both,
                             random_seed::Union{Int, Nothing} = nothing,
                             show_progress::Bool = true)

Perform cluster-based permutation test on prepared ERP data.

# Arguments
- `prepared::PermutationTestData`: Prepared data from `prepare_permutation_data`
- `n_permutations::Int`: Number of permutations (default: 1000)
- `threshold::Float64`: P-value threshold (default: 0.05)
- `threshold_method::Symbol`: Threshold method - `:parametric` (default), `:nonparametric_individual`, or `:nonparametric_common`
- `cluster_type::Symbol`: Type of clustering - `:spatial`, `:temporal`, or `:spatiotemporal` (default)
- `cluster_statistic::Symbol`: Cluster statistic - `:sum` (default), `:max`, `:size`, or `:wcm`
- `min_cluster_size::Int`: Minimum number of points in cluster (default: 0)
- `min_num_neighbors::Int`: Minimum number of neighboring significant channels required (FieldTrip's minNumChannels, default: 0). Points with fewer neighbors are removed before clustering.
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`
- `random_seed::Union{Int, Nothing}`: Random seed for reproducibility (default: nothing)
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `ClusterPermutationResult`: Complete results structure

# Examples
```julia
# Prepare data
prepared = prepare_permutation_data("erps_good", 1, 2, design=:paired, input_dir="data/")

# Run permutation test
result = cluster_permutation_test(prepared, n_permutations=1000, threshold=0.05)

# Access results
println("Found ", length(result.positive_clusters), " positive clusters")
println("Found ", length(result.negative_clusters), " negative clusters")
```
"""
function cluster_permutation_test(
    prepared::PermutationTestData;
    n_permutations::Int = 1000,
    threshold::Float64 = 0.05,
    threshold_method::Symbol = :parametric,
    cluster_type::Symbol = :spatiotemporal,
    cluster_statistic::Symbol = :sum,
    min_cluster_size::Int = 0,
    min_num_neighbors::Int = 0,
    tail::Symbol = :both,
    random_seed::Union{Int, Nothing} = nothing,
    show_progress::Bool = true
)
    # Validate inputs
    validate_permutation_inputs(
        prepared, n_permutations, threshold, cluster_type, cluster_statistic, tail, min_cluster_size
    )
    
    # Only parametric thresholding is implemented for now
    if threshold_method != :parametric
        error("Only :parametric thresholding is currently implemented. " *
              "Got :$threshold_method")
    end
    
    # Compute observed t-matrix and df-matrix
    @info "Computing t-statistics..."
    t_matrix, df_matrix = compute_t_matrix(prepared)
    
    # Compute critical t-values
    @info "Computing critical t-values..."
    critical_t_values = compute_critical_t_values(df_matrix, threshold, tail)
    
    # Threshold observed data
    @info "Thresholding observed data..."
    mask_positive, mask_negative = threshold_t_matrix_parametric(
        t_matrix, df_matrix, critical_t_values, tail
    )
    
    # Build connectivity matrix
    @info "Building connectivity matrix..."
    spatial_connectivity, _, _ = build_connectivity_matrix(
        prepared.electrodes, prepared.layout, cluster_type
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
        mask_positive, mask_negative, prepared.electrodes, prepared.time_points,
        spatial_connectivity, cluster_type, min_cluster_size
    )
    
    # Compute observed cluster statistics
    cluster_stats_positive = Float64[]
    cluster_stats_negative = Float64[]
    
    if !isempty(positive_clusters)
        @info "Computing cluster statistics for $(length(positive_clusters)) positive clusters..."
        positive_clusters, cluster_stats_positive = compute_cluster_statistics(
            positive_clusters, t_matrix, prepared.electrodes, cluster_statistic
        )
    end
    
    if !isempty(negative_clusters)
        @info "Computing cluster statistics for $(length(negative_clusters)) negative clusters..."
        negative_clusters, cluster_stats_negative = compute_cluster_statistics(
            negative_clusters, t_matrix, prepared.electrodes, cluster_statistic
        )
    end
    
    # Run permutations
    @info "Running $n_permutations permutations..."
    permutation_max_positive, permutation_max_negative = run_permutations(
        prepared, n_permutations, threshold, critical_t_values, spatial_connectivity,
        cluster_type, cluster_statistic, tail, min_cluster_size, min_num_neighbors,
        random_seed, show_progress
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
    
    # Assemble results
    result = ClusterPermutationResult(
        prepared.design,
        prepared.design,
        t_matrix,
        df_matrix,
        mask_positive,
        mask_negative,
        positive_clusters,
        negative_clusters,
        cluster_stats_positive,
        cluster_stats_negative,
        n_permutations,
        permutation_max_positive,
        permutation_max_negative,
        random_seed,
        prepared.electrodes,
        prepared.time_points,
        threshold,
        threshold_method,
        cluster_type,
        cluster_statistic,
        min_cluster_size,
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
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]
- `df_matrix::Array{Float64, 2}`: Degrees of freedom [electrodes × time]
- `significant_mask_positive::BitArray{2}`: Positive significant points [electrodes × time]
- `significant_mask_negative::BitArray{2}`: Negative significant points [electrodes × time]
- `correction_method::Symbol`: Multiple comparison correction method used
- `alpha::Float64`: Significance threshold used
- `tail::Symbol`: Test tail used
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points in seconds
"""
struct AnalyticTTestResult
    test_type::Symbol
    t_matrix::Array{Float64, 2}
    p_matrix::Array{Float64, 2}
    df_matrix::Array{Float64, 2}
    significant_mask_positive::BitArray{2}
    significant_mask_negative::BitArray{2}
    correction_method::Symbol
    alpha::Float64
    tail::Symbol
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
end

"""
    compute_p_matrix(t_matrix::Array{Float64, 2}, df_matrix::Array{Float64, 2}, tail::Symbol = :both)

Compute p-values from t-statistics.

# Arguments
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df_matrix::Array{Float64, 2}`: Degrees of freedom [electrodes × time]
- `tail::Symbol`: Test tail - `:both` (two-tailed), `:left`, or `:right` (default: `:both`)

# Returns
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]
"""
function compute_p_matrix(t_matrix::Array{Float64, 2}, df_matrix::Array{Float64, 2}, tail::Symbol = :both)
    p_matrix = Array{Float64, 2}(undef, size(t_matrix))
    
    for i in eachindex(t_matrix)
        t_val = t_matrix[i]
        df = df_matrix[i]
        
        if isnan(t_val) || isnan(df) || isinf(t_val) || isinf(df) || df <= 0
            p_matrix[i] = NaN
            continue
        end
        
        dist = TDist(df)
        if tail == :both
            p_matrix[i] = 2 * (1 - cdf(dist, abs(t_val)))
        elseif tail == :left
            p_matrix[i] = cdf(dist, t_val)
        elseif tail == :right
            p_matrix[i] = 1 - cdf(dist, t_val)
        else
            error("tail must be :both, :left, or :right, got :$tail")
        end
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
            for i in eachindex(p_matrix)
                corrected_mask[i] = !isnan(p_matrix[i]) && p_matrix[i] <= bonferroni_alpha
            end
        end
        
    else
        error("Unsupported correction method: $method. Use :no or :bonferroni")
    end
    
    return corrected_mask
end

"""
    analytic_ttest(prepared::PermutationTestData;
                  alpha::Float64 = 0.05,
                  tail::Symbol = :both,
                  correction_method::Symbol = :no)

Perform analytic (parametric) t-test without permutation (FieldTrip's 'analytic' method).

# Arguments
- `prepared::PermutationTestData`: Prepared data from `prepare_permutation_data`
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
    prepared::PermutationTestData;
    alpha::Float64 = 0.05,
    tail::Symbol = :both,
    correction_method::Symbol = :no
)
    # Validate correction method
    if !(correction_method in (:no, :bonferroni))
        error("correction_method must be :no or :bonferroni. Got :$correction_method")
    end
    
    # Compute t-statistics and degrees of freedom
    @info "Computing t-statistics..."
    t_matrix, df_matrix = compute_t_matrix(prepared)
    
    # Compute p-values
    @info "Computing p-values..."
    p_matrix = compute_p_matrix(t_matrix, df_matrix, tail)
    
    # Apply multiple comparison correction
    @info "Applying $(correction_method) correction for multiple comparisons..."
    corrected_mask = apply_multiple_comparison_correction(p_matrix, alpha, correction_method)
    
    # Create positive and negative masks based on t-values
    mask_positive = BitArray{2}(undef, size(t_matrix))
    mask_negative = BitArray{2}(undef, size(t_matrix))
    
    for i in eachindex(t_matrix)
        t_val = t_matrix[i]
        is_sig = corrected_mask[i]
        
        if isnan(t_val) || !is_sig
            mask_positive[i] = false
            mask_negative[i] = false
        else
            if tail == :both
                mask_positive[i] = t_val > 0
                mask_negative[i] = t_val < 0
            elseif tail == :right
                mask_positive[i] = true
                mask_negative[i] = false
            elseif tail == :left
                mask_positive[i] = false
                mask_negative[i] = true
            end
        end
    end
    
    result = AnalyticTTestResult(
        prepared.design,
        t_matrix,
        p_matrix,
        df_matrix,
        mask_positive,
        mask_negative,
        correction_method,
        alpha,
        tail,
        prepared.electrodes,
        prepared.time_points
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
    
    n_sig_pos = count(result.significant_mask_positive)
    n_sig_neg = count(result.significant_mask_negative)
    
    println(io, "AnalyticTTestResult")
    println(io, "├─ Test type: $(result.test_type)")
    println(io, "├─ Alpha: $(result.alpha)")
    println(io, "├─ Tail: $(result.tail)")
    println(io, "├─ Correction method: $(result.correction_method)")
    println(io, "├─ Data dimensions: $n_electrodes electrodes × $n_time_points time points ($time_range)")
    println(io, "└─ Significant points: $n_sig_pos positive, $n_sig_neg negative")
end
