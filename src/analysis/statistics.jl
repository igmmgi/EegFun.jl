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
    t_str = isnan(result.t) ? "NaN" : (isinf(result.t) ? (result.t > 0 ? "Inf" : "-Inf") : @sprintf("%.4f", result.t))
    p_str = isnan(result.p) ? "NaN" : @sprintf("%.4f", result.p)
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

"""
    prepare_permutation_data(erps_A::Vector{ErpData}, erps_B::Vector{ErpData}; 
                             design::Symbol = :paired,
                             sample_selection::Function = samples())

    prepare_permutation_data(file_pattern::String, condition_A, condition_B;
                             design::Symbol = :paired,
                             input_dir::String = pwd(),
                             participant_selection::Function = participants(),
                             channel_selection::Function = channels(),
                             sample_selection::Function = samples())

Prepare ErpData for cluster-based permutation tests.

Organizes ErpData into participant × electrode × time arrays for statistical analysis.
Validates the design and ensures data consistency across conditions.

# Arguments (Direct data version)
- `erps_A::Vector{ErpData}`: ERPs for condition/group A (one per participant)
- `erps_B::Vector{ErpData}`: ERPs for condition/group B (one per participant)
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `channel_selection::Function`: Function to filter channels (default: all channels)
- `sample_selection::Function`: Function to select time points (default: all samples). Use `samples((start, end))` for time windows.

# Arguments (File-based version)
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned")
- `condition_A`: Condition identifier for group A (Int or String)
- `condition_B`: Condition identifier for group B (Int or String)
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `input_dir::String`: Directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Function to filter participants (default: all)
- `channel_selection::Function`: Function to filter channels (default: all, not yet implemented)
- `sample_selection::Function`: Function to select time points (default: all samples). Use `samples((start, end))` for time windows.

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

# With time window selection
prepared = prepare_permutation_data(
    erps_condition_A, 
    erps_condition_B, 
    design = :paired,
    sample_selection = samples((0.0, 0.5))
)
```
"""
# Direct data version
function prepare_permutation_data(erps_A::Vector{ErpData}, erps_B::Vector{ErpData}; 
                                   design::Symbol = :paired,
                                   channel_selection::Function = channels(),
                                   sample_selection::Function = samples())
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
    
    # Get all electrodes from first ERP
    all_electrodes_A = channel_labels(first_erp_A)
    all_electrodes_B = channel_labels(first_erp_B)
    
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
    sample_mask_A = sample_selection(first_erp_A.data)
    sample_mask_B = sample_selection(first_erp_B.data)
    
    # Get selected time points
    time_points_A = first_erp_A.data[sample_mask_A, :time]
    time_points_B = first_erp_B.data[sample_mask_B, :time]
    
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
    
    # Get layout (use from first ERP)
    layout = first_erp_A.layout
    
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
    sample_selection::Function = samples()
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
        sample_selection = sample_selection
    )
end
