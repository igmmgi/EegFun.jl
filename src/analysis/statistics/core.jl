# This file contains the core data preparation and t-test computation functions
# for statistical testing of EEG/ERP data.

"""
    prepare_stats(erps::Vector{ErpData}; design::Symbol = :paired, condition_selection::Function = conditions([1, 2]), channel_selection::Function = channels(), interval_selection::Interval = times(), baseline_interval::Interval = nothing, analysis_interval::Interval = times())

Prepare ErpData for comparing two conditions in statistical tests (permutation and analytic tests).

Organizes ErpData into participant × electrode × time arrays for statistical analysis.
Validates the design and ensures data consistency across conditions.

# Arguments
- `erps::Vector{ErpData}`: ERPs containing data for multiple conditions/participants
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `condition_selection::Function`: Predicate to select exactly 2 conditions for comparison (default: `conditions([1, 2])`)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()` - all channels)
- `interval_selection::Interval`: Time interval as tuple (e.g., (0.0, 1.0)) or interval object for initial data selection (default: nothing - all samples)
- `baseline_interval::Interval`: Baseline interval as tuple (e.g., (-0.2, 0.0)) or interval object (default: nothing - no baseline)
- `analysis_interval::Interval`: Analysis interval as tuple (e.g., (0.3, 0.5)) or interval object for statistical testing (default: nothing - use interval_selection)

# Returns
- `StatisticalData`: Prepared data structure ready for statistical testing
"""
function prepare_stats(
    erps::Vector{ErpData};
    design::Symbol = :paired,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    analysis_interval::Interval = times(),
)
    # Group all ERPs by condition first
    erps_by_condition = group_by_condition(erps)

    # Apply condition selection to the sorted condition numbers
    all_cond_nums = collect(keys(erps_by_condition))  # Already sorted by group_by_condition
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    # Validate exactly 2 conditions
    length(selected_cond_nums) == 2 ||
        @minimal_error "Statistical tests require exactly 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select exactly 2 conditions."

    condition1 = erps_by_condition[selected_cond_nums[1]]
    condition2 = erps_by_condition[selected_cond_nums[2]]

    # Validate design
    design ∉ (:paired, :independent) && @minimal_error "design must be :paired or :independent, got :$design"

    # Extract participant IDs from filenames (using utility from batch.jl)
    vps1 = [_extract_participant_id(basename(data.file)) for data in condition1]
    vps2 = [_extract_participant_id(basename(data.file)) for data in condition2]

    # Validate design
    if design == :paired # Paired design: same participants in both conditions, in the same order
        vps1 != vps2 && @minimal_error "Paired design requires same participants in both conditions"
    elseif design == :independent # Independent design: different participants (or allow overlap)
        (length(vps1) < 2 || length(vps2) < 2) && @minimal_error "Independent design requires at least 2 participants per group"
    end

    # Validate all ERPs have same structure within each condition
    _have_same_structure(condition1) || @minimal_error("Condition 1: ERPs have inconsistent structure")
    _have_same_structure(condition2) || @minimal_error("Condition 2: ERPs have inconsistent structure")
    _have_same_structure(condition1[1], condition2[1]) || @minimal_error("Condition 1 vs. 2: ERPs have inconsistent structure")

    # Convert intervals to samples() predicates for subset()
    sample_sel = _interval_to_samples(interval_selection)

    condition1 = subset(condition1; channel_selection = channel_selection, sample_selection = sample_sel)
    isempty(condition1) && @minimal_error "No data matched the selection criteria!"

    condition2 = subset(condition2; channel_selection = channel_selection, sample_selection = sample_sel)
    isempty(condition2) && @minimal_error "No data matched the selection criteria!"

    # baseline (if no interval specified, use the full time range)
    if isnothing(baseline_interval)
        baseline!.(condition1)
        baseline!.(condition2)
    else
        baseline!.(condition1, Ref(baseline_interval))
        baseline!.(condition2, Ref(baseline_interval))
    end

    # create grand averages for ease of use in plotting results
    condition1_avg = _create_grand_average(condition1, selected_cond_nums[1])
    condition2_avg = _create_grand_average(condition2, selected_cond_nums[2])

    # Compute per-condition SE over the full display interval (before analysis subset narrows it)
    # SE = std(across participants) / sqrt(n) at each electrode × time point
    display_electrodes = channel_labels(condition1[1])
    n_display_electrodes = length(display_electrodes)
    n_display_time = nrow(condition1[1].data)
    n_cond1 = length(condition1)
    n_cond2 = length(condition2)

    # Stack per-participant data: [participants × electrodes × time]
    cond1_stacked = _extract_erp_array(condition1, display_electrodes, n_display_time)
    cond2_stacked = _extract_erp_array(condition2, display_electrodes, n_display_time)

    # SE = std across participants (dim 1) / sqrt(n), result is [electrodes × time]
    se_cond1 = dropdims(std(cond1_stacked, dims = 1), dims = 1) ./ sqrt(n_cond1)
    se_cond2 = dropdims(std(cond2_stacked, dims = 1), dims = 1) ./ sqrt(n_cond2)

    # SE of the mean difference (for paired: std of per-participant differences / sqrt(n))
    if design == :paired
        diff_stacked = cond1_stacked .- cond2_stacked  # [participants × electrodes × time]
        se_diff = dropdims(std(diff_stacked, dims = 1), dims = 1) ./ sqrt(n_cond1)
    else
        # Independent: SE of difference = sqrt(SE_A^2 + SE_B^2)
        se_diff = sqrt.(se_cond1 .^ 2 .+ se_cond2 .^ 2)
    end

    # create second subset with analysis_interval for statistical tests
    # Use analysis_interval if provided, otherwise use interval_selection
    analysis_sel = if !isnothing(analysis_interval)
        _interval_to_samples(analysis_interval)
    else
        sample_sel
    end

    condition1 = subset(condition1; channel_selection = channel_selection, sample_selection = analysis_sel)
    isempty(condition1) && @minimal_error "No data matched the analysis interval criteria!"

    condition2 = subset(condition2; channel_selection = channel_selection, sample_selection = analysis_sel)
    isempty(condition2) && @minimal_error "No data matched the analysis interval criteria!"

    # Get dimensions and metadata from analysis subset
    electrodes = channel_labels(condition1[1])
    n_electrodes = length(electrodes)
    time_points = condition1[1].data[!, :time]
    n_time = length(time_points)

    # Extract data arrays: [participants × electrodes × time]
    condition1_array = _extract_erp_array(condition1, electrodes, n_time)
    condition2_array = _extract_erp_array(condition2, electrodes, n_time)

    return StatisticalData(
        [condition1_avg, condition2_avg],
        AnalysisData(design, [condition1_array, condition2_array], time_points),
        se_cond1,
        se_cond2,
        se_diff,
    )

end

"""Load ERP/TF files by pattern and call the typed `prepare_stats` method."""
function prepare_stats(
    file_pattern::String,
    design::Symbol;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    kwargs...,
)
    # Load data and auto-detect type
    all_data = read_all_data(joinpath(input_dir, file_pattern), participant_selection)
    isempty(all_data) && @minimal_error "No valid data found matching pattern '$file_pattern' in $input_dir"

    first_item = first(all_data)
    if first_item isa TimeFreqData
        return prepare_stats(Vector{TimeFreqData}(all_data); design = design, kwargs...)
    elseif first_item isa ErpData
        return prepare_stats(Vector{ErpData}(all_data); design = design, kwargs...)
    else
        @minimal_error "Unsupported data type: $(typeof(first_item)). Expected ErpData or TimeFreqData."
    end
end

"""
    _extract_erp_array(erps::Vector{ErpData}, electrodes::Vector{Symbol}, n_time::Int)

Extract ERP data into a 3D array [participants × electrodes × time].
Pre-allocates the target array to avoid slow Splatting/Catting operations.
"""
function _extract_erp_array(erps::Vector{ErpData}, electrodes::Vector{Symbol}, n_time::Int)
    n_participants = length(erps)
    n_electrodes = length(electrodes)
    data = Array{Float64,3}(undef, n_participants, n_electrodes, n_time)

    for (p, erp) in enumerate(erps)
        mat = Matrix{Float64}(erp.data[!, electrodes])
        for e = 1:n_electrodes
            @inbounds for t = 1:n_time
                data[p, e, t] = mat[t, e]
            end
        end
    end

    return data
end


"""
    _validate_permutation_inputs(prepared::StatisticalData, 
                                n_permutations::Int,
                                threshold::Float64,
                                cluster_type::Symbol,
                                tail::Symbol)

Validate inputs for cluster permutation test.

# Arguments
- `prepared::StatisticalData`: Prepared data
- `n_permutations::Int`: Number of permutations
- `threshold::Float64`: P-value threshold
- `cluster_type::Symbol`: Type of clustering
- `tail::Symbol`: Test tail
"""
function _validate_permutation_inputs(
    prepared::StatisticalData,
    n_permutations::Int,
    threshold::Float64,
    cluster_type::Symbol,
    tail::Symbol,
)
    # Some basic validations 
    n_permutations < 1 && @minimal_error "n_permutations must be >= 1, got $n_permutations"
    (threshold <= 0.0 || threshold >= 1.0) && @minimal_error "threshold must be between 0 and 1, got $threshold"
    cluster_type ∉ (:spatial, :temporal, :spatiotemporal) &&
        @minimal_error "cluster_type must be :spatial, :temporal, or :spatiotemporal, got :$cluster_type"
    tail ∉ (:both, :left, :right) && @minimal_error "tail must be :both, :left, or :right, got :$tail"

    # Auto-compute neighbours if missing (only when needed for spatial clustering)
    if cluster_type ∈ (:spatial, :spatiotemporal)
        if isnothing(prepared.data[1].layout.neighbours)
            @minimal_warning "Layout.neighbours is not set. Computing with default distance criterion (0.35)."
            # Need to update layout.neighbours for both conditions
            get_neighbours_xy!(prepared.data[1].layout, 0.35)
            if length(prepared.data) > 1
                # Ensure they match if they have separate layouts
                # They should share the layout anyway!
                get_neighbours_xy!(prepared.data[2].layout, 0.35)
            end
        end
    end

    return nothing
end

# === T-TEST COMPUTATION ===

"""
    compute_t_matrix(prepared::StatisticalData; tail::Symbol = :both)

Compute t-statistics and p-values for all electrode × time points.

# Arguments
- `prepared::StatisticalData`: Prepared data for statistical test
- `tail::Symbol`: Test tail - `:both` (default), `:left`, or `:right`

# Returns
- `t_matrix::Array{Float64, 2}`: T-statistics [electrodes × time]
- `df::Float64`: Degrees of freedom (constant across all electrodes/time points)
- `p_matrix::Array{Float64, 2}`: P-values [electrodes × time]
- `se_matrix::Array{Float64, 2}`: Standard error of the mean difference [electrodes × time]

# Examples
```julia
t_matrix, df, p_matrix = compute_t_matrix(prepared, tail=:both)
```
"""
# Optimized version that accepts raw arrays (for permutation loop - avoids StatisticalData creation)
"""Compute paired or independent t-statistics across all electrode × time points (raw-array optimised path)."""
function _compute_t_matrix(
    data1::Array{Float64,3},
    data2::Array{Float64,3},
    design::Symbol;
    tail::Symbol = :both,
    mean1_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    mean2_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    mean_diff_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    std_diff_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    t_matrix_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    compute_p_values::Bool = true,
)
    n_participants, n_electrodes, n_time = size(data1)

    if !isnothing(t_matrix_buffer)
        t_matrix = t_matrix_buffer
    else
        t_matrix = Array{Float64,2}(undef, n_electrodes, n_time)
    end

    p_matrix = compute_p_values ? Array{Float64,2}(undef, n_electrodes, n_time) : Array{Float64,2}(undef, 0, 0)

    if design == :paired

        # Compute mean of differences: mean(data1 - data2) = mean(data1) - mean(data2)
        # Use pre-allocated buffers if provided, otherwise allocate
        if !isnothing(mean1_buffer)
            mean1 = mean1_buffer
            mean2 = mean2_buffer
            mean_diff = mean_diff_buffer
        else
            mean1 = Array{Float64,2}(undef, n_electrodes, n_time)
            mean2 = Array{Float64,2}(undef, n_electrodes, n_time)
            mean_diff = Array{Float64,2}(undef, n_electrodes, n_time)
        end

        # Compute means and std of differences in one pass (avoids diff array allocation)
        if !isnothing(std_diff_buffer)
            std_diff = std_diff_buffer
        else
            std_diff = Array{Float64,2}(undef, n_electrodes, n_time)
        end

        @inbounds for t_idx = 1:n_time
            @inbounds for e_idx = 1:n_electrodes
                sum1 = 0.0
                sum2 = 0.0
                sum_diff = 0.0
                sum_diff_sq = 0.0
                for p_idx = 1:n_participants
                    val1 = data1[p_idx, e_idx, t_idx]
                    val2 = data2[p_idx, e_idx, t_idx]
                    diff_val = val1 - val2
                    sum1 += val1
                    sum2 += val2
                    sum_diff += diff_val
                    sum_diff_sq += diff_val * diff_val
                end
                mean1_val = sum1 / n_participants
                mean2_val = sum2 / n_participants
                mean_diff_val = sum_diff / n_participants
                mean_diff[e_idx, t_idx] = mean_diff_val

                # Compute variance: var = mean(x^2) - mean(x)^2, then std = sqrt(var * n/(n-1))
                variance = (sum_diff_sq / n_participants) - (mean_diff_val * mean_diff_val)
                std_diff_val = sqrt(variance * n_participants / (n_participants - 1))
                std_diff[e_idx, t_idx] = std_diff_val

                # Only store mean1/mean2 if buffers were provided
                if !isnothing(mean1_buffer)
                    mean1[e_idx, t_idx] = mean1_val
                    mean2[e_idx, t_idx] = mean2_val
                end

                # Compute SE and t-statistic directly in loop (avoids multiple allocations)
                if std_diff_val == 0.0
                    t_matrix[e_idx, t_idx] = mean_diff_val == 0.0 ? NaN : Inf
                else
                    t_matrix[e_idx, t_idx] = mean_diff_val / (std_diff_val / sqrt(n_participants))
                end
            end
        end

        # Compute SE for returned array
        se_matrix = std_diff ./ sqrt(n_participants)

        # Degrees of freedom (same for all points in paired design)
        df = Float64(n_participants - 1)

        # Compute p-values using internal function (avoids code duplication)
        # Use pre-allocated p_matrix buffer
        if compute_p_values
            p_matrix = _compute_p_matrix(t_matrix, df, tail, p_matrix)
        end

    else
        # Independent design: need to loop (but df is constant across all points)
        se_matrix = Array{Float64,2}(undef, n_electrodes, n_time)
        n_A = size(data1, 1)
        n_B = size(data2, 1)
        df = Float64(n_A + n_B - 2)

        @inbounds for t_idx = 1:n_time
            @inbounds for e_idx = 1:n_electrodes
                sum1 = 0.0;
                sum2 = 0.0;
                sum1_sq = 0.0;
                sum2_sq = 0.0

                @simd for p_idx = 1:n_A
                    val = data1[p_idx, e_idx, t_idx]
                    sum1 += val
                    sum1_sq += val * val
                end

                @simd for p_idx = 1:n_B
                    val = data2[p_idx, e_idx, t_idx]
                    sum2 += val
                    sum2_sq += val * val
                end

                mean1 = sum1 / n_A
                mean2 = sum2 / n_B
                var1 = (sum1_sq / n_A - mean1 * mean1) * n_A / (n_A - 1)
                var2 = (sum2_sq / n_B - mean2 * mean2) * n_B / (n_B - 1)

                # Pooled standard error (assuming equal variances as documented)
                pooled_var = ((n_A - 1) * var1 + (n_B - 1) * var2) / df
                se = sqrt(pooled_var * (1.0 / n_A + 1.0 / n_B))

                t_matrix[e_idx, t_idx] = (mean1 - mean2) / se
                se_matrix[e_idx, t_idx] = se
            end
        end

        if compute_p_values
            p_matrix = _compute_p_matrix(t_matrix, df, tail, p_matrix)
        end
    end

    return t_matrix, df, p_matrix, se_matrix
end

"""Convenience wrapper: extract arrays from `StatisticalData` and call the raw-array `_compute_t_matrix`."""
function _compute_t_matrix(prepared::StatisticalData; tail::Symbol = :both)
    t_matrix, df, p_matrix, se_matrix =
        _compute_t_matrix(prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, tail = tail)
    return t_matrix, df, p_matrix, se_matrix
end

"""
    _compute_critical_t_values(df::Float64, 
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
- `critical_t_values::Array{Float64, 2}`: Critical t-values [electrodes × time]. All values are uniform.

# Examples
```julia
critical_t = _compute_critical_t_values(df, size(t_matrix), 0.05, :both)
```
"""
function _compute_critical_t_values(df::Float64, matrix_size::Tuple{Int,Int}, alpha::Float64 = 0.05, tail::Symbol = :both)
    critical_t_values = Array{Float64,2}(undef, matrix_size)

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

"""
    _compute_p_matrix(t_matrix, df, tail, p_matrix_buffer)

Compute p-values from t-statistics.

Internal helper function used by compute_t_matrix.
"""
function _compute_p_matrix(
    t_matrix::Array{Float64,2},
    df::Float64,
    tail::Symbol,
    p_matrix_buffer::Union{Nothing,Array{Float64,2}} = nothing,
)
    if !isnothing(p_matrix_buffer)
        p_matrix = p_matrix_buffer
    else
        p_matrix = Array{Float64,2}(undef, size(t_matrix))
    end
    dist = TDist(df)
    @inbounds for i in eachindex(t_matrix)
        t_val = t_matrix[i]
        p_matrix[i] = if isnan(t_val) || isinf(t_val)
            NaN
        elseif tail == :both
            2 * (1 - cdf(dist, abs(t_val)))
        elseif tail == :left
            cdf(dist, t_val)
        else  # :right
            1 - cdf(dist, t_val)
        end
    end
    return p_matrix
end
