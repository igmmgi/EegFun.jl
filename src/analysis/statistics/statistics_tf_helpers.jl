# ===================
# TF DATA PREPARATION (4D: participants × electrodes × frequencies × time)
# ===================

"""
    prepare_stats(tfs::Vector{TimeFreqData}; design, condition_selection, channel_selection, frequency_selection, interval_selection, baseline_interval, baseline_method)

Prepare TimeFreqData for comparing two conditions in statistical tests.

Organizes TimeFreqData power values into participant × electrode × frequency × time arrays
for statistical analysis. Validates the design and ensures data consistency across conditions.

# Arguments
- `tfs::Vector{TimeFreqData}`: TF data containing data for multiple conditions/participants
- `design::Symbol`: Design type - `:paired` or `:independent`
- `condition_selection::Function`: Predicate to select exactly 2 conditions (default: `conditions()` - all conditions)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()` - all channels)
- `frequency_selection::Interval`: Frequency range as tuple (e.g., `(4.0, 30.0)`) or `nothing` for all frequencies
- `interval_selection::Interval`: Time interval as tuple (e.g., `(0.0, 1.0)`) or `nothing` for all time points
- `baseline_interval::Union{Nothing,Tuple{Real,Real}}`: Baseline window in seconds (e.g., `(-0.5, -0.1)`). Default: `nothing` (no baseline)
- `baseline_method::Symbol`: Baseline method (default: `:db`). Options: `:db`, `:absolute`, `:relative`, `:relchange`, `:percent`, `:zscore`

# Returns
- `TFStatisticalData`: Prepared data structure ready for statistical testing

# Examples
```julia
# Prepare with baseline correction
prepared = prepare_stats(tfs; design=:paired,
    baseline_interval=(-0.5, -0.1), baseline_method=:db,
    frequency_selection=(4.0, 30.0))

# Run permutation test
result = permutation_test(prepared; n_permutations=1000)
```
"""
function prepare_stats(
    tfs::Vector{TimeFreqData};
    design::Symbol = :paired,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    frequency_selection::Interval = nothing,
    interval_selection::Interval = nothing,
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    baseline_method::Symbol = :db,
)
    # Group all TF data by condition
    tfs_by_condition = group_by_condition(tfs)

    # Apply condition selection to the sorted condition numbers
    all_cond_nums = collect(keys(tfs_by_condition))
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    # Validate exactly 2 conditions
    length(selected_cond_nums) == 2 ||
        @minimal_error "Statistical tests require exactly 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select exactly 2 conditions."

    condition1 = tfs_by_condition[selected_cond_nums[1]]
    condition2 = tfs_by_condition[selected_cond_nums[2]]

    # Apply baseline correction if requested (copy first to avoid mutating original data)
    if !isnothing(baseline_interval)
        condition1 = [copy(tf) for tf in condition1]
        condition2 = [copy(tf) for tf in condition2]
        for tf in condition1
            tf_baseline!(tf, baseline_interval; method = baseline_method)
        end
        for tf in condition2
            tf_baseline!(tf, baseline_interval; method = baseline_method)
        end
    end

    # Validate design
    design ∉ (:paired, :independent) && @minimal_error "design must be :paired or :independent, got :$design"

    # Extract participant IDs from filenames
    vps1 = [_extract_participant_id(basename(data.file)) for data in condition1]
    vps2 = [_extract_participant_id(basename(data.file)) for data in condition2]

    if design == :paired
        vps1 != vps2 && @minimal_error "Paired design requires same participants in both conditions"
    elseif design == :independent
        (length(vps1) < 2 || length(vps2) < 2) && @minimal_error "Independent design requires at least 2 participants per group"
    end

    # Validate all TF data have same structure (same channels, same frequencies, same time points)
    _have_same_structure(condition1) || @minimal_error("Condition 1: TF data have inconsistent structure")
    _have_same_structure(condition2) || @minimal_error("Condition 2: TF data have inconsistent structure")
    _have_same_structure(condition1[1], condition2[1]) ||
        @minimal_error("Conditions have inconsistent structure (different channels, frequencies, or time points)")

    # Get metadata from first TF dataset 
    ref_tf = condition1[1]
    all_electrodes = channel_labels(ref_tf)
    all_freqs = Float64.(sort(unique(ref_tf.data_power.freq)))
    all_times = Float64.(sort(unique(ref_tf.data_power.time)))

    # Apply channel selection
    selected_ch_mask = channel_selection(1:length(all_electrodes))
    electrodes = all_electrodes[selected_ch_mask]
    isempty(electrodes) && @minimal_error "No channels matched the channel selection."

    # Apply frequency selection
    if !isnothing(frequency_selection) && frequency_selection isa Tuple
        freq_mask = (all_freqs .>= frequency_selection[1]) .& (all_freqs .<= frequency_selection[2])
        frequencies = all_freqs[freq_mask]
    else
        frequencies = all_freqs
    end
    isempty(frequencies) && @minimal_error "No frequencies matched the frequency selection."

    # Apply time (interval) selection
    if !isnothing(interval_selection) && interval_selection isa Tuple
        time_mask = (all_times .>= interval_selection[1]) .& (all_times .<= interval_selection[2])
        time_points = all_times[time_mask]
    else
        time_points = all_times
    end
    isempty(time_points) && @minimal_error "No time points matched the interval selection."

    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    # Extract 4D arrays: [participants × electrodes × frequencies × time]
    data1 = _extract_tf_array(condition1, electrodes, frequencies, time_points)
    data2 = _extract_tf_array(condition2, electrodes, frequencies, time_points)

    # Create grand averages from pre-extracted arrays (avoids re-extracting)
    condition1_avg = _create_tf_grand_average(condition1[1], data1, electrodes, frequencies, time_points, selected_cond_nums[1])
    condition2_avg = _create_tf_grand_average(condition2[1], data2, electrodes, frequencies, time_points, selected_cond_nums[2])

    return TFStatisticalData([condition1_avg, condition2_avg], TFAnalysisData(design, [data1, data2], frequencies, time_points))
end
function prepare_stats(
    ::Type{TimeFreqData},
    file_pattern::String,
    design::Symbol;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    kwargs...,
)
    all_tfs = read_all_data(TimeFreqData, joinpath(input_dir, file_pattern), participant_selection)
    isempty(all_tfs) && @minimal_error "No valid TF data found matching pattern '$file_pattern' in $input_dir"

    return prepare_stats(all_tfs; design = design, kwargs...)
end


# ===================
# TF HELPER FUNCTIONS
# ===================

"""
    _extract_tf_array(tfs, electrodes, frequencies, time_points)

Extract power values from TimeFreqData into a 4D array.

# Returns
- `Array{Float64, 4}`: [participants × electrodes × frequencies × time]
"""
function _extract_tf_array(
    tfs::Vector{TimeFreqData},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
)
    n_participants = length(tfs)
    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    data = Array{Float64,4}(undef, n_participants, n_electrodes, n_freqs, n_time)

    for (p_idx, tf) in enumerate(tfs)
        # Build row index once per participant, reuse for all channels
        row_index = _tf_build_row_index(tf.data_power)

        for (e_idx, elec) in enumerate(electrodes)
            data[p_idx, e_idx, :, :] = _tf_df_to_matrix(tf.data_power, elec, frequencies, time_points, row_index)
        end
    end

    return data
end


"""
    _create_tf_grand_average(ref_tf, data_4d, electrodes, frequencies, time_points, cond_num)

Create a grand average TimeFreqData by averaging a pre-extracted 4D array over participants.

# Arguments
- `ref_tf::TimeFreqData`: Reference TF data for metadata (condition_name, layout, etc.)
- `data_4d::Array{Float64,4}`: Pre-extracted [participants × electrodes × freqs × time] array
- `electrodes`, `frequencies`, `time_points`: Axis labels
- `cond_num::Int`: Condition number for the grand average
"""
function _create_tf_grand_average(
    ref_tf::TimeFreqData,
    data_4d::Array{Float64,4},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    cond_num::Int,
)
    # Average over participants (dim 1)
    avg_power = dropdims(mean(data_4d, dims = 1), dims = 1)  # [electrodes × freqs × time]

    # Build grand average DataFrame with same structure as input
    n_freqs = length(frequencies)
    n_time = length(time_points)

    # Create time and freq columns (time varies fastest within each freq)
    time_col = repeat(time_points, outer = n_freqs)
    freq_col = repeat(frequencies, inner = n_time)

    # Build DataFrame with zero-initialized channel columns
    power_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
    for elec in electrodes
        power_df[!, elec] = zeros(n_freqs * n_time)
    end

    # Fill using shared helper (avoids fragile transpose+vec pattern)
    row_index = _tf_build_row_index(power_df)
    for (e_idx, elec) in enumerate(electrodes)
        _tf_matrix_to_df!(power_df, elec, avg_power[e_idx, :, :], frequencies, time_points, row_index)
    end

    # Create a matching phase DataFrame (all zeros for grand average since phase averages to zero)
    phase_df = copy(power_df)
    for elec in electrodes
        phase_df[!, elec] .= 0.0
    end

    return TimeFreqData(
        "grand_average",
        cond_num,
        ref_tf.condition_name,
        power_df,
        phase_df,
        copy(ref_tf.layout),
        ref_tf.sample_rate,
        ref_tf.method,
        ref_tf.baseline,
        ref_tf.analysis_info,
    )
end


# ===================
# TF T-TEST COMPUTATION (4D)
# ===================

"""
    _compute_t_matrix_tf(data1, data2, design; tail)

Compute t-statistics and p-values for all electrode × frequency × time points.

# Arguments
- `data1::Array{Float64, 4}`: Condition 1 data [participants × electrodes × freqs × time]
- `data2::Array{Float64, 4}`: Condition 2 data [participants × electrodes × freqs × time]
- `design::Symbol`: `:paired` or `:independent`
- `tail::Symbol`: `:both`, `:left`, or `:right`

# Returns
- `t_matrix::Array{Float64, 3}`: T-statistics [electrodes × freqs × time]
- `df::Float64`: Degrees of freedom
- `p_matrix::Array{Float64, 3}`: P-values [electrodes × freqs × time]
"""
function _compute_t_matrix_tf(data1::Array{Float64,4}, data2::Array{Float64,4}, design::Symbol; tail::Symbol = :both, compute_p_values::Bool = true, t_matrix_buffer::Union{Nothing,Array{Float64,3}} = nothing)
    n_participants, n_electrodes, n_freqs, n_time = size(data1)
    
    if !isnothing(t_matrix_buffer)
        t_matrix = t_matrix_buffer
    else
        t_matrix = Array{Float64,3}(undef, n_electrodes, n_freqs, n_time)
    end
    
    p_matrix = compute_p_values ? Array{Float64,3}(undef, n_electrodes, n_freqs, n_time) : Array{Float64,3}(undef, 0, 0, 0)

    if design == :paired
        @inbounds for t_idx = 1:n_time
            for f_idx = 1:n_freqs
                for e_idx = 1:n_electrodes
                    sum_diff = 0.0
                    sum_diff_sq = 0.0
                    for p_idx = 1:n_participants
                        diff_val = data1[p_idx, e_idx, f_idx, t_idx] - data2[p_idx, e_idx, f_idx, t_idx]
                        sum_diff += diff_val
                        sum_diff_sq += diff_val * diff_val
                    end
                    mean_diff_val = sum_diff / n_participants
                    variance = (sum_diff_sq / n_participants) - (mean_diff_val * mean_diff_val)
                    std_diff = sqrt(variance * n_participants / (n_participants - 1))

                    if std_diff == 0.0
                        t_matrix[e_idx, f_idx, t_idx] = mean_diff_val == 0.0 ? NaN : Inf
                    else
                        t_matrix[e_idx, f_idx, t_idx] = mean_diff_val / (std_diff / sqrt(n_participants))
                    end
                end
            end
        end

        df = Float64(n_participants - 1)
        if compute_p_values
            _compute_p_matrix_tf!(p_matrix, t_matrix, df, tail)
        end

    else
        # Independent design
        n_A = size(data1, 1)
        n_B = size(data2, 1)
        df = Float64(n_A + n_B - 2)
        @inbounds for t_idx = 1:n_time
            for f_idx = 1:n_freqs
                for e_idx = 1:n_electrodes
                    sum1 = 0.0; sum2 = 0.0; sum1_sq = 0.0; sum2_sq = 0.0
                    
                    @simd for p_idx = 1:n_A
                        val = data1[p_idx, e_idx, f_idx, t_idx]
                        sum1 += val
                        sum1_sq += val * val
                    end
                    
                    @simd for p_idx = 1:n_B
                        val = data2[p_idx, e_idx, f_idx, t_idx]
                        sum2 += val
                        sum2_sq += val * val
                    end
                    
                    mean1 = sum1 / n_A
                    mean2 = sum2 / n_B
                    var1 = (sum1_sq / n_A - mean1 * mean1) * n_A / (n_A - 1)
                    var2 = (sum2_sq / n_B - mean2 * mean2) * n_B / (n_B - 1)
                    
                    pooled_var = ((n_A - 1) * var1 + (n_B - 1) * var2) / df
                    se = sqrt(pooled_var * (1.0 / n_A + 1.0 / n_B))
                    
                    t_matrix[e_idx, f_idx, t_idx] = (mean1 - mean2) / se
                end
            end
        end
        
        if compute_p_values
            _compute_p_matrix_tf!(p_matrix, t_matrix, df, tail)
        end
    end

    return t_matrix, df, p_matrix
end

function _compute_t_matrix_tf(prepared::TFStatisticalData; tail::Symbol = :both)
    return _compute_t_matrix_tf(prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, tail = tail)
end


"""
    _compute_p_matrix_tf!(p_matrix, t_matrix, df, tail)

In-place computation of p-values from 3D t-statistics matrix.
"""
function _compute_p_matrix_tf!(p_matrix::Array{Float64,3}, t_matrix::Array{Float64,3}, df::Float64, tail::Symbol)
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


"""
    _compute_critical_t_values_tf(df, matrix_size, alpha, tail)

Compute critical t-values for parametric thresholding (3D version).

# Returns
- `critical_t_values::Array{Float64, 3}`: Uniform critical t-values [electrodes × freqs × time]
"""
function _compute_critical_t_values_tf(df::Float64, matrix_size::Tuple{Int,Int,Int}, alpha::Float64 = 0.05, tail::Symbol = :both)
    critical_t_values = Array{Float64,3}(undef, matrix_size)

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
