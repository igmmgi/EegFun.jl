# This file contains the core data preparation and t-test computation functions
# for statistical testing of EEG/ERP data.

"""
    prepare_stats(erps::Vector{ErpData}; design::Symbol = :paired, condition_selection::Function = conditions([1, 2]), channel_selection::Function = channels(), interval_selection::Interval = times(), baseline_interval::Interval = times(), analysis_interval::Interval = times())

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
    baseline_interval::Interval = times(),
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
        @minimal_error_throw "Statistical tests require exactly 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select exactly 2 conditions."

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
        length(vps1) < 2 || length(vps2) < 2 && @minimal_error "Independent design requires at least 2 participants per group"
    end

    # Validate all ERPs have same structure within each condition
    _have_same_structure(condition1) || @minimal_error("Condition 1: ERPs have inconsistent structure")
    _have_same_structure(condition2) || @minimal_error("Condition 2: ERPs have inconsistent structure")
    _have_same_structure(condition1[1], condition2[1]) || @minimal_error("Condition 1 vs. 2: ERPs have inconsistent structure")

    # Convert intervals to samples() predicates for subset()
    sample_sel = _interval_to_samples(interval_selection)

    condition1 = subset(condition1; channel_selection = channel_selection, sample_selection = sample_sel)
    isempty(condition1) && @minimal_error_throw "No data matched the selection criteria!"

    condition2 = subset(condition2; channel_selection = channel_selection, sample_selection = sample_sel)
    isempty(condition2) && @minimal_error_throw "No data matched the selection criteria!"

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

    # create second subset with analysis_interval for statistical tests
    # Use analysis_interval if provided, otherwise use interval_selection
    analysis_sel = if !isnothing(analysis_interval)
        _interval_to_samples(analysis_interval)
    else
        sample_sel
    end

    condition1 = subset(condition1; channel_selection = channel_selection, sample_selection = analysis_sel)
    isempty(condition1) && @minimal_error_throw "No data matched the analysis interval criteria!"

    condition2 = subset(condition2; channel_selection = channel_selection, sample_selection = analysis_sel)
    isempty(condition2) && @minimal_error_throw "No data matched the analysis interval criteria!"

    # Get dimensions and metadata from analysis subset
    electrodes = channel_labels(condition1[1])
    n_electrodes = length(electrodes)
    time_points = condition1[1].data[!, :time]
    n_time = length(time_points)

    # Extract data arrays: [participants × electrodes × time]
    condition1 = cat([reshape(Matrix(erp.data[!, electrodes])', 1, n_electrodes, n_time) for erp in condition1]..., dims = 1)
    condition2 = cat([reshape(Matrix(erp.data[!, electrodes])', 1, n_electrodes, n_time) for erp in condition2]..., dims = 1)

    return StatisticalData([condition1_avg, condition2_avg], AnalysisData(design, [condition1, condition2], time_points))

end

"""
    prepare_stats(file_pattern::String, design::Symbol; input_dir::String = pwd(), participant_selection::Function = participants(), condition_selection::Function = conditions([1, 2]), channel_selection::Function = channels(), interval_selection::Interval = times(), baseline_interval::Interval = times(), analysis_interval::Interval = times())

Prepare ErpData for statistical tests from JLD2 files (convenience wrapper).

Loads ErpData from JLD2 files matching the pattern and prepares them for statistical testing.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned")
- `design::Symbol`: Design type - `:paired` (same participants in both conditions) or `:independent` (different participants)
- `input_dir::String`: Directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Predicate to filter participants (default: `participants()` - all participants)
- `condition_selection::Function`: Predicate to select exactly 2 conditions for comparison (default: `conditions([1, 2])`)
- `channel_selection::Function`: Predicate to filter channels (default: `channels()` - all channels)
- `interval_selection::Interval`: Time interval as tuple (e.g., (0.0, 1.0)) or interval object for initial data selection (default: nothing - all samples)
- `baseline_interval::Interval`: Baseline interval as tuple (e.g., (-0.2, 0.0)) or interval object (default: nothing - no baseline)
- `analysis_interval::Interval`: Analysis interval as tuple (e.g., (0.3, 0.5)) or interval object for statistical testing (default: nothing - use interval_selection)

# Returns
- `StatisticalData`: Prepared data structure ready for statistical testing
"""
function prepare_stats(
    file_pattern::String,
    design::Symbol;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions([1, 2]),
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = times(),
    analysis_interval::Interval = times(),
)
    # just load all appropriate data and call the main preparation function
    all_erps = read_all_data(ErpData, file_pattern, input_dir, participant_selection)
    isempty(all_erps) && @minimal_error_throw "No valid ERP data found matching pattern '$file_pattern' in $input_dir"

    return prepare_stats(
        all_erps;
        design = design,
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        analysis_interval = analysis_interval,
    )
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
    threshold <= 0.0 || threshold >= 1.0 && @minimal_error "threshold must be between 0 and 1, got $threshold"
    cluster_type ∉ (:spatial, :temporal, :spatiotemporal) &&
        @minimal_error "cluster_type must be :spatial, :temporal, or :spatiotemporal, got :$cluster_type"
    tail ∉ (:both, :left, :right) && @minimal_error "tail must be :both, :left, or :right, got :$tail"

    # Auto-compute neighbours if missing (only when needed for spatial clustering)
    if cluster_type ∈ (:spatial, :spatiotemporal)
        if isnothing(prepared.data[1].layout.neighbours)
            @minimal_warning "Layout.neighbours is not set. Computing with default distance criterion (0.25)."
            # Compute neighbours if missing (using default criterion)
            get_neighbours_xy!(prepared.data[1].layout, 0.25)
            # Also update the other condition's layout if it's a different layout object
            if length(prepared.data) > 1 && prepared.data[2].layout !== prepared.data[1].layout
                # Different layout object, update it too
                get_neighbours_xy!(prepared.data[2].layout, 0.25)
            end
        end
    end

    return nothing
end

# ===================
# T-TEST COMPUTATION
# ===================

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

# Examples
```julia
t_matrix, df, p_matrix = compute_t_matrix(prepared, tail=:both)
```
"""
# Optimized version that accepts raw arrays (for permutation loop - avoids StatisticalData creation)
# Vectorized computation to avoid function call overhead and allocations
# Can accept pre-allocated buffers to avoid allocations
function _compute_t_matrix(
    data1::Array{Float64,3},
    data2::Array{Float64,3},
    design::Symbol;
    tail::Symbol = :both,
    mean1_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    mean2_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    mean_diff_buffer::Union{Nothing,Array{Float64,2}} = nothing,
    std_diff_buffer::Union{Nothing,Array{Float64,2}} = nothing,
)
    n_participants, n_electrodes, n_time = size(data1)
    t_matrix = Array{Float64,2}(undef, n_electrodes, n_time)
    p_matrix = Array{Float64,2}(undef, n_electrodes, n_time)

    if design == :paired

        # Compute mean of differences: mean(data1 - data2) = mean(data1) - mean(data2)
        # Use pre-allocated buffers if provided, otherwise allocate
        if mean1_buffer !== nothing
            mean1 = mean1_buffer
            mean2 = mean2_buffer
            mean_diff = mean_diff_buffer
        else
            mean1 = Array{Float64,2}(undef, n_electrodes, n_time)
            mean2 = Array{Float64,2}(undef, n_electrodes, n_time)
            mean_diff = Array{Float64,2}(undef, n_electrodes, n_time)
        end

        # Compute means and std of differences in one pass (avoids diff array allocation)
        if std_diff_buffer !== nothing
            std_diff = std_diff_buffer
        else
            std_diff = Array{Float64,2}(undef, n_electrodes, n_time)
        end

        @inbounds for e_idx = 1:n_electrodes
            @inbounds for t_idx = 1:n_time
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
                std_diff[e_idx, t_idx] = sqrt(variance * n_participants / (n_participants - 1))

                # Only store mean1/mean2 if buffers were provided
                if mean1_buffer !== nothing
                    mean1[e_idx, t_idx] = mean1_val
                    mean2[e_idx, t_idx] = mean2_val
                end
            end
        end

        # Compute t-statistics: t = mean_diff / (std_diff / sqrt(n))
        # Use in-place assignment to fill pre-allocated t_matrix
        # Handle division by zero
        zero_std_mask = std_diff .== 0.0
        zero_mean_mask = mean_diff .== 0.0

        # Fill pre-allocated t_matrix in-place
        t_matrix .= mean_diff ./ (std_diff ./ sqrt(n_participants))
        # Where std is zero: NaN if mean is also zero, Inf otherwise
        t_matrix[zero_std_mask.&zero_mean_mask] .= NaN
        t_matrix[zero_std_mask.&.!zero_mean_mask] .= Inf

        # Degrees of freedom (same for all points in paired design)
        df = Float64(n_participants - 1)

        # Compute p-values using internal function (avoids code duplication)
        # Use pre-allocated p_matrix buffer
        p_matrix = _compute_p_matrix(t_matrix, df, tail, p_matrix)

    else
        # Independent design: need to loop (but df is constant across all points)
        result = nothing
        @inbounds for e_idx = 1:n_electrodes
            @inbounds for t_idx = 1:n_time
                data_A = view(data1, :, e_idx, t_idx)
                data_B = view(data2, :, e_idx, t_idx)
                result = independent_ttest(data_A, data_B, tail = tail)
                t_matrix[e_idx, t_idx] = result.t
                p_matrix[e_idx, t_idx] = result.p
            end
        end
        df = result.df
    end

    return t_matrix, df, p_matrix
end

function _compute_t_matrix(prepared::StatisticalData; tail::Symbol = :both)
    return _compute_t_matrix(prepared.analysis.data[1], prepared.analysis.data[2], prepared.analysis.design, tail = tail)
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
    if p_matrix_buffer !== nothing
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
        @minimal_error_throw "Statistical tests require exactly 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select exactly 2 conditions."

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
        length(vps1) < 2 || length(vps2) < 2 && @minimal_error "Independent design requires at least 2 participants per group"
    end

    # Validate all TF data have same structure (same channels, same frequencies, same time points)
    _have_same_structure(condition1) || @minimal_error("Condition 1: TF data have inconsistent structure")
    _have_same_structure(condition2) || @minimal_error("Condition 2: TF data have inconsistent structure")
    _have_same_structure(condition1[1], condition2[1]) || @minimal_error("Condition 1 vs. 2: TF data have inconsistent structure")

    # Get metadata from first TF dataset 
    ref_tf = condition1[1]
    all_electrodes = channel_labels(ref_tf)
    all_freqs = sort(unique(ref_tf.data_power.freq))
    all_times = sort(unique(ref_tf.data_power.time))

    # Apply channel selection
    selected_ch_mask = channel_selection(1:length(all_electrodes))
    electrodes = all_electrodes[selected_ch_mask]
    isempty(electrodes) && @minimal_error_throw "No channels matched the channel selection."

    # Apply frequency selection
    if !isnothing(frequency_selection) && frequency_selection isa Tuple
        freq_mask = (all_freqs .>= frequency_selection[1]) .& (all_freqs .<= frequency_selection[2])
        frequencies = all_freqs[freq_mask]
    else
        frequencies = all_freqs
    end
    isempty(frequencies) && @minimal_error_throw "No frequencies matched the frequency selection."

    # Apply time (interval) selection
    if !isnothing(interval_selection) && interval_selection isa Tuple
        time_mask = (all_times .>= interval_selection[1]) .& (all_times .<= interval_selection[2])
        time_points = all_times[time_mask]
    else
        time_points = all_times
    end
    isempty(time_points) && @minimal_error_throw "No time points matched the interval selection."

    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    # Extract 4D arrays: [participants × electrodes × frequencies × time]
    data1 = _extract_tf_array(condition1, electrodes, frequencies, time_points)
    data2 = _extract_tf_array(condition2, electrodes, frequencies, time_points)

    # Create grand averages (average power across participants)
    condition1_avg = _create_tf_grand_average(condition1, electrodes, frequencies, time_points, selected_cond_nums[1])
    condition2_avg = _create_tf_grand_average(condition2, electrodes, frequencies, time_points, selected_cond_nums[2])

    return TFStatisticalData([condition1_avg, condition2_avg], TFAnalysisData(design, [data1, data2], frequencies, time_points))
end

"""
    prepare_stats(::Type{TimeFreqData}, file_pattern, design; input_dir, participant_selection, kwargs...)

Prepare TF data for statistical tests from JLD2 files (convenience wrapper).
"""
function prepare_stats(
    ::Type{TimeFreqData},
    file_pattern::String,
    design::Symbol;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    kwargs...,
)
    all_tfs = read_all_data(TimeFreqData, file_pattern, input_dir, participant_selection)
    isempty(all_tfs) && @minimal_error_throw "No valid TF data found matching pattern '$file_pattern' in $input_dir"

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
        df = tf.data_power

        # Build index vectors for this TF dataset
        tf_freqs = df.freq
        tf_times = df.time

        for (f_idx, freq) in enumerate(frequencies)
            for (t_idx, tp) in enumerate(time_points)
                # Find the row matching this (freq, time) pair
                row_mask = (tf_freqs .== freq) .& (tf_times .== tp)
                row_idx = findfirst(row_mask)

                isnothing(row_idx) && @minimal_error "Missing data point for freq=$freq, time=$tp in participant $(tf.file)"

                for (e_idx, elec) in enumerate(electrodes)
                    data[p_idx, e_idx, f_idx, t_idx] = df[row_idx, elec]
                end
            end
        end
    end

    return data
end


"""
    _create_tf_grand_average(tfs, electrodes, frequencies, time_points, cond_num)

Create a grand average TimeFreqData by averaging power across participants.
"""
function _create_tf_grand_average(
    tfs::Vector{TimeFreqData},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    cond_num::Int,
)
    # Extract 4D data and average over participants (dim 1)
    data_4d = _extract_tf_array(tfs, electrodes, frequencies, time_points)
    avg_power = dropdims(mean(data_4d, dims = 1), dims = 1)  # [electrodes × freqs × time]

    # Build grand average DataFrame with same structure as input
    n_freqs = length(frequencies)
    n_time = length(time_points)
    n_rows = n_freqs * n_time

    # Create time and freq columns (time varies fastest within each freq)
    time_col = repeat(time_points, outer = n_freqs)
    freq_col = repeat(frequencies, inner = n_time)

    # Build DataFrame
    power_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
    for (e_idx, elec) in enumerate(electrodes)
        power_df[!, elec] = vec(avg_power[e_idx, :, :]')  # [freqs × time]' → [time × freqs], vec gives time-fastest
    end

    # Create a matching phase DataFrame (all zeros for grand average since phase averages to zero)
    phase_df = copy(power_df)
    for elec in electrodes
        phase_df[!, elec] .= 0.0
    end

    ref_tf = tfs[1]
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
function _compute_t_matrix_tf(data1::Array{Float64,4}, data2::Array{Float64,4}, design::Symbol; tail::Symbol = :both)
    n_participants, n_electrodes, n_freqs, n_time = size(data1)
    t_matrix = Array{Float64,3}(undef, n_electrodes, n_freqs, n_time)
    p_matrix = Array{Float64,3}(undef, n_electrodes, n_freqs, n_time)

    if design == :paired
        @inbounds for e_idx = 1:n_electrodes
            for f_idx = 1:n_freqs
                for t_idx = 1:n_time
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
        _compute_p_matrix_tf!(p_matrix, t_matrix, df, tail)

    else
        # Independent design
        result = nothing
        @inbounds for e_idx = 1:n_electrodes
            for f_idx = 1:n_freqs
                for t_idx = 1:n_time
                    data_A = view(data1, :, e_idx, f_idx, t_idx)
                    data_B = view(data2, :, e_idx, f_idx, t_idx)
                    result = independent_ttest(data_A, data_B, tail = tail)
                    t_matrix[e_idx, f_idx, t_idx] = result.t
                    p_matrix[e_idx, f_idx, t_idx] = result.p
                end
            end
        end
        df = result.df
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
