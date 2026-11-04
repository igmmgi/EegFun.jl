# Core data preparation and t-statistic kernels for EEG/ERP statistical testing.

function _select_two_conditions(items::AbstractVector, design::Symbol, condition_selection, label::String)
    design ∈ (:paired, :independent) || @minimal_error "design must be :paired or :independent, got :$design"

    by_cond = group_by_condition(items)
    all_cond_nums = collect(keys(by_cond))
    cond_func = _to_condition_func(condition_selection)
    selected_cond_nums = all_cond_nums[cond_func(1:length(all_cond_nums))]

    length(selected_cond_nums) == 2 ||
        @minimal_error "Statistical tests require exactly 2 conditions, got $(length(selected_cond_nums)): $selected_cond_nums. Use condition_selection to select exactly 2 conditions."

    cond1 = by_cond[selected_cond_nums[1]]
    cond2 = by_cond[selected_cond_nums[2]]

    vps1 = [_extract_participant_id(basename(d.file)) for d in cond1]
    vps2 = [_extract_participant_id(basename(d.file)) for d in cond2]

    if design == :paired
        vps1 == vps2 || @minimal_error "Paired design requires same participants in both conditions"
    else
        (length(vps1) >= 2 && length(vps2) >= 2) || @minimal_error "Independent design requires at least 2 participants per group"
    end

    _have_same_structure(cond1) || @minimal_error("Condition 1: $label have inconsistent structure")
    _have_same_structure(cond2) || @minimal_error("Condition 2: $label have inconsistent structure")
    _have_same_structure(cond1[1], cond2[1]) || @minimal_error("Condition 1 vs. 2: $label have inconsistent structure")

    return cond1, cond2, selected_cond_nums
end

@inline function _fill_baselined_slice!(
    data::Array{Float64,3},
    p::Int,
    e::Int,
    col::Vector{Float64},
    rows::Vector{Int},
    bl_start::Int,
    bl_stop::Int,
)
    bl_sum = 0.0
    @inbounds @simd for idx = bl_start:bl_stop
        bl_sum += col[rows[idx]]
    end
    bl_mean = bl_sum / (bl_stop - bl_start + 1)
    @inbounds @simd for t in eachindex(rows)
        data[p, e, t] = col[rows[t]] - bl_mean
    end
end

@inline function _fill_baselined_slice!(
    data::Array{Float64,3},
    p::Int,
    e::Int,
    col::AbstractVector,
    rows::Vector{Int},
    bl_start::Int,
    bl_stop::Int,
)
    bl_sum = 0.0
    @inbounds for idx = bl_start:bl_stop
        bl_sum += Float64(col[rows[idx]])
    end
    bl_mean = bl_sum / (bl_stop - bl_start + 1)
    @inbounds for t in eachindex(rows)
        data[p, e, t] = Float64(col[rows[t]]) - bl_mean
    end
end

function _extract_baselined_erp_array(
    erps::AbstractVector{<:ErpData},
    electrodes::Vector{Symbol},
    rows::Vector{Int},
    bl_range::Tuple{Int,Int},
)
    data = Array{Float64,3}(undef, length(erps), length(electrodes), length(rows))
    bl_start, bl_stop = bl_range
    for (p, erp) in enumerate(erps)
        df = erp.data
        for (e, el) in enumerate(electrodes)
            _fill_baselined_slice!(data, p, e, df[!, el], rows, bl_start, bl_stop)
        end
    end
    return data
end

function _grand_avg_from_stacked!(ref::ErpData, stacked::Array{Float64,3}, electrodes::Vector{Symbol}, cond_num::Int, total_epochs::Int)
    n_p, _, n_t = size(stacked)
    inv_n = 1.0 / n_p
    df = ref.data
    for (e, el) in enumerate(electrodes)
        col = Vector{Float64}(undef, n_t)
        @inbounds for t = 1:n_t
            s = 0.0
            @simd for p = 1:n_p
                s += stacked[p, e, t]
            end
            col[t] = s * inv_n
        end
        df[!, el] = col
    end
    return ErpData(
        "grand_avg",
        cond_num,
        "grand_avg_$(ref.condition_name)",
        df,
        ref.layout,
        ref.sample_rate,
        ref.analysis_info,
        total_epochs,
    )
end

function _compute_display_se(cond1::Array{Float64,3}, cond2::Array{Float64,3}, design::Symbol)
    n1, n_e, n_t = size(cond1)
    n2 = size(cond2, 1)
    se1 = Matrix{Float64}(undef, n_e, n_t)
    se2 = Matrix{Float64}(undef, n_e, n_t)
    se_diff = Matrix{Float64}(undef, n_e, n_t)
    inv_n1_var = 1.0 / ((n1 - 1) * n1)
    inv_n2_var = 1.0 / ((n2 - 1) * n2)

    @inbounds for t = 1:n_t, e = 1:n_e
        s1 = 0.0
        @simd for p = 1:n1
            s1 += cond1[p, e, t]
        end
        m1 = s1 / n1
        ss1 = 0.0
        @simd for p = 1:n1
            d = cond1[p, e, t] - m1
            ss1 += d * d
        end
        v1 = max(ss1 * inv_n1_var, 0.0)
        se1[e, t] = sqrt(v1)

        s2 = 0.0
        @simd for p = 1:n2
            s2 += cond2[p, e, t]
        end
        m2 = s2 / n2
        ss2 = 0.0
        @simd for p = 1:n2
            d = cond2[p, e, t] - m2
            ss2 += d * d
        end
        v2 = max(ss2 * inv_n2_var, 0.0)
        se2[e, t] = sqrt(v2)

        if design == :paired
            md = m1 - m2
            ssd = 0.0
            @simd for p = 1:n1
                dd = (cond1[p, e, t] - cond2[p, e, t]) - md
                ssd += dd * dd
            end
            se_diff[e, t] = sqrt(max(ssd * inv_n1_var, 0.0))
        else
            se_diff[e, t] = sqrt(v1 + v2)
        end
    end
    return se1, se2, se_diff
end

"""
    prepare_stats(erps::AbstractVector{<:ErpData}; design::Symbol = :paired, condition_selection = conditions(), channel_selection = channels(), interval_selection::Interval = times(), baseline_interval::Interval = nothing, analysis_interval::Interval = times())

Prepare `ErpData` for comparing two conditions in statistical tests (permutation and analytic tests).

Organizes `ErpData` into `participant × electrode × time` arrays, computes grand averages and
standard errors over the display interval, and optionally narrows the statistical testing window via `analysis_interval`.

# Returns
- `StatisticalData`: Prepared data structure ready for statistical testing
"""
function prepare_stats(
    erps::AbstractVector{<:ErpData};
    design::Symbol = :paired,
    condition_selection = conditions(),
    channel_selection = channels(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    analysis_interval::Interval = times(),
)
    condition1, condition2, selected_cond_nums = _select_two_conditions(erps, design, condition_selection, "ERPs")

    sample_sel = _interval_to_samples(interval_selection)
    ref1 = subset(condition1[1]; channel_selection = channel_selection, sample_selection = sample_sel)
    ref2 = subset(condition2[1]; channel_selection = channel_selection, sample_selection = sample_sel)

    display_electrodes = channel_labels(ref1)
    display_rows = findall(_to_sample_func(sample_sel)(condition1[1].data))
    (isempty(display_electrodes) || isempty(display_rows)) && @minimal_error "No data matched the selection criteria!"

    bl_range = _validate_baseline_interval(ref1, isnothing(baseline_interval) ? times() : baseline_interval)

    cond1_stacked = _extract_baselined_erp_array(condition1, display_electrodes, display_rows, bl_range)
    cond2_stacked = _extract_baselined_erp_array(condition2, display_electrodes, display_rows, bl_range)

    epochs1 = sum(erp.n_epochs for erp in condition1)
    epochs2 = sum(erp.n_epochs for erp in condition2)
    condition1_avg = _grand_avg_from_stacked!(ref1, cond1_stacked, display_electrodes, selected_cond_nums[1], epochs1)
    condition2_avg = _grand_avg_from_stacked!(ref2, cond2_stacked, display_electrodes, selected_cond_nums[2], epochs2)

    se_cond1, se_cond2, se_diff = _compute_display_se(cond1_stacked, cond2_stacked, design)

    if isnothing(analysis_interval) || analysis_interval isa AllSelection
        time_points = copy(condition1_avg.data[!, :time])
        condition1_array = cond1_stacked
        condition2_array = cond2_stacked
    else
        analysis_mask = _to_sample_func(_interval_to_samples(analysis_interval))(condition1_avg.data)
        any(analysis_mask) || @minimal_error "No data matched the analysis interval criteria!"
        time_points = condition1_avg.data[analysis_mask, :time]
        condition1_array = cond1_stacked[:, :, analysis_mask]
        condition2_array = cond2_stacked[:, :, analysis_mask]
    end

    return StatisticalData(
        [condition1_avg, condition2_avg],
        AnalysisData(design, [condition1_array, condition2_array], time_points),
        se_cond1,
        se_cond2,
        se_diff,
    )
end

"""Load ERP/TF files by pattern and call the typed `prepare_stats` method."""
function prepare_stats(file_pattern::String, design::Symbol; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
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

function _validate_permutation_inputs(
    prepared::StatisticalData,
    n_permutations::Int,
    alpha::Float64,
    cluster_type::Symbol,
    tail::Symbol;
    cluster_alpha::Float64 = alpha,
)
    n_permutations < 1 && @minimal_error "n_permutations must be >= 1, got $n_permutations"
    (alpha <= 0.0 || alpha >= 1.0) && @minimal_error "alpha must be between 0 and 1, got $alpha"
    (cluster_alpha <= 0.0 || cluster_alpha >= 1.0) && @minimal_error "cluster_alpha must be between 0 and 1, got $cluster_alpha"
    cluster_type ∉ (:spatial, :temporal, :spatiotemporal) &&
        @minimal_error "cluster_type must be :spatial, :temporal, or :spatiotemporal, got :$cluster_type"
    tail ∉ (:both, :left, :right) && @minimal_error "tail must be :both, :left, or :right, got :$tail"

    if cluster_type ∈ (:spatial, :spatiotemporal) && isnothing(prepared.data[1].layout.neighbors)
        @minimal_warning "Layout.neighbors is not set. Computing with default distance criterion (0.35)."
    end
    return nothing
end

# === UNIFIED N-D T-STATISTIC KERNELS ===

function _t_stat_kernel!(
    t_flat::AbstractVector{Float64},
    se_flat::Union{Nothing,AbstractVector{Float64}},
    d1_flat::AbstractMatrix{Float64},
    d2_flat::AbstractMatrix{Float64},
    design::Symbol,
)
    n_voxels = size(d1_flat, 2)
    store_se = !isnothing(se_flat)

    if design == :paired
        n_p = size(d1_flat, 1)
        df = Float64(n_p - 1)
        inv_n = 1.0 / n_p
        inv_var_se = 1.0 / (df * n_p)

        @inbounds for v = 1:n_voxels
            sum_diff = 0.0
            @simd for p = 1:n_p
                sum_diff += d1_flat[p, v] - d2_flat[p, v]
            end
            mean_diff = sum_diff * inv_n

            ss_diff = 0.0
            @simd for p = 1:n_p
                d = (d1_flat[p, v] - d2_flat[p, v]) - mean_diff
                ss_diff += d * d
            end

            se = sqrt(max(ss_diff * inv_var_se, 0.0))
            t_flat[v] = se == 0.0 ? (mean_diff == 0.0 ? NaN : copysign(Inf, mean_diff)) : (mean_diff / se)
            if store_se
                se_flat[v] = se
            end
        end
        return df
    else
        n_A = size(d1_flat, 1)
        n_B = size(d2_flat, 1)
        df = Float64(n_A + n_B - 2)
        inv_A = 1.0 / n_A
        inv_B = 1.0 / n_B
        se_scale = (inv_A + inv_B) / df

        @inbounds for v = 1:n_voxels
            sum1 = 0.0
            @simd for p = 1:n_A
                sum1 += d1_flat[p, v]
            end
            sum2 = 0.0
            @simd for p = 1:n_B
                sum2 += d2_flat[p, v]
            end
            mean1 = sum1 * inv_A
            mean2 = sum2 * inv_B

            ss1 = 0.0
            @simd for p = 1:n_A
                d1 = d1_flat[p, v] - mean1
                ss1 += d1 * d1
            end
            ss2 = 0.0
            @simd for p = 1:n_B
                d2 = d2_flat[p, v] - mean2
                ss2 += d2 * d2
            end

            se = sqrt(max((ss1 + ss2) * se_scale, 0.0))
            diff_mean = mean1 - mean2
            t_flat[v] = se == 0.0 ? (diff_mean == 0.0 ? NaN : copysign(Inf, diff_mean)) : (diff_mean / se)
            if store_se
                se_flat[v] = se
            end
        end
        return df
    end
end

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
    compute_se::Bool = true,
)
    _, n_electrodes, n_time = size(data1)
    n_voxels = n_electrodes * n_time

    t_matrix = isnothing(t_matrix_buffer) ? Matrix{Float64}(undef, n_electrodes, n_time) : t_matrix_buffer
    se_matrix = compute_se ? Matrix{Float64}(undef, n_electrodes, n_time) : Matrix{Float64}(undef, 0, 0)

    df = _t_stat_kernel!(
        vec(t_matrix),
        compute_se ? vec(se_matrix) : nothing,
        reshape(data1, size(data1, 1), n_voxels),
        reshape(data2, size(data2, 1), n_voxels),
        design,
    )

    p_matrix = compute_p_values ? _compute_p_matrix(t_matrix, df, tail) : Matrix{Float64}(undef, 0, 0)
    return t_matrix, df, p_matrix, se_matrix
end

function _compute_t_matrix(prepared::StatisticalData; tail::Symbol = :both, compute_p_values::Bool = true, compute_se::Bool = true)
    return _compute_t_matrix(
        prepared.analysis.data[1],
        prepared.analysis.data[2],
        prepared.analysis.design;
        tail = tail,
        compute_p_values = compute_p_values,
        compute_se = compute_se,
    )
end

function _critical_t_scalar(df::Real, alpha::Real, tail::Symbol)
    dist = TDist(df)
    tail == :both && return quantile(dist, 1.0 - alpha / 2.0)
    tail == :right && return quantile(dist, 1.0 - alpha)
    tail == :left && return quantile(dist, alpha)
    @minimal_error("tail must be :both, :left, or :right, got :$tail")
end

function _compute_critical_t_values(df::Real, matrix_size::NTuple{N,Int}, alpha::Real = 0.05, tail::Symbol = :both) where {N}
    critical_t_values = Array{Float64,N}(undef, matrix_size)
    if !isfinite(df) || df <= 0
        fill!(critical_t_values, NaN)
        return critical_t_values
    end
    fill!(critical_t_values, _critical_t_scalar(df, alpha, tail))
    return critical_t_values
end

function _compute_p_matrix(
    t_matrix::Array{Float64,N},
    df::Float64,
    tail::Symbol,
    p_matrix_buffer::Union{Nothing,Array{Float64,N}} = nothing,
) where {N}
    p_matrix = isnothing(p_matrix_buffer) ? Array{Float64,N}(undef, size(t_matrix)) : p_matrix_buffer
    dist = TDist(df)
    if tail == :both
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            p_matrix[i] = isnan(t_val) ? NaN : 2.0 * ccdf(dist, abs(t_val))
        end
    elseif tail == :left
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            p_matrix[i] = isnan(t_val) ? NaN : cdf(dist, t_val)
        end
    else
        @inbounds for i in eachindex(t_matrix)
            t_val = t_matrix[i]
            p_matrix[i] = isnan(t_val) ? NaN : ccdf(dist, t_val)
        end
    end
    return p_matrix
end
