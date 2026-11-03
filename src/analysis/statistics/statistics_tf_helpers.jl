# Time-frequency (4D: participants × electrodes × frequencies × time) statistical preparation and helpers.

function _filter_axis_interval(values::Vector{Float64}, sel::Interval, empty_msg::String)
    (isnothing(sel) || sel isa AllSelection) && return values, nothing
    lo, hi = sel isa TimeSelection ? (sel.start, sel.stop) : (sel[1], sel[2])
    mask = (values .>= lo) .& (values .<= hi)
    filtered = values[mask]
    isempty(filtered) && @minimal_error(empty_msg)
    return filtered, mask
end

function _baseline_tf_4d!(data::Array{Float64,4}, base_range::UnitRange{Int}, method::Symbol)
    n_p, n_e, _, _ = size(data)
    _foreach_chunk(n_p * n_e) do idx_range
        for idx in idx_range
            p = mod1(idx, n_p)
            e = cld(idx, n_p)
            _tf_baseline_grid!(@view(data[p, e, :, :]), base_range, method)
        end
    end
    return data
end

"""
    prepare_stats(tfs::AbstractVector{<:TimeFreqData}; design = :paired, condition_selection = conditions(), channel_selection = channels(), frequency_selection = nothing, interval_selection = nothing, baseline_interval = nothing, baseline_method = :db, analysis_interval = nothing)

Prepare `TimeFreqData` for comparing two conditions in statistical tests.

Organizes `TimeFreqData` power values into `participant × electrode × frequency × time` arrays,
computes grand averages over the display interval, and optionally applies baseline normalization and an `analysis_interval`.

# Returns
- `TFStatisticalData`: Prepared data structure ready for statistical testing
"""
function prepare_stats(
    tfs::AbstractVector{<:TimeFreqData};
    design::Symbol = :paired,
    condition_selection = conditions(),
    channel_selection = channels(),
    frequency_selection::Interval = nothing,
    interval_selection::Interval = nothing,
    baseline_interval::Interval = nothing,
    baseline_method::Symbol = :db,
    analysis_interval::Interval = nothing,
)
    condition1, condition2, selected_cond_nums = _select_two_conditions(tfs, design, condition_selection, "TF data")

    ref_tf = condition1[1]
    all_electrodes = channel_labels(ref_tf)
    all_freqs = Float64.(sort(unique(ref_tf.data_power.freq)))
    all_times = Float64.(sort(unique(ref_tf.data_power.time)))

    electrodes = all_electrodes[_to_channel_func(channel_selection)(all_electrodes)]
    isempty(electrodes) && @minimal_error "No channels matched the channel selection."

    frequencies, _ = _filter_axis_interval(all_freqs, frequency_selection, "No frequencies matched the frequency selection.")
    display_times, time_mask = _filter_axis_interval(all_times, interval_selection, "No time points matched the interval selection.")
    time_points, analysis_mask =
        _filter_axis_interval(display_times, analysis_interval, "No time points matched the analysis interval criteria!")

    if !isnothing(baseline_interval)
        for tf in condition1
            _check_tf_not_baselined(tf.baseline)
        end
        for tf in condition2
            _check_tf_not_baselined(tf.baseline)
        end
        _validate_tf_baseline_method(baseline_method)
        bl_tuple = baseline_interval isa AllSelection ? _full_tf_interval(ref_tf) : _to_tf_interval(baseline_interval)
        b_start, b_end = _validate_and_resolve_tf_baseline_times(all_times, bl_tuple)
        base_range = b_start:b_end

        full_data1 = _baseline_tf_4d!(_extract_tf_array(condition1, electrodes, frequencies, all_times), base_range, baseline_method)
        full_data2 = _baseline_tf_4d!(_extract_tf_array(condition2, electrodes, frequencies, all_times), base_range, baseline_method)

        display_data1 = isnothing(time_mask) ? full_data1 : full_data1[:, :, :, time_mask]
        display_data2 = isnothing(time_mask) ? full_data2 : full_data2[:, :, :, time_mask]
        bl_info = BaselineInfo(method = baseline_method, window = bl_tuple)
    else
        display_data1 = _extract_tf_array(condition1, electrodes, frequencies, display_times)
        display_data2 = _extract_tf_array(condition2, electrodes, frequencies, display_times)
        bl_info = ref_tf.baseline
    end

    condition1_avg = _create_tf_grand_average(
        condition1[1],
        display_data1,
        electrodes,
        frequencies,
        display_times,
        selected_cond_nums[1];
        baseline = bl_info,
    )
    condition2_avg = _create_tf_grand_average(
        condition2[1],
        display_data2,
        electrodes,
        frequencies,
        display_times,
        selected_cond_nums[2];
        baseline = bl_info,
    )

    data1 = isnothing(analysis_mask) ? display_data1 : display_data1[:, :, :, analysis_mask]
    data2 = isnothing(analysis_mask) ? display_data2 : display_data2[:, :, :, analysis_mask]

    return TFStatisticalData([condition1_avg, condition2_avg], TFAnalysisData(design, [data1, data2], frequencies, time_points))
end

function prepare_stats(
    ::Type{TimeFreqData},
    file_pattern::String,
    design::Symbol;
    input_dir::String = pwd(),
    participant_selection = participants(),
    kwargs...,
)
    all_tfs = read_all_data(TimeFreqData, joinpath(input_dir, file_pattern), participant_selection)
    isempty(all_tfs) && @minimal_error "No valid TF data found matching pattern '$file_pattern' in $input_dir"
    return prepare_stats(all_tfs; design = design, kwargs...)
end

function _extract_tf_array(
    tfs::AbstractVector{<:TimeFreqData},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
)
    n_participants = length(tfs)
    n_electrodes = length(electrodes)
    n_freqs = length(frequencies)
    n_time = length(time_points)

    first_df = tfs[1].data_power
    if n_freqs > 0 &&
       n_time > 0 &&
       nrow(first_df) == n_freqs * n_time &&
       first_df.time[1] == first_df.time[min(n_freqs, end)] &&
       isapprox(first_df.time[1], time_points[1]; atol = 1e-6) &&
       isapprox(first_df.time[end], time_points[end]; atol = 1e-6) &&
       isapprox(first_df.freq[1], frequencies[1]; atol = 1e-6) &&
       isapprox(first_df.freq[end], frequencies[end]; atol = 1e-6) &&
       all(tf -> nrow(tf.data_power) == n_freqs * n_time, tfs)
        data = Array{Float64,4}(undef, n_participants, n_electrodes, n_freqs, n_time)
        for (p_idx, tf) in enumerate(tfs), (e_idx, elec) in enumerate(electrodes)
            copyto!(@view(data[p_idx, e_idx, :, :]), tf.data_power[!, elec]::Vector{Float64})
        end
        return data
    end

    all_freqs = unique(first_df.freq)
    all_times = unique(first_df.time)
    n_all_f = length(all_freqs)
    n_all_t = length(all_times)

    is_regular =
        n_all_f > 0 &&
        n_all_t > 0 &&
        nrow(first_df) == n_all_f * n_all_t &&
        first_df.time[1] == first_df.time[min(n_all_f, end)] &&
        all(
            tf -> nrow(tf.data_power) == n_all_f * n_all_t && tf.data_power.time == first_df.time && tf.data_power.freq == first_df.freq,
            tfs,
        )

    if is_regular
        rf_map = Dict{Float64,Int}(round(f, digits = 6) => i for (i, f) in enumerate(all_freqs))
        rt_map = Dict{Float64,Int}(round(t, digits = 6) => i for (i, t) in enumerate(all_times))
        f_idxs = Int[get(rf_map, round(f, digits = 6), 0) for f in frequencies]
        t_idxs = Int[get(rt_map, round(t, digits = 6), 0) for t in time_points]

        if all(>(0), f_idxs) && all(>(0), t_idxs)
            data = Array{Float64,4}(undef, n_participants, n_electrodes, n_freqs, n_time)
            for (p_idx, tf) in enumerate(tfs), (e_idx, elec) in enumerate(electrodes)
                src = reshape(tf.data_power[!, elec]::Vector{Float64}, n_all_f, n_all_t)
                @inbounds for (ti, st) in enumerate(t_idxs), (fi, sf) in enumerate(f_idxs)
                    data[p_idx, e_idx, fi, ti] = src[sf, st]
                end
            end
            return data
        end
    end

    data = fill(NaN, n_participants, n_electrodes, n_freqs, n_time)
    rounded_freqs = round.(frequencies, digits = 6)
    rounded_times = round.(time_points, digits = 6)

    for (p_idx, tf) in enumerate(tfs)
        row_index = _tf_build_row_index(tf.data_power)
        for (e_idx, elec) in enumerate(electrodes)
            col_data = tf.data_power[!, elec]::Vector{Float64}
            @inbounds for (ti, rt) in enumerate(rounded_times), (fi, rf) in enumerate(rounded_freqs)
                row = get(row_index, (rf, rt), 0)
                if row > 0
                    data[p_idx, e_idx, fi, ti] = col_data[row]
                end
            end
        end
    end
    return data
end

function _create_tf_grand_average(
    ref_tf::TimeFreqData,
    data_4d::Array{Float64,4},
    electrodes::Vector{Symbol},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    cond_num::Int;
    baseline::Union{BaselineInfo,Nothing} = ref_tf.baseline,
)
    n_p, _, n_freqs, n_time = size(data_4d)
    n_rows = n_freqs * n_time
    inv_n = 1.0 / n_p

    time_col = repeat(time_points, inner = n_freqs)
    freq_col = repeat(frequencies, outer = n_time)

    power_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
    phase_df = DataFrame(time = copy(time_col), freq = copy(freq_col), copycols = false)

    for (e_idx, elec) in enumerate(electrodes)
        col = Vector{Float64}(undef, n_rows)
        idx = 0
        @inbounds for t = 1:n_time, f = 1:n_freqs
            s = 0.0
            @simd for p = 1:n_p
                s += data_4d[p, e_idx, f, t]
            end
            idx += 1
            col[idx] = s * inv_n
        end
        power_df[!, elec] = col
        phase_df[!, elec] = fill(NaN, n_rows)
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
        baseline,
        copy(ref_tf.analysis_info),
    )
end

function _compute_t_matrix_tf(
    data1::Array{Float64,4},
    data2::Array{Float64,4},
    design::Symbol;
    tail::Symbol = :both,
    compute_p_values::Bool = true,
    t_matrix_buffer::Union{Nothing,Array{Float64,3}} = nothing,
)
    _, n_electrodes, n_freqs, n_time = size(data1)
    n_voxels = n_electrodes * n_freqs * n_time

    t_matrix = isnothing(t_matrix_buffer) ? Array{Float64,3}(undef, n_electrodes, n_freqs, n_time) : t_matrix_buffer
    df = _t_stat_kernel!(vec(t_matrix), nothing, reshape(data1, size(data1, 1), n_voxels), reshape(data2, size(data2, 1), n_voxels), design)

    p_matrix = compute_p_values ? _compute_p_matrix(t_matrix, df, tail) : Array{Float64,3}(undef, 0, 0, 0)
    return t_matrix, df, p_matrix
end

function _compute_t_matrix_tf(prepared::TFStatisticalData; tail::Symbol = :both, compute_p_values::Bool = true)
    return _compute_t_matrix_tf(
        prepared.analysis.data[1],
        prepared.analysis.data[2],
        prepared.analysis.design;
        tail = tail,
        compute_p_values = compute_p_values,
    )
end

_compute_p_matrix_tf!(p_matrix::Array{Float64,3}, t_matrix::Array{Float64,3}, df::Float64, tail::Symbol) =
    _compute_p_matrix(t_matrix, df, tail, p_matrix)

const _compute_critical_t_values_tf = _compute_critical_t_values
