# Decibel, percent-change, relative, and z-score baseline correction for TF data.

######################################
######### Types & Constants ##########
######################################

const _VALID_TF_BASELINE_METHODS = (:absolute, :relative, :relchange, :normchange, :db, :zscore, :percent)

const _TFDataObject = Union{TimeFreqData,TimeFreqEpochData}

######################################
########### User Functions ###########
######################################

"""
    tf_baseline(tf_data::Union{TimeFreqData,TimeFreqEpochData}, baseline_interval = times(); method = nothing)
    tf_baseline(tf_data::AbstractVector{<:Union{TimeFreqData,TimeFreqEpochData}}, baseline_interval = times(); method = nothing)
    tf_baseline(file_pattern::String, baseline_interval; method = nothing, input_dir = pwd(), output_dir = nothing,
                participant_selection = participants(), condition_selection = conditions())

Apply baseline correction to `TimeFreqData` or single-trial `TimeFreqEpochData` without modifying `tf_data`
(or batch-process matching `.jld2` files when given `file_pattern::String`).

Defaults to `:db` for power decompositions and `:absolute` (subtraction) for ITC data.
Supported `method` values: `:db` (`10 * log10(P / B)`), `:absolute` (`P - B`),
`:relative` (`P / B`), `:relchange` (`(P - B) / B`), `:percent` (`100 * (P - B) / B`),
`:normchange` (`(P - B) / (P + B)`), and `:zscore` (`(P - B) / SD_B`).

# Examples
```julia
tf_corr = tf_baseline(tf_data, (-0.3, 0.0); method = :db)
tf_baseline("tf_epochs_wavelet", (-0.3, 0.0); method = :db)
```
"""
tf_baseline(tf_data::_TFDataObject, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    tf_baseline!(copy(tf_data), baseline_interval; method = method)

tf_baseline(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    [tf_baseline(tf, baseline_interval; method = method) for tf in tf_data]

"""
    tf_baseline!(tf_data::TimeFreqData, baseline_interval = times(); method = nothing)
    tf_baseline!(tf_data::TimeFreqEpochData, baseline_interval = times(); method = nothing)

In-place mutating version of [`tf_baseline`](@ref).
"""
function tf_baseline!(tf_data::TimeFreqData, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing)
    win = _resolve_tf_baseline_window(tf_data, baseline_interval)
    isnothing(win) && return tf_data
    _check_tf_not_baselined(tf_data.baseline)
    eff_method = isnothing(method) ? (_is_itc(tf_data) ? :absolute : :db) : method
    _validate_tf_baseline_method(eff_method)
    if _is_itc(tf_data) && eff_method ∉ (:absolute, :zscore)
        @minimal_warning "Baseline method :$eff_method is a power ratio and is unusual for ITC values; consider method = :absolute (subtraction)."
    end
    ch_labels = channel_labels(tf_data)
    @info "Applying $(eff_method) baseline correction ($(win[1])s to $(win[2])s)"
    _tf_baseline_df!(tf_data.data_power, ch_labels, win, eff_method; parallel_channels = true)
    tf_data.baseline = BaselineInfo(method = eff_method, window = win)
    return tf_data
end

function tf_baseline!(tf_data::TimeFreqEpochData, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing)
    win = _resolve_tf_baseline_window(tf_data, baseline_interval)
    isnothing(win) && return tf_data
    _check_tf_not_baselined(tf_data.baseline)
    eff_method = isnothing(method) ? (_is_itc(tf_data) ? :absolute : :db) : method
    _validate_tf_baseline_method(eff_method)
    isempty(tf_data.data_power) && @minimal_error("Cannot baseline-correct empty TimeFreqEpochData")

    df1 = tf_data.data_power[1]
    times_vec = unique(df1.time)
    freqs_unique = unique(df1.freq)
    n_freqs, n_times = length(freqs_unique), length(times_vec)
    b1, b2 = _validate_and_resolve_tf_baseline_times(times_vec, win)
    base_range = b1:b2
    n_grid = n_freqs * n_times
    regular_ref = (nrow(df1) == n_grid) && (n_grid > 0) && (df1.time[1] == df1.time[min(n_freqs, end)])

    ch_labels = channel_labels(tf_data)
    @info "Applying $(eff_method) baseline correction to $(length(tf_data.data_power)) TF epochs ($(win[1])s to $(win[2])s)"

    _foreach_chunk(length(tf_data.data_power)) do ep_range
        bm_buf = Vector{Float64}(undef, n_freqs)
        cnt_buf = Vector{Int}(undef, n_freqs)
        work_buf = Vector{Float64}(undef, n_freqs)
        for ep_idx in ep_range
            df = tf_data.data_power[ep_idx]
            if regular_ref && nrow(df) == n_grid && df.time[1] == times_vec[1] && df.time[end] == times_vec[end]
                for col in _ensure_f64_channel_cols!(df, ch_labels)
                    _tf_baseline_grid!(reshape(col, n_freqs, n_times), base_range, eff_method, bm_buf, cnt_buf, work_buf)
                end
            else
                _tf_baseline_df!(df, ch_labels, win, eff_method; parallel_channels = false)
            end
        end
    end
    tf_data.baseline = BaselineInfo(method = eff_method, window = win)
    return tf_data
end

tf_baseline!(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    (foreach(tf -> tf_baseline!(tf, baseline_interval; method = method), tf_data); tf_data)

function tf_baseline(
    file_pattern::String,
    baseline_interval::Interval;
    method::Union{Symbol,Nothing} = nothing,
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection = participants(),
    condition_selection = conditions(),
)
    _with_batch_logging("TF baseline", "tf_baseline.log", input_dir, output_dir) do ctx
        _ = output_dir
        @log_call "tf_baseline"
        @info "baseline: $baseline_interval, method: $(something(method, :auto))"
        default_dir = joinpath(input_dir, "tf_baseline_$(_clean_pattern(file_pattern))")
        return _batch_jld2_files(
            data -> (tf_baseline!(data, baseline_interval; method = method), "Baseline corrected"),
            Vector{<:_TFDataObject},
            ctx,
            default_dir,
            file_pattern,
            input_dir,
            participant_selection,
            "TF baseline";
            condition_selection = condition_selection,
        )
    end
end

"""
    baseline(tf_data::Union{TimeFreqData,TimeFreqEpochData}, baseline_interval = times(); method = nothing)

Alias for [`tf_baseline`](@ref) on time-frequency containers.
"""
baseline(tf_data::_TFDataObject, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    tf_baseline(tf_data, baseline_interval; method = method)

baseline(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    tf_baseline(tf_data, baseline_interval; method = method)

baseline!(tf_data::_TFDataObject, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    tf_baseline!(tf_data, baseline_interval; method = method)

baseline!(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::Interval = times(); method::Union{Symbol,Nothing} = nothing) =
    tf_baseline!(tf_data, baseline_interval; method = method)

######################################
######### Internal Functions #########
######################################

function _validate_tf_baseline_method(method::Symbol)
    method ∈ _VALID_TF_BASELINE_METHODS ||
        @minimal_error("Unknown baseline method: $method. Use :absolute, :relative, :relchange, :normchange, :db, :zscore, or :percent")
    return nothing
end

function _check_tf_not_baselined(baseline_info::Union{BaselineInfo,Nothing})
    if !isnothing(baseline_info)
        @minimal_error(
            "Baseline correction has already been applied to this data (method: $(baseline_info.method), window: $(baseline_info.window)). " *
            "Baseline corrections are non-linear and cannot be chained. Use the original data to apply a different baseline.",
        )
    end
    return nothing
end

function _resolve_tf_baseline_window(tf_data::_TFDataObject, interval::Interval)::Union{Nothing,Tuple{Float64,Float64}}
    isnothing(interval) && return nothing
    interval isa Tuple{Real,Real} && return (Float64(interval[1]), Float64(interval[2]))
    interval isa TimeSelection && return (interval.start, interval.stop)
    if interval isa AbstractVector{<:Real}
        length(interval) == 2 || @minimal_error("Baseline interval vector must have 2 elements [start, stop], got: $interval")
        return (Float64(interval[1]), Float64(interval[2]))
    end
    if tf_data isa TimeFreqEpochData && isempty(tf_data.data_power)
        @minimal_error(
            interval isa AllSelection ? "TimeFreqEpochData contains no time points" : "Cannot baseline-correct empty TimeFreqEpochData"
        )
    end
    df = tf_data isa TimeFreqData ? tf_data.data_power : tf_data.data_power[1]
    if interval isa AllSelection
        isempty(df.time) && @minimal_error(
            tf_data isa TimeFreqData ? "Time-frequency data contains no time points" : "TimeFreqEpochData contains no time points"
        )
        return Float64.(extrema(df.time))
    end
    idxs = findall(interval(df))
    isempty(idxs) && @minimal_error("Baseline window selection returned no samples!")
    return Float64.(extrema(@view df.time[idxs]))
end

function _validate_and_resolve_tf_baseline_times(times::Vector{Float64}, baseline_interval::Tuple{Real,Real})
    b_start, b_stop = Float64(baseline_interval[1]), Float64(baseline_interval[2])
    (isnan(b_start) || isnan(b_stop) || b_start > b_stop) &&
        @minimal_error("Invalid baseline interval $(baseline_interval): start must be <= stop")
    isempty(times) && @minimal_error("Time-frequency data contains no time points")

    t_min, t_max = extrema(times)
    (b_stop < t_min || b_start > t_max) &&
        @minimal_error("Baseline window $(baseline_interval) does not overlap with data times ($(t_min) to $(t_max) s)")

    i1 = _find_closest_time_index(times, b_start)
    i2 = _find_closest_time_index(times, b_stop)
    return minmax(i1, i2)
end

function _ensure_f64_channel_cols!(df::DataFrame, ch_labels::Vector{Symbol})::Vector{Vector{Float64}}
    n_ch = length(ch_labels)
    cols = Vector{Vector{Float64}}(undef, n_ch)
    for (ci, ch) in enumerate(ch_labels)
        col = df[!, ch]
        if col isa Vector{Float64}
            cols[ci] = col
        else
            col_f64 = Float64.(col)
            df[!, ch] = col_f64
            cols[ci] = col_f64
        end
    end
    return cols
end

function _tf_baseline_df!(
    df::DataFrame,
    ch_labels::Vector{Symbol},
    baseline_interval::Tuple{Real,Real},
    method::Symbol;
    parallel_channels::Bool = true,
)
    times = unique(df.time)
    freqs_unique = unique(df.freq)
    n_freqs, n_times, n_channels = length(freqs_unique), length(times), length(ch_labels)

    b1, b2 = _validate_and_resolve_tf_baseline_times(times, baseline_interval)
    base_range = b1:b2
    is_regular = (nrow(df) == n_freqs * n_times) && (n_freqs > 0 && n_times > 0) && (df.time[1] == df.time[min(n_freqs, end)])
    cols = _ensure_f64_channel_cols!(df, ch_labels)

    if is_regular
        if parallel_channels && n_channels > 1
            _foreach_chunk(n_channels) do ch_range
                local bm_buf = Vector{Float64}(undef, n_freqs)
                local cnt_buf = Vector{Int}(undef, n_freqs)
                local work_buf = Vector{Float64}(undef, n_freqs)
                for ci in ch_range
                    _tf_baseline_grid!(reshape(cols[ci], n_freqs, n_times), base_range, method, bm_buf, cnt_buf, work_buf)
                end
            end
        else
            bm_seq = Vector{Float64}(undef, n_freqs)
            cnt_seq = Vector{Int}(undef, n_freqs)
            work_seq = Vector{Float64}(undef, n_freqs)
            for col in cols
                _tf_baseline_grid!(reshape(col, n_freqs, n_times), base_range, method, bm_seq, cnt_seq, work_seq)
            end
        end
    else
        row_index = _tf_build_row_index(df)
        power_3d = fill(NaN, n_channels, n_freqs, n_times)
        @inbounds for (fi, f) in enumerate(freqs_unique)
            rf = round(f, digits = 6)
            for (ti, t) in enumerate(times)
                row = get(row_index, (rf, round(t, digits = 6)), nothing)
                if !isnothing(row)
                    for ci = 1:n_channels
                        power_3d[ci, fi, ti] = cols[ci][row]
                    end
                end
            end
        end
        bm_irr = Vector{Float64}(undef, n_freqs)
        cnt_irr = Vector{Int}(undef, n_freqs)
        work_irr = Vector{Float64}(undef, n_freqs)
        for ci = 1:n_channels
            _tf_baseline_grid!(@view(power_3d[ci, :, :]), base_range, method, bm_irr, cnt_irr, work_irr)
        end
        @inbounds for (fi, f) in enumerate(freqs_unique)
            rf = round(f, digits = 6)
            for (ti, t) in enumerate(times)
                row = get(row_index, (rf, round(t, digits = 6)), nothing)
                if !isnothing(row)
                    for ci = 1:n_channels
                        cols[ci][row] = power_3d[ci, fi, ti]
                    end
                end
            end
        end
    end
    return nothing
end

function _tf_baseline_grid!(
    mat::AbstractMatrix{Float64},
    base_range::UnitRange{Int},
    method::Symbol,
    bm_buf::Vector{Float64} = Vector{Float64}(undef, size(mat, 1)),
    cnt_buf::Vector{Int} = Vector{Int}(undef, size(mat, 1)),
    work_buf::Vector{Float64} = Vector{Float64}(undef, size(mat, 1)),
)
    n_freqs, n_times = size(mat)
    fill!(bm_buf, 0.0)
    fill!(cnt_buf, 0)
    @inbounds for ti in base_range, fi = 1:n_freqs
        v = mat[fi, ti]
        if !isnan(v)
            bm_buf[fi] += v
            cnt_buf[fi] += 1
        end
    end
    @inbounds for fi = 1:n_freqs
        c = cnt_buf[fi]
        bm_buf[fi] = c > 0 ? (bm_buf[fi] / c) : NaN
    end

    if method === :db
        @inbounds for fi = 1:n_freqs
            bm = bm_buf[fi]
            work_buf[fi] = (isnan(bm) || bm <= 0.0) ? NaN : (1.0 / bm)
        end
        @inbounds for ti = 1:n_times, fi = 1:n_freqs
            inv_bm = work_buf[fi]
            p = mat[fi, ti]
            mat[fi, ti] = (isnan(inv_bm) || isnan(p) || p <= 0.0) ? NaN : 10.0 * log10(p * inv_bm)
        end
    elseif method === :absolute
        @inbounds for ti = 1:n_times
            @simd for fi = 1:n_freqs
                mat[fi, ti] -= bm_buf[fi]
            end
        end
    elseif method === :relative
        @inbounds for fi = 1:n_freqs
            bm = bm_buf[fi]
            work_buf[fi] = (isnan(bm) || bm == 0.0) ? NaN : (1.0 / bm)
        end
        @inbounds for ti = 1:n_times
            @simd for fi = 1:n_freqs
                mat[fi, ti] *= work_buf[fi]
            end
        end
    elseif method === :relchange || method === :percent
        num = method === :percent ? 100.0 : 1.0
        @inbounds for fi = 1:n_freqs
            bm = bm_buf[fi]
            work_buf[fi] = (isnan(bm) || bm == 0.0) ? NaN : (num / bm)
        end
        @inbounds for ti = 1:n_times
            @simd for fi = 1:n_freqs
                mat[fi, ti] = (mat[fi, ti] - bm_buf[fi]) * work_buf[fi]
            end
        end
    elseif method === :normchange
        @inbounds for ti = 1:n_times, fi = 1:n_freqs
            bm = bm_buf[fi]
            p = mat[fi, ti]
            denom = p + bm
            mat[fi, ti] = denom == 0.0 ? NaN : (p - bm) / denom
        end
    elseif method === :zscore
        fill!(work_buf, 0.0)
        @inbounds for ti in base_range, fi = 1:n_freqs
            v = mat[fi, ti]
            if !isnan(v)
                d = v - bm_buf[fi]
                work_buf[fi] += d * d
            end
        end
        @inbounds for fi = 1:n_freqs
            c = cnt_buf[fi]
            bstd = c >= 2 ? sqrt(work_buf[fi] / (c - 1)) : NaN
            work_buf[fi] = (isnan(bstd) || bstd == 0.0) ? NaN : (1.0 / bstd)
        end
        @inbounds for ti = 1:n_times
            @simd for fi = 1:n_freqs
                mat[fi, ti] = (mat[fi, ti] - bm_buf[fi]) * work_buf[fi]
            end
        end
    end
    return mat
end
