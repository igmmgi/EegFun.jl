const _VALID_TF_BASELINE_METHODS = (:absolute, :relative, :relchange, :normchange, :db, :zscore, :percent)
const _TFDataObject = Union{TimeFreqData,TimeFreqEpochData}
const _TFDataOrVec = Union{TimeFreqData,TimeFreqEpochData,AbstractVector{<:TimeFreqData},AbstractVector{<:TimeFreqEpochData}}
const _TFAltInterval = Union{AbstractVector{<:Real},TimeSelection}

function _validate_tf_baseline_method(method::Symbol)
    if method ∉ _VALID_TF_BASELINE_METHODS
        @minimal_error("Unknown baseline method: $method. Use :absolute, :relative, :relchange, :normchange, :db, :zscore, or :percent")
    end
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

function _validate_and_resolve_tf_baseline_times(times::Vector{Float64}, baseline_interval::Tuple{Real,Real})
    b_start, b_stop = Float64(baseline_interval[1]), Float64(baseline_interval[2])
    (isnan(b_start) || isnan(b_stop) || b_start > b_stop) &&
        @minimal_error("Invalid baseline interval $(baseline_interval): start must be <= stop")
    isempty(times) && @minimal_error("Time-frequency data contains no time points")

    t_min, t_max = extrema(times)
    (b_stop < t_min || b_start > t_max) &&
        @minimal_error("Baseline window $(baseline_interval) does not overlap with data times ($(t_min) to $(t_max) s)")

    i1 = find_closest_time_index(times, b_start)
    i2 = find_closest_time_index(times, b_stop)
    return minmax(i1, i2)
end

@inline function _baseline_df_col!(df::DataFrame, ch::Symbol, n_freqs::Int, n_times::Int, base_range::UnitRange{Int}, method::Symbol)
    col = df[!, ch]
    if col isa Vector{Float64}
        _tf_baseline_grid!(reshape(col, n_freqs, n_times), base_range, method)
    else
        col_f64 = Float64.(col)
        _tf_baseline_grid!(reshape(col_f64, n_freqs, n_times), base_range, method)
        df[!, ch] = col_f64
    end
    return nothing
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

    if is_regular
        if parallel_channels && n_channels > 1
            Threads.@threads for ch in ch_labels
                _baseline_df_col!(df, ch, n_freqs, n_times, base_range, method)
            end
        else
            for ch in ch_labels
                _baseline_df_col!(df, ch, n_freqs, n_times, base_range, method)
            end
        end
    else
        row_index = _tf_build_row_index(df)
        power_3d = fill(NaN, n_channels, n_freqs, n_times)
        for (fi, f) in enumerate(freqs_unique)
            rf = round(f, digits = 6)
            for (ti, t) in enumerate(times)
                row = get(row_index, (rf, round(t, digits = 6)), nothing)
                if !isnothing(row)
                    for (ci, ch) in enumerate(ch_labels)
                        power_3d[ci, fi, ti] = Float64(df[row, ch])
                    end
                end
            end
        end
        for ci = 1:n_channels
            _tf_baseline_grid!(@view(power_3d[ci, :, :]), base_range, method)
        end
        for (fi, f) in enumerate(freqs_unique)
            rf = round(f, digits = 6)
            for (ti, t) in enumerate(times)
                row = get(row_index, (rf, round(t, digits = 6)), nothing)
                if !isnothing(row)
                    for (ci, ch) in enumerate(ch_labels)
                        df[row, ch] = power_3d[ci, fi, ti]
                    end
                end
            end
        end
    end
    return nothing
end

function _tf_baseline_grid!(mat::AbstractMatrix{Float64}, base_range::UnitRange{Int}, method::Symbol)
    n_freqs, n_times = size(mat)
    n_base = length(base_range)

    @inbounds for fi = 1:n_freqs
        s = 0.0
        cnt = 0
        for ti in base_range
            v = mat[fi, ti]
            if !isnan(v)
                s += v
                cnt += 1
            end
        end
        bm = cnt > 0 ? s / cnt : NaN
        if isnan(bm)
            for ti = 1:n_times
                mat[fi, ti] = NaN
            end
            continue
        end

        if method === :absolute
            @simd for ti = 1:n_times
                mat[fi, ti] -= bm
            end
        elseif method === :relative
            if bm == 0.0
                for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_bm = 1.0 / bm
                @simd for ti = 1:n_times
                    mat[fi, ti] *= inv_bm
                end
            end
        elseif method === :relchange || method === :percent
            if bm == 0.0
                for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                scale = (method === :percent ? 100.0 : 1.0) / bm
                @simd for ti = 1:n_times
                    mat[fi, ti] = (mat[fi, ti] - bm) * scale
                end
            end
        elseif method === :normchange
            for ti = 1:n_times
                p = mat[fi, ti]
                denom = p + bm
                mat[fi, ti] = denom == 0.0 ? NaN : (p - bm) / denom
            end
        elseif method === :db
            if bm <= 0.0
                for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_bm = 1.0 / bm
                for ti = 1:n_times
                    p = mat[fi, ti]
                    mat[fi, ti] = (isnan(p) || p <= 0.0) ? NaN : 10.0 * log10(p * inv_bm)
                end
            end
        elseif method === :zscore
            ssq = 0.0
            for ti in base_range
                v = mat[fi, ti]
                if !isnan(v)
                    d = v - bm
                    ssq += d * d
                end
            end
            bstd = cnt >= 2 ? sqrt(ssq / (cnt - 1)) : NaN
            if isnan(bstd) || bstd == 0.0
                for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_std = 1.0 / bstd
                @simd for ti = 1:n_times
                    mat[fi, ti] = (mat[fi, ti] - bm) * inv_std
                end
            end
        end
    end
    return mat
end

"""
    tf_baseline!(tf_data::TimeFreqData, baseline_interval; method=:db)
    tf_baseline!(tf_data::TimeFreqEpochData, baseline_interval; method=:db)

Apply baseline correction to `TimeFreqData` or single-trial `TimeFreqEpochData` in-place.

Supported `method` values: `:db` (`10 * log10(P / B)`), `:absolute` (`P - B`),
`:relative` (`P / B`), `:relchange` (`(P - B) / B`), `:percent` (`100 * (P - B) / B`),
`:normchange` (`(P - B) / (P + B)`), and `:zscore` (`(P - B) / SD_B`).

# Example
```julia
tf_baseline!(tf_data, (-0.3, 0.0); method=:db)
```
"""
function tf_baseline!(tf_data::TimeFreqData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    _check_tf_not_baselined(tf_data.baseline)
    _validate_tf_baseline_method(method)
    ch_labels = channel_labels(tf_data)
    @info "Applying $(method) baseline correction ($(baseline_interval[1])s to $(baseline_interval[2])s)"
    _tf_baseline_df!(tf_data.data_power, ch_labels, baseline_interval, method; parallel_channels = true)
    tf_data.baseline = BaselineInfo(method = method, window = (Float64(baseline_interval[1]), Float64(baseline_interval[2])))
    return tf_data
end

function tf_baseline!(tf_data::TimeFreqEpochData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    _check_tf_not_baselined(tf_data.baseline)
    _validate_tf_baseline_method(method)
    isempty(tf_data.data_power) && @minimal_error("Cannot baseline-correct empty TimeFreqEpochData")

    df1 = tf_data.data_power[1]
    times = unique(df1.time)
    freqs_unique = unique(df1.freq)
    n_freqs, n_times = length(freqs_unique), length(times)
    b1, b2 = _validate_and_resolve_tf_baseline_times(times, baseline_interval)
    base_range = b1:b2
    n_grid = n_freqs * n_times
    regular_ref = (nrow(df1) == n_grid) && (n_grid > 0) && (df1.time[1] == df1.time[min(n_freqs, end)])

    ch_labels = channel_labels(tf_data)
    @info "Applying $(method) baseline correction to $(length(tf_data.data_power)) TF epochs ($(baseline_interval[1])s to $(baseline_interval[2])s)"

    Threads.@threads for ep_idx in eachindex(tf_data.data_power)
        df = tf_data.data_power[ep_idx]
        if regular_ref && nrow(df) == n_grid && df.time[1] == times[1] && df.time[end] == times[end]
            for ch in ch_labels
                _baseline_df_col!(df, ch, n_freqs, n_times, base_range, method)
            end
        else
            _tf_baseline_df!(df, ch_labels, baseline_interval, method; parallel_channels = false)
        end
    end
    tf_data.baseline = BaselineInfo(method = method, window = (Float64(baseline_interval[1]), Float64(baseline_interval[2])))
    return tf_data
end

# Non-mutating, vector, interval-conversion, and `baseline`/`baseline!` dispatches
tf_baseline(tf_data::_TFDataObject, baseline_interval::Tuple{Real,Real}; method::Symbol = :db) =
    tf_baseline!(copy(tf_data), baseline_interval; method = method)

function tf_baseline!(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    for tf in tf_data
        tf_baseline!(tf, baseline_interval; method = method)
    end
    return tf_data
end

tf_baseline(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db) =
    [tf_baseline(tf, baseline_interval; method = method) for tf in tf_data]

_to_tf_interval(interval::Tuple{Real,Real}) = (Float64(interval[1]), Float64(interval[2]))
_to_tf_interval(interval::TimeSelection) = (interval.start, interval.stop)
function _to_tf_interval(interval::AbstractVector{<:Real})
    length(interval) == 2 || @minimal_error("Baseline interval vector must have 2 elements [start, stop], got: $interval")
    return (Float64(interval[1]), Float64(interval[2]))
end

function _full_tf_interval(tf_data::TimeFreqData)
    isempty(tf_data.data_power.time) && @minimal_error("Time-frequency data contains no time points")
    return extrema(tf_data.data_power.time)
end

function _full_tf_interval(tf_data::TimeFreqEpochData)
    (isempty(tf_data.data_power) || isempty(tf_data.data_power[1].time)) && @minimal_error("TimeFreqEpochData contains no time points")
    return extrema(tf_data.data_power[1].time)
end

tf_baseline!(tf_data::_TFDataOrVec, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, _to_tf_interval(baseline_interval); method = method)
tf_baseline(tf_data::_TFDataOrVec, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline(tf_data, _to_tf_interval(baseline_interval); method = method)

tf_baseline!(tf_data::_TFDataObject, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline!(tf_data, _full_tf_interval(tf_data); method = method)
function tf_baseline!(tf_data::AbstractVector{<:_TFDataObject}, ::AllSelection = times(); method::Symbol = :db)
    for tf in tf_data
        tf_baseline!(tf, _full_tf_interval(tf); method = method)
    end
    return tf_data
end
tf_baseline(tf_data::_TFDataObject, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline(tf_data, _full_tf_interval(tf_data); method = method)
tf_baseline(tf_data::AbstractVector{<:_TFDataObject}, ::AllSelection = times(); method::Symbol = :db) =
    [tf_baseline(tf, _full_tf_interval(tf); method = method) for tf in tf_data]

const _TFInterval = Union{Tuple{Real,Real},_TFAltInterval}

baseline!(tf_data::_TFDataObject, baseline_interval::_TFInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, baseline_interval; method = method)
baseline!(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::_TFInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, baseline_interval; method = method)
baseline(tf_data::_TFDataObject, baseline_interval::_TFInterval; method::Symbol = :db) =
    tf_baseline(tf_data, baseline_interval; method = method)
baseline(tf_data::AbstractVector{<:_TFDataObject}, baseline_interval::_TFInterval; method::Symbol = :db) =
    tf_baseline(tf_data, baseline_interval; method = method)
baseline!(tf_data::_TFDataObject, sel::AllSelection = times(); method::Symbol = :db) = tf_baseline!(tf_data, sel; method = method)
baseline!(tf_data::AbstractVector{<:_TFDataObject}, sel::AllSelection = times(); method::Symbol = :db) =
    tf_baseline!(tf_data, sel; method = method)
baseline(tf_data::_TFDataObject, sel::AllSelection = times(); method::Symbol = :db) = tf_baseline(tf_data, sel; method = method)
baseline(tf_data::AbstractVector{<:_TFDataObject}, sel::AllSelection = times(); method::Symbol = :db) =
    tf_baseline(tf_data, sel; method = method)

"""
    tf_baseline(file_pattern::String, baseline_interval;
                method=:db, input_dir=pwd(), output_dir=nothing,
                participant_selection=participants(), condition_selection=conditions())

Apply baseline correction to `TimeFreqData` or `TimeFreqEpochData` `.jld2` files in batch.

# Example
```julia
tf_baseline("tf_epochs_wavelet", (-0.3, 0.0); method=:db)
```
"""
function tf_baseline(
    file_pattern::String,
    baseline_interval::Tuple{Real,Real};
    method::Symbol = :db,
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection = participants(),
    condition_selection = conditions(),
)
    _with_batch_logging("TF baseline", "tf_baseline.log", input_dir, output_dir) do ctx
        @log_call "tf_baseline"
        output_dir = _prepare_output_dir!(ctx, joinpath(input_dir, "tf_baseline_$(clean_pattern(file_pattern))"))
        @info "baseline: $baseline_interval, method: $method"

        process_fn = (input_path, output_path) -> begin
            filename = basename(input_path)
            data = read_data(input_path)
            if isnothing(data) || !(data isa Vector{<:Union{TimeFreqData,TimeFreqEpochData}})
                return BatchResult(false, filename, "Invalid data type")
            end
            data = _condition_select(data, condition_selection)
            tf_baseline!(data, baseline_interval; method = method)
            jldsave(output_path; data = data, format_version = EEGFUN_FORMAT_VERSION)
            return BatchResult(true, filename, "Baseline corrected")
        end

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "TF baseline")
    end
end
