const _VALID_TF_BASELINE_METHODS = (:absolute, :relative, :relchange, :normchange, :vssum, :db, :zscore, :percent)

function _validate_tf_baseline_method(method::Symbol)
    if method ∉ _VALID_TF_BASELINE_METHODS
        @minimal_error(
            "Unknown baseline method: $method. Use :absolute, :relative, :relchange, :normchange, :db, :vssum (alias for :normchange), :zscore, or :percent",
        )
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
    b_start = Float64(baseline_interval[1])
    b_stop = Float64(baseline_interval[2])
    if isnan(b_start) || isnan(b_stop) || b_start > b_stop
        @minimal_error("Invalid baseline interval $(baseline_interval): start must be <= stop")
    end
    if isempty(times)
        @minimal_error("Time-frequency data contains no time points")
    end
    t_min, t_max = extrema(times)
    if b_stop < t_min || b_start > t_max
        @minimal_error("Baseline window $(baseline_interval) does not overlap with data times ($(t_min) to $(t_max) s)")
    end

    baseline_start_idx = find_closest_time_index(times, b_start)
    baseline_end_idx = find_closest_time_index(times, b_stop)
    if baseline_start_idx > baseline_end_idx
        baseline_start_idx, baseline_end_idx = baseline_end_idx, baseline_start_idx
    end
    return baseline_start_idx, baseline_end_idx
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
    n_freqs = length(freqs_unique)
    n_times = length(times)
    n_channels = length(ch_labels)

    baseline_start_idx, baseline_end_idx = _validate_and_resolve_tf_baseline_times(times, baseline_interval)
    base_range = baseline_start_idx:baseline_end_idx

    is_regular_grid = (nrow(df) == n_freqs * n_times) && (n_freqs > 0 && n_times > 0) && (df.time[1] == df.time[min(n_freqs, end)])

    if is_regular_grid
        if parallel_channels && n_channels > 1
            Threads.@threads for ch in ch_labels
                col = df[!, ch]
                if col isa Vector{Float64}
                    _tf_baseline_grid!(reshape(col, n_freqs, n_times), base_range, method)
                else
                    col_f64 = Float64.(col)
                    _tf_baseline_grid!(reshape(col_f64, n_freqs, n_times), base_range, method)
                    df[!, ch] = col_f64
                end
            end
        else
            for ch in ch_labels
                col = df[!, ch]
                if col isa Vector{Float64}
                    _tf_baseline_grid!(reshape(col, n_freqs, n_times), base_range, method)
                else
                    col_f64 = Float64.(col)
                    _tf_baseline_grid!(reshape(col_f64, n_freqs, n_times), base_range, method)
                    df[!, ch] = col_f64
                end
            end
        end
    else
        base_mask = falses(n_times)
        base_mask[base_range] .= true

        row_index = _tf_build_row_index(df)
        power_3d = Array{Float64,3}(undef, n_channels, n_freqs, n_times)
        for (fi, f) in enumerate(freqs_unique)
            rf = round(f, digits = 6)
            for (ti, t) in enumerate(times)
                row = get(row_index, (rf, round(t, digits = 6)), nothing)
                if !isnothing(row)
                    for (ci, ch) in enumerate(ch_labels)
                        power_3d[ci, fi, ti] = Float64(df[row, ch])
                    end
                else
                    power_3d[:, fi, ti] .= NaN
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

"""
    tf_baseline!(tf_data::TimeFreqData, baseline_interval; method=:db)
    tf_baseline!(tf_data::TimeFreqEpochData, baseline_interval; method=:db)

Apply baseline correction to `TimeFreqData` or single-trial `TimeFreqEpochData` in-place.

# Arguments
- `tf_data`: Time-frequency data (`TimeFreqData` or `TimeFreqEpochData`) to baseline correct
- `baseline_interval::Tuple{Real,Real}`: Time window for baseline `(start, stop)` in seconds

# Keyword Arguments
- `method::Symbol=:db`: Baseline method
  - `:absolute`: Absolute change (`power - baseline_mean`) - simple subtraction, no normalization
  - `:relative`: Relative power (`power / baseline_mean`) - ratio
  - `:relchange`: Relative change (`(power - baseline_mean) / baseline_mean`) - fractional change
  - `:normchange`: Normalized change (`(power - baseline) / (power + baseline)`) - symmetric normalization
  - `:db`: Decibel change (`10 * log10(power / baseline_mean)`)
  - `:vssum`: Alias for `:normchange`
  - `:zscore`: Z-score normalization (`(power - baseline_mean) / baseline_std`, sample SD)
  - `:percent`: Percent change (`100 * (power - baseline) / baseline`)

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

    # Validate interval on the first epoch up-front before spawning threads
    _validate_and_resolve_tf_baseline_times(unique(tf_data.data_power[1].time), baseline_interval)

    ch_labels = channel_labels(tf_data)
    @info "Applying $(method) baseline correction to $(length(tf_data.data_power)) TF epochs ($(baseline_interval[1])s to $(baseline_interval[2])s)"

    Threads.@threads for ep_idx in eachindex(tf_data.data_power)
        _tf_baseline_df!(tf_data.data_power[ep_idx], ch_labels, baseline_interval, method; parallel_channels = false)
    end
    tf_data.baseline = BaselineInfo(method = method, window = (Float64(baseline_interval[1]), Float64(baseline_interval[2])))
    return tf_data
end

"""
    _tf_baseline_grid!(mat::AbstractMatrix{Float64}, base_range, method::Symbol)

Baseline-correct one channel stored as a regular `[freq × time]` grid, in place.
Kept as a separate function so that each thread gets its own local variables.
"""
function _tf_baseline_grid!(mat::AbstractMatrix{Float64}, base_range::UnitRange{Int}, method::Symbol)
    n_freqs, n_times = size(mat)
    n_base = length(base_range)
    inv_n_base = n_base > 0 ? 1.0 / n_base : NaN

    for fi = 1:n_freqs
        base_vals = @view mat[fi, base_range]
        has_nan = any(isnan, base_vals)
        bm = if has_nan
            valid = filter(!isnan, base_vals)
            isempty(valid) ? NaN : mean(valid)
        elseif n_base == 0
            NaN
        else
            s = 0.0
            @inbounds @simd for ti in base_range
                s += mat[fi, ti]
            end
            s * inv_n_base
        end

        if isnan(bm)
            @inbounds for ti = 1:n_times
                mat[fi, ti] = NaN
            end
            continue
        end

        if method == :absolute
            @inbounds @simd for ti = 1:n_times
                mat[fi, ti] -= bm
            end
        elseif method == :relative
            if bm == 0.0
                @inbounds for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_bm = 1.0 / bm
                @inbounds @simd for ti = 1:n_times
                    mat[fi, ti] *= inv_bm
                end
            end
        elseif method == :relchange
            if bm == 0.0
                @inbounds for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_bm = 1.0 / bm
                @inbounds @simd for ti = 1:n_times
                    mat[fi, ti] = (mat[fi, ti] - bm) * inv_bm
                end
            end
        elseif method == :percent
            if bm == 0.0
                @inbounds for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_bm = 100.0 / bm
                @inbounds @simd for ti = 1:n_times
                    mat[fi, ti] = (mat[fi, ti] - bm) * inv_bm
                end
            end
        elseif method == :normchange || method == :vssum
            @inbounds for ti = 1:n_times
                p = mat[fi, ti]
                denom = p + bm
                mat[fi, ti] = denom == 0.0 ? NaN : (p - bm) / denom
            end
        elseif method == :db
            if bm <= 0.0
                @inbounds for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_bm = 1.0 / bm
                @inbounds for ti = 1:n_times
                    p = mat[fi, ti]
                    mat[fi, ti] = (isnan(p) || p <= 0.0) ? NaN : 10.0 * log10(p * inv_bm)
                end
            end
        elseif method == :zscore
            bstd = if has_nan
                valid = filter(!isnan, base_vals)
                length(valid) < 2 ? NaN : std(valid, corrected = true)
            else
                n_base < 2 ? NaN : std(base_vals, corrected = true)
            end
            if isnan(bstd) || bstd == 0.0
                @inbounds for ti = 1:n_times
                    mat[fi, ti] = NaN
                end
            else
                inv_std = 1.0 / bstd
                @inbounds @simd for ti = 1:n_times
                    mat[fi, ti] = (mat[fi, ti] - bm) * inv_std
                end
            end
        end
    end
    return mat
end

# Non-mutating versions
function tf_baseline(tf_data::TimeFreqData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline!(copy(tf_data), baseline_interval; method = method)
end

function tf_baseline(tf_data::TimeFreqEpochData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline!(copy(tf_data), baseline_interval; method = method)
end

# Vector versions
function tf_baseline!(tf_data::AbstractVector{<:TimeFreqData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    for tf in tf_data
        tf_baseline!(tf, baseline_interval; method = method)
    end
    return tf_data
end

function tf_baseline!(tf_data::AbstractVector{<:TimeFreqEpochData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    for tf in tf_data
        tf_baseline!(tf, baseline_interval; method = method)
    end
    return tf_data
end

function tf_baseline(tf_data::AbstractVector{<:TimeFreqData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return [tf_baseline(tf, baseline_interval; method = method) for tf in tf_data]
end

function tf_baseline(tf_data::AbstractVector{<:TimeFreqEpochData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return [tf_baseline(tf, baseline_interval; method = method) for tf in tf_data]
end

# Overloads for baseline! and baseline for TimeFreqData & TimeFreqEpochData
function baseline!(tf_data::TimeFreqData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline!(tf_data, baseline_interval; method = method)
end

function baseline!(tf_data::TimeFreqEpochData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline!(tf_data, baseline_interval; method = method)
end

function baseline!(tf_data::AbstractVector{<:TimeFreqData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline!(tf_data, baseline_interval; method = method)
end

function baseline!(tf_data::AbstractVector{<:TimeFreqEpochData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline!(tf_data, baseline_interval; method = method)
end

function baseline(tf_data::TimeFreqData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline(tf_data, baseline_interval; method = method)
end

function baseline(tf_data::TimeFreqEpochData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline(tf_data, baseline_interval; method = method)
end

function baseline(tf_data::AbstractVector{<:TimeFreqData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline(tf_data, baseline_interval; method = method)
end

function baseline(tf_data::AbstractVector{<:TimeFreqEpochData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return tf_baseline(tf_data, baseline_interval; method = method)
end

# Interval normalization helpers for TF baseline
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

const _TFDataObject = Union{TimeFreqData,TimeFreqEpochData}
const _TFDataOrVec = Union{TimeFreqData,TimeFreqEpochData,AbstractVector{<:TimeFreqData},AbstractVector{<:TimeFreqEpochData}}
const _TFAltInterval = Union{AbstractVector{<:Real},TimeSelection}

tf_baseline!(tf_data::_TFDataObject, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, _to_tf_interval(baseline_interval); method = method)

tf_baseline!(tf_data::AbstractVector{<:TimeFreqData}, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, _to_tf_interval(baseline_interval); method = method)

tf_baseline!(tf_data::AbstractVector{<:TimeFreqEpochData}, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, _to_tf_interval(baseline_interval); method = method)

tf_baseline(tf_data::_TFDataOrVec, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline(tf_data, _to_tf_interval(baseline_interval); method = method)

baseline!(tf_data::_TFDataOrVec, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline!(tf_data, _to_tf_interval(baseline_interval); method = method)

baseline(tf_data::_TFDataOrVec, baseline_interval::_TFAltInterval; method::Symbol = :db) =
    tf_baseline(tf_data, _to_tf_interval(baseline_interval); method = method)

# Full-epoch / default interval overloads so 1-arg `baseline!(tf_epoch)` or `baseline!(tf_epoch, times())`
# never silently dispatches to ERP `baseline!(::EegData, ::Interval = times())`
tf_baseline!(tf_data::_TFDataObject, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline!(tf_data, _full_tf_interval(tf_data); method = method)

function tf_baseline!(tf_data::AbstractVector{<:TimeFreqData}, ::AllSelection = times(); method::Symbol = :db)
    for tf in tf_data
        tf_baseline!(tf, _full_tf_interval(tf); method = method)
    end
    return tf_data
end

function tf_baseline!(tf_data::AbstractVector{<:TimeFreqEpochData}, ::AllSelection = times(); method::Symbol = :db)
    for tf in tf_data
        tf_baseline!(tf, _full_tf_interval(tf); method = method)
    end
    return tf_data
end

tf_baseline(tf_data::_TFDataObject, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline(tf_data, _full_tf_interval(tf_data); method = method)

tf_baseline(tf_data::AbstractVector{<:TimeFreqData}, ::AllSelection = times(); method::Symbol = :db) =
    [tf_baseline(tf, _full_tf_interval(tf); method = method) for tf in tf_data]

tf_baseline(tf_data::AbstractVector{<:TimeFreqEpochData}, ::AllSelection = times(); method::Symbol = :db) =
    [tf_baseline(tf, _full_tf_interval(tf); method = method) for tf in tf_data]

baseline!(tf_data::_TFDataObject, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline!(tf_data, _full_tf_interval(tf_data); method = method)

baseline!(tf_data::AbstractVector{<:TimeFreqData}, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline!(tf_data, times(); method = method)

baseline!(tf_data::AbstractVector{<:TimeFreqEpochData}, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline!(tf_data, times(); method = method)

baseline(tf_data::_TFDataObject, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline(tf_data, _full_tf_interval(tf_data); method = method)

baseline(tf_data::AbstractVector{<:TimeFreqData}, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline(tf_data, times(); method = method)

baseline(tf_data::AbstractVector{<:TimeFreqEpochData}, ::AllSelection = times(); method::Symbol = :db) =
    tf_baseline(tf_data, times(); method = method)

"""
    tf_baseline(file_pattern::String, baseline_interval;
                method=:db, input_dir=pwd(), output_dir=nothing,
                participant_selection=participants(), condition_selection=conditions())

Apply baseline correction to TimeFreqData or TimeFreqEpochData files.

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
            jldsave(output_path; data = data)
            return BatchResult(true, filename, "Baseline corrected")
        end

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "TF baseline")
    end
end
