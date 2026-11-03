"""Human-readable labels paired with measurement type strings, for use in GUI menus."""
const MEASUREMENT_TYPE_LABELS = [
    "Mean Amplitude" => "mean_amplitude",
    "Max Peak Amplitude" => "max_peak_amplitude",
    "Min Peak Amplitude" => "min_peak_amplitude",
    "Max Peak Latency" => "max_peak_latency",
    "Min Peak Latency" => "min_peak_latency",
    "Peak-to-Peak Amplitude" => "peak_to_peak_amplitude",
    "Peak-to-Peak Latency" => "peak_to_peak_latency",
    "Rectified Area" => "rectified_area",
    "Integral" => "integral",
    "Positive Area" => "positive_area",
    "Negative Area" => "negative_area",
    "Fractional Area Latency" => "fractional_area_latency",
    "Fractional Peak Latency" => "fractional_peak_latency",
]

"""Valid ERP measurement type strings."""
const VALID_MEASUREMENT_TYPES = last.(MEASUREMENT_TYPE_LABELS)

const ERP_MEASUREMENTS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :local_interval => (
        3,
        "Number of samples on each side of peak (total interval = 2*local_interval + 1). Peak must be larger than neighbors and local averages within this window.",
    ),
    :fractional_area_fraction =>
        (0.5, "Fraction for fractional area latency (0.0-1.0). Finds latency where this fraction of area is to the left."),
    :area_boundary =>
        (:rectified, "Area boundary mode for fractional area latency: :rectified (default), :positive, :negative, or :signed"),
    :fractional_peak_fraction =>
        (0.5, "Fraction for fractional peak latency (0.0-1.0). Finds latency where amplitude is this fraction of peak."),
    :fractional_peak_direction => (:onset, "Direction for fractional peak latency: :onset (before peak) or :offset (after peak)"),
    :peak_polarity =>
        (:abs, "Peak polarity for fractional peak latency: :abs (default, largest absolute peak), :positive, or :negative"),
)

struct _MeasurementSpec
    mode::Symbol
    local_interval::Int
    fractional_area_fraction::Float64
    area_boundary::Symbol
    fractional_peak_fraction::Float64
    fractional_peak_direction::Symbol
    peak_polarity::Symbol
end

function _MeasurementSpec(analysis_type::Union{AbstractString,Symbol}, measurement_kwargs::Dict{Symbol,Any})
    t_str = string(analysis_type)
    t_str ∈ VALID_MEASUREMENT_TYPES || @minimal_error "Analysis type must be one of: $(join(VALID_MEASUREMENT_TYPES, ", ")). Got: $t_str"

    local_interval = measurement_kwargs[:local_interval]
    (local_interval isa Integer && local_interval >= 1) || @minimal_error "local_interval must be >= 1, got: $local_interval"

    faf = measurement_kwargs[:fractional_area_fraction]
    (faf isa Real && 0.0 <= faf <= 1.0) || @minimal_error "fractional_area_fraction must be in [0.0, 1.0], got: $faf"

    ab = measurement_kwargs[:area_boundary]
    ab ∈ (:rectified, :positive, :negative, :signed) ||
        @minimal_error "area_boundary must be :rectified, :positive, :negative, or :signed, got: $ab"

    fpf = measurement_kwargs[:fractional_peak_fraction]
    (fpf isa Real && 0.0 <= fpf <= 1.0) || @minimal_error "fractional_peak_fraction must be in [0.0, 1.0], got: $fpf"

    fpd = measurement_kwargs[:fractional_peak_direction]
    (fpd === :onset || fpd === :offset) || @minimal_error "fractional_peak_direction must be :onset or :offset, got: $fpd"

    pp = measurement_kwargs[:peak_polarity]
    pp ∈ (:abs, :positive, :negative) || @minimal_error "peak_polarity must be :abs, :positive, or :negative, got: $pp"

    return _MeasurementSpec(Symbol(t_str), Int(local_interval), Float64(faf), ab, Float64(fpf), fpd, pp)
end

function _validate_measurement_params(
    analysis_type::Union{AbstractString,Symbol},
    analysis_interval::Interval,
    measurement_kwargs::Dict{Symbol,Any},
)
    spec = _MeasurementSpec(analysis_type, measurement_kwargs)
    if !isnothing(analysis_interval) && !(analysis_interval isa AllSelection) && !(analysis_interval isa Function)
        iv = _normalize_interval(analysis_interval)
        iv[1] > iv[2] && @minimal_error "Analysis interval start must be <= stop. Got: ($(iv[1]), $(iv[2]))"
    end
    return spec
end

_default_measurements_output_dir(input_dir::String, analysis_type::Union{AbstractString,Symbol}) =
    joinpath(input_dir, "measurements_$(string(analysis_type))")

function _normalize_interval(interval::Interval)
    isnothing(interval) && return nothing
    if interval isa AbstractRange
        return (Float64(first(interval)), Float64(last(interval)))
    elseif interval isa AbstractVector{<:Real}
        length(interval) != 2 && @minimal_error "Interval vector must contain exactly 2 elements (start, stop), got $(length(interval))"
        return (Float64(interval[1]), Float64(interval[2]))
    elseif interval isa TimeSelection
        return (Float64(interval.start), Float64(interval.stop))
    elseif interval isa Tuple{Real,Real}
        return (Float64(interval[1]), Float64(interval[2]))
    else
        return interval
    end
end

# === SCALAR MEASUREMENT KERNELS ===

function _find_robust_peak(data::AbstractVector, peak_type::Symbol; local_interval::Int = 3)
    n = length(data)
    total_window = 2 * local_interval + 1
    if total_window > n
        @minimal_warning "local_interval ($local_interval) requires $(total_window) samples but data has only $n, cannot detect robust peak"
        return (nothing, nothing)
    elseif local_interval > n ÷ 2
        @minimal_warning "local_interval ($local_interval) is > 50% of data length ($n), may be too large for robust peak detection"
    end

    found = false
    best_val = 0.0
    best_idx = 0
    inv_w = 1.0 / local_interval

    @inbounds for i = (local_interval+1):(n-local_interval)
        val = Float64(data[i])
        if peak_type === :max
            (val > data[i-1] && val > data[i+1]) || continue
            (!found || val > best_val) || continue
            left_sum = 0.0
            right_sum = 0.0
            for k = 1:local_interval
                left_sum += data[i-k]
                right_sum += data[i+k]
            end
            if val > left_sum * inv_w && val > right_sum * inv_w
                best_val = val
                best_idx = i
                found = true
            end
        else
            (val < data[i-1] && val < data[i+1]) || continue
            (!found || val < best_val) || continue
            left_sum = 0.0
            right_sum = 0.0
            for k = 1:local_interval
                left_sum += data[i-k]
                right_sum += data[i+k]
            end
            if val < left_sum * inv_w && val < right_sum * inv_w
                best_val = val
                best_idx = i
                found = true
            end
        end
    end

    return found ? (best_val, best_idx) : (nothing, nothing)
end

function _compute_peak_measurement(chan_data::AbstractVector, peak_type::Symbol, local_interval::Int, channel_name::Symbol)
    total_window = 2 * local_interval + 1
    if total_window > length(chan_data)
        peak_name = peak_type === :max ? "maximum" : "minimum"
        @minimal_warning "Channel $channel_name: local_interval ($local_interval, total interval $total_window) > analysis interval length ($(length(chan_data)) samples). Cannot detect robust peak, using simple $peak_name."
    end

    peak_val, peak_idx = _find_robust_peak(chan_data, peak_type; local_interval = local_interval)
    if isnothing(peak_val)
        peak_name = peak_type === :max ? "maximum" : "minimum"
        @minimal_warning "Channel $channel_name: No robust $peak_name peak found, using simple $peak_name"
        v, idx = peak_type === :max ? findmax(chan_data) : findmin(chan_data)
        return (Float64(v), idx)
    end
    return (peak_val, peak_idx)
end

@inline function _area_sample(x::Float64, mode::Symbol)::Float64
    mode === :rectified && return abs(x)
    mode === :positive && return max(x, 0.0)
    mode === :negative && return abs(min(x, 0.0))
    mode === :signed && return x
    @minimal_error "area_boundary must be :rectified, :positive, :negative, or :signed, got: $mode"
end

function _trapezoidal_area(data::AbstractVector, time_col::AbstractVector, mode::Symbol)::Float64
    n = min(length(data), length(time_col))
    n <= 1 && return 0.0
    total = 0.0
    prev_y = _area_sample(Float64(data[1]), mode)
    prev_t = Float64(time_col[1])
    @inbounds for i = 2:n
        cur_y = _area_sample(Float64(data[i]), mode)
        cur_t = Float64(time_col[i])
        total += 0.5 * (prev_y + cur_y) * (cur_t - prev_t)
        prev_y = cur_y
        prev_t = cur_t
    end
    return total
end

function _fractional_area_latency(data::AbstractVector, time_col::AbstractVector, fraction::Real; area_boundary::Symbol = :rectified)
    (isempty(data) || isempty(time_col)) && return NaN
    n = min(length(data), length(time_col))
    n == 1 && return Float64(time_col[1])
    t_first = Float64(time_col[1])
    t_last = Float64(time_col[n])
    t_last - t_first <= 0.0 && return t_first

    total_area = _trapezoidal_area(data, time_col, area_boundary)
    abs(total_area) < 1e-12 && return t_first

    frac = Float64(fraction)
    frac <= 0.0 && return t_first
    frac >= 1.0 && return t_last

    target_area = total_area * frac
    cum = 0.0
    prev_y = _area_sample(Float64(data[1]), area_boundary)
    prev_t = t_first

    @inbounds for i = 2:n
        cur_y = _area_sample(Float64(data[i]), area_boundary)
        cur_t = Float64(time_col[i])
        y0 = cum
        cum = y0 + 0.5 * (prev_y + cur_y) * (cur_t - prev_t)
        crossed = total_area > 0.0 ? (cum >= target_area) : (cum <= target_area)
        if crossed
            abs(cum - y0) < 1e-14 && return prev_t
            return prev_t + (target_area - y0) * (cur_t - prev_t) / (cum - y0)
        end
        prev_y = cur_y
        prev_t = cur_t
    end
    return t_last
end

function _fractional_peak_latency(data::AbstractVector, time_col::AbstractVector, peak_idx::Int, fraction::Float64, direction::Symbol)
    (isempty(data) || isempty(time_col) || peak_idx < 1 || peak_idx > length(data)) && return NaN
    peak_val = Float64(data[peak_idx])
    abs(peak_val) < 1e-12 && return Float64(time_col[peak_idx])
    target_val = peak_val * fraction

    if direction === :onset
        @inbounds for i = (peak_idx-1):-1:1
            yi = Float64(data[i])
            if (peak_val > 0.0 && yi <= target_val) || (peak_val < 0.0 && yi >= target_val)
                y1 = Float64(data[i+1])
                x0, x1 = Float64(time_col[i]), Float64(time_col[i+1])
                abs(y1 - yi) < 1e-14 && return x0
                return x0 + (target_val - yi) * (x1 - x0) / (y1 - yi)
            end
        end
        return Float64(time_col[1])
    else
        @inbounds for i = (peak_idx+1):length(data)
            yi = Float64(data[i])
            if (peak_val > 0.0 && yi <= target_val) || (peak_val < 0.0 && yi >= target_val)
                y0 = Float64(data[i-1])
                x0, x1 = Float64(time_col[i-1]), Float64(time_col[i])
                abs(yi - y0) < 1e-14 && return x1
                return x0 + (target_val - y0) * (x1 - x0) / (yi - y0)
            end
        end
        return Float64(time_col[end])
    end
end

function _eval_channel_measurement(
    chan_data::AbstractVector,
    selected_times::AbstractVector,
    spec::_MeasurementSpec,
    channel_name::Symbol,
)::Float64
    m = spec.mode
    if m === :mean_amplitude
        s = 0.0
        @inbounds @simd for i in eachindex(chan_data)
            s += Float64(chan_data[i])
        end
        return s / length(chan_data)
    elseif m === :max_peak_amplitude || m === :min_peak_amplitude || m === :max_peak_latency || m === :min_peak_latency
        peak_type = (m === :max_peak_amplitude || m === :max_peak_latency) ? :max : :min
        peak_val, peak_idx = _compute_peak_measurement(chan_data, peak_type, spec.local_interval, channel_name)
        return (m === :max_peak_amplitude || m === :min_peak_amplitude) ? peak_val : Float64(selected_times[peak_idx])
    elseif m === :peak_to_peak_amplitude || m === :peak_to_peak_latency
        max_val, max_idx = _compute_peak_measurement(chan_data, :max, spec.local_interval, channel_name)
        min_val, min_idx = _compute_peak_measurement(chan_data, :min, spec.local_interval, channel_name)
        return m === :peak_to_peak_amplitude ? (max_val - min_val) :
               abs(Float64(selected_times[max_idx]) - Float64(selected_times[min_idx]))
    elseif m === :rectified_area || m === :integral || m === :positive_area || m === :negative_area
        (length(selected_times) <= 1 || length(chan_data) <= 1) && return 0.0
        dt = (Float64(selected_times[end]) - Float64(selected_times[1])) / (length(selected_times) - 1)
        if dt <= 0.0
            @minimal_warning "Channel $channel_name: time step is <= 0, returning 0.0 for area measurement"
            return 0.0
        end
        bmode = m === :rectified_area ? :rectified : m === :integral ? :signed : m === :positive_area ? :positive : :negative
        return _trapezoidal_area(chan_data, selected_times, bmode)
    elseif m === :fractional_area_latency
        return _fractional_area_latency(chan_data, selected_times, spec.fractional_area_fraction; area_boundary = spec.area_boundary)
    elseif m === :fractional_peak_latency
        pol = spec.peak_polarity
        peak_idx = if pol === :positive
            _, idx = _compute_peak_measurement(chan_data, :max, spec.local_interval, channel_name)
            idx
        elseif pol === :negative
            _, idx = _compute_peak_measurement(chan_data, :min, spec.local_interval, channel_name)
            idx
        else
            max_v, max_i = _compute_peak_measurement(chan_data, :max, spec.local_interval, channel_name)
            min_v, min_i = _compute_peak_measurement(chan_data, :min, spec.local_interval, channel_name)
            abs(max_v) >= abs(min_v) ? max_i : min_i
        end
        return _fractional_peak_latency(chan_data, selected_times, peak_idx, spec.fractional_peak_fraction, spec.fractional_peak_direction)
    end
    @minimal_error "Unknown analysis_type: '$(spec.mode)'. Valid types: $(join(VALID_MEASUREMENT_TYPES, ", "))"
end

function _compute_measurement(
    chan_data::AbstractVector,
    selected_times::AbstractVector,
    analysis_type::String,
    measurement_kwargs::Dict{Symbol,Any},
    channel_name::Symbol,
)
    return _eval_channel_measurement(chan_data, selected_times, _MeasurementSpec(analysis_type, measurement_kwargs), channel_name)
end

# === TABLE EXTRACTION ENGINE ===

function _resolve_in_memory_participant(dat::Union{ErpData,EpochData}, participant::Union{Int,Nothing})::Int
    !isnothing(participant) && return participant
    first_df = dat isa ErpData ? dat.data : (isempty(dat.data) ? nothing : dat.data[1])
    if !isnothing(first_df) && hasproperty(first_df, :participant) && !isempty(first_df)
        return Int(first_df[1, :participant])
    end
    stem, _ = splitext(basename(dat.file))
    return any(isdigit, stem) ? _extract_participant_id(basename(dat.file)) : 0
end

function _apply_baseline_correction!(dfs::AbstractVector{DataFrame}, baseline_interval::Interval, all_channels::AbstractVector{Symbol})
    (isnothing(baseline_interval) || isempty(dfs)) && return dfs
    first_df = first(dfs)
    eeg_channels = Symbol[ch for ch in all_channels if hasproperty(first_df, ch)]
    if isempty(eeg_channels)
        @minimal_warning "No EEG channels found for baseline correction"
        return dfs
    end
    baseline_int = _normalize_interval(baseline_interval)
    start_idx, stop_idx = _find_idx_start_end(first_df.time, baseline_int[1], baseline_int[2])
    if start_idx <= stop_idx
        @info "Applying baseline correction to $(length(eeg_channels)) channels using interval: $baseline_int"
        _apply_baseline!(dfs, eeg_channels, (start_idx, stop_idx))
    else
        @minimal_warning "No samples found in baseline interval $baseline_int; skipping baseline correction"
    end
    return dfs
end

function _fill_channel_measurements!(
    out_col::Vector{Float64},
    dfs,
    valid_idxs::Vector{Int},
    time_vecs::Vector{Vector{Float64}},
    time_ranges::Vector{UnitRange{Int}},
    bl_ranges::Vector{UnitRange{Int}},
    ch::Symbol,
    col_idx::Int,
    spec::_MeasurementSpec,
    apply_bl::Bool,
    chan_buf::Vector{Float64},
)
    for r in eachindex(valid_idxs)
        df = dfs[valid_idxs[r]]
        t_range = time_ranges[r]
        raw_col = (col_idx > 0 && ncol(df) >= col_idx) ? df[!, col_idx] : df[!, ch]
        t_col = time_vecs[r]
        sel_times = @view t_col[t_range]
        bl_r = bl_ranges[r]

        if apply_bl && !isempty(bl_r)
            n_w = length(t_range)
            resize!(chan_buf, n_w)
            offset = first(t_range) - 1
            if raw_col isa Vector{Float64}
                rc = raw_col::Vector{Float64}
                b_sum = 0.0
                @inbounds @simd for i in bl_r
                    b_sum += rc[i]
                end
                b_mean = b_sum / length(bl_r)
                @inbounds @simd for j = 1:n_w
                    chan_buf[j] = rc[offset+j] - b_mean
                end
            else
                b_mean = mean(@view raw_col[bl_r])
                @inbounds for j = 1:n_w
                    chan_buf[j] = Float64(raw_col[offset+j]) - b_mean
                end
            end
            out_col[r] = _eval_channel_measurement(chan_buf, sel_times, spec, ch)
        else
            if raw_col isa Vector{Float64}
                rc = raw_col::Vector{Float64}
                out_col[r] = _eval_channel_measurement(@view(rc[t_range]), sel_times, spec, ch)
            else
                out_col[r] = _eval_channel_measurement(@view(raw_col[t_range]), sel_times, spec, ch)
            end
        end
    end
end

function _extract_measurements_df(
    dat::Union{ErpData,EpochData},
    analysis_type_str::String,
    spec::_MeasurementSpec;
    analysis_interval::Interval = times(),
    baseline_interval::Interval = nothing,
    channel_selection = channels(),
    participant::Int = 0,
    average_channels::Bool = false,
    include_extra::Bool = false,
    apply_baseline::Bool = true,
)::Union{DataFrame,Nothing}
    dfs = dat isa ErpData ? (dat.data,) : dat.data
    isempty(dfs) && return nothing

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = include_extra)
    if isempty(selected_channels)
        @minimal_warning "No channels selected"
        return nothing
    end

    use_all_times = isnothing(analysis_interval) || analysis_interval isa AllSelection
    ai_norm = use_all_times ? nothing : _normalize_interval(analysis_interval)
    bl_norm = isnothing(baseline_interval) ? nothing : _normalize_interval(baseline_interval)

    n_dfs = length(dfs)
    valid_idxs = Int[]
    sizehint!(valid_idxs, n_dfs)
    time_ranges = UnitRange{Int}[]
    sizehint!(time_ranges, n_dfs)
    bl_ranges = UnitRange{Int}[]
    sizehint!(bl_ranges, n_dfs)
    ai_starts = Float64[]
    sizehint!(ai_starts, n_dfs)
    ai_ends = Float64[]
    sizehint!(ai_ends, n_dfs)

    for k = 1:n_dfs
        df = dfs[k]
        if isempty(df)
            @minimal_warning "DataFrame is empty"
            continue
        end
        hasproperty(df, :time) || @minimal_error "DataFrame must have a :time column"
        t_col = df.time

        t_range = if use_all_times
            1:nrow(df)
        else
            s_i, e_i = _find_idx_start_end(t_col, ai_norm[1], ai_norm[2])
            s_i:e_i
        end
        if isempty(t_range)
            @minimal_warning "No time points found matching analysis interval"
            continue
        end

        bl_r = if apply_baseline && !isnothing(bl_norm)
            bs_i, be_i = _find_idx_start_end(t_col, bl_norm[1], bl_norm[2])
            bs_i:be_i
        else
            1:0
        end

        push!(valid_idxs, k)
        push!(time_ranges, t_range)
        push!(bl_ranges, bl_r)
        push!(ai_starts, use_all_times ? Float64(t_col[first(t_range)]) : ai_norm[1])
        push!(ai_ends, use_all_times ? Float64(t_col[last(t_range)]) : ai_norm[2])
    end

    n_valid = length(valid_idxs)
    n_valid == 0 && return nothing

    out_df = DataFrame(; copycols = false)
    out_df[!, :file] = fill(dat.file, n_valid)
    out_df[!, :participant] = fill(participant, n_valid)
    out_df[!, :condition] = [hasproperty(dfs[k], :condition) ? dfs[k][1, :condition] : dat.condition for k in valid_idxs]
    out_df[!, :condition_name] =
        String[hasproperty(dfs[k], :condition_name) ? String(dfs[k][1, :condition_name]) : dat.condition_name for k in valid_idxs]
    if hasproperty(dfs[valid_idxs[1]], :epoch)
        out_df[!, :epoch] = [dfs[k][1, :epoch] for k in valid_idxs]
    end
    out_df[!, :analysis_type] = fill(analysis_type_str, n_valid)
    out_df[!, :baseline_interval_start] = isnothing(bl_norm) ? fill(missing, n_valid) : fill(bl_norm[1], n_valid)
    out_df[!, :baseline_interval_end] = isnothing(bl_norm) ? fill(missing, n_valid) : fill(bl_norm[2], n_valid)
    out_df[!, :analysis_interval_start] = ai_starts
    out_df[!, :analysis_interval_end] = ai_ends

    first_valid_df = dfs[valid_idxs[1]]
    time_vecs = Vector{Float64}[
        let t = dfs[k].time
            t isa Vector{Float64} ? t : Float64.(t)
        end for k in valid_idxs
    ]
    chan_buf = Float64[]

    if average_channels
        dv_col = zeros(Float64, n_valid)
        tmp_col = Vector{Float64}(undef, n_valid)
        for ch in selected_channels
            col_idx = columnindex(first_valid_df, ch)
            _fill_channel_measurements!(
                tmp_col,
                dfs,
                valid_idxs,
                time_vecs,
                time_ranges,
                bl_ranges,
                ch,
                col_idx,
                spec,
                apply_baseline,
                chan_buf,
            )
            dv_col .+= tmp_col
        end
        dv_col ./= length(selected_channels)
        out_df[!, :dv_channels] = dv_col
    else
        for ch in selected_channels
            col_idx = columnindex(first_valid_df, ch)
            col_out = Vector{Float64}(undef, n_valid)
            _fill_channel_measurements!(
                col_out,
                dfs,
                valid_idxs,
                time_vecs,
                time_ranges,
                bl_ranges,
                ch,
                col_idx,
                spec,
                apply_baseline,
                chan_buf,
            )
            out_df[!, ch] = col_out
        end
    end

    return out_df
end

# === PUBLIC API ===

"""
    erp_measurements!(dat::ErpData, analysis_type; analysis_interval = times(), baseline_interval = nothing, channel_selection = channels(), participant = 0, kwargs...)
    erp_measurements!(dat::EpochData, analysis_type; kwargs...)
    erp_measurements!(data::AbstractVector{<:Union{ErpData,EpochData}}, analysis_type; kwargs...)
    erp_measurements(data::Union{ErpData,EpochData,AbstractVector{<:Union{ErpData,EpochData}}}, analysis_type; kwargs...) -> Union{ErpMeasurementsResult,Nothing}
    erp_measurements(file_pattern::String, analysis_type; analysis_interval = times(), baseline_interval = nothing,
                     participant_selection = participants(), condition_selection = conditions(), channel_selection = channels(),
                     input_dir = pwd(), output_dir = nothing, output_file = "erp_measurements", kwargs...)

Extract ERP amplitude, latency, or area measurements from EEG data.
- `erp_measurements` on in-memory `ErpData`, `EpochData`, or `Vector` returns an [`ErpMeasurementsResult`](@ref) without mutating input data.
- `erp_measurements!` mutates input data in-place when applying `baseline_interval` / condition columns and returns raw `NamedTuple` row(s).
- `erp_measurements(file_pattern, ...)` batch-processes JLD2 files and saves a CSV.

# `analysis_type` Options
- `"mean_amplitude"` — mean amplitude in `analysis_interval`
- `"max_peak_amplitude"` / `"min_peak_amplitude"` — robust local peak amplitude (with simple extremum fallback)
- `"max_peak_latency"` / `"min_peak_latency"` — peak latency in seconds
- `"peak_to_peak_amplitude"` / `"peak_to_peak_latency"` — difference between max and min peaks
- `"integral"` / `"rectified_area"` / `"positive_area"` / `"negative_area"` — trapezoidal area in µV·s
- `"fractional_area_latency"` / `"fractional_peak_latency"` — interpolated onset/offset latencies

$(_generate_kwargs_doc(ERP_MEASUREMENTS_KWARGS))
"""
function erp_measurements!(
    dat::Union{ErpData,EpochData},
    analysis_type::Union{AbstractString,Symbol};
    analysis_interval::Interval = times(),
    baseline_interval::Interval = nothing,
    channel_selection = channels(),
    participant::Int = 0,
    average_channels::Bool = false,
    include_extra::Bool = false,
    kwargs...,
)
    analysis_type_str = string(analysis_type)
    measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)
    spec = _validate_measurement_params(analysis_type_str, analysis_interval, measurement_kwargs)

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = include_extra)
    if isempty(selected_channels)
        @minimal_warning "No channels selected"
        return nothing
    end
    all_channels = include_extra ? union(channel_labels(dat), selected_channels) : channel_labels(dat)

    dfs = dat isa ErpData ? [dat.data] : dat.data
    _apply_baseline_correction!(dfs, baseline_interval, all_channels)
    for df in dfs
        hasproperty(df, :condition) || insertcols!(df, 1, :condition => dat.condition)
        hasproperty(df, :condition_name) || insertcols!(df, 2, :condition_name => dat.condition_name)
    end

    res_df = _extract_measurements_df(
        dat,
        analysis_type_str,
        spec;
        analysis_interval = analysis_interval,
        baseline_interval = baseline_interval,
        channel_selection = channel_selection,
        participant = participant,
        average_channels = average_channels,
        include_extra = include_extra,
        apply_baseline = false,
    )
    isnothing(res_df) && return nothing
    rows = Tables.rowtable(res_df)
    return dat isa ErpData ? only(rows) : rows
end

function erp_measurements!(
    data::AbstractVector{<:Union{ErpData,EpochData}},
    analysis_type::Union{AbstractString,Symbol};
    average_channels::Bool = false,
    kwargs...,
)
    if isempty(data)
        measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)
        _validate_measurement_params(analysis_type, get(kwargs, :analysis_interval, times()), measurement_kwargs)
        return nothing
    end

    all_results = NamedTuple[]
    for dat in data
        result = erp_measurements!(dat, analysis_type; average_channels = average_channels, kwargs...)
        isnothing(result) && continue
        if result isa AbstractVector
            append!(all_results, result)
        else
            push!(all_results, result)
        end
    end
    return isempty(all_results) ? nothing : all_results
end

"""
    erp_measurements(data::Union{ErpData,EpochData}, analysis_type; kwargs...) -> Union{ErpMeasurementsResult,Nothing}
    erp_measurements(data::AbstractVector{<:Union{ErpData,EpochData}}, analysis_type; condition_selection = conditions(), kwargs...) -> Union{ErpMeasurementsResult,Nothing}

Non-mutating in-memory ERP measurement extraction returning an [`ErpMeasurementsResult`](@ref).
"""
function erp_measurements(
    dat::Union{ErpData,EpochData},
    analysis_type::Union{AbstractString,Symbol};
    analysis_interval::Interval = times(),
    baseline_interval::Interval = nothing,
    channel_selection = channels(),
    participant::Union{Int,Nothing} = nothing,
    average_channels::Bool = false,
    include_extra::Bool = false,
    kwargs...,
)::Union{ErpMeasurementsResult,Nothing}
    analysis_type_str = string(analysis_type)
    measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)
    spec = _validate_measurement_params(analysis_type_str, analysis_interval, measurement_kwargs)
    part_id = _resolve_in_memory_participant(dat, participant)

    res_df = _extract_measurements_df(
        dat,
        analysis_type_str,
        spec;
        analysis_interval = analysis_interval,
        baseline_interval = baseline_interval,
        channel_selection = channel_selection,
        participant = part_id,
        average_channels = average_channels,
        include_extra = include_extra,
        apply_baseline = true,
    )
    isnothing(res_df) && return nothing
    return ErpMeasurementsResult(res_df, analysis_type_str, analysis_interval, baseline_interval)
end

function erp_measurements(
    data::AbstractVector{<:Union{ErpData,EpochData}},
    analysis_type::Union{AbstractString,Symbol};
    analysis_interval::Interval = times(),
    baseline_interval::Interval = nothing,
    condition_selection = conditions(),
    channel_selection = channels(),
    participant::Union{Int,Nothing} = nothing,
    average_channels::Bool = false,
    include_extra::Bool = false,
    kwargs...,
)::Union{ErpMeasurementsResult,Nothing}
    analysis_type_str = string(analysis_type)
    measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)
    spec = _validate_measurement_params(analysis_type_str, analysis_interval, measurement_kwargs)

    isempty(data) && return nothing
    selected_indices = get_selected_conditions(data, condition_selection)
    if isempty(selected_indices)
        @minimal_warning "No conditions selected for ERP measurements"
        return nothing
    end

    dfs = DataFrame[]
    for idx in selected_indices
        dat = data[idx]
        part_id = _resolve_in_memory_participant(dat, participant)
        df = _extract_measurements_df(
            dat,
            analysis_type_str,
            spec;
            analysis_interval = analysis_interval,
            baseline_interval = baseline_interval,
            channel_selection = channel_selection,
            participant = part_id,
            average_channels = average_channels,
            include_extra = include_extra,
            apply_baseline = true,
        )
        !isnothing(df) && push!(dfs, df)
    end

    isempty(dfs) && return nothing
    combined = length(dfs) == 1 ? dfs[1] : vcat(dfs...; cols = :union)
    return ErpMeasurementsResult(combined, analysis_type_str, analysis_interval, baseline_interval)
end

function _process_measurements_file(
    filepath::String,
    analysis_interval::Interval,
    analysis_type_str::String,
    spec::_MeasurementSpec,
    baseline_interval::Interval,
    average_channels::Bool,
    condition_selection,
    channel_selection;
    include_extra::Bool = false,
)::Union{DataFrame,Nothing}
    filename = basename(filepath)
    participant = _extract_participant_id(filename)

    data_var = read_data(filepath)
    if isnothing(data_var)
        @minimal_warning "No data variables found in $filename"
        return nothing
    end

    items = if data_var isa Union{ErpData,EpochData}
        [data_var]
    elseif data_var isa AbstractVector{<:Union{ErpData,EpochData}}
        data_var
    else
        @minimal_warning "Invalid data type in $filename: expected Vector{ErpData} or Vector{EpochData}, got $(typeof(data_var))"
        return nothing
    end

    selected_indices = get_selected_conditions(items, condition_selection)
    isempty(selected_indices) && @minimal_error "No conditions left following condition selection in $filename"

    dfs = DataFrame[]
    for idx in selected_indices
        df = _extract_measurements_df(
            items[idx],
            analysis_type_str,
            spec;
            analysis_interval = analysis_interval,
            baseline_interval = baseline_interval,
            channel_selection = channel_selection,
            participant = participant,
            average_channels = average_channels,
            include_extra = include_extra,
            apply_baseline = true,
        )
        !isnothing(df) && push!(dfs, df)
    end
    isempty(dfs) && return nothing
    return length(dfs) == 1 ? dfs[1] : vcat(dfs...; cols = :union)
end

function erp_measurements(
    file_pattern::String,
    analysis_type::Union{AbstractString,Symbol};
    analysis_interval::Interval = times(),
    baseline_interval::Interval = nothing,
    participant_selection = participants(),
    condition_selection = conditions(),
    channel_selection = channels(),
    average_channels::Bool = false,
    include_extra::Bool = false,
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    output_file::String = "erp_measurements",
    kwargs...,
)
    _with_batch_logging("ERP measurement analysis", "$(output_file).log", input_dir, output_dir) do ctx
        @log_call "erp_measurements"

        analysis_type_str = string(analysis_type)
        measurement_kwargs = _merge_plot_kwargs(ERP_MEASUREMENTS_KWARGS, kwargs)
        spec = _validate_measurement_params(analysis_type_str, analysis_interval, measurement_kwargs)

        output_dir = _prepare_output_dir!(ctx, _default_measurements_output_dir(input_dir, analysis_type_str))
        @info "Analysis type: $analysis_type_str"

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files"

        file_dfs = DataFrame[]
        processed_count = 0
        error_count = 0

        p = Progress(length(files); desc = "Measuring ERPs...")
        for (file_idx, file) in enumerate(files)
            input_path = joinpath(input_dir, file)
            @info "Processing: $file ($file_idx/$(length(files)))"
            try
                file_df = _process_measurements_file(
                    input_path,
                    analysis_interval,
                    analysis_type_str,
                    spec,
                    baseline_interval,
                    average_channels,
                    condition_selection,
                    channel_selection;
                    include_extra = include_extra,
                )
                if !isnothing(file_df) && !isempty(file_df)
                    push!(file_dfs, file_df)
                    processed_count += 1
                    @info "  ✓ Extracted $(nrow(file_df)) measurement(s)"
                else
                    error_count += 1
                    @minimal_warning "  ✗ No measurements extracted from $file"
                end
            catch e
                @minimal_warning "Error processing $file: $(sprint(showerror, e))"
                error_count += 1
            end
            next!(p)
        end

        if isempty(file_dfs)
            @minimal_warning "No results to save"
            @info "Analysis complete! Processed $processed_count files successfully, $error_count errors"
            return nothing
        end

        results_df = length(file_dfs) == 1 ? file_dfs[1] : vcat(file_dfs...; cols = :union)
        sort_cols = hasproperty(results_df, :epoch) ? [:participant, :condition, :epoch] : [:participant, :condition]
        sort!(results_df, sort_cols)

        output_csv = joinpath(output_dir, "$(output_file).csv")
        CSV.write(output_csv, results_df)

        @info "Saved $(nrow(results_df)) measurement(s) to: $output_csv"
        @info "Analysis complete! Processed $processed_count files successfully, $error_count errors"
        return ErpMeasurementsResult(results_df, analysis_type_str, analysis_interval, baseline_interval)
    end
end
