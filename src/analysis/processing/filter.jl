"""
    FilterInfo

Metadata and DSP filter object for a configured digital EEG filter.

# Fields
- `filter_type::Symbol`: `:hp` (highpass) or `:lp` (lowpass)
- `filter_object::Union{ZeroPoleGain,Vector{Float64}}`: `ZeroPoleGain` for IIR or tap vector for FIR
- `filter_method::Symbol`: `:iir` or `:fir`
- `cutoff_freq::Float64`: Cutoff frequency in Hz
- `sample_rate::Float64`: Sampling rate in Hz
- `order::Int`: Filter order (for IIR filters)
- `n_taps::Union{Nothing,Int}`: Tap count for FIR filters (`nothing` for IIR)
- `transition_band::Float64`: Transition band width in Hz
"""
struct FilterInfo
    filter_type::Symbol
    filter_object::Union{ZeroPoleGain,Vector{Float64}}
    filter_method::Symbol
    cutoff_freq::Float64
    sample_rate::Float64
    order::Int
    n_taps::Union{Nothing,Int}
    transition_band::Float64
end

"""
    create_lowpass_filter(cutoff_freq::Real, sample_rate::Real; filter_method::Symbol = :iir,
                          order::Integer = 2, transition_width = nothing,
                          transition_band_hz = nothing, min_taps::Int = 101) -> FilterInfo

Create a lowpass digital filter specification (`:iir` Butterworth or `:fir` Hamming-windowed).
When applied with zero-phase `filtfilt` (default), effective filter order is doubled.
"""
create_lowpass_filter(
    cutoff_freq::Real,
    sample_rate::Real;
    filter_method::Symbol = :iir,
    order::Integer = 2,
    transition_width::Union{Nothing,Real} = nothing,
    transition_band_hz::Union{Nothing,Real} = nothing,
    min_taps::Int = 101,
) = _create_filter(
    :lp,
    cutoff_freq,
    sample_rate;
    filter_method = filter_method,
    order = order,
    transition_width = transition_width,
    transition_band_hz = transition_band_hz,
    min_taps = min_taps,
)

"""
    create_highpass_filter(cutoff_freq::Real, sample_rate::Real; filter_method::Symbol = :iir,
                           order::Integer = 2, transition_width = nothing,
                           transition_band_hz = nothing, min_taps::Int = 101) -> FilterInfo

Create a highpass digital filter specification (`:iir` Butterworth or `:fir` Hamming-windowed).
When applied with zero-phase `filtfilt` (default), effective filter order is doubled.
"""
create_highpass_filter(
    cutoff_freq::Real,
    sample_rate::Real;
    filter_method::Symbol = :iir,
    order::Integer = 2,
    transition_width::Union{Nothing,Real} = nothing,
    transition_band_hz::Union{Nothing,Real} = nothing,
    min_taps::Int = 101,
) = _create_filter(
    :hp,
    cutoff_freq,
    sample_rate;
    filter_method = filter_method,
    order = order,
    transition_width = transition_width,
    transition_band_hz = transition_band_hz,
    min_taps = min_taps,
)

function _resolve_transition_band_hz(
    cutoff_freq::Real,
    sample_rate::Real,
    filter_type::Symbol,
    transition_width::Union{Nothing,Real},
)::Float64
    if isnothing(transition_width)
        bound = filter_type == :hp ? Float64(cutoff_freq) : (sample_rate / 2.0 - cutoff_freq)
        return min(max(0.25 * cutoff_freq, 2.0), bound)
    end
    transition_width <= 0 && @minimal_error("transition_width must be positive: $transition_width")
    return transition_width <= 1.0 ? Float64(cutoff_freq * transition_width) : Float64(transition_width)
end

function _create_filter(
    filter_type::Symbol,
    cutoff_freq::Real,
    sample_rate::Real;
    filter_method::Symbol = :iir,
    order::Integer = 2,
    transition_width::Union{Nothing,Real} = nothing,
    transition_band_hz::Union{Nothing,Real} = nothing,
    min_taps::Int = 101,
)
    filter_type in (:hp, :lp) || @minimal_error("filter_type '$filter_type' must be one of: (:hp, :lp)")
    filter_method in (:iir, :fir) || @minimal_error("filter_method '$filter_method' must be one of: (:iir, :fir)")
    order > 0 || @minimal_error("filter order must be positive: $order")
    sample_rate > 0 || @minimal_error("sample_rate must be positive: $sample_rate")
    (0 < cutoff_freq < sample_rate / 2) ||
        @minimal_error("cutoff_freq must be between 0 and Nyquist frequency ($(sample_rate/2) Hz): $cutoff_freq")

    prototype = filter_type == :hp ? Highpass(cutoff_freq) : Lowpass(cutoff_freq)
    transition_band = if !isnothing(transition_band_hz)
        transition_band_hz > 0 || @minimal_error("transition_band_hz must be positive: $transition_band_hz")
        Float64(transition_band_hz)
    else
        _resolve_transition_band_hz(cutoff_freq, sample_rate, filter_type, transition_width)
    end

    if filter_method == :iir
        fobj = digitalfilter(prototype, Butterworth(order); fs = sample_rate)
        return FilterInfo(filter_type, fobj, filter_method, cutoff_freq, sample_rate, order, nothing, transition_band)
    else
        n_taps = max(ceil(Int, 4.0 * sample_rate / transition_band), max(3, min_taps))
        iseven(n_taps) && (n_taps += 1)
        fobj = digitalfilter(prototype, FIRWindow(hamming(n_taps)); fs = sample_rate)
        return FilterInfo(filter_type, fobj, filter_method, cutoff_freq, sample_rate, order, n_taps, transition_band)
    end
end

function _ensure_fir_transition_width(
    dat::EegData,
    transition_width::Union{Nothing,Real},
    cutoff_freq::Real,
    filter_type::Symbol = :hp,
)::Float64
    min_length = dat.data isa Vector ? minimum(nrow, dat.data) : n_samples(dat)
    min_length < 10 && @minimal_error("Data length ($min_length) is too short for zero-phase FIR filtering (minimum 10 samples required)")
    max_taps_raw = (min_length - 1) ÷ 3
    max_odd_taps = max(3, iseven(max_taps_raw) ? max_taps_raw - 1 : max_taps_raw)
    min_req_tb = (4.0 * dat.sample_rate / max_odd_taps) * (1.0 + 1e-12)
    actual_tb = _resolve_transition_band_hz(cutoff_freq, dat.sample_rate, filter_type, transition_width)
    if actual_tb < min_req_tb
        @minimal_warning "FIR transition band $(round(actual_tb, digits=3)) Hz too narrow for data length $(min_length). Widening to $(round(min_req_tb, digits=3)) Hz to prevent crash."
        return min_req_tb
    end
    return actual_tb
end

@inline _prepared_filter_object(fobj::ZeroPoleGain) = convert(SecondOrderSections, fobj)
@inline _prepared_filter_object(fobj::Vector{Float64}) = fobj

@inline function _filter_column!(df::DataFrame, col_idx::Int, fobj, use_filtfilt::Bool)
    col = df[!, col_idx]
    if col isa Vector{Float64}
        col_f = col::Vector{Float64}
        col_f .= use_filtfilt ? filtfilt(fobj, col_f) : filt(fobj, col_f)
    else
        col .= use_filtfilt ? filtfilt(fobj, col) : filt(fobj, col)
    end
    return nothing
end

function _apply_filter!(dat::DataFrame, channels::AbstractVector{Symbol}, filter::FilterInfo; filter_func::Symbol = :filtfilt)::Nothing
    isempty(channels) && return nothing
    fobj = _prepared_filter_object(filter.filter_object)
    use_filtfilt = filter_func == :filtfilt
    col_indices = Int[columnindex(dat, ch) for ch in channels]
    if length(col_indices) > 1 && length(col_indices) * nrow(dat) >= 8_192
        Threads.@threads :static for idx in eachindex(col_indices)
            _filter_column!(dat, col_indices[idx], fobj, use_filtfilt)
        end
    else
        for col_idx in col_indices
            _filter_column!(dat, col_idx, fobj, use_filtfilt)
        end
    end
    return nothing
end

function _apply_filter!(
    dat::AbstractVector{DataFrame},
    channels::AbstractVector{Symbol},
    filter::FilterInfo;
    filter_func::Symbol = :filtfilt,
)::Nothing
    (isempty(dat) || isempty(channels)) && return nothing
    fobj = _prepared_filter_object(filter.filter_object)
    use_filtfilt = filter_func == :filtfilt
    col_indices = Int[columnindex(dat[1], ch) for ch in channels]
    if length(dat) > 1 && length(dat) * length(col_indices) * nrow(dat[1]) >= 8_192
        Threads.@threads :static for ep_idx in eachindex(dat)
            epoch = dat[ep_idx]
            for col_idx in col_indices
                _filter_column!(epoch, col_idx, fobj, use_filtfilt)
            end
        end
    else
        for epoch in dat, col_idx in col_indices
            _filter_column!(epoch, col_idx, fobj, use_filtfilt)
        end
    end
    return nothing
end

function _update_filter_info!(dat::EegData, filter_info::FilterInfo)::Nothing
    if filter_info.filter_type == :hp
        dat.analysis_info.hp_filter = filter_info.cutoff_freq
    elseif filter_info.filter_type == :lp
        dat.analysis_info.lp_filter = filter_info.cutoff_freq
    end
    return nothing
end

function _filter_eeg!(
    dat::EegData,
    cutoff_freq::Real,
    filter_type::Symbol;
    order::Int,
    transition_width::Union{Nothing,Real} = nothing,
    filter_method::Symbol = :iir,
    channel_selection = channels(),
    filter_func::Symbol = :filtfilt,
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for $(filter_type == :lp ? "lowpass" : "highpass") filtering"
        return dat
    end

    transition_band_hz = nothing
    min_taps = 101
    if filter_method == :fir && filter_func == :filtfilt
        transition_band_hz = _ensure_fir_transition_width(dat, transition_width, cutoff_freq, filter_type)
        min_length = dat.data isa Vector ? minimum(nrow, dat.data) : n_samples(dat)
        min_taps = min_length >= 304 ? 101 : 3
    end

    filter_info = _create_filter(
        filter_type,
        cutoff_freq,
        dat.sample_rate;
        filter_method = filter_method,
        order = order,
        transition_width = transition_width,
        transition_band_hz = transition_band_hz,
        min_taps = min_taps,
    )
    _update_filter_info!(dat, filter_info)
    _apply_filter!(dat.data, selected_channels, filter_info; filter_func = filter_func)
    return dat
end

"""
    lowpass_filter!(dat::EegData, cutoff_freq::Real; order::Int = 3, transition_width = nothing,
                    filter_method::Symbol = :iir, channel_selection = channels(), filter_func::Symbol = :filtfilt)
    lowpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :lowpass)

Apply a digital lowpass filter to `dat` (or a `Vector` of `EegData`) in-place.
"""
lowpass_filter!(
    dat::EegData,
    cutoff_freq::Real;
    order::Int = 3,
    transition_width::Union{Nothing,Real} = nothing,
    filter_method::Symbol = :iir,
    channel_selection = channels(),
    filter_func::Symbol = :filtfilt,
) = _filter_eeg!(
    dat,
    cutoff_freq,
    :lp;
    order = order,
    transition_width = transition_width,
    filter_method = filter_method,
    channel_selection = channel_selection,
    filter_func = filter_func,
)

"""
    highpass_filter!(dat::EegData, cutoff_freq::Real; order::Int = 1, transition_width = nothing,
                     filter_method::Symbol = :iir, channel_selection = channels(), filter_func::Symbol = :filtfilt)
    highpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :highpass)

Apply a digital highpass filter to `dat` (or a `Vector` of `EegData`) in-place.
"""
highpass_filter!(
    dat::EegData,
    cutoff_freq::Real;
    order::Int = 1,
    transition_width::Union{Nothing,Real} = nothing,
    filter_method::Symbol = :iir,
    channel_selection = channels(),
    filter_func::Symbol = :filtfilt,
) = _filter_eeg!(
    dat,
    cutoff_freq,
    :hp;
    order = order,
    transition_width = transition_width,
    filter_method = filter_method,
    channel_selection = channel_selection,
    filter_func = filter_func,
)

"""
    lowpass_filter(dat, cutoff_freq::Real; kwargs...)
    lowpass_filter(file_pattern::String, cutoff_freq::Real; input_dir = pwd(),
                   participant_selection = participants(), condition_selection = conditions(), output_dir = nothing)

Non-mutating or batch JLD2 version of [`lowpass_filter!`](@ref).
"""
lowpass_filter(dat::EegData, args...; kwargs...) = lowpass_filter!(copy(dat), args...; kwargs...)

"""
    highpass_filter(dat, cutoff_freq::Real; kwargs...)
    highpass_filter(file_pattern::String, cutoff_freq::Real; input_dir = pwd(),
                    participant_selection = participants(), condition_selection = conditions(), output_dir = nothing)

Non-mutating or batch JLD2 version of [`highpass_filter!`](@ref).
"""
highpass_filter(dat::EegData, args...; kwargs...) = highpass_filter!(copy(dat), args...; kwargs...)

function lowpass_filter!(dat::AbstractVector{<:EegData}, cutoff_freq::Real; kwargs...)
    for d in dat
        lowpass_filter!(d, cutoff_freq; kwargs...)
    end
    return dat
end

function highpass_filter!(dat::AbstractVector{<:EegData}, cutoff_freq::Real; kwargs...)
    for d in dat
        highpass_filter!(d, cutoff_freq; kwargs...)
    end
    return dat
end

lowpass_filter(dat::AbstractVector{<:EegData}, cutoff_freq::Real; kwargs...) = [lowpass_filter(d, cutoff_freq; kwargs...) for d in dat]
highpass_filter(dat::AbstractVector{<:EegData}, cutoff_freq::Real; kwargs...) = [highpass_filter(d, cutoff_freq; kwargs...) for d in dat]

function lowpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :lowpass)
    sec = getfield(filter_cfg, section)
    sec.apply && lowpass_filter!(dat, sec.freq; order = sec.order, filter_method = sec.method, filter_func = sec.func)
    return dat
end

function lowpass_filter!(dat::AbstractVector{<:EegData}, filter_cfg::FilterConfig; section::Symbol = :lowpass)
    for d in dat
        lowpass_filter!(d, filter_cfg; section = section)
    end
    return dat
end

function highpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :highpass)
    sec = getfield(filter_cfg, section)
    sec.apply && highpass_filter!(dat, sec.freq; order = sec.order, filter_method = sec.method, filter_func = sec.func)
    return dat
end

function highpass_filter!(dat::AbstractVector{<:EegData}, filter_cfg::FilterConfig; section::Symbol = :highpass)
    for d in dat
        highpass_filter!(d, filter_cfg; section = section)
    end
    return dat
end

lowpass_filter(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :lowpass) =
    lowpass_filter!(copy(dat), filter_cfg; section = section)
highpass_filter(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :highpass) =
    highpass_filter!(copy(dat), filter_cfg; section = section)
lowpass_filter(dat::AbstractVector{<:EegData}, filter_cfg::FilterConfig; section::Symbol = :lowpass) =
    [lowpass_filter(d, filter_cfg; section = section) for d in dat]
highpass_filter(dat::AbstractVector{<:EegData}, filter_cfg::FilterConfig; section::Symbol = :highpass) =
    [highpass_filter(d, filter_cfg; section = section) for d in dat]

function _fir_freqresp(h::Vector{Float64}, w::AbstractVector{Float64})::Vector{ComplexF64}
    n_pts = length(w)
    m = length(h)
    resp = Vector{ComplexF64}(undef, n_pts)
    @inbounds for idx = 1:n_pts
        u = cis(-w[idx])
        acc = ComplexF64(h[m])
        for k = (m-1):-1:1
            acc = muladd(acc, u, h[k])
        end
        resp[idx] = acc
    end
    return resp
end

function _find_db_crossings(mag_db::Vector{Float64}, target_db::Float64, threshold::Float64 = 0.5)::Vector{Int}
    n = length(mag_db)
    crossings = Int[]
    best_idx = 1
    best_diff = abs(mag_db[1] - target_db)
    @inbounds for i = 2:(n-1)
        d = abs(mag_db[i] - target_db)
        if d < best_diff
            best_diff = d
            best_idx = i
        end
        if d < threshold && d <= abs(mag_db[i-1] - target_db) && d <= abs(mag_db[i+1] - target_db)
            push!(crossings, i)
        end
    end
    if n >= 2
        d_end = abs(mag_db[n] - target_db)
        d_end < best_diff && (best_idx = n)
    end
    isempty(crossings) && push!(crossings, best_idx)
    return crossings
end

"""
    get_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000, filter_func = :filtfilt)

Compute key frequency-response characteristics (`cutoff_freq_3db`, `cutoff_freq_6db`, `passband_ripple`,
`stopband_atten`, `phase_delay`, `group_delay`, `transition_band`, `filter_type`) for `filter_info`.
When `filter_func = :filtfilt` (default), dB attenuation reflects two-pass zero-phase filtering.
"""
function get_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000, filter_func::Union{Symbol,Function} = :filtfilt)
    freqs = exp10.(range(log10(0.01), log10(filter_info.sample_rate / 2), length = npoints))
    w = (2π / filter_info.sample_rate) .* freqs

    resp = if filter_info.filter_method == :fir
        _fir_freqresp(filter_info.filter_object::Vector{Float64}, w)
    else
        freqresp(filter_info.filter_object, w)
    end

    is_zero_phase = (filter_func == :filtfilt || filter_func === filtfilt)
    db_scale = is_zero_phase ? 40.0 : 20.0
    mag_db = Vector{Float64}(undef, npoints)
    phase = Vector{Float64}(undef, npoints)
    @inbounds for i = 1:npoints
        mag_db[i] = db_scale * log10(abs(resp[i]))
        phase[i] = angle(resp[i])
    end

    cutoff_freq_3db = freqs[_find_db_crossings(mag_db, -3.0)]
    cutoff_freq_6db = freqs[_find_db_crossings(mag_db, -6.0)]
    fc = cutoff_freq_3db[1]

    if filter_info.filter_type == :lp
        passband_mask = freqs .<= fc
        stopband_mask = freqs .>= (fc + filter_info.transition_band)
    else
        passband_mask = freqs .>= fc
        stopband_mask = freqs .<= (fc - filter_info.transition_band)
    end

    passband_ripple = any(passband_mask) ? maximum(abs, @view mag_db[passband_mask]) : NaN
    stopband_atten = if any(stopband_mask)
        inf_floor = is_zero_phase ? -200.0 : -100.0
        s = 0.0
        cnt = 0
        @inbounds for i = 1:npoints
            if stopband_mask[i]
                v = mag_db[i]
                s += (v == -Inf ? inf_floor : v)
                cnt += 1
            end
        end
        s / cnt
    else
        NaN
    end

    unwrapped = DSP.unwrap(phase)
    group_delay = Vector{Float64}(undef, max(0, npoints - 1))
    @inbounds for i = 1:(npoints-1)
        group_delay[i] = -(unwrapped[i+1] - unwrapped[i]) / (w[i+1] - w[i])
    end
    phase_delay = -phase ./ w
    length(phase_delay) >= 2 && (phase_delay[1] = phase_delay[2])
    gd_mask = @view passband_mask[1:(end-1)]

    return (
        cutoff_freq_3db = cutoff_freq_3db,
        cutoff_freq_6db = cutoff_freq_6db,
        passband_ripple = passband_ripple,
        stopband_atten = stopband_atten,
        phase_delay = any(passband_mask) ? mean(@view phase_delay[passband_mask]) : NaN,
        group_delay = any(gd_mask) ? mean(@view group_delay[gd_mask]) : NaN,
        transition_band = filter_info.transition_band,
        filter_type = filter_info.filter_type,
    )
end

"""
    print_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000, filter_func = :filtfilt)

Log a formatted summary of filter cutoff, transition band, stopband attenuation, and delay characteristics.
"""
function print_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000, filter_func::Union{Symbol,Function} = :filtfilt)
    chars = get_filter_characteristics(filter_info; npoints = npoints, filter_func = filter_func)
    @info "--------------------------------------------------------------------------------"
    @info "Filter Characteristics: $(uppercase(string(chars.filter_type))) $(filter_info.cutoff_freq) Hz"
    @info "Cutoff: -3 dB = $(round.(chars.cutoff_freq_3db, digits=2)) Hz, -6 dB = $(round.(chars.cutoff_freq_6db, digits=2)) Hz"
    @info "Transition band: $(round(chars.transition_band, digits=1)) Hz"
    @info "Stopband attenuation: $(round(chars.stopband_atten, digits=1)) dB"
    @info "Single-pass characteristics: Group delay = $(round(chars.group_delay * 1000/filter_info.sample_rate, digits=1)) ms, Phase delay = $(round(chars.phase_delay * 1000/filter_info.sample_rate, digits=1)) ms "
    return nothing
end

_default_filter_output_dir(input_dir::String, pattern::String, filter_type::Symbol, freq::Real) =
    joinpath(input_dir, "filtered_$(clean_pattern(pattern))_$(filter_type)_$(freq)hz")

function _process_filter_file(filepath::String, output_path::String, filter_type::Symbol, cutoff_freq::Real, condition_selection)
    filename = basename(filepath)
    data = read_data(filepath)
    isnothing(data) && return BatchResult(false, filename, "No data variables found")
    data isa AbstractVector{<:Union{ErpData,EpochData}} ||
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")

    data = _condition_select(data, condition_selection)
    filter_type == :lp ? lowpass_filter!(data, cutoff_freq) : highpass_filter!(data, cutoff_freq)
    jldsave(output_path; data = data, format_version = EEGFUN_FORMAT_VERSION)
    return BatchResult(true, filename, "Filtered successfully")
end

lowpass_filter(
    file_pattern::String,
    cutoff_freq::Real;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
) = _run_filter_batch(
    file_pattern,
    cutoff_freq,
    :lp;
    input_dir = input_dir,
    participant_selection = participant_selection,
    condition_selection = condition_selection,
    output_dir = output_dir,
)

highpass_filter(
    file_pattern::String,
    cutoff_freq::Real;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
) = _run_filter_batch(
    file_pattern,
    cutoff_freq,
    :hp;
    input_dir = input_dir,
    participant_selection = participant_selection,
    condition_selection = condition_selection,
    output_dir = output_dir,
)

function _run_filter_batch(
    file_pattern::String,
    cutoff_freq::Real,
    filter_type::Symbol;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging(
        "$(filter_type == :lp ? "lowpass" : "highpass") filtering",
        "filter_$(filter_type).log",
        input_dir,
        output_dir,
    ) do ctx
        @info "  cutoff: $cutoff_freq Hz"
        cutoff_freq <= 0 && @minimal_error("Cutoff frequency must be positive, got: $cutoff_freq")
        output_dir = _prepare_output_dir!(ctx, _default_filter_output_dir(input_dir, file_pattern, filter_type, cutoff_freq))
        process_fn =
            (input_path, output_path) -> _process_filter_file(input_path, output_path, filter_type, cutoff_freq, condition_selection)
        return batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Filtering")
    end
end
