# Moving-window and per-epoch peak-to-peak amplitude artifact detection.

######################################
########### User Functions ###########
######################################

"""
    is_peak_to_peak(dat, threshold::Real, window_size::Real; channel_selection = channels(),
                    sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined)

Detect moving-window peak-to-peak artifacts without modifying `dat`. Returns a `Vector{Bool}`
(or per-epoch vector) in `:combined` mode, or a `DataFrame` (or per-epoch vector) in `:separate` mode.
"""
function is_peak_to_peak(dat::_DetectableEeg, threshold::Real, window_size::Real; mode::Symbol = :combined, kwargs...)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    return _detect_artifacts(_make_p2p_detector(threshold, window_samples), dat, ch -> Symbol("is_peak_to_peak_$(ch)"); mode, kwargs...)
end

is_peak_to_peak(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...) =
    [is_peak_to_peak(d, threshold, window_size; kwargs...) for d in dat]

"""
    is_peak_to_peak!(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                     channel_selection = channels(), sample_selection = samples(),
                     interval_selection = times(), mode::Symbol = :combined, channel_out = nothing)
    is_peak_to_peak!(dat::MultiDataFrameEeg, threshold::Real, window_size::Real; kwargs...)
    is_peak_to_peak!(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...)

Detect moving-window peak-to-peak artifacts (`max - min > threshold` within `window_size` seconds),
storing boolean flag columns in-place (`:is_peak_to_peak` or `:is_peak_to_peak_<ch>` by default).
"""
function is_peak_to_peak!(
    dat::_DetectableEeg,
    threshold::Real,
    window_size::Real;
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
    kwargs...,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    combined_col = something(channel_out, :is_peak_to_peak)
    sep_fn = ch -> Symbol("$(combined_col)_$(ch)")
    return _flag_artifacts!(_make_p2p_detector(threshold, window_samples), dat, combined_col, sep_fn; mode, kwargs...)
end

is_peak_to_peak!(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...) =
    (foreach(d -> is_peak_to_peak!(d, threshold, window_size; kwargs...), dat); dat)

"""
    n_peak_to_peak(dat, threshold::Real, window_size::Real; channel_selection = channels(),
                   sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined)

Count peak-to-peak artifact samples across selected channels, returning a total `Int` (`:combined`)
or a per-channel `DataFrame` (`:separate`).
"""
function n_peak_to_peak(dat::_DetectableEeg, threshold::Real, window_size::Real; mode::Symbol = :combined, kwargs...)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    return _count_artifacts(_make_p2p_detector(threshold, window_samples), dat, :n_peak_to_peak; mode, kwargs...)
end

n_peak_to_peak(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...) =
    [n_peak_to_peak(d, threshold, window_size; kwargs...) for d in dat]

######################################
######### Internal Functions #########
######################################

@inline function _make_p2p_detector(threshold::Real, window_samples::Int)
    w = max(1, window_samples)
    thresh_f = Float64(threshold)
    return function ()
        ws = ntuple(_ -> Vector{Float64}(undef, w), 6)
        return (buf, sig) -> _is_peak_to_peak!(buf, sig, thresh_f, w, ws)
    end
end

function _is_peak_to_peak!(
    mask::AbstractVector{Bool},
    signal::AbstractVector{<:Real},
    threshold::Real,
    window_samples::Int,
    ws::Union{Nothing,NTuple{6,Vector{Float64}}} = nothing,
)
    n = length(signal)
    fill!(mask, false)
    n == 0 && return mask
    w = max(1, window_samples)
    if w >= n
        (maximum(signal) - minimum(signal) > threshold) && fill!(mask, true)
        return mask
    end
    w == 1 && return mask

    pref_max, pref_min, suff_max_prev, suff_min_prev, suff_max_curr, suff_min_curr =
        isnothing(ws) ? ntuple(_ -> Vector{Float64}(undef, w), 6) : ws

    n_blocks = cld(n, w)
    last_marked = 0
    thresh_f = Float64(threshold)

    @inbounds for b = 1:n_blocks
        b_start = (b - 1) * w + 1
        b_end = min(b_start + w - 1, n)
        b_len = b_end - b_start + 1

        p_max = p_min = Float64(signal[b_start])
        pref_max[1] = p_max
        pref_min[1] = p_min
        for k = 2:b_len
            s_val = Float64(signal[b_start+k-1])
            p_max = max(p_max, s_val)
            p_min = min(p_min, s_val)
            pref_max[k] = p_max
            pref_min[k] = p_min
        end

        s_max = s_min = Float64(signal[b_end])
        suff_max_curr[b_len] = s_max
        suff_min_curr[b_len] = s_min
        for k = (b_len-1):-1:1
            s_val = Float64(signal[b_start+k-1])
            s_max = max(s_max, s_val)
            s_min = min(s_min, s_val)
            suff_max_curr[k] = s_max
            suff_min_curr[k] = s_min
        end

        if b == 1
            if b_len == w && (pref_max[w] - pref_min[w]) > thresh_f
                mask[1:w] .= true
                last_marked = w
            end
        else
            for k = 2:w
                j_global = (b - 1) * w + k - 1
                j_global > n && break
                i_global = j_global - w + 1
                if (max(suff_max_prev[k], pref_max[k-1]) - min(suff_min_prev[k], pref_min[k-1])) > thresh_f
                    s_idx = max(last_marked + 1, i_global)
                    if s_idx <= j_global
                        mask[s_idx:j_global] .= true
                        last_marked = j_global
                    end
                end
            end
            if b_len == w && (pref_max[w] - pref_min[w]) > thresh_f
                s_idx = max(last_marked + 1, b_start)
                if s_idx <= b_end
                    mask[s_idx:b_end] .= true
                    last_marked = b_end
                end
            end
        end

        suff_max_prev, suff_max_curr = suff_max_curr, suff_max_prev
        suff_min_prev, suff_min_curr = suff_min_curr, suff_min_prev
    end

    return mask
end

_is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int) =
    _is_peak_to_peak!(falses(length(signal)), signal, threshold, window_samples)
