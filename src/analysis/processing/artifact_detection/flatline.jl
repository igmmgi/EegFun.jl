function _is_flatline!(mask::AbstractVector{Bool}, signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)
    n = length(signal)
    fill!(mask, false)
    n == 0 && return mask
    w = max(1, window_samples)
    if w >= n
        std(signal) < threshold && fill!(mask, true)
        return mask
    end
    w <= 1 && return mask

    w_view = view(signal, 1:w)
    μ = mean(w_view)
    M2 = var(w_view) * (w - 1)

    inv_w = 1.0 / w
    inv_w_minus_1 = 1.0 / (w - 1)
    thresh_sq = Float64(threshold)^2

    last_marked = 0
    if (M2 * inv_w_minus_1) < thresh_sq
        mask[1:w] .= true
        last_marked = w
    end

    @inbounds for i = 2:(n-w+1)
        x_out = Float64(signal[i-1])
        x_in = Float64(signal[i+w-1])

        if (i & 255) == 0
            cur_view = view(signal, i:(i+w-1))
            μ = mean(cur_view)
            M2 = var(cur_view) * (w - 1)
        else
            d = x_in - x_out
            old_μ = μ
            μ = old_μ + d * inv_w
            M2 = max(0.0, M2 + d * ((x_in - μ) + (x_out - old_μ)))
        end

        if (M2 * inv_w_minus_1) < thresh_sq
            s_idx = max(last_marked + 1, i)
            e_idx = i + w - 1
            if s_idx <= e_idx
                mask[s_idx:e_idx] .= true
                last_marked = e_idx
            end
        end
    end
    return mask
end

_is_flatline(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int) =
    _is_flatline!(falses(length(signal)), signal, threshold, window_samples)

"""
    is_flatline!(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                 channel_selection = channels(), sample_selection = samples(),
                 interval_selection = times(), mode::Symbol = :combined, channel_out = nothing)
    is_flatline!(dat::MultiDataFrameEeg, threshold::Real, window_size::Real; kwargs...)
    is_flatline!(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...)

Detect flatline artifacts in sliding windows of `window_size` seconds where `std < threshold`,
storing boolean flag columns in-place (`:is_flatline` or `:is_flatline_<ch>` by default).
"""
function is_flatline!(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    combined_col = something(channel_out, :is_flatline)
    sep_fn = ch -> Symbol("$(combined_col)_$(ch)")
    return _flag_artifacts!(
        (buf, sig) -> _is_flatline!(buf, sig, threshold, window_samples),
        dat,
        combined_col,
        sep_fn;
        channel_selection,
        sample_selection,
        interval_selection,
        mode,
    )
end

function is_flatline!(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    combined_col = something(channel_out, :is_flatline)
    sep_fn = ch -> Symbol("$(combined_col)_$(ch)")
    return _flag_artifacts!(
        (buf, sig) -> _is_flatline!(buf, sig, threshold, window_samples),
        dat,
        combined_col,
        sep_fn;
        channel_selection,
        sample_selection,
        interval_selection,
        epoch_selection,
        mode,
    )
end

function is_flatline!(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...)
    for d in dat
        is_flatline!(d, threshold, window_size; kwargs...)
    end
    return dat
end

"""
    is_flatline(dat, threshold::Real, window_size::Real; channel_selection = channels(),
                sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined)

Detect flatline artifacts without modifying `dat`. Returns a `Vector{Bool}` (or per-epoch vector)
in `:combined` mode, or a `DataFrame` (or per-epoch vector) in `:separate` mode.
"""
function is_flatline(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    return _detect_artifacts(
        (buf, sig) -> _is_flatline!(buf, sig, threshold, window_samples),
        dat,
        ch -> Symbol("is_flatline_$(ch)");
        channel_selection,
        sample_selection,
        interval_selection,
        mode,
    )
end

function is_flatline(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    return _detect_artifacts(
        (buf, sig) -> _is_flatline!(buf, sig, threshold, window_samples),
        dat,
        ch -> Symbol("is_flatline_$(ch)");
        channel_selection,
        sample_selection,
        interval_selection,
        epoch_selection,
        mode,
    )
end

is_flatline(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...) =
    [is_flatline(d, threshold, window_size; kwargs...) for d in dat]

"""
    n_flatline(dat, threshold::Real, window_size::Real; channel_selection = channels(),
               sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined)

Count flatline samples across selected channels, returning a total `Int` (`:combined`) or a per-channel `DataFrame` (`:separate`).
"""
function n_flatline(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    return _count_artifacts(
        (buf, sig) -> _is_flatline!(buf, sig, threshold, window_samples),
        dat,
        :n_flatline;
        channel_selection,
        sample_selection,
        interval_selection,
        mode,
    )
end

function n_flatline(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold, window_size)
    window_samples = round(Int, window_size * dat.sample_rate)
    return _count_artifacts(
        (buf, sig) -> _is_flatline!(buf, sig, threshold, window_samples),
        dat,
        :n_flatline;
        channel_selection,
        sample_selection,
        interval_selection,
        epoch_selection,
        mode,
    )
end

n_flatline(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...) =
    [n_flatline(d, threshold, window_size; kwargs...) for d in dat]
