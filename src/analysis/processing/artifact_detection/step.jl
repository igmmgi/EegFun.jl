# Sample-to-sample voltage step discontinuity detection.

######################################
########### User Functions ###########
######################################

"""
    is_step_value(dat, threshold::Real; channel_selection = channels(), sample_selection = samples(),
                  interval_selection = times(), mode::Symbol = :combined)

Detect step artifacts (`|x[i] - x[i-1]| > threshold`) without modifying `dat`. Returns a `Vector{Bool}`
(or per-epoch vector) in `:combined` mode, or a `DataFrame` (or per-epoch vector) in `:separate` mode.
"""
function is_step_value(dat::_DetectableEeg, threshold::Real; mode::Symbol = :combined, kwargs...)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    sep_fn = ch -> Symbol("is_step_value_$(ch)_$(threshold)")
    return _detect_artifacts(() -> (buf, sig) -> _is_step_value!(buf, sig, thresh_f), dat, sep_fn; mode, kwargs...)
end

is_step_value(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...) = [is_step_value(d, threshold; kwargs...) for d in dat]

"""
    is_step_value!(dat::SingleDataFrameEeg, threshold::Real; channel_selection = channels(),
                   sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined, channel_out = nothing)
    is_step_value!(dat::MultiDataFrameEeg, threshold::Real; channel_selection = channels(),
                   sample_selection = samples(), interval_selection = times(), epoch_selection = epochs(),
                   mode::Symbol = :combined, channel_out = nothing)
    is_step_value!(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...)

Detect step artifacts (`|x[i] - x[i-1]| > threshold`) and store boolean flag columns in-place
(`:is_step_value_<threshold>` in `:combined` mode, or `:is_step_value_<ch>_<threshold>` in `:separate` mode).
"""
function is_step_value!(
    dat::_DetectableEeg,
    threshold::Real;
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
    kwargs...,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    combined_col = something(channel_out, Symbol("is_step_value_$(threshold)"))
    sep_fn = ch -> Symbol("is_step_value_$(ch)_$(threshold)")
    return _flag_artifacts!(() -> (buf, sig) -> _is_step_value!(buf, sig, thresh_f), dat, combined_col, sep_fn; mode, kwargs...)
end

is_step_value!(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...) =
    (foreach(d -> is_step_value!(d, threshold; kwargs...), dat); dat)

"""
    n_step_value(dat, threshold::Real; channel_selection = channels(), sample_selection = samples(),
                 interval_selection = times(), mode::Symbol = :combined)

Count step artifacts (`|x[i] - x[i-1]| > threshold`) across selected channels, returning a total `Int`
(`:combined`) or a per-channel `DataFrame` (`:separate`).
"""
function n_step_value(dat::_DetectableEeg, threshold::Real; mode::Symbol = :combined, kwargs...)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    return _count_artifacts(() -> (buf, sig) -> _is_step_value!(buf, sig, thresh_f), dat, :n_step; mode, kwargs...)
end

n_step_value(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...) = [n_step_value(d, threshold; kwargs...) for d in dat]

######################################
######### Internal Functions #########
######################################

function _is_step_value!(mask::AbstractVector{Bool}, signal::AbstractVector{<:Real}, threshold::Real)
    n = length(signal)
    length(mask) == n || @minimal_error("Mask and signal must have the same length (mask: $(length(mask)), signal: $n)")
    n == 0 && return nothing
    mask[1] = false
    @views mask[2:n] .= abs.(signal[2:n] .- signal[1:(n-1)]) .> threshold
    return nothing
end

function _is_step_value(signal::AbstractVector{<:Real}, threshold::Real)
    mask = Vector{Bool}(undef, length(signal))
    _is_step_value!(mask, signal, threshold)
    return mask
end
