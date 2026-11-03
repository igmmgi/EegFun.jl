_is_extreme_value(signal::AbstractVector{<:Real}, threshold::Real) = abs.(signal) .> threshold

function _is_extreme_value!(mask::AbstractVector{Bool}, signal::AbstractVector{<:Real}, threshold::Real)
    length(mask) == length(signal) ||
        @minimal_error("Mask and signal must have the same length (mask: $(length(mask)), signal: $(length(signal)))")
    mask .= abs.(signal) .> threshold
    return nothing
end

function _detect_extreme_values(dat::SingleDataFrameEeg, threshold::Real; channel_selection = channels(), sample_selection = samples())
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for extreme value detection")
    selected_samples = get_selected_samples(dat, sample_selection)
    isempty(selected_samples) && @minimal_error("No samples selected for extreme value detection")

    n_samp = n_samples(dat)
    sample_mask = falses(n_samp)
    sample_mask[selected_samples] .= true
    thresh_f = Float64(threshold)
    return Dict{Symbol,Vector{Bool}}(
        ch => Vector{Bool}(_is_extreme_value(dat.data[!, ch], thresh_f) .& sample_mask) for ch in selected_channels
    )
end

_n_extreme_value(df::DataFrame, channels::AbstractVector{Symbol}, threshold::Real) =
    [count(x -> abs(x) > threshold, df[!, ch]) for ch in channels]

"""
    is_extreme_value!(dat::SingleDataFrameEeg, threshold::Real; channel_selection = channels(),
                      sample_selection = samples(), interval_selection = times(), mode::Symbol = :combined, channel_out = nothing)
    is_extreme_value!(dat::MultiDataFrameEeg, threshold::Real; channel_selection = channels(),
                      sample_selection = samples(), interval_selection = times(), epoch_selection = epochs(),
                      mode::Symbol = :combined, channel_out = nothing)
    is_extreme_value!(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...)

Detect extreme values (`|x| > threshold`) and store boolean flag columns in-place (`:is_extreme_value_<threshold>`
in `:combined` mode, or `:is_extreme_value_<ch>_<threshold>` in `:separate` mode).
"""
function is_extreme_value!(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    combined_col = something(channel_out, Symbol("is_extreme_value_$(threshold)"))
    sep_fn = ch -> Symbol("is_extreme_value_$(ch)_$(threshold)")
    return _flag_artifacts!(
        (buf, sig) -> _is_extreme_value!(buf, sig, thresh_f),
        dat,
        combined_col,
        sep_fn;
        channel_selection,
        sample_selection,
        interval_selection,
        mode,
    )
end

function is_extreme_value!(
    dat::MultiDataFrameEeg,
    threshold::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    combined_col = something(channel_out, Symbol("is_extreme_value_$(threshold)"))
    sep_fn = ch -> Symbol("is_extreme_value_$(ch)_$(threshold)")
    return _flag_artifacts!(
        (buf, sig) -> _is_extreme_value!(buf, sig, thresh_f),
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

function is_extreme_value!(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...)
    for d in dat
        is_extreme_value!(d, threshold; kwargs...)
    end
    return dat
end

"""
    is_extreme_value(dat, threshold::Real; channel_selection = channels(), sample_selection = samples(),
                     interval_selection = times(), mode::Symbol = :combined)

Detect extreme values (`|x| > threshold`) without modifying `dat`. Returns a `Vector{Bool}` (or per-epoch vector)
in `:combined` mode, or a `DataFrame` (or per-epoch vector of `DataFrame`s) in `:separate` mode.
"""
function is_extreme_value(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    sep_fn = ch -> Symbol("is_extreme_value_$(ch)_$(threshold)")
    return _detect_artifacts(
        (buf, sig) -> _is_extreme_value!(buf, sig, thresh_f),
        dat,
        sep_fn;
        channel_selection,
        sample_selection,
        interval_selection,
        mode,
    )
end

function is_extreme_value(
    dat::MultiDataFrameEeg,
    threshold::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    sep_fn = ch -> Symbol("is_extreme_value_$(ch)_$(threshold)")
    return _detect_artifacts(
        (buf, sig) -> _is_extreme_value!(buf, sig, thresh_f),
        dat,
        sep_fn;
        channel_selection,
        sample_selection,
        interval_selection,
        epoch_selection,
        mode,
    )
end

is_extreme_value(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...) = [is_extreme_value(d, threshold; kwargs...) for d in dat]

"""
    n_extreme_value(dat, threshold::Real; channel_selection = channels(), sample_selection = samples(),
                    interval_selection = times(), mode::Symbol = :combined)

Count extreme values (`|x| > threshold`) across selected channels, returning a total `Int` (`:combined`)
or a per-channel `DataFrame` (`:separate`).
"""
function n_extreme_value(
    dat::SingleDataFrameEeg,
    threshold::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    return _count_artifacts(
        (buf, sig) -> _is_extreme_value!(buf, sig, thresh_f),
        dat,
        :n_extreme;
        channel_selection,
        sample_selection,
        interval_selection,
        mode,
    )
end

function n_extreme_value(
    dat::MultiDataFrameEeg,
    threshold::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
)
    _validate_detector_args(mode, threshold)
    thresh_f = Float64(threshold)
    return _count_artifacts(
        (buf, sig) -> _is_extreme_value!(buf, sig, thresh_f),
        dat,
        :n_extreme;
        channel_selection,
        sample_selection,
        interval_selection,
        epoch_selection,
        mode,
    )
end

n_extreme_value(dat::AbstractVector{<:EpochData}, threshold::Real; kwargs...) = [n_extreme_value(d, threshold; kwargs...) for d in dat]

"""
    n_values(df::DataFrame, column::Symbol)
    n_values(dat::SingleDataFrameEeg, column::Symbol)
    n_values(dat::MultiDataFrameEeg, column::Symbol)
    n_values(dat::AbstractVector{<:EegData}, column::Symbol)

Count the number of `true` values in boolean `column` (returning `0` if `column` is absent).
"""
n_values(df::DataFrame, column::Symbol) = hasproperty(df, column) ? count(df[!, column]) : 0
n_values(dat::SingleDataFrameEeg, column::Symbol) = n_values(dat.data, column)
n_values(dat::MultiDataFrameEeg, column::Symbol) = sum(n_values(epoch, column) for epoch in dat; init = 0)
n_values(dat::AbstractVector{<:EegData}, column::Symbol) = sum(n_values(d, column) for d in dat; init = 0)
