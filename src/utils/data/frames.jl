# === CHANNEL-DATA CONTAINERS ===
#
# Uniform access to the DataFrame(s) holding channel values in each EEG data type, so that
# channel/condition operations can be written once for both EEG and time-frequency data.
#
# - EEG types (`ContinuousData`, `ErpData`, `EpochData`) hold one container: `data`.
# - TF types (`TimeFreqData`, `TimeFreqEpochData`) hold two: `data_power` and `data_phase`,
#   which must always be transformed together to stay in sync.
#
# A "container" is a `DataFrame` (single-DataFrame types) or a `Vector{DataFrame}` (epoch types).

"""
    _data_frames(dat::EegData) -> Tuple

Return every channel-data container of `dat` as a tuple: `(dat.data,)` for EEG data and
`(dat.data_power, dat.data_phase)` for time-frequency data.
"""
_data_frames(dat::EegData) = (dat.data,)
_data_frames(dat::Union{TimeFreqData,TimeFreqEpochData}) = (dat.data_power, dat.data_phase)

"""
    _set_data_frames!(dat::EegData, frames::Tuple) -> dat

Replace the channel-data containers of `dat` (inverse of [`_data_frames`](@ref)).
"""
function _set_data_frames!(dat::EegData, frames::Tuple{Any})
    dat.data = frames[1]
    return dat
end
function _set_data_frames!(dat::Union{TimeFreqData,TimeFreqEpochData}, frames::Tuple{Any,Any})
    dat.data_power, dat.data_phase = frames
    return dat
end

"""Apply `f(df)` to each DataFrame of a container (`DataFrame` or `Vector{DataFrame}`)."""
_foreach_df(f, df::DataFrame) = (f(df); nothing)
_foreach_df(f, dfs::AbstractVector{<:DataFrame}) = (foreach(f, dfs); nothing)

"""Map `f(df) -> DataFrame` over a container, preserving its shape."""
_map_df(f, df::DataFrame) = f(df)
_map_df(f, dfs::AbstractVector{<:DataFrame}) = DataFrame[f(df) for df in dfs]

"""Short tag used in log messages to distinguish TF from EEG operations."""
_kind_tag(::Type{<:Union{TimeFreqData,TimeFreqEpochData}}) = "TF "
_kind_tag(::Type) = ""

# === DERIVED CONDITIONS ===

"""Data types that represent a single (averaged) condition and support condition arithmetic."""
const ConditionAverageData = Union{ErpData,TimeFreqData}

"""
    _derive_condition(template, condition, condition_name, frames; n_epochs) -> typeof(template)

Build a new condition object of the same type as `template` holding the given channel-data
`frames` (as returned by [`_data_frames`](@ref)). File, sample rate and (copies of) layout and
analysis info are taken from `template`. `n_epochs` is only used by `ErpData`.
"""
function _derive_condition(t::ErpData, condition::Int, condition_name::String, frames::Tuple{DataFrame}; n_epochs::Int)
    return ErpData(t.file, condition, condition_name, frames[1], copy(t.layout), t.sample_rate, copy(t.analysis_info), n_epochs)
end

function _derive_condition(t::TimeFreqData, condition::Int, condition_name::String, frames::Tuple{DataFrame,DataFrame}; n_epochs::Int = 0)
    phase_df = frames[2]
    if startswith(condition_name, "difference_")
        phase_df = copy(phase_df)
        for ch in channel_labels(t)
            if hasproperty(phase_df, ch)
                phase_df[!, ch] = rem2pi.(phase_df[!, ch], RoundNearest)
            end
        end
    end
    return TimeFreqData(
        t.file,
        condition,
        condition_name,
        frames[1],
        phase_df,
        copy(t.layout),
        t.sample_rate,
        t.method,
        t.baseline,
        copy(t.analysis_info),
    )
end

"""Channels (electrode columns only, excluding metadata and extra columns) that condition arithmetic operates on."""
_condition_channels(dat::EegData) = channel_labels(dat)

"""Wrap phase differences in radians onto `[-π, π]` using circular arithmetic."""
_circular_phase_diff(p1::Real, p2::Real) = rem2pi(Float64(p1) - Float64(p2), RoundNearest)
_circular_phase_diff(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}) = rem2pi.(Float64.(p1) .- Float64.(p2), RoundNearest)

"""Compute circular mean of phase angles in radians via `atan(mean(sin), mean(cos))`."""
function _circular_phase_mean(phases::AbstractVector{<:AbstractVector{<:Real}})::Vector{Float64}
    isempty(phases) && return Float64[]
    n_rows = length(phases[1])
    n_conds = length(phases)
    sin_sum = zeros(Float64, n_rows)
    cos_sum = zeros(Float64, n_rows)
    for p in phases
        @inbounds @simd for i = 1:n_rows
            val = Float64(p[i])
            sin_sum[i] += sin(val)
            cos_sum[i] += cos(val)
        end
    end
    inv_n = 1.0 / n_conds
    return [atan(sin_sum[i] * inv_n, cos_sum[i] * inv_n) for i = 1:n_rows]
end

"""Compute condition-difference channel-data frames, using circular difference for `TimeFreqData` phase."""
function _difference_frames(dat1::ErpData, dat2::ErpData)
    diff_df = copy(dat1.data)
    for ch in _condition_channels(dat1)
        hasproperty(dat2.data, ch) && (diff_df[!, ch] .-= dat2.data[!, ch])
    end
    return (diff_df,)
end

function _difference_frames(dat1::TimeFreqData, dat2::TimeFreqData)
    diff_power = copy(dat1.data_power)
    diff_phase = copy(dat1.data_phase)
    for ch in _condition_channels(dat1)
        if hasproperty(dat2.data_power, ch)
            diff_power[!, ch] .-= dat2.data_power[!, ch]
        end
        if hasproperty(dat2.data_phase, ch)
            diff_phase[!, ch] = _circular_phase_diff(dat1.data_phase[!, ch], dat2.data_phase[!, ch])
        end
    end
    return (diff_power, diff_phase)
end

"""Compute condition-average channel-data frames, using circular mean for `TimeFreqData` phase."""
function _average_frames(dats::AbstractVector{<:ErpData})
    avg_df = copy(dats[1].data)
    n = length(dats)
    for ch in _condition_channels(dats[1])
        avg_df[!, ch] .= 0.0
        for d in dats
            avg_df[!, ch] .+= d.data[!, ch]
        end
        avg_df[!, ch] ./= n
    end
    return (avg_df,)
end

function _average_frames(dats::AbstractVector{<:TimeFreqData})
    avg_power = copy(dats[1].data_power)
    avg_phase = copy(dats[1].data_phase)
    n = length(dats)
    for ch in _condition_channels(dats[1])
        avg_power[!, ch] .= 0.0
        for d in dats
            avg_power[!, ch] .+= d.data_power[!, ch]
        end
        avg_power[!, ch] ./= n
        avg_phase[!, ch] = _circular_phase_mean([d.data_phase[!, ch] for d in dats])
    end
    return (avg_power, avg_phase)
end
