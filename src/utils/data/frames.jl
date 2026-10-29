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
    data = frames[1]
    # condition/condition_name/n_epochs live in the struct, not as columns
    legacy_cols = filter(col -> hasproperty(data, col), [:condition, :condition_name, :n_epochs])
    isempty(legacy_cols) || select!(data, Not(legacy_cols))
    return ErpData(t.file, condition, condition_name, data, copy(t.layout), t.sample_rate, copy(t.analysis_info), n_epochs)
end

function _derive_condition(t::TimeFreqData, condition::Int, condition_name::String, frames::Tuple{DataFrame,DataFrame}; n_epochs::Int = 0)
    return TimeFreqData(
        t.file,
        condition,
        condition_name,
        frames[1],
        frames[2],
        copy(t.layout),
        t.sample_rate,
        t.method,
        t.baseline,
        copy(t.analysis_info),
    )
end

"""Channels (non-metadata columns) of `dat` that condition arithmetic operates on."""
_condition_channels(dat::EegData) = setdiff(propertynames(dat.data), meta_labels(dat))
