_data_frames(dat::EegData) = (dat.data,)
_data_frames(dat::Union{TimeFreqData,TimeFreqEpochData}) = (dat.data_power, dat.data_phase)

function _set_data_frames!(dat::EegData, frames::Tuple{Any})
    dat.data = frames[1]
    return dat
end
function _set_data_frames!(dat::Union{TimeFreqData,TimeFreqEpochData}, frames::Tuple{Any,Any})
    dat.data_power, dat.data_phase = frames
    return dat
end

_foreach_df(f, df::DataFrame) = (f(df); nothing)
_foreach_df(f, dfs::AbstractVector{<:DataFrame}) = (foreach(f, dfs); nothing)

_map_df(f, df::DataFrame) = f(df)
_map_df(f, dfs::AbstractVector{<:DataFrame}) = DataFrame[f(df) for df in dfs]

_kind_tag(::Type{<:Union{TimeFreqData,TimeFreqEpochData}}) = "TF "
_kind_tag(::Type) = ""

const ConditionAverageData = Union{ErpData,TimeFreqData}

@inline function _accum_col!(dst::Vector{Float64}, src::Vector{Float64})
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] += src[i]
    end
end
@inline function _accum_col!(dst::AbstractVector, src::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] += src[i]
    end
end

@inline function _rem2pi_col!(col::AbstractVector)
    @inbounds @simd for i in eachindex(col)
        col[i] = rem2pi(col[i], RoundNearest)
    end
end

function _derive_condition(t::ErpData, condition::Int, condition_name::String, frames::Tuple{DataFrame}; n_epochs::Int)
    return ErpData(t.file, condition, condition_name, frames[1], copy(t.layout), t.sample_rate, copy(t.analysis_info), n_epochs)
end

function _derive_condition(t::TimeFreqData, condition::Int, condition_name::String, frames::Tuple{DataFrame,DataFrame}; n_epochs::Int = 0)
    phase_df = frames[2]
    if startswith(condition_name, "difference_")
        phase_df = copy(phase_df)
        for ch in channel_labels(t)
            hasproperty(phase_df, ch) && _rem2pi_col!(phase_df[!, ch])
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

_condition_channels(dat::EegData) = channel_labels(dat)
