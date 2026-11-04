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

@inline function _sub_col!(dst::Vector{Float64}, src::Vector{Float64})
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] -= src[i]
    end
end
@inline function _sub_col!(dst::AbstractVector, src::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] -= src[i]
    end
end

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

@inline function _scale_col!(dst::AbstractVector, scale::Float64)
    @inbounds @simd for i in eachindex(dst)
        dst[i] *= scale
    end
end

@inline function _rem2pi_col!(col::AbstractVector)
    @inbounds @simd for i in eachindex(col)
        col[i] = rem2pi(col[i], RoundNearest)
    end
end

@inline function _rem2pi_diff_col!(dst::Vector{Float64}, src::Vector{Float64})
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] = rem2pi(dst[i] - src[i], RoundNearest)
    end
end
@inline function _rem2pi_diff_col!(dst::AbstractVector, src::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] = rem2pi(Float64(dst[i]) - Float64(src[i]), RoundNearest)
    end
end

@inline function _accum_sincos!(sin_sum::Vector{Float64}, cos_sum::Vector{Float64}, p::Vector{Float64})
    @inbounds for i in eachindex(sin_sum, cos_sum, p)
        s, c = sincos(p[i])
        sin_sum[i] += s
        cos_sum[i] += c
    end
end
@inline function _accum_sincos!(sin_sum::Vector{Float64}, cos_sum::Vector{Float64}, p::AbstractVector)
    @inbounds for i in eachindex(sin_sum, cos_sum, p)
        s, c = sincos(Float64(p[i]))
        sin_sum[i] += s
        cos_sum[i] += c
    end
end

@inline function _write_atan!(out::AbstractVector, sin_sum::Vector{Float64}, cos_sum::Vector{Float64}, inv_n::Float64)
    @inbounds for i in eachindex(out, sin_sum, cos_sum)
        out[i] = atan(sin_sum[i] * inv_n, cos_sum[i] * inv_n)
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

function _circular_phase_mean!(
    out::AbstractVector,
    sin_sum::Vector{Float64},
    cos_sum::Vector{Float64},
    dats::AbstractVector{<:TimeFreqData},
    ch::Symbol,
)
    fill!(sin_sum, 0.0)
    fill!(cos_sum, 0.0)
    for d in dats
        _accum_sincos!(sin_sum, cos_sum, d.data_phase[!, ch])
    end
    _write_atan!(out, sin_sum, cos_sum, 1.0 / length(dats))
    return out
end

function _difference_frames(dat1::ErpData, dat2::ErpData)
    diff_df = copy(dat1.data)
    for ch in _condition_channels(dat1)
        hasproperty(dat2.data, ch) && _sub_col!(diff_df[!, ch], dat2.data[!, ch])
    end
    return (diff_df,)
end

function _difference_frames(dat1::TimeFreqData, dat2::TimeFreqData)
    diff_power = copy(dat1.data_power)
    diff_phase = copy(dat1.data_phase)
    for ch in _condition_channels(dat1)
        hasproperty(dat2.data_power, ch) && _sub_col!(diff_power[!, ch], dat2.data_power[!, ch])
        hasproperty(dat2.data_phase, ch) && _rem2pi_diff_col!(diff_phase[!, ch], dat2.data_phase[!, ch])
    end
    return (diff_power, diff_phase)
end

function _average_frames(dats::AbstractVector{<:ErpData})
    avg_df = copy(dats[1].data)
    inv_n = 1.0 / length(dats)
    for ch in _condition_channels(dats[1])
        dst = avg_df[!, ch]
        for k = 2:length(dats)
            _accum_col!(dst, dats[k].data[!, ch])
        end
        _scale_col!(dst, inv_n)
    end
    return (avg_df,)
end

function _average_frames(dats::AbstractVector{<:TimeFreqData})
    avg_power = copy(dats[1].data_power)
    avg_phase = copy(dats[1].data_phase)
    inv_n = 1.0 / length(dats)
    n_rows = nrow(avg_phase)
    sin_sum = Vector{Float64}(undef, n_rows)
    cos_sum = Vector{Float64}(undef, n_rows)
    for ch in _condition_channels(dats[1])
        dst = avg_power[!, ch]
        for k = 2:length(dats)
            _accum_col!(dst, dats[k].data_power[!, ch])
        end
        _scale_col!(dst, inv_n)
        _circular_phase_mean!(avg_phase[!, ch], sin_sum, cos_sum, dats, ch)
    end
    return (avg_power, avg_phase)
end
