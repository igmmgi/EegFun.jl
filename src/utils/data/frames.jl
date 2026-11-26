# Conversion between EegFun data containers and tabular `DataFrame` representations.

######################################
######### Types & Constants ##########
######################################

const ConditionAverageData = Union{ErpData,TimeFreqData}

######################################
######### Internal Functions #########
######################################

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

@inline function _copy_f64!(dst::AbstractVector{Float64}, src::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] = Float64(src[i])
    end
    return dst
end

@inline function _accum_col!(dst::AbstractVector{Float64}, src::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] += Float64(src[i])
    end
    return dst
end

@inline function _accum_col!(dst::AbstractVector, src::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] += src[i]
    end
    return dst
end

@inline function _add_scaled_col!(dst::Vector{Float64}, src::AbstractVector, scale::Float64)
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] = muladd(Float64(src[i]), scale, dst[i])
    end
    return dst
end

@inline function _sub_col!(dst::Vector{Float64}, src1::AbstractVector, src2::AbstractVector)
    @inbounds @simd for i in eachindex(dst, src1, src2)
        dst[i] = Float64(src1[i]) - Float64(src2[i])
    end
    return dst
end

function _fill_colmeans_idx!(out::Vector{Float64}, df::DataFrame, idxs::AbstractVector{<:Integer})::Vector{Float64}
    isempty(idxs) && return fill!(out, 0.0)
    _copy_f64!(out, df[!, idxs[1]])
    for m = 2:length(idxs)
        _accum_col!(out, df[!, idxs[m]])
    end
    if length(idxs) > 1
        k_f64 = Float64(length(idxs))
        @inbounds @simd for i in eachindex(out)
            out[i] /= k_f64
        end
    end
    return out
end

@inline function _rem2pi_col!(col::AbstractVector)
    @inbounds @simd for i in eachindex(col)
        col[i] = rem2pi(col[i], RoundNearest)
    end
end

function _sync_condition_cols!(df::DataFrame, condition::Int, condition_name::String)
    hasproperty(df, :condition) && (df[!, :condition] = fill(condition, nrow(df)))
    hasproperty(df, :condition_name) && (df[!, :condition_name] = fill(condition_name, nrow(df)))
    return df
end

function _derive_condition(t::ErpData, condition::Int, condition_name::String, frames::Tuple{DataFrame}; n_epochs::Int)
    _sync_condition_cols!(frames[1], condition, condition_name)
    return ErpData(t.file, condition, condition_name, frames[1], copy(t.layout), t.sample_rate, copy(t.analysis_info), n_epochs)
end

function _derive_condition(t::TimeFreqData, condition::Int, condition_name::String, frames::Tuple{DataFrame,DataFrame}; n_epochs::Int = 0)
    power_df = _sync_condition_cols!(frames[1], condition, condition_name)
    phase_df = frames[2]
    if startswith(condition_name, "difference_")
        phase_df = copy(phase_df)
        for ch in _condition_channels(t)
            hasproperty(phase_df, ch) && _rem2pi_col!(phase_df[!, ch])
        end
    end
    _sync_condition_cols!(phase_df, condition, condition_name)
    return TimeFreqData(
        t.file,
        condition,
        condition_name,
        power_df,
        phase_df,
        copy(t.layout),
        t.sample_rate,
        t.method,
        t.baseline,
        copy(t.analysis_info),
    )
end

function _condition_channels(dat::EegData)::Vector{Symbol}
    chs = copy(channel_labels(dat))
    df = _first_df(dat)
    isnothing(df) && return chs
    for ex in extra_labels(dat)
        if hasproperty(df, ex) && eltype(df[!, ex]) <: AbstractFloat && ex ∉ chs
            push!(chs, ex)
        end
    end
    return chs
end
