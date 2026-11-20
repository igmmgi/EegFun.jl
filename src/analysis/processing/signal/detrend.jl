function _detrend_precompute(time::AbstractVector, idx_interval::Tuple{Int,Int})
    start_idx, end_idx = idx_interval
    k = end_idx - start_idx + 1
    inv_k = 1.0 / k
    sum_t = 0.0
    @inbounds @simd for i = start_idx:end_idx
        sum_t += Float64(time[i])
    end
    mean_t = sum_t * inv_k
    stt = 0.0
    @inbounds @simd for i = start_idx:end_idx
        dt = Float64(time[i]) - mean_t
        stt = muladd(dt, dt, stt)
    end
    return mean_t, inv_k, stt
end

function _apply_detrend_df!(
    df::DataFrame,
    ch_indices::Vector{Int},
    idx_interval::Tuple{Int,Int},
    type::Symbol,
    mean_t::Float64,
    inv_k::Float64,
    stt::Float64,
)
    start_idx, end_idx = idx_interval
    if type !== :linear || stt < 1e-15
        _apply_baseline!(df, ch_indices, idx_interval)
        return nothing
    end

    inv_stt = 1.0 / stt
    t_col = df[!, :time]
    for col_idx in ch_indices
        col = df[!, col_idx]
        if col isa Vector{Float64} && t_col isa Vector{Float64}
            col_f64 = col::Vector{Float64}
            t_f64 = t_col::Vector{Float64}
            sum_y = 0.0
            sty = 0.0
            @inbounds @simd for i = start_idx:end_idx
                yi = col_f64[i]
                sum_y += yi
                sty = muladd(t_f64[i] - mean_t, yi, sty)
            end
            mean_y = sum_y * inv_k
            m = sty * inv_stt
            b = mean_y - m * mean_t
            @inbounds @simd for i in eachindex(col_f64, t_f64)
                col_f64[i] -= muladd(m, t_f64[i], b)
            end
        else
            sum_y = 0.0
            sty = 0.0
            @inbounds for i = start_idx:end_idx
                yi = Float64(col[i])
                sum_y += yi
                sty = muladd(Float64(t_col[i]) - mean_t, yi, sty)
            end
            mean_y = sum_y * inv_k
            m = sty * inv_stt
            b = mean_y - m * mean_t
            df[!, col_idx] = Float64[Float64(col[i]) - muladd(m, Float64(t_col[i]), b) for i in eachindex(col, t_col)]
        end
    end
    return nothing
end

function _apply_detrend!(dat::DataFrame, channels::AbstractVector{Symbol}, idx_interval::Tuple{Int,Int}, type::Symbol)
    ch_indices = Int[columnindex(dat, ch) for ch in channels]
    mean_t, inv_k, stt = _detrend_precompute(dat[!, :time], idx_interval)
    _apply_detrend_df!(dat, ch_indices, idx_interval, type, mean_t, inv_k, stt)
    return nothing
end

function _apply_detrend!(dat::AbstractVector{DataFrame}, channels::AbstractVector{Symbol}, idx_interval::Tuple{Int,Int}, type::Symbol)
    isempty(dat) && return nothing
    ch_indices = Int[columnindex(dat[1], ch) for ch in channels]
    mean_t, inv_k, stt = _detrend_precompute(dat[1][!, :time], idx_interval)
    if length(dat) > 1 && Threads.nthreads() > 1 && length(dat) * length(ch_indices) * nrow(dat[1]) >= 65_536
        Threads.@threads for df in dat
            _apply_detrend_df!(df, ch_indices, idx_interval, type, mean_t, inv_k, stt)
        end
    else
        for df in dat
            _apply_detrend_df!(df, ch_indices, idx_interval, type, mean_t, inv_k, stt)
        end
    end
    return nothing
end

"""
    detrend!(dat::Union{ContinuousData,ErpData,EpochData}; channel_selection = channels(), interval_selection = times(), type = :linear)
    detrend!(dat::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}; channel_selection = channels(), interval_selection = times(), type = :linear)

Remove the least-squares linear trend (`type = :linear`) or constant offset (`type = :constant`) estimated over
`interval_selection` from each selected channel in-place. Passing `interval_selection = nothing` leaves `dat` unchanged.
"""
function detrend!(dat::_TimeDomainEegData; channel_selection = channels(), interval_selection::Interval = times(), type::Symbol = :linear)
    isnothing(interval_selection) && return dat
    type in (:linear, :constant, :mean) || @minimal_error("Unknown detrend type :$type. Use :linear or :constant.")
    dat isa EpochData && isempty(dat.data) && return dat
    n_samples(dat) == 0 && return dat

    idx_interval = _validate_baseline_interval(dat, interval_selection)
    if type === :linear && (idx_interval[2] - idx_interval[1] + 1) < 2
        @minimal_error("Linear detrending requires at least 2 samples in interval_selection")
    end

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for detrending"
        return dat
    end

    _apply_detrend!(dat.data, selected_channels, idx_interval, type)
    return dat
end

function detrend!(
    dat::AbstractVector{<:_TimeDomainEegData};
    channel_selection = channels(),
    interval_selection::Interval = times(),
    type::Symbol = :linear,
)
    isnothing(interval_selection) && return dat
    for d in dat
        detrend!(d; channel_selection = channel_selection, interval_selection = interval_selection, type = type)
    end
    return dat
end

"""
    detrend(dat::Union{ContinuousData,ErpData,EpochData}; channel_selection = channels(), interval_selection = times(), type = :linear)
    detrend(dat::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}; channel_selection = channels(), interval_selection = times(), type = :linear)

Non-mutating version of [`detrend!`](@ref).
"""
detrend(dat::_TimeDomainEegData; kwargs...) = detrend!(copy(dat); kwargs...)
detrend(dat::AbstractVector{<:_TimeDomainEegData}; kwargs...) = [detrend(d; kwargs...) for d in dat]
