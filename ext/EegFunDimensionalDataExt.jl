module EegFunDimensionalDataExt

using EegFun
using DimensionalData

function _try_extract_standard_times(df_freq::AbstractVector, df_time::AbstractVector, freqs::Vector{Float64})
    n_f = length(freqs)
    n_rows = length(df_time)
    (n_f == 0 || n_rows == 0 || n_rows % n_f != 0 || length(df_freq) != n_rows) && return nothing
    n_t = n_rows ÷ n_f
    times = Vector{Float64}(undef, n_t)
    idx = 0
    @inbounds for ti = 1:n_t
        t_val = Float64(df_time[idx+1])
        (ti > 1 && !(t_val > times[ti-1])) && return nothing
        times[ti] = t_val
        for fi = 1:n_f
            idx += 1
            (Float64(df_freq[idx]) == freqs[fi] && Float64(df_time[idx]) == t_val) || return nothing
        end
    end
    return times
end

function _is_standard_tf_grid(df_freq::AbstractVector, df_time::AbstractVector, freqs::Vector{Float64}, times::Vector{Float64})
    n_f = length(freqs)
    n_t = length(times)
    n_rows = n_f * n_t
    (length(df_freq) == n_rows && length(df_time) == n_rows) || return false
    idx = 0
    @inbounds for ti = 1:n_t
        t_val = times[ti]
        for fi = 1:n_f
            idx += 1
            (Float64(df_freq[idx]) == freqs[fi] && Float64(df_time[idx]) == t_val) || return false
        end
    end
    return true
end

function _map_tf_rows(df_freq::AbstractVector, df_time::AbstractVector, f_idx::Dict{Float64,Int}, t_idx::Dict{Float64,Int})
    n_rows = length(df_freq)
    fi_map = Vector{Int}(undef, n_rows)
    ti_map = Vector{Int}(undef, n_rows)
    @inbounds for r = 1:n_rows
        fi_map[r] = get(f_idx, Float64(df_freq[r]), 0)
        ti_map[r] = get(t_idx, Float64(df_time[r]), 0)
    end
    return fi_map, ti_map
end

function _fill_spec_col!(A::Matrix{Float64}, ci::Int, col::AbstractVector)
    @inbounds for r in axes(A, 1)
        A[r, ci] = Float64(col[r])
    end
    return nothing
end

function _fill_tf_channel_3d!(A::Array{Float64,3}, ci::Int, col::AbstractVector)
    n_f, n_t, _ = size(A)
    idx = 0
    @inbounds for ti = 1:n_t, fi = 1:n_f
        idx += 1
        A[fi, ti, ci] = Float64(col[idx])
    end
    return nothing
end

function _fill_tf_channel_3d_mapped!(A::Array{Float64,3}, ci::Int, col::AbstractVector, fi_map::Vector{Int}, ti_map::Vector{Int})
    @inbounds for r in eachindex(fi_map)
        fi = fi_map[r]
        ti = ti_map[r]
        if fi > 0 && ti > 0
            A[fi, ti, ci] = Float64(col[r])
        end
    end
    return nothing
end

function _fill_tf_channel_4d!(A::Array{Float64,4}, ep::Int, ci::Int, col::AbstractVector)
    _, n_f, n_t, _ = size(A)
    idx = 0
    @inbounds for ti = 1:n_t, fi = 1:n_f
        idx += 1
        A[ep, fi, ti, ci] = Float64(col[idx])
    end
    return nothing
end

function _fill_tf_channel_4d_mapped!(A::Array{Float64,4}, ep::Int, ci::Int, col::AbstractVector, fi_map::Vector{Int}, ti_map::Vector{Int})
    @inbounds for r in eachindex(fi_map)
        fi = fi_map[r]
        ti = ti_map[r]
        if fi > 0 && ti > 0
            A[ep, fi, ti, ci] = Float64(col[r])
        end
    end
    return nothing
end

function _populate_tf_3d!(A::Array{Float64,3}, df, freqs::Vector{Float64}, times::Vector{Float64}, channels::Vector{Symbol}, is_std::Bool)
    (isempty(freqs) || isempty(times) || isempty(channels)) && return fill!(A, NaN)
    if is_std
        for (ci, ch) in enumerate(channels)
            _fill_tf_channel_3d!(A, ci, df[!, ch])
        end
    else
        fill!(A, NaN)
        f_idx = Dict{Float64,Int}(f => i for (i, f) in enumerate(freqs))
        t_idx = Dict{Float64,Int}(t => i for (i, t) in enumerate(times))
        fi_map, ti_map = _map_tf_rows(df.freq, df.time, f_idx, t_idx)
        for (ci, ch) in enumerate(channels)
            _fill_tf_channel_3d_mapped!(A, ci, df[!, ch], fi_map, ti_map)
        end
    end
    return A
end

function DimensionalData.DimArray(dat::EegFun.SingleDataFrameEeg)
    return DimArray(Matrix(dat), (Ti(EegFun.time_vector(dat)), Dim{:Channel}(EegFun.channel_labels(dat))))
end

function DimensionalData.DimArray(dat::EegFun.MultiDataFrameEeg)
    return DimArray(
        Array(dat),
        (Dim{:Epoch}(1:EegFun.n_epochs(dat)), Ti(EegFun.time_vector(dat)), Dim{:Channel}(EegFun.channel_labels(dat))),
    )
end

function DimensionalData.DimArray(dat::EegFun.SpectrumData)
    c = EegFun.channel_labels(dat)
    f = Float64.(dat.data.freq)
    A = Matrix{Float64}(undef, length(f), length(c))
    for (ci, ch) in enumerate(c)
        _fill_spec_col!(A, ci, dat.data[!, ch])
    end
    return DimArray(A, (Dim{:Freq}(f), Dim{:Channel}(c)))
end

function DimensionalData.DimArray(dat::EegFun.TimeFreqData)
    df = dat.data_power
    freqs = EegFun.frequency_vector(dat)
    c = EegFun.channel_labels(dat)
    std_times = _try_extract_standard_times(df.freq, df.time, freqs)
    is_std = std_times !== nothing
    times = is_std ? std_times : sort!(Float64.(unique(df.time)))
    A = Array{Float64,3}(undef, length(freqs), length(times), length(c))
    _populate_tf_3d!(A, df, freqs, times, c, is_std)
    return DimArray(A, (Dim{:Freq}(freqs), Ti(times), Dim{:Channel}(c)))
end

function DimensionalData.DimArray(dat::EegFun.TimeFreqEpochData)
    n_ep = EegFun.n_epochs(dat)
    c = EegFun.channel_labels(dat)
    if n_ep == 0
        A0 = Array{Float64,4}(undef, 0, 0, 0, length(c))
        return DimArray(A0, (Dim{:Epoch}(1:0), Dim{:Freq}(Float64[]), Ti(Float64[]), Dim{:Channel}(c)))
    end

    freqs = EegFun.frequency_vector(dat)
    df1 = dat.data_power[1]
    std_times1 = _try_extract_standard_times(df1.freq, df1.time, freqs)
    times = std_times1 !== nothing ? std_times1 : sort!(Float64.(unique(df1.time)))
    n_f, n_t, n_c = length(freqs), length(times), length(c)
    A = Array{Float64,4}(undef, n_ep, n_f, n_t, n_c)

    if n_f == 0 || n_t == 0 || n_c == 0
        fill!(A, NaN)
        return DimArray(A, (Dim{:Epoch}(1:n_ep), Dim{:Freq}(freqs), Ti(times), Dim{:Channel}(c)))
    end

    f_idx = nothing
    t_idx = nothing
    for ep = 1:n_ep
        df = dat.data_power[ep]
        is_std = (ep == 1 && std_times1 !== nothing) || _is_standard_tf_grid(df.freq, df.time, freqs, times)
        if is_std
            for (ci, ch) in enumerate(c)
                _fill_tf_channel_4d!(A, ep, ci, df[!, ch])
            end
        else
            @views fill!(A[ep, :, :, :], NaN)
            if f_idx === nothing
                f_idx = Dict{Float64,Int}(f => i for (i, f) in enumerate(freqs))
                t_idx = Dict{Float64,Int}(t => i for (i, t) in enumerate(times))
            end
            fi_map, ti_map = _map_tf_rows(df.freq, df.time, f_idx, t_idx)
            for (ci, ch) in enumerate(c)
                _fill_tf_channel_4d_mapped!(A, ep, ci, df[!, ch], fi_map, ti_map)
            end
        end
    end

    return DimArray(A, (Dim{:Epoch}(1:n_ep), Dim{:Freq}(freqs), Ti(times), Dim{:Channel}(c)))
end

end
