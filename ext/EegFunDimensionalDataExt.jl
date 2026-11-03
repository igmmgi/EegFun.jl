module EegFunDimensionalDataExt

using EegFun
using DimensionalData

function _tf_power_cube(df, freqs::Vector{Float64}, times::Vector{Float64}, channels::Vector{Symbol})
    n_f = length(freqs)
    n_t = length(times)
    n_c = length(channels)
    A = fill(NaN, n_f, n_t, n_c)
    if n_f == 0 || n_t == 0 || n_c == 0
        return A
    end

    # Fast path when rows follow standard EegFun layout (freq inner, time outer)
    if size(df, 1) == n_f * n_t && df.freq == repeat(freqs, outer = n_t) && df.time == repeat(times, inner = n_f)
        for (ci, ch) in enumerate(channels)
            A[:, :, ci] = reshape(Float64.(df[!, ch]), n_f, n_t)
        end
    else
        f_idx = Dict{Float64,Int}(f => i for (i, f) in enumerate(freqs))
        t_idx = Dict{Float64,Int}(t => i for (i, t) in enumerate(times))
        df_freqs = Float64.(df.freq)
        df_times = Float64.(df.time)
        for r in eachindex(df_freqs)
            fi = get(f_idx, df_freqs[r], 0)
            ti = get(t_idx, df_times[r], 0)
            if fi > 0 && ti > 0
                for (ci, ch) in enumerate(channels)
                    A[fi, ti, ci] = Float64(df[r, ch])
                end
            end
        end
    end
    return A
end

function DimensionalData.DimArray(dat::EegFun.SingleDataFrameEeg)
    A = Matrix(dat)
    t = EegFun.time_vector(dat)
    c = EegFun.channel_labels(dat)

    # Dimensions: Time × Channel
    return DimArray(A, (Ti(t), Dim{:Channel}(c)))
end

function DimensionalData.DimArray(dat::EegFun.MultiDataFrameEeg)
    A = Array(dat)
    ep = 1:EegFun.n_epochs(dat)
    t = EegFun.time_vector(dat)
    c = EegFun.channel_labels(dat)

    # Dimensions: Epoch × Time × Channel
    return DimArray(A, (Dim{:Epoch}(ep), Ti(t), Dim{:Channel}(c)))
end

function DimensionalData.DimArray(dat::EegFun.SpectrumData)
    c = EegFun.channel_labels(dat)
    f = Float64.(dat.data.freq)
    A = isempty(c) ? Array{Float64,2}(undef, length(f), 0) : Matrix{Float64}(dat.data[!, c])

    # Dimensions: Freq × Channel
    return DimArray(A, (Dim{:Freq}(f), Dim{:Channel}(c)))
end

function DimensionalData.DimArray(dat::EegFun.TimeFreqData)
    freqs = EegFun.frequency_vector(dat)
    times = Float64.(sort(unique(dat.data_power.time)))
    c = EegFun.channel_labels(dat)
    A = _tf_power_cube(dat.data_power, freqs, times, c)

    # Dimensions: Freq × Time × Channel
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
    times = Float64.(sort(unique(dat.data_power[1].time)))
    A = Array{Float64,4}(undef, n_ep, length(freqs), length(times), length(c))
    for ep = 1:n_ep
        A[ep, :, :, :] = _tf_power_cube(dat.data_power[ep], freqs, times, c)
    end

    # Dimensions: Epoch × Freq × Time × Channel
    return DimArray(A, (Dim{:Epoch}(1:n_ep), Dim{:Freq}(freqs), Ti(times), Dim{:Channel}(c)))
end

end
