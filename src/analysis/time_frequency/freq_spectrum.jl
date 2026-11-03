function _welch_channel_epochs!(
    avg_power::Vector{Float64},
    dfs::Vector{DataFrame},
    channel::Symbol,
    cfg::DSP.Periodograms.WelchConfig,
    sig_buf::Vector{Float64},
)
    fill!(avg_power, 0.0)
    n_ep = length(dfs)
    for ep_idx = 1:n_ep
        col = dfs[ep_idx][!, channel]
        pgram = if col isa Vector{Float64}
            DSP.welch_pgram(col, cfg)
        else
            @inbounds for i in eachindex(sig_buf)
                sig_buf[i] = Float64(col[i])
            end
            DSP.welch_pgram(sig_buf, cfg)
        end
        p = DSP.power(pgram)
        @inbounds @simd for i in eachindex(avg_power)
            avg_power[i] += p[i]
        end
    end
    avg_power .*= (1.0 / n_ep)
    return avg_power
end

"""
    freq_spectrum(dat::EegData;
                  channel_selection=channels(),
                  window_size::Union{Int,Nothing}=nothing,
                  overlap::Real=0.5,
                  window_function::Function=DSP.hanning,
                  max_freq::Union{Nothing,Real}=nothing)

Compute the Welch power spectral density (`SpectrumData`, in μV²/Hz) for `ContinuousData`,
`ErpData`, or `EpochData` (computed per epoch and averaged across epochs).

# Keyword Arguments
- `channel_selection=channels()`: Channel selection predicate.
- `window_size::Union{Int,Nothing}=nothing`: FFT window size in samples (defaults to `min(2048, n_samples(dat))`).
- `overlap::Real=0.5`: Window overlap fraction in `[0, 1)`.
- `window_function::Function=DSP.hanning`: Spectral window function.
- `max_freq::Union{Nothing,Real}=nothing`: Optional upper frequency cutoff in Hz.

# Examples
```julia
spectrum = freq_spectrum(epochs)
spectrum = freq_spectrum(epochs; channel_selection=channels(:Cz), window_size=2048, max_freq=100.0)
```
"""
function freq_spectrum(
    dat::EegData;
    channel_selection = channels(),
    window_size::Union{Int,Nothing} = nothing,
    overlap::Real = 0.5,
    window_function::Function = DSP.hanning,
    max_freq::Union{Nothing,Real} = nothing,
)
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected. Available channels: $(channel_labels(dat))")
    (dat isa MultiDataFrameEeg && isempty(dat.data)) && @minimal_error("Cannot compute frequency spectrum on empty epoch data")

    n_samp = n_samples(dat)
    n_samp < 1 && @minimal_error("Cannot compute frequency spectrum on data with $n_samp samples")

    actual_window_size = isnothing(window_size) ? min(2048, n_samp) : window_size
    (1 <= actual_window_size <= n_samp) ||
        @minimal_error("`window_size` must satisfy 1 <= window_size <= n_samples ($n_samp), got $actual_window_size")
    (0 <= overlap < 1) || @minimal_error("`overlap` must be in range [0, 1), got $overlap")

    noverlap = min(actual_window_size - 1, round(Int, actual_window_size * overlap))
    n_ch = length(selected_channels)
    n_workers = min(Threads.nthreads(), n_ch)

    chunks = collect(Iterators.partition(1:n_ch, cld(n_ch, n_workers)))
    n_chunks = length(chunks)

    # Pre-plan Welch configurations sequentially on main thread (FFTW planning is not thread-safe)
    cfgs = [
        DSP.Periodograms.WelchConfig(
            n_samp,
            Float64;
            n = actual_window_size,
            noverlap = noverlap,
            fs = dat.sample_rate,
            window = window_function,
        ) for _ = 1:n_chunks
    ]
    full_freqs = collect(Float64, cfgs[1].freq)
    n_full = length(full_freqs)
    keep_len = isnothing(max_freq) ? n_full : something(findlast(<=(max_freq), full_freqs), 0)
    freqs_out = keep_len == n_full ? full_freqs : full_freqs[1:keep_len]

    ch_cols = Vector{Vector{Float64}}(undef, n_ch)
    if dat isa MultiDataFrameEeg
        dfs = dat.data
        _foreach_chunk(n_chunks; ntasks = n_chunks) do tid_range
            for tid in tid_range
                cfg = cfgs[tid]
                avg_buf = Vector{Float64}(undef, n_full)
                sig_buf = Vector{Float64}(undef, n_samp)
                for ci in chunks[tid]
                    _welch_channel_epochs!(avg_buf, dfs, selected_channels[ci], cfg, sig_buf)
                    ch_cols[ci] = avg_buf[1:keep_len]
                end
            end
        end
    else
        _foreach_chunk(n_chunks; ntasks = n_chunks) do tid_range
            for tid in tid_range
                cfg = cfgs[tid]
                for ci in chunks[tid]
                    sig = channel_data(dat, selected_channels[ci])
                    p = DSP.power(DSP.welch_pgram(sig, cfg))
                    ch_cols[ci] = keep_len == n_full ? copy(p) : p[1:keep_len]
                end
            end
        end
    end

    cols = Vector{Any}(undef, n_ch + 1)
    cols[1] = freqs_out
    for ci = 1:n_ch
        cols[ci+1] = ch_cols[ci]
    end
    spectrum_df = DataFrame(cols, vcat([:freq], selected_channels); copycols = false)

    cond_num = hasproperty(dat, :condition) ? dat.condition : 0
    return SpectrumData(
        dat.file,
        cond_num,
        condition_name(dat),
        spectrum_df,
        copy(dat.layout),
        dat.sample_rate,
        :welch,
        copy(dat.analysis_info),
    )
end
