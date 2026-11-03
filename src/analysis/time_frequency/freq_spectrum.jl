"""
    freq_spectrum(dat::EegData;
                  channel_selection=channels(),
                  window_size::Union{Int,Nothing}=nothing,
                  overlap::Real=0.5,
                  window_function::Function=DSP.hanning,
                  max_freq::Union{Nothing,Real}=nothing)

Compute power spectrum using Welch's method for EEG data.

Compute the frequency-domain power spectrum (no time dimension) using
Welch's method with overlapping windows. For `EpochData`, Welch's periodogram is
computed separately for each epoch (never concatenating across trial discontinuities)
and averaged across epochs. For `ErpData` and `ContinuousData` (`SingleDataFrameEeg`),
power is computed directly from the single DataFrame.

# Arguments
- `dat::EegData`: EEG data (EpochData, ErpData, or ContinuousData)

# Keyword Arguments
- `channel_selection=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `window_size::Union{Int,Nothing}=nothing`: Size of the FFT window for spectral estimation (in samples).
  - If `nothing` (default), dynamically set to `min(2048, n_samples(dat))` for optimal high-resolution defaults that fit within each epoch/signal.
  - Must satisfy `1 <= window_size <= n_samples(dat)`.
- `overlap::Real=0.5`: Overlap fraction between windows (`0.0 <= overlap < 1.0`). Default is `0.5` (50% overlap)
- `window_function::Function=DSP.hanning`: Window function for spectral estimation
  - Options: `DSP.hanning`, `DSP.hamming`, `DSP.blackman`, etc.
- `max_freq::Union{Nothing,Real}=nothing`: Maximum frequency to return in Hz.
  - If `nothing`, returns all frequencies up to Nyquist
  - If specified, filters results to frequencies <= max_freq

# Returns
- `SpectrumData`: Power spectrum data structure with:
  - `data::DataFrame`: DataFrame with columns: freq, [electrode channels...] containing power spectral density (μV²/Hz)
  - Other metadata: file, condition, condition_name, layout, sample_rate, method, analysis_info

# Example
```julia
# Compute power spectrum for all channels
spectrum = freq_spectrum(epochs)

# Single channel with custom parameters
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
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected. Available channels: $(channel_labels(dat))")

    if dat isa MultiDataFrameEeg && isempty(dat.data)
        @minimal_error("Cannot compute frequency spectrum on empty epoch data")
    end

    n_samp = n_samples(dat)
    n_samp < 1 && @minimal_error("Cannot compute frequency spectrum on data with $n_samp samples")

    # Dynamically determine window size if not provided (clamped to n_samples per epoch/signal)
    actual_window_size = isnothing(window_size) ? min(2048, n_samp) : window_size
    if actual_window_size < 1 || actual_window_size > n_samp
        @minimal_error("`window_size` must satisfy 1 <= window_size <= n_samples ($n_samp), got $actual_window_size")
    end

    # Validate overlap
    (overlap < 0 || overlap >= 1) && @minimal_error("`overlap` must be in range [0, 1), got $overlap")

    noverlap = min(actual_window_size - 1, Int(round(actual_window_size * overlap)))

    first_channel = selected_channels[1]
    if dat isa MultiDataFrameEeg
        n_ep = length(dat.data)
        first_pgram = DSP.welch_pgram(
            Float64.(dat.data[1][!, first_channel]),
            actual_window_size,
            noverlap;
            fs = dat.sample_rate,
            window = window_function,
        )
        freqs = DSP.freq(first_pgram)
        spectrum_df = DataFrame(freq = freqs)
        inv_n_ep = 1.0 / n_ep

        for channel in selected_channels
            avg_power = zeros(Float64, length(freqs))
            for ep_idx = 1:n_ep
                ep_col = dat.data[ep_idx][!, channel]
                ep_sig = ep_col isa Vector{Float64} ? ep_col : Float64.(ep_col)
                pgram = DSP.welch_pgram(ep_sig, actual_window_size, noverlap; fs = dat.sample_rate, window = window_function)
                p = DSP.power(pgram)
                @inbounds @simd for i in eachindex(avg_power)
                    avg_power[i] += p[i]
                end
            end
            avg_power .*= inv_n_ep
            spectrum_df[!, channel] = avg_power
        end
    else
        first_signal = channel_data(dat, first_channel)
        pgram = DSP.welch_pgram(first_signal, actual_window_size, noverlap; fs = dat.sample_rate, window = window_function)
        freqs = DSP.freq(pgram)
        spectrum_df = DataFrame(freq = freqs)

        for channel in selected_channels
            signal = channel_data(dat, channel)
            pgram = DSP.welch_pgram(signal, actual_window_size, noverlap; fs = dat.sample_rate, window = window_function)
            spectrum_df[!, channel] = DSP.power(pgram)
        end
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = spectrum_df.freq .<= max_freq
        spectrum_df = spectrum_df[mask, :]
    end

    # Return SpectrumData type
    cond_num = hasproperty(dat, :condition) ? dat.condition : 0
    cond_name = condition_name(dat)

    return SpectrumData(dat.file, cond_num, cond_name, spectrum_df, copy(dat.layout), dat.sample_rate, :welch, copy(dat.analysis_info))
end
