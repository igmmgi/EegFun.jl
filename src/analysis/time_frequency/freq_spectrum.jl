"""
    freq_spectrum(dat::EegData;
                  channel_selection::Function=channels(),
                  window_size::Int=256,
                  overlap::Real=0.5,
                  window_function::Function=DSP.hanning,
                  max_freq::Union{Nothing,Real}=nothing)

Compute power spectrum using Welch's method for EEG data.

This function computes the frequency-domain power spectrum (no time dimension) using
Welch's method with overlapping windows. For EpochData, power is averaged across epochs.
For ErpData and ContinuousData (SingleDataFrameEeg), power is computed directly from the single DataFrame.

# Arguments
- `dat::EegData`: EEG data (EpochData, ErpData, or ContinuousData)

# Keyword Arguments
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `window_size::Union{Int,Nothing}=nothing`: Size of the FFT window for spectral estimation (in samples).
  - If `nothing` (default), dynamically set to `min(2048, n_samples)` for optimal high-resolution defaults that don't crash on short epochs.
- `overlap::Real=0.5`: Overlap fraction between windows (0.0 to 1.0). Default is 0.5 (50% overlap)
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
    channel_selection::Function = channels(),
    window_size::Union{Int,Nothing} = nothing,
    overlap::Real = 0.5,
    window_function::Function = DSP.hanning,
    max_freq::Union{Nothing,Real} = nothing,
)
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Get frequency vector from first channel's first signal
    first_channel = selected_channels[1]
    first_signal = channel_data(dat, first_channel)

    # Dynamically determine window size if not provided
    actual_window_size = isnothing(window_size) ? min(2048, length(first_signal)) : window_size

    # Validate overlap
    (overlap < 0 || overlap >= 1) && error("`overlap` must be in range [0, 1), got $overlap")

    noverlap = Int(round(actual_window_size * overlap))

    pgram = DSP.welch_pgram(first_signal, actual_window_size, noverlap; fs = dat.sample_rate, window = window_function)
    freqs = DSP.freq(pgram)

    # Initialize output DataFrame with frequency column
    spectrum_df = DataFrame(freq = freqs)

    # Process each channel
    for channel in selected_channels
        signal = channel_data(dat, channel)
        pgram = DSP.welch_pgram(signal, actual_window_size, noverlap; fs = dat.sample_rate, window = window_function)
        spectrum_df[!, channel] = DSP.power(pgram)
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
