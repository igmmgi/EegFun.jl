using FFTW
using LinearAlgebra
using Statistics
using DataFrames
using DSP

export morlet_wavelet, tf_wavelet, tf_analysis

# ============================================================================
# Workspace Structs for Performance
# ============================================================================

"""
    TFMorletWorkspace

Workspace for batch Morlet wavelet processing.
Pre-allocates FFT buffers and plans for frequency-domain convolution.
"""
struct TFMorletWorkspace
    padded_signal::Matrix{ComplexF64} # (n_conv_batch x n_channels)
    signal_fft::Matrix{ComplexF64}
    conv_fft::Matrix{ComplexF64}
    eegconv::Matrix{ComplexF64}
    wavelet_padded::Vector{ComplexF64}
    wavelet_fft::Vector{ComplexF64}
    p_fft::Any
    p_ifft::Any
    p_wavelet_fft::Any
end

function TFMorletWorkspace(n_conv_batch, n_channels, n_trials, n_times)
    ps = zeros(ComplexF64, n_conv_batch, n_channels)
    sf = zeros(ComplexF64, n_conv_batch, n_channels)
    cf = zeros(ComplexF64, n_conv_batch, n_channels)
    ec = zeros(ComplexF64, n_conv_batch, n_channels)
    wp = zeros(ComplexF64, n_conv_batch)
    wf = zeros(ComplexF64, n_conv_batch)

    p_fft = plan_fft(ps, 1, flags = FFTW.ESTIMATE)
    p_ifft = plan_ifft(cf, 1, flags = FFTW.ESTIMATE)
    p_wavelet_fft = plan_fft(wp, flags = FFTW.ESTIMATE)

    return TFMorletWorkspace(ps, sf, cf, ec, wp, wf, p_fft, p_ifft, p_wavelet_fft)
end

"""
    TFHanningFixedWorkspace

Workspace for Hanning taper with fixed window length (Short-Time FFT).
Pre-allocates buffers for real FFT (rfft) which is faster for real input.
"""
struct TFHanningFixedWorkspace
    tmpdat::Matrix{Float64}  # Windowed signal buffer (real)
    fdat::Matrix{ComplexF64}  # FFT output buffer
    p_rfft::Any              # Real FFT plan
end

function TFHanningFixedWorkspace(timewinidx, n_trials)
    tmpdat = zeros(Float64, timewinidx, n_trials)
    n_rfft_out = div(timewinidx, 2) + 1
    fdat = zeros(ComplexF64, n_rfft_out, n_trials)
    p_rfft = plan_rfft(tmpdat, 1, flags = FFTW.ESTIMATE)
    return TFHanningFixedWorkspace(tmpdat, fdat, p_rfft)
end

"""
    TFHanningAdaptiveWorkspace

Workspace for Hanning taper with frequency-dependent window length.
Pre-allocates buffers using the maximum window size needed across all frequencies.
"""
struct TFHanningAdaptiveWorkspace
    tmpdat::Matrix{Float64}   # Windowed signal buffer (max size, real)
    fdat::Matrix{ComplexF64}  # FFT output buffer (max size)
    p_rfft::Any               # Real FFT plan
end

function TFHanningAdaptiveWorkspace(max_timewinidx::Int, n_trials::Int)
    tmpdat = zeros(Float64, max_timewinidx, n_trials)
    n_rfft_out = div(max_timewinidx, 2) + 1
    fdat = zeros(ComplexF64, n_rfft_out, n_trials)
    p_rfft = plan_rfft(tmpdat, 1, flags = FFTW.ESTIMATE)
    return TFHanningAdaptiveWorkspace(tmpdat, fdat, p_rfft)
end

"""
    TFMultitaperWorkspace

Workspace for multi-taper method (DPSS tapers) - Cohen Chapter 16.
Pre-allocates buffers for windowed data, tapered data, and FFT results.
"""
struct TFMultitaperWorkspace
    tmpdat::Matrix{Float64}      # Windowed signal buffer (real)
    tapered_buf::Matrix{Float64} # Tapered signal buffer (real)
    fdat::Matrix{ComplexF64}     # FFT output buffer
    taper_power::Matrix{Float64} # Accumulated power across tapers
    p_rfft::Any                  # Real FFT plan
end

function TFMultitaperWorkspace(timewinidx, n_trials, n_freqs)
    tmpdat = zeros(Float64, timewinidx, n_trials)
    tapered_buf = zeros(Float64, timewinidx, n_trials)
    n_rfft_out = div(timewinidx, 2) + 1
    fdat = zeros(ComplexF64, n_rfft_out, n_trials)
    taper_power = zeros(n_freqs, n_trials)
    p_rfft = plan_rfft(tapered_buf, 1, flags = FFTW.ESTIMATE)
    return TFMultitaperWorkspace(tmpdat, tapered_buf, fdat, taper_power, p_rfft)
end

"""
    morlet_wavelet(f, fs; n_cycles=7)

Create a complex Morlet wavelet for frequency `f` at sampling rate `fs`.
`n_cycles` defines the width of the Gaussian (number of cycles).
"""
function morlet_wavelet(f, fs; n_cycles = 7)
    s = n_cycles / (2 * pi * f)
    wavelet_time = -2:(1/fs):2

    # Normalization factor
    # A = 1 / sqrt(s * sqrt(pi))
    # This normalization ensures that the peak of the wavelet's FFT spectrum is 1.0 (approx)
    # allowing for direct comparison of power across frequencies.
    A = 1 / sqrt(s * sqrt(pi))

    sine_wave = exp.(2 * im * pi * f .* wavelet_time)
    gaussian_win = exp.(-wavelet_time .^ 2 ./ (2 * s^2))

    return A .* sine_wave .* gaussian_win
end

"""
    tf_wavelet(epochs::EpochData, channel, freqs; n_cycles=7)

Perform time-frequency analysis using Morlet wavelets on epoched EEG data.
Extracts data for the specified channel from all epochs and processes them together.

# Arguments
- `epochs::EpochData`: Epoched EEG data
- `channel::Symbol`: Channel name to analyze (e.g., `:Channel1`)
- `freqs`: Vector of frequencies to analyze (Hz)

# Keyword Arguments
- `n_cycles::Number=7`: Number of wavelet cycles (can be a single number or a tuple for frequency-dependent cycles)

# Returns
- Power matrix (freqs × samples) averaged across epochs
"""
function tf_wavelet(epochs::EpochData, channel::Symbol, freqs; n_cycles = 7)
    # Extract channel data from all epochs
    n_epochs = length(epochs.data)
    n_samples = nrow(epochs.data[1])
    
    # Stack epochs into matrix (samples × trials)
    data_matrix = zeros(Float64, n_samples, n_epochs)
    for (ei, epoch_df) in enumerate(epochs.data)
        if !hasproperty(epoch_df, channel)
            error("Channel $channel not found in epoch data. Available channels: $(filter(c -> c != :time && c != :epoch, propertynames(epoch_df)))")
        end
        data_matrix[:, ei] .= epoch_df[!, channel]
    end
    
    # Use existing tf_wavelet implementation
    return tf_wavelet(data_matrix, Float64(epochs.sample_rate), freqs; n_cycles = n_cycles)
end

"""
    tf_wavelet(data, fs, freqs; n_cycles=7)

Perform time-frequency analysis using Morlet wavelets and FFT convolution.
`data` can be a vector (samples,) or a matrix (samples x trials).
`freqs` is a vector of frequencies to analyze.
Returns power matrix (freqs x samples).
"""
function tf_wavelet(data, fs, freqs; n_cycles = 7)
    n_samples = size(data, 1)
    if data isa Vector
        n_trials = 1
        data = reshape(data, :, 1)
    else
        n_trials = size(data, 2)
    end

    n_freqs = length(freqs)
    tf = zeros(n_freqs, n_samples)

    # Setup wavelet parameters
    wavelet_time = -2:(1/fs):2
    n_wavelet = length(wavelet_time)
    half_wavelet = floor(Int, n_wavelet / 2)

    n_conv = n_samples + n_wavelet - 1
    n_fft = nextpow(2, n_conv)

    # FFT data
    data_padded = [data; zeros(n_fft - n_samples, n_trials)]
    data_fft = fft(data_padded, 1)

    for (fi, f) in enumerate(freqs)
        # Create wavelet
        s = n_cycles / (2 * pi * f)
        A = 1 / sqrt(s * sqrt(pi))

        # We can optimize by only creating the wavelet of necessary length
        # but +/- 2s is safe and usually sufficient.
        wavelet = A .* exp.(2 * im * pi * f .* wavelet_time) .* exp.(-wavelet_time .^ 2 ./ (2 * s^2))

        wavelet_padded = [wavelet; zeros(n_fft - n_wavelet)]
        wavelet_fft = fft(wavelet_padded)

        # Batch convolution for all trials at once
        # Broadcoast multiplication along trials dimension
        conv_fft = data_fft .* wavelet_fft
        conv_res_full = ifft(conv_fft, 1)

        # Remove padding and extract relevant part (samples x trials)
        conv_res = conv_res_full[half_wavelet+1:half_wavelet+n_samples, :]

        # Compute power: |conv_res|^2
        # Average power across trials
        # sum |z|^2 / n_trials
        # Note: abs2.(conv_res) works on the whole matrix
        power_all_trials = abs2.(conv_res)
        tf[fi, :] .= vec(mean(power_all_trials, dims = 2))
    end

    return tf
end

# ============================================================================
# Short-Time FFT (Hanning Taper) - Cohen Chapter 15
# ============================================================================

"""
    tf_hanning(signal, times, sample_rate, frequencies, time_steps;
               window_length=nothing, cycles=nothing, keeptrials=false, workspace=nothing)

Perform Short-Time FFT analysis using Hanning taper (Cohen Chapter 15).

This implements the Short-Time FFT method where a fixed or frequency-dependent
window is applied at each time point, followed by FFT to extract power.

# Arguments
- `signal`: Signal data (samples × trials) or (samples,)
- `times`: Time vector for the signal
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Time points of interest (seconds)

# Keyword Arguments
- `window_length`: Fixed window length in seconds (for fixed window method)
- `cycles`: Number of cycles per frequency (for adaptive window method)
- `keeptrials`: Keep individual trials (true) or average (false)
- `workspace`: Pre-allocated workspace (optional, for performance)

# Returns
- `power`: Power matrix (freqs × times) or (freqs × times × trials) if keeptrials=true
- `times_out`: Output time vector
- `freqs_out`: Output frequency vector

# Example
```julia
# Fixed window length (100ms)
power, times, freqs = tf_hanning(signal, times, 256, 2:1:80, -0.5:0.01:1.0; 
                                  window_length=0.1)

# Frequency-dependent window (5 cycles per frequency)
power, times, freqs = tf_hanning(signal, times, 256, 2:1:80, -0.5:0.01:1.0; 
                                  cycles=5)
```
"""
function tf_hanning(
    signal,
    times,
    sample_rate,
    frequencies,
    time_steps;
    window_length = nothing,
    cycles = nothing,
    keeptrials::Bool = false,
    workspace = nothing,
)
    if isnothing(window_length) && isnothing(cycles)
        throw(ArgumentError("Must specify either window_length or cycles"))
    end

    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)

    time_of_interest = collect(Float64, time_steps)
    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_of_interest]
    n_frex = length(frequencies)
    n_timepoints = length(time_of_interest)

    # Fixed window length (Cohen's Short-Time FFT - Figure 15.2)
    if isnothing(cycles)
        timewinidx = round(Int, window_length * sample_rate)
        timewinidx += iseven(timewinidx)  # Make odd for symmetric window
        
        half_win = timewinidx ÷ 2

        # Create Hanning window (Cohen eq. 15.1)
        hann_win = 0.5 * (1 .- cos.(2π * (0:(timewinidx-1)) / (timewinidx - 1)))

        # Frequency indices for requested frequencies
        # FFT frequencies: linspace(0, sample_rate/2, floor(timewinidx/2)+1)
        freq_indices = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
        # Clamp to valid range
        max_freq_idx = div(timewinidx, 2) + 1
        freq_indices = clamp.(freq_indices, 1, max_freq_idx)

        # Create workspace
        ws = isnothing(workspace) ? TFHanningFixedWorkspace(timewinidx, n_trials) : workspace

        # Allocate output (power and complex)
        if keeptrials
            tf_trials = zeros(n_frex, n_timepoints, n_trials)
            tf_complex = zeros(ComplexF64, n_frex, n_timepoints, n_trials)
        else
            tf_avg = zeros(n_frex, n_timepoints)
            tf_complex_avg = zeros(ComplexF64, n_frex, n_timepoints)
        end
        inv_n_trials = 1.0 / n_trials

        @inbounds for (ti, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win

            # Extract signal and apply window (Cohen's approach)
            if start_idx >= 1 && end_idx <= n_samples
                for trial = 1:n_trials
                    @simd for s = 1:timewinidx
                        ws.tmpdat[s, trial] = signal[start_idx+s-1, trial] * hann_win[s]
                    end
                end
            else
                # Handle edge cases (zero-padding)
                fill!(ws.tmpdat, 0.0)
                v_start = max(1, start_idx)
                v_end = min(n_samples, end_idx)
                offset = v_start - start_idx
                len = v_end - v_start + 1
                for trial = 1:n_trials
                    @simd for s = 1:len
                        ws.tmpdat[offset+s, trial] = signal[v_start+s-1, trial] * hann_win[offset+s]
                    end
                end
            end

            # Real-to-Complex FFT (Cohen: fft(taperdat,[],1)/timewinidx)
            mul!(ws.fdat, ws.p_rfft, ws.tmpdat)
            ws.fdat ./= timewinidx  # Normalize by window length

            # Extract power and complex at requested frequencies
            if keeptrials
                for trial = 1:n_trials
                    for (fi, f_idx) in enumerate(freq_indices)
                        z = ws.fdat[f_idx, trial]
                        tf_trials[fi, ti, trial] = abs2(z)
                        tf_complex[fi, ti, trial] = z
                    end
                end
            else
                # Average over trials
                for trial = 1:n_trials
                    for (fi, f_idx) in enumerate(freq_indices)
                        z = ws.fdat[f_idx, trial]
                        tf_avg[fi, ti] += abs2(z) * inv_n_trials
                        tf_complex_avg[fi, ti] += z * inv_n_trials
                    end
                end
            end
        end

        if !keeptrials
            return tf_avg, tf_complex_avg, time_of_interest, frequencies
        else
            return tf_trials, tf_complex, time_of_interest, frequencies
        end
    end

    # Frequency-dependent window length (cycles specified)
    # Calculate maximum window size for workspace allocation
    max_timewinidx = round(Int, cycles / minimum(frequencies) * sample_rate)
    max_timewinidx += iseven(max_timewinidx)
    
    # Create workspace with maximum window size
    ws = isnothing(workspace) ? TFHanningAdaptiveWorkspace(max_timewinidx, n_trials) : workspace

    # Allocate output (power and complex)
    if keeptrials
        tf_trials = zeros(n_frex, n_timepoints, n_trials)
        tf_complex = zeros(ComplexF64, n_frex, n_timepoints, n_trials)
    else
        tf_avg = zeros(n_frex, n_timepoints)
        tf_complex_avg = zeros(ComplexF64, n_frex, n_timepoints)
    end
    inv_n_trials = 1.0 / n_trials

    # Loop over frequencies (each has different window length)
    for (fi, freq) in enumerate(frequencies)
        timewinidx = round(Int, cycles / freq * sample_rate)
        timewinidx += iseven(timewinidx)
        half_win = timewinidx ÷ 2

        # Create Hanning window for this frequency
        hann_win = 0.5 * (1 .- cos.(2π * (0:(timewinidx-1)) / (timewinidx - 1)))

        # Frequency index for this frequency
        f_idx = round(Int, freq * timewinidx / sample_rate) + 1
        f_idx = clamp(f_idx, 1, div(timewinidx, 2) + 1)

        @inbounds for (ti, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win

            # Extract signal and apply window
            if start_idx >= 1 && end_idx <= n_samples
                for trial = 1:n_trials
                    @simd for s = 1:timewinidx
                        ws.tmpdat[s, trial] = signal[start_idx+s-1, trial] * hann_win[s]
                    end
                end
            else
                fill!(view(ws.tmpdat, 1:timewinidx, :), 0.0)
                v_start = max(1, start_idx)
                v_end = min(n_samples, end_idx)
                offset = v_start - start_idx
                len = v_end - v_start + 1
                for trial = 1:n_trials
                    @simd for s = 1:len
                        ws.tmpdat[offset+s, trial] = signal[v_start+s-1, trial] * hann_win[offset+s]
                    end
                end
            end

            # Real-to-Complex FFT
            # Create a view and compute FFT (rfft will handle the plan internally)
            tmpdat_view = view(ws.tmpdat, 1:timewinidx, :)
            fdat = rfft(tmpdat_view, 1) ./ timewinidx

            # Extract power and complex at this frequency
            if keeptrials
                for trial = 1:n_trials
                    z = fdat[f_idx, trial]
                    tf_trials[fi, ti, trial] = abs2(z)
                    tf_complex[fi, ti, trial] = z
                end
            else
                for trial = 1:n_trials
                    z = fdat[f_idx, trial]
                    tf_avg[fi, ti] += abs2(z) * inv_n_trials
                    tf_complex_avg[fi, ti] += z * inv_n_trials
                end
            end
        end
    end

    if !keeptrials
        return tf_avg, tf_complex_avg, time_of_interest, frequencies
    else
        return tf_trials, tf_complex, time_of_interest, frequencies
    end
end

# ============================================================================
# Multi-Taper Method (DPSS) - Cohen Chapter 16
# ============================================================================

"""
    tf_multitaper(signal, times, sample_rate, frequencies, time_steps;
                  time_window_length=0.5, time_bandwidth_product=3.0, keeptrials=false, workspace=nothing)

Perform time-frequency analysis using multi-taper method with DPSS (Slepian) tapers (Cohen Chapter 16).

This implements the multi-taper method where multiple orthogonal tapers (DPSS) are applied
to the same data window, and power is averaged across tapers to reduce variance.

# Arguments
- `signal`: Signal data (samples × trials) or (samples,)
- `times`: Time vector for the signal
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Time points of interest (seconds)

# Keyword Arguments
- `time_window_length`: Time window length in seconds (default: 0.5)
- `time_bandwidth_product`: Time-bandwidth product (nw) - controls number of tapers (default: 3.0)
  - Higher values = more tapers = lower variance but more frequency smoothing
  - Cohen uses nw_product = 3, which gives 2*nw-1 = 5 tapers
- `keeptrials`: Keep individual trials (true) or average (false)
- `workspace`: Pre-allocated workspace (optional, for performance)

# Returns
- `power`: Power matrix (freqs × times) or (freqs × times × trials) if keeptrials=true
- `times_out`: Output time vector
- `freqs_out`: Output frequency vector

# Example
```julia
# Multi-taper with 400ms window and nw=3 (Cohen's Figure 16.2)
power, times, freqs = tf_multitaper(signal, times, 256, 2:1:80, -0.3:0.05:1.0; 
                                    time_window_length=0.4, time_bandwidth_product=3.0)
```
"""
function tf_multitaper(
    signal,
    times,
    sample_rate,
    frequencies,
    time_steps;
    time_window_length::Float64 = 0.5,
    time_bandwidth_product::Float64 = 3.0,
    keeptrials::Bool = false,
    workspace = nothing,
)
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)

    time_of_interest = collect(Float64, time_steps)
    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_of_interest]
    n_frex = length(frequencies)
    n_timepoints = length(time_of_interest)

    # Window parameters (Cohen: timewinidx = round(timewin/(1000/EEG.srate)))
    timewinidx = round(Int, time_window_length * sample_rate)
    timewinidx += iseven(timewinidx)  # Make odd for symmetric window
    half_win = timewinidx ÷ 2

    # Generate DPSS tapers (Cohen: tapers = dpss(timewinidx, nw_product))
    # nw = time_bandwidth_product (Cohen's nw_product)
    nw = time_bandwidth_product
    n_tapers = max(1, floor(Int, 2 * nw - 1))  # Number of tapers = 2*nw - 1
    tapers = DSP.dpss(timewinidx, nw, n_tapers)

    # Frequency indices (Cohen: f = linspace(0, sample_rate/2, floor(timewinidx/2)+1))
    # We extract power at requested frequencies
    freq_indices = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
    max_freq_idx = div(timewinidx, 2) + 1
    freq_indices = clamp.(freq_indices, 1, max_freq_idx)

    # Create workspace
    ws = isnothing(workspace) ? TFMultitaperWorkspace(timewinidx, n_trials, n_frex) : workspace

    # Allocate output (power and complex)
    if keeptrials
        tf_trials = zeros(n_frex, n_timepoints, n_trials)
        tf_complex = zeros(ComplexF64, n_frex, n_timepoints, n_trials)
    else
        tf_avg = zeros(n_frex, n_timepoints)
        tf_complex_avg = zeros(ComplexF64, n_frex, n_timepoints)
    end
    inv_n_trials = 1.0 / n_trials
    inv_n_tapers = 1.0 / n_tapers

    @inbounds for (ti, center_idx) in enumerate(tois_idx)
        start_idx = center_idx - half_win
        end_idx = center_idx + half_win

        # Extract signal window (Cohen: data = EEG.data(..., times2saveidx(ti)-floor(timewinidx/2)+1:times2saveidx(ti)+ceil(timewinidx/2), :))
        if start_idx >= 1 && end_idx <= n_samples
            for trial = 1:n_trials
                @simd for s = 1:timewinidx
                    ws.tmpdat[s, trial] = signal[start_idx+s-1, trial]
                end
            end
        else
            # Handle edge cases (zero-padding)
            fill!(ws.tmpdat, 0.0)
            v_start = max(1, start_idx)
            v_end = min(n_samples, end_idx)
            offset = v_start - start_idx
            len = v_end - v_start + 1
            for trial = 1:n_trials
                @simd for s = 1:len
                    ws.tmpdat[offset+s, trial] = signal[v_start+s-1, trial]
                end
            end
        end

        # Reset taper power accumulator (Cohen: taperpow = zeros(...))
        fill!(ws.taper_power, 0.0)
        
        # Accumulate complex values across tapers
        taper_complex = zeros(ComplexF64, n_frex, n_trials)

        # Loop through tapers (Cohen: for tapi = 1:size(tapers,2)-1)
        # Note: Cohen uses size(tapers,2)-1, but we'll use all tapers
        for tapi = 1:n_tapers
            # Apply taper to data (Cohen: data = bsxfun(@times, ..., tapers(:,tapi)))
            for trial = 1:n_trials
                @simd for s = 1:timewinidx
                    ws.tapered_buf[s, trial] = ws.tmpdat[s, trial] * tapers[s, tapi]
                end
            end

            # Real-to-Complex FFT (Cohen: pow = fft(data, timewinidx)/timewinidx)
            mul!(ws.fdat, ws.p_rfft, ws.tapered_buf)
            ws.fdat ./= timewinidx  # Normalize by window length

            # Extract power and complex, accumulate across tapers
            # For power: average abs2 (Cohen: taperpow = taperpow + mean(pow.*conj(pow), 2))
            # For complex: average the complex values
            for trial = 1:n_trials
                for (fi, f_idx) in enumerate(freq_indices)
                    z = ws.fdat[f_idx, trial]
                    ws.taper_power[fi, trial] += abs2(z)
                    taper_complex[fi, trial] += z
                end
            end
        end

        # Average across tapers (Cohen: taperpow / tapi)
        ws.taper_power .*= inv_n_tapers
        taper_complex .*= inv_n_tapers

        # Store results
        if keeptrials
            for trial = 1:n_trials
                for fi = 1:n_frex
                    tf_trials[fi, ti, trial] = ws.taper_power[fi, trial]
                    tf_complex[fi, ti, trial] = taper_complex[fi, trial]
                end
            end
        else
            # Average across trials
            for trial = 1:n_trials
                for fi = 1:n_frex
                    tf_avg[fi, ti] += ws.taper_power[fi, trial] * inv_n_trials
                    tf_complex_avg[fi, ti] += taper_complex[fi, trial] * inv_n_trials
                end
            end
        end
    end

    if !keeptrials
        return tf_avg, tf_complex_avg, time_of_interest, frequencies
    else
        return tf_trials, tf_complex, time_of_interest, frequencies
    end
end

# ============================================================================
# Validation Functions
# ============================================================================

"""
    _validate_tf_method_options(method::Symbol, kwargs)

Internal function: Validate that kwargs are appropriate for the selected method.

Raises `ArgumentError` if invalid options are provided.
"""
function _validate_tf_method_options(method::Symbol, kwargs)
    # Define valid options for each method
    valid_options = Dict(
        :wavelet => Set([:width, :log_freqs]),
        :superlet => Set([:order, :width, :combine, :log_freqs]),
        :multitaper => Set([:time_window_length, :time_bandwidth_product]),
        :hanning_fixed => Set([:time_window_length]),
        :hanning_adaptive => Set([:cycles]),
    )

    # Check if method is valid
    if !haskey(valid_options, method)
        valid_methods = join(sort(collect(keys(valid_options))), ", ")
        error("Unknown method `:$method`. Valid methods are: $valid_methods")
    end

    # Common options that are always valid
    common_options = Set([:keep_epochs, :channel_selection, :method])

    # Get valid options for this method
    method_valid = valid_options[method]
    all_valid = union(method_valid, common_options)

    # Check for invalid options
    invalid_options = Set{Symbol}()
    for (key, _) in kwargs
        if !(key in all_valid)
            push!(invalid_options, key)
        end
    end

    if !isempty(invalid_options)
        invalid_str = join(sort(collect(invalid_options)), ", ")
        valid_str = join(sort(collect(method_valid)), ", ")
        error(
            "Invalid option(s) for method `:$method`: $invalid_str. " * "Valid options for `:$method` are: $valid_str",
        )
    end

    # Method-specific validation warnings
    if !(method in [:wavelet, :superlet]) && haskey(kwargs, :log_freqs)
        @minimal_warning "Option `log_freqs` is only used by `:wavelet` and `:superlet` methods. " *
                         "It will be ignored for method `:$method`."
    end
end

# ============================================================================
# Batch Processing Functions
# ============================================================================

"""
    _tf_batch_wavelet!(power_results, ws, n_samples, n_trials, n_channels, tois_idx, freqs_out, cycles, sr, keeptrials)

Internal function: Batch process all frequencies using pre-computed FFT of signal.
This is much more efficient than computing FFT per frequency.
"""
function _tf_batch_wavelet!(
    power_results,      # Matrix/Array to fill with power
    complex_results,    # Matrix/Array to fill with complex values (optional)
    ws::TFMorletWorkspace,
    n_samples,
    n_trials,
    n_channels,
    tois_idx,
    freqs_out,
    cycles,
    sr,
    keeptrials,
)
    n_conv_batch = size(ws.padded_signal, 1)
    inv_n_trials = 1.0 / n_trials

    # FFT of signal is already done in ws.signal_fft

    for (fi, freq) in enumerate(freqs_out)
        sigma = cycles[fi] / (2π * freq)
        hw = ceil(Int, 6 * sigma * sr) ÷ 2
        wl = hw * 2 + 1
        w_time = range(-wl / 2, wl / 2, length = wl) ./ sr

        fill!(ws.wavelet_padded, 0.0)
        @. ws.wavelet_padded[1:wl] =
            sqrt(1 / (sigma * sqrt(pi))) * exp(2im * pi * freq * w_time) * exp(-w_time^2 / (2 * sigma^2))
        mul!(ws.wavelet_fft, ws.p_wavelet_fft, ws.wavelet_padded)

        for ci = 1:n_channels
            @simd for i = 1:n_conv_batch
                ws.conv_fft[i, ci] = ws.signal_fft[i, ci] * ws.wavelet_fft[i]
            end
        end
        mul!(ws.eegconv, ws.p_ifft, ws.conv_fft)

        valid_start = hw + 1
        n_tois = length(tois_idx)
        
        for ci = 1:n_channels
            for trial = 1:n_trials
                t_offset = (trial - 1) * n_samples
                for ti = 1:n_tois
                    t_idx = tois_idx[ti]
                    # t_idx is 1-indexed position in original time array
                    # After concatenation: trial starts at t_offset+1
                    # After convolution: valid data starts at valid_start
                    # So sample t_idx from trial is at: valid_start + t_offset + t_idx - 1
                    @inbounds z = ws.eegconv[valid_start + t_offset + t_idx - 1, ci]
                    @inbounds val = abs2(z)
                    
                    if keeptrials
                        # Handle List of Matrices (channel_results)
                        if power_results isa Vector
                            @inbounds power_results[ci][fi, ti, trial] = val
                        else
                            @inbounds power_results[fi, ti, trial, ci] = val
                        end
                        if !isnothing(complex_results)
                            if complex_results isa Vector
                                @inbounds complex_results[ci][fi, ti, trial] = z
                            else
                                @inbounds complex_results[fi, ti, trial, ci] = z
                            end
                        end
                    else
                        # Average across trials
                        if power_results isa Vector
                            @inbounds power_results[ci][fi, ti] += val * inv_n_trials
                        elseif ndims(power_results) == 3
                            @inbounds power_results[fi, ti, ci] += val * inv_n_trials
                        else
                            @inbounds power_results[fi, ti] += val * inv_n_trials
                        end
                        if !isnothing(complex_results)
                            if complex_results isa Vector
                                @inbounds complex_results[ci][fi, ti] += z * inv_n_trials
                            elseif ndims(complex_results) == 3
                                @inbounds complex_results[fi, ti, ci] += z * inv_n_trials
                            else
                                @inbounds complex_results[fi, ti] += z * inv_n_trials
                            end
                        end
                    end
                end
            end
        end
    end
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
    _power_to_tf_dataframe(power, times, freqs, channel_labels) -> DataFrame
    _power_to_tf_dataframe(complex_data, times, freqs, channel_labels) -> (DataFrame, DataFrame)

Convert 3D power array (freqs × times × channels) or complex array to DataFrame format.

If complex data is provided, returns a tuple of (power_df, phase_df).
If power (real) is provided, returns a single power DataFrame.
"""
function _power_to_tf_dataframe(
    power::AbstractArray{<:Real,3},
    times::AbstractVector,
    freqs::AbstractVector,
    channel_labels::Vector{Symbol},
)
    n_freqs, n_times, n_channels = size(power)
    n_rows = n_freqs * n_times

    # Pre-allocate all columns at once
    cols = Dict{Symbol,Vector{Float64}}()
    cols[:time] = Vector{Float64}(undef, n_rows)
    cols[:freq] = Vector{Float64}(undef, n_rows)
    for ch in channel_labels
        cols[ch] = Vector{Float64}(undef, n_rows)
    end

    # Fill columns in a single pass
    idx = 1
    for ti = 1:n_times
        for fi = 1:n_freqs
            cols[:time][idx] = times[ti]
            cols[:freq][idx] = freqs[fi]
            for (ci, ch) in enumerate(channel_labels)
                cols[ch][idx] = power[fi, ti, ci]
            end
            idx += 1
        end
    end

    # Create DataFrame with explicit column order: time, freq, then channels in order
    column_order = [:time, :freq, channel_labels...]
    return DataFrame(cols, copycols = false)[!, column_order]
end

function _power_to_tf_dataframe(
    complex_data::AbstractArray{<:Complex,3},
    times::AbstractVector,
    freqs::AbstractVector,
    channel_labels::Vector{Symbol},
)
    n_freqs, n_times, n_channels = size(complex_data)
    n_rows = n_freqs * n_times

    # Pre-allocate columns for power and phase DataFrames
    cols_power = Dict{Symbol,Vector{Float64}}()
    cols_phase = Dict{Symbol,Vector{Float64}}()
    cols_power[:time] = Vector{Float64}(undef, n_rows)
    cols_power[:freq] = Vector{Float64}(undef, n_rows)
    cols_phase[:time] = Vector{Float64}(undef, n_rows)
    cols_phase[:freq] = Vector{Float64}(undef, n_rows)
    for ch in channel_labels
        cols_power[ch] = Vector{Float64}(undef, n_rows)
        cols_phase[ch] = Vector{Float64}(undef, n_rows)
    end

    # Fill columns in a single pass
    idx = 1
    for ti = 1:n_times
        for fi = 1:n_freqs
            cols_power[:time][idx] = times[ti]
            cols_power[:freq][idx] = freqs[fi]
            cols_phase[:time][idx] = times[ti]
            cols_phase[:freq][idx] = freqs[fi]
            for (ci, ch) in enumerate(channel_labels)
                z = complex_data[fi, ti, ci]
                cols_power[ch][idx] = abs2(z)  # Power: |z|^2
                cols_phase[ch][idx] = angle(z)  # Phase: arg(z) in radians
            end
            idx += 1
        end
    end

    # Create DataFrames with explicit column order: time, freq, then channels
    column_order = [:time, :freq, channel_labels...]
    power_df = DataFrame(cols_power, copycols = false)[!, column_order]
    phase_df = DataFrame(cols_phase, copycols = false)[!, column_order]
    return (power_df, phase_df)
end

# ============================================================================
# Main tf_analysis Function
# ============================================================================

"""
    tf_analysis(epochs::EpochData, frequencies, time_steps; 
                 method=:wavelet, keep_epochs=false, kwargs...) -> TimeFreqData or TimeFreqEpochData

Perform time-frequency analysis on epoched EEG data.

This is the main time-frequency analysis function for end users.

# Arguments
- `epochs::EpochData`: Epoched EEG data
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (seconds) for TF decomposition

# Keyword Arguments

## Common Options (all methods)
- `method::Symbol=:wavelet`: Analysis method (see Method-Specific Options below)
- `keep_epochs::Bool=false`: Keep individual epochs (true) or average (false)
- `channel_selection::Function=channels()`: Channel selection predicate
- `log_freqs::Bool=true`: Use logarithmic frequency spacing (only used by `:wavelet` and `:superlet`)

## Method-Specific Options

### `:wavelet` - Morlet Wavelet Convolution
- `width`: Number of wavelet cycles
  - `Int` (e.g., `7`): Fixed cycles for all frequencies (default: `7`)
  - `Tuple{Number,Number}`: Frequency-dependent cycles (e.g., `(3, 10)` for log-spaced)
  - `Vector{Number}`: Explicit cycles per frequency

### `:hanning_fixed` - Short-Time FFT with Fixed Window (Cohen Chapter 15)
- `time_window_length::Float64=0.5`: Fixed window length in seconds for all frequencies (default: `0.5`)
- Uses a Hanning taper with the same window length at all frequencies
- Corresponds to FieldTrip's `'mtmconvol'` with `cfg.taper = 'hanning'` and fixed `cfg.t_ftimwin`

### `:hanning_adaptive` - Short-Time FFT with Frequency-Dependent Window (Cohen Chapter 15)
- `cycles::Int=5`: Number of cycles per frequency (default: `5`)
- Window length = cycles / frequency (adaptive)
- Uses a Hanning taper with frequency-dependent window length
- Corresponds to FieldTrip's `'mtmconvol'` with `cfg.taper = 'hanning'` and frequency-dependent `cfg.t_ftimwin`

### `:multitaper` - Multi-Taper Method with DPSS Tapers (Cohen Chapter 16)
- `time_window_length::Float64=0.5`: Time window length in seconds (default: `0.5`)
- `time_bandwidth_product::Float64=3.0`: Time-bandwidth product (nw) - controls number of tapers (default: `3.0`)
  - Higher values = more tapers = lower variance but more frequency smoothing
  - Number of tapers = 2*nw - 1 (e.g., nw=3 gives 5 tapers)
  - Cohen uses nw_product = 3 in Figure 16.2
- Uses DPSS (Slepian) tapers to reduce variance in power estimates
- Corresponds to FieldTrip's `'mtmconvol'` with `cfg.taper = 'dpss'`

# Returns
- `TimeFreqData`: If `keep_epochs=false` (default), returns averaged time-frequency data
- `TimeFreqEpochData`: If `keep_epochs=true`, returns trial-level time-frequency data

# Examples
```julia
# Wavelet method (default)
tf_data = tf_analysis(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, width=7)

# Keep individual epochs
tf_epochs = tf_analysis(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, keep_epochs=true)
```

# Validation
Invalid options for the selected method will raise an `ArgumentError` with a helpful message.
"""
function tf_analysis(
    epochs::EpochData,
    frequencies,
    time_steps;
    method::Symbol = :wavelet,
    keep_epochs::Bool = false,
    channel_selection::Function = channels(),
    kwargs...,
)
    # Validate method-specific options
    _validate_tf_method_options(method, kwargs)

    # Get channels to process
    selected_channels = get_selected_channels(epochs, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected")

    n_channels = length(selected_channels)
    n_trials = length(epochs.data)
    n_samples = nrow(epochs.data[1])
    times = collect(epochs.data[1].time)
    sr = epochs.sample_rate

    frequency_of_interest = collect(Float64, frequencies)
    time_of_interest = collect(Float64, time_steps)

    # Pre-compute common parameters
    log_freqs = get(kwargs, :log_freqs, true)
    n_freqs = length(frequency_of_interest)
    n_times = length(time_of_interest)

    # Output time and frequency vectors
    freqs_out =
        log_freqs && method in [:wavelet, :superlet] ?
        exp.(range(log(float(frequency_of_interest[1])), log(float(frequency_of_interest[end])), length = n_freqs)) :
        collect(Float64, frequency_of_interest)
    # Find indices of time points of interest
    tois_idx = [argmin(abs.(times .- t)) for t in time_of_interest]
    times_out = time_of_interest

    # For now, only implement :wavelet method
    if method == :wavelet
        # Get width parameter
        width_val = get(kwargs, :width, 7)
        if width_val isa Tuple{<:Number,<:Number}
            cycles = exp.(range(log(float(width_val[1])), log(float(width_val[2])), length = n_freqs))
        elseif width_val isa Number
            cycles = fill(float(width_val), n_freqs)
        else
            cycles = collect(Float64, width_val)
        end

        # 1. CONVERT DATA TO 3D ARRAY (Optimized Layout: Samples × Trials × Channels)
        data_3d = zeros(Float64, n_samples, n_trials, n_channels)
        for (ti, trial_df) in enumerate(epochs.data)
            for (ci, ch) in enumerate(selected_channels)
                @inbounds data_3d[:, ti, ci] .= trial_df[!, ch]
            end
        end

        # 2. BATCH PROCESSING (Cohen's Trick)
        n_total = n_samples * n_trials
        data_batch = reshape(data_3d, n_total, n_channels)

        # Determine convolution length for the concatenated signal
        max_sig = maximum(cycles ./ (2π .* freqs_out))
        n_convolution = n_total + ceil(Int, 6 * max_sig * sr)
        n_conv_batch = DSP.nextfastfft(n_convolution)

        # Setup Workspace
        ws = TFMorletWorkspace(n_conv_batch, n_channels, n_trials, n_times)
        # Zero padding and fill signal
        @views ws.padded_signal[1:n_total, :] .= data_batch
        @views ws.padded_signal[n_total+1:end, :] .= 0.0

        # Batch FFT (compute once for all channels)
        mul!(ws.signal_fft, ws.p_fft, ws.padded_signal)

        # 3. ALLOCATE OUTPUT ARRAYS (power and complex)
        if keep_epochs
            channel_results = [zeros(Float64, n_freqs, n_times, n_trials) for _ = 1:n_channels]
            channel_complex = [zeros(ComplexF64, n_freqs, n_times, n_trials) for _ = 1:n_channels]
        else
            full_power = zeros(Float64, n_freqs, n_times, n_channels)
            full_complex = zeros(ComplexF64, n_freqs, n_times, n_channels)
        end

        # 4. PROCESS ALL FREQUENCIES (batch processing)
        _tf_batch_wavelet!(
            keep_epochs ? channel_results : full_power,
            keep_epochs ? channel_complex : full_complex,
            ws,
            n_samples,
            n_trials,
            n_channels,
            tois_idx,
            freqs_out,
            cycles,
            sr,
            keep_epochs,
        )

        # Assemble output data structures
        if keep_epochs
            trial_dfs_power = Vector{DataFrame}(undef, n_trials)
            trial_dfs_phase = Vector{DataFrame}(undef, n_trials)
            for trial_idx = 1:n_trials
                # Temporary matrices for this trial (freqs x times x channels)
                trial_power = zeros(n_freqs, n_times, n_channels)
                trial_complex = zeros(ComplexF64, n_freqs, n_times, n_channels)
                for ci = 1:n_channels
                    @views trial_power[:, :, ci] .= channel_results[ci][:, :, trial_idx]
                    @views trial_complex[:, :, ci] .= channel_complex[ci][:, :, trial_idx]
                end
                power_df, phase_df = _power_to_tf_dataframe(trial_complex, times_out, freqs_out, selected_channels)
                trial_dfs_power[trial_idx] = power_df
                trial_dfs_phase[trial_idx] = phase_df
            end

            return TimeFreqEpochData(
                epochs.file,
                epochs.condition,
                epochs.condition_name,
                trial_dfs_power,
                trial_dfs_phase,
                epochs.layout,
                epochs.sample_rate,
                method,
                nothing,
                epochs.analysis_info,
            )
        else
            power_df, phase_df = _power_to_tf_dataframe(full_complex, times_out, freqs_out, selected_channels)

            return TimeFreqData(
                epochs.file,
                epochs.condition,
                epochs.condition_name,
                power_df,
                phase_df,
                epochs.layout,
                epochs.sample_rate,
                method,
                nothing,
                epochs.analysis_info,
            )
        end
    elseif method in [:hanning_fixed, :hanning_adaptive]
        # Short-Time FFT (Cohen Chapter 15)
        # Process each channel separately
        if keep_epochs
            channel_results = [zeros(Float64, n_freqs, n_times, n_trials) for _ = 1:n_channels]
            channel_complex = [zeros(ComplexF64, n_freqs, n_times, n_trials) for _ = 1:n_channels]
        else
            full_power = zeros(Float64, n_freqs, n_times, n_channels)
            full_complex = zeros(ComplexF64, n_freqs, n_times, n_channels)
        end

        # Get method-specific parameters
        h_win_len = method == :hanning_fixed ? get(kwargs, :time_window_length, 0.5) : nothing
        h_cycles = method == :hanning_adaptive ? get(kwargs, :cycles, 5) : nothing

        # Pre-allocate workspace if using fixed window
        h_ws = nothing
        if method == :hanning_fixed && !isnothing(h_win_len)
            timewinidx = round(Int, h_win_len * sr)
            timewinidx += iseven(timewinidx)
            h_ws = TFHanningFixedWorkspace(timewinidx, n_trials)
        end

        for ci = 1:n_channels
            # Extract channel data (samples × trials)
            signal_ci = zeros(Float64, n_samples, n_trials)
            for (ti, trial_df) in enumerate(epochs.data)
                @inbounds signal_ci[:, ti] .= trial_df[!, selected_channels[ci]]
            end

            # Run Short-Time FFT
            power, complex_data, _, _ = tf_hanning(
                signal_ci,
                times,
                sr,
                freqs_out,
                time_of_interest;
                window_length = h_win_len,
                cycles = h_cycles,
                keeptrials = keep_epochs,
                workspace = h_ws,
            )

            if keep_epochs
                channel_results[ci] .= ndims(power) == 2 ? reshape(power, size(power)..., 1) : power
                channel_complex[ci] .= ndims(complex_data) == 2 ? reshape(complex_data, size(complex_data)..., 1) : complex_data
            else
                full_power[:, :, ci] .= power
                full_complex[:, :, ci] .= complex_data
            end
        end

        # Assemble output data structures
        if keep_epochs
            trial_dfs_power = Vector{DataFrame}(undef, n_trials)
            trial_dfs_phase = Vector{DataFrame}(undef, n_trials)
            for trial_idx = 1:n_trials
                # Temporary matrices for this trial (freqs x times x channels)
                trial_power = zeros(n_freqs, n_times, n_channels)
                trial_complex = zeros(ComplexF64, n_freqs, n_times, n_channels)
                for ci = 1:n_channels
                    @views trial_power[:, :, ci] .= channel_results[ci][:, :, trial_idx]
                    @views trial_complex[:, :, ci] .= channel_complex[ci][:, :, trial_idx]
                end
                power_df, phase_df = _power_to_tf_dataframe(trial_complex, times_out, freqs_out, selected_channels)
                trial_dfs_power[trial_idx] = power_df
                trial_dfs_phase[trial_idx] = phase_df
            end

            return TimeFreqEpochData(
                epochs.file,
                epochs.condition,
                epochs.condition_name,
                trial_dfs_power,
                trial_dfs_phase,
                epochs.layout,
                epochs.sample_rate,
                method,
                nothing,
                epochs.analysis_info,
            )
        else
            power_df, phase_df = _power_to_tf_dataframe(full_complex, times_out, freqs_out, selected_channels)

            return TimeFreqData(
                epochs.file,
                epochs.condition,
                epochs.condition_name,
                power_df,
                phase_df,
                epochs.layout,
                epochs.sample_rate,
                method,
                nothing,
                epochs.analysis_info,
            )
        end
    elseif method == :multitaper
        # Multi-Taper Method (Cohen Chapter 16)
        # Process each channel separately
        if keep_epochs
            channel_results = [zeros(Float64, n_freqs, n_times, n_trials) for _ = 1:n_channels]
            channel_complex = [zeros(ComplexF64, n_freqs, n_times, n_trials) for _ = 1:n_channels]
        else
            full_power = zeros(Float64, n_freqs, n_times, n_channels)
            full_complex = zeros(ComplexF64, n_freqs, n_times, n_channels)
        end

        # Get method-specific parameters
        mt_win_len = get(kwargs, :time_window_length, 0.5)
        mt_nw = get(kwargs, :time_bandwidth_product, 3.0)

        # Pre-allocate workspace
        timewinidx = round(Int, mt_win_len * sr)
        timewinidx += iseven(timewinidx)
        mt_ws = TFMultitaperWorkspace(timewinidx, n_trials, n_freqs)

        for ci = 1:n_channels
            # Extract channel data (samples × trials)
            signal_ci = zeros(Float64, n_samples, n_trials)
            for (ti, trial_df) in enumerate(epochs.data)
                @inbounds signal_ci[:, ti] .= trial_df[!, selected_channels[ci]]
            end

            # Run Multi-Taper analysis
            power, complex_data, _, _ = tf_multitaper(
                signal_ci,
                times,
                sr,
                freqs_out,
                time_of_interest;
                time_window_length = mt_win_len,
                time_bandwidth_product = mt_nw,
                keeptrials = keep_epochs,
                workspace = mt_ws,
            )

            if keep_epochs
                channel_results[ci] .= ndims(power) == 2 ? reshape(power, size(power)..., 1) : power
                channel_complex[ci] .= ndims(complex_data) == 2 ? reshape(complex_data, size(complex_data)..., 1) : complex_data
            else
                full_power[:, :, ci] .= power
                full_complex[:, :, ci] .= complex_data
            end
        end

        # Assemble output data structures
        if keep_epochs
            trial_dfs_power = Vector{DataFrame}(undef, n_trials)
            trial_dfs_phase = Vector{DataFrame}(undef, n_trials)
            for trial_idx = 1:n_trials
                # Temporary matrices for this trial (freqs x times x channels)
                trial_power = zeros(n_freqs, n_times, n_channels)
                trial_complex = zeros(ComplexF64, n_freqs, n_times, n_channels)
                for ci = 1:n_channels
                    @views trial_power[:, :, ci] .= channel_results[ci][:, :, trial_idx]
                    @views trial_complex[:, :, ci] .= channel_complex[ci][:, :, trial_idx]
                end
                power_df, phase_df = _power_to_tf_dataframe(trial_complex, times_out, freqs_out, selected_channels)
                trial_dfs_power[trial_idx] = power_df
                trial_dfs_phase[trial_idx] = phase_df
            end

            return TimeFreqEpochData(
                epochs.file,
                epochs.condition,
                epochs.condition_name,
                trial_dfs_power,
                trial_dfs_phase,
                epochs.layout,
                epochs.sample_rate,
                method,
                nothing,
                epochs.analysis_info,
            )
        else
            power_df, phase_df = _power_to_tf_dataframe(full_complex, times_out, freqs_out, selected_channels)

            return TimeFreqData(
                epochs.file,
                epochs.condition,
                epochs.condition_name,
                power_df,
                phase_df,
                epochs.layout,
                epochs.sample_rate,
                method,
                nothing,
                epochs.analysis_info,
            )
        end
    else
        error("Method `:$method` not yet implemented. Currently supported: `:wavelet`, `:hanning_fixed`, `:hanning_adaptive`, `:multitaper`")
    end
end
