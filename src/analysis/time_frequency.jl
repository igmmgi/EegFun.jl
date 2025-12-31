"""
Time-Frequency Analysis Module

Provides FieldTrip-style time-frequency analysis with multiple methods.

Main function for end users:
- `tf_analysis(epochs, ...)`: High-level function for EpochData → TimeFreqData

Low-level functions for raw signals:
- `tf_wavelet`: Morlet wavelet convolution
- `tf_multitaper`: Multi-taper sliding window analysis (Hanning or DPSS)
- `tf_superlet`: Superlet transform (super-resolution)
- `tf_spectrum`: Multi-taper FFT (non-time-resolved power spectrum)

# Examples
```julia
# High-level: EpochData → TimeFreqData (recommended)
tf_data = tf_analysis(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, width=7)

# Low-level: raw signals → power arrays
power, t, f = tf_wavelet(signal, times, sr, 2:30, -0.5:0.05:1.5; width=7)
power, t, f = tf_multitaper(signal, times, sr, 2:30, -0.5:0.05:1.5; taper=:dpss)
power, f = tf_spectrum(signal, sr, 2:40; taper=:dpss)
```
"""

# ============================================================================
# DPSS Taper (Slepian sequences)
# ============================================================================

"""
    dpss_tapers(n::Int, nw::Real, k::Int)

Compute Discrete Prolate Spheroidal Sequences (DPSS/Slepian tapers).

# Arguments
- `n::Int`: Window length
- `nw::Real`: Time-bandwidth product (typically 2-4)
- `k::Int`: Number of tapers (typically 2*nw - 1)

# Returns
- Matrix of size (n, k) containing the tapers
"""
function dpss_tapers(n::Int, nw::Real, k::Int)
    return DSP.dpss(n, nw, k)
end

# ============================================================================
# Shared Helper Functions
# ============================================================================

function _prepare_morlet_wavelets(frequencies, width, sample_rate, n_samples; log_freqs=true)
    freqs = collect(Float64, frequencies)
    n_freqs = length(freqs)
    
    # 1. Frequency vector
    freqs_out = log_freqs ? exp.(range(log(float(freqs[1])), log(float(freqs[end])), length = n_freqs)) : freqs
    
    # 2. Cycles (width)
    if width isa Tuple{<:Number, <:Number}
        # Variable cycles: interpolate
        cycles = exp.(range(log(float(width[1])), log(float(width[2])), length = n_freqs))
    elseif width isa Number
        cycles = fill(float(width), n_freqs)
    else
        cycles = collect(Float64, width)
    end
    
    # 3. Sigmas and Window lengths
    sigmas = cycles ./ (2 * pi .* freqs_out)
    half_wins = Int[]
    for sig in sigmas
        wl = ceil(Int, 6 * sig * sample_rate)
        wl += iseven(wl)
        push!(half_wins, wl ÷ 2)
    end
    
    # 4. Convolution length
    # Use nextfastfft for optimal performance (handles factors 2, 3, 5, 7)
    # instead of restrictive nextpow(2)
    n_convolution = n_samples + ceil(Int, 6 * maximum(sigmas) * sample_rate)
    n_conv = DSP.nextfastfft(n_convolution)
    
    # 5. Pre-compute wavelet FFTs
    wavelet_ffts = zeros(ComplexF64, n_conv, n_freqs)
    tmp_w_padded = zeros(ComplexF64, n_conv)
    
    for fi in 1:n_freqs
        freq = freqs_out[fi]
        sigma = sigmas[fi]
        hw = half_wins[fi]
        wl = hw * 2 + 1 
        
        fill!(tmp_w_padded, 0.0)
        w_time = range(-wl / 2, wl / 2, length = wl) ./ sample_rate
        @. tmp_w_padded[1:wl] = sqrt(1 / (sigma * sqrt(pi))) * exp(2im * pi * freq * w_time) * exp(-w_time^2 / (2 * sigma^2))
        wavelet_ffts[:, fi] = fft(tmp_w_padded)
    end
    
    return freqs_out, half_wins, n_conv, wavelet_ffts
end

# ============================================================================
# Option 1: Individual TF functions with keyword arguments
# ============================================================================

"""
    tf_wavelet(signal, times, sample_rate, frequencies, time_steps; 
               width=7, keeptrials=false, log_freqs=true) -> (power, times_out, freqs_out)

Time-frequency analysis using Morlet wavelets.

# Arguments
- `signal`: Signal matrix (samples × trials) or vector
- `times`: Time vector corresponding to signal samples
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (seconds)

# Keyword Arguments
- `width`: Number of wavelet cycles. Can be:
  - `Int` (e.g., `7`): Fixed cycles for all frequencies
  - `Tuple{Int,Int}` (e.g., `(3, 10)`): Variable cycles, linearly interpolated from low to high freq
- `keeptrials::Bool=false`: Keep individual trials (true) or average (false)
- `log_freqs::Bool=true`: Use logarithmically spaced frequencies

# Returns
- `power`: Time-frequency power (freqs × times) or (freqs × times × trials)
- `times_out`: Output time vector
- `freqs_out`: Output frequency vector

# Example
```julia
# Fixed 7 cycles for all frequencies
power, t, f = tf_wavelet(signal, times, 256, 2:30, -0.5:0.05:1.5; width=7)

# Variable cycles: 3 at lowest freq → 10 at highest freq (better time-freq tradeoff)
power, t, f = tf_wavelet(signal, times, 256, 2:30, -0.5:0.05:1.5; width=(3, 10))
```
"""
function tf_wavelet(signal, times, sample_rate, frequencies, time_steps; 
                    width::Union{Int, Tuple{Int,Int}, Vector{Int}}=7, 
                    keeptrials::Bool=false, log_freqs::Bool=true)
    
    # Ensure signal is 2D
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)
    
    # Prepare parameters and wavelets
    freqs_out, half_wins, n_conv, wavelet_ffts = _prepare_morlet_wavelets(frequencies, width, sample_rate, n_samples; log_freqs=log_freqs)
    
    # Setup output time indices
    toi = collect(Float64, time_steps)
    tois_idx = [argmin(abs.(times .- t)) for t in toi]
    
    # Create workspace
    ws = TFMorletWorkspace(n_conv, n_trials)
    
    # Convert signal to 3D for kernel (1 channel)
    data_3d = reshape(signal, n_samples, 1, n_trials)
    
    # Allocate output
    n_freqs = length(freqs_out)
    n_times = length(toi)
    
    if keeptrials
        power = zeros(n_freqs, n_times, n_trials)
        tf_morlet_optimized!(power, ws, data_3d, 1, n_samples, n_trials, tois_idx, half_wins, wavelet_ffts, true)
    else
        power = zeros(n_freqs, n_times)
        tf_morlet_optimized!(power, ws, data_3d, 1, n_samples, n_trials, tois_idx, half_wins, wavelet_ffts, false)
    end
    
    return power, toi, freqs_out
end

"""
    tf_superlet(signal, times, sample_rate, frequencies, time_steps;
                order=5, width=3, combine=:multiplicative, keeptrials=false)
                -> (power, times_out, freqs_out)

Time-frequency analysis using Superlets (Moca et al. 2021, Nature Communications).

Superlets combine multiple Morlet wavelets with progressively increasing cycle widths 
via geometric mean, providing super-resolution in time-frequency space by overcoming 
the Heisenberg uncertainty trade-off.

# Arguments
- `signal`: Signal matrix (samples × trials) or vector
- `times`: Time vector corresponding to signal samples
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (seconds)

# Keyword Arguments
- `order`: Superlet order (number of wavelets to combine). Can be:
  - `Int` (e.g., `5`): Same order for all frequencies
  - `Vector{Int}`: Different order per frequency (higher for higher frequencies)
- `width::Int=3`: Base wavelet width (cycles) - typically lower than standard wavelet (3 vs 7)
- `combine::Symbol=:multiplicative`: How cycle widths progress
  - `:multiplicative`: cycles = width * [1, 2, 3, ..., order]
  - `:additive`: cycles = [width, width+1, width+2, ..., width+order-1]
- `keeptrials::Bool=false`: Keep individual trials (true) or average (false)

# Returns
- `power`: Time-frequency power (freqs × times) or (freqs × times × trials)
- `times_out`: Output time vector
- `freqs_out`: Output frequency vector

# Example
```julia
# Basic superlet with order 5 (default)
power, t, f = tf_superlet(signal, times, 256, 2:30, -0.5:0.05:1.5)

# Higher order for better resolution (more wavelets combined)
power, t, f = tf_superlet(signal, times, 256, 2:30, -0.5:0.05:1.5; order=10)

# Frequency-adaptive order: higher order for higher frequencies
power, t, f = tf_superlet(signal, times, 256, 2:30, -0.5:0.05:1.5; order=collect(3:31))
```

# Reference
Moca VV, Bârzan H, Nagy-Dăbâcan A, Mureșan RC (2021). Time-frequency super-resolution 
with superlets. Nature Communications 12:337. https://doi.org/10.1038/s41467-020-20539-9
"""


function tf_superlet(signal, times, sample_rate, frequencies, time_steps;
                     order::Union{Int, Vector{Int}}=5, width::Int=3,
                     combine::Symbol=:multiplicative, keeptrials::Bool=false)
    
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)
    
    freqs = collect(Float64, frequencies)
    toi = collect(Float64, time_steps)
    n_freqs = length(freqs)
    
    orders = order isa Int ? fill(order, n_freqs) : order
    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in toi]
    
    max_sigma = (width * maximum(orders)) / (2π * minimum(freqs))
    n_conv = nextpow(2, n_samples + ceil(Int, 6 * max_sigma * sample_rate))
    
    ws = TFMorletWorkspace(n_conv, n_trials)
    
    if keeptrials
        power = zeros(n_freqs, length(toi), n_trials)
        tf_superlet_optimized!(power, ws, signal, 1, n_samples, n_trials, tois_idx, freqs, orders, width, combine, true, sample_rate)
        return power, toi, freqs
    else
        power = zeros(n_freqs, length(toi))
        tf_superlet_optimized!(power, ws, signal, 1, n_samples, n_trials, tois_idx, freqs, orders, width, combine, false, sample_rate)
        return power, toi, freqs
    end
end




"""
    tf_multitaper(signal, times, sample_rate, frequencies, time_steps;

                  taper=:hanning, t_ftimwin=0.5, tapsmofrq=4.0, keeptrials=false) 
                  -> (power, times_out, freqs_out)

Time-frequency analysis using sliding window with tapers (Hanning or DPSS).

# Arguments
- `signal`: Signal matrix (samples × trials) or vector
- `times`: Time vector corresponding to signal samples
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (seconds)

# Keyword Arguments
- `taper::Symbol=:hanning`: Taper type (`:hanning` or `:dpss`)
- `t_ftimwin::Float64=0.5`: Time window length in seconds
- `tapsmofrq::Float64=4.0`: Frequency smoothing in Hz (for `:dpss` only)
- `keeptrials::Bool=false`: Keep individual trials (true) or average (false)

# Returns
- `power`: Time-frequency power (freqs × times) or (freqs × times × trials)
- `times_out`: Output time vector
- `freqs_out`: Output frequency vector

# Example
```julia
# With Hanning taper
power, t, f = tf_multitaper(signal, times, 256, 2:30, -0.5:0.05:1.5; taper=:hanning)

# With DPSS tapers (multi-taper)
power, t, f = tf_multitaper(signal, times, 256, 2:30, -0.5:0.05:1.5; taper=:dpss, tapsmofrq=4.0)
```
"""
function tf_multitaper(signal, times, sample_rate, frequencies, time_steps;
                       taper::Symbol=:hanning, t_ftimwin::Float64=0.5, 
                       tapsmofrq::Float64=4.0, keeptrials::Bool=false)
    freqs = collect(Float64, frequencies)
    toi = collect(Float64, time_steps)
    
    if taper == :hanning
        power, times_out, freqs_out = tf_hanning(signal, times, sample_rate, freqs, toi; 
                                                  window_length=t_ftimwin, keeptrials=keeptrials)
    elseif taper == :dpss
        power, times_out, freqs_out = tf_mtm(signal, times, sample_rate, freqs, toi;
                                              t_ftimwin=t_ftimwin, tapsmofrq=tapsmofrq, keeptrials=keeptrials)
        if !keeptrials && ndims(power) == 3
             power = dropdims(mean(power, dims=3), dims=3)
        end
    else
        throw(ArgumentError("Unknown taper: $taper. Use :hanning or :dpss"))
    end
    
    return power, times_out, freqs_out
end

"""
    tf_spectrum(signal, sample_rate, frequencies; 
                taper=:dpss, tapsmofrq=4.0, keeptrials=false) -> (power, freqs_out)

Non-time-resolved power spectrum using multi-taper method.

# Arguments
- `signal`: Signal matrix (samples × trials) or vector
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)

# Keyword Arguments
- `taper::Symbol=:dpss`: Taper type (`:hanning` or `:dpss`)
- `tapsmofrq::Float64=4.0`: Frequency smoothing in Hz (for `:dpss`)
- `keeptrials::Bool=false`: Keep individual trials (true) or average (false)

# Returns
- `power`: Power spectrum (freqs,) or (freqs × trials)
- `freqs_out`: Frequency vector

# Example
```julia
power, freqs = tf_spectrum(signal, 256, 2:40; taper=:dpss, tapsmofrq=2.0)
```
"""
function tf_spectrum(signal, sample_rate, frequencies; 
                     taper::Symbol=:dpss, tapsmofrq::Float64=4.0, keeptrials::Bool=false)
    freqs = collect(Float64, frequencies)
    power, freqs_out = tf_mtmfft(signal, sample_rate, freqs; taper=taper, tapsmofrq=tapsmofrq)
    
    if !keeptrials && ndims(power) == 2
        power = vec(mean(power, dims=2))
    end
    
    return power, freqs_out
end





# ============================================================================
# Multi-taper FFT (non-time-resolved)
# ============================================================================

"""
    tf_mtmfft(signal, sample_rate, frequencies; taper=:dpss, tapsmofrq=4.0)

Multi-taper frequency analysis (non-time-resolved power spectrum).

# Arguments
- `signal`: Signal matrix (samples × trials) or vector
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest
- `taper`: Taper type (:hanning or :dpss)
- `tapsmofrq`: Frequency smoothing in Hz (for :dpss)

# Returns
- `power`: Power spectrum (freqs × trials)
- `frequencies`: Frequency vector
"""
function tf_mtmfft(signal, sample_rate, frequencies; taper::Symbol=:dpss, tapsmofrq::Float64=4.0)
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)
    n_freqs = length(frequencies)
    
    # Compute frequency indices
    freq_resolution = sample_rate / n_samples
    freq_indices = round.(Int, frequencies ./ freq_resolution) .+ 1
    
    # Generate tapers
    if taper == :hanning
        tapers = reshape(0.5 * (1 .- cos.(2π * (0:n_samples-1) / (n_samples - 1))), :, 1)
        n_tapers = 1
    elseif taper == :dpss
        # Time-bandwidth product: nw = tapsmofrq * n_samples / sample_rate / 2
        nw = tapsmofrq * n_samples / sample_rate / 2
        n_tapers = max(1, floor(Int, 2 * nw - 1))
        tapers = dpss_tapers(n_samples, nw, n_tapers)
    else
        throw(ArgumentError("Unknown taper: $taper"))
    end
    
    power = zeros(n_freqs, n_trials)
    
    @inbounds for trial in 1:n_trials
        taper_power = zeros(n_freqs)
        for t in 1:n_tapers
            tapered_signal = @view(signal[:, trial]) .* @view(tapers[:, t])
            spectrum = fft(tapered_signal)
            taper_power .+= abs2.(@view(spectrum[freq_indices]))
        end
        power[:, trial] = taper_power ./ n_tapers
    end
    
    return power, collect(Float64, frequencies)
end


# ============================================================================
# Signal generation for testing
# ============================================================================

"""
    generate_signal(n_trials, time_window, sample_rate, signal_freqs, signal_amps, signal_times, noise)

Generate synthetic signals with known frequency components for testing TF analysis.

# Arguments
- `n_trials`: Number of trials to generate
- `time_window`: [start_time, end_time] in seconds
- `sample_rate`: Sampling rate in Hz
- `signal_freqs`: Vector of frequencies for each component
- `signal_amps`: Vector of amplitudes for each component
- `signal_times`: Vector of [start_time, end_time] for each component
- `noise`: Amplitude of random noise

# Returns
- `times`: Time vector
- `signal`: Matrix (samples × trials)
"""
function generate_signal(

    n_trials::Int,
    time_window::Vector{Float64},
    sample_rate::Real,
    signal_freqs::Vector{<:Real},
    signal_amps::Vector{<:Real},
    signal_times::Vector{Vector{Float64}},
    noise::Real,
)
    # Pre-compute time vector
    dt = 1 / sample_rate
    x_time = range(time_window[1], time_window[2] - dt, step = dt)
    n_samples = length(x_time)

    # Pre-allocate output matrix with noise
    signal_trials = zeros(n_samples, n_trials)
    if noise > 0
        @views for trial = 1:n_trials
            signal_trials[:, trial] .= randn(n_samples) .* noise
        end
    end

    # Pre-allocate temporary arrays
    temp_signal = zeros(n_samples)

    # Add frequency components to all trials
    @inbounds for trial = 1:n_trials
        # Reset temporary signal array
        fill!(temp_signal, 0.0)

        # Add each frequency component
        for (freq, amp, times) in zip(signal_freqs, signal_amps, signal_times)
            # Find time window indices
            time_mask = (x_time .>= times[1]) .& (x_time .< times[2])
            window_times = @view x_time[time_mask]

            # Generate signal for this component
            @. temp_signal[time_mask] += amp * sin(2π * freq * window_times)
        end

        # Add to output
        signal_trials[:, trial] .+= temp_signal
    end

    return x_time, signal_trials
end


function average_over_trials(x::Array{<:AbstractFloat,3})
    return dropdims(mean(x, dims = 3), dims = 3)  # Returns a freq × time matrix
end

"""
    signal_to_epochs(times, signal, channel_name::Symbol, sample_rate::Int;
                     file::String="synthetic", condition::Int=1, condition_name::String="test") -> EpochData

Convert signal data from `generate_signal` to `EpochData` format.

This function takes the output of `generate_signal` (a time vector and a signal matrix)
and converts it to `EpochData`, which can then be used with `tf_analysis` and other
high-level analysis functions.

# Arguments
- `times`: Time vector (from `generate_signal`)
- `signal`: Signal matrix (samples × trials) from `generate_signal`
- `channel_name::Symbol`: Name for the channel (e.g., `:Channel1`)
- `sample_rate::Int`: Sampling rate in Hz

# Keyword Arguments
- `file::String="synthetic"`: Source filename
- `condition::Int=1`: Condition number
- `condition_name::String="test"`: Condition name

# Returns
- `EpochData`: Epoch data ready for time-frequency analysis

# Example
```julia
times, signal = generate_signal(100, [-1.0, 2.0], 256, [10.0], [1.0], [[0.0, 1.0]], 0.5)
epochs = signal_to_epochs(times, signal, :Channel1, 256)
tf_data = tf_analysis(epochs, 2:40, -0.5:0.02:1.5; method=:wavelet)
```
"""
function signal_to_epochs(times, signal, channel_name::Symbol, sample_rate::Int;
                          file::String="synthetic", condition::Int=1, condition_name::String="test")
    # Ensure signal is 2D (samples × trials)
    if ndims(signal) == 1
        signal = reshape(signal, :, 1)
    end
    
    n_samples, n_trials = size(signal)
    
    # Convert times to vector if needed
    times_vec = collect(Float64, times)
    @assert length(times_vec) == n_samples "Time vector length ($(length(times_vec))) must match signal samples ($n_samples)"
    
    # Create Vector{DataFrame} - one DataFrame per trial
    trial_dfs = Vector{DataFrame}(undef, n_trials)
    
    for trial_idx in 1:n_trials
        # Create DataFrame for this trial with time and channel columns
        trial_dfs[trial_idx] = DataFrame(
            :time => times_vec,
            channel_name => signal[:, trial_idx],
            :epoch => trial_idx
        )
    end
    
    # Create minimal layout for the channel
    layout_df = DataFrame(:label => [channel_name])
    layout = Layout(layout_df, nothing, nothing)
    
    # Create AnalysisInfo
    analysis_info = AnalysisInfo()
    
    # Create and return EpochData
    return EpochData(
        file,
        condition,
        condition_name,
        trial_dfs,
        layout,
        sample_rate,
        analysis_info
    )
end

# baseline operations
function apply_tf_baseline_db(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(@view(tf_data[:, base_idx]), dims = 2)
    return 10 .* log10.(tf_data ./ baseline_power)
end

function apply_tf_baseline_perchange(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(@view(tf_data[:, base_idx]), dims = 2)
    return 100 .* (tf_data .- baseline_power) ./ baseline_power
end

function apply_tf_baseline_relchange(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(@view(tf_data[:, base_idx]), dims = 2)
    return tf_data ./ baseline_power
end

function apply_tf_baseline_ztransforms(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = @view tf_data[:, base_idx]
    baseline_mean = mean(baseline_power, dims = 2)
    baseline_std = std(baseline_power, dims = 2)
    return (tf_data .- baseline_mean) ./ baseline_std
end


# Optimized internal workspace for Morlet transform to avoid allocations
struct TFMorletWorkspace
    padded_signal::Matrix{ComplexF64}
    signal_fft::Matrix{ComplexF64}
    conv_fft::Matrix{ComplexF64}
    eegconv::Matrix{ComplexF64}
    p_fft::Any
    p_ifft::Any
end

function TFMorletWorkspace(n_conv, n_trials)
    ps = zeros(ComplexF64, n_conv, n_trials)
    sf = zeros(ComplexF64, n_conv, n_trials)
    cf = zeros(ComplexF64, n_conv, n_trials)
    ec = zeros(ComplexF64, n_conv, n_trials)
    # Use FFTW.ESTIMATE for faster planning (MEASURE is too slow for per-thread workspace creation)
    p_fft = plan_fft(ps, 1, flags=FFTW.ESTIMATE)
    p_ifft = plan_ifft(cf, 1, flags=FFTW.ESTIMATE)
    return TFMorletWorkspace(ps, sf, cf, ec, p_fft, p_ifft)
end

function tf_morlet_optimized!(output_power, ws::TFMorletWorkspace, data_3d, ci, 
                              n_samples, n_trials, tois_idx, half_wins, 
                              wavelet_ffts, keeptrials)

    n_freq = length(half_wins)
    n_times = length(tois_idx)
    n_conv = size(ws.padded_signal, 1)

    # 1. Prepare signal and compute signal FFT once per channel
    # This is often the memory bottleneck if not careful
    for trial in 1:n_trials
        @inbounds for s in 1:n_samples
            ws.padded_signal[s, trial] = ComplexF64(data_3d[s, ci, trial], 0.0)
        end
        @inbounds for s in (n_samples+1):n_conv
            ws.padded_signal[s, trial] = ComplexF64(0.0, 0.0)
        end
    end
    mul!(ws.signal_fft, ws.p_fft, ws.padded_signal)

    # 2. Iterate frequencies (Reverting Mega-Batch for better cache/memory profile)
    inv_n_trials = 1.0 / n_trials
    
    @inbounds for fi in 1:n_freq
        half_win = half_wins[fi]
        valid_start = half_win + 1
        
        # Get precomputed wavelet FFT
        @views w_fft = wavelet_ffts[:, fi]
        
        # Convolve: conv_fft = signal_fft * wavelet_fft
        # Manual loop for absolute zero allocations and SIMD
        for trial in 1:n_trials
            @simd for i in 1:n_conv
                ws.conv_fft[i, trial] = ws.signal_fft[i, trial] * w_fft[i]
            end
        end
        
        # IFFT back to time domain
        mul!(ws.eegconv, ws.p_ifft, ws.conv_fft)
        
        # 4. Extract power at TOIs
        if keeptrials
            for trial in 1:n_trials
                for ti in 1:n_times
                    t_idx = tois_idx[ti]
                    output_power[fi, ti, trial] = abs2(ws.eegconv[valid_start + t_idx - 1, trial])
                end
            end
        else
            # Average directly
            for trial in 1:n_trials
                for ti in 1:n_times
                    t_idx = tois_idx[ti]
                    output_power[fi, ti] += abs2(ws.eegconv[valid_start + t_idx - 1, trial]) * inv_n_trials
                end
            end
        end
    end
end






function tf_superlet_optimized!(output_power, ws::TFMorletWorkspace, signal, ci, # signal can be 2D (samples x trials) or 3D slice if we adapt
                                n_samples, n_trials, tois_idx, freqs, orders, width, combine, keeptrials, sample_rate)

    n_freqs = length(freqs)
    n_times = length(tois_idx)
    n_conv = size(ws.padded_signal, 1)
    
    # 1. Prepare signal
    # If signal is 3D (samples x channels x trials), ci is index.
    # If signal is 2D (samples x trials), ci is ignored or used as offset?
    # Let's enforce signal to be the 2D matrix for simplicity in this helper
    for trial in 1:n_trials
        @inbounds for s in 1:n_samples
            ws.padded_signal[s, trial] = ComplexF64(signal[s, trial], 0.0)
        end
        @inbounds for s in (n_samples+1):n_conv
            ws.padded_signal[s, trial] = ComplexF64(0.0, 0.0)
        end
    end
    mul!(ws.signal_fft, ws.p_fft, ws.padded_signal)

    wavelet_padded = zeros(ComplexF64, n_conv)
    # Buffers for geometric mean accumulation
    prod_power = zeros(Float64, n_times, n_trials)
    
    inv_n_trials = 1.0 / n_trials
    
    @inbounds for fi in 1:n_freqs
        freq = freqs[fi]
        ord = orders[fi]
        cycles_list = combine == :multiplicative ? [width * i for i in 1:ord] : [width + i - 1 for i in 1:ord]
        
        fill!(prod_power, 1.0)
        
        for cyc in cycles_list
            sigma = cyc / (2π * freq)
            win_len = ceil(Int, 6 * sigma * sample_rate)
            win_len += iseven(win_len)
            half_win = win_len ÷ 2
            
            fill!(wavelet_padded, 0.0)
            w_time = range(-win_len / 2, win_len / 2, length = win_len) ./ sample_rate
            @. wavelet_padded[1:win_len] = sqrt(1 / (sigma * sqrt(pi))) * exp(2im * pi * freq * w_time) * exp(-w_time^2 / (2 * sigma^2))
            
            # FFT of wavelet (allocating a bit here, but negligible vs IFFT)
            wavelet_fft = fft(wavelet_padded)
            
            # Convolution
            for trial in 1:n_trials
                @simd for i in 1:n_conv
                    ws.conv_fft[i, trial] = ws.signal_fft[i, trial] * wavelet_fft[i]
                end
            end
            mul!(ws.eegconv, ws.p_ifft, ws.conv_fft)
            
            valid_start = half_win + 1
            for trial in 1:n_trials
                for ti in 1:n_times
                    t_idx = tois_idx[ti]
                    prod_power[ti, trial] *= abs2(ws.eegconv[valid_start + t_idx - 1, trial])
                end
            end
        end
        
        # Power calculation (Geometric Mean)
        inv_ord = 1.0 / ord
        if keeptrials
             for trial in 1:n_trials
                @simd for ti in 1:n_times
                    output_power[fi, ti, trial] = prod_power[ti, trial] ^ inv_ord
                end
            end
        else
            # Average directly
            for trial in 1:n_trials
                @simd for ti in 1:n_times
                    output_power[fi, ti] += (prod_power[ti, trial] ^ inv_ord) * inv_n_trials
                end
            end
        end
    end
end



function tf_hanning(signal, times, sample_rate, frequencies, time_steps; window_length = nothing, cycles = nothing, keeptrials::Bool=false)

    if isnothing(window_length) && isnothing(cycles)
        throw(ArgumentError("Must specify either window_length or cycles"))
    end

    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)

    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_steps]
    n_frex = length(frequencies)
    n_timepoints = length(time_steps)
    
    # Optimization: If window_length is fixed (cycles=nothing), we can do 1 FFT per time point
    # and extract all frequencies, similar to MTM. This avoids n_freqs * n_times FFTs.
    if isnothing(cycles)
        timewinidx = round(Int, window_length * sample_rate)
        timewinidx += iseven(timewinidx)
        half_win = timewinidx ÷ 2
        
        hann_win = 0.5 * (1 .- cos.(2π * (0:(timewinidx-1)) / (timewinidx - 1)))
        
        # Frequency indices for all frequencies
        freq_indices = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
        
        # Pre-allocate buffers (Real input for rfft)
        tmpdat = zeros(Float64, timewinidx, n_trials)
        p_rfft = plan_rfft(tmpdat, 1) # Plan once
        
        # Output size of rfft is div(n, 2) + 1
        n_rfft_out = div(timewinidx, 2) + 1
        fdat = zeros(ComplexF64, n_rfft_out, n_trials)
        
        # Allocate output
        if keeptrials
            tf_trials = zeros(n_frex, n_timepoints, n_trials)
        else
            tf_avg = zeros(n_frex, n_timepoints)
        end
        inv_n_trials = 1.0 / n_trials
        
        @inbounds for (ti, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win

            # Extract signal and apply window in one pass (Loop Fusion)
            if start_idx >= 1 && end_idx <= n_samples
                for trial in 1:n_trials
                    @simd for s in 1:timewinidx
                        tmpdat[s, trial] = signal[start_idx + s - 1, trial] * hann_win[s]
                    end
                end
            else
                fill!(tmpdat, 0.0)
                v_start = max(1, start_idx)
                v_end = min(n_samples, end_idx)
                offset = v_start - start_idx
                len = v_end - v_start + 1
                for trial in 1:n_trials
                    @simd for s in 1:len
                        tmpdat[offset + s, trial] = signal[v_start + s - 1, trial] * hann_win[offset + s]
                    end
                end
            end

            # Real-to-Complex FFT (2x faster than Complex FFT)
            mul!(fdat, p_rfft, tmpdat)
            
            # Extract power
            if keeptrials
                 for trial in 1:n_trials
                    for (fi, f_idx) in enumerate(freq_indices)
                        tf_trials[fi, ti, trial] = abs2(fdat[f_idx, trial])
                    end
                end
            else
                # Average directly
                for trial in 1:n_trials
                    for (fi, f_idx) in enumerate(freq_indices)
                        tf_avg[fi, ti] += abs2(fdat[f_idx, trial])
                    end
                end
                # Normalize (can also be done at end, but here keeps numbers smaller?)
                # Actually, better to normalize at very end? No, `tf_avg` is Float64.
                # Accumulate then multiply.
            end
        end
        
        if !keeptrials
            tf_avg .*= inv_n_trials
            return tf_avg, times[tois_idx], frequencies
        else
            return tf_trials, times[tois_idx], frequencies
        end
    end

    # Variable window length (cycles specified) - must loop over frequencies
    # Also handle keeptrials here
    if keeptrials
        tf_trials = zeros(n_frex, n_timepoints, n_trials)
    else
        tf_avg = zeros(n_frex, n_timepoints)
    end
    
    @inbounds for (fi, freq) in enumerate(frequencies)
        timewinidx = round(Int, cycles / freq * sample_rate)
        timewinidx += iseven(timewinidx)
        half_win = timewinidx ÷ 2
        
        hann_win = 0.5 * (1 .- cos.(2π * (0:(timewinidx-1)) / (timewinidx - 1)))
        frex_idx = round(Int, freq * timewinidx / sample_rate) + 1
        
        # Pre-allocate buffers for this frequency
        tmpdat = zeros(ComplexF64, timewinidx, n_trials)
        fdat = zeros(ComplexF64, timewinidx, n_trials)
        p_fft = plan_fft(tmpdat, 1)

        for (ti, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win

            if start_idx >= 1 && end_idx <= n_samples
                # Fast path: no boundary issues
                for trial in 1:n_trials
                    @simd for s in 1:timewinidx
                        tmpdat[s, trial] = signal[start_idx + s - 1, trial]
                    end
                end
            else
                # Slow path: handle boundaries
                fill!(tmpdat, 0.0)
                v_start = max(1, start_idx)
                v_end = min(n_samples, end_idx)
                offset = v_start - start_idx
                len = v_end - v_start + 1
                for trial in 1:n_trials
                    @simd for s in 1:len
                        tmpdat[offset + s, trial] = signal[v_start + s - 1, trial]
                    end
                end
            end

            # Apply Hanning window
            for trial in 1:n_trials
                @simd for s in 1:timewinidx
                    tmpdat[s, trial] *= hann_win[s]
                end
            end
            
            # FFT using in-place operation
            mul!(fdat, p_fft, tmpdat)
            
            # Extract power
            if keeptrials
                @simd for trial in 1:n_trials
                    tf_trials[fi, ti, trial] = abs2(fdat[frex_idx, trial])
                end
            else
                sum_p = 0.0
                @simd for trial in 1:n_trials
                    sum_p += abs2(fdat[frex_idx, trial])
                end
                tf_avg[fi, ti] = sum_p / n_trials
            end
        end
    end

    if !keeptrials
        return tf_avg, times[tois_idx], frequencies
    else
        return tf_trials, times[tois_idx], frequencies
    end
end

function tf_mtm(signal, times, sample_rate, frequencies, time_steps; t_ftimwin::Float64=0.5, tapsmofrq::Float64=4.0, keeptrials::Bool=false)
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)
    
    tois_idx = [findfirst(≈(t, atol=(1000/sample_rate)/1000), times) for t in time_steps]
    
    n_freqs = length(frequencies)
    n_timepoints = length(time_steps)
    
    timewinidx = round(Int, t_ftimwin * sample_rate)
    timewinidx += iseven(timewinidx)
    half_win = timewinidx ÷ 2
    
    # DPSS tapers
    nw = tapsmofrq * timewinidx / sample_rate / 2
    n_tapers = max(1, floor(Int, 2 * nw - 1))
    tapers = dpss_tapers(timewinidx, nw, n_tapers)
    
    # Frequency indices
    freq_indices = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
    
    tf_trials = zeros(n_freqs, n_timepoints, n_trials)
    
    # Pre-allocate all buffers outside the time loop to avoid repeated allocations
    tmpdat = zeros(ComplexF64, timewinidx, n_trials)
    tapered_buf = zeros(ComplexF64, timewinidx, n_trials)
    fdat = zeros(ComplexF64, timewinidx, n_trials)
    taper_power = zeros(n_freqs, n_trials)
    
    # Create FFT plan for in-place operations
    p_fft = plan_fft(tapered_buf, 1)

    @inbounds for (ti, center_idx) in enumerate(tois_idx)
        start_idx = center_idx - half_win
        end_idx = center_idx + half_win
        
        if start_idx < 1 || end_idx > n_samples
            continue
        end
        
        # Extract window for all trials
        for trial in 1:n_trials
            @simd for s in 1:timewinidx
                tmpdat[s, trial] = signal[start_idx + s - 1, trial]
            end
        end
        
        # Reset taper_power for this time point
        fill!(taper_power, 0.0)
        
        # Multi-taper power accumulation
        for t in 1:n_tapers
            # Apply taper - manual loop for better performance
            for trial in 1:n_trials
                @simd for s in 1:timewinidx
                    tapered_buf[s, trial] = tmpdat[s, trial] * tapers[s, t]
                end
            end
            
            # FFT all trials at once using in-place operation
            mul!(fdat, p_fft, tapered_buf)
            
            # Accumulate power at frequency indices
            for trial in 1:n_trials
                @simd for fi in 1:n_freqs
                    taper_power[fi, trial] += abs2(fdat[freq_indices[fi], trial])
                end
            end
        end
        
        # Store averaged power
        inv_n_tapers = 1.0 / n_tapers
        for trial in 1:n_trials
            @simd for fi in 1:n_freqs
                tf_trials[fi, ti, trial] = taper_power[fi, trial] * inv_n_tapers
            end
        end
    end
    
    return tf_trials, collect(Float64, time_steps), collect(Float64, frequencies)
end
function power_to_tf_dataframe(power::AbstractArray{<:Real,3}, 
                               times::AbstractVector, 
                               freqs::AbstractVector, 
                               channel_labels::Vector{Symbol})
    n_freqs, n_times, n_channels = size(power)
    n_rows = n_freqs * n_times
    
    # Pre-allocate all columns at once
    cols = Dict{Symbol, Vector{Float64}}()
    cols[:time] = Vector{Float64}(undef, n_rows)
    cols[:freq] = Vector{Float64}(undef, n_rows)
    for ch in channel_labels
        cols[ch] = Vector{Float64}(undef, n_rows)
    end
    
    # Fill columns in a single pass
    idx = 1
    for ti in 1:n_times
        for fi in 1:n_freqs
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
    return DataFrame(cols, copycols=false)[!, column_order]
end

"""
    tf_analysis(epochs::EpochData, frequencies, time_steps; 
                 method=:wavelet, keeptrials=false, kwargs...) -> TimeFreqData or TimeFreqEpochData

Perform time-frequency analysis on epoched EEG data.

This is the main time-frequency analysis function for end users.

# Arguments
- `epochs::EpochData`: Epoched EEG data
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (seconds) for TF decomposition

# Keyword Arguments
- `method::Symbol=:wavelet`: Analysis method (`:wavelet`, `:superlet`, `:multitaper`)
- `keeptrials::Bool=false`: Keep individual trials (true) or average (false)
- `channel_selection::Function=channels()`: Channel selection predicate
- Additional kwargs passed to the underlying TF method

# Returns
- `TimeFreqData`: Averaged TF result (if keeptrials=false)
- `TimeFreqEpochData`: Trial-level TF result (if keeptrials=true)

# Examples
```julia
# Wavelet with default settings
tf_data = tf_analysis(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, width=7)

# Superlet with higher order
tf_data = tf_analysis(epochs, 2:40, -0.3:0.02:0.8; method=:superlet, order=10)

# Keep individual trials
tf_epochs = tf_analysis(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, keeptrials=true)
```
"""
function tf_analysis(epochs::EpochData, frequencies, time_steps;
                      method::Symbol=:wavelet, keeptrials::Bool=false,
                      channel_selection::Function=channels(), kwargs...)
    
    # Get channels to process
    selected_channels = get_selected_channels(epochs, channel_selection; include_meta=false, include_extra=false)
    isempty(selected_channels) && error("No channels selected")
    
    n_channels = length(selected_channels)
    n_trials = length(epochs.data)
    n_samples = nrow(epochs.data[1])
    times = collect(epochs.data[1].time)
    sr = epochs.sample_rate
    
    foi = collect(Float64, frequencies)
    toi = collect(Float64, time_steps)
    
    # Pre-compute common parameters
    log_freqs = get(kwargs, :log_freqs, true)
    
    # We need to know the output dimensions (freqs and times)
    # The first channel can serve as a template
    # Determine output dimensions
    n_freqs = length(foi)
    n_times = length(time_steps)
    freqs_out = log_freqs ? exp.(range(log(float(foi[1])), log(float(foi[end])), length = n_freqs)) : collect(Float64, foi)
    times_out = collect(Float64, time_steps)
    tois_idx = [argmin(abs.(times .- t)) for t in time_steps]

    # Pre-compute wavelet parameters (only needed for wavelet method)
    # Pre-compute wavelet parameters (only needed for wavelet method)
    if method == :wavelet
        width_val = get(kwargs, :width, 7)
        # Re-calc freqs_out to ensure consistency with helper (though logic is same)
        _, half_wins, n_conv, wavelet_ffts = _prepare_morlet_wavelets(foi, width_val, sr, n_samples; log_freqs=log_freqs)
        
        # Convert epochs.data to 3D array (only needed for wavelet)
        data_3d = zeros(Float64, n_samples, n_channels, n_trials)
        for (ti, trial_df) in enumerate(epochs.data)
            for (ci, ch) in enumerate(selected_channels)
                @inbounds data_3d[:, ci, ti] .= trial_df[!, ch]
            end
        end

        ws = TFMorletWorkspace(n_conv, n_trials)
    elseif method == :superlet
        # Superlet workspace setup
        # Superlet needs convolution buffer similar to wavelet
        # Max sigma calculation
        order = get(kwargs, :order, 5)
        width = get(kwargs, :width, 3)
        orders = order isa Int ? fill(order, n_freqs) : order
        max_sigma = (width * maximum(orders)) / (2π * minimum(freqs_out))
        n_convolution = n_samples + ceil(Int, 6 * max_sigma * sr)
        n_conv_pow2 = nextpow(2, n_convolution)
        
        ws = TFMorletWorkspace(n_conv_pow2, n_trials)
    end

    if !keeptrials
        full_power = zeros(Float64, n_freqs, n_times, n_channels)
        
        for ci in 1:n_channels
            # Extract signal for this channel (samples × trials)
            signal = zeros(n_samples, n_trials)
            for (ti, trial_df) in enumerate(epochs.data)
                signal[:, ti] = trial_df[!, selected_channels[ci]]
            end
            
            if method == :wavelet
                tf_morlet_optimized!(@view(full_power[:, :, ci]), ws, data_3d, ci, 
                                     n_samples, n_trials, tois_idx, half_wins, 
                                     wavelet_ffts, false)
            elseif method == :multitaper
                # Filter kwargs for multitaper
                multitaper_kwargs = Base.filter(p -> p.first in (:taper, :t_ftimwin, :tapsmofrq), kwargs)
                power, _, _ = tf_multitaper(signal, times, sr, foi, toi; keeptrials=false, multitaper_kwargs...)
                full_power[:, :, ci] = power
            elseif method == :superlet
                # Filter kwargs for superlet
                superlet_kwargs = Base.filter(p -> p.first in (:order, :width, :combine), kwargs)
                # Parse defaults
                sl_order = get(kwargs, :order, 5)
                sl_width = get(kwargs, :width, 3)
                sl_combine = get(kwargs, :combine, :multiplicative)
                orders = sl_order isa Int ? fill(sl_order, n_freqs) : sl_order
                
                tf_superlet_optimized!(@view(full_power[:, :, ci]), ws, signal, ci,
                                     n_samples, n_trials, tois_idx, freqs_out, orders, sl_width, sl_combine, false, sr)
            else
                error("Unknown method: $method. Use :wavelet, :multitaper, or :superlet")
            end
        end
        
        df = power_to_tf_dataframe(full_power, times_out, freqs_out, selected_channels)
        
        return TimeFreqData(epochs.file, epochs.condition, epochs.condition_name, df, 
                            epochs.layout, epochs.sample_rate, method, nothing, epochs.analysis_info)
    else
        # Keep trials version
        channel_results = Vector{Array{Float64, 3}}(undef, n_channels)
        
        for ci in 1:n_channels
            # Extract signal for this channel (samples × trials)
            signal = zeros(n_samples, n_trials)
            for (ti, trial_df) in enumerate(epochs.data)
                signal[:, ti] = trial_df[!, selected_channels[ci]]
            end
            
            pwr = zeros(Float64, n_freqs, n_times, n_trials)
            if method == :wavelet
                tf_morlet_optimized!(pwr, ws, data_3d, ci, 
                                     n_samples, n_trials, tois_idx, half_wins, 
                                     wavelet_ffts, true)
            elseif method == :multitaper
                # Filter kwargs for multitaper
                multitaper_kwargs = Base.filter(p -> p.first in (:taper, :t_ftimwin, :tapsmofrq), kwargs)
                power, _, _ = tf_multitaper(signal, times, sr, foi, toi; keeptrials=true, multitaper_kwargs...)
                # Ensure 3D (freqs × times × trials)
                if ndims(power) == 2
                    pwr = reshape(power, size(power)..., 1)
                else
                    pwr = power
                end
            elseif method == :superlet
                # Filter kwargs for superlet
                superlet_kwargs = Base.filter(p -> p.first in (:order, :width, :combine), kwargs)
                # Parse defaults
                sl_order = get(kwargs, :order, 5)
                sl_width = get(kwargs, :width, 3)
                sl_combine = get(kwargs, :combine, :multiplicative)
                orders = sl_order isa Int ? fill(sl_order, n_freqs) : sl_order
                
                tf_superlet_optimized!(pwr, ws, signal, ci,
                                     n_samples, n_trials, tois_idx, freqs_out, orders, sl_width, sl_combine, true, sr)
                # Ensure 3D (freqs × times × trials)
                if ndims(power) == 2
                    pwr = reshape(power, size(power)..., 1)
                else
                    pwr = power
                end
            else
                error("Unknown method: $method. Use :wavelet, :multitaper, or :superlet")
            end
            channel_results[ci] = pwr
        end
        
        # Reorganize into Vector{DataFrame} (one per trial)
        trial_dfs = Vector{DataFrame}(undef, n_trials)
        for trial_idx in 1:n_trials
            # Temporary matrix for this trial (freqs × times × channels)
            trial_power = zeros(n_freqs, n_times, n_channels)
            for ci in 1:n_channels
                @views trial_power[:, :, ci] .= channel_results[ci][:, :, trial_idx]
            end
            trial_dfs[trial_idx] = power_to_tf_dataframe(trial_power, times_out, freqs_out, selected_channels)
        end
        
        return TimeFreqEpochData(
            epochs.file,
            epochs.condition,
            epochs.condition_name,
            trial_dfs,
            epochs.layout,
            epochs.sample_rate,
            method,
            nothing,
            epochs.analysis_info
        )
    end
end

# ============================================================================
# Batch TF transform: Process multiple files

"""Generate default output directory name for TF transform operation."""
function _default_tf_output_dir(input_dir::String, pattern::String, method::Symbol)
    joinpath(input_dir, "tf_$(pattern)_$(method)")
end

"""
Process a single file through TF transform pipeline.
Returns BatchResult with success/failure info.
"""
function _process_tf_file(
    filepath::String,
    output_path::String,
    frequencies,
    time_steps,
    method::Symbol,
    keeptrials::Bool,
    condition_selection::Function,
    tf_kwargs::NamedTuple
)
    filename = basename(filepath)

    # Load data
    data = load_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is EpochData
    if !(data isa Vector{EpochData})
        return BatchResult(false, filename, "Invalid data type: expected Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Apply TF transform to each condition
    tf_results = Vector{Union{TimeFreqData,TimeFreqEpochData}}()
    for epochs in data
        tf_result = tf_analysis(epochs, frequencies, time_steps; 
                                  method=method, keeptrials=keeptrials, tf_kwargs...)
        push!(tf_results, tf_result)
    end

    # Save results
    jldsave(output_path; data=tf_results)

    return BatchResult(true, filename, "TF transform completed ($(length(tf_results)) conditions)")
end

"""
    tf_analysis(file_pattern::String, frequencies, time_steps;
                 method=:wavelet, keeptrials=false,
                 input_dir=pwd(), output_dir=nothing,
                 participant_selection=participants(), 
                 condition_selection=conditions(),
                 kwargs...) -> Nothing

Apply time-frequency analysis to EEG data from JLD2 files and save results.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "cleaned", or custom pattern)
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (seconds) for TF decomposition

# Keyword Arguments
- `method::Symbol=:wavelet`: Analysis method (`:wavelet`, `:superlet`, `:multitaper`)
- `keeptrials::Bool=false`: Keep individual trials (true) or average (false)
- `input_dir::String=pwd()`: Input directory containing JLD2 files
- `output_dir::Union{String,Nothing}=nothing`: Output directory (default: creates subdirectory)
- `participant_selection::Function=participants()`: Participant selection predicate
- `condition_selection::Function=conditions()`: Condition selection predicate
- Additional kwargs passed to underlying TF method (e.g., `width`, `order`, `taper`)

# Example
```julia
# Wavelet TF on all epoch files
tf_analysis("epochs", 2:40, -0.3:0.02:0.8; method=:wavelet, width=7)

# Superlet with specific participants
tf_analysis("epochs", 2:40, -0.3:0.02:0.8; 
             method=:superlet, order=5,
             participant_selection=participants([1, 2, 3]))

# Keep individual trials
tf_analysis("epochs", 2:40, -0.3:0.02:0.8; 
             method=:wavelet, keeptrials=true)
```
"""
function tf_analysis(
    file_pattern::String,
    frequencies,
    time_steps;
    method::Symbol=:wavelet,
    keeptrials::Bool=false,
    input_dir::String=pwd(),
    output_dir::Union{String,Nothing}=nothing,
    participant_selection::Function=participants(),
    condition_selection::Function=conditions(),
    kwargs...
)
    # Setup logging
    log_file = "tf_analysis.log"
    setup_global_logging(log_file)

    try
        @info ""
        @info "Batch TF transform started at $(now())"
        @log_call "tf_analysis"

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_tf_output_dir(input_dir, file_pattern, method))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "TF method: $method, Frequencies: $(first(frequencies))-$(last(frequencies)) Hz"
        @info "Time steps: $(first(time_steps))-$(last(time_steps)) s, keeptrials: $keeptrials"

        # Convert kwargs to NamedTuple for passing through
        tf_kwargs = NamedTuple(kwargs)

        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> _process_tf_file(
            input_path, output_path, frequencies, time_steps, 
            method, keeptrials, condition_selection, tf_kwargs
        )

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; 
                                       operation_name="TF transform")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

