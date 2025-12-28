"""
Time-Frequency Analysis Module

Provides FieldTrip-style time-frequency analysis with multiple methods:
- `tf_wavelet`: Morlet wavelet convolution
- `tf_multitaper`: Multi-taper sliding window analysis (Hanning or DPSS)
- `tf_spectrum`: Multi-taper FFT (non-time-resolved power spectrum)

Or use the unified dispatcher:
- `tf_analysis(method, ...)`: Dispatches to the appropriate method

# Examples
```julia
# Direct function calls (Option 1)
power, t, f = tf_wavelet(signal, times, sr, 2:30, -0.5:0.05:1.5; width=7)
power, t, f = tf_multitaper(signal, times, sr, 2:30, -0.5:0.05:1.5; taper=:dpss)
power, f = tf_spectrum(signal, sr, 2:40; taper=:dpss)

# Unified dispatcher (Option 2)
power, t, f = tf_analysis(:wavelet, signal, times, sr, 2:30, -0.5:0.05:1.5; width=7)
power, t, f = tf_analysis(:multitaper, signal, times, sr, 2:30, -0.5:0.05:1.5)
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
    freqs = collect(Float64, frequencies)
    toi = collect(Float64, time_steps)
    
    # Convert width to vector format expected by tf_morlet
    cycles = width isa Int ? [width] : collect(width)
    
    power, times_out, freqs_out = tf_morlet(signal, times, sample_rate, freqs, cycles; 
                                             tois=toi, log_freqs=log_freqs)
    
    if !keeptrials && ndims(power) == 3
        power = dropdims(mean(power, dims=3), dims=3)
    end
    
    return power, times_out, freqs_out
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
    n_times = length(toi)
    
    # Expand order to per-frequency if scalar
    orders = order isa Int ? fill(order, n_freqs) : order
    if length(orders) != n_freqs
        throw(ArgumentError("order vector must have same length as frequencies"))
    end
    
    # Pre-allocate output
    power = zeros(n_freqs, n_times, n_trials)
    
    # Process each frequency
    for (fi, freq) in enumerate(freqs)
        ord = orders[fi]
        
        # Generate cycle widths based on combine mode
        if combine == :multiplicative
            cycles = [width * i for i in 1:ord]
        elseif combine == :additive
            cycles = [width + i - 1 for i in 1:ord]
        else
            throw(ArgumentError("combine must be :multiplicative or :additive"))
        end
        
        # Compute wavelet response for each cycle width
        wavelet_powers = zeros(ord, n_times, n_trials)
        for (ci, cyc) in enumerate(cycles)
            # Use tf_morlet for each single frequency with specific cycle width
            pwr, _, _ = tf_morlet(signal, times, sample_rate, [freq], [cyc]; 
                                   tois=toi, log_freqs=false)
            wavelet_powers[ci, :, :] = pwr[1, :, :]
        end
        
        # Geometric mean across wavelets (per time point and trial)
        for ti in 1:n_times
            for trial in 1:n_trials
                # Product of powers ^ (1/order) = geometric mean
                power[fi, ti, trial] = prod(wavelet_powers[:, ti, trial]) ^ (1.0 / ord)
            end
        end
    end
    
    # Average over trials if requested
    if !keeptrials && n_trials > 1
        power = dropdims(mean(power, dims=3), dims=3)
    elseif n_trials == 1
        power = dropdims(power, dims=3)
    end
    
    return power, toi, freqs
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
                                                  window_length=t_ftimwin)
    elseif taper == :dpss
        power, times_out, freqs_out = tf_mtm(signal, times, sample_rate, freqs, toi;
                                              t_ftimwin=t_ftimwin, tapsmofrq=tapsmofrq)
    else
        throw(ArgumentError("Unknown taper: $taper. Use :hanning or :dpss"))
    end
    
    if !keeptrials && ndims(power) == 3
        power = dropdims(mean(power, dims=3), dims=3)
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
# Option 2: Unified tf_analysis dispatcher
# ============================================================================

"""
    tf_analysis(method, signal, times, sample_rate, frequencies, time_steps; kwargs...)

Unified time-frequency analysis dispatcher.

# Arguments
- `method::Symbol`: Analysis method (`:wavelet`, `:multitaper`, or `:spectrum`)
- `signal`: Signal matrix (samples × trials) or vector
- `times`: Time vector (not used for `:spectrum`)
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest (Hz)
- `time_steps`: Times of interest (not used for `:spectrum`)

# Keyword Arguments
Passed to the underlying method. See `tf_wavelet`, `tf_multitaper`, `tf_spectrum`.

# Returns
- For `:wavelet`, `:multitaper`: (power, times_out, freqs_out)
- For `:spectrum`: (power, freqs_out)

# Example
```julia
# Wavelet analysis
power, t, f = tf_analysis(:wavelet, signal, times, 256, 2:30, -0.5:0.05:1.5; width=7)

# Multi-taper with DPSS
power, t, f = tf_analysis(:multitaper, signal, times, 256, 2:30, -0.5:0.05:1.5; taper=:dpss)

# Power spectrum (times/time_steps ignored)
power, f = tf_analysis(:spectrum, signal, nothing, 256, 2:40, nothing; taper=:dpss)
```
"""
function tf_analysis(method::Symbol, signal, times, sample_rate, frequencies, time_steps; kwargs...)
    if method == :wavelet
        return tf_wavelet(signal, times, sample_rate, frequencies, time_steps; kwargs...)
    elseif method == :superlet
        return tf_superlet(signal, times, sample_rate, frequencies, time_steps; kwargs...)
    elseif method == :multitaper
        return tf_multitaper(signal, times, sample_rate, frequencies, time_steps; kwargs...)
    elseif method == :spectrum
        return tf_spectrum(signal, sample_rate, frequencies; kwargs...)
    else
        throw(ArgumentError("Unknown method: $method. Use :wavelet, :superlet, :multitaper, or :spectrum"))
    end
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
# Multi-taper sliding window (mtmconvol with DPSS)
# ============================================================================

"""
    tf_mtm(signal, times, sample_rate, frequencies, time_steps; t_ftimwin=0.5, tapsmofrq=4.0)

Multi-taper time-frequency analysis with DPSS tapers.

# Arguments
- `signal`: Signal matrix (samples × trials) or vector
- `times`: Time vector
- `sample_rate`: Sampling rate in Hz
- `frequencies`: Frequencies of interest
- `time_steps`: Times of interest
- `t_ftimwin`: Window length in seconds
- `tapsmofrq`: Frequency smoothing in Hz

# Returns
- `power`: Time-frequency power (freqs × times × trials)
- `times_out`: Output time vector
- `frequencies`: Frequency vector
"""
function tf_mtm(signal, times, sample_rate, frequencies, time_steps; t_ftimwin::Float64=0.5, tapsmofrq::Float64=4.0)
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)
    
    # Convert time_steps to indices
    tois_idx = [findfirst(≈(t, atol=(1000/sample_rate)/1000), times) for t in time_steps]
    
    n_freqs = length(frequencies)
    n_timepoints = length(time_steps)
    
    # Window parameters
    timewinidx = round(Int, t_ftimwin * sample_rate)
    timewinidx += iseven(timewinidx)
    half_win = fld(timewinidx, 2)
    
    # DPSS tapers
    nw = tapsmofrq * timewinidx / sample_rate / 2
    n_tapers = max(1, floor(Int, 2 * nw - 1))
    tapers = dpss_tapers(timewinidx, nw, n_tapers)
    
    # Frequency indices for this window size
    freq_indices = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
    
    tf_trials = zeros(n_freqs, n_timepoints, n_trials)
    tmpdat = zeros(timewinidx)
    
    @inbounds for trial in 1:n_trials
        for (ti, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win
            
            # Handle edge cases
            if start_idx < 1 || end_idx > n_samples
                continue
            end
            
            # Extract window
            @views tmpdat .= signal[start_idx:end_idx, trial]
            
            # Multi-taper power
            taper_power = zeros(n_freqs)
            for t in 1:n_tapers
                tapered = tmpdat .* @view(tapers[:, t])
                spectrum = fft(tapered)
                taper_power .+= abs2.(@view(spectrum[freq_indices]))
            end
            tf_trials[:, ti, trial] = taper_power ./ n_tapers
        end
    end
    
    return tf_trials, collect(Float64, time_steps), collect(Float64, frequencies)
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


function tf_morlet(signal, times, sample_rate, freqs, cycles; log_freqs = true, tois = nothing)

    # if vector
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal

    # data dimensions
    n_samples, n_trials = size(signal)
    n_freq = length(freqs)

    # time of interest
    tois_idx = 1:n_samples
    if !isnothing(tois)
        tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in tois]
    end

    # frequency steps
    if log_freqs
        freqs = exp.(range(log(freqs[1]), log(freqs[end]), length = n_freq))
    end

    # variable or fixed cycles/frequency
    if length(cycles) == 1
        cycles = cycles[1] ./ (2 * pi .* freqs) # variable cycles/frequency
    else
        cycles = exp.(range(log(cycles[1]), log(cycles[end]), length = n_freq)) ./ (2 * pi .* freqs)
    end

    # pre-compute signal FFT once maintining trial dimension
    # use lowest freq for max padding
    n_convolution = n_samples + ceil(Int, 6 * (maximum(cycles) / (2π * freqs[1])) * sample_rate)
    n_conv_pow2 = nextpow(2, n_convolution)
    signal_fft = fft([signal; zeros(n_conv_pow2 - n_samples, n_trials)], 1)  # FFT along time dimension

    tf_trials = zeros(n_freq, length(tois_idx), n_trials)
    @inbounds for idx_freq in eachindex(freqs)
        # frequency-dependent window length
        window_length = min(n_convolution, ceil(Int, 6 * cycles[idx_freq] * sample_rate))  # Cap window length
        window_length += iseven(window_length)  # ensure odd length

        #  wavelet
        wavelet_time = range(-window_length / 2, window_length / 2, length = window_length) ./ sample_rate
        wavelet =
            sqrt(1 ./ (cycles[idx_freq] .* sqrt(pi))) .* exp.(2im * pi * freqs[idx_freq] .* wavelet_time) .*
            exp.(-wavelet_time .^ 2 ./ (2 * (cycles[idx_freq] .^ 2)))

        # zero pad wavelet to match signal FFT length
        wavelet_fft = fft([wavelet; zeros(n_conv_pow2 - length(wavelet))])

        # convolution across trials
        eegconv = ifft(wavelet_fft .* signal_fft, 1)  # IFFT along time dimension

        # trim to original size and handle edge effects
        half_window = window_length ÷ 2
        valid_idx = (half_window+1):(half_window+n_samples)

        # store result
        @views tf_trials[idx_freq, :, :] = abs2.(eegconv[valid_idx, :])[tois_idx, :]

    end

    return tf_trials, times[tois_idx], freqs

end


function tf_hanning(signal, times, sample_rate, frequencies, time_steps; window_length = nothing, cycles = nothing)

    if isnothing(window_length) && isnothing(cycles)
        throw(ArgumentError("Must specify either window_length or cycles"))
    end

    # if vector
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal

    # convert time_steps to index values
    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_steps]

    # Pre-allocate output and temporary arrays
    n_frex = length(frequencies)
    n_timepoints = length(time_steps)
    n_trials = size(signal, 2)
    tf_trials = zeros(n_frex, n_timepoints, n_trials)

    # Pre-compute maximum window length
    max_timewinidx = if !isnothing(cycles)
        round(Int, cycles / minimum(frequencies) * sample_rate)
    else
        round(Int, window_length * sample_rate)
    end
    max_timewinidx += iseven(max_timewinidx)

    # Pre-allocate reusable arrays
    tmpdat = zeros(max_timewinidx, n_trials)
    fdat = zeros(ComplexF64, max_timewinidx, n_trials)

    # Pre-compute window lengths and Hanning windows for each frequency
    timewin_indices = Vector{Int}(undef, n_frex)
    hann_windows = Vector{Vector{Float64}}(undef, n_frex)
    frex_indices = Vector{Int}(undef, n_frex)

    for (fi, freq) in enumerate(frequencies)
        timewinidx = if !isnothing(cycles)
            round(Int, cycles / freq * sample_rate)
        else
            round(Int, window_length * sample_rate)
        end
        timewinidx += iseven(timewinidx)
        timewin_indices[fi] = timewinidx

        # Pre-compute Hanning window
        hann_windows[fi] = 0.5 * (1 .- cos.(2π * (0:(timewinidx-1)) / (timewinidx - 1)))
        frex_indices[fi] = round(Int, freq * timewinidx / sample_rate) + 1
    end

    # Main computation loop
    @inbounds for fi = 1:n_frex
        timewinidx = timewin_indices[fi]
        hann_win = hann_windows[fi]
        frex_idx = frex_indices[fi]
        half_win = fld(timewinidx, 2)

        for (timepointi, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win

            # Handle edge cases and data extraction
            fill!(view(tmpdat, 1:timewinidx, :), 0)
            if start_idx < 1 || end_idx > size(signal, 1)
                pad_left = max(1 - start_idx, 0)
                valid_start = max(start_idx, 1)
                valid_end = min(end_idx, size(signal, 1))
                valid_length = valid_end - valid_start + 1
                @views tmpdat[(pad_left+1):(pad_left+valid_length), :] .= signal[valid_start:valid_end, :]
            else
                @views tmpdat[1:timewinidx, :] .= signal[start_idx:end_idx, :]
            end

            # Apply Hanning taper and compute FFT
            @views tmpdat[1:timewinidx, :] .*= hann_win
            @views fdat[1:timewinidx, :] .= tmpdat[1:timewinidx, :]
            fft!(view(fdat, 1:timewinidx, :), 1)  # FFT along time dimension

            # Store power
            @views tf_trials[fi, timepointi, :] .= abs2.(fdat[frex_idx, :])
        end
    end

    return tf_trials, times[tois_idx], frequencies

end

# ============================================================================
# High-level TF transform: EpochData → TimeFreqData / TimeFreqEpochData
# ============================================================================

"""
    power_to_tf_dataframe(power, times, freqs, channel_labels) -> DataFrame

Convert a 3D power array (freqs × times × channels) to a long-format DataFrame.

# Arguments
- `power`: 3D array (freqs × times × channels)
- `times`: Time vector
- `freqs`: Frequency vector
- `channel_labels`: Vector of channel symbols

# Returns
- `DataFrame`: With columns [:time, :freq, channel1, channel2, ...]
"""
function power_to_tf_dataframe(power::AbstractArray{<:Real,3}, 
                               times::AbstractVector, 
                               freqs::AbstractVector, 
                               channel_labels::Vector{Symbol})
    n_freqs, n_times, n_channels = size(power)
    n_rows = n_freqs * n_times
    
    # Pre-allocate columns
    time_col = Vector{Float64}(undef, n_rows)
    freq_col = Vector{Float64}(undef, n_rows)
    
    # Build time × freq grid (all freqs for each time point)
    idx = 1
    for ti in 1:n_times
        for fi in 1:n_freqs
            time_col[idx] = times[ti]
            freq_col[idx] = freqs[fi]
            idx += 1
        end
    end
    
    # Start DataFrame with time and freq columns
    df = DataFrame(:time => time_col, :freq => freq_col)
    
    # Add each channel's power values
    for (ci, ch) in enumerate(channel_labels)
        ch_power = Vector{Float64}(undef, n_rows)
        idx = 1
        for ti in 1:n_times
            for fi in 1:n_freqs
                ch_power[idx] = power[fi, ti, ci]
                idx += 1
            end
        end
        df[!, ch] = ch_power
    end
    
    return df
end

"""
    tf_transform(epochs::EpochData, frequencies, time_steps; 
                 method=:wavelet, keeptrials=false, kwargs...) -> TimeFreqData or TimeFreqEpochData

Perform time-frequency analysis on epoched EEG data.

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
tf_data = tf_transform(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, width=7)

# Superlet with higher order
tf_data = tf_transform(epochs, 2:40, -0.3:0.02:0.8; method=:superlet, order=10)

# Keep individual trials
tf_epochs = tf_transform(epochs, 2:40, -0.3:0.02:0.8; method=:wavelet, keeptrials=true)
```
"""
function tf_transform(epochs::EpochData, frequencies, time_steps;
                      method::Symbol=:wavelet, keeptrials::Bool=false,
                      channel_selection::Function=channels(), kwargs...)
    
    # Get channels to process
    ch_labels = channel_labels(epochs)
    selected_channels = filter(ch -> channel_selection(ch), ch_labels)
    isempty(selected_channels) && error("No channels selected")
    
    n_channels = length(selected_channels)
    n_trials = length(epochs.data)
    n_samples = nrow(epochs.data[1])
    times = collect(epochs.data[1].time)
    sr = epochs.sample_rate
    
    foi = collect(Float64, frequencies)
    toi = collect(Float64, time_steps)
    
    # First pass: determine output dimensions
    # Process first channel to get dimensions
    test_signal = zeros(n_samples, n_trials)
    for (i, trial_df) in enumerate(epochs.data)
        test_signal[:, i] = trial_df[!, selected_channels[1]]
    end
    
    if method == :wavelet
        test_power, times_out, freqs_out = tf_morlet(test_signal, times, sr, foi, 
            get(kwargs, :width, 7) isa Int ? [get(kwargs, :width, 7)] : collect(get(kwargs, :width, 7)); 
            tois=toi, log_freqs=get(kwargs, :log_freqs, true))
    elseif method == :superlet
        test_power, times_out, freqs_out = tf_superlet(test_signal, times, sr, foi, toi; 
            keeptrials=true, filter(p -> p.first in (:order, :width, :combine), kwargs)...)
        test_power = ndims(test_power) == 2 ? reshape(test_power, size(test_power)..., 1) : test_power
    elseif method == :multitaper
        test_power, times_out, freqs_out = tf_multitaper(test_signal, times, sr, foi, toi; 
            keeptrials=true, filter(p -> p.first in (:taper, :t_ftimwin, :tapsmofrq), kwargs)...)
        test_power = ndims(test_power) == 2 ? reshape(test_power, size(test_power)..., 1) : test_power
    else
        error("Unknown method: $method. Use :wavelet, :superlet, or :multitaper")
    end
    
    n_freqs = length(freqs_out)
    n_times = length(times_out)
    
    if keeptrials
        # Process each trial separately, store in Vector{DataFrame}
        trial_dfs = Vector{DataFrame}(undef, n_trials)
        
        for trial_idx in 1:n_trials
            # Collect power for all channels in this trial
            trial_power = zeros(n_freqs, n_times, n_channels)
            
            for (ci, ch) in enumerate(selected_channels)
                signal = reshape(epochs.data[trial_idx][!, ch], :, 1)
                
                if method == :wavelet
                    width_arg = get(kwargs, :width, 7) isa Int ? [get(kwargs, :width, 7)] : collect(get(kwargs, :width, 7))
                    power, _, _ = tf_morlet(signal, times, sr, foi, width_arg; 
                                            tois=toi, log_freqs=get(kwargs, :log_freqs, true))
                elseif method == :superlet
                    power, _, _ = tf_superlet(signal, times, sr, foi, toi; keeptrials=false,
                        filter(p -> p.first in (:order, :width, :combine), kwargs)...)
                elseif method == :multitaper
                    power, _, _ = tf_multitaper(signal, times, sr, foi, toi; keeptrials=false,
                        filter(p -> p.first in (:taper, :t_ftimwin, :tapsmofrq), kwargs)...)
                end
                
                trial_power[:, :, ci] = ndims(power) == 3 ? dropdims(power, dims=3) : power
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
            epochs.analysis_info
        )
    else
        # Average across trials: collect power for all channels
        avg_power = zeros(n_freqs, n_times, n_channels)
        
        for (ci, ch) in enumerate(selected_channels)
            # Extract signal for this channel (samples × trials)
            signal = zeros(n_samples, n_trials)
            for (i, trial_df) in enumerate(epochs.data)
                signal[:, i] = trial_df[!, ch]
            end
            
            # Compute TF and average across trials
            if method == :wavelet
                width_arg = get(kwargs, :width, 7) isa Int ? [get(kwargs, :width, 7)] : collect(get(kwargs, :width, 7))
                power, _, _ = tf_morlet(signal, times, sr, foi, width_arg; 
                                        tois=toi, log_freqs=get(kwargs, :log_freqs, true))
            elseif method == :superlet
                power, _, _ = tf_superlet(signal, times, sr, foi, toi; keeptrials=true,
                    filter(p -> p.first in (:order, :width, :combine), kwargs)...)
            elseif method == :multitaper
                power, _, _ = tf_multitaper(signal, times, sr, foi, toi; keeptrials=true,
                    filter(p -> p.first in (:taper, :t_ftimwin, :tapsmofrq), kwargs)...)
            end
            
            # Average over trials
            if ndims(power) == 3
                avg_power[:, :, ci] = dropdims(mean(power, dims=3), dims=3)
            else
                avg_power[:, :, ci] = power
            end
        end
        
        df = power_to_tf_dataframe(avg_power, times_out, freqs_out, selected_channels)
        
        return TimeFreqData(
            epochs.file,
            epochs.condition,
            epochs.condition_name,
            df,
            epochs.layout,
            epochs.sample_rate,
            method,
            epochs.analysis_info
        )
    end
end

# ============================================================================
# Batch TF transform: Process multiple files
# ============================================================================

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
        tf_result = tf_transform(epochs, frequencies, time_steps; 
                                  method=method, keeptrials=keeptrials, tf_kwargs...)
        push!(tf_results, tf_result)
    end

    # Save results
    jldsave(output_path; data=tf_results)

    return BatchResult(true, filename, "TF transform completed ($(length(tf_results)) conditions)")
end

"""
    tf_transform(file_pattern::String, frequencies, time_steps;
                 method=:wavelet, keeptrials=false,
                 input_dir=pwd(), output_dir=nothing,
                 participant_selection=participants(), 
                 condition_selection=conditions(),
                 kwargs...) -> Nothing

Apply time-frequency transform to EEG data from JLD2 files and save results.

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
tf_transform("epochs", 2:40, -0.3:0.02:0.8; method=:wavelet, width=7)

# Superlet with specific participants
tf_transform("epochs", 2:40, -0.3:0.02:0.8; 
             method=:superlet, order=5,
             participant_selection=participants([1, 2, 3]))

# Keep individual trials
tf_transform("epochs", 2:40, -0.3:0.02:0.8; 
             method=:wavelet, keeptrials=true)
```
"""
function tf_transform(
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
    log_file = "tf_transform.log"
    setup_global_logging(log_file)

    try
        @info ""
        @info "Batch TF transform started at $(now())"
        @log_call "tf_transform"

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

