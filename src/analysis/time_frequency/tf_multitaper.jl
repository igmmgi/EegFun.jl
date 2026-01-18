"""
    tf_multitaper(dat::EpochData; 
                  channel_selection::Function=channels(),
                  sample_selection::Function=samples(),
                  frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40),
                  cycles::Real,
                  frequency_smoothing::Union{Nothing,Real}=nothing,
                  time_steps::Real=0.05,
                  pad::Union{Nothing,Symbol}=nothing,
                  return_trials::Bool=false,
                  filter_edges::Bool=true)

Multitaper time-frequency analysis using DPSS (Discrete Prolate Spheroidal Sequences) tapers (Cohen Chapter 16, equivalent to FieldTrip's 'mtmconvol' method).

Uses multiple orthogonal tapers (Slepian sequences) to reduce variance in spectral estimates compared to single-taper methods. Uses adaptive window lengths (cycles per frequency) to match the time-frequency trade-off of Morlet wavelets and adaptive STFT.

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `sample_selection::Function=samples()`: Sample selection predicate. See `samples()` for options.
  - Example: `sample_selection=samples((-0.5, 2.0))` for time window from -0.5 to 2.0 seconds
  - Example: `sample_selection=samples()` for all time points (default)
  - Default: all samples
- `frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40)`: Frequency specification.
  - Can be any range or vector of frequencies in Hz
  - For linear spacing: `frequencies=1:1:40` or `frequencies=range(1, 40, length=40)`
  - For logarithmic spacing: `frequencies=logrange(1, 40, length=30)`
  - Default: `range(1, 40, length=40)` (40 linearly-spaced frequencies from 1 to 40 Hz)
- `cycles::Real`: Number of cycles per frequency. Window length = cycles / frequency (in seconds).
  - Example: `cycles=5` uses 5 cycles for all frequencies (FieldTrip: `cfg.t_ftimwin = 5./cfg.foi`)
  - Lower frequencies will have longer windows, higher frequencies will have shorter windows
- `frequency_smoothing::Union{Nothing,Real}=nothing`: Frequency smoothing parameter (FieldTrip's `tapsmofrq`).
  - If `nothing`, uses `frequency_smoothing = 0.4 * frequency` (FieldTrip default: `cfg.tapsmofrq = 0.4 * cfg.foi`)
  - If a number, uses that value multiplied by frequency: `tapsmofrq = frequency_smoothing * frequency`
  - Example: `frequency_smoothing=0.4` matches FieldTrip's default
  - Controls time-bandwidth product: `NW = tapsmofrq * window_length / 2`
  - Number of tapers used: `K = 2*NW - 1` (rounded down)
- `time_steps::Real=0.05`: Step size for extracting time points in seconds.
  - Creates time points from the selected time range with the specified step size
  - Default: `0.05` (50 ms)
  - Example: `time_steps=0.01` creates time points every 0.01 seconds within the selected time window
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default)
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.
- `filter_edges::Bool=true`: If `true` (default), filters out edge regions where the window extends beyond the data
  (FieldTrip's "cone of influence" - marks edges as NaN). If `false`, uses all convolution results including edges.

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Examples
```julia
# Default: linear frequencies 1-40 Hz, 40 points
tf_data = tf_multitaper(epochs; cycles=5)

# Log-spaced frequencies with default frequency smoothing (0.4 * frequency)
tf_data = tf_multitaper(epochs; frequencies=logrange(2, 80, length=30), cycles=5)

# Custom frequency smoothing
# FieldTrip equivalent: cfg.t_ftimwin = 5./cfg.foi, cfg.tapsmofrq = 0.4 * cfg.foi
tf_data = tf_multitaper(epochs; frequencies=1:2:30, cycles=5, frequency_smoothing=0.4)
```
"""
function tf_multitaper(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    frequencies::Union{AbstractRange,AbstractVector{<:Real}} = range(1, 40, length=40),
    time_steps::Real = 0.05,
    cycles::Real,
    frequency_smoothing::Union{Nothing,Real} = nothing,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
    filter_edges::Bool = true,
)
    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end

    # Subset data with channel and sample selection
    dat = subset(dat; channel_selection = channel_selection, sample_selection = sample_selection)
    isempty(dat.data) && error("No data remaining after subsetting")

    # Get selected channels (after subsetting)
    selected_channels = channel_labels(dat)

    # Use frequency input directly (ranges and vectors both work)
    num_frex = length(frequencies)
    
    # Validate frequencies
    if num_frex == 0
        error("`frequencies` must contain at least one frequency")
    end
    if any(f -> f <= 0, frequencies)
        error("All frequencies in `frequencies` must be positive")
    end

    # Validate cycles parameter (adaptive window only)
    if cycles <= 0
        error("`cycles` must be positive, got $cycles")
    end

    # Default frequency smoothing (FieldTrip: cfg.tapsmofrq = 0.4 * cfg.foi)
    if isnothing(frequency_smoothing)
        frequency_smoothing = 0.4
    end
    if frequency_smoothing <= 0
        error("`frequency_smoothing` must be positive, got $frequency_smoothing")
    end

    # Get original data time range (before padding) - these are the time points we want in output
    n_samples_original_unpadded = n_samples(dat)  # Store original unpadded length for edge filtering
    times_original = time(dat)

    # Apply padding if requested (mutating dat directly for consistency with tf_morlet)
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get sample rate and time vector from processed data
    times_processed = time(dat)
    n_samples_processed = n_samples(dat)  # Number of samples per epoch (may be padded)

    # Handle time_steps parameter - determine which time points to extract from results
    # After padding, processed data has extended time range - validate against processed data
    # Create time points with specified step size within the selected time range
    time_min = minimum(times_original)
    time_max = maximum(times_original)
    time_steps_range = time_min:time_steps:time_max
        time_indices, times_out = find_times(times_processed, time_steps_range)
        if isempty(time_indices)
        error("No valid time points found with step size $time_steps in range ($time_min to $time_max seconds)")
    end

    # Use processed data dimensions
    n_samples_original = n_samples_processed

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples_original

    # Use frequencies directly (convert to vector if needed for indexing)
    freqs = collect(frequencies)

    # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
    cycles_vec = fill(Float64(cycles), num_frex)
    window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
    if any(n -> n < 2, n_window_samples_per_freq)
        min_samples = minimum(n_window_samples_per_freq)
        error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end

    n_times = length(times_out)

    # Pre-compute DPSS tapers for each frequency
    tapers_per_freq = Vector{Matrix{Float64}}(undef, num_frex)  # Each element is (n_window_samples × n_tapers)
    n_tapers_per_freq = zeros(Int, num_frex)
    inv_n_tapers_per_freq = Vector{Float64}(undef, num_frex)
    
    # Determine padding length: pad to at least the data length and largest window
    max_window_samples = maximum(n_window_samples_per_freq)
    n_samples_padded = max(n_samples_per_epoch, max_window_samples)
    
    # Pre-compute tapered wavelets and their FFTs for each frequency-taper combination
    # For multitaper: tapered_wavelet = taper * (cos + i*sin) at target frequency
    tapered_wavelet_ffts = Vector{Vector{Vector{ComplexF64}}}(undef, num_frex)  # [frequency][taper] = FFT of tapered wavelet
    
    # Pre-compute FFT plans for padded data (batch process all trials)
    template_padded_batch = zeros(ComplexF64, n_samples_padded, n_trials)
    fft_plan_padded_batch = plan_fft(template_padded_batch, 1, flags = FFTW.MEASURE)
    
    # Pre-compute IFFT plan (per trial)
    template_complex = zeros(ComplexF64, n_samples_padded)
    ifft_plan_padded = plan_ifft(template_complex, flags = FFTW.MEASURE)
    
    for (fi, freq) in enumerate(freqs)
        n_window_samples = n_window_samples_per_freq[fi]
        window_length_sec = window_lengths_sec[fi]
        
        # Calculate frequency smoothing (FieldTrip: cfg.tapsmofrq = frequency_smoothing * cfg.foi)
        tapsmofrq = frequency_smoothing * freq  # Frequency smoothing in Hz
        
        # Time-bandwidth product: NW = tapsmofrq * window_length / 2
        NW = tapsmofrq * window_length_sec / 2
        
        # Ensure NW is valid (must be > 0 and < n_window_samples/2)
        if NW <= 0 || NW >= n_window_samples / 2
            error("Invalid time-bandwidth product NW=$NW for frequency $freq Hz. NW must be > 0 and < $(n_window_samples/2). Window length: $window_length_sec s, tapsmofrq: $tapsmofrq Hz")
        end
        
        # Number of tapers: K = 2*NW - 1 (Shannon number)
        K = max(1, Int(floor(2 * NW - 1)))
        K = min(K, n_window_samples ÷ 2)  # Cap at half window length
        
        # Generate DPSS tapers
        dpss_result = DSP.dpss(n_window_samples, NW, K)
        
        # Handle different return types from dpss
        if dpss_result isa Tuple
            tapers = dpss_result[1]
        else
            tapers = dpss_result
        end
        
        # Convert to matrix format (n_window_samples × K)
        if tapers isa AbstractMatrix
            tapers_per_freq[fi] = Matrix{Float64}(tapers)
        elseif tapers isa AbstractVector
            # If it's a vector (happens when K=1), reshape to matrix
            tapers_per_freq[fi] = reshape(Vector{Float64}(tapers), length(tapers), 1)
        else
            error("DSP.dpss returned unexpected type: $(typeof(tapers)) for freq=$freq Hz, NW=$NW, K=$K")
        end
        n_tapers_per_freq[fi] = K
        inv_n_tapers_per_freq[fi] = 1.0 / K
        
        # Pre-compute tapered wavelets and their FFTs for this frequency
        # Tapered wavelet = taper * (cos + i*sin) at target frequency
        tapered_wavelet_ffts[fi] = Vector{Vector{ComplexF64}}(undef, K)
        
        # Create complex exponential at target frequency
        angle_in = range(-(n_window_samples-1)/2, (n_window_samples-1)/2, length = n_window_samples) .* (2π * freq / dat.sample_rate)
        cos_wav = cos.(angle_in)
        sin_wav = sin.(angle_in)
        
        # For each taper, create tapered wavelet and compute FFT
        for taper_idx = 1:K
            taper_col = @view tapers_per_freq[fi][:, taper_idx]
            # Tapered wavelet: taper * (cos + i*sin)
            tapered_wavelet = taper_col .* (cos_wav .+ im .* sin_wav)
            
            # Pad to n_samples_padded and compute FFT
            tapered_wavelet_padded = zeros(ComplexF64, n_samples_padded)
            tapered_wavelet_padded[1:n_window_samples] = tapered_wavelet
            tapered_wavelet_ffts[fi][taper_idx] = fft(tapered_wavelet_padded)
        end
    end

    # Pre-allocate reusable buffers for frequency-domain convolution
    data_padded = Matrix{Float64}(undef, n_samples_padded, n_trials)
    data_fft = Matrix{ComplexF64}(undef, n_samples_padded, n_trials)
    conv_result = Matrix{ComplexF64}(undef, n_samples_padded, n_trials)
    ifft_temp = Vector{ComplexF64}(undef, n_samples_padded)
    trial_signals_matrix = Matrix{Float64}(undef, n_samples_per_epoch, n_trials)

    # Initialize output structures - allocate appropriate type based on return_trials
    # Pre-compute shared time and freq columns (same for power and phase)
    time_col = repeat(times_out, inner = num_frex)
    freq_col = repeat(freqs, outer = n_times)
    
    if return_trials
        power_df = [DataFrame(time = time_col, freq = freq_col, copycols = false) for _ = 1:n_trials]
        phase_df = [DataFrame(time = time_col, freq = freq_col, copycols = false) for _ = 1:n_trials]
    else
        power_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
        phase_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
    end

    # Pre-allocate reusable output buffers (reused across all channels)
    if return_trials
        eegpower = zeros(Float64, num_frex, n_times, n_trials)
        eegconv = zeros(ComplexF64, num_frex, n_times, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times)
        eegconv = zeros(ComplexF64, num_frex, n_times)
    end

    # Process each selected channel
    for channel in selected_channels
        # Pre-extract all trial data for this channel into a matrix (n_samples × n_trials) for better cache locality
        trial_signals_matrix = Matrix{Float64}(undef, n_samples_per_epoch, n_trials)
        for trial_idx = 1:n_trials
            col = dat.data[trial_idx][!, channel]
            # Copy directly without intermediate Vector allocation
            @inbounds @simd for i = 1:n_samples_per_epoch
                trial_signals_matrix[i, trial_idx] = Float64(col[i])
            end
        end
        
        # Clear/initialize output buffers for this channel
        if return_trials
            fill!(eegpower, 0.0)
            fill!(eegconv, 0.0im)
        else
            fill!(eegpower, 0.0)
            fill!(eegconv, 0.0im)
        end

        # FieldTrip frequency-domain convolution approach for multitaper
        # FFT entire data once, then multiply by tapered wavelet FFTs and IFFT
        
        # Pad data to n_samples_padded (zero-padding at the end)
        fill!(data_padded, 0.0)
        data_padded[1:n_samples_per_epoch, :] = trial_signals_matrix
        
        # FFT entire padded data (batch process all trials at once)
        # Convert real to complex and perform batch FFT (out-of-place)
        data_fft .= complex.(data_padded)
        data_fft .= fft_plan_padded_batch * data_fft
        
        # Process each frequency
        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            n_tapers = n_tapers_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]
            
            # Initialize accumulation buffers for this frequency (accumulate across tapers)
            if return_trials
                power_accum = zeros(Float64, n_times, n_trials)
                complex_accum = zeros(ComplexF64, n_times, n_trials)
            else
                power_accum = zeros(Float64, n_times)
                complex_accum = zeros(ComplexF64, n_times)
            end
            
            # Process each taper and accumulate results
            for taper_idx = 1:n_tapers
                tapered_wavelet_fft = tapered_wavelet_ffts[fi][taper_idx]
                
                # Frequency-domain convolution: multiply data FFT by tapered wavelet FFT (broadcast across all trials)
                conv_result .= data_fft .* tapered_wavelet_fft
                
                # IFFT to get time-domain result (FFTW's IFFT doesn't normalize)
                @inbounds for trial = 1:n_trials
                    mul!(ifft_temp, ifft_plan_padded, view(conv_result, :, trial))
                    @simd for i = 1:n_samples_padded
                        conv_result[i, trial] = ifft_temp[i]
                    end
                end
                
                # Apply FieldTrip normalization: sqrt(2 ./ timwinsample)
                norm_factor = sqrt(2.0 / n_window_samples)  # FieldTrip: sqrt(2 ./ timwinsample)
                @inbounds @simd for i in eachindex(conv_result)
                    conv_result[i] *= norm_factor
                end
                
                # Extract requested time points and accumulate across tapers
                # Shift by half window length to account for FieldTrip's fftshift centering
                half_window = n_window_samples ÷ 2
                @inbounds for ti_idx = 1:n_times
                    sample_idx = time_indices[ti_idx]
                    # Adjust for shift: FieldTrip centers the result with fftshift
                    adjusted_idx = sample_idx + half_window
                    # Clamp to valid range
                    if adjusted_idx < 1
                        adjusted_idx = 1
                    elseif adjusted_idx > n_samples_padded
                        adjusted_idx = n_samples_padded
                    end
                    conv_vals = @view conv_result[adjusted_idx, :]
                    if return_trials
                        power_accum[ti_idx, :] .+= abs2.(conv_vals)
                        complex_accum[ti_idx, :] .+= conv_vals
                    else
                        power_accum[ti_idx] += sum(abs2, conv_vals)
                        complex_accum[ti_idx] += sum(conv_vals)
                    end
                end
            end
            
            # Average across tapers and store results
            if return_trials
                @inbounds for ti_idx = 1:n_times
                    eegpower[fi, ti_idx, :] .= power_accum[ti_idx, :] .* inv_n_tapers
                    eegconv[fi, ti_idx, :] .= complex_accum[ti_idx, :] .* inv_n_tapers
                end
            else
                inv_n_trials = 1.0 / n_trials
                @inbounds for ti_idx = 1:n_times
                    eegpower[fi, ti_idx] = power_accum[ti_idx] * inv_n_tapers * inv_n_trials
                    eegconv[fi, ti_idx] = complex_accum[ti_idx] * inv_n_tapers * inv_n_trials
                end
            end
        end

        # Apply edge filtering if requested (FieldTrip's "cone of influence")
        if filter_edges
            # Compute exact window lengths in samples (floating point) for edge filtering
            # FieldTrip uses: nsamplefreqoi = timwin(ifreqoi) .* fsample (exact floating point)
            # For multitaper: timwin = cycles / foi (adaptive window)
            window_lengths_samples_exact = [(cycles / freqs[fi]) * dat.sample_rate for fi = 1:num_frex]
            # Use unpadded data length for edge filtering (FieldTrip uses ndatsample = size(dat, 2), unpadded)
            _filter_edges!(eegpower, eegconv, num_frex, time_indices, window_lengths_samples_exact, n_samples_per_epoch) 
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = vec(eegpower[:, :, trial_idx])
                phase_df[trial_idx][!, channel] = vec(angle.(eegconv[:, :, trial_idx]))
            end
        else 
            eegpower ./= n_trials
            eegconv ./= n_trials
            power_df[!, channel] = vec(eegpower)
            phase_df[!, channel] = vec(angle.(eegconv))
        end
    end

    # Create and return appropriate data type
    return_type = return_trials ? TimeFreqEpochData : TimeFreqData
    return return_type(
        dat.file,
        dat.condition,
        dat.condition_name,
        power_df,
        phase_df,
        dat.layout,
        dat.sample_rate,
        :multitaper,
        nothing,  # baseline
        dat.analysis_info,
    )
end




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
- `window_size::Int=1024`: Size of the FFT window for spectral estimation (in samples)
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
    window_size::Int = 256,
    overlap::Real = 0.5,
    window_function::Function = DSP.hanning,
    max_freq::Union{Nothing,Real} = nothing,
)
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate overlap
    overlap < 0 || overlap >= 1 && error("`overlap` must be in range [0, 1), got $overlap")

    noverlap = Int(round(window_size * overlap))

    # Get frequency vector from first channel's first signal
    first_channel = selected_channels[1]
    first_signal = channel_data(dat, first_channel)
    pgram = DSP.welch_pgram(first_signal, window_size, noverlap; fs = dat.sample_rate, window = window_function)
    freqs = DSP.freq(pgram)

    # Initialize output DataFrame with frequency column
    spectrum_df = DataFrame(freq = freqs)

    # Process each channel
    for channel in selected_channels
        signal = channel_data(dat, channel)
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = dat.sample_rate, window = window_function)
        spectrum_df[!, channel] = DSP.power(pgram)
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = spectrum_df.freq .<= max_freq
        spectrum_df = spectrum_df[mask, :]
    end

    # Return SpectrumData type (dispatch handles condition info)
    condition_num, condition_name = condition_info(dat)
    return SpectrumData(
        dat.file,
        condition_num,
        condition_name,
        spectrum_df,
        dat.layout,
        dat.sample_rate,
        :welch,
        dat.analysis_info,
    )
end

