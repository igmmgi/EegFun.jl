"""
    tf_morlet(dat::EpochData; 
              lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
              log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
              cycles::Union{Real,Tuple{Real,Real}}=(3, 10),
              time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
              channel_selection::Function=channels(),
              pad::Union{Nothing,Symbol}=nothing,
              return_trials::Bool=false)

Time-frequency analysis using Morlet wavelets (Cohen Chapter 13).

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Linear frequency spacing as (start, stop, step).
  - Example: `lin_freqs=(2, 80, 2)` creates frequencies [2, 4, 6, ..., 80] Hz
- `log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing`: Logarithmic frequency spacing as (start, stop, number).
  - Example: `log_freqs=(2, 80, 30)` creates 30 log-spaced frequencies from 2 to 80 Hz
- **Exactly one of `lin_freqs` or `log_freqs` must be specified.**
- `cycles::Union{Real,Tuple{Real,Real}}=(3, 10)`: Number of cycles. Can be:
  - Single number: fixed cycles for all frequencies
  - Tuple (min, max): log-spaced cycles from min to max
- `time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Time points of interest as (start, stop, step) in seconds.
  - If `nothing`, uses all time points from the data
  - Example: `time_steps=(-0.5, 2.0, 0.01)` creates time points from -0.5 to 2.0 with 0.01s steps
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default)
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Example
```julia
# Log-spaced frequencies (30 frequencies from 2 to 80 Hz), single channel, averaged
tf_data = tf_morlet(epochs; log_freqs=(2, 80, 30), channel_selection=channels(:Cz))

# Linear-spaced frequencies with padding and individual trials
tf_epochs = tf_morlet(epochs; lin_freqs=(2, 80, 2), time_steps=(-0.5, 2.0, 0.01), 
    pad=:both, return_trials=true)
```
"""
function tf_morlet(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    cycles::Union{Real,Tuple{Real,Real}} = (3, 10),
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
)

    # Get selected channels using channel selection predicate
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    if isempty(selected_channels)
        error("No channels selected. Available channels: $(channel_labels(dat))")
    end

    # Validate frequency specification - exactly one must be provided
    if isnothing(lin_freqs) && isnothing(log_freqs)
        error("Either `lin_freqs` or `log_freqs` must be specified")
    end
    if !isnothing(lin_freqs) && !isnothing(log_freqs) 
        error("Only one of `lin_freqs` or `log_freqs` can be specified, not both")
    end

    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end

    # Get original data time range 
    times_original = time(dat)

    # Apply padding if requested 
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get sample rate and time vector from processed data
    times_processed = time(dat)
    n_samples_processed = n_samples(dat)  # Number of samples per epoch (may be padded)

    # Handle time_steps parameter - determine which time points to extract from results
    if isnothing(time_steps)
        time_indices, times_out = find_times(times_processed, times_original)
    else
        start_time, stop_time, step_time = time_steps
        
        # Check if requested range extends beyond processed data range and warn
        time_min_processed = minimum(times_processed)
        time_max_processed = maximum(times_processed)
        if start_time < time_min_processed || stop_time > time_max_processed
            @minimal_warning "Requested time range ($start_time to $stop_time seconds) extends beyond processed data range ($time_min_processed to $time_max_processed seconds). Clipping to available range."
        end
        
        time_steps_range = Float64(start_time):Float64(step_time):Float64(stop_time)
        time_indices, times_out = find_times(times_processed, time_steps_range)
        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Use processed data dimensions for convolution
    n_samples_original = n_samples_processed

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples_original  # Use full signal for convolution

    # Define frequencies based on user specification
    if !isnothing(log_freqs) # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
    else # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = range(Float64(min_freq), Float64(max_freq), step = Float64(step))
    end
    num_frex = length(freqs)  # Update num_frex for use in rest of function

    # Define cycles/sigma
    if cycles isa Tuple
        cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
    else
        cycles_vec = fill(Float64(cycles), num_frex)
    end

    # Initialize output structures based on return_trials flag
    n_times = length(times_out)

    # Initialize output structures 
    time_col = repeat(times_out, inner = num_frex)
    freq_col = repeat(freqs, outer = n_times)
    if return_trials
        power_df = [DataFrame(time = time_col, freq = freq_col, copycols = false) for _ = 1:n_trials]
        phase_df = [DataFrame(time = time_col, freq = freq_col, copycols = false) for _ = 1:n_trials]
    else
        power_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
        phase_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
    end

    # Pre-compute convolution length for batch processing (all trials concatenated)
    max_sigma = maximum(cycles_vec ./ (2 * pi .* freqs))
    max_hw = ceil(Int, 6 * max_sigma * dat.sample_rate) ÷ 2
    max_wl = max_hw * 2 + 1
    n_total = n_samples_per_epoch * n_trials  # Total samples across all trials
    n_conv_batch = max_wl + n_total - 1
    n_conv_pow2_batch = nextpow(2, n_conv_batch)
    
    # Pre-allocate reusable buffers for batch processing (same size for all channels)
    signal_padded_batch = zeros(ComplexF64, n_conv_pow2_batch)
    eegfft_batch = zeros(ComplexF64, n_conv_pow2_batch)
    wavelet_padded = zeros(ComplexF64, n_conv_pow2_batch)
    conv_buffer = zeros(ComplexF64, n_conv_pow2_batch)
    eegconv_batch = zeros(ComplexF64, n_conv_pow2_batch)
    
    # Create FFT plans using the actual buffers (plans keep references, so buffers must persist)
    p_fft_batch = plan_fft(signal_padded_batch, flags=FFTW.MEASURE)
    p_fft_wavelet_batch = plan_fft(wavelet_padded, flags=FFTW.MEASURE)
    p_ifft_batch = plan_ifft(conv_buffer, flags=FFTW.MEASURE)

    # Pre-compute constants (same for all channels and frequencies)
    inv_sr = 1.0 / dat.sample_rate 
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    hw_per_freq = Vector{Int}(undef, num_frex)
    valid_start_per_freq = Vector{Int}(undef, num_frex)
    for fi = 1:num_frex
        sigma = cycles_vec[fi] / (two_pi * freqs[fi])
        hw = ceil(Int, 6 * sigma * dat.sample_rate) ÷ 2
        wl = hw * 2 + 1
        hw_per_freq[fi] = hw
        valid_start_per_freq[fi] = hw + 1
        
        # Create wavelet directly in padded buffer
        fill!(wavelet_padded, 0)
        A = sqrt(1 / (sigma * sqrt_pi))
        two_pi_freq = two_pi * freqs[fi]
        inv_2sigma2 = 1.0 / (2 * sigma^2)
        @inbounds @simd for i = 1:wl
            t = (-wl / 2 + i - 1) * inv_sr
            t2 = t * t
            wavelet_padded[i] = A * exp(im * two_pi_freq * t) * exp(-t2 * inv_2sigma2)
        end

        # FFT of wavelet - compute once, reuse for all channels and trials
        wavelet_fft_freq = zeros(ComplexF64, n_conv_pow2_batch)
        mul!(wavelet_fft_freq, p_fft_wavelet_batch, wavelet_padded)
        wavelet_ffts[fi] = wavelet_fft_freq
    end

    # Pre-allocate reusable output buffers (reused across all channels)
    if return_trials
        eegpower = zeros(Float64, num_frex, n_times, n_trials)
        eegconv = zeros(ComplexF64, num_frex, n_times, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times)
        eegconv = zeros(ComplexF64, num_frex, n_times)
    end

    # Process each selected channel - batch process all trials together
    for channel in selected_channels

        # Concatenate all trials for this channel into one long signal
        fill!(signal_padded_batch, 0)
        offset = 0
        for trial_idx = 1:n_trials
            signal_padded_batch[offset+1:offset+n_samples_per_epoch] .= dat.data[trial_idx][!, channel]
            offset += n_samples_per_epoch
        end
        
        # Single FFT for all trials (much faster than n_trials separate FFTs)
        mul!(eegfft_batch, p_fft_batch, signal_padded_batch)

        # Loop through frequencies - reuse pre-computed wavelet FFTs
        for fi = 1:num_frex
            wavelet_fft = wavelet_ffts[fi]
            valid_start = valid_start_per_freq[fi]

            # Convolution (MATLAB: ifft(wavelet.*eegfft)) - use @simd for faster multiplication
            @inbounds @simd for i = 1:n_conv_pow2_batch
                conv_buffer[i] = wavelet_fft[i] * eegfft_batch[i]
            end
            mul!(eegconv_batch, p_ifft_batch, conv_buffer)

            # Extract results per trial using trial offsets
            for trial_idx = 1:n_trials
                t_offset = (trial_idx - 1) * n_samples_per_epoch
                @inbounds for ti in eachindex(time_indices)
                    conv_idx = valid_start + t_offset + time_indices[ti] - 1
                    val = eegconv_batch[conv_idx]
                    if return_trials
                        eegpower[fi, ti, trial_idx] = abs2(val)
                        eegconv[fi, ti, trial_idx] = val
                    else
                        eegpower[fi, ti] += abs2(val) 
                        eegconv[fi, ti] += val 
                    end
                end
            end
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                @views power_trial = eegpower[:, :, trial_idx]
                @views conv_trial = eegconv[:, :, trial_idx]
                power_df[trial_idx][!, channel] = vec(power_trial)
                phase_df[trial_idx][!, channel] = vec(angle.(conv_trial))
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
        :wavelet,
        nothing,  # baseline
        dat.analysis_info,
    )
end



"""
    tf_stft_fixed(dat::EpochData; 
            lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
            window_length::Real,
            overlap::Real=0.5,
            time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            channel_selection::Function=channels(),
            pad::Union{Nothing,Symbol}=nothing,
            return_trials::Bool=false)

Short-Time Fourier Transform (STFT) time-frequency analysis using a fixed-length sliding Hanning window (Cohen Chapter 15, equivalent to FieldTrip's 'hanning_fixed' method).

Uses a fixed window length for all frequencies, which provides consistent time resolution across the frequency range.

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Linear frequency spacing as (start, stop, step).
  - Example: `lin_freqs=(2, 80, 2)` creates frequencies [2, 4, 6, ..., 80] Hz
- `log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing`: Logarithmic frequency spacing as (start, stop, number).
  - Example: `log_freqs=(2, 80, 30)` creates 30 log-spaced frequencies from 2 to 80 Hz
- **Exactly one of `lin_freqs` or `log_freqs` must be specified.**
- `window_length::Real`: Fixed window length in seconds (same for all frequencies).
  - Example: `window_length=0.3` uses a fixed 0.3 second window for all frequencies
- `overlap::Real=0.5`: Overlap fraction between windows (0 to 1). Default is 0.5 (50% overlap).
- `time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Time points of interest as (start, stop, step) in seconds.
  - If `nothing`, uses all time points from the data
  - Example: `time_steps=(-0.5, 2.0, 0.01)` creates time points from -0.5 to 2.0 with 0.01s steps
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default). Time points where windows extend beyond data boundaries are excluded with a warning.
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Example
```julia
# Fixed window: Log-spaced frequencies with 0.3s fixed window
tf_data = tf_stft_fixed(epochs; log_freqs=(2, 80, 30), channel_selection=channels(:Cz), window_length=0.3)

# With padding and individual trials
tf_epochs = tf_stft_fixed(epochs; lin_freqs=(2, 80, 2), time_steps=(-0.5, 2.0, 0.01), 
                          window_length=0.3, overlap=0.5, pad=:both, return_trials=true)
```
"""
function tf_stft_fixed(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    window_length::Real,
    overlap::Real = 0.5,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
)
    # Get selected channels using channel selection predicate
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate frequency specification - exactly one must be provided
    if isnothing(lin_freqs) && isnothing(log_freqs) 
        error("Either `lin_freqs` or `log_freqs` must be specified")
    end
    if !isnothing(lin_freqs) && !isnothing(log_freqs) 
        error("Only one of `lin_freqs` or `log_freqs` can be specified, not both")
    end

    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end

    # Validate overlap
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    # Validate window_length
    if window_length <= 0
        error("`window_length` must be positive, got $window_length")
    end

    # Get original data time range (for when time_steps is nothing - use all original points)
    times_original = time(dat)

    # Apply padding if requested (mutating dat directly for consistency with tf_morlet)
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get sample rate and time vector from processed data
    times_processed = time(dat)

    # Handle time_steps parameter - determine which time points to extract from results
    # find_times already filters out-of-range times, so no need to check separately
    if isnothing(time_steps)
        time_indices, times_out = find_times(times_processed, times_original)
    else
        start_time, stop_time, step_time = time_steps
        time_steps_range = start_time:step_time:stop_time
        time_indices, times_out = find_times(times_processed, time_steps_range)
        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples(dat)

    # Define frequencies based on user specification
    if !isnothing(log_freqs) # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
    else # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = range(Float64(min_freq), Float64(max_freq), step = Float64(step))
    end
    num_frex = length(freqs)

    # Fixed window: same window length for all frequencies
    n_window_samples = Int(round(window_length * dat.sample_rate))
    if n_window_samples < 2
        error("Window length is too short. Minimum window size is 2 samples, got $n_window_samples samples.")
    end

    # Compute window positions for each requested time point and filter to valid ones
    half_window = n_window_samples ÷ 2
    window_starts = Vector{Int}(undef, length(times_out))
    valid_mask = BitVector(undef, length(times_out))
    
    for (ti_idx, t_requested) in enumerate(times_out)
        sample_idx = searchsortedfirst(times_processed, t_requested)
        window_start = sample_idx - half_window
        window_end = window_start + n_window_samples - 1
        window_starts[ti_idx] = window_start
        valid_mask[ti_idx] = (window_start >= 1 && window_end <= n_samples_per_epoch)
    end
    
    # Filter to only valid time points
    n_valid = count(valid_mask)
    if n_valid < length(times_out)
        n_excluded = length(times_out) - n_valid
        @minimal_warning "Excluding $n_excluded time point(s) where windows extend beyond data boundaries. Consider using `pad=:both` to avoid this."
        valid_indices = findall(valid_mask)
        times_out = times_out[valid_indices]
        time_indices = time_indices[valid_indices]
        window_starts = window_starts[valid_indices]
    end
    n_times = length(times_out)

    # Fixed window: all frequencies use the same window size and FFT size
    n_fft = nextpow(2, n_window_samples)
    hanning_window = DSP.hanning(n_window_samples)
    inv_n_window_sq = 1.0 / (n_window_samples * n_window_samples)
    
    # Pre-compute FFT plan for batch processing (2D matrix: n_fft × n_trials)
    template_2d = zeros(Float64, n_fft, n_trials)
    fft_plan_batch = plan_rfft(template_2d, 1, flags = FFTW.MEASURE)
    
    # Pre-allocate reusable buffers 
    signal_padded_batch = zeros(Float64, n_fft, n_trials)
    signal_fft_batch = zeros(ComplexF64, n_fft÷2+1, n_trials)

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
    
    # Pre-compute constants for bin index calculation
    max_bin = n_fft ÷ 2 + 1  # rFFT output size
    freq_to_bin_factor = n_fft / dat.sample_rate

    # Process each selected channel
    for channel in selected_channels

        # extract all trial data for this channel into a matrix (n_samples × n_trials) for batch processing
        trial_signals_matrix = Matrix{Float64}(undef, n_samples_per_epoch, n_trials)
        for trial_idx = 1:n_trials
            trial_signals_matrix[:, trial_idx] = dat.data[trial_idx][!, channel]
        end

        # Clear/initialize output buffers for this channel
        fill!(eegpower, 0.0)
        fill!(eegconv, 0.0im)

        # Fixed window - do 1 FFT per time point, batch all trials together
        @inbounds for ti_idx = 1:n_times
            window_start = window_starts[ti_idx]
            
            # Extract windowed signal, apply Hanning window, and pad for FFT
            fill!(signal_padded_batch, 0.0)  # Zero pad for FFT
            @inbounds @simd for trial = 1:n_trials
                for i = 1:n_window_samples
                    signal_padded_batch[i, trial] = trial_signals_matrix[window_start + i - 1, trial] * hanning_window[i]
                end
            end
            
            # Batch FFT: process all trials at once 
            mul!(signal_fft_batch, fft_plan_batch, signal_padded_batch)
            
            # Extract power and complex values for all frequencies from this single FFT
            @inbounds for (fi, freq) in enumerate(freqs)
                # Compute FFT bin index directly: round(freq * n_fft / sample_rate) + 1
                freq_idx = clamp(round(Int, freq * freq_to_bin_factor) + 1, 1, max_bin)
                complex_vals = @view signal_fft_batch[freq_idx, :]
                if return_trials
                    eegpower[fi, ti_idx, :] .= abs2.(complex_vals) .* inv_n_window_sq
                    eegconv[fi, ti_idx, :] .= complex_vals
                else
                    eegpower[fi, ti_idx] += sum(abs2, complex_vals) * inv_n_window_sq
                    eegconv[fi, ti_idx] += sum(complex_vals)
                end
            end
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
        :taper_fixed,
        nothing,  # baseline
        dat.analysis_info,
    )
end



"""
    tf_stft_adaptive(dat::EpochData; 
            lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
            cycles::Real,
            overlap::Real=0.5,
            time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            channel_selection::Function=channels(),
            pad::Union{Nothing,Symbol}=nothing,
            return_trials::Bool=false)

Short-Time Fourier Transform (STFT) time-frequency analysis using an adaptive-length sliding Hanning window (Cohen Chapter 15, equivalent to FieldTrip's 'hanning_adaptive' method).

Uses frequency-dependent window lengths (window_length = cycles / frequency), which provides better frequency resolution at lower frequencies and better time resolution at higher frequencies. This matches the trade-off used in Morlet wavelet analysis.

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Linear frequency spacing as (start, stop, step).
  - Example: `lin_freqs=(2, 80, 2)` creates frequencies [2, 4, 6, ..., 80] Hz
- `log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing`: Logarithmic frequency spacing as (start, stop, number).
  - Example: `log_freqs=(2, 80, 30)` creates 30 log-spaced frequencies from 2 to 80 Hz
- **Exactly one of `lin_freqs` or `log_freqs` must be specified.**
- `cycles::Real`: Number of cycles per frequency. Window length = cycles / frequency (in seconds).
  - Example: `cycles=5` uses 5 cycles for all frequencies (FieldTrip: `cfg.t_ftimwin = 5./cfg.foi`)
  - Lower frequencies will have longer windows, higher frequencies will have shorter windows
- `overlap::Real=0.5`: Overlap fraction between windows (0 to 1). Default is 0.5 (50% overlap).
  - Overlap is calculated per frequency based on that frequency's window length.
- `time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Time points of interest as (start, stop, step) in seconds.
  - If `nothing`, uses all time points from the data
  - Example: `time_steps=(-0.5, 2.0, 0.01)` creates time points from -0.5 to 2.0 with 0.01s steps
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default)
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Example
```julia
# Adaptive window: Linear-spaced frequencies with 5 cycles per frequency
tf_data = tf_stft_adaptive(epochs; lin_freqs=(2, 80, 2), cycles=5)

# With padding and individual trials
tf_epochs = tf_stft_adaptive(epochs; log_freqs=(2, 80, 30), time_steps=(-0.5, 2.0, 0.01), 
                              cycles=5, overlap=0.5, pad=:both, return_trials=true)
```
"""
function tf_stft_adaptive(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    cycles::Real,
    overlap::Real = 0.5,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
)

    # Get selected channels using channel selection predicate
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate frequency specification - exactly one must be provided
    if isnothing(lin_freqs) && isnothing(log_freqs) 
        error("Either `lin_freqs` or `log_freqs` must be specified")
    end
    if !isnothing(lin_freqs) && !isnothing(log_freqs) 
        error("Only one of `lin_freqs` or `log_freqs` can be specified, not both")
    end

    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end

    # Validate overlap
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    # Validate cycles parameter (adaptive window only)
    if cycles <= 0
        error("`cycles` must be positive, got $cycles")
    end

    # Get original data time range (for when time_steps is nothing - use all original points)
    times_original = time(dat)

    # Apply padding if requested (mutating dat directly for consistency with tf_morlet)
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get sample rate and time vector from processed data
    times_processed = time(dat)

    # Handle time_steps parameter - determine which time points to extract from results
    # After padding, processed data has extended time range - validate against processed data
    if isnothing(time_steps)
        time_indices, times_out = find_times(times_processed, times_original)
    else
        start_time, stop_time, step_time = time_steps
        # Check if requested range extends beyond processed data range and warn
        time_min_processed = minimum(times_processed)
        time_max_processed = maximum(times_processed)
        if start_time < time_min_processed || stop_time > time_max_processed
            @minimal_warning "Requested time range ($start_time to $stop_time seconds) extends beyond processed data range ($time_min_processed to $time_max_processed seconds). Clipping to available range."
        end
        time_steps_range = start_time:step_time:stop_time
        time_indices, times_out = find_times(times_processed, time_steps_range)
        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples(dat)

    # Define frequencies based on user specification
    if !isnothing(log_freqs) # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
    else # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = range(Float64(min_freq), Float64(max_freq), step = Float64(step))
    end
    num_frex = length(freqs)

    # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
    cycles_vec = fill(Float64(cycles), num_frex)
    window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds (matches FieldTrip: t_ftimwin = cycles./foi)

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
    if any(n -> n < 2, n_window_samples_per_freq)
        min_samples = minimum(n_window_samples_per_freq)
        error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end

    n_times = length(times_out)

    # Pre-compute Hanning windows for each frequency (adaptive windows - these vary)
    hanning_windows = [DSP.hanning(n_win) for n_win in n_window_samples_per_freq]
    
    # Use frequency-domain convolution (FieldTrip approach)
    # Determine padding length: pad to at least the data length and largest window
    # (FFTW is efficient for many sizes, not just powers of 2)
    max_window_samples = maximum(n_window_samples_per_freq)
    n_samples_padded = max(n_samples_per_epoch, max_window_samples)
    
    # Pre-compute complex wavelets and their FFTs for each frequency
    # Wavelet = Hanning window * complex exponential (cos + i*sin)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    for fi = 1:num_frex
        n_window_samples = n_window_samples_per_freq[fi]
        freq = freqs[fi]
        hanning_window = hanning_windows[fi]
        
        # Create complex wavelet: Hanning * (cos + i*sin) at frequency
        # Phase: center of wavelet has angle = 0 (FieldTrip convention)
        angle_in = range(-(n_window_samples-1)/2, (n_window_samples-1)/2, length = n_window_samples) .* (2π * freq / dat.sample_rate)
        cos_wav = hanning_window .* cos.(angle_in)
        sin_wav = hanning_window .* sin.(angle_in)
        wavelet = cos_wav .+ im .* sin_wav
        
        # Pad wavelet to match data FFT length and compute FFT
        # Place wavelet at the beginning (not centered) - centering happens through convolution
        wavelet_padded = zeros(ComplexF64, n_samples_padded)
        wavelet_padded[1:n_window_samples] = wavelet
        wavelet_ffts[fi] = fft(wavelet_padded)
    end
    
    # Plan for batch FFT of entire padded data (all trials at once)
    # Convert real to complex for batch FFT
    template_padded_batch = zeros(ComplexF64, n_samples_padded, n_trials)
    fft_plan_padded_batch = plan_fft(template_padded_batch, 1, flags = FFTW.MEASURE)  # FFT along first dimension
    
    # Plan for IFFT (same size, per trial)
    template_complex = zeros(ComplexF64, n_samples_padded)
    ifft_plan_padded = plan_ifft(template_complex, flags = FFTW.MEASURE)


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

    # For adaptive window: time_indices are already sample indices in processed data
    # Since we pad at the end, these indices work directly for the first n_samples_per_epoch samples

    # Pre-allocate reusable output buffers (reused across all channels)
    if return_trials
        eegpower = zeros(Float64, num_frex, n_times, n_trials)
        eegconv = zeros(ComplexF64, num_frex, n_times, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times)
        eegconv = zeros(ComplexF64, num_frex, n_times)
    end

    # Pre-allocate reusable buffers (reused across all channels and frequencies)
    data_padded = Matrix{Float64}(undef, n_samples_padded, n_trials)
    data_fft = Matrix{ComplexF64}(undef, n_samples_padded, n_trials)
    conv_result = Matrix{ComplexF64}(undef, n_samples_padded, n_trials)
    ifft_temp = Vector{ComplexF64}(undef, n_samples_padded)
    trial_signals_matrix = Matrix{Float64}(undef, n_samples_per_epoch, n_trials)

    # Process each selected channel
    for channel in selected_channels

        # Pre-extract all trial data for this channel into a matrix (n_samples × n_trials) for batch processing
        for trial_idx = 1:n_trials
            trial_signals_matrix[:, trial_idx] = dat.data[trial_idx][!, channel]
        end

        # Clear/initialize output buffers for this channel
        if return_trials
            fill!(eegpower, 0.0)
            fill!(eegconv, 0.0im)
        else
            fill!(eegpower, 0.0)
            fill!(eegconv, 0.0im)
        end

        # Adaptive window: FieldTrip frequency-domain convolution approach
        # FFT entire data once, then multiply by wavelet FFTs and IFFT
        
        # Pad data to n_samples_padded (zero-padding at the end)
        fill!(data_padded, 0.0)
        data_padded[1:n_samples_per_epoch, :] = trial_signals_matrix
        
        # FFT entire padded data (batch process all trials at once)
        # Convert real to complex and perform batch FFT (out-of-place)
        data_fft .= complex.(data_padded)
        data_fft .= fft_plan_padded_batch * data_fft
        
        # Process each frequency
        inv_n_samples_padded = 1.0 / n_samples_padded
        for fi = 1:num_frex
            # n_window_samples = n_window_samples_per_freq[fi]
            # wavelet_fft = wavelet_ffts[fi]
            
            # Frequency-domain convolution: multiply data FFT by wavelet FFT (broadcast across all trials)
            conv_result .= data_fft .* wavelet_ffts[fi]
            
            # IFFT to get time-domain result (FFTW's IFFT doesn't normalize, so we need to divide by n_samples_padded)
            @inbounds for trial = 1:n_trials
                # Direct column access - Julia optimizes this and FFTW is fast with contiguous arrays
                mul!(ifft_temp, ifft_plan_padded, view(conv_result, :, trial))
                @simd for i = 1:n_samples_padded
                    conv_result[i, trial] = ifft_temp[i] * inv_n_samples_padded
                end
            end
            
            # Extract requested time points (time_indices are sample indices in processed data)
            # The convolution result is already aligned with the data (no fftshift needed)
            # Normalization: for frequency-domain convolution, power = abs2(conv_result) / window_length
            inv_n_window = 1.0 / n_window_samples_per_freq[fi]
            @inbounds for ti_idx = 1:n_times
                sample_idx = time_indices[ti_idx]
                conv_vals = @view conv_result[sample_idx, :]
                if return_trials
                    eegpower[fi, ti_idx, :] .= abs2.(conv_vals) .* inv_n_window
                    eegconv[fi, ti_idx, :] .= conv_vals
                else
                    eegpower[fi, ti_idx] = sum(abs2, conv_vals) * inv_n_window
                    eegconv[fi, ti_idx] = sum(conv_vals)
                end
            end
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
        :hanning_adaptive,
        nothing,  # baseline
        dat.analysis_info,
    )
end








"""
    tf_multitaper(dat::EpochData; 
                  lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
                  log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
                  cycles::Real,
                  frequency_smoothing::Union{Nothing,Real}=nothing,
                  time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
                  channel_selection::Function=channels(),
                  pad::Union{Nothing,Symbol}=nothing,
                  return_trials::Bool=false)

Multitaper time-frequency analysis using DPSS (Discrete Prolate Spheroidal Sequences) tapers (Cohen Chapter 16, equivalent to FieldTrip's 'mtmconvol' method).

Uses multiple orthogonal tapers (Slepian sequences) to reduce variance in spectral estimates compared to single-taper methods. Uses adaptive window lengths (cycles per frequency) to match the time-frequency trade-off of Morlet wavelets and adaptive STFT.

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Linear frequency spacing as (start, stop, step).
  - Example: `lin_freqs=(2, 80, 2)` creates frequencies [2, 4, 6, ..., 80] Hz
- `log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing`: Logarithmic frequency spacing as (start, stop, number).
  - Example: `log_freqs=(2, 80, 30)` creates 30 log-spaced frequencies from 2 to 80 Hz
- **Exactly one of `lin_freqs` or `log_freqs` must be specified.**
- `cycles::Real`: Number of cycles per frequency. Window length = cycles / frequency (in seconds).
  - Example: `cycles=5` uses 5 cycles for all frequencies (FieldTrip: `cfg.t_ftimwin = 5./cfg.foi`)
  - Lower frequencies will have longer windows, higher frequencies will have shorter windows
- `frequency_smoothing::Union{Nothing,Real}=nothing`: Frequency smoothing parameter (FieldTrip's `tapsmofrq`).
  - If `nothing`, uses `frequency_smoothing = 0.4 * frequency` (FieldTrip default: `cfg.tapsmofrq = 0.4 * cfg.foi`)
  - If a number, uses that value multiplied by frequency: `tapsmofrq = frequency_smoothing * frequency`
  - Example: `frequency_smoothing=0.4` matches FieldTrip's default
  - Controls time-bandwidth product: `NW = tapsmofrq * window_length / 2`
  - Number of tapers used: `K = 2*NW - 1` (rounded down)
- `time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Time points of interest as (start, stop, step) in seconds.
  - If `nothing`, uses all time points from the data
  - Example: `time_steps=(-0.5, 2.0, 0.01)` creates time points from -0.5 to 2.0 with 0.01s steps
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default)
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Example
```julia
# Adaptive window with default frequency smoothing (0.4 * frequency)
tf_data = tf_multitaper(epochs; log_freqs=(2, 80, 30), cycles=5)

# Custom frequency smoothing
# FieldTrip equivalent: cfg.t_ftimwin = 5./cfg.foi, cfg.tapsmofrq = 0.4 * cfg.foi
tf_data = tf_multitaper(epochs; lin_freqs=(1, 30, 2), cycles=5, frequency_smoothing=0.4)
```
"""
function tf_multitaper(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    cycles::Real,
    frequency_smoothing::Union{Nothing,Real} = nothing,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
)
    # Get selected channels using channel selection predicate
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate frequency specification - exactly one must be provided
    isnothing(lin_freqs) && isnothing(log_freqs) && error("Either `lin_freqs` or `log_freqs` must be specified")
    !isnothing(lin_freqs) &&
        !isnothing(log_freqs) &&
        error("Only one of `lin_freqs` or `log_freqs` can be specified, not both")

    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
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

    # Apply padding if requested (mutating dat directly for consistency with tf_morlet)
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get sample rate and time vector from processed data
    times_processed = time(dat)
    n_samples_processed = n_samples(dat)  # Number of samples per epoch (may be padded)

    # Get original data time range (for when time_steps is nothing - use all original points)
    times_original = time(dat)

    # Handle time_steps parameter - determine which time points to extract from results
    # After padding, processed data has extended time range - validate against processed data
    time_min_processed = minimum(times_processed)
    time_max_processed = maximum(times_processed)

    if isnothing(time_steps)
        # Use all original time points (all in seconds)
        # Find indices in processed data that correspond to original time points
        time_indices, times_out = find_times(times_processed, times_original)
    else
        start_time, stop_time, step_time = time_steps
        
        # Check if requested range extends beyond processed data range and warn
        if start_time < time_min_processed || stop_time > time_max_processed
            @minimal_warning "Requested time range ($start_time to $stop_time seconds) extends beyond processed data range ($time_min_processed to $time_max_processed seconds). Clipping to available range."
        end
        
        time_steps_range = Float64(start_time):Float64(step_time):Float64(stop_time)
        time_indices, times_out = find_times(times_processed, time_steps_range)
        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Use processed data dimensions
    n_samples_original = n_samples_processed

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples_original

    # Define frequencies based on user specification
    if !isnothing(log_freqs) # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
    else # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = range(Float64(min_freq), Float64(max_freq), step = Float64(step))
    end
    num_frex = length(freqs)

    # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
    cycles_vec = fill(Float64(cycles), num_frex)
    window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
    if any(n -> n < 2, n_window_samples_per_freq)
        min_samples = minimum(n_window_samples_per_freq)
        error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end

    # Compute window center positions for each requested time point
    n_times = length(times_out)
    window_center_indices = Vector{Int}(undef, n_times)
    for (ti_idx, t_requested) in enumerate(times_out)
        sample_idx = searchsortedfirst(times_processed, t_requested)
        window_center_indices[ti_idx] = sample_idx
    end

    # Pre-compute DPSS tapers and FFT parameters for each frequency
    # Use actual window size for FFT (FFTW is efficient for many sizes, not just powers of 2)
    n_fft_per_freq = n_window_samples_per_freq
    tapers_per_freq = Vector{Matrix{Float64}}(undef, num_frex)  # Each element is (n_window_samples × n_tapers)
    n_tapers_per_freq = zeros(Int, num_frex)
    # Pre-compute normalization factors per frequency (avoid repeated computation)
    inv_n_window_samples_sq_per_freq = Vector{Float64}(undef, num_frex)
    inv_n_tapers_per_freq = Vector{Float64}(undef, num_frex)
    
    # For each frequency, find the FFT bin that corresponds to it
    freq_indices = zeros(Int, num_frex)
    
    # Pre-compute FFT plans for each unique FFT size (reuse for same sizes)
    unique_fft_sizes = unique(n_fft_per_freq)
    fft_plans_batch = Dict{Int, Any}()  # 2D plan for batch (n_fft × n_trials) - use Any for 2D plan type
    for n_fft in unique_fft_sizes
        # Plan for 2D matrix batch (n_fft × n_trials), FFT along dimension 1
        template_2d = zeros(ComplexF64, n_fft, n_trials)
        fft_plans_batch[n_fft] = plan_fft(template_2d, 1, flags = FFTW.MEASURE)
    end
    
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
        
        # Pre-compute normalization factors for this frequency
        inv_n_window_samples_sq_per_freq[fi] = 1.0 / (n_window_samples^2)
        inv_n_tapers_per_freq[fi] = 1.0 / n_tapers_per_freq[fi]
        
        # Compute FFT bin index directly (faster than searching)
        n_fft = n_fft_per_freq[fi]
        freq_to_bin_factor = n_fft / dat.sample_rate
        max_bin = n_fft ÷ 2 + 1  # rFFT output size
        freq_indices[fi] = clamp(round(Int, freq * freq_to_bin_factor) + 1, 1, max_bin)
    end

    # Pre-allocate reusable buffers per unique FFT size (reuse across frequencies with same sizes)
    unique_fft_sizes = unique(n_fft_per_freq)
    unique_window_sizes = unique(n_window_samples_per_freq)
    windowed_signal_buffers = Dict{Int, Matrix{Float64}}()  # Batch buffers: (n_window × n_trials)
    tapered_signal_buffers = Dict{Int, Matrix{Float64}}()   # Batch buffers: (n_window × n_trials)
    signal_padded_buffers = Dict{Int, Matrix{ComplexF64}}()  # Batch buffers: (n_fft × n_trials)
    signal_fft_buffers = Dict{Int, Matrix{ComplexF64}}()      # Batch buffers: (n_fft × n_trials)
    for n_window in unique_window_sizes
        windowed_signal_buffers[n_window] = zeros(Float64, n_window, n_trials)
        tapered_signal_buffers[n_window] = zeros(Float64, n_window, n_trials)
    end
    for n_fft in unique_fft_sizes
        signal_padded_buffers[n_fft] = zeros(ComplexF64, n_fft, n_trials)
        signal_fft_buffers[n_fft] = zeros(ComplexF64, n_fft, n_trials)
    end

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

        # Process each frequency (for adaptive windows, window size varies per frequency)
        # Batch process all trials together for each frequency-time-taper combination
        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            n_fft = n_fft_per_freq[fi]
            tapers = tapers_per_freq[fi]  # (n_window_samples × n_tapers)
            n_tapers = n_tapers_per_freq[fi]
            freq_idx = freq_indices[fi]
            p_fft_batch = fft_plans_batch[n_fft]  # Use batch plan for matrix operations

            # Reuse pre-allocated batch buffers for this frequency
            windowed_signal_batch = windowed_signal_buffers[n_window_samples]  # (n_window × n_trials)
            tapered_signal_batch = tapered_signal_buffers[n_window_samples]    # (n_window × n_trials)
            signal_padded_batch = signal_padded_buffers[n_fft]                  # (n_fft × n_trials)
            signal_fft_batch = signal_fft_buffers[n_fft]                        # (n_fft × n_trials)
            
            inv_n_window_samples_sq = inv_n_window_samples_sq_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]

            # Process each requested time point
            half_window = n_window_samples ÷ 2
            for ti_idx = 1:n_times
                window_center = window_center_indices[ti_idx]
                window_start = window_center - half_window
                window_end = window_start + n_window_samples - 1

                # Extract windowed signal for ALL trials at once
                # Check if window is valid (should be after filtering, but check for safety)
                if window_start >= 1 && window_end <= n_samples_per_epoch
                    @inbounds @simd for trial = 1:n_trials
                        for i = 1:n_window_samples
                            windowed_signal_batch[i, trial] = trial_signals_matrix[window_start + i - 1, trial]
                        end
                    end
                else
                    # Edge case: pad with zeros (shouldn't happen after filtering, but handle gracefully)
                    fill!(windowed_signal_batch, 0.0)
                    actual_start = max(1, window_start)
                    actual_end = min(n_samples_per_epoch, window_end)
                    if actual_start <= actual_end
                        copy_start = max(1, 1 - window_start + 1)
                        copy_end = copy_start + (actual_end - actual_start)
                        if copy_end <= n_window_samples && copy_start >= 1
                            @inbounds @simd for trial = 1:n_trials
                                for i = copy_start:copy_end
                                    windowed_signal_batch[i, trial] = trial_signals_matrix[actual_start + i - copy_start, trial]
                                end
                            end
                        end
                    end
                end

                # Reset power accumulation buffers for this time point
                if return_trials
                    fill!(view(eegpower, fi, ti_idx, :), 0.0)
                    fill!(view(eegconv, fi, ti_idx, :), 0.0im)
                else
                    power_sum_trials = zeros(Float64, n_trials)
                    complex_sum_trials = zeros(ComplexF64, n_trials)
                end
                
                # Apply each taper and compute FFT for all trials, then average across tapers
                for taper_idx = 1:n_tapers
                    # Apply taper to ALL trials at once
                    taper_col = @view tapers[:, taper_idx]
                    @inbounds @simd for trial = 1:n_trials
                        for i = 1:n_window_samples
                            tapered_signal_batch[i, trial] = windowed_signal_batch[i, trial] * taper_col[i]
                        end
                    end
                    
                    # Pad to FFT length for all trials (zero-pad if needed)
                    fill!(signal_padded_batch, 0.0)
                    @inbounds @simd for trial = 1:n_trials
                        for i = 1:n_window_samples
                            signal_padded_batch[i, trial] = tapered_signal_batch[i, trial]
                        end
                    end
                    
                    # Batch FFT: process all trials at once (ONE FFT instead of n_trials FFTs!)
                    mul!(signal_fft_batch, p_fft_batch, signal_padded_batch)
                    
                    # Extract value at requested frequency for all trials
                    if return_trials
                        @inbounds for trial = 1:n_trials
                            complex_val = signal_fft_batch[freq_idx, trial]
                            eegpower[fi, ti_idx, trial] += abs2(complex_val) * inv_n_window_samples_sq
                            eegconv[fi, ti_idx, trial] += complex_val
                        end
                    else
                        @inbounds for trial = 1:n_trials
                            complex_val = signal_fft_batch[freq_idx, trial]
                            power_sum_trials[trial] += abs2(complex_val) * inv_n_window_samples_sq
                            complex_sum_trials[trial] += complex_val
                        end
                    end
                end
                
                # Average across tapers and accumulate across trials
                if return_trials
                    @inbounds for trial = 1:n_trials
                        eegpower[fi, ti_idx, trial] *= inv_n_tapers
                        eegconv[fi, ti_idx, trial] *= inv_n_tapers
                    end
                else
                    inv_n_trials = 1.0 / n_trials
                    @inbounds for trial = 1:n_trials
                        power_avg = power_sum_trials[trial] * inv_n_tapers
                        complex_avg = complex_sum_trials[trial] * inv_n_tapers
                        eegpower[fi, ti_idx] += power_avg * inv_n_trials
                        eegconv[fi, ti_idx] += complex_avg * inv_n_trials
                    end
                end
            end
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

