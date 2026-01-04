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

    # Apply padding if requested (non-mutating)
    dat_processed = isnothing(pad) ? dat : mirror(dat, pad)

    # Get sample rate and time vector from processed data
    sr = Float64(dat_processed.sample_rate)
    times_processed = time(dat_processed)
    n_samples_processed = n_samples(dat_processed)  # Number of samples per epoch (may be padded)

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
        # Generate time points from tuple (start, stop, step) - all in SECONDS
        start_time, stop_time, step_time = time_steps
        # Generate range of time values in seconds
        time_steps_vec = collect(Float64, range(Float64(start_time), Float64(stop_time), step = Float64(step_time)))
        # If the last value is significantly less than stop_time, add it
        if !isempty(time_steps_vec) && (stop_time - time_steps_vec[end]) > step_time / 2
            push!(time_steps_vec, Float64(stop_time))
        end

        # Validate that requested times (in seconds) are within the processed data range (which may be padded)
        if start_time < time_min_processed || stop_time > time_max_processed
            @minimal_warning "Requested time range ($start_time to $stop_time seconds) extends beyond processed data range ($time_min_processed to $time_max_processed seconds). Clipping to available range."
        end

        # Filter time_steps_vec to only include points within processed data range
        time_steps_vec_filtered = Base.filter(t -> time_min_processed <= t <= time_max_processed, time_steps_vec)

        if isempty(time_steps_vec_filtered)
            error(
                "No valid time points found in requested range ($start_time to $stop_time seconds) within processed data range ($time_min_processed to $time_max_processed seconds)",
            )
        end

        # Find nearest time points in processed data
        time_indices, times_out = find_times(times_processed, time_steps_vec_filtered)

        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Use processed data dimensions for convolution
    n_samples_original = n_samples_processed

    # Get number of trials/epochs
    n_trials = n_epochs(dat_processed)
    n_samples_per_epoch = n_samples_original  # Use full signal for convolution

    # Define frequencies based on user specification
    if !isnothing(log_freqs)
        # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
        freqs = collect(Float64, freqs)
    else
        # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = collect(Float64, range(Float64(min_freq), Float64(max_freq), step = Float64(step)))
    end
    num_frex = length(freqs)  # Update num_frex for use in rest of function

    # Define wavelet time window (MATLAB: -1:1/EEG.srate:1)
    wavelet_time = -1:(1/sr):1
    n_wavelet = length(wavelet_time)
    half_of_wavelet_size = (n_wavelet - 1) / 2

    # Define cycles/sigma
    # MATLAB: s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex)
    if cycles isa Tuple
        cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
    else
        cycles_vec = fill(Float64(cycles), num_frex)
    end
    sigma = cycles_vec ./ (2 * pi .* freqs)

    # Define convolution parameters (same for all channels)
    n_data = n_samples_per_epoch * n_trials
    n_convolution = n_wavelet + n_data - 1
    n_conv_pow2 = nextpow(2, n_convolution)

    # Initialize output structures based on return_trials flag
    n_times = length(times_out)

    if return_trials
        # Initialize per-trial DataFrames
        power_dfs = [DataFrame() for _ = 1:n_trials]
        phase_dfs = [DataFrame() for _ = 1:n_trials]
        for trial_idx = 1:n_trials
            power_dfs[trial_idx].time = repeat(times_out, inner = num_frex)
            power_dfs[trial_idx].freq = repeat(freqs, outer = n_times)
            phase_dfs[trial_idx].time = copy(power_dfs[trial_idx].time)
            phase_dfs[trial_idx].freq = copy(power_dfs[trial_idx].freq)
        end
    else
        # Initialize single averaged DataFrame
        power_df = DataFrame()
        power_df.time = repeat(times_out, inner = num_frex)  # Each time repeated for all freqs
        power_df.freq = repeat(freqs, outer = n_times)       # All freqs repeated for each time
        phase_df = DataFrame()
        phase_df.time = copy(power_df.time)
        phase_df.freq = copy(power_df.freq)
    end

    # Process each selected channel
    for channel in selected_channels
        # Extract full signal data for convolution (use all time points)
        signal = zeros(Float64, n_samples_original * n_trials)
        for (trial_idx, epoch_df) in enumerate(dat_processed.data)
            signal[(trial_idx-1)*n_samples_original+1:trial_idx*n_samples_original] = epoch_df[!, channel]
        end

        # Get FFT of data 
        signal_padded = zeros(ComplexF64, n_conv_pow2)
        signal_padded[1:n_data] = signal
        eegfft = fft(signal_padded)

        # Initialize output power and phase matrices (frequencies × time × trials) - use full time resolution
        # Store complex values for proper phase averaging
        eegpower_full = zeros(Float64, num_frex, n_samples_per_epoch, n_trials)
        eegconv_full = zeros(ComplexF64, num_frex, n_samples_per_epoch, n_trials)  # Store complex for phase

        # Loop through frequencies
        for fi = 1:num_frex
            # Create wavelet (MATLAB: sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))))
            A = sqrt(1 / (sigma[fi] * sqrt(pi)))
            wavelet = A .* exp.(2im * pi * freqs[fi] .* wavelet_time) .* exp.(-wavelet_time .^ 2 ./ (2 * (sigma[fi]^2)))

            # FFT of wavelet
            wavelet_padded = zeros(ComplexF64, n_conv_pow2)
            wavelet_padded[1:n_wavelet] = wavelet
            wavelet_fft = fft(wavelet_padded)

            # Convolution (MATLAB: ifft(wavelet.*eegfft))
            eegconv = ifft(wavelet_fft .* eegfft)

            # Extract valid part (MATLAB: eegconv(1:n_convolution) then eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size))
            eegconv = eegconv[1:n_convolution]
            eegconv = eegconv[Int(half_of_wavelet_size)+1:end-Int(half_of_wavelet_size)]

            # Reshape to (n_samples_per_epoch × n_trials)
            eegconv_reshaped = reshape(eegconv, n_samples_per_epoch, n_trials)

            # Store power and complex values for each trial
            eegpower_full[fi, :, :] = abs2.(eegconv_reshaped)
            eegconv_full[fi, :, :] = eegconv_reshaped
        end

        # Extract only requested time points from full power and complex matrices
        eegpower = eegpower_full[:, time_indices, :]  # (num_frex × n_times × n_trials)
        eegconv = eegconv_full[:, time_indices, :]    # (num_frex × n_times × n_trials)

        if return_trials
            # Store each trial separately
            for trial_idx = 1:n_trials
                # Extract trial data: (num_frex × n_times)
                power_trial = eegpower[:, :, trial_idx]
                phase_trial = angle.(eegconv[:, :, trial_idx])  # Compute phase from complex values

                # Reshape to long format: [freq1_t1, freq2_t1, ..., freqN_t1, freq1_t2, ...]
                power_values = vec(power_trial)
                phase_values = vec(phase_trial)

                power_dfs[trial_idx][!, channel] = power_values
                phase_dfs[trial_idx][!, channel] = phase_values
            end
        else
            # Average across trials
            eegpower_avg = mean(eegpower, dims = 3)  # (num_frex × n_times × 1)
            eegphase_avg = angle.(mean(eegconv, dims = 3))  # Mean of complex values, then angle

            # Reshape to long format
            power_values = vec(eegpower_avg)
            phase_values = vec(eegphase_avg)

            power_df[!, channel] = power_values
            phase_df[!, channel] = phase_values
        end
    end

    # No need to unmirror DataFrames - we already extracted only the original time range

    # Create and return appropriate data type
    if return_trials
        return TimeFreqEpochData(
            dat.file,
            dat.condition,
            dat.condition_name,
            power_dfs,
            phase_dfs,
            dat.layout,
            dat.sample_rate,
            :wavelet,
            nothing,  # baseline
            dat.analysis_info,
        )
    else
        return TimeFreqData(
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
end

"""
    tf_stft(dat::EpochData; 
            lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
            window_length::Union{Nothing,Real}=nothing,
            cycles::Union{Nothing,Real,Tuple{Real,Real}}=nothing,
            overlap::Real=0.5,
            time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            channel_selection::Function=channels(),
            pad::Union{Nothing,Symbol}=nothing,
            return_trials::Bool=false)

Short-Time Fourier Transform (STFT) time-frequency analysis using a sliding Hanning window (Cohen Chapter 15, equivalent to FieldTrip's 'hanning' method).

Supports both fixed and adaptive window lengths:
- **Fixed window** (FieldTrip's 'hanning_fixed'): Use `window_length` to specify a fixed window size in seconds
- **Adaptive window** (FieldTrip's 'hanning_adaptive'): Use `cycles` to specify frequency-dependent window sizes

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Linear frequency spacing as (start, stop, step).
  - Example: `lin_freqs=(2, 80, 2)` creates frequencies [2, 4, 6, ..., 80] Hz
- `log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing`: Logarithmic frequency spacing as (start, stop, number).
  - Example: `log_freqs=(2, 80, 30)` creates 30 log-spaced frequencies from 2 to 80 Hz
- **Exactly one of `lin_freqs` or `log_freqs` must be specified.**
- `window_length::Union{Nothing,Real}=nothing`: Fixed window length in seconds (for fixed window mode).
  - If `nothing`, uses adaptive window mode (requires `cycles`).
  - Example: `window_length=0.3` uses a fixed 0.3 second window for all frequencies
- `cycles::Union{Nothing,Real,Tuple{Real,Real}}=nothing`: Number of cycles for adaptive window mode.
  - Can be a single number (fixed cycles for all frequencies) or a tuple (min, max) for log-spaced cycles.
  - Window length = cycles / frequency (in seconds)
  - Example: `cycles=5` uses 5 cycles for all frequencies
  - Example: `cycles=(3, 10)` uses log-spaced cycles from 3 to 10 across frequencies
  - **Exactly one of `window_length` or `cycles` must be specified.**
- `overlap::Real=0.5`: Overlap fraction between windows (0 to 1). Default is 0.5 (50% overlap).
  - Note: For adaptive windows, overlap is calculated per frequency based on that frequency's window length.
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
# Fixed window: Log-spaced frequencies with 0.3s fixed window
tf_data = tf_stft(epochs; log_freqs=(2, 80, 30), channel_selection=channels(:Cz), window_length=0.3)

# Adaptive window: Linear-spaced frequencies with 5 cycles per frequency
tf_data = tf_stft(epochs; lin_freqs=(2, 80, 2), cycles=5)

# Adaptive window with variable cycles and padding
tf_epochs = tf_stft(epochs; lin_freqs=(2, 80, 2), time_steps=(-0.5, 2.0, 0.01), 
                    cycles=(3, 10), overlap=0.5, pad=:both, return_trials=true)
```
"""
function tf_stft(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    window_length::Union{Nothing,Real} = nothing,
    cycles::Union{Nothing,Real,Tuple{Real,Real}} = nothing,
    overlap::Real = 0.5,
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

    # Validate overlap
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    # Validate window specification - exactly one must be provided
    isnothing(window_length) && isnothing(cycles) && error("Either `window_length` (fixed) or `cycles` (adaptive) must be specified")
    !isnothing(window_length) &&
        !isnothing(cycles) &&
        error("Only one of `window_length` or `cycles` can be specified, not both")

    # Determine mode
    use_fixed_window = !isnothing(window_length)
    if use_fixed_window
        if window_length <= 0
            error("`window_length` must be positive, got $window_length")
        end
    else
        if cycles isa Tuple
            if cycles[1] <= 0 || cycles[2] <= 0 || cycles[1] >= cycles[2]
                error("`cycles` tuple must be (min, max) with 0 < min < max, got $cycles")
            end
        elseif cycles <= 0
            error("`cycles` must be positive, got $cycles")
        end
    end

    # Apply padding if requested (non-mutating)
    dat_processed = isnothing(pad) ? dat : mirror(dat, pad)

    # Get sample rate and time vector from processed data
    sr = Float64(dat_processed.sample_rate)
    times_processed = time(dat_processed)
    n_samples_processed = n_samples(dat_processed)  # Number of samples per epoch (may be padded)

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
        # Generate time points from tuple (start, stop, step) - all in SECONDS
        start_time, stop_time, step_time = time_steps
        # Generate range of time values in seconds
        time_steps_vec = collect(Float64, range(Float64(start_time), Float64(stop_time), step = Float64(step_time)))
        # If the last value is significantly less than stop_time, add it
        if !isempty(time_steps_vec) && (stop_time - time_steps_vec[end]) > step_time / 2
            push!(time_steps_vec, Float64(stop_time))
        end

        # Validate that requested times (in seconds) are within the processed data range (which may be padded)
        if start_time < time_min_processed || stop_time > time_max_processed
            @minimal_warning "Requested time range ($start_time to $stop_time seconds) extends beyond processed data range ($time_min_processed to $time_max_processed seconds). Clipping to available range."
        end

        # Filter time_steps_vec to only include points within processed data range
        time_steps_vec_filtered = Base.filter(t -> time_min_processed <= t <= time_max_processed, time_steps_vec)

        if isempty(time_steps_vec_filtered)
            error(
                "No valid time points found in requested range ($start_time to $stop_time seconds) within processed data range ($time_min_processed to $time_max_processed seconds)",
            )
        end

        # Find nearest time points in processed data
        time_indices, times_out = find_times(times_processed, time_steps_vec_filtered)

        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Use processed data dimensions
    n_samples_original = n_samples_processed

    # Get number of trials/epochs
    n_trials = n_epochs(dat_processed)
    n_samples_per_epoch = n_samples_original

    # Define frequencies based on user specification
    if !isnothing(log_freqs)
        # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
        freqs = collect(Float64, freqs)
    else
        # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = collect(Float64, range(Float64(min_freq), Float64(max_freq), step = Float64(step)))
    end
    num_frex = length(freqs)

    # Determine window lengths (fixed or adaptive)
    # FieldTrip: cfg.t_ftimwin = 7./cfg.foi means window_length = cycles / frequency (in seconds)
    if use_fixed_window
        # Fixed window mode: same window length for all frequencies (FieldTrip: fixed t_ftimwin)
        window_lengths_sec = fill(Float64(window_length), num_frex)
        method_symbol = :hanning_fixed
    else
        # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
        if cycles isa Tuple
            cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
        else
            cycles_vec = fill(Float64(cycles), num_frex)
        end
        window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds (matches FieldTrip: t_ftimwin = cycles./foi)
        method_symbol = :hanning_adaptive
    end

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * sr)) for wl in window_lengths_sec]
    if any(n -> n < 2, n_window_samples_per_freq)
        min_samples = minimum(n_window_samples_per_freq)
        error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end

    # Determine window center positions (in samples) for each requested time point
    # For each requested time, find the window that centers on it
    n_times = length(times_out)
    window_center_indices = zeros(Int, n_times)
    for (ti_idx, t_requested) in enumerate(times_out)
        # Find the sample index in processed data that corresponds to this time
        sample_idx = searchsortedfirst(times_processed, t_requested)
        if sample_idx > 1 && abs(times_processed[sample_idx-1] - t_requested) < abs(times_processed[min(sample_idx, length(times_processed))] - t_requested)
            sample_idx = sample_idx - 1
        end
        sample_idx = min(sample_idx, length(times_processed))
        window_center_indices[ti_idx] = sample_idx
    end

    # Pre-compute FFT parameters for each frequency (for adaptive windows, these vary)
    n_fft_per_freq = [nextpow(2, n_win) for n_win in n_window_samples_per_freq]
    hanning_windows = [DSP.hanning(n_win) for n_win in n_window_samples_per_freq]
    
    # For each frequency, find the FFT bin that corresponds to it
    freq_indices = zeros(Int, num_frex)
    for (fi, freq) in enumerate(freqs)
        n_fft = n_fft_per_freq[fi]
        freqs_fft = range(0.0, sr, length = n_fft + 1)[1:(n_fft÷2+1)]  # FFT frequencies (0 to Nyquist)
        # Find nearest FFT frequency bin
        freq_idx = argmin(abs.(freqs_fft .- freq))
        freq_indices[fi] = freq_idx
    end

    # Initialize output structures based on return_trials flag
    if return_trials
        # Initialize per-trial DataFrames
        power_dfs = [DataFrame() for _ = 1:n_trials]
        phase_dfs = [DataFrame() for _ = 1:n_trials]
        for trial_idx = 1:n_trials
            power_dfs[trial_idx].time = repeat(times_out, inner = num_frex)
            power_dfs[trial_idx].freq = repeat(freqs, outer = n_times)
            phase_dfs[trial_idx].time = copy(power_dfs[trial_idx].time)
            phase_dfs[trial_idx].freq = copy(power_dfs[trial_idx].freq)
        end
    else
        # Initialize single averaged DataFrame
        power_df = DataFrame()
        power_df.time = repeat(times_out, inner = num_frex)  # Each time repeated for all freqs
        power_df.freq = repeat(freqs, outer = n_times)       # All freqs repeated for each time
        phase_df = DataFrame()
        phase_df.time = copy(power_df.time)
        phase_df.freq = copy(power_df.freq)
    end

    # Process each selected channel
    for channel in selected_channels
        # Initialize output power and complex matrices (frequencies × time × trials)
        # Store complex values for proper phase averaging
        eegpower_full = zeros(Float64, num_frex, n_times, n_trials)
        eegconv_full = zeros(ComplexF64, num_frex, n_times, n_trials)

        # Process each trial
        for trial_idx = 1:n_trials
            # Extract signal for this trial
            epoch_df = dat_processed.data[trial_idx]
            signal = Float64.(epoch_df[!, channel])

            # Process each frequency (for adaptive windows, window size varies per frequency)
            for fi = 1:num_frex
                n_window_samples = n_window_samples_per_freq[fi]
                n_fft = n_fft_per_freq[fi]
                hanning_window = hanning_windows[fi]
                freq_idx = freq_indices[fi]

                # Process each requested time point
                for (ti_idx, window_center) in enumerate(window_center_indices)
                    # Calculate window start and end indices
                    window_start = window_center - n_window_samples ÷ 2
                    window_end = window_start + n_window_samples - 1

                    # Check bounds and handle edge cases
                    if window_start < 1 || window_end > n_samples_per_epoch
                        # Edge case: pad with zeros or use available data
                        windowed_signal = zeros(Float64, n_window_samples)
                        actual_start = max(1, window_start)
                        actual_end = min(n_samples_per_epoch, window_end)
                        copy_start = max(1, 1 - window_start + 1)
                        copy_end = copy_start + (actual_end - actual_start)
                        if copy_end <= n_window_samples && copy_start >= 1
                            windowed_signal[copy_start:copy_end] = signal[actual_start:actual_end]
                        end
                    else
                        # Normal case: extract windowed signal
                        windowed_signal = signal[window_start:window_end]
                    end

                    # Apply Hanning window
                    windowed_signal .*= hanning_window

                    # Pad to FFT length and compute FFT
                    signal_padded = zeros(ComplexF64, n_fft)
                    signal_padded[1:n_window_samples] = windowed_signal
                    signal_fft = fft(signal_padded)

                    # Extract power and complex values at requested frequency
                    # FFT output is (n_fft÷2+1) frequencies from 0 to Nyquist
                    # Normalize power by window length squared (matches old tf_hanning implementation)
                    if freq_idx <= length(signal_fft)
                        complex_val = signal_fft[freq_idx]
                        eegpower_full[fi, ti_idx, trial_idx] = abs2(complex_val) / (n_window_samples^2)
                        eegconv_full[fi, ti_idx, trial_idx] = complex_val
                    end
                end
            end
        end

        if return_trials
            # Store each trial separately
            for trial_idx = 1:n_trials
                # Extract trial data: (num_frex × n_times)
                power_trial = eegpower_full[:, :, trial_idx]
                phase_trial = angle.(eegconv_full[:, :, trial_idx])  # Compute phase from complex values

                # Reshape to long format: [freq1_t1, freq2_t1, ..., freqN_t1, freq1_t2, ...]
                power_values = vec(power_trial)
                phase_values = vec(phase_trial)

                power_dfs[trial_idx][!, channel] = power_values
                phase_dfs[trial_idx][!, channel] = phase_values
            end
        else
            # Average across trials
            eegpower_avg = mean(eegpower_full, dims = 3)  # (num_frex × n_times × 1)
            eegphase_avg = angle.(mean(eegconv_full, dims = 3))  # Mean of complex values, then angle

            # Reshape to long format
            power_values = vec(eegpower_avg)
            phase_values = vec(eegphase_avg)

            power_df[!, channel] = power_values
            phase_df[!, channel] = phase_values
        end
    end

    # Create and return appropriate data type
    if return_trials
        return TimeFreqEpochData(
            dat.file,
            dat.condition,
            dat.condition_name,
            power_dfs,
            phase_dfs,
            dat.layout,
            dat.sample_rate,
            method_symbol,
            nothing,  # baseline
            dat.analysis_info,
        )
    else
        return TimeFreqData(
            dat.file,
            dat.condition,
            dat.condition_name,
            power_df,
            phase_df,
            dat.layout,
            dat.sample_rate,
            method_symbol,
            nothing,  # baseline
            dat.analysis_info,
        )
    end
end


