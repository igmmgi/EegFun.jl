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

    # Define cycles/sigma
    # MATLAB: s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex)
    if cycles isa Tuple
        cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
    else
        cycles_vec = fill(Float64(cycles), num_frex)
    end

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

    # Pre-compute convolution length and FFT plans for single trial (same for all channels)
    max_sigma = maximum(cycles_vec ./ (2 * pi .* freqs))
    max_hw = ceil(Int, 6 * max_sigma * sr) ÷ 2
    max_wl = max_hw * 2 + 1
    n_conv_trial = max_wl + n_samples_per_epoch - 1
    n_conv_pow2_trial = nextpow(2, n_conv_trial)
    
    # Pre-allocate FFT plans for single trial (smaller, faster) - same for all channels
    signal_padded_trial_template = zeros(ComplexF64, n_conv_pow2_trial)
    wavelet_padded_template = zeros(ComplexF64, n_conv_pow2_trial)
    conv_buffer_template = zeros(ComplexF64, n_conv_pow2_trial)
    p_fft_trial = plan_fft(signal_padded_trial_template, flags=FFTW.ESTIMATE)
    p_fft_wavelet_trial = plan_fft(wavelet_padded_template, flags=FFTW.ESTIMATE)
    p_ifft_trial = plan_ifft(conv_buffer_template, flags=FFTW.ESTIMATE)

    # Pre-allocate reusable buffers (same size for all channels, can be reused)
    signal_padded_trial = zeros(ComplexF64, n_conv_pow2_trial)
    eegfft_trial = zeros(ComplexF64, n_conv_pow2_trial)
    wavelet_padded = zeros(ComplexF64, n_conv_pow2_trial)
    wavelet_fft = zeros(ComplexF64, n_conv_pow2_trial)
    conv_buffer = zeros(ComplexF64, n_conv_pow2_trial)
    eegconv_trial = zeros(ComplexF64, n_conv_pow2_trial)

    # Pre-compute constants (same for all channels and frequencies)
    inv_sr = 1.0 / sr
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    hw_per_freq = Vector{Int}(undef, num_frex)
    valid_start_per_freq = Vector{Int}(undef, num_frex)
    for fi = 1:num_frex
        sigma = cycles_vec[fi] / (two_pi * freqs[fi])
        hw = ceil(Int, 6 * sigma * sr) ÷ 2
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
        wavelet_fft_freq = zeros(ComplexF64, n_conv_pow2_trial)
        mul!(wavelet_fft_freq, p_fft_wavelet_trial, wavelet_padded)
        wavelet_ffts[fi] = wavelet_fft_freq
    end

    # Process each selected channel - process trials separately for better performance
    for channel in selected_channels
        # Initialize output arrays - only for requested time points!
        if return_trials
            eegpower = zeros(Float64, num_frex, n_times, n_trials)
            eegconv = zeros(ComplexF64, num_frex, n_times, n_trials)
        else
            eegpower = zeros(Float64, num_frex, n_times)
            eegconv = zeros(ComplexF64, num_frex, n_times)
        end
        inv_n_trials = 1.0 / n_trials

        # Process each trial separately
        for trial_idx = 1:n_trials
            # Extract signal for this trial
            # signal_trial = dat_processed.data[trial_idx][!, channel]
            
            # FFT of trial signal - zero pad first, then copy signal
            fill!(signal_padded_trial, 0)
            signal_padded_trial[1:n_samples_per_epoch] .= dat_processed.data[trial_idx][!, channel]
            mul!(eegfft_trial, p_fft_trial, signal_padded_trial)

            # Loop through frequencies - reuse pre-computed wavelet FFTs
            for fi = 1:num_frex
                wavelet_fft = wavelet_ffts[fi]
                valid_start = valid_start_per_freq[fi]

                # Convolution (MATLAB: ifft(wavelet.*eegfft)) - use @simd for faster multiplication
                @inbounds @simd for i = 1:n_conv_pow2_trial
                    conv_buffer[i] = wavelet_fft[i] * eegfft_trial[i]
                end
                mul!(eegconv_trial, p_ifft_trial, conv_buffer)

                # Directly extract only requested time points from convolution result
                # Pre-compute conv_idx offsets for better performance
                @inbounds for ti in eachindex(time_indices)
                    conv_idx = valid_start + time_indices[ti] - 1
                    val = eegconv_trial[conv_idx]
                    if return_trials
                        eegpower[fi, ti, trial_idx] = abs2(val)
                        eegconv[fi, ti, trial_idx] = val
                    else
                        eegpower[fi, ti] += abs2(val) * inv_n_trials
                        eegconv[fi, ti] += val * inv_n_trials
                    end
                end
            end
        end

        if return_trials
            # Store each trial separately
            for trial_idx = 1:n_trials
                # Extract trial data and assign directly (vec flattens 2D to match long format)
                @views power_trial = eegpower[:, :, trial_idx]
                @views conv_trial = eegconv[:, :, trial_idx]
                power_dfs[trial_idx][!, channel] = vec(power_trial)
                phase_dfs[trial_idx][!, channel] = vec(angle.(conv_trial))
            end
        else
            # Already averaged during accumulation - just compute phase and assign directly
            power_df[!, channel] = vec(eegpower)
            phase_df[!, channel] = vec(angle.(eegconv))
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

"""
    tf_multitaper(dat::EpochData; 
                  lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
                  log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
                  window_length::Union{Nothing,Real}=nothing,
                  cycles::Union{Nothing,Real,Tuple{Real,Real}}=nothing,
                  frequency_smoothing::Union{Nothing,Real}=nothing,
                  time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
                  channel_selection::Function=channels(),
                  pad::Union{Nothing,Symbol}=nothing,
                  return_trials::Bool=false)

Multitaper time-frequency analysis using DPSS (Discrete Prolate Spheroidal Sequences) tapers (Cohen Chapter 16, equivalent to FieldTrip's 'mtmconvol' method).

Uses multiple orthogonal tapers (Slepian sequences) to reduce variance in spectral estimates. Supports both fixed and adaptive window lengths.

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
  - Example: `window_length=0.5` uses a fixed 0.5 second window for all frequencies
- `cycles::Union{Nothing,Real,Tuple{Real,Real}}=nothing`: Number of cycles for adaptive window mode.
  - Window length = cycles / frequency (in seconds)
  - Example: `cycles=5` uses 5 cycles for all frequencies (FieldTrip: `cfg.t_ftimwin = 5./cfg.foi`)
  - Example: `cycles=(3, 10)` uses log-spaced cycles from 3 to 10 across frequencies
  - **Exactly one of `window_length` or `cycles` must be specified.**
- `frequency_smoothing::Union{Nothing,Real}=nothing`: Frequency smoothing parameter (FieldTrip's `tapsmofrq`).
  - If `nothing`, uses `frequency_smoothing = 0.4 * frequency` (FieldTrip default: `cfg.tapsmofrq = 0.4 * cfg.foi`)
  - If a number, uses that value multiplied by frequency: `tapsmofrq = frequency_smoothing * frequency`
  - Example: `frequency_smoothing=0.4` matches FieldTrip's default
  - Controls time-bandwidth product: `NW = tapsmofrq * window_length / 2`
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

# Fixed window with custom frequency smoothing
tf_data = tf_multitaper(epochs; lin_freqs=(1, 30, 2), window_length=0.5, frequency_smoothing=0.4)

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
    window_length::Union{Nothing,Real} = nothing,
    cycles::Union{Nothing,Real,Tuple{Real,Real}} = nothing,
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

    # Default frequency smoothing (FieldTrip: cfg.tapsmofrq = 0.4 * cfg.foi)
    if isnothing(frequency_smoothing)
        frequency_smoothing = 0.4
    end
    if frequency_smoothing <= 0
        error("`frequency_smoothing` must be positive, got $frequency_smoothing")
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
    # FieldTrip: cfg.t_ftimwin = 5./cfg.foi means window_length = cycles / frequency (in seconds)
    if use_fixed_window
        # Fixed window mode: same window length for all frequencies
        window_lengths_sec = fill(Float64(window_length), num_frex)
    else
        # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
        if cycles isa Tuple
            cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
        else
            cycles_vec = fill(Float64(cycles), num_frex)
        end
        window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds
    end

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * sr)) for wl in window_lengths_sec]
    if any(n -> n < 2, n_window_samples_per_freq)
        min_samples = minimum(n_window_samples_per_freq)
        error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end

    # Determine window center positions (in samples) for each requested time point
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

    # Pre-compute DPSS tapers and FFT parameters for each frequency
    n_fft_per_freq = [nextpow(2, n_win) for n_win in n_window_samples_per_freq]
    tapers_per_freq = Vector{Matrix{Float64}}(undef, num_frex)  # Each element is (n_window_samples × n_tapers)
    n_tapers_per_freq = zeros(Int, num_frex)
    
    # For each frequency, find the FFT bin that corresponds to it
    freq_indices = zeros(Int, num_frex)
    
    for (fi, freq) in enumerate(freqs)
        n_window_samples = n_window_samples_per_freq[fi]
        window_length_sec = window_lengths_sec[fi]
        
        # Calculate frequency smoothing (FieldTrip: cfg.tapsmofrq = frequency_smoothing * cfg.foi)
        tapsmofrq = frequency_smoothing * freq  # Frequency smoothing in Hz
        
        # Time-bandwidth product: NW = tapsmofrq * window_length / 2
        # (matches old code: nw = frequency_smoothing * timewinidx / sample_rate / 2)
        NW = tapsmofrq * window_length_sec / 2
        
        # Ensure NW is valid (must be > 0 and < n_window_samples/2)
        if NW <= 0 || NW >= n_window_samples / 2
            error("Invalid time-bandwidth product NW=$NW for frequency $freq Hz. NW must be > 0 and < $(n_window_samples/2). Window length: $window_length_sec s, tapsmofrq: $tapsmofrq Hz")
        end
        
        # Number of tapers: K = 2*NW - 1 (but ensure at least 1, and cap at reasonable value)
        # (matches old code: n_tapers = max(1, floor(Int, 2 * nw - 1)))
        K = max(1, Int(floor(2 * NW - 1)))
        K = min(K, n_window_samples ÷ 2)  # Cap at half window length to avoid issues
        
        # If K=1 and NW is very small, dpss may have issues - try to get at least 2 tapers
        # by slightly increasing NW if possible
        if K == 1 && NW < 1.0
            # Try to get K=2 by adjusting NW slightly
            NW_adjusted = 1.5  # This will give K = floor(2*1.5 - 1) = 2
            K = 2
        else
            NW_adjusted = NW
        end
        
        # Generate DPSS tapers (dpss returns (tapers, eigenvalues) where tapers is a matrix)
        # Match old code: tapers = DSP.dpss(timewinidx, nw, n_tapers)
        # Note: When K=1, dpss may return a vector instead of a matrix
        dpss_result = DSP.dpss(n_window_samples, NW_adjusted, K)
        
        # Handle different return types from dpss
        if dpss_result isa Tuple
            tapers = dpss_result[1]
        else
            # If not a tuple, use directly (unlikely)
            tapers = dpss_result
        end
        
        # Convert to matrix format (n_window_samples × K)
        if tapers isa AbstractMatrix
            tapers_per_freq[fi] = Matrix{Float64}(tapers)
        elseif tapers isa AbstractVector
            # If it's a vector (happens when K=1), reshape to matrix (n_window_samples × 1)
            tapers_per_freq[fi] = reshape(Vector{Float64}(tapers), length(tapers), 1)
        elseif tapers isa Number
            # If it's a scalar, something went wrong - but try to create a matrix from it
            # This shouldn't happen, but if NW is very small, dpss might misbehave
            error("DSP.dpss returned a scalar ($tapers) instead of a matrix/vector. NW=$NW, K=$K, n_window_samples=$n_window_samples, freq=$freq Hz. Try increasing frequency_smoothing or window_length.")
        else
            error("DSP.dpss returned unexpected type: $(typeof(tapers)), value: $tapers. NW=$NW, K=$K, n_window_samples=$n_window_samples, freq=$freq Hz")
        end
        n_tapers_per_freq[fi] = K
        
        # Find FFT bin for this frequency
        n_fft = n_fft_per_freq[fi]
        freqs_fft = range(0.0, sr, length = n_fft + 1)[1:(n_fft÷2+1)]  # FFT frequencies (0 to Nyquist)
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
                tapers = tapers_per_freq[fi]  # (n_window_samples × n_tapers)
                n_tapers = n_tapers_per_freq[fi]
                freq_idx = freq_indices[fi]

                # Process each requested time point
                for (ti_idx, window_center) in enumerate(window_center_indices)
                    # Calculate window start and end indices
                    window_start = window_center - n_window_samples ÷ 2
                    window_end = window_start + n_window_samples - 1

                    # Extract windowed signal (handle edge cases)
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

                    # Apply each taper and compute FFT, then average across tapers
                    power_sum = 0.0
                    complex_sum = 0.0im
                    
                    for taper_idx = 1:n_tapers
                        # Apply taper
                        tapered_signal = windowed_signal .* tapers[:, taper_idx]
                        
                        # Pad to FFT length and compute FFT
                        signal_padded = zeros(ComplexF64, n_fft)
                        signal_padded[1:n_window_samples] = tapered_signal
                        signal_fft = fft(signal_padded)
                        
                        # Extract value at requested frequency
                        if freq_idx <= length(signal_fft)
                            complex_val = signal_fft[freq_idx]
                            # Normalize by window length squared (consistent with STFT)
                            power_sum += abs2(complex_val) / (n_window_samples^2)
                            complex_sum += complex_val
                        end
                    end
                    
                    # Average across tapers
                    eegpower_full[fi, ti_idx, trial_idx] = power_sum / n_tapers
                    eegconv_full[fi, ti_idx, trial_idx] = complex_sum / n_tapers
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
            :multitaper,
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
            :multitaper,
            nothing,  # baseline
            dat.analysis_info,
        )
    end
end

"""
    freq_spectrum(dat::EegData;
                  channel_selection::Function=channels(),
                  window_size::Int=1024,
                  overlap::Real=0.5,
                  window_function::Function=DSP.hanning,
                  max_freq::Union{Nothing,Real}=nothing)

Compute power spectrum using Welch's method for EEG data.

This function computes the frequency-domain power spectrum (no time dimension) using
Welch's method with overlapping windows. For EpochData, power is averaged across epochs.

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

# Plot
using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Frequency (Hz)", ylabel="Power (μV²/Hz)")
for ch in names(spectrum)[2:end]  # Skip 'freq' column
    lines!(ax, spectrum.freq, spectrum[!, ch], label=string(ch))
end
axislegend(ax)
```
"""
function freq_spectrum(
    dat::EpochData;
    channel_selection::Function = channels(),
    window_size::Int = 1024,
    overlap::Real = 0.5,
    window_function::Function = DSP.hanning,
    max_freq::Union{Nothing,Real} = nothing,
)
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate overlap
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    sr = Float64(dat.sample_rate)
    n_trials = n_epochs(dat)
    noverlap = Int(round(window_size * overlap))

    # Initialize output DataFrame
    spectrum_df = DataFrame()
    freqs_initialized = false

    # Process each channel
    for channel in selected_channels
        # Collect power spectra from all epochs
        power_spectra = Vector{Vector{Float64}}()

        for trial_idx = 1:n_trials
            epoch_df = dat.data[trial_idx]
            signal = Float64.(epoch_df[!, channel])

            # Compute power spectrum using Welch's method
            pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = sr, window = window_function)
            freqs, power = DSP.freq(pgram), DSP.power(pgram)

            # Initialize frequency column on first iteration
            if !freqs_initialized
                spectrum_df.freq = collect(Float64, freqs)
                freqs_initialized = true
            end

            push!(power_spectra, collect(Float64, power))
        end

        # Average power across epochs
        power_avg = mean(hcat(power_spectra...), dims = 2)[:, 1]
        spectrum_df[!, channel] = power_avg
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = spectrum_df.freq .<= max_freq
        spectrum_df = spectrum_df[mask, :]
    end

    # Return SpectrumData type
    return SpectrumData(
        dat.file,
        dat.condition,
        dat.condition_name,
        spectrum_df,
        dat.layout,
        dat.sample_rate,
        :welch,
        dat.analysis_info,
    )
end

function freq_spectrum(
    dat::ErpData;
    channel_selection::Function = channels(),
    window_size::Int = 1024,
    overlap::Real = 0.5,
    window_function::Function = DSP.hanning,
    max_freq::Union{Nothing,Real} = nothing,
)
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate overlap
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    sr = Float64(dat.sample_rate)
    noverlap = Int(round(window_size * overlap))

    # Initialize output DataFrame
    spectrum_df = DataFrame()
    freqs_initialized = false

    # Process each channel
    for channel in selected_channels
        signal = Float64.(dat.data[!, channel])

        # Compute power spectrum using Welch's method
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = sr, window = window_function)
        freqs, power = DSP.freq(pgram), DSP.power(pgram)

        # Initialize frequency column on first iteration
        if !freqs_initialized
            spectrum_df.freq = collect(Float64, freqs)
            freqs_initialized = true
        end

        spectrum_df[!, channel] = collect(Float64, power)
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = spectrum_df.freq .<= max_freq
        spectrum_df = spectrum_df[mask, :]
    end

    # Return SpectrumData type
    return SpectrumData(
        dat.file,
        dat.condition,
        dat.condition_name,
        spectrum_df,
        dat.layout,
        dat.sample_rate,
        :welch,
        dat.analysis_info,
    )
end

function freq_spectrum(
    dat::ContinuousData;
    channel_selection::Function = channels(),
    window_size::Int = 1024,
    overlap::Real = 0.5,
    window_function::Function = DSP.hanning,
    max_freq::Union{Nothing,Real} = nothing,
)
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate overlap
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    sr = Float64(dat.sample_rate)
    noverlap = Int(round(window_size * overlap))

    # Initialize output DataFrame
    spectrum_df = DataFrame()
    freqs_initialized = false

    # Process each channel
    for channel in selected_channels
        signal = Float64.(dat.data[!, channel])

        # Compute power spectrum using Welch's method
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = sr, window = window_function)
        freqs, power = DSP.freq(pgram), DSP.power(pgram)

        # Initialize frequency column on first iteration
        if !freqs_initialized
            spectrum_df.freq = collect(Float64, freqs)
            freqs_initialized = true
        end

        spectrum_df[!, channel] = collect(Float64, power)
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = spectrum_df.freq .<= max_freq
        spectrum_df = spectrum_df[mask, :]
    end

    # Return SpectrumData type
    # ContinuousData doesn't have condition/condition_name fields, use defaults
    return SpectrumData(
        dat.file,
        1,  # Default condition number
        "Continuous",  # Default condition name
        spectrum_df,
        dat.layout,
        dat.sample_rate,
        :welch,
        dat.analysis_info,
    )
end

