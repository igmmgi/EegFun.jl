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

    # Get original data time range (for when time_steps is nothing - use all original points)
    # Extract this BEFORE potentially mutating dat
    times_original = time(dat)

    # Apply padding if requested (mutating dat directly since we extract time points later anyway)
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
    if !isnothing(log_freqs)
        # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
        # Broadcasting already allocates an array, no collect() needed
    else
        # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = range(Float64(min_freq), Float64(max_freq), step = Float64(step))
        # Range works fine for indexing, repeat(), and broadcasting - no collect() needed
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

    # Pre-compute convolution length and FFT plans for single trial (same for all channels)
    max_sigma = maximum(cycles_vec ./ (2 * pi .* freqs))
    max_hw = ceil(Int, 6 * max_sigma * dat.sample_rate) ÷ 2
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
        wavelet_fft_freq = zeros(ComplexF64, n_conv_pow2_trial)
        mul!(wavelet_fft_freq, p_fft_wavelet_trial, wavelet_padded)
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

    # Process each selected channel - process trials separately for better performance
    for channel in selected_channels

        # Process each trial separately
        for trial_idx = 1:n_trials
            
            # FFT of trial signal - zero pad first, then copy signal
            fill!(signal_padded_trial, 0)
            signal_padded_trial[1:n_samples_per_epoch] .= dat.data[trial_idx][!, channel]
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
                        eegpower[fi, ti] += abs2(val)
                        eegconv[fi, ti] += val
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

    # Validate window specification - exactly one must be provided
    if isnothing(window_length) && isnothing(cycles) 
        error("Either `window_length` (fixed) or `cycles` (adaptive) must be specified")
    end
    if !isnothing(window_length) && !isnothing(cycles) 
        error("Only one of `window_length` or `cycles` can be specified, not both")
    end

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

    # Determine window lengths (fixed or adaptive)
    # FieldTrip: cfg.t_ftimwin = 7./cfg.foi means window_length = cycles / frequency (in seconds)
    if use_fixed_window
        window_lengths_sec = fill(Float64(window_length), num_frex)
        # lets just match FieldTrip naming from https://www.fieldtriptoolbox.org/tutorial/sensor/timefrequencyanalysis/ 
        method_symbol = :hanning_fixed 
    else
        # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
        if cycles isa Tuple
            cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
        else
            cycles_vec = fill(Float64(cycles), num_frex)
        end
        window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds (matches FieldTrip: t_ftimwin = cycles./foi)
        # lets just match FieldTrip naming from https://www.fieldtriptoolbox.org/tutorial/sensor/timefrequencyanalysis/ 
        method_symbol = :hanning_adaptive
    end

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
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
    inv_n_window_sq_per_freq = [1.0 / (n_win * n_win) for n_win in n_window_samples_per_freq]
    
    # Pre-compute FFT plans and frequency bins for each unique FFT size (reuse for same sizes)
    # Use real-to-complex FFT (rfft) which is ~2x faster than complex FFT for real input data
    unique_fft_sizes = unique(n_fft_per_freq)
    fft_plans = Dict{Int, FFTW.rFFTWPlan{Float64, -1, false, 1}}()  # Real FFT plan
    fft_freqs_cache = Dict{Int, Vector{Float64}}()  # Cache FFT frequency vectors per FFT size
    for n_fft in unique_fft_sizes
        template = zeros(Float64, n_fft)  # Real input (not complex)
        # PATIENT flag: even more thorough than MEASURE, finds the absolute fastest algorithm
        # Real FFT is ~2x faster than complex FFT for real input data
        fft_plans[n_fft] = plan_rfft(template, flags = FFTW.PATIENT)
        # Pre-compute FFT frequency bins once per unique FFT size
        # rfft output size is n_fft÷2+1 (frequencies from 0 to Nyquist)
        fft_freqs_cache[n_fft] = collect(range(0.0, dat.sample_rate, length = n_fft + 1)[1:(n_fft÷2+1)])
    end
    
    # For each frequency, find the FFT bin that corresponds to it (use cached frequency vectors)
    freq_indices = Vector{Int}(undef, num_frex)
    for (fi, freq) in enumerate(freqs)
        n_fft = n_fft_per_freq[fi]
        freqs_fft = fft_freqs_cache[n_fft]  # Use cached frequency vector
        # Find nearest FFT frequency bin (avoid allocations from broadcasting)
        freq_idx = 1
        min_diff = abs(freqs_fft[1] - freq)
        @inbounds for i = 2:length(freqs_fft)
            diff = abs(freqs_fft[i] - freq)
            if diff < min_diff
                min_diff = diff
                freq_idx = i
            end
        end
        freq_indices[fi] = freq_idx
    end

    # Pre-allocate reusable buffers for each unique FFT size (match FFT plan sizes)
    max_window_samples = maximum(n_window_samples_per_freq)
    # Create buffers for each unique FFT size (allows direct mul! without views)
    # Use real buffers for input (rfft takes real input)
    signal_padded_buffers = Dict{Int, Vector{Float64}}()  # Real input buffer
    signal_fft_buffers = Dict{Int, Vector{ComplexF64}}()  # Complex output buffer
    for n_fft in unique_fft_sizes
        signal_padded_buffers[n_fft] = zeros(Float64, n_fft)  # Real input
        # rfft output size is n_fft÷2+1 (frequencies from 0 to Nyquist)
        signal_fft_buffers[n_fft] = zeros(ComplexF64, n_fft÷2+1)  # Complex output
    end
    signal_buffer = zeros(Float64, n_samples_per_epoch)  # Reusable buffer for trial signals

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

    # Pre-compute window start/end indices for each frequency and time point (avoid repeated calculations)
    window_starts_per_freq = Vector{Vector{Int}}(undef, num_frex)
    window_ends_per_freq = Vector{Vector{Int}}(undef, num_frex)
    needs_edge_per_freq = Vector{BitVector}(undef, num_frex)
    for fi = 1:num_frex
        n_window_samples = n_window_samples_per_freq[fi]
        half_window = n_window_samples ÷ 2
        window_starts = Vector{Int}(undef, n_times)
        window_ends = Vector{Int}(undef, n_times)
        needs_edge = BitVector(undef, n_times)
        @inbounds for (ti_idx, window_center) in enumerate(window_center_indices)
            window_start = window_center - half_window
            window_end = window_start + n_window_samples - 1
            window_starts[ti_idx] = window_start
            window_ends[ti_idx] = window_end
            needs_edge[ti_idx] = (window_start < 1 || window_end > n_samples_per_epoch)
        end
        window_starts_per_freq[fi] = window_starts
        window_ends_per_freq[fi] = window_ends
        needs_edge_per_freq[fi] = needs_edge
    end

    # Pre-allocate reusable output buffers (reused across all channels)
    if return_trials
        eegpower_full = zeros(Float64, num_frex, n_times, n_trials)
        eegconv_full = zeros(ComplexF64, num_frex, n_times, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times)
        eegconv = zeros(ComplexF64, num_frex, n_times)
    end

    # Process each selected channel
    for channel in selected_channels
        # Pre-extract all trial data for this channel (avoid repeated DataFrame lookups)
        trial_signals = Vector{Vector{Float64}}(undef, n_trials)
        for trial_idx = 1:n_trials
            trial_signals[trial_idx] = Vector{Float64}(dat.data[trial_idx][!, channel])
        end

        # Clear/initialize output buffers for this channel
        if return_trials
            fill!(eegpower_full, 0.0)
            fill!(eegconv_full, 0.0im)
        else
            fill!(eegpower, 0.0)
            fill!(eegconv, 0.0im)
        end

        # Process each frequency (for adaptive windows, window size varies per frequency)
        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            n_fft = n_fft_per_freq[fi]
            hanning_window = hanning_windows[fi]
            freq_idx = freq_indices[fi]
            inv_n_window_sq_val = inv_n_window_sq_per_freq[fi]
            fft_plan = fft_plans[n_fft]
            window_starts = window_starts_per_freq[fi]
            window_ends = window_ends_per_freq[fi]
            needs_edge = needs_edge_per_freq[fi]
            
            # Get buffers for this FFT size (pre-allocated, exact size match)
            signal_padded = signal_padded_buffers[n_fft]
            signal_fft = signal_fft_buffers[n_fft]

            # Process each trial
            for trial_idx = 1:n_trials
                signal_buffer = trial_signals[trial_idx]  # Use pre-extracted signal

                # Process each requested time point
                @inbounds for ti_idx = 1:n_times
                    window_start = window_starts[ti_idx]
                    window_end = window_ends[ti_idx]

                    # Extract windowed signal, apply Hanning window, and pad to FFT length in one pass
                    if needs_edge[ti_idx]
                        # Edge case: pad with zeros or use available data
                        # Zero out the entire padded buffer first
                        @inbounds @simd for i = 1:n_fft
                            signal_padded[i] = 0.0
                        end
                        actual_start = max(1, window_start)
                        actual_end = min(n_samples_per_epoch, window_end)
                        copy_start = max(1, 1 - window_start + 1)
                        copy_end = copy_start + (actual_end - actual_start)
                        if copy_end <= n_window_samples && copy_start >= 1
                            # Copy and apply Hanning window directly into padded buffer
                            @inbounds @simd for i = copy_start:copy_end
                                signal_padded[i] = signal_buffer[actual_start + i - copy_start] * hanning_window[i]
                            end
                        end
                    else
                        # Normal case: copy windowed signal, apply Hanning window, and pad in one pass
                        # Zero out tail first (if needed)
                        if n_fft > n_window_samples
                            @inbounds @simd for i = (n_window_samples+1):n_fft
                                signal_padded[i] = 0.0
                            end
                        end
                        # Copy and apply Hanning window directly into padded buffer
                        @inbounds @simd for i = 1:n_window_samples
                            signal_padded[i] = signal_buffer[window_start + i - 1] * hanning_window[i]
                        end
                    end

                    # Compute real-to-complex FFT using pre-computed plan (~2x faster than complex FFT)
                    # rfft takes real input and produces complex output of size n_fft÷2+1
                    mul!(signal_fft, fft_plan, signal_padded)

                    # Extract power and complex values at requested frequency
                    # FFT output is (n_fft÷2+1) frequencies from 0 to Nyquist
                    # Normalize power by window length squared (matches old tf_hanning implementation)
                    if freq_idx <= length(signal_fft)
                        @inbounds complex_val = signal_fft[freq_idx]
                        power_val = abs2(complex_val) * inv_n_window_sq_val
                        if return_trials
                            eegpower_full[fi, ti_idx, trial_idx] = power_val
                            eegconv_full[fi, ti_idx, trial_idx] = complex_val
                        else
                            eegpower[fi, ti_idx] += power_val
                            eegconv[fi, ti_idx] += complex_val
                        end
                    end
                end
            end
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = vec(eegpower_full[:, :, trial_idx])
                phase_df[trial_idx][!, channel] = vec(angle.(eegconv_full[:, :, trial_idx]))
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
        method_symbol,
        nothing,  # baseline
        dat.analysis_info,
    )
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
    if !isnothing(log_freqs)
        # Logarithmic spacing: (start, stop, number)
        min_freq, max_freq, num_frex = log_freqs
        freqs = exp.(range(log(Float64(min_freq)), log(Float64(max_freq)), length = num_frex))
        # Broadcasting already allocates an array, no collect() needed
    else
        # Linear spacing: (start, stop, step)
        min_freq, max_freq, step = lin_freqs
        freqs = range(Float64(min_freq), Float64(max_freq), step = Float64(step))
        # Range works fine for indexing, repeat(), and broadcasting - no collect() needed
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
    n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
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
    
    # Pre-compute FFT plans for each unique FFT size (reuse plans for same sizes)
    # Use FFTW.MEASURE for better performance (takes longer to create but FFTs are faster)
    unique_fft_sizes = unique(n_fft_per_freq)
    fft_plans = Dict{Int, FFTW.cFFTWPlan{ComplexF64, -1, false, 1}}()
    for n_fft in unique_fft_sizes
        template = zeros(ComplexF64, n_fft)
        fft_plans[n_fft] = plan_fft(template, flags = FFTW.MEASURE)
    end
    
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
        freqs_fft = range(0.0, dat.sample_rate, length = n_fft + 1)[1:(n_fft÷2+1)]  # FFT frequencies (0 to Nyquist)
        freq_idx = argmin(abs.(freqs_fft .- freq))
        freq_indices[fi] = freq_idx
    end

    # Pre-allocate reusable buffers (use maximum sizes needed, computed after tapers)
    max_window_samples = maximum(n_window_samples_per_freq)
    max_fft = maximum(n_fft_per_freq)
    windowed_signal_buffer = zeros(Float64, max_window_samples)
    tapered_signal_buffer = zeros(Float64, max_window_samples)
    signal_padded_buffer = zeros(ComplexF64, max_fft)
    signal_fft_buffer = zeros(ComplexF64, max_fft)
    signal_buffer = zeros(Float64, n_samples_per_epoch)  # Reusable buffer for trial signals

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
        eegpower_full = zeros(Float64, num_frex, n_times, n_trials)
        eegconv_full = zeros(ComplexF64, num_frex, n_times, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times)
        eegconv = zeros(ComplexF64, num_frex, n_times)
    end

    # Process each selected channel
    for channel in selected_channels
        # Clear/initialize output buffers for this channel
        if return_trials
            fill!(eegpower_full, 0.0)
            fill!(eegconv_full, 0.0im)
        else
            fill!(eegpower, 0.0)
            fill!(eegconv, 0.0im)
        end

        # Process each trial
        for trial_idx = 1:n_trials
            # Extract signal for this trial into pre-allocated buffer (avoid allocation)
            # Get column once to avoid repeated DataFrame lookups
            # epoch_df = dat.data[trial_idx]
            # col = epoch_df[!, channel]
            signal_buffer .= dat.data[trial_idx][!, channel]
            # Copy with type conversion in one pass (more efficient than Float64.(col) which allocates)
            # @inbounds @simd for i = 1:n_samples_per_epoch
            #     signal_buffer[i] = Float64(col[i])
            # end

            # Process each frequency (for adaptive windows, window size varies per frequency)
            for fi = 1:num_frex
                n_window_samples = n_window_samples_per_freq[fi]
                n_fft = n_fft_per_freq[fi]
                tapers = tapers_per_freq[fi]  # (n_window_samples × n_tapers)
                n_tapers = n_tapers_per_freq[fi]
                freq_idx = freq_indices[fi]
                p_fft = fft_plans[n_fft]

                # Get views into pre-allocated buffers for this frequency
                windowed_signal = @view windowed_signal_buffer[1:n_window_samples]
                tapered_signal = @view tapered_signal_buffer[1:n_window_samples]
                signal_padded = @view signal_padded_buffer[1:n_fft]
                signal_fft = @view signal_fft_buffer[1:n_fft]
                
                # Pre-compute whether padding is needed (same for all time points and tapers at this frequency)
                needs_padding = n_window_samples < n_fft
                padding_start = n_window_samples + 1

                # Process each requested time point
                for (ti_idx, window_center) in enumerate(window_center_indices)
                    # Calculate window start and end indices
                    window_start = window_center - n_window_samples ÷ 2
                    window_end = window_start + n_window_samples - 1

                    # Extract windowed signal into pre-allocated buffer
                    if window_start < 1 || window_end > n_samples_per_epoch
                        # Edge case: pad with zeros or use available data
                        fill!(windowed_signal, 0.0)
                        actual_start = max(1, window_start)
                        actual_end = min(n_samples_per_epoch, window_end)
                        copy_start = max(1, 1 - window_start + 1)
                        copy_end = copy_start + (actual_end - actual_start)
                        if copy_end <= n_window_samples && copy_start >= 1
                            @inbounds @simd for i = copy_start:copy_end
                                windowed_signal[i] = signal_buffer[actual_start + i - copy_start]
                            end
                        end
                    else
                        # Normal case: copy windowed signal (use @inbounds for speed)
                        @inbounds @simd for i = 1:n_window_samples
                            windowed_signal[i] = signal_buffer[window_start + i - 1]
                        end
                    end

                    # Apply each taper and compute FFT, then average across tapers
                    # Pre-compute normalization factor to avoid repeated division
                    inv_n_window_samples_sq = 1.0 / (n_window_samples^2)
                    inv_n_tapers = 1.0 / n_tapers
                    power_sum = 0.0
                    complex_sum = 0.0im
                    
                    for taper_idx = 1:n_tapers
                        # Apply taper (in-place multiplication into pre-allocated buffer)
                        taper_col = @view tapers[:, taper_idx]  # Views are cheap, create on the fly
                        @inbounds @simd for i = 1:n_window_samples
                            tapered_signal[i] = windowed_signal[i] * taper_col[i]
                        end
                        
                        # Pad to FFT length (zero out padding part, then copy tapered signal)
                        # Only zero the padding part if needed (pre-computed check)
                        if needs_padding
                            @inbounds @simd for i = padding_start:n_fft
                                signal_padded[i] = 0.0
                            end
                        end
                        @inbounds @simd for i = 1:n_window_samples
                            signal_padded[i] = tapered_signal[i]
                        end
                        
                        # Compute FFT using pre-computed plan (reuse buffer)
                        mul!(signal_fft, p_fft, signal_padded)
                        
                        # Extract value at requested frequency
                        if freq_idx <= length(signal_fft)
                            @inbounds complex_val = signal_fft[freq_idx]
                            # Normalize by window length squared (consistent with STFT)
                            power_sum += abs2(complex_val) * inv_n_window_samples_sq
                            complex_sum += complex_val
                        end
                    end
                    
                    # Average across tapers (use pre-computed inverse)
                    power_avg = power_sum * inv_n_tapers
                    complex_avg = complex_sum * inv_n_tapers
                    if return_trials
                        eegpower_full[fi, ti_idx, trial_idx] = power_avg
                        eegconv_full[fi, ti_idx, trial_idx] = complex_avg
                    else
                        eegpower[fi, ti_idx] += power_avg
                        eegconv[fi, ti_idx] += complex_avg
                    end
                end
            end
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = vec(eegpower_full[:, :, trial_idx])
                phase_df[trial_idx][!, channel] = vec(angle.(eegconv_full[:, :, trial_idx]))
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
                  window_size::Int=1024,
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
    dat::EpochData;
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

    n_trials = n_epochs(dat)
    noverlap = Int(round(window_size * overlap))

    # Use first epoch of first channel to get frequency vector
    first_channel = selected_channels[1]
    first_epoch_df = dat.data[1]
    first_signal = Float64.(first_epoch_df[!, first_channel])
    pgram_first = DSP.welch_pgram(first_signal, window_size, noverlap; fs = dat.sample_rate, window = window_function)

    # Initialize output DataFrame with frequency column
    freqs = DSP.freq(pgram_first)
    spectrum_df = DataFrame(freq = freqs)

    # Process each channel
    for channel in selected_channels
        # Accumulate power across epochs (avoid storing all spectra)
        power_sum = zeros(Float64, length(freqs))

        for trial_idx = 1:n_trials
            signal = dat.data[trial_idx][!, channel]

            # Compute power spectrum using Welch's method
            pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = dat.sample_rate, window = window_function)
            power_sum .+= DSP.power(pgram)
        end

        # Average power across epochs
        spectrum_df[!, channel] = power_sum ./ n_trials
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
    dat::SingleDataFrameEeg;
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
    if overlap < 0 || overlap >= 1
        error("`overlap` must be in range [0, 1), got $overlap")
    end

    noverlap = Int(round(window_size * overlap))

    # Compute frequencies once (same for all channels since they use same parameters)
    # Use first channel to get frequency vector
    first_channel = selected_channels[1]
    first_signal = Float64.(dat.data[!, first_channel])
    pgram_first = DSP.welch_pgram(first_signal, window_size, noverlap; fs = dat.sample_rate, window = window_function)
    freqs = collect(Float64, DSP.freq(pgram_first))

    # Initialize output DataFrame with frequency column
    spectrum_df = DataFrame(freq = freqs)

    # Process each channel
    for channel in selected_channels
        signal = dat.data[!, channel]

        # Compute power spectrum using Welch's method
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = dat.sample_rate, window = window_function)
        power = DSP.power(pgram)

        spectrum_df[!, channel] = power
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = spectrum_df.freq .<= max_freq
        spectrum_df = spectrum_df[mask, :]
    end

    # Return SpectrumData type
    # Handle condition/condition_name: ErpData has them, ContinuousData doesn't
    condition_num = hasproperty(dat, :condition) ? dat.condition : 1
    condition_name = hasproperty(dat, :condition_name) ? dat.condition_name : "Continuous"
    
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

