"""
    tf_morlet(dat::EpochData; 
              lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
              log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
              cycles::Union{Real,Tuple{Real,Real}}=(3, 10),
              time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
              channel_selection::Function=channels(),
              pad::Union{Nothing,Symbol}=nothing,
              return_trials::Bool=false,
              filter_edges::Bool=true)

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
- `filter_edges::Bool=true`: If `true` (default), filters out edge regions where the wavelet extends beyond the data

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
    filter_edges::Bool = true,
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
    n_samples_original_unpadded = n_samples(dat)  # Store original unpadded length for edge filtering

    # Apply padding if requested 
    if !isnothing(pad)
        dat = mirror(dat, pad)
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
        
        time_steps_range = start_time:step_time:stop_time
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

    # Pre-compute convolution length for single trial processing
    max_sigma = maximum(cycles_vec ./ (2 * pi .* freqs))
    max_hw = ceil(Int, 6 * max_sigma * dat.sample_rate) ÷ 2
    max_wl = max_hw * 2 + 1
    n_conv = max_wl + n_samples_per_epoch - 1
    n_conv_pow2 = nextpow(2, n_conv)
    
    # Pre-allocate reusable buffers for single trial processing
    signal_padded = zeros(ComplexF64, n_conv_pow2)
    eegfft = zeros(ComplexF64, n_conv_pow2)
    wavelet_padded = zeros(ComplexF64, n_conv_pow2)
    conv_buffer = zeros(ComplexF64, n_conv_pow2)
    eegconv_buffer = zeros(ComplexF64, n_conv_pow2)
    
    # Create FFT plans using the actual buffers (plans keep references, so buffers must persist)
    p_fft = plan_fft(signal_padded, flags=FFTW.MEASURE)
    p_fft_wavelet = plan_fft(wavelet_padded, flags=FFTW.MEASURE)
    p_ifft = plan_ifft(eegconv_buffer, flags=FFTW.MEASURE)

    # Pre-compute constants (same for all channels and frequencies)
    inv_sr = 1.0 / dat.sample_rate 
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    hw_per_freq = Vector{Int}(undef, num_frex)
    wl_per_freq = Vector{Int}(undef, num_frex)  # Store actual wavelet length for edge filtering
    valid_start_per_freq = Vector{Int}(undef, num_frex)
    # Pre-compute convolution indices for each frequency to avoid computation in inner loop
    conv_indices_per_freq = Vector{Vector{Int}}(undef, num_frex)

    for fi = 1:num_frex

        sigma = cycles_vec[fi] / (two_pi * freqs[fi])
        hw = ceil(Int, 6 * sigma * dat.sample_rate) ÷ 2
        wl = hw * 2 + 1
        hw_per_freq[fi] = hw
        wl_per_freq[fi] = wl
        valid_start = hw + 1
        valid_start_per_freq[fi] = valid_start
        
        # Pre-compute convolution indices for this frequency
        conv_indices_per_freq[fi] = [valid_start + sample_idx - 1 for sample_idx in time_indices]
        
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
        wavelet_fft_freq = zeros(ComplexF64, n_conv_pow2)
        mul!(wavelet_fft_freq, p_fft_wavelet, wavelet_padded)
        wavelet_ffts[fi] = wavelet_fft_freq
    end

    # Pre-allocate reusable output buffers (reused across all channels)
    # Always initialize to zeros for accumulation (we'll set invalid regions to NaN later if filtering edges)
    if return_trials
        eegpower = zeros(Float64, num_frex, n_times, n_trials)
        eegconv = zeros(ComplexF64, num_frex, n_times, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times)
        eegconv = zeros(ComplexF64, num_frex, n_times)
    end

    norm_factor = sqrt(2.0 / dat.sample_rate)
    
    # Pre-compute taper lengths for edge filtering (if needed) - avoid recomputing
    taper_lengths_samples_exact = filter_edges ? [Float64(wl_per_freq[fi]) for fi = 1:num_frex] : nothing

    # Process each selected channel and each trial separately
    for channel in selected_channels
        # Reset output arrays for this channel
        fill!(eegpower, 0)
        fill!(eegconv, 0)
        
        for trial_idx = 1:n_trials
            
            # Copy single trial data to padded buffer
            fill!(signal_padded, 0)
            signal_padded[1:n_samples_per_epoch] .= dat.data[trial_idx][!, channel]
            
            # FFT for this trial
            mul!(eegfft, p_fft, signal_padded)

            # Loop through frequencies - reuse pre-computed wavelet FFTs
            for fi = 1:num_frex

                # Convolution (MATLAB: ifft(wavelet.*eegfft)) - use @simd for faster multiplication
                @inbounds @simd for i = 1:n_conv_pow2
                    conv_buffer[i] = wavelet_ffts[fi][i] * eegfft[i]
                end
                mul!(eegconv_buffer, p_ifft, conv_buffer)
                
                # Apply norm_factor and extract in one pass - use pre-computed indices
                conv_indices = conv_indices_per_freq[fi]
                @inbounds for ti in eachindex(time_indices)
                    val = eegconv_buffer[conv_indices[ti]] * norm_factor
                    
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
            if filter_edges
                _filter_edges!(eegpower, eegconv, num_frex, time_indices, taper_lengths_samples_exact, n_samples_original_unpadded)
            end
            # Pre-allocate phase vectors to avoid repeated allocations
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = vec(@view eegpower[:, :, trial_idx])
                # Compute angle directly without intermediate view
                phase_df[trial_idx][!, channel] = vec(angle.(@view eegconv[:, :, trial_idx]))
            end
        else
            eegpower ./= n_trials
            eegconv ./= n_trials
            if filter_edges
                _filter_edges!(eegpower, eegconv, num_frex, time_indices, taper_lengths_samples_exact, n_samples_original_unpadded)
            end
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


# helper function to filter edge regions 
function _filter_edges!(
    eegpower::AbstractArray,
    eegconv::AbstractArray,
    num_frex::Int,
    time_indices::AbstractVector{Int},
    window_lengths_per_freq::Union{Vector{Int}, Vector{Float64}},
    n_samples_per_epoch::Int,
)
    for fi = 1:num_frex

        half_nsamplefreqoi = window_lengths_per_freq[fi] / 2.0
        min_valid_threshold = half_nsamplefreqoi
        max_valid_threshold = n_samples_per_epoch - half_nsamplefreqoi
        
        for ti in eachindex(time_indices)
            sample_idx = time_indices[ti]
            if !(sample_idx >= min_valid_threshold && sample_idx < max_valid_threshold)
                if ndims(eegpower) == 3  # return_trials = true
                    @views eegpower[fi, ti, :] .= NaN
                    @views eegconv[fi, ti, :] .= NaN * im
                else  # return_trials = false
                    eegpower[fi, ti] = NaN
                    eegconv[fi, ti] = NaN * im
                end
            end
        end
    end
end





"""
    tf_stft(dat::EpochData; 
            lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
            window_length::Union{Nothing,Real}=nothing,
            cycles::Union{Nothing,Real}=nothing,
            overlap::Real=0.5,
            time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
            channel_selection::Function=channels(),
            pad::Union{Nothing,Symbol}=nothing,
            return_trials::Bool=false,
            filter_edges::Bool=true)

Short-Time Fourier Transform (STFT) time-frequency analysis using a sliding Hanning window (equivalent to FieldTrip's 'mtmconvol' method).

Supports both fixed-length windows (consistent time resolution) and adaptive windows (constant cycles per frequency).

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing`: Linear frequency spacing as (start, stop, step).
  - Example: `lin_freqs=(2, 80, 2)` creates frequencies [2, 4, 6, ..., 80] Hz
- `log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing`: Logarithmic frequency spacing as (start, stop, number).
  - Example: `log_freqs=(2, 80, 30)` creates 30 log-spaced frequencies from 2 to 80 Hz
- **Exactly one of `lin_freqs` or `log_freqs` must be specified.**
- `window_length::Union{Nothing,Real}=nothing`: Fixed window length in seconds (same for all frequencies).
  - Example: `window_length=0.3` uses a fixed 0.3 second window for all frequencies
  - **Exactly one of `window_length` or `cycles` must be specified.**
- `cycles::Union{Nothing,Real}=nothing`: Number of cycles per frequency (adaptive window).
  - Example: `cycles=7` uses 7 cycles per frequency (window length = 7/frequency)
  - **Exactly one of `window_length` or `cycles` must be specified.**
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
- `filter_edges::Bool=true`: If `true` (default), filters out edge regions where the window extends beyond the data
  (FieldTrip's "cone of influence" - marks edges as NaN). If `false`, uses all convolution results including edges.

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Example
```julia
# Fixed window: Log-spaced frequencies with 0.3s fixed window
tf_data = tf_stft(epochs; log_freqs=(2, 80, 30), channel_selection=channels(:Cz), window_length=0.3)

# Adaptive window: 7 cycles per frequency (FieldTrip equivalent)
tf_data = tf_stft(epochs; lin_freqs=(2, 30, 1), channel_selection=channels(:Cz), cycles=7, time_steps=(-0.5, 1.5, 0.05))
```
"""
function tf_stft(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    window_length::Union{Nothing,Real} = nothing,
    cycles::Union{Nothing,Real} = nothing,
    overlap::Real = 0.5,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
    filter_edges::Bool = true,
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

    # Validate window_length and cycles - exactly one must be provided
    if isnothing(window_length) && isnothing(cycles)
        error("Either `window_length` (fixed window) or `cycles` (adaptive window) must be specified")
    end
    if !isnothing(window_length) && !isnothing(cycles)
        error("Only one of `window_length` or `cycles` can be specified, not both")
    end
    if !isnothing(window_length) && window_length <= 0
        error("`window_length` must be positive, got $window_length")
    end
    if !isnothing(cycles) && cycles <= 0
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

    # Determine window mode: fixed or adaptive
    if !isnothing(window_length)
        # Fixed window mode: same window length for all frequencies
        n_window_samples = Int(round(window_length * dat.sample_rate))
        if n_window_samples < 2
            error("Window length is too short. Minimum window size is 2 samples, got $n_window_samples samples.")
        end
        n_window_samples_per_freq = fill(n_window_samples, num_frex)
        is_fixed = true
    else
        # Adaptive window mode: window length = cycles / frequency (FieldTrip: cfg.t_ftimwin = cycles./cfg.foi)
        cycles_vec = fill(Float64(cycles), num_frex)
        window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds (matches FieldTrip: t_ftimwin = cycles./foi)
        n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
        if any(n -> n < 2, n_window_samples_per_freq)
            min_samples = minimum(n_window_samples_per_freq)
            error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
        end
        is_fixed = false
    end

    n_times = length(times_out)

    # Pre-compute Hanning windows for each frequency
    # FieldTrip normalizes Hanning window by Frobenius norm: tap = tap./norm(tap, 'fro')
    # For fixed window, all frequencies use the same window, so create it once
    if is_fixed
        n_win = n_window_samples_per_freq[1]  # All the same for fixed window
        hanning_window_normalized = DSP.hanning(n_win) ./ norm(DSP.hanning(n_win), 2)  # L2 norm for vector = Frobenius norm
        hanning_windows = [hanning_window_normalized for _ = 1:num_frex]  # Reuse same window for all frequencies
    else
        # Adaptive window: each frequency has different window size
        hanning_windows = Vector{Vector{Float64}}(undef, num_frex)
        for fi = 1:num_frex
            n_win = n_window_samples_per_freq[fi]
            hanning_windows[fi] = DSP.hanning(n_win) ./ norm(DSP.hanning(n_win), 2)  # L2 norm for vector = Frobenius norm
        end
    end
    
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

        # FieldTrip frequency-domain convolution approach (works for both fixed and adaptive windows)
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
            n_window_samples = n_window_samples_per_freq[fi]
            
            # Frequency-domain convolution: multiply data FFT by wavelet FFT (broadcast across all trials)
            conv_result .= data_fft .* wavelet_ffts[fi]
            
            # IFFT to get time-domain result (FFTW's IFFT doesn't normalize)
            @inbounds for trial = 1:n_trials
                # Direct column access - Julia optimizes this and FFTW is fast with contiguous arrays
                mul!(ifft_temp, ifft_plan_padded, view(conv_result, :, trial))
                @simd for i = 1:n_samples_padded
                    conv_result[i, trial] = ifft_temp[i]  # No IFFT normalization (matching fixed window)
                end
            end
            
            # Apply FieldTrip normalization: sqrt(2 ./ timwinsample)
            norm_factor = sqrt(2.0 / n_window_samples)  # FieldTrip: sqrt(2 ./ timwinsample)
            @inbounds @simd for i in eachindex(conv_result)
                conv_result[i] *= norm_factor
            end
            
            # Extract requested time points (time_indices are sample indices in processed data)
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
                    eegpower[fi, ti_idx, :] .= abs2.(conv_vals)
                    eegconv[fi, ti_idx, :] .= conv_vals
                else
                    eegpower[fi, ti_idx] = sum(abs2, conv_vals)
                    eegconv[fi, ti_idx] = sum(conv_vals)
                end
            end
        end

        # Apply edge filtering if requested (FieldTrip's "cone of influence")
        if filter_edges
            # Compute exact window lengths in samples (floating point) for edge filtering
            # FieldTrip uses: nsamplefreqoi = timwin(ifreqoi) .* fsample (exact floating point)
            # We need to recompute from window_lengths_sec to get exact values, not rounded integers
            if is_fixed
                # Fixed window: same for all frequencies
                window_lengths_samples_exact = fill(Float64(window_length * dat.sample_rate), num_frex)
            else
                # Adaptive window: cycles / frequency
                window_lengths_samples_exact = [(cycles / freqs[fi]) * dat.sample_rate for fi = 1:num_frex]
            end
            # Use unpadded data length for edge filtering (FieldTrip uses ndatsample = size(dat, 2), unpadded)
            # Padding is only for FFT efficiency, but edge filtering should be based on actual data length
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
        is_fixed ? :taper_fixed : :taper_adaptive,
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
- `filter_edges::Bool=true`: If `true` (default), filters out edge regions where the window extends beyond the data
  (FieldTrip's "cone of influence" - marks edges as NaN). If `false`, uses all convolution results including edges.

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
    filter_edges::Bool = true,
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

