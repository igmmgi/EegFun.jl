"""
    tf_morlet(dat::EpochData; 
              lin_freqs::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
              log_freqs::Union{Nothing,Tuple{Real,Real,Int}}=nothing,
              cycles::Union{Real,Tuple{Real,Real}}=(3, 10),
              time_steps::Union{Nothing,Tuple{Real,Real,Real}}=nothing,
              channel_selection::Function=channels())

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

# Returns
- `TimeFreqData`: Time-frequency data object containing power (and phase, currently empty) for selected channels

# Example
```julia
# Log-spaced frequencies (30 frequencies from 2 to 80 Hz), single channel
tf_data = tf_morlet(epochs; log_freqs=(2, 80, 30), channel_selection=channels(:Cz))

# Linear-spaced frequencies (2 Hz steps from 2 to 80 Hz), specific time resolution, multiple channels
tf_data = tf_morlet(epochs; lin_freqs=(2, 80, 2), time_steps=(-0.5, 2.0, 0.01), channel_selection=channels([:Cz, :Pz]))
```
"""
function tf_morlet(
    dat::EpochData;
    channel_selection::Function = channels(),
    time_steps::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    lin_freqs::Union{Nothing,Tuple{Real,Real,Real}} = nothing,
    log_freqs::Union{Nothing,Tuple{Real,Real,Int}} = nothing,
    cycles::Union{Real,Tuple{Real,Real}} = (3, 10),
)

    # Get selected channels using channel selection predicate
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(dat))")

    # Validate frequency specification - exactly one must be provided
    isnothing(lin_freqs) && isnothing(log_freqs) && error("Either `lin_freqs` or `log_freqs` must be specified")
    !isnothing(lin_freqs) && !isnothing(log_freqs) && error("Only one of `lin_freqs` or `log_freqs` can be specified, not both")

    # Get sample rate
    sr = Float64(dat.sample_rate)

    # Get time vector from first epoch (all epochs should have same time points)
    times_original = dat.data[1][!, :time]
    n_samples_original = length(times_original)

    # Handle time_steps parameter - determine which time points to extract from results
    # time_steps is in SECONDS: (start, stop, step) where all values are in seconds
    if isnothing(time_steps)
        # Use all original time points (all in seconds)
        time_indices = 1:n_samples_original
        times_out = times_original  # times_original is in seconds
    else
        # Generate time points from tuple (start, stop, step) - all in SECONDS
        start_time, stop_time, step_time = time_steps
        # Generate range of time values in seconds
        time_steps_vec = collect(Float64, range(Float64(start_time), Float64(stop_time), step = Float64(step_time)))
        # If the last value is significantly less than stop_time, add it
        if !isempty(time_steps_vec) && (stop_time - time_steps_vec[end]) > step_time / 2
            push!(time_steps_vec, Float64(stop_time))
        end
        
        # Validate that requested times (in seconds) are within the original data range (in seconds)
        time_min = minimum(times_original)  # seconds
        time_max = maximum(times_original)  # seconds
        if start_time < time_min || stop_time > time_max
            @minimal_warning "Requested time range ($start_time to $stop_time seconds) extends beyond data range ($time_min to $time_max seconds). Clipping to available range."
        end
        
        # Find nearest time points using find_times utility function
        time_indices, times_out = find_times(times_original, time_steps_vec)
        
        if isempty(time_indices)
            error("No valid time points found in requested range ($start_time to $stop_time seconds)")
        end
    end

    # Get number of trials/epochs
    n_trials = length(dat.data)
    n_samples = n_samples_original  # Use full signal for convolution

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
    n_data = n_samples * n_trials
    n_convolution = n_wavelet + n_data - 1
    n_conv_pow2 = nextpow(2, n_convolution)

    # Initialize output DataFrames (will accumulate data from all channels)
    n_times = length(times_out)
    power_df = DataFrame()
    power_df.time = repeat(times_out, inner = num_frex)  # Each time repeated for all freqs
    power_df.freq = repeat(freqs, outer = n_times)        # All freqs repeated for each time
    
    phase_df = DataFrame()
    phase_df.time = copy(power_df.time)
    phase_df.freq = copy(power_df.freq)

    # Process each selected channel
    for channel in selected_channels
        # Extract full signal data for convolution (use all time points)
        # MATLAB: reshape(EEG.data(chan,:,:), 1, EEG.pnts*EEG.trials)
        signal = zeros(Float64, n_samples_original * n_trials)
        for (trial_idx, epoch_df) in enumerate(dat.data)
            signal[(trial_idx-1)*n_samples_original+1:trial_idx*n_samples_original] = epoch_df[!, channel]
        end

        # Get FFT of data (MATLAB: fft(reshape(...), n_conv_pow2))
        signal_padded = zeros(ComplexF64, n_conv_pow2)
        signal_padded[1:n_data] = signal
        eegfft = fft(signal_padded)

        # Initialize output power matrix (frequencies × time) - use full time resolution
        eegpower_full = zeros(Float64, num_frex, n_samples)

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

            # Reshape to (n_samples × n_trials) and average power over trials
            # MATLAB: mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2)
            eegconv_reshaped = reshape(eegconv, n_samples, n_trials)
            eegpower_full[fi, :] = vec(mean(abs2.(eegconv_reshaped), dims = 2))
        end

        # Extract only requested time points from full power matrix
        eegpower = eegpower_full[:, time_indices]

        # Reshape power matrix to match expected order: [freq1_t1, freq2_t1, ..., freqN_t1, freq1_t2, ...]
        # eegpower is (num_frex × n_times), we need column-major order
        # vec(eegpower) gives: [freq1_time1, freq2_time1, ..., freqN_time1, freq1_time2, ...]
        power_values = vec(eegpower)
        power_df[!, channel] = power_values

        # Add phase column (NaN for now)
        phase_df[!, channel] = fill(NaN, nrow(power_df))
    end

    # Create TimeFreqData object
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