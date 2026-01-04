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


