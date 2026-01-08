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
