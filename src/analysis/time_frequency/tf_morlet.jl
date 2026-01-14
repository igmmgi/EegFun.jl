"""
    tf_morlet(dat::EpochData; 
              channel_selection::Function=channels(),
              sample_selection::Function=samples(),
              frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40),
              cycles::Union{Real,Tuple{Real,Real}}=7,
              pad::Union{Nothing,Symbol}=nothing,
              return_trials::Bool=false,
              filter_edges::Bool=true)

Time-frequency analysis using Morlet wavelets.

Performs continuous wavelet transform using complex Morlet wavelets. The function supports both 
linear and logarithmic frequency spacing, with optional padding to reduce edge artifacts.

# Arguments
- `dat::EpochData`: Epoched EEG data to analyze

# Keyword Arguments
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
  - Default: all channels
- `sample_selection::Function=samples()`: Sample selection predicate. See `samples()` for options.
  - Example: `sample_selection=samples((-0.5, 2.0))` for time window from -0.5 to 2.0 seconds
  - Example: `sample_selection=samples()` for all time points (default)
  - Default: all samples
- `frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40)`: Frequency specification.
  - Can be any range or vector of frequencies in Hz
  - For linear spacing: `frequencies=1:1:40` or `frequencies=range(1, 40, length=40)`
  - For logarithmic spacing: `frequencies=logrange(1, 40, length=30)`
  - Default: `range(1, 40, length=40)` (40 linearly-spaced frequencies from 1 to 40 Hz)
- `cycles::Union{Real,Tuple{Real,Real}}=7`: Number of cycles in the wavelet. Controls time-frequency trade-off:
  - Single number: fixed cycles for all frequencies (e.g., `cycles=5`)
  - Tuple `(min, max)`: log-spaced cycles from `min` to `max` across frequencies
  - More cycles = better frequency resolution, worse time resolution
  - Default: `7` cycles
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default)
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended for best edge artifact reduction)
  - Padding extends the data for convolution, then results are automatically unpadded
- `return_trials::Bool=false`: Whether to preserve individual trials:
  - `false` (default): Returns `TimeFreqData` with trials averaged
  - `true`: Returns `TimeFreqEpochData` with individual trials preserved
- `filter_edges::Bool=true`: Whether to filter edge regions where the wavelet extends beyond the data:
  - `true` (default): Sets edge regions to `NaN` where the wavelet window extends beyond data boundaries
  - `false`: Keeps all computed values (may include edge artifacts)
  - Edge filtering accounts for padding if applied

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Examples
```julia
# Default: linear frequencies 1-40 Hz, 40 points, averaged across trials
tf_data = tf_morlet(epochs)

# Log-spaced frequencies (30 frequencies from 1 to 40 Hz)
tf_data = tf_morlet(epochs; frequencies=logrange(1, 40, length=30))

# Log-spaced frequencies, single channel
tf_data = tf_morlet(epochs; channel_selection=channels(:Cz), frequencies=logrange(2, 80, length=30))

# Linear frequencies with padding and individual trials preserved
tf_epochs = tf_morlet(epochs; sample_selection=samples((-0.5, 2.0)), frequencies=2:2:80, pad=:both, return_trials=true)
```
"""
function tf_morlet(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    frequencies::Union{AbstractRange,AbstractVector{<:Real}} = range(1, 40, length = 40),
    cycles::Union{Real,Tuple{Real,Real}} = 7,
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

    # Get original data time range (before padding) - these are the time points we want in output
    n_samples_original_unpadded = n_samples(dat)  # Store original unpadded length for edge filtering

    # Apply padding if requested 
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples(dat) # Use full signal for convolution (may be padded)

    # Get all time points from processed (possibly padded) data
    times_all = time(dat)

    # Use frequency input directly (ranges and vectors both work)
    num_frex = length(frequencies)

    # Validate frequencies
    if num_frex == 0
        error("`frequencies` must contain at least one frequency")
    end
    if any(f -> f <= 0, frequencies)
        error("All frequencies in `frequencies` must be positive")
    end

    # Define cycles/sigma
    if cycles isa Tuple
        cycles_vec = logrange(cycles[1], cycles[2], length = num_frex)
    else
        cycles_vec = fill(cycles, num_frex)
    end

    # Calculate which time indices to keep (original unpadded samples)
    if !isnothing(pad)
        if pad == :pre || pad == :both
            n_pre_pad = n_samples_original_unpadded - 1
        else  # :post
            n_pre_pad = 0
        end
        time_indices_out = (n_pre_pad+1):(n_pre_pad+n_samples_original_unpadded)
    else
        time_indices_out = 1:n_samples_original_unpadded
    end
    n_times_out = length(time_indices_out)

    # Initialize output structures with only original time points as we unpad!
    times_out = times_all[time_indices_out]
    time_col = repeat(times_out, inner = num_frex)
    freq_col = repeat(frequencies, outer = n_times_out)
    if return_trials
        power_df = [DataFrame(time = time_col, freq = freq_col, copycols = false) for _ = 1:n_trials]
        phase_df = [DataFrame(time = time_col, freq = freq_col, copycols = false) for _ = 1:n_trials]
    else
        power_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
        phase_df = DataFrame(time = time_col, freq = freq_col, copycols = false)
    end

    # Pre-compute convolution length for single trial processing
    max_sigma = maximum(cycles_vec ./ (2 * pi .* frequencies))
    max_hw = ceil(Int, 6 * max_sigma * dat.sample_rate) ÷ 2
    max_wl = max_hw * 2 + 1
    n_conv_pow2 = nextpow(2, max_wl + n_samples_per_epoch - 1)

    # Pre-allocate reusable buffers for single trial processing
    signal_padded = zeros(ComplexF64, n_conv_pow2)
    eegfft = zeros(ComplexF64, n_conv_pow2)
    wavelet_padded = zeros(ComplexF64, n_conv_pow2)
    conv_buffer = zeros(ComplexF64, n_conv_pow2)
    eegconv_buffer = zeros(ComplexF64, n_conv_pow2)

    # Create FFT plans using the actual buffers 
    p_fft_signal = plan_fft(signal_padded, flags = FFTW.MEASURE)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)
    p_ifft_conv = plan_ifft(eegconv_buffer, flags = FFTW.MEASURE)

    # Pre-compute some constants 
    inv_sr = 1.0 / dat.sample_rate
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    wl_per_freq = Vector{Int}(undef, num_frex)
    conv_indices_per_freq = Vector{Vector{Int}}(undef, num_frex)

    for fi = 1:num_frex
        freq_val = frequencies[fi]
        sigma = cycles_vec[fi] / (two_pi * freq_val)
        hw = ceil(Int, 6 * sigma * dat.sample_rate) ÷ 2
        wl = hw * 2 + 1
        wl_per_freq[fi] = wl
        valid_start = hw + 1

        # Pre-compute convolution indices for all samples in padded data
        conv_indices_per_freq[fi] = [valid_start + sample_idx - 1 for sample_idx = 1:n_samples_per_epoch]

        # Create wavelet directly in padded buffer
        fill!(wavelet_padded, 0)
        A = sqrt(1 / (sigma * sqrt_pi))
        two_pi_freq = two_pi * freq_val
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
    # Only allocate for original unpadded samples
    if return_trials
        eegpower = zeros(Float64, num_frex, n_times_out, n_trials)
        eegconv = zeros(ComplexF64, num_frex, n_times_out, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times_out)
        eegconv = zeros(ComplexF64, num_frex, n_times_out)
    end

    # TODO: initial testing shows this is just as fast as concatenating 
    # all trials and/or channels into a single matrix and then processing that
    # but it still seems too slow!!!

    # Process each selected channel
    for channel in channel_labels(dat)

        # Reset output arrays for this channel
        fill!(eegpower, 0)
        fill!(eegconv, 0)

        # ==============================================================
        # BATCH PROCESSING: Collect all trials into matrix
        # ==============================================================
        # Allocate matrix for all trials [n_samples × n_trials]
        signal_matrix = zeros(Float64, n_samples_per_epoch, n_trials)
        @inbounds for trial_idx = 1:n_trials
            signal_matrix[:, trial_idx] .= dat.data[trial_idx][!, channel]
        end

        # ==============================================================
        # FFT all trials at once (1 call instead of n_trials calls)
        # ==============================================================
        # Pad signal matrix for FFT
        signal_matrix_padded = zeros(ComplexF64, n_conv_pow2, n_trials)
        signal_matrix_padded[1:n_samples_per_epoch, :] .= signal_matrix

        # FFT along first dimension for all trials
        signal_fft_all = zeros(ComplexF64, n_conv_pow2, n_trials)
        for trial_idx = 1:n_trials
            mul!(@view(signal_fft_all[:, trial_idx]), p_fft_signal, @view(signal_matrix_padded[:, trial_idx]))
        end

        # ==============================================================
        # Process each frequency with all trials together
        # ==============================================================
        for fi = 1:num_frex
            wavelet_fft = wavelet_ffts[fi]
            conv_indices = conv_indices_per_freq[fi]

            # ==============================================================
            # Multiply wavelet FFT with all trial FFTs and IFFT
            # (40 IFFTs instead of 16,000)
            # ==============================================================
            conv_result_all = zeros(ComplexF64, n_conv_pow2, n_trials)

            for trial_idx = 1:n_trials
                # Element-wise multiply
                @inbounds @simd for i = 1:n_conv_pow2
                    conv_buffer[i] = wavelet_fft[i] * signal_fft_all[i, trial_idx]
                end

                # IFFT for this trial
                mul!(@view(conv_result_all[:, trial_idx]), p_ifft_conv, conv_buffer)
            end

            # ==============================================================
            # Extract results for all trials
            # ==============================================================
            @inbounds for trial_idx = 1:n_trials
                for (ti_out, ti_padded) in enumerate(time_indices_out)
                    val = conv_result_all[conv_indices[ti_padded], trial_idx]

                    if return_trials
                        eegpower[fi, ti_out, trial_idx] = abs2(val)
                        eegconv[fi, ti_out, trial_idx] = val
                    else
                        eegpower[fi, ti_out] += abs2(val)
                        eegconv[fi, ti_out] += val
                    end
                end
            end
        end

        if !return_trials
            eegpower ./= n_trials
            eegconv ./= n_trials
        end

        if filter_edges
            # Use padded indices for edge filtering - padding extends valid region
            _filter_edges!(eegpower, eegconv, num_frex, time_indices_out, wl_per_freq, n_samples_per_epoch)
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = vec(@view eegpower[:, :, trial_idx])
                phase_df[trial_idx][!, channel] = vec(angle.(@view eegconv[:, :, trial_idx]))
            end
        else
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
