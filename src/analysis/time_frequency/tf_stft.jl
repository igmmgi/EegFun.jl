"""
    tf_stft(dat::EpochData; 
            channel_selection::Function=channels(),
            interval_selection::Interval=samples(),
            frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40),
            window_length::Union{Nothing,Real}=nothing,
            cycles::Union{Nothing,Real}=nothing,
            time_steps::Real=0.05,
            pad::Union{Nothing,Symbol}=nothing,
            return_trials::Bool=false,
            filter_edges::Bool=true)

Short-Time Fourier Transform (STFT) time-frequency analysis using a sliding Hanning window 

Supports both fixed-length windows (consistent time resolution) and adaptive windows (constant cycles per frequency).

# Arguments
- `dat::EpochData`: Epoched EEG data

# Keyword Arguments
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `interval_selection::Interval=samples()`: Sample selection predicate. See `samples()` for options.
  - Example: `sample_selection=samples((-0.5, 2.0))` for time interval from -0.5 to 2.0 seconds
  - Example: `sample_selection=samples()` for all time points (default)
  - Default: all samples
- `frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40)`: Frequency specification.
  - Can be any range or vector of frequencies in Hz
  - For linear spacing: `frequencies=1:1:40` or `frequencies=range(1, 40, length=40)`
  - For logarithmic spacing: `frequencies=logrange(1, 40, length=30)`
  - Default: `range(1, 40, length=40)` (40 linearly-spaced frequencies from 1 to 40 Hz)
- `window_length::Union{Nothing,Real}=nothing`: Fixed window length in seconds (same for all frequencies).
  - Example: `window_length=0.3` uses a fixed 0.3 second window for all frequencies
  - **Exactly one of `window_length` or `cycles` must be specified.**
- `cycles::Union{Nothing,Real}=nothing`: Number of cycles per frequency (adaptive window).
  - Example: `cycles=7` uses 7 cycles per frequency (window length = 7/frequency)
  - **Exactly one of `window_length` or `cycles` must be specified.**
- `time_steps::Real=0.05`: Step size for extracting time points in seconds.
  - Creates time points from the selected time range with the specified step size
  - Default: `0.05` (50 ms)
  - Example: `time_steps=0.01` creates time points every 0.01 seconds within the selected time interval
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default). Time points where intervals extend beyond data boundaries are excluded with a warning.
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.
- `filter_edges::Bool=true`: If `true` (default), filters out edge regions where the interval extends beyond the data

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Examples
```julia
# Default: linear frequencies 1-40 Hz, 40 points, fixed window
tf_data = tf_stft(epochs; window_length=0.3)

# Log-spaced frequencies with fixed window
tf_data = tf_stft(epochs; frequencies=logrange(2, 80, length=30), window_length=0.3)

# Adaptive window: 7 cycles per frequency 
tf_data = tf_stft(epochs; frequencies=2:1:30, cycles=7, sample_selection=samples((-0.5, 1.5)), time_steps=0.05)
```
"""
function tf_stft(
    dat::EpochData;
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    frequencies::Union{AbstractRange,AbstractVector{<:Real}} = range(1, 40, length = 40),
    time_steps::Real = 0.05,
    window_length::Union{Nothing,Real} = nothing,
    cycles::Union{Nothing,Real} = nothing,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
    filter_edges::Bool = true,
)

    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end

    # Subset data with channel and interval selection
    dat = subset(dat; channel_selection = channel_selection, interval_selection = interval_selection)
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

    # Get original data time range (before padding) - these are the time points we want in output
    n_samples_original_unpadded = n_samples(dat)  # Store original unpadded length for edge filtering
    times_original = time_vector(dat)

    # Apply padding if requested (mutating dat directly for consistency with tf_morlet)
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get sample rate and time vector from processed data
    times_processed = time_vector(dat)

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

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples(dat)

    # Use frequencies directly (convert to vector if needed for indexing)
    freqs = collect(frequencies)

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
        # Adaptive window mode: window length = cycles / frequency 
        cycles_vec = fill(Float64(cycles), num_frex)
        window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds 
        n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
        if any(n -> n < 2, n_window_samples_per_freq)
            min_samples = minimum(n_window_samples_per_freq)
            error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
        end
        is_fixed = false
    end

    n_times = length(times_out)

    # Pre-compute Hanning windows for each frequency
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

        angle_in =
            range(-(n_window_samples - 1) / 2, (n_window_samples - 1) / 2, length = n_window_samples) .* (2π * freq / dat.sample_rate)
        cos_wav = hanning_window .* cos.(angle_in)
        sin_wav = hanning_window .* sin.(angle_in)
        wavelet = cos_wav .+ im .* sin_wav

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
        fill!(eegpower, 0.0)
        fill!(eegconv, 0.0im)

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

            norm_factor = sqrt(2.0 / n_window_samples)
            @inbounds @simd for i in eachindex(conv_result)
                conv_result[i] *= norm_factor
            end

            # Extract requested time points (time_indices are sample indices in processed data)
            half_window = n_window_samples ÷ 2
            @inbounds for ti_idx = 1:n_times
                sample_idx = time_indices[ti_idx]
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

        if filter_edges
            # Compute exact window lengths in samples (floating point) for edge filtering
            if is_fixed
                # Fixed window: same for all frequencies
                window_lengths_samples_exact = fill(Float64(window_length * dat.sample_rate), num_frex)
            else
                # Adaptive window: cycles / frequency
                window_lengths_samples_exact = [(cycles / freqs[fi]) * dat.sample_rate for fi = 1:num_frex]
            end
            # Padding is only for FFT efficiency, but edge filtering should be based on actual data length
            _filter_edges!(eegpower, eegconv, num_frex, time_indices, window_lengths_samples_exact, n_samples_per_epoch)
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = vec(@view eegpower[:, :, trial_idx])
                phase_df[trial_idx][!, channel] = vec(angle.(@view eegconv[:, :, trial_idx]))
            end
        else
            eegpower ./= n_trials
            eegconv ./= n_trials
            power_df[!, channel] = copy(vec(eegpower))
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

# Vector version: process each EpochData element
function tf_stft(data_vec::Vector{<:EpochData}; kwargs...)
    return [tf_stft(dat; kwargs...) for dat in data_vec]
end


# =============================================================================
#     BATCH PROCESSING FUNCTIONS
# =============================================================================
function tf_stft(
    file_pattern::String;
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    kwargs...,
)
    log_file = "tf_stft.log"
    setup_global_logging(log_file)

    try
        @info "Batch tf_stft started at $(now())"
        @log_call "tf_stft"

        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_stft_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files for tf_stft analysis"

        process_fn =
            (input_path, output_path) -> begin
                filename = basename(input_path)
                data = read_data(input_path)
                if isnothing(data) || !(data isa Vector{<:EpochData})
                    return BatchResult(false, filename, "Invalid data type")
                end
                data = _condition_select(data, condition_selection)
                tf_results = tf_stft(data; kwargs...)
                jldsave(output_path; data = tf_results)
                return BatchResult(true, filename, "TF STFT analysis complete ($(length(tf_results)) conditions)")
            end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "TF STFT")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
