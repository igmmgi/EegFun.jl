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
- `interval_selection::Interval=times()`: Interval selection predicate. See `times()` for options.
  - Example: `interval_selection=times((-0.5, 2.0))` for time interval from -0.5 to 2.0 seconds
  - Example: `interval_selection=times()` for all time points (default)
  - Default: all time points
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
    return_phase::Bool = false,
    filter_edges::Bool = true,
    use_gpu::Bool = false,
)

    # Validate padding parameter
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end

    gpu_active = use_gpu && is_gpu_available()
    if use_gpu && !gpu_active
        @warn "GPU requested but CUDA is not available or functional. Falling back to CPU."
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

    # Get time points and determine output time grid
    times = time_vector(dat)
    time_min = minimum(times)
    time_max = maximum(times)
    time_steps_range = time_min:time_steps:time_max
    time_indices, times_out = find_times(times, time_steps_range)
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
    n_pre_pad = (!isnothing(pad) && (pad == :both || pad == :pre)) ? n_samples_per_epoch - 1 : 0
    n_post_pad = (!isnothing(pad) && (pad == :both || pad == :post)) ? n_samples_per_epoch - 1 : 0
    n_padded_samples = n_pre_pad + n_samples_per_epoch + n_post_pad

    max_window_samples = maximum(n_window_samples_per_freq)
    n_samples_padded = max(n_padded_samples, max_window_samples)

    # Pre-compute complex wavelets and their FFTs for each frequency using a pre-planned FFT
    wavelet_padded = zeros(ComplexF64, n_samples_padded)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)
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

        fill!(wavelet_padded, 0)
        wavelet_padded[1:n_window_samples] = wavelet
        w_fft = zeros(ComplexF64, n_samples_padded)
        mul!(w_fft, p_fft_wavelet, wavelet_padded)
        wavelet_ffts[fi] = w_fft
    end

    if gpu_active
        template_padded_batch_gpu = gpu_array(zeros(ComplexF64, n_samples_padded, n_trials))
        fft_plan_padded_batch = plan_fft!(template_padded_batch_gpu, 1)
        ifft_plan_padded_batch = plan_ifft!(template_padded_batch_gpu, 1)

        gpu_buffers = (
            local_data_padded_gpu = gpu_array(zeros(Float64, n_samples_padded, n_trials)),
            local_data_fft_gpu = gpu_array(zeros(ComplexF64, n_samples_padded, n_trials)),
            local_conv_result_gpu = gpu_array(zeros(ComplexF64, n_samples_padded, n_trials)),
            curr_wavelet_gpu = gpu_array(zeros(ComplexF64, n_samples_padded)),
            eegpower_trials_gpu = return_trials ? gpu_array(zeros(Float64, n_trials, num_frex, n_times)) : nothing,
            eegconv_trials_gpu = (return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_trials, num_frex, n_times)) : nothing,
            eegpower_avg_gpu = !return_trials ? gpu_array(zeros(Float64, num_frex, n_times)) : nothing,
            eegconv_avg_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, num_frex, n_times)) : nothing,
            adjusted_indices_gpu = gpu_array(zeros(Int, n_times)),
            extracted_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_times, n_trials)) : nothing,
            sum_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_times, 1)) : nothing,
            sum_complex_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_times, 1)) : nothing,
        )
    else
        # Plan for batch FFT of entire padded data (all trials at once)
        template_padded_batch = zeros(ComplexF64, n_samples_padded, n_trials)
        fft_plan_padded_batch = plan_fft!(template_padded_batch, 1, flags = FFTW.MEASURE)
        ifft_plan_padded_batch = plan_ifft!(template_padded_batch, 1, flags = FFTW.MEASURE)
        gpu_buffers = nothing
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

    # A lock to safely write into the DataFrames
    df_lock = ReentrantLock()

    function _stft_cpu_inner!(
        fi::Int,
        n_trials::Int,
        n_samples_padded::Int,
        n_pre_pad::Int,
        half_window::Int,
        norm_factor::Float64,
        n_times::Int,
        time_indices::Vector{Int},
        return_trials::Bool,
        return_phase::Bool,
        local_data_fft::Matrix{ComplexF64},
        curr_wavelet::Vector{ComplexF64},
        ifft_plan_padded_batch::Any,
        local_conv_result::Matrix{ComplexF64},
        eegpower_trials::Array{Float64,3},
        eegconv_trials::Array{ComplexF64,3},
        eegpower_avg::Matrix{Float64},
        eegconv_avg::Matrix{ComplexF64},
    )
        @inbounds for trial_idx = 1:n_trials
            @simd for i = 1:n_samples_padded
                local_conv_result[i, trial_idx] = local_data_fft[i, trial_idx] * curr_wavelet[i]
            end
        end

        # IFFT to get time-domain result (in-place)
        mul!(local_conv_result, ifft_plan_padded_batch, local_conv_result)

        @inbounds @simd for i in eachindex(local_conv_result)
            local_conv_result[i] *= norm_factor
        end

        # Extract requested time points
        @inbounds for ti_idx = 1:n_times
            sample_idx = time_indices[ti_idx]
            # Shift by n_pre_pad because local_conv_result contains virtually padded signal
            adjusted_idx = sample_idx + n_pre_pad + half_window
            # Clamp to valid range
            if adjusted_idx < 1
                adjusted_idx = 1
            elseif adjusted_idx > n_samples_padded
                adjusted_idx = n_samples_padded
            end

            if return_trials
                @inbounds @simd for trial_idx = 1:n_trials
                    val = local_conv_result[adjusted_idx, trial_idx]
                    eegpower_trials[trial_idx, fi, ti_idx] = abs2(val)
                    if return_phase
                        eegconv_trials[trial_idx, fi, ti_idx] = val
                    end
                end
            else
                power_sum::Float64 = 0.0
                conv_sum::ComplexF64 = 0.0im
                @inbounds @simd for trial_idx = 1:n_trials
                    val = local_conv_result[adjusted_idx, trial_idx]
                    power_sum += abs2(val)
                    if return_phase
                        conv_sum += val
                    end
                end
                eegpower_avg[fi, ti_idx] = power_sum
                if return_phase
                    eegconv_avg[fi, ti_idx] = conv_sum
                end
            end
        end
    end

    function _process_stft_channel!(
        channel::Symbol,
        dat::EpochData,
        n_trials::Int,
        n_samples_per_epoch::Int,
        n_padded_samples::Int,
        n_samples_padded::Int,
        num_frex::Int,
        return_trials::Bool,
        return_phase::Bool,
        filter_edges::Bool,
        pad::Union{Nothing,Symbol},
        n_times::Int,
        time_indices::Vector{Int},
        n_window_samples_per_freq::Vector{Int},
        fft_plan_padded_batch,
        ifft_plan_padded_batch,
        wavelet_ffts::Vector{Vector{ComplexF64}},
        df_lock::ReentrantLock,
        power_df,
        phase_df,
        gpu_buffers,
    )
        # Thread-local reusable buffers (CPU)
        local_data_padded = Matrix{Float64}(undef, n_samples_padded, n_trials)
        local_data_fft = Matrix{ComplexF64}(undef, n_samples_padded, n_trials)
        local_conv_result = Matrix{ComplexF64}(undef, n_samples_padded, n_trials)
        local_trial_signals = Matrix{Float64}(undef, n_samples_per_epoch, n_trials)

        gpu_active = !isnothing(gpu_buffers)
        if gpu_active
            local_data_padded_gpu = gpu_buffers.local_data_padded_gpu
            local_data_fft_gpu = gpu_buffers.local_data_fft_gpu
            local_conv_result_gpu = gpu_buffers.local_conv_result_gpu
            eegpower_trials_gpu = gpu_buffers.eegpower_trials_gpu
            eegconv_trials_gpu = gpu_buffers.eegconv_trials_gpu
            eegpower_avg_gpu = gpu_buffers.eegpower_avg_gpu
            eegconv_avg_gpu = gpu_buffers.eegconv_avg_gpu

            adjusted_indices_gpu = gpu_buffers.adjusted_indices_gpu
            extracted_power_gpu = gpu_buffers.extracted_power_gpu
            sum_power_gpu = gpu_buffers.sum_power_gpu
            sum_complex_gpu = gpu_buffers.sum_complex_gpu
            adjusted_indices_cpu = Vector{Int}(undef, n_times)

            # Reset output arrays
            if return_trials
                fill!(eegpower_trials_gpu, 0.0)
                if return_phase
                    fill!(eegconv_trials_gpu, 0.0im)
                end
            else
                fill!(eegpower_avg_gpu, 0.0)
                if return_phase
                    fill!(eegconv_avg_gpu, 0.0im)
                end
            end
        else
            # Pre-allocate reusable output buffers (reused across trials for this channel)
            local eegpower_trials::Array{Float64,3} =
                return_trials ? zeros(Float64, n_trials, num_frex, n_times) : Array{Float64,3}(undef, 0, 0, 0)
            local eegconv_trials::Array{ComplexF64,3} =
                (return_trials && return_phase) ? zeros(ComplexF64, n_trials, num_frex, n_times) : Array{ComplexF64,3}(undef, 0, 0, 0)
            local eegpower_avg::Matrix{Float64} = return_trials ? Matrix{Float64}(undef, 0, 0) : zeros(Float64, num_frex, n_times)
            local eegconv_avg::Matrix{ComplexF64} =
                (!return_trials && return_phase) ? zeros(ComplexF64, num_frex, n_times) : Matrix{ComplexF64}(undef, 0, 0)
        end

        # Pre-extract all trial data for this channel into a matrix for batch processing
        for trial_idx = 1:n_trials
            col = dat.data[trial_idx][!, channel]
            if col isa Vector{Float64}
                col_f64 = col::Vector{Float64}
                @inbounds @simd for i = 1:n_samples_per_epoch
                    local_trial_signals[i, trial_idx] = col_f64[i]
                end
            else
                @inbounds @simd for i = 1:n_samples_per_epoch
                    local_trial_signals[i, trial_idx] = Float64(col[i])
                end
            end
        end

        # Clear/initialize output buffers for this channel
        if !gpu_active
            if return_trials
                fill!(eegpower_trials, 0.0)
                return_phase && fill!(eegconv_trials, 0.0im)
            else
                fill!(eegpower_avg, 0.0)
                return_phase && fill!(eegconv_avg, 0.0im)
            end
        end

        # Pad data to n_samples_padded (zero-padding at the end)
        fill!(local_data_padded, 0.0)

        # Explicit loop for padding to avoid copy allocations (Virtual Padding)
        @inbounds for trial_idx = 1:n_trials
            for j = 1:n_padded_samples
                local src_j
                if j <= n_pre_pad
                    src_j = n_samples_per_epoch - j + 1
                elseif j > n_pre_pad + n_samples_per_epoch
                    idx_in_post = j - (n_pre_pad + n_samples_per_epoch)
                    src_j = n_samples_per_epoch - idx_in_post
                else
                    src_j = j - n_pre_pad
                end

                local_data_padded[j, trial_idx] = local_trial_signals[src_j, trial_idx]
            end
        end

        # FFT entire padded data (batch process all trials at once)
        @inbounds @simd for i in eachindex(local_data_padded)
            local_data_fft[i] = ComplexF64(local_data_padded[i])
        end
        if gpu_active
            copyto!(local_data_fft_gpu, local_data_fft)
            mul!(local_data_fft_gpu, fft_plan_padded_batch, local_data_fft_gpu)
        else
            mul!(local_data_fft, fft_plan_padded_batch, local_data_fft)
        end

        # Process each frequency
        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]

            # Frequency-domain convolution: multiply data FFT by wavelet FFT
            # Explicit loops to avoid broadcast allocations on large matrices
            curr_wavelet = wavelet_ffts[fi]
            norm_factor = sqrt(2.0 / n_window_samples)
            half_window = n_window_samples ÷ 2

            if gpu_active
                # GPU inner loop broadcast
                copyto!(gpu_buffers.curr_wavelet_gpu, curr_wavelet)
                local_conv_result_gpu .= local_data_fft_gpu .* gpu_buffers.curr_wavelet_gpu
                mul!(local_conv_result_gpu, ifft_plan_padded_batch, local_conv_result_gpu)
                local_conv_result_gpu .*= norm_factor

                # Precompute adjusted indices
                for ti_idx = 1:n_times
                    adjusted_indices_cpu[ti_idx] = min(max(time_indices[ti_idx] + n_pre_pad + half_window, 1), n_samples_padded)
                end
                copyto!(adjusted_indices_gpu, adjusted_indices_cpu)

                extracted_gpu = @view local_conv_result_gpu[adjusted_indices_gpu, :]

                if return_trials
                    @views eegpower_trials_gpu[:, fi, :] .= transpose(abs2.(extracted_gpu))
                    if return_phase
                        @views eegconv_trials_gpu[:, fi, :] .= transpose(extracted_gpu)
                    end
                else
                    extracted_power_gpu .= abs2.(extracted_gpu)
                    sum!(sum_power_gpu, extracted_power_gpu)
                    @views eegpower_avg_gpu[fi, :] .= vec(sum_power_gpu)

                    if return_phase
                        sum!(sum_complex_gpu, extracted_gpu)
                        @views eegconv_avg_gpu[fi, :] .= vec(sum_complex_gpu)
                    end
                end
            else
                _stft_cpu_inner!(
                    fi,
                    n_trials,
                    n_samples_padded,
                    n_pre_pad,
                    half_window,
                    norm_factor,
                    n_times,
                    time_indices,
                    return_trials,
                    return_phase,
                    local_data_fft,
                    curr_wavelet,
                    ifft_plan_padded_batch,
                    local_conv_result,
                    eegpower_trials,
                    eegconv_trials,
                    eegpower_avg,
                    eegconv_avg,
                )
            end
        end

        if gpu_active
            if !return_trials
                eegpower_avg_gpu ./= n_trials
                if return_phase
                    eegconv_avg_gpu ./= n_trials
                end
            end

            if return_trials
                eegpower_trials = Array(eegpower_trials_gpu)
                if return_phase
                    eegconv_trials = Array(eegconv_trials_gpu)
                end
            else
                eegpower_avg = Array(eegpower_avg_gpu)
                if return_phase
                    eegconv_avg = Array(eegconv_avg_gpu)
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
            adjusted_time_indices = time_indices .+ n_pre_pad
            _filter_edges!(
                return_trials ? eegpower_trials : eegpower_avg,
                return_trials ? (return_phase ? eegconv_trials : nothing) : (return_phase ? eegconv_avg : nothing),
                num_frex,
                adjusted_time_indices,
                window_lengths_samples_exact,
                n_padded_samples,
            )
        end

        # Normalise outside the lock to minimise lock hold time
        if !gpu_active
            if !return_trials
                eegpower_avg ./= n_trials
                if return_phase
                    eegconv_avg ./= n_trials
                end
            end
        end

        lock(df_lock) do
            if return_trials # Store each trial separately
                for trial_idx = 1:n_trials
                    power_df[trial_idx][!, channel] = copy(vec(@view eegpower_trials[trial_idx, :, :]))
                    if return_phase
                        phase_df[trial_idx][!, channel] = copy(vec(angle.(@view eegconv_trials[trial_idx, :, :])))
                    else
                        phase_df[trial_idx][!, channel] = fill(NaN, num_frex * n_times)
                    end
                end
            else
                power_df[!, channel] = copy(vec(eegpower_avg))
                if return_phase
                    phase_vec = Vector{Float64}(undef, num_frex * n_times)
                    @inbounds @simd for i in eachindex(eegconv_avg)
                        phase_vec[i] = angle(eegconv_avg[i])
                    end
                    phase_df[!, channel] = phase_vec
                else
                    phase_df[!, channel] = fill(NaN, num_frex * n_times)
                end
            end
        end
    end

    if gpu_active
        for channel in selected_channels
            _process_stft_channel!(
                channel,
                dat,
                n_trials,
                n_samples_per_epoch,
                n_padded_samples,
                n_samples_padded,
                num_frex,
                return_trials,
                return_phase,
                filter_edges,
                pad,
                n_times,
                time_indices,
                n_window_samples_per_freq,
                fft_plan_padded_batch,
                ifft_plan_padded_batch,
                wavelet_ffts,
                df_lock,
                power_df,
                phase_df,
                gpu_buffers,
            )
        end
    else
        Threads.@threads for channel in selected_channels
            _process_stft_channel!(
                channel,
                dat,
                n_trials,
                n_samples_per_epoch,
                n_padded_samples,
                n_samples_padded,
                num_frex,
                return_trials,
                return_phase,
                filter_edges,
                pad,
                n_times,
                time_indices,
                n_window_samples_per_freq,
                fft_plan_padded_batch,
                ifft_plan_padded_batch,
                wavelet_ffts,
                df_lock,
                power_df,
                phase_df,
                gpu_buffers,
            )
        end
    end

    # Restore original column order which might be shuffled by multithreading
    expected_cols = vcat([:time, :freq], selected_channels)
    if return_trials
        for i = 1:n_trials
            select!(power_df[i], expected_cols)
            select!(phase_df[i], expected_cols)
        end
    else
        select!(power_df, expected_cols)
        select!(phase_df, expected_cols)
    end

    # Create and return appropriate data type
    return_type = return_trials ? TimeFreqEpochData : TimeFreqData
    return return_type(
        dat.file,
        dat.condition,
        dat.condition_name,
        power_df,
        phase_df,
        copy(dat.layout),
        dat.sample_rate,
        is_fixed ? :taper_fixed : :taper_adaptive,
        nothing,  # baseline
        copy(dat.analysis_info),
    )
end

# Vector version: process each EpochData element
function tf_stft(data_vec::Vector{<:EpochData}; kwargs...)
    return [tf_stft(dat; kwargs...) for dat in data_vec]
end


# === BATCH PROCESSING FUNCTIONS ===
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

        output_dir = something(output_dir, joinpath(input_dir, "tf_stft_$(clean_pattern(file_pattern))"))
        mkpath(output_dir)

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

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "TF STFT")

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
