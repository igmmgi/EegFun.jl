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
    channel_selection = channels(),
    interval_selection::Interval = times(),
    frequencies::AbstractVector{<:Real} = range(1, 40, length = 40),
    time_steps::Real = 0.05,
    window_length::Union{Nothing,Real} = nothing,
    cycles::Union{Nothing,Real} = nothing,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
    return_phase::Bool = false,
    filter_edges::Bool = true,
    use_gpu::Bool = false,
)

    # Validate arguments and subset data with channel and interval selection
    sdat, gpu_active = _tf_prepare(dat, pad, use_gpu, channel_selection, interval_selection, frequencies)

    # Get selected channels (after subsetting)
    selected_channels = channel_labels(sdat)

    # Use frequency input directly (ranges and vectors both work)
    num_frex = length(frequencies)

    # Validate window_length and cycles - exactly one must be provided
    if isnothing(window_length) && isnothing(cycles)
        @minimal_error("Either `window_length` (fixed window) or `cycles` (adaptive window) must be specified")
    end
    if !isnothing(window_length) && !isnothing(cycles)
        @minimal_error("Only one of `window_length` or `cycles` can be specified, not both")
    end
    if !isnothing(window_length) && window_length <= 0
        @minimal_error("`window_length` must be positive, got $window_length")
    end
    if !isnothing(cycles) && cycles <= 0
        @minimal_error("`cycles` must be positive, got $cycles")
    end

    # Get time points and determine output time grid
    time_indices, times_out = _tf_time_grid(sdat, time_steps)

    # Get number of trials/epochs
    n_trials = n_epochs(sdat)
    n_samples_per_epoch = n_samples(sdat)

    # Use frequencies directly (convert to vector if needed for indexing)
    freqs = collect(frequencies)

    # Determine window mode: fixed or adaptive
    # (single destructuring assignment so the channel closure below captures these without boxing)
    is_fixed, n_window_samples_per_freq, window_lengths_samples_exact = if !isnothing(window_length)
        # Fixed window mode: same window length for all frequencies
        n_window_samples_fixed = Int(round(window_length * sdat.sample_rate))
        _validate_tf_window_samples(n_window_samples_fixed)
        (true, fill(n_window_samples_fixed, num_frex), fill(Float64(window_length * sdat.sample_rate), num_frex))
    else
        # Adaptive window mode: window length = cycles / frequency 
        cycles_vec = fill(Float64(cycles), num_frex)
        window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds 
        n_window_samples_adaptive = [Int(round(wl * sdat.sample_rate)) for wl in window_lengths_sec]
        _validate_tf_window_samples(n_window_samples_adaptive)
        (false, n_window_samples_adaptive, [(cycles / freqs[fi]) * sdat.sample_rate for fi = 1:num_frex])
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
    n_pre_pad, _, n_padded_samples = _tf_padding(pad, n_samples_per_epoch)
    adjusted_time_indices = time_indices .+ n_pre_pad

    max_window_samples = maximum(n_window_samples_per_freq)
    n_samples_padded = max(n_padded_samples, max_window_samples)

    g = TFGeometry(
        n_trials,
        n_samples_per_epoch,
        n_pre_pad,
        n_padded_samples,
        n_samples_padded,
        num_frex,
        n_times,
        return_trials,
        return_phase,
        filter_edges,
    )

    # Pre-compute complex wavelets and their FFTs for each frequency using a pre-planned FFT
    wavelet_padded = zeros(ComplexF64, n_samples_padded)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)

    for fi = 1:num_frex
        n_window_samples = n_window_samples_per_freq[fi]
        freq = freqs[fi]
        hanning_window = hanning_windows[fi]

        angle_in =
            range(-(n_window_samples - 1) / 2, (n_window_samples - 1) / 2, length = n_window_samples) .* (2π * freq / sdat.sample_rate)
        cos_wav = hanning_window .* cos.(angle_in)
        sin_wav = hanning_window .* sin.(angle_in)
        wavelet = cos_wav .+ im .* sin_wav

        wavelet_ffts[fi] = _tf_padded_fft(wavelet, p_fft_wavelet, wavelet_padded)
    end

    fft_plan, ifft_plan = _tf_batch_fft_plans(n_samples_padded, n_trials, gpu_active)
    gpu_buffers =
        gpu_active ?
        _tf_gpu_buffers(
            g;
            extracted_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_times, n_trials)) : nothing,
            sum_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_times, 1)) : nothing,
            sum_complex_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_times, 1)) : nothing,
        ) : nothing
    out = TFOutput(times_out, freqs, g)

    _tf_foreach_channel(selected_channels, gpu_active, host_outputs -> TFBuffers(g; host_outputs)) do channel, buf
        _process_stft_channel!(
            channel,
            buf,
            sdat,
            g,
            time_indices,
            adjusted_time_indices,
            n_window_samples_per_freq,
            window_lengths_samples_exact,
            fft_plan,
            ifft_plan,
            wavelet_ffts,
            out,
            gpu_buffers,
        )
    end

    return _tf_result(sdat, out, selected_channels, is_fixed ? :taper_fixed : :taper_adaptive, g)
end

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
    ifft_plan_padded_batch::P,
    local_conv_result::Matrix{ComplexF64},
    eegpower_trials::Array{Float64,3},
    eegconv_trials::Array{ComplexF64,3},
    eegpower_avg::Matrix{Float64},
    eegconv_avg::Matrix{ComplexF64},
) where {P}
    @inbounds for trial_idx = 1:n_trials
        @simd for i = 1:n_samples_padded
            local_conv_result[i, trial_idx] = local_data_fft[i, trial_idx] * curr_wavelet[i]
        end
    end

    mul!(local_conv_result, ifft_plan_padded_batch, local_conv_result)

    @inbounds @simd for i in eachindex(local_conv_result)
        local_conv_result[i] *= norm_factor
    end

    @inbounds for ti_idx = 1:n_times
        sample_idx = time_indices[ti_idx]
        adjusted_idx = sample_idx + n_pre_pad + half_window
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
    buf::TFBuffers,
    dat::EpochData,
    g::TFGeometry,
    time_indices::Vector{Int},
    adjusted_time_indices::Vector{Int},
    n_window_samples_per_freq::Vector{Int},
    window_lengths_samples_exact::Vector{Float64},
    fft_plan_padded_batch::P1,
    ifft_plan_padded_batch::P2,
    wavelet_ffts::Vector{Vector{ComplexF64}},
    out::TFOutput,
    gpu_buffers,
) where {P1,P2}
    n_trials, n_samples_padded, num_frex, n_times, n_pre_pad = g.n_trials, g.n_fft, g.num_frex, g.n_times, g.n_pre_pad
    return_trials, return_phase = g.return_trials, g.return_phase

    _tf_reset_outputs!(isnothing(gpu_buffers) ? buf : gpu_buffers, g)

    # Extract trials, mirror-pad and batch FFT (into buf.data_fft or the GPU buffer)
    _tf_load_channel_fft!(buf, dat, channel, g, fft_plan_padded_batch, gpu_buffers)

    if !isnothing(gpu_buffers)
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

        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            norm_factor = sqrt(2.0 / n_window_samples)
            half_window = n_window_samples ÷ 2

            copyto!(gpu_buffers.curr_wavelet_gpu, wavelet_ffts[fi])
            local_conv_result_gpu .= local_data_fft_gpu .* gpu_buffers.curr_wavelet_gpu
            mul!(local_conv_result_gpu, ifft_plan_padded_batch, local_conv_result_gpu)
            local_conv_result_gpu .*= norm_factor

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
        end

        if !return_trials
            eegpower_avg_gpu ./= n_trials
            if return_phase
                eegconv_avg_gpu ./= n_trials
            end
        end
    else
        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            norm_factor = sqrt(2.0 / n_window_samples)
            half_window = n_window_samples ÷ 2
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
                buf.data_fft,
                wavelet_ffts[fi],
                ifft_plan_padded_batch,
                buf.conv_result,
                buf.power_trials,
                buf.conv_trials,
                buf.power_avg,
                buf.conv_avg,
            )
        end
    end

    # Final [trial × freq × time] or [freq × time] arrays (downloaded from the GPU if needed)
    power, conv = _tf_channel_results(buf, gpu_buffers, g)

    if g.filter_edges
        _filter_edges!(power, conv, num_frex, adjusted_time_indices, window_lengths_samples_exact, g.n_padded)
    end

    # CPU path averages over trials after edge filtering (GPU path averaged on the device above)
    if isnothing(gpu_buffers) && !return_trials
        power ./= n_trials
        if return_phase
            conv ./= n_trials
        end
    end

    _tf_store_channel!(out, channel, power, conv, g)
    return nothing
end

# Vector version: process each EpochData element
function tf_stft(data_vec::AbstractVector{<:EpochData}; kwargs...)
    return [tf_stft(dat; kwargs...) for dat in data_vec]
end


# === BATCH PROCESSING FUNCTIONS ===
function tf_stft(
    file_pattern::String;
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection = participants(),
    condition_selection = conditions(),
    kwargs...,
)
    _with_batch_logging("tf_stft", "tf_stft.log", input_dir, output_dir) do ctx
        @log_call "tf_stft"

        output_dir = _prepare_output_dir!(ctx, joinpath(input_dir, "tf_stft_$(clean_pattern(file_pattern))"))

        process_fn = _tf_batch_worker(tf_stft, "STFT", condition_selection, kwargs)

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "TF STFT")
    end
end
