"""
    tf_multitaper(dat::EpochData; 
                  channel_selection::Function=channels(),
                  interval_selection::Interval=samples(),
                  frequencies::Union{AbstractRange,AbstractVector{<:Real}}=range(1, 40, length=40),
                  cycles::Real,
                  frequency_smoothing::Union{Nothing,Real}=nothing,
                  time_steps::Real=0.05,
                  pad::Union{Nothing,Symbol}=nothing,
                  return_trials::Bool=false,
                  filter_edges::Bool=true)

Multitaper time-frequency analysis using DPSS (Discrete Prolate Spheroidal Sequences) tapers (Cohen Chapter 16, equivalent to FieldTrip's 'mtmconvol' method).

Uses multiple orthogonal tapers (Slepian sequences) to reduce variance in spectral estimates compared to single-taper methods. Uses adaptive window lengths (cycles per frequency) to match the time-frequency trade-off of Morlet wavelets and adaptive STFT.

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
- `cycles::Real`: Number of cycles per frequency. Window length = cycles / frequency (in seconds).
  - Example: `cycles=5` uses 5 cycles for all frequencies 
  - Lower frequencies will have longer windows, higher frequencies will have shorter windows
- `frequency_smoothing::Union{Nothing,Real}=nothing`: Frequency smoothing parameter 
  - If `nothing`, uses `frequency_smoothing = 0.4 * frequency` 
  - If a number, uses that value multiplied by frequency: `tapsmofrq = frequency_smoothing * frequency`
  - Example: `frequency_smoothing=0.4` 
  - Controls time-bandwidth product: `NW = tapsmofrq * window_length / 2`
  - Number of tapers used: `K = 2*NW - 1` (rounded down)
- `time_steps::Real=0.05`: Step size for extracting time points in seconds.
  - Creates time points from the selected time range with the specified step size
  - Default: `0.05` (50 ms)
  - Example: `time_steps=0.01` creates time points every 0.01 seconds within the selected time interval
- `pad::Union{Nothing,Symbol}=nothing`: Padding method to reduce edge artifacts. Options:
  - `nothing`: No padding (default)
  - `:pre`: Mirror data before each epoch
  - `:post`: Mirror data after each epoch
  - `:both`: Mirror data on both sides (recommended)
- `return_trials::Bool=false`: If `true`, returns `TimeFreqEpochData` with individual trials preserved.
  - If `false` (default), returns `TimeFreqData` with trials averaged.
- `filter_edges::Bool=true`: If `true` (default), filters out edge regions where the window extends beyond the data

# Returns
- `TimeFreqData` (if `return_trials=false`): Time-frequency data with trials averaged
- `TimeFreqEpochData` (if `return_trials=true`): Time-frequency data with individual trials preserved

# Examples
```julia
# Default: linear frequencies 1-40 Hz, 40 points
tf_data = tf_multitaper(epochs; cycles=5)

# Log-spaced frequencies with default frequency smoothing (0.4 * frequency)
tf_data = tf_multitaper(epochs; frequencies=logrange(2, 80, length=30), cycles=5)

# Custom frequency smoothing
tf_data = tf_multitaper(epochs; frequencies=1:2:30, cycles=5, frequency_smoothing=0.4)
```
"""
function tf_multitaper(
    dat::EpochData;
    channel_selection = channels(),
    interval_selection::Interval = times(),
    frequencies::AbstractVector{<:Real} = range(1, 40, length = 40),
    time_steps::Real = 0.05,
    cycles::Real,
    frequency_smoothing::Union{Nothing,Real} = nothing,
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

    # Validate cycles parameter (adaptive window only)
    if cycles <= 0
        @minimal_error("`cycles` must be positive, got $cycles")
    end

    # Default frequency smoothing 
    if isnothing(frequency_smoothing)
        frequency_smoothing = 0.4
    end
    if frequency_smoothing <= 0
        @minimal_error("`frequency_smoothing` must be positive, got $frequency_smoothing")
    end

    # Get time points and determine output time grid
    n_samples_per_epoch = n_samples(sdat)
    time_indices, times_out = _tf_time_grid(sdat, time_steps)

    # Get number of trials/epochs
    n_trials = n_epochs(sdat)

    # Use frequencies directly (convert to vector if needed for indexing)
    freqs = collect(frequencies)

    # Adaptive window mode: window length = cycles / frequency 
    cycles_vec = fill(Float64(cycles), num_frex)
    window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * sdat.sample_rate)) for wl in window_lengths_sec]
    _validate_tf_window_samples(n_window_samples_per_freq)

    n_times = length(times_out)

    # Pre-compute DPSS tapers for each frequency
    tapers_per_freq = Vector{Matrix{Float64}}(undef, num_frex)  # Each element is (n_window_samples × n_tapers)
    n_tapers_per_freq = zeros(Int, num_frex)
    inv_n_tapers_per_freq = Vector{Float64}(undef, num_frex)

    # Pre-calculate padding amounts for virtual padding (no new dataframe allocations)
    n_pre_pad, _, n_padded_samples = _tf_padding(pad, n_samples_per_epoch)

    # Maximum window dictates if we need additional padding for FFT beyond our explicit pad
    max_window_samples = maximum(n_window_samples_per_freq)
    n_samples_padded = nextpow(2, max(n_padded_samples, max_window_samples))

    adjusted_time_indices = time_indices .+ n_pre_pad
    window_lengths_samples_exact = [(cycles / freqs[fi]) * sdat.sample_rate for fi = 1:num_frex]

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

    # Pre-compute tapered wavelets and their FFTs for each frequency-taper combination
    # For multitaper: tapered_wavelet = taper * (cos + i*sin) at target frequency
    tapered_wavelet_ffts = Vector{Vector{Vector{ComplexF64}}}(undef, num_frex)  # [frequency][taper] = FFT of tapered wavelet

    # Pre-allocate wavelet buffer and FFT plan
    wavelet_padded = zeros(ComplexF64, n_samples_padded)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)

    for fi = 1:num_frex
        freq = freqs[fi]
        n_window_samples = n_window_samples_per_freq[fi]
        window_length_sec = window_lengths_sec[fi]

        # Calculate frequency smoothing 
        tapsmofrq = frequency_smoothing * freq  # Frequency smoothing in Hz

        # Time-bandwidth product: NW = tapsmofrq * window_length / 2
        NW = tapsmofrq * window_length_sec / 2

        # Ensure NW is valid (must be > 0 and < n_window_samples/2)
        if NW <= 0 || NW >= n_window_samples / 2
            @minimal_error(
                "Invalid time-bandwidth product NW=$NW for frequency $freq Hz. NW must be > 0 and < $(n_window_samples/2). Window length: $window_length_sec s, tapsmofrq: $tapsmofrq Hz",
            )
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
        angle_in =
            range(-(n_window_samples - 1) / 2, (n_window_samples - 1) / 2, length = n_window_samples) .* (2π * freq / sdat.sample_rate)
        cos_wav = cos.(angle_in)
        sin_wav = sin.(angle_in)

        # For each taper, create tapered wavelet and compute FFT using pre-planned FFT
        for taper_idx = 1:K
            taper_col = @view tapers_per_freq[fi][:, taper_idx]
            # Tapered wavelet: taper * (cos + i*sin)
            tapered_wavelet = taper_col .* (cos_wav .+ im .* sin_wav)
            tapered_wavelet_ffts[fi][taper_idx] = _tf_padded_fft(tapered_wavelet, p_fft_wavelet, wavelet_padded)
        end
    end

    fft_plan, ifft_plan = _tf_batch_fft_plans(n_samples_padded, n_trials, gpu_active)
    gpu_buffers =
        gpu_active ?
        _tf_gpu_buffers(
            g;
            power_accum_trials_gpu = return_trials ? gpu_array(zeros(Float64, n_trials, n_times)) : nothing,
            complex_accum_trials_gpu = (return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_trials, n_times)) : nothing,
            power_accum_avg_gpu = !return_trials ? gpu_array(zeros(Float64, n_times)) : nothing,
            complex_accum_avg_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_times)) : nothing,
            extracted_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_times, n_trials)) : nothing,
            sum_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_times, 1)) : nothing,
            sum_complex_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_times, 1)) : nothing,
        ) : nothing
    out = TFOutput(times_out, freqs, g)

    make_buffers = host_outputs -> (TFBuffers(g; host_outputs), MultitaperAccumulators(g; host_outputs))
    _tf_foreach_channel(selected_channels, gpu_active, make_buffers) do channel, (buf, acc)
        _process_multitaper_channel!(
            channel,
            buf,
            acc,
            sdat,
            g,
            time_indices,
            adjusted_time_indices,
            n_window_samples_per_freq,
            window_lengths_samples_exact,
            n_tapers_per_freq,
            inv_n_tapers_per_freq,
            fft_plan,
            ifft_plan,
            tapered_wavelet_ffts,
            out,
            gpu_buffers,
        )
    end

    return _tf_result(sdat, out, selected_channels, :multitaper, g)
end

"""Per-thread host accumulators summing over tapers (empty when not needed or on the GPU)."""
struct MultitaperAccumulators
    power_trials::Matrix{Float64}      # n_trials × n_times
    complex_trials::Matrix{ComplexF64}
    power_avg::Vector{Float64}         # n_times
    complex_avg::Vector{ComplexF64}
end

function MultitaperAccumulators(g::TFGeometry; host_outputs::Bool = true)
    rt, rp = g.return_trials, g.return_phase
    return MultitaperAccumulators(
        (host_outputs && rt) ? Matrix{Float64}(undef, g.n_trials, g.n_times) : Matrix{Float64}(undef, 0, 0),
        (host_outputs && rt && rp) ? Matrix{ComplexF64}(undef, g.n_trials, g.n_times) : Matrix{ComplexF64}(undef, 0, 0),
        (host_outputs && !rt) ? Vector{Float64}(undef, g.n_times) : Vector{Float64}(undef, 0),
        (host_outputs && !rt && rp) ? Vector{ComplexF64}(undef, g.n_times) : Vector{ComplexF64}(undef, 0),
    )
end

function _process_multitaper_channel!(
    channel::Symbol,
    buf::TFBuffers,
    acc::MultitaperAccumulators,
    dat::EpochData,
    g::TFGeometry,
    time_indices::Vector{Int},
    adjusted_time_indices::Vector{Int},
    n_window_samples_per_freq::Vector{Int},
    window_lengths_samples_exact::Vector{Float64},
    n_tapers_per_freq::Vector{Int},
    inv_n_tapers_per_freq::Vector{Float64},
    fft_plan_padded_batch::P1,
    ifft_plan_padded_batch::P2,
    tapered_wavelet_ffts::Vector{Vector{Vector{ComplexF64}}},
    out::TFOutput,
    gpu_buffers,
) where {P1,P2}
    n_trials, n_samples_padded, num_frex, n_times, n_pre_pad = g.n_trials, g.n_fft, g.num_frex, g.n_times, g.n_pre_pad
    return_trials, return_phase = g.return_trials, g.return_phase

    _tf_reset_outputs!(isnothing(gpu_buffers) ? buf : gpu_buffers, g)

    # Extract trials, mirror-pad and batch FFT (into buf.data_fft or the GPU buffer)
    _tf_load_channel_fft!(buf, dat, channel, g, fft_plan_padded_batch, gpu_buffers)

    inv_n_trials = 1.0 / n_trials

    if !isnothing(gpu_buffers)
        local_data_fft_gpu = gpu_buffers.local_data_fft_gpu
        local_conv_result_gpu = gpu_buffers.local_conv_result_gpu
        eegpower_trials_gpu = gpu_buffers.eegpower_trials_gpu
        eegconv_trials_gpu = gpu_buffers.eegconv_trials_gpu
        eegpower_avg_gpu = gpu_buffers.eegpower_avg_gpu
        eegconv_avg_gpu = gpu_buffers.eegconv_avg_gpu
        power_accum_trials_gpu = gpu_buffers.power_accum_trials_gpu
        complex_accum_trials_gpu = gpu_buffers.complex_accum_trials_gpu
        power_accum_avg_gpu = gpu_buffers.power_accum_avg_gpu
        complex_accum_avg_gpu = gpu_buffers.complex_accum_avg_gpu
        adjusted_indices_gpu = gpu_buffers.adjusted_indices_gpu
        extracted_power_gpu = gpu_buffers.extracted_power_gpu
        sum_power_gpu = gpu_buffers.sum_power_gpu
        sum_complex_gpu = gpu_buffers.sum_complex_gpu
        adjusted_indices_cpu = Vector{Int}(undef, n_times)

        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            n_tapers = n_tapers_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]
            norm_factor = sqrt(2.0 / n_window_samples)
            half_window = n_window_samples ÷ 2

            for ti_idx = 1:n_times
                adjusted_indices_cpu[ti_idx] = min(max(time_indices[ti_idx] + n_pre_pad + half_window, 1), n_samples_padded)
            end
            copyto!(adjusted_indices_gpu, adjusted_indices_cpu)

            if return_trials
                fill!(power_accum_trials_gpu, 0.0)
                return_phase && fill!(complex_accum_trials_gpu, 0.0im)
            else
                fill!(power_accum_avg_gpu, 0.0)
                return_phase && fill!(complex_accum_avg_gpu, 0.0im)
            end

            for taper_idx = 1:n_tapers
                copyto!(gpu_buffers.curr_wavelet_gpu, tapered_wavelet_ffts[fi][taper_idx])
                local_conv_result_gpu .= local_data_fft_gpu .* gpu_buffers.curr_wavelet_gpu
                mul!(local_conv_result_gpu, ifft_plan_padded_batch, local_conv_result_gpu)
                local_conv_result_gpu .*= norm_factor

                extracted_gpu = @view local_conv_result_gpu[adjusted_indices_gpu, :]

                if return_trials
                    power_accum_trials_gpu .+= transpose(abs2.(extracted_gpu))
                    if return_phase
                        complex_accum_trials_gpu .+= transpose(extracted_gpu)
                    end
                else
                    extracted_power_gpu .= abs2.(extracted_gpu)
                    sum!(sum_power_gpu, extracted_power_gpu)
                    power_accum_avg_gpu .+= vec(sum_power_gpu)

                    if return_phase
                        sum!(sum_complex_gpu, extracted_gpu)
                        complex_accum_avg_gpu .+= vec(sum_complex_gpu)
                    end
                end
            end

            if return_trials
                @views eegpower_trials_gpu[:, fi, :] .= power_accum_trials_gpu .* inv_n_tapers
                if return_phase
                    @views eegconv_trials_gpu[:, fi, :] .= complex_accum_trials_gpu .* inv_n_tapers
                end
            else
                @views eegpower_avg_gpu[fi, :] .= power_accum_avg_gpu .* (inv_n_tapers * inv_n_trials)
                if return_phase
                    @views eegconv_avg_gpu[fi, :] .= complex_accum_avg_gpu .* (inv_n_tapers * inv_n_trials)
                end
            end
        end
    else
        local_data_fft = buf.data_fft
        local_conv_result = buf.conv_result
        eegpower_trials = buf.power_trials
        eegconv_trials = buf.conv_trials
        eegpower_avg = buf.power_avg
        eegconv_avg = buf.conv_avg
        power_accum_trials = acc.power_trials
        complex_accum_trials = acc.complex_trials
        power_accum_avg = acc.power_avg
        complex_accum_avg = acc.complex_avg

        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            n_tapers = n_tapers_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]
            norm_factor = sqrt(2.0 / n_window_samples)
            half_window = n_window_samples ÷ 2

            if return_trials
                fill!(power_accum_trials, 0.0)
                return_phase && fill!(complex_accum_trials, 0.0im)
            else
                fill!(power_accum_avg, 0.0)
                return_phase && fill!(complex_accum_avg, 0.0im)
            end

            for taper_idx = 1:n_tapers
                tapered_wavelet_fft = tapered_wavelet_ffts[fi][taper_idx]

                @inbounds for trial_idx = 1:n_trials
                    @simd for i = 1:n_samples_padded
                        local_conv_result[i, trial_idx] = local_data_fft[i, trial_idx] * tapered_wavelet_fft[i]
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
                            power_accum_trials[trial_idx, ti_idx] += abs2(val)
                            if return_phase
                                complex_accum_trials[trial_idx, ti_idx] += val
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
                        power_accum_avg[ti_idx] += power_sum
                        if return_phase
                            complex_accum_avg[ti_idx] += conv_sum
                        end
                    end
                end
            end

            if return_trials
                @inbounds for ti_idx = 1:n_times
                    @simd for trial_idx = 1:n_trials
                        eegpower_trials[trial_idx, fi, ti_idx] = power_accum_trials[trial_idx, ti_idx] * inv_n_tapers
                    end
                    if return_phase
                        @simd for trial_idx = 1:n_trials
                            eegconv_trials[trial_idx, fi, ti_idx] = complex_accum_trials[trial_idx, ti_idx] * inv_n_tapers
                        end
                    end
                end
            else
                factor = inv_n_tapers * inv_n_trials
                @inbounds for ti_idx = 1:n_times
                    eegpower_avg[fi, ti_idx] = power_accum_avg[ti_idx] * factor
                    if return_phase
                        eegconv_avg[fi, ti_idx] = complex_accum_avg[ti_idx] * factor
                    end
                end
            end
        end
    end

    # Final [trial × freq × time] or [freq × time] arrays (downloaded from the GPU if needed)
    power, conv = _tf_channel_results(buf, gpu_buffers, g)

    if g.filter_edges
        _filter_edges!(power, conv, num_frex, adjusted_time_indices, window_lengths_samples_exact, g.n_padded)
    end

    _tf_store_channel!(out, channel, power, conv, g)
    return nothing
end

# Vector version: process each EpochData element
function tf_multitaper(data_vec::AbstractVector{<:EpochData}; kwargs...)
    return [tf_multitaper(dat; kwargs...) for dat in data_vec]
end


# === BATCH PROCESSING FUNCTIONS ===
function tf_multitaper(
    file_pattern::String;
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection = participants(),
    condition_selection = conditions(),
    kwargs...,
)
    _with_batch_logging("tf_multitaper", "tf_multitaper.log", input_dir, output_dir) do ctx
        @log_call "tf_multitaper"

        output_dir = _prepare_output_dir!(ctx, joinpath(input_dir, "tf_multitaper_$(clean_pattern(file_pattern))"))

        process_fn = _tf_batch_worker(tf_multitaper, "multitaper", condition_selection, kwargs)

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "TF multitaper")
    end
end
