function tf_morlet(
    dat::EpochData;
    channel_selection = channels(),
    interval_selection::Interval = times(),
    frequencies::AbstractVector{<:Real} = range(1, 40, length = 40),
    cycles::Union{Real,Tuple{Real,Real}} = 7,
    pad::Union{Nothing,Symbol} = nothing,
    return_trials::Bool = false,
    return_phase::Bool = false,
    filter_edges::Bool = true,
    use_gpu::Bool = false,
)

    # Validate arguments and subset data with channel and interval selection
    sdat, gpu_active = _tf_prepare(dat, pad, use_gpu, channel_selection, interval_selection, frequencies)

    # Original (unpadded) samples are the time points we want in the output
    n_samples_unpadded = n_samples(sdat)
    n_trials = n_epochs(sdat)
    times_out = time_vector(sdat)

    # Use frequency input directly (ranges and vectors both work)
    num_frex = length(frequencies)

    # Define cycles/sigma
    cycles_vec = cycles isa Tuple ? logrange(cycles[1], cycles[2], length = num_frex) : fill(cycles, num_frex)

    # Calculate which time indices to keep (original unpadded samples)
    n_pre_pad, _, n_padded_samples = _tf_padding(pad, n_samples_unpadded)
    time_indices_out = (n_pre_pad+1):(n_pre_pad+n_samples_unpadded)

    # Pre-compute convolution length for single trial processing
    max_sigma = maximum(cycles_vec ./ (2 * pi .* frequencies))
    max_hw = ceil(Int, 6 * max_sigma * sdat.sample_rate) ÷ 2
    max_wl = max_hw * 2 + 1
    n_conv_pow2 = nextpow(2, max_wl + n_padded_samples - 1)

    g = TFGeometry(
        n_trials,
        n_samples_unpadded,
        n_pre_pad,
        n_padded_samples,
        n_conv_pow2,
        num_frex,
        length(time_indices_out),
        return_trials,
        return_phase,
        filter_edges,
    )

    # Pre-allocate wavelet buffer
    wavelet_padded = zeros(ComplexF64, n_conv_pow2)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)

    # Pre-compute some constants 
    inv_sr = 1.0 / sdat.sample_rate
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    wl_per_freq = Vector{Int}(undef, num_frex)
    valid_starts_per_freq = Vector{Int}(undef, num_frex)

    for fi = 1:num_frex
        freq_val = frequencies[fi]
        sigma = cycles_vec[fi] / (two_pi * freq_val)
        hw = ceil(Int, 6 * sigma * sdat.sample_rate) ÷ 2
        wl = hw * 2 + 1
        wl_per_freq[fi] = wl
        valid_start = hw + 1
        valid_starts_per_freq[fi] = valid_start

        # Create wavelet directly in padded buffer
        fill!(wavelet_padded, 0)
        A = sqrt(1 / (sigma * sqrt_pi))
        two_pi_freq = two_pi * freq_val
        inv_2sigma2 = 1.0 / (2 * sigma^2)
        @inbounds @simd for i = 1:wl
            t = (i - 1 - hw) * inv_sr
            t2 = t * t
            wavelet_padded[i] = A * exp(im * two_pi_freq * t) * exp(-t2 * inv_2sigma2)
        end

        # FFT of wavelet - compute once, reuse for all channels and trials
        wavelet_fft_freq = zeros(ComplexF64, n_conv_pow2)
        mul!(wavelet_fft_freq, p_fft_wavelet, wavelet_padded)
        wavelet_ffts[fi] = wavelet_fft_freq
    end

    fft_plan, ifft_plan = _tf_batch_fft_plans(n_conv_pow2, n_trials, gpu_active)
    gpu_buffers = gpu_active ? _tf_gpu_buffers(g) : nothing
    out = TFOutput(times_out, frequencies, g)

    # Process each selected channel (all trials batched)
    channel_names = channel_labels(sdat)
    _tf_foreach_channel(channel_names, gpu_active, host_outputs -> TFBuffers(g; host_outputs)) do channel, buf
        _process_morlet_channel!(
            channel,
            buf,
            sdat,
            g,
            valid_starts_per_freq,
            time_indices_out,
            wl_per_freq,
            fft_plan,
            ifft_plan,
            wavelet_ffts,
            out,
            gpu_buffers,
        )
    end

    return _tf_result(sdat, out, channel_names, :wavelet, g)
end

function _process_morlet_channel!(
    channel::Symbol,
    buf::TFBuffers,
    dat::EpochData,
    g::TFGeometry,
    valid_starts_per_freq::Vector{Int},
    time_indices_out,
    wl_per_freq::Vector{Int},
    fft_plan_padded_batch,
    ifft_plan_padded_batch,
    wavelet_ffts::Vector{Vector{ComplexF64}},
    out::TFOutput,
    gpu_buffers,
)
    n_trials, n_conv_pow2, num_frex, n_times_out = g.n_trials, g.n_fft, g.num_frex, g.n_times
    return_trials, return_phase = g.return_trials, g.return_phase

    # Reset output arrays (morlet accumulates into them)
    _tf_reset_outputs!(isnothing(gpu_buffers) ? buf : gpu_buffers, g)

    # Extract trials, mirror-pad and batch FFT (into buf.data_fft or the GPU buffer)
    _tf_load_channel_fft!(buf, dat, channel, g, fft_plan_padded_batch, gpu_buffers)

    # Process each frequency
    sqrt2 = sqrt(2)
    if !isnothing(gpu_buffers)
        local_data_fft_gpu = gpu_buffers.local_data_fft_gpu
        local_conv_result_gpu = gpu_buffers.local_conv_result_gpu
        curr_wavelet_gpu = gpu_buffers.curr_wavelet_gpu
        eegpower_trials_gpu = gpu_buffers.eegpower_trials_gpu
        eegconv_trials_gpu = gpu_buffers.eegconv_trials_gpu
        eegpower_avg_gpu = gpu_buffers.eegpower_avg_gpu
        eegconv_avg_gpu = gpu_buffers.eegconv_avg_gpu
        adjusted_indices_gpu = gpu_buffers.adjusted_indices_gpu
        adjusted_indices_cpu = Vector{Int}(undef, n_times_out)

        for fi = 1:num_frex
            # GPU inner loop broadcast
            copyto!(curr_wavelet_gpu, wavelet_ffts[fi])
            local_conv_result_gpu .= local_data_fft_gpu .* curr_wavelet_gpu
            mul!(local_conv_result_gpu, ifft_plan_padded_batch, local_conv_result_gpu)
            local_conv_result_gpu .*= sqrt2

            # Extract requested time points (precompute valid indices)
            valid_start = valid_starts_per_freq[fi]
            for ti_out = 1:n_times_out
                adjusted_indices_cpu[ti_out] = valid_start + time_indices_out[ti_out] - 1
            end
            copyto!(adjusted_indices_gpu, adjusted_indices_cpu)
            extracted_gpu = @view local_conv_result_gpu[adjusted_indices_gpu, :]

            if return_trials
                @views eegpower_trials_gpu[:, fi, :] .= transpose(abs2.(extracted_gpu))
                if return_phase
                    @views eegconv_trials_gpu[:, fi, :] .= transpose(extracted_gpu)
                end
            else
                @views eegpower_avg_gpu[fi, :] .= vec(sum(abs2.(extracted_gpu), dims = 2))
                if return_phase
                    @views eegconv_avg_gpu[fi, :] .= vec(sum(extracted_gpu, dims = 2))
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
        local_data_fft = buf.data_fft
        local_conv_result = buf.conv_result
        eegpower_trials = buf.power_trials
        eegconv_trials = buf.conv_trials
        eegpower_avg = buf.power_avg
        eegconv_avg = buf.conv_avg

        for fi = 1:num_frex
            # Frequency-domain convolution: multiply data FFT by wavelet FFT (broadcast across all trials)
            # We use explicit loops to avoid .* broadcast allocations on huge matrices
            curr_wavelet = wavelet_ffts[fi]
            @inbounds for trial_idx = 1:n_trials
                @simd for i = 1:n_conv_pow2
                    local_conv_result[i, trial_idx] = local_data_fft[i, trial_idx] * curr_wavelet[i]
                end
            end

            # IFFT to get time-domain result (in-place along dim 2)
            mul!(local_conv_result, ifft_plan_padded_batch, local_conv_result)

            # Scale to approximate MNE's normalization
            @inbounds @simd for i in eachindex(local_conv_result)
                local_conv_result[i] *= sqrt2
            end

            # Extract requested time points (time_indices are sample indices in processed data)
            valid_start = valid_starts_per_freq[fi]
            @inbounds for ti_out = 1:n_times_out
                ti_padded = time_indices_out[ti_out]
                idx = valid_start + ti_padded - 1

                if return_trials
                    @inbounds @simd for trial_idx = 1:n_trials
                        val = local_conv_result[idx, trial_idx]
                        eegpower_trials[trial_idx, fi, ti_out] = abs2(val)
                        if return_phase
                            eegconv_trials[trial_idx, fi, ti_out] = val
                        end
                    end
                else
                    power_sum = 0.0
                    conv_sum = 0.0im
                    @inbounds @simd for trial_idx = 1:n_trials
                        val = local_conv_result[idx, trial_idx]
                        power_sum += abs2(val)
                        if return_phase
                            conv_sum += val
                        end
                    end
                    eegpower_avg[fi, ti_out] += power_sum
                    if return_phase
                        eegconv_avg[fi, ti_out] += conv_sum
                    end
                end
            end
        end

        if !return_trials
            eegpower_avg ./= n_trials
            if return_phase
                eegconv_avg ./= n_trials
            end
        end
    end

    # Final [trial × freq × time] or [freq × time] arrays (downloaded from the GPU if needed)
    power, conv = _tf_channel_results(buf, gpu_buffers, g)

    if g.filter_edges
        # Use padded indices for edge filtering - padding extends valid region
        _filter_edges!(power, conv, num_frex, time_indices_out, wl_per_freq, g.n_padded)
    end

    _tf_store_channel!(out, channel, power, conv, g)
    return nothing
end

# Vector version: process each EpochData element
function tf_morlet(data_vec::AbstractVector{<:EpochData}; kwargs...)
    return [tf_morlet(dat; kwargs...) for dat in data_vec]
end

# === BATCH PROCESSING FUNCTIONS ===

"""
    tf_morlet(file_pattern::String;
              input_dir::String=pwd(),
              output_dir::Union{String,Nothing}=nothing,
              participant_selection::Function=participants(),
              condition_selection::Function=conditions(),
              kwargs...)

Batch Morlet wavelet time-frequency analysis on EpochData files.

Loads `.jld2` files containing `EpochData` (one or more conditions per file),
runs `tf_morlet` on each condition, and saves the results as `Vector{TimeFreqData}`.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "epochs", "epochs_unrejected")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- All other keyword arguments are forwarded to `tf_morlet` (e.g., `frequencies`, `cycles`, `pad`, `filter_edges`)

# Examples
```julia
# Run Morlet wavelet TF on all epoch files
tf_morlet("epochs_unrejected"; cycles=7)

# Specific participants, custom frequencies
tf_morlet("epochs_unrejected";
          participant_selection=participants([1, 2, 3]),
          frequencies=logrange(2, 40, length=30),
          cycles=(3, 10))

# Custom output directory
tf_morlet("epochs"; cycles=7, output_dir="/data/tf_results")
```
"""
function tf_morlet(
    file_pattern::String;
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection = participants(),
    condition_selection = conditions(),
    kwargs...,
)
    _with_batch_logging("tf_morlet", "tf_morlet.log", input_dir, output_dir) do ctx
        @log_call "tf_morlet"

        output_dir = _prepare_output_dir!(ctx, joinpath(input_dir, "tf_morlet_$(clean_pattern(file_pattern))"))

        process_fn = _tf_batch_worker(tf_morlet, "morlet", condition_selection, kwargs)

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "TF morlet")
    end
end
