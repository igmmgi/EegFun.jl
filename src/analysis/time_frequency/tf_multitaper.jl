"""
    tf_multitaper(dat::EpochData;
                  channel_selection=channels(),
                  interval_selection::Interval=times(),
                  frequencies::AbstractVector{<:Real}=range(1, 40, length=40),
                  time_steps::Real=0.05,
                  cycles::Real,
                  frequency_smoothing::Union{Nothing,Real}=nothing,
                  pad::Union{Nothing,Symbol}=nothing,
                  return_trials::Bool=false,
                  return_phase::Bool=false,
                  filter_edges::Bool=true,
                  use_gpu::Bool=false)

Multitaper time-frequency decomposition using orthogonal DPSS (Slepian) tapers (`mtmconvol`).
Uses adaptive window lengths (`cycles / frequency`) and `frequency_smoothing * frequency` bandwidth
(defaulting to `0.4 * frequency` when `frequency_smoothing = nothing`).

Returns `TimeFreqData` when `return_trials=false` or `TimeFreqEpochData` when `return_trials=true`.

# Examples
```julia
tf_data = tf_multitaper(epochs; cycles=5)
tf_log = tf_multitaper(epochs; frequencies=logrange(2, 80, length=30), cycles=5, frequency_smoothing=0.4)
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
    sdat, gpu_active = _tf_prepare(dat, pad, use_gpu, channel_selection, interval_selection, frequencies)
    selected_channels = channel_labels(sdat)
    num_frex = length(frequencies)

    cycles <= 0 && @minimal_error("`cycles` must be positive, got $cycles")
    freq_smooth = isnothing(frequency_smoothing) ? 0.4 : Float64(frequency_smoothing)
    freq_smooth <= 0 && @minimal_error("`frequency_smoothing` must be positive, got $freq_smooth")

    n_samples_per_epoch = n_samples(sdat)
    time_indices, times_out = _tf_time_grid(sdat, time_steps)
    n_trials = n_epochs(sdat)
    n_times = length(times_out)
    freqs = collect(frequencies)

    window_lengths_sec = [Float64(cycles) / freqs[fi] for fi = 1:num_frex]
    window_lengths_samples_exact = [wl * sdat.sample_rate for wl in window_lengths_sec]
    n_window_samples_per_freq = [round(Int, w) for w in window_lengths_samples_exact]
    _validate_tf_window_samples(n_window_samples_per_freq)

    n_pre_pad, _, n_padded_samples = _tf_padding(pad, n_samples_per_epoch)
    n_samples_padded = nextpow(2, n_padded_samples + maximum(n_window_samples_per_freq) - 1)
    adjusted_time_indices = time_indices .+ n_pre_pad

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

    n_tapers_per_freq = Vector{Int}(undef, num_frex)
    inv_n_tapers_per_freq = Vector{Float64}(undef, num_frex)
    tapered_wavelet_ffts = Vector{Vector{Vector{ComplexF64}}}(undef, num_frex)

    wavelet_padded = Vector{ComplexF64}(undef, n_samples_padded)
    p_fft_wavelet = plan_fft(wavelet_padded; flags = FFTW.MEASURE)
    carrier = Vector{ComplexF64}(undef, maximum(n_window_samples_per_freq))
    tapered_buf = similar(carrier)

    for fi = 1:num_frex
        freq = freqs[fi]
        n_win = n_window_samples_per_freq[fi]
        wl_sec = window_lengths_sec[fi]
        tapsmofrq = freq_smooth * freq
        NW = tapsmofrq * wl_sec / 2
        if NW <= 0 || NW >= n_win / 2
            @minimal_error(
                "Invalid time-bandwidth product NW=$NW for frequency $freq Hz. NW must be > 0 and < $(n_win/2). Window length: $wl_sec s, tapsmofrq: $tapsmofrq Hz",
            )
        end

        K = min(max(1, floor(Int, 2 * NW - 1)), n_win ÷ 2)
        dpss_res = DSP.dpss(n_win, NW, K)
        raw_tapers = dpss_res isa Tuple ? dpss_res[1] : dpss_res
        tapers = if raw_tapers isa AbstractMatrix
            Matrix{Float64}(raw_tapers)
        elseif raw_tapers isa AbstractVector
            reshape(Vector{Float64}(raw_tapers), length(raw_tapers), 1)
        else
            error("DSP.dpss returned unexpected type: $(typeof(raw_tapers)) for freq=$freq Hz, NW=$NW, K=$K")
        end

        n_tapers_per_freq[fi] = K
        inv_n_tapers_per_freq[fi] = 1.0 / K
        tapered_wavelet_ffts[fi] = Vector{Vector{ComplexF64}}(undef, K)

        step_angle = 2π * freq / sdat.sample_rate
        half = (n_win - 1) * 0.5
        @inbounds @simd for i = 1:n_win
            carrier[i] = cis(((i - 1) - half) * step_angle)
        end

        tapered_view = @view tapered_buf[1:n_win]
        for k = 1:K
            @inbounds @simd for i = 1:n_win
                tapered_view[i] = tapers[i, k] * carrier[i]
            end
            tapered_view .*= sqrt(2.0) / norm(tapered_view, 2)
            tapered_wavelet_ffts[fi][k] = _tf_padded_fft(tapered_view, p_fft_wavelet, wavelet_padded)
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
        ) : nothing
    out = TFOutput(times_out, freqs, length(selected_channels), g)

    make_buffers = host_outputs -> (TFBuffers(g; host_outputs), MultitaperAccumulators(g; host_outputs))
    _tf_foreach_channel(selected_channels, gpu_active, make_buffers) do ch_idx, channel, (buf, acc)
        _process_multitaper_channel!(
            ch_idx,
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

struct MultitaperAccumulators
    power_trials::Matrix{Float64}
    complex_trials::Matrix{ComplexF64}
    power_avg::Vector{Float64}
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
    ch_idx::Int,
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
    fft_plan::P1,
    ifft_plan::P2,
    tapered_wavelet_ffts::Vector{Vector{Vector{ComplexF64}}},
    out::TFOutput,
    gpu_buffers,
) where {P1,P2}
    n_trials, n_fft, num_frex, n_times, n_pre_pad = g.n_trials, g.n_fft, g.num_frex, g.n_times, g.n_pre_pad
    return_trials, return_phase = g.return_trials, g.return_phase
    inv_n_trials = 1.0 / n_trials

    _tf_load_channel_fft!(buf, dat, channel, g, fft_plan, gpu_buffers)

    if !isnothing(gpu_buffers)
        _tf_reset_outputs!(gpu_buffers, g)
        gb = gpu_buffers
        adjusted_indices_cpu = Vector{Int}(undef, n_times)

        for fi = 1:num_frex
            n_tapers = n_tapers_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]
            off = n_pre_pad + (n_window_samples_per_freq[fi] ÷ 2)

            @inbounds for ti = 1:n_times
                adjusted_indices_cpu[ti] = time_indices[ti] + off
            end
            copyto!(gb.adjusted_indices_gpu, adjusted_indices_cpu)

            if return_trials
                fill!(gb.power_accum_trials_gpu, 0.0)
                return_phase && fill!(gb.complex_accum_trials_gpu, 0.0im)
            else
                fill!(gb.power_accum_avg_gpu, 0.0)
                return_phase && fill!(gb.complex_accum_avg_gpu, 0.0im)
            end

            for taper_idx = 1:n_tapers
                copyto!(gb.curr_wavelet_gpu, tapered_wavelet_ffts[fi][taper_idx])
                gb.local_conv_result_gpu .= gb.local_data_fft_gpu .* gb.curr_wavelet_gpu
                mul!(gb.local_conv_result_gpu, ifft_plan, gb.local_conv_result_gpu)

                extracted_gpu = @view gb.local_conv_result_gpu[gb.adjusted_indices_gpu, :]
                if return_trials
                    gb.power_accum_trials_gpu .+= transpose(abs2.(extracted_gpu))
                    (return_phase && taper_idx == 1) && (gb.complex_accum_trials_gpu .+= transpose(extracted_gpu))
                else
                    gb.extracted_power_gpu .= abs2.(extracted_gpu)
                    sum!(gb.sum_power_gpu, gb.extracted_power_gpu)
                    gb.power_accum_avg_gpu .+= vec(gb.sum_power_gpu)
                    if return_phase && taper_idx == 1
                        sum!(gb.sum_complex_gpu, extracted_gpu)
                        gb.complex_accum_avg_gpu .+= vec(gb.sum_complex_gpu)
                    end
                end
            end

            if return_trials
                @views gb.eegpower_trials_gpu[:, fi, :] .= gb.power_accum_trials_gpu .* inv_n_tapers
                return_phase && (@views gb.eegconv_trials_gpu[:, fi, :] .= gb.complex_accum_trials_gpu)
            else
                @views gb.eegpower_avg_gpu[fi, :] .= gb.power_accum_avg_gpu .* (inv_n_tapers * inv_n_trials)
                return_phase && (@views gb.eegconv_avg_gpu[fi, :] .= gb.complex_accum_avg_gpu .* inv_n_trials)
            end
        end
    else
        local_data_fft, local_conv_result = buf.data_fft, buf.conv_result
        eegpower_trials, eegconv_trials = buf.power_trials, buf.conv_trials
        eegpower_avg, eegconv_avg = buf.power_avg, buf.conv_avg
        power_accum_trials, complex_accum_trials = acc.power_trials, acc.complex_trials
        power_accum_avg, complex_accum_avg = acc.power_avg, acc.complex_avg

        for fi = 1:num_frex
            n_tapers = n_tapers_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]
            off = n_pre_pad + (n_window_samples_per_freq[fi] ÷ 2)

            if return_trials
                fill!(power_accum_trials, 0.0)
                return_phase && fill!(complex_accum_trials, 0.0im)
            else
                fill!(power_accum_avg, 0.0)
                return_phase && fill!(complex_accum_avg, 0.0im)
            end

            for taper_idx = 1:n_tapers
                tapered_wavelet_fft = tapered_wavelet_ffts[fi][taper_idx]
                accum_phase = return_phase && taper_idx == 1

                @inbounds for tr = 1:n_trials
                    @simd for i = 1:n_fft
                        local_conv_result[i, tr] = local_data_fft[i, tr] * tapered_wavelet_fft[i]
                    end
                end
                mul!(local_conv_result, ifft_plan, local_conv_result)

                if return_trials
                    if accum_phase
                        @inbounds for ti = 1:n_times
                            adjusted_idx = time_indices[ti] + off
                            @simd for tr = 1:n_trials
                                val = local_conv_result[adjusted_idx, tr]
                                power_accum_trials[tr, ti] += abs2(val)
                                complex_accum_trials[tr, ti] += val
                            end
                        end
                    else
                        @inbounds for ti = 1:n_times
                            adjusted_idx = time_indices[ti] + off
                            @simd for tr = 1:n_trials
                                power_accum_trials[tr, ti] += abs2(local_conv_result[adjusted_idx, tr])
                            end
                        end
                    end
                elseif accum_phase
                    @inbounds for ti = 1:n_times
                        adjusted_idx = time_indices[ti] + off
                        p_sum = 0.0
                        c_sum = 0.0im
                        @simd for tr = 1:n_trials
                            val = local_conv_result[adjusted_idx, tr]
                            p_sum += abs2(val)
                            c_sum += val
                        end
                        power_accum_avg[ti] += p_sum
                        complex_accum_avg[ti] += c_sum
                    end
                else
                    @inbounds for ti = 1:n_times
                        adjusted_idx = time_indices[ti] + off
                        p_sum = 0.0
                        @simd for tr = 1:n_trials
                            p_sum += abs2(local_conv_result[adjusted_idx, tr])
                        end
                        power_accum_avg[ti] += p_sum
                    end
                end
            end

            if return_trials
                @inbounds for ti = 1:n_times
                    @simd for tr = 1:n_trials
                        eegpower_trials[tr, fi, ti] = power_accum_trials[tr, ti] * inv_n_tapers
                    end
                    if return_phase
                        @simd for tr = 1:n_trials
                            eegconv_trials[tr, fi, ti] = complex_accum_trials[tr, ti]
                        end
                    end
                end
            else
                factor = inv_n_tapers * inv_n_trials
                @inbounds for ti = 1:n_times
                    eegpower_avg[fi, ti] = power_accum_avg[ti] * factor
                    if return_phase
                        eegconv_avg[fi, ti] = complex_accum_avg[ti] * inv_n_trials
                    end
                end
            end
        end
    end

    power, conv = _tf_channel_results(buf, gpu_buffers, g)
    g.filter_edges && _filter_edges!(power, conv, num_frex, adjusted_time_indices, window_lengths_samples_exact, g.n_padded)
    _tf_store_channel!(out, ch_idx, power, conv, g)
    return nothing
end

tf_multitaper(data_vec::AbstractVector{<:EpochData}; kwargs...) = [tf_multitaper(dat; kwargs...) for dat in data_vec]

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
