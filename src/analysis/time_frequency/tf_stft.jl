"""
    tf_stft(dat::EpochData;
            channel_selection=channels(),
            interval_selection::Interval=times(),
            frequencies::AbstractVector{<:Real}=range(1, 40, length=40),
            time_steps::Real=0.05,
            window_length::Union{Nothing,Real}=nothing,
            cycles::Union{Nothing,Real}=nothing,
            pad::Union{Nothing,Symbol}=nothing,
            return_trials::Bool=false,
            return_phase::Bool=false,
            filter_edges::Bool=true,
            use_gpu::Bool=false)

Short-Time Fourier Transform (STFT) time-frequency decomposition using a sliding Hanning window.
Specify **exactly one** of `window_length` (fixed window in seconds) or `cycles` (adaptive window `cycles / freq`).

Returns `TimeFreqData` when `return_trials=false` or `TimeFreqEpochData` when `return_trials=true`.

# Examples
```julia
tf_data = tf_stft(epochs; window_length=0.3)
tf_adaptive = tf_stft(epochs; frequencies=2:1:30, cycles=7, interval_selection=times(-0.5, 1.5), time_steps=0.05)
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
    sdat, gpu_active = _tf_prepare(dat, pad, use_gpu, channel_selection, interval_selection, frequencies)
    selected_channels = channel_labels(sdat)
    num_frex = length(frequencies)

    isnothing(window_length) &&
        isnothing(cycles) &&
        @minimal_error("Either `window_length` (fixed window) or `cycles` (adaptive window) must be specified")
    !isnothing(window_length) && !isnothing(cycles) && @minimal_error("Only one of `window_length` or `cycles` can be specified, not both")
    !isnothing(window_length) && window_length <= 0 && @minimal_error("`window_length` must be positive, got $window_length")
    !isnothing(cycles) && cycles <= 0 && @minimal_error("`cycles` must be positive, got $cycles")

    time_indices, times_out = _tf_time_grid(sdat, time_steps)
    n_trials = n_epochs(sdat)
    n_samples_per_epoch = n_samples(sdat)
    freqs = collect(frequencies)

    is_fixed, n_window_samples_per_freq, window_lengths_samples_exact = if !isnothing(window_length)
        n_win_fixed = round(Int, window_length * sdat.sample_rate)
        _validate_tf_window_samples(n_win_fixed)
        (true, fill(n_win_fixed, num_frex), fill(Float64(window_length * sdat.sample_rate), num_frex))
    else
        exact = [(cycles / freqs[fi]) * sdat.sample_rate for fi = 1:num_frex]
        n_win_adapt = [round(Int, w) for w in exact]
        _validate_tf_window_samples(n_win_adapt)
        (false, n_win_adapt, exact)
    end

    n_times = length(times_out)
    n_pre_pad, _, n_padded_samples = _tf_padding(pad, n_samples_per_epoch)
    adjusted_time_indices = time_indices .+ n_pre_pad
    sample_offsets = [n_pre_pad + (n_window_samples_per_freq[fi] ÷ 2) for fi = 1:num_frex]
    n_samples_padded = nextpow(2, n_padded_samples + maximum(n_window_samples_per_freq) - 1)

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

    wavelet_padded = Vector{ComplexF64}(undef, n_samples_padded)
    p_fft_wavelet = plan_fft(wavelet_padded; flags = FFTW.MEASURE)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    fixed_hann = if is_fixed
        w = DSP.hanning(n_window_samples_per_freq[1])
        w ./ norm(w, 2)
    else
        Float64[]
    end

    for fi = 1:num_frex
        n_win = n_window_samples_per_freq[fi]
        hann = is_fixed ? fixed_hann : (w = DSP.hanning(n_win); w ./ norm(w, 2))
        step_angle = 2π * freqs[fi] / sdat.sample_rate
        half = (n_win - 1) * 0.5
        sqrt2 = sqrt(2.0)
        @inbounds @simd for i = 1:n_win
            wavelet_padded[i] = (hann[i] * sqrt2) * cis(((i - 1) - half) * step_angle)
        end
        @inbounds @simd for i = (n_win+1):n_samples_padded
            wavelet_padded[i] = zero(ComplexF64)
        end
        out_fft = Vector{ComplexF64}(undef, n_samples_padded)
        mul!(out_fft, p_fft_wavelet, wavelet_padded)
        wavelet_ffts[fi] = out_fft
    end

    fft_plan, ifft_plan = _tf_batch_fft_plans(n_samples_padded, n_trials, gpu_active)
    gpu_buffers = gpu_active ? _tf_gpu_buffers(g) : nothing
    out = TFOutput(times_out, freqs, length(selected_channels), g)

    _tf_foreach_channel(selected_channels, gpu_active, host_outputs -> TFBuffers(g; host_outputs)) do ch_idx, channel, buf
        _tf_single_taper_channel!(
            ch_idx,
            channel,
            buf,
            sdat,
            g,
            time_indices,
            sample_offsets,
            adjusted_time_indices,
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

tf_stft(data_vec::AbstractVector{<:EpochData}; kwargs...) = [tf_stft(dat; kwargs...) for dat in data_vec]

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
