"""
    tf_morlet(dat::EpochData;
              channel_selection=channels(),
              interval_selection::Interval=times(),
              frequencies::AbstractVector{<:Real}=range(1, 40, length=40),
              cycles::Union{Real,Tuple{Real,Real}}=7,
              pad::Union{Nothing,Symbol}=nothing,
              return_trials::Bool=false,
              return_phase::Bool=false,
              filter_edges::Bool=true,
              use_gpu::Bool=false)

Morlet wavelet time-frequency decomposition using MNE-style discrete ``L^2`` wavelet
normalization (``\\|W\\|_2 = \\sqrt{2}``).

Returns `TimeFreqData` when `return_trials=false` or `TimeFreqEpochData` when `return_trials=true`.

# Examples
```julia
tf_data = tf_morlet(epochs; frequencies=range(2, 40, length=30), cycles=7)
tf_trials = tf_morlet(epochs; frequencies=logrange(2, 40, length=20), cycles=(3, 10), return_trials=true, return_phase=true)
```
"""
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
    sdat, gpu_active = _tf_prepare(dat, pad, use_gpu, channel_selection, interval_selection, frequencies)
    n_samples_unpadded = n_samples(sdat)
    n_trials = n_epochs(sdat)
    times_out = time_vector(sdat)
    num_frex = length(frequencies)

    cycles_vec = cycles isa Tuple ? logrange(cycles[1], cycles[2], length = num_frex) : fill(cycles, num_frex)
    n_pre_pad, _, n_padded_samples = _tf_padding(pad, n_samples_unpadded)
    time_indices_out = (n_pre_pad+1):(n_pre_pad+n_samples_unpadded)

    max_sigma = maximum(cycles_vec ./ (2 * pi .* frequencies))
    max_hw = ceil(Int, 6 * max_sigma * sdat.sample_rate) ÷ 2
    n_conv_pow2 = nextpow(2, (max_hw * 2 + 1) + n_padded_samples - 1)

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

    wavelet_padded = Vector{ComplexF64}(undef, n_conv_pow2)
    p_fft_wavelet = plan_fft(wavelet_padded; flags = FFTW.MEASURE)
    inv_sr = 1.0 / sdat.sample_rate
    two_pi = 2 * pi

    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    wl_per_freq = Vector{Int}(undef, num_frex)
    sample_offsets = Vector{Int}(undef, num_frex)

    for fi = 1:num_frex
        freq_val = frequencies[fi]
        sigma = cycles_vec[fi] / (two_pi * freq_val)
        hw = ceil(Int, 6 * sigma * sdat.sample_rate) ÷ 2
        wl = hw * 2 + 1
        wl_per_freq[fi] = wl
        sample_offsets[fi] = hw

        two_pi_freq = two_pi * freq_val
        inv_2sigma2 = 1.0 / (2 * sigma^2)
        @inbounds @simd for i = 1:wl
            t = (i - 1 - hw) * inv_sr
            wavelet_padded[i] = cis(two_pi_freq * t) * exp(-(t * t) * inv_2sigma2)
        end
        @inbounds @simd for i = (wl+1):n_conv_pow2
            wavelet_padded[i] = zero(ComplexF64)
        end
        w_view = @view wavelet_padded[1:wl]
        w_view .*= sqrt(2.0) / norm(w_view, 2)

        wavelet_fft_freq = Vector{ComplexF64}(undef, n_conv_pow2)
        mul!(wavelet_fft_freq, p_fft_wavelet, wavelet_padded)
        wavelet_ffts[fi] = wavelet_fft_freq
    end

    fft_plan, ifft_plan = _tf_batch_fft_plans(n_conv_pow2, n_trials, gpu_active)
    gpu_buffers = gpu_active ? _tf_gpu_buffers(g) : nothing
    channel_names = channel_labels(sdat)
    out = TFOutput(times_out, frequencies, length(channel_names), g)

    _tf_foreach_channel(channel_names, gpu_active, host_outputs -> TFBuffers(g; host_outputs)) do ch_idx, channel, buf
        _tf_single_taper_channel!(
            ch_idx,
            channel,
            buf,
            sdat,
            g,
            time_indices_out,
            sample_offsets,
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

tf_morlet(data_vec::AbstractVector{<:EpochData}; kwargs...) = [tf_morlet(dat; kwargs...) for dat in data_vec]

"""
    tf_morlet(file_pattern::String;
              input_dir::String=pwd(),
              output_dir::Union{String,Nothing}=nothing,
              participant_selection=participants(),
              condition_selection=conditions(),
              kwargs...)

Batch Morlet wavelet time-frequency analysis on `EpochData` `.jld2` files.
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
