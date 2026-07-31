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
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    frequencies::Union{AbstractRange,AbstractVector{<:Real}} = range(1, 40, length = 40),
    time_steps::Real = 0.05,
    cycles::Real,
    frequency_smoothing::Union{Nothing,Real} = nothing,
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

    # Validate cycles parameter (adaptive window only)
    if cycles <= 0
        error("`cycles` must be positive, got $cycles")
    end

    # Default frequency smoothing 
    if isnothing(frequency_smoothing)
        frequency_smoothing = 0.4
    end
    if frequency_smoothing <= 0
        error("`frequency_smoothing` must be positive, got $frequency_smoothing")
    end

    # Get time points and determine output time grid
    times = time_vector(dat)
    n_samples_per_epoch = n_samples(dat)
    time_min = minimum(times)
    time_max = maximum(times)
    time_steps_range = time_min:time_steps:time_max
    time_indices, times_out = find_times(times, time_steps_range)
    if isempty(time_indices)
        error("No valid time points found with step size $time_steps in range ($time_min to $time_max seconds)")
    end

    # Get number of trials/epochs
    n_trials = n_epochs(dat)

    # Use frequencies directly (convert to vector if needed for indexing)
    freqs = collect(frequencies)

    # Adaptive window mode: window length = cycles / frequency 
    cycles_vec = fill(Float64(cycles), num_frex)
    window_lengths_sec = cycles_vec ./ freqs  # Window length in seconds

    # Convert window lengths from seconds to samples (per frequency)
    n_window_samples_per_freq = [Int(round(wl * dat.sample_rate)) for wl in window_lengths_sec]
    if any(n -> n < 2, n_window_samples_per_freq)
        min_samples = minimum(n_window_samples_per_freq)
        error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end

    n_times = length(times_out)

    # Pre-compute DPSS tapers for each frequency
    tapers_per_freq = Vector{Matrix{Float64}}(undef, num_frex)  # Each element is (n_window_samples × n_tapers)
    n_tapers_per_freq = zeros(Int, num_frex)
    inv_n_tapers_per_freq = Vector{Float64}(undef, num_frex)

    # Pre-calculate padding amounts for virtual padding (no new dataframe allocations)
    n_pre_pad = (!isnothing(pad) && (pad == :both || pad == :pre)) ? n_samples_per_epoch - 1 : 0
    n_post_pad = (!isnothing(pad) && (pad == :both || pad == :post)) ? n_samples_per_epoch - 1 : 0
    n_padded_samples = n_pre_pad + n_samples_per_epoch + n_post_pad

    # Maximum window dictates if we need additional padding for FFT beyond our explicit pad
    max_window_samples = maximum(n_window_samples_per_freq)
    n_samples_padded = nextpow(2, max(n_padded_samples, max_window_samples))

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
            error(
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
            range(-(n_window_samples - 1) / 2, (n_window_samples - 1) / 2, length = n_window_samples) .* (2π * freq / dat.sample_rate)
        cos_wav = cos.(angle_in)
        sin_wav = sin.(angle_in)

        # For each taper, create tapered wavelet and compute FFT using pre-planned FFT
        for taper_idx = 1:K
            taper_col = @view tapers_per_freq[fi][:, taper_idx]
            # Tapered wavelet: taper * (cos + i*sin)
            tapered_wavelet = taper_col .* (cos_wav .+ im .* sin_wav)

            fill!(wavelet_padded, 0)
            wavelet_padded[1:n_window_samples] = tapered_wavelet
            w_fft = zeros(ComplexF64, n_samples_padded)
            mul!(w_fft, p_fft_wavelet, wavelet_padded)
            tapered_wavelet_ffts[fi][taper_idx] = w_fft
        end
    end

    if gpu_active
        template_padded_batch_gpu = gpu_array(zeros(ComplexF64, n_trials, n_samples_padded))
        fft_plan_padded_batch = plan_fft!(template_padded_batch_gpu, 2)
        ifft_plan_padded_batch = plan_ifft!(template_padded_batch_gpu, 2)
        
        gpu_buffers = (
            local_data_padded_gpu = gpu_array(zeros(Float64, n_trials, n_samples_padded)),
            local_data_fft_gpu = gpu_array(zeros(ComplexF64, n_trials, n_samples_padded)),
            local_conv_result_gpu = gpu_array(zeros(ComplexF64, n_trials, n_samples_padded)),
            curr_wavelet_gpu = gpu_array(zeros(ComplexF64, n_samples_padded)),
            eegpower_trials_gpu = return_trials ? gpu_array(zeros(Float64, n_trials, num_frex, n_times)) : nothing,
            eegconv_trials_gpu = (return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_trials, num_frex, n_times)) : nothing,
            eegpower_avg_gpu = !return_trials ? gpu_array(zeros(Float64, num_frex, n_times)) : nothing,
            eegconv_avg_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, num_frex, n_times)) : nothing,
            power_accum_trials_gpu = return_trials ? gpu_array(zeros(Float64, n_trials, n_times)) : nothing,
            complex_accum_trials_gpu = (return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_trials, n_times)) : nothing,
            power_accum_avg_gpu = !return_trials ? gpu_array(zeros(Float64, n_times)) : nothing,
            complex_accum_avg_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_times)) : nothing,
            adjusted_indices_gpu = gpu_array(zeros(Int, n_times)),
            extracted_power_gpu = !return_trials ? gpu_array(zeros(Float64, n_trials, n_times)) : nothing,
            sum_power_gpu = !return_trials ? gpu_array(zeros(Float64, 1, n_times)) : nothing,
            sum_complex_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, 1, n_times)) : nothing
        )
    else
        # Plan for batch FFT of entire padded data (all trials at once)
        template_padded_batch = zeros(ComplexF64, n_trials, n_samples_padded)
        fft_plan_padded_batch = plan_fft!(template_padded_batch, 2, flags = FFTW.MEASURE)
        ifft_plan_padded_batch = plan_ifft!(template_padded_batch, 2, flags = FFTW.MEASURE)
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

    # Process each selected channel
    function _process_multitaper_channel!(
        channel, dat, n_trials, n_samples_per_epoch, n_padded_samples, n_samples_padded, num_frex,
        return_trials, return_phase, filter_edges, pad, n_times, time_indices,
        n_window_samples_per_freq, n_tapers_per_freq, inv_n_tapers_per_freq,
        fft_plan_padded_batch, ifft_plan_padded_batch, tapered_wavelet_ffts, df_lock, power_df, phase_df, gpu_buffers
    )
        gpu_active = !isnothing(gpu_buffers)
        if gpu_active
            local_data_padded_gpu = gpu_buffers.local_data_padded_gpu
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
            
            if return_trials
                fill!(eegpower_trials_gpu, 0.0)
                return_phase && fill!(eegconv_trials_gpu, 0.0im)
            else
                fill!(eegpower_avg_gpu, 0.0)
                return_phase && fill!(eegconv_avg_gpu, 0.0im)
            end
        else
            local eegpower_trials::Array{Float64,3} = return_trials ? zeros(Float64, n_trials, num_frex, n_times) : Array{Float64,3}(undef, 0, 0, 0)
            local eegconv_trials::Array{ComplexF64,3} = (return_trials && return_phase) ? zeros(ComplexF64, n_trials, num_frex, n_times) : Array{ComplexF64,3}(undef, 0, 0, 0)
            local eegpower_avg::Matrix{Float64} = return_trials ? Matrix{Float64}(undef, 0, 0) : zeros(Float64, num_frex, n_times)
            local eegconv_avg::Matrix{ComplexF64} = (!return_trials && return_phase) ? zeros(ComplexF64, num_frex, n_times) : Matrix{ComplexF64}(undef, 0, 0)
            
            local power_accum_trials::Matrix{Float64} = return_trials ? Matrix{Float64}(undef, n_trials, n_times) : Matrix{Float64}(undef, 0, 0)
            local complex_accum_trials::Matrix{ComplexF64} = (return_trials && return_phase) ? Matrix{ComplexF64}(undef, n_trials, n_times) : Matrix{ComplexF64}(undef, 0, 0)
            local power_accum_avg::Vector{Float64} = return_trials ? Vector{Float64}(undef, 0) : Vector{Float64}(undef, n_times)
            local complex_accum_avg::Vector{ComplexF64} = (!return_trials && return_phase) ? Vector{ComplexF64}(undef, n_times) : Vector{ComplexF64}(undef, 0)
        end
        
        local_data_padded = Matrix{Float64}(undef, n_trials, n_samples_padded)
        local_data_fft = Matrix{ComplexF64}(undef, n_trials, n_samples_padded)
        local_conv_result = Matrix{ComplexF64}(undef, n_trials, n_samples_padded)
        local_trial_signals = Matrix{Float64}(undef, n_trials, n_samples_per_epoch)
        
        for trial_idx = 1:n_trials
            col = dat.data[trial_idx][!, channel]
            if col isa Vector{Float64}
                col_f64 = col::Vector{Float64}
                @inbounds @simd for i = 1:n_samples_per_epoch
                    local_trial_signals[trial_idx, i] = col_f64[i]
                end
            else
                @inbounds @simd for i = 1:n_samples_per_epoch
                    local_trial_signals[trial_idx, i] = Float64(col[i])
                end
            end
        end
        
        fill!(local_data_padded, 0.0)
        
        @inbounds for j = 1:n_padded_samples
            local src_j
            if j <= n_pre_pad
                src_j = n_samples_per_epoch - j + 1
            elseif j > n_pre_pad + n_samples_per_epoch
                idx_in_post = j - (n_pre_pad + n_samples_per_epoch)
                src_j = n_samples_per_epoch - idx_in_post
            else
                src_j = j - n_pre_pad
            end

            for i = 1:n_trials
                local_data_padded[i, j] = local_trial_signals[i, src_j]
            end
        end
        
        @inbounds @simd for i in eachindex(local_data_padded)
            local_data_fft[i] = ComplexF64(local_data_padded[i])
        end
        
        if gpu_active
            copyto!(local_data_fft_gpu, local_data_fft)
            fft_plan_padded_batch * local_data_fft_gpu
        else
            fft_plan_padded_batch * local_data_fft
        end
        

        for fi = 1:num_frex
            n_window_samples = n_window_samples_per_freq[fi]
            n_tapers = n_tapers_per_freq[fi]
            inv_n_tapers = inv_n_tapers_per_freq[fi]
            
            if gpu_active
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
            else
                if return_trials
                    fill!(power_accum_trials, 0.0)
                    return_phase && fill!(complex_accum_trials, 0.0im)
                else
                    fill!(power_accum_avg, 0.0)
                    return_phase && fill!(complex_accum_avg, 0.0im)
                end
            end
            
            for taper_idx = 1:n_tapers
                tapered_wavelet_fft = tapered_wavelet_ffts[fi][taper_idx]
                norm_factor = sqrt(2.0 / n_window_samples)
                half_window = n_window_samples ÷ 2
                
                if gpu_active
                    copyto!(gpu_buffers.curr_wavelet_gpu, tapered_wavelet_fft)
                    local_conv_result_gpu .= local_data_fft_gpu .* transpose(gpu_buffers.curr_wavelet_gpu)
                    ifft_plan_padded_batch * local_conv_result_gpu
                    local_conv_result_gpu .*= norm_factor
                    
                    extracted_gpu = @view local_conv_result_gpu[:, adjusted_indices_gpu]
                    
                    if return_trials
                        power_accum_trials_gpu .+= abs2.(extracted_gpu)
                        if return_phase
                            complex_accum_trials_gpu .+= extracted_gpu
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
                else
                    @inbounds for i = 1:n_samples_padded
                        wv = tapered_wavelet_fft[i]
                        @simd for trial_idx = 1:n_trials
                            local_conv_result[trial_idx, i] = local_data_fft[trial_idx, i] * wv
                        end
                    end
                    ifft_plan_padded_batch * local_conv_result
                    
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
                        conv_vals = @view local_conv_result[:, adjusted_idx]
                        
                        if return_trials
                            power_accum_trials[:, ti_idx] .+= abs2.(conv_vals)
                            if return_phase
                                complex_accum_trials[:, ti_idx] .+= conv_vals
                            end
                        else
                            power_accum_avg[ti_idx] += sum(abs2, conv_vals)
                            if return_phase
                                complex_accum_avg[ti_idx] += sum(conv_vals)
                            end
                        end
                    end
                end
            end
            
            if gpu_active
                if return_trials
                    @views eegpower_trials_gpu[:, fi, :] .= power_accum_trials_gpu .* inv_n_tapers
                    if return_phase
                        @views eegconv_trials_gpu[:, fi, :] .= complex_accum_trials_gpu .* inv_n_tapers
                    end
                else
                    @views eegpower_avg_gpu[fi, :] .= power_accum_avg_gpu .* inv_n_tapers
                    if return_phase
                        @views eegconv_avg_gpu[fi, :] .= complex_accum_avg_gpu .* inv_n_tapers
                    end
                end
            else
                if return_trials
                    @inbounds for ti_idx = 1:n_times
                        eegpower_trials[:, fi, ti_idx] .= power_accum_trials[:, ti_idx] .* inv_n_tapers
                        if return_phase
                            eegconv_trials[:, fi, ti_idx] .= complex_accum_trials[:, ti_idx] .* inv_n_tapers
                        end
                    end
                else
                    @inbounds for ti_idx = 1:n_times
                        eegpower_avg[fi, ti_idx] = power_accum_avg[ti_idx] * inv_n_tapers
                        if return_phase
                            eegconv_avg[fi, ti_idx] = complex_accum_avg[ti_idx] * inv_n_tapers
                        end
                    end
                end
            end
        end
        if gpu_active
            if return_trials
                eegpower_trials_cpu = Array(eegpower_trials_gpu)
                if return_phase
                    eegconv_trials_cpu = Array(eegconv_trials_gpu)
                end
            else
                eegpower_avg_cpu = Array(eegpower_avg_gpu)
                if return_phase
                    eegconv_avg_cpu = Array(eegconv_avg_gpu)
                end
            end
        else
            eegpower_trials_cpu = eegpower_trials
            eegconv_trials_cpu = eegconv_trials
            eegpower_avg_cpu = eegpower_avg
            eegconv_avg_cpu = eegconv_avg
        end

        if filter_edges
            # Compute exact window lengths in samples (floating point) for edge filtering
            window_lengths_samples_exact = [(cycles / freqs[fi]) * dat.sample_rate for fi = 1:num_frex]
            adjusted_time_indices = time_indices .+ n_pre_pad
            _filter_edges!(
                return_trials ? eegpower_trials_cpu : eegpower_avg_cpu,
                return_trials ? (return_phase ? eegconv_trials_cpu : nothing) : (return_phase ? eegconv_avg_cpu : nothing),
                num_frex,
                adjusted_time_indices,
                window_lengths_samples_exact,
                n_padded_samples,
            )
        end

        lock(df_lock) do
            if return_trials # Store each trial separately
                for trial_idx = 1:n_trials
                    power_df[trial_idx][!, channel] = copy(vec(@view eegpower_trials_cpu[trial_idx, :, :]))
                    if return_phase
                        phase_df[trial_idx][!, channel] = copy(vec(angle.(@view eegconv_trials_cpu[trial_idx, :, :])))
                    else
                        phase_df[trial_idx][!, channel] = fill(NaN, num_frex * n_times)
                    end
                end
            else
                power_df[!, channel] = copy(vec(eegpower_avg_cpu))
                if return_phase
                    phase_vec = Vector{Float64}(undef, num_frex * n_times)
                    @inbounds @simd for i in eachindex(eegconv_avg_cpu)
                        phase_vec[i] = angle(eegconv_avg_cpu[i])
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
            _process_multitaper_channel!(
                channel, dat, n_trials, n_samples_per_epoch, n_padded_samples, n_samples_padded, num_frex,
                return_trials, return_phase, filter_edges, pad, n_times, time_indices,
                n_window_samples_per_freq, n_tapers_per_freq, inv_n_tapers_per_freq,
                fft_plan_padded_batch, ifft_plan_padded_batch, tapered_wavelet_ffts, df_lock, power_df, phase_df, gpu_buffers
            )
        end
    else
        Threads.@threads for channel in selected_channels
            _process_multitaper_channel!(
                channel, dat, n_trials, n_samples_per_epoch, n_padded_samples, n_samples_padded, num_frex,
                return_trials, return_phase, filter_edges, pad, n_times, time_indices,
                n_window_samples_per_freq, n_tapers_per_freq, inv_n_tapers_per_freq,
                fft_plan_padded_batch, ifft_plan_padded_batch, tapered_wavelet_ffts, df_lock, power_df, phase_df, gpu_buffers
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
        :multitaper,
        nothing,  # baseline
        copy(dat.analysis_info),
    )
end

# Vector version: process each EpochData element
function tf_multitaper(data_vec::Vector{<:EpochData}; kwargs...)
    return [tf_multitaper(dat; kwargs...) for dat in data_vec]
end


# =============================================================================
#     BATCH PROCESSING FUNCTIONS
# =============================================================================
function tf_multitaper(
    file_pattern::String;
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    kwargs...,
)
    log_file = "tf_multitaper.log"
    setup_global_logging(log_file)

    try
        @info "Batch tf_multitaper started at $(now())"
        @log_call "tf_multitaper"

        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_multitaper_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files for tf_multitaper analysis"

        process_fn =
            (input_path, output_path) -> begin
                filename = basename(input_path)
                data = read_data(input_path)
                if isnothing(data) || !(data isa Vector{<:EpochData})
                    return BatchResult(false, filename, "Invalid data type")
                end
                data = _condition_select(data, condition_selection)
                tf_results = tf_multitaper(data; kwargs...)
                jldsave(output_path; data = tf_results)
                return BatchResult(true, filename, "TF multitaper analysis complete ($(length(tf_results)) conditions)")
            end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "TF multitaper")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
