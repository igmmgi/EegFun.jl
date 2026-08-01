function tf_morlet(
    dat::EpochData;
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    frequencies::Union{AbstractRange,AbstractVector{<:Real}} = range(1, 40, length = 40),
    cycles::Union{Real,Tuple{Real,Real}} = 7,
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

    # Get original data time range (before padding) - these are the time points we want in output
    n_samples_original_unpadded = n_samples(dat)  # Store original unpadded length for edge filtering

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples(dat) # Use full signal for convolution (may be padded)

    times_out = time_vector(dat)

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
        n_post_pad = (pad == :post || pad == :both) ? n_samples_original_unpadded - 1 : 0
    else
        n_pre_pad = 0
        n_post_pad = 0
    end
    n_padded_samples = n_pre_pad + n_samples_original_unpadded + n_post_pad

    time_indices_out = (n_pre_pad+1):(n_pre_pad+n_samples_original_unpadded)
    n_times_out = length(time_indices_out)

    # Initialize output structures with only original time points
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
    n_conv_pow2 = nextpow(2, max_wl + n_padded_samples - 1)

    # GPU buffers
    if gpu_active
        template_padded_batch_gpu = gpu_array(zeros(ComplexF64, n_conv_pow2, n_trials))
        fft_plan_padded_batch = plan_fft!(template_padded_batch_gpu, 1)
        ifft_plan_padded_batch = plan_ifft!(template_padded_batch_gpu, 1)

        gpu_buffers = (
            local_data_padded_gpu = gpu_array(zeros(Float64, n_conv_pow2, n_trials)),
            local_data_fft_gpu = gpu_array(zeros(ComplexF64, n_conv_pow2, n_trials)),
            local_conv_result_gpu = gpu_array(zeros(ComplexF64, n_conv_pow2, n_trials)),
            curr_wavelet_gpu = gpu_array(zeros(ComplexF64, n_conv_pow2)),
            eegpower_trials_gpu = return_trials ? gpu_array(zeros(Float64, n_trials, num_frex, n_times_out)) : nothing,
            eegconv_trials_gpu = (return_trials && return_phase) ? gpu_array(zeros(ComplexF64, n_trials, num_frex, n_times_out)) : nothing,
            eegpower_avg_gpu = !return_trials ? gpu_array(zeros(Float64, num_frex, n_times_out)) : nothing,
            eegconv_avg_gpu = (!return_trials && return_phase) ? gpu_array(zeros(ComplexF64, num_frex, n_times_out)) : nothing,
            adjusted_indices_gpu = gpu_array(zeros(Int, n_times_out)),
        )
    else
        template_padded_batch = zeros(ComplexF64, n_conv_pow2, n_trials)
        fft_plan_padded_batch = plan_fft!(template_padded_batch, 1, flags = FFTW.MEASURE)
        ifft_plan_padded_batch = plan_ifft!(template_padded_batch, 1, flags = FFTW.MEASURE)
        gpu_buffers = nothing
    end

    # Pre-allocate wavelet buffer
    wavelet_padded = zeros(ComplexF64, n_conv_pow2)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)

    # Pre-compute some constants 
    inv_sr = 1.0 / dat.sample_rate
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    wl_per_freq = Vector{Int}(undef, num_frex)
    valid_starts_per_freq = Vector{Int}(undef, num_frex)

    for fi = 1:num_frex
        freq_val = frequencies[fi]
        sigma = cycles_vec[fi] / (two_pi * freq_val)
        hw = ceil(Int, 6 * sigma * dat.sample_rate) ÷ 2
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

    # A lock to safely write into the DataFrames
    df_lock = ReentrantLock()

    # Process each selected channel and each trial separately
    channel_names = channel_labels(dat)

    process_channel =
        (channel) -> begin
            _process_tf_channel!(
                channel,
                dat,
                n_trials,
                n_samples_per_epoch,
                n_conv_pow2,
                num_frex,
                return_trials,
                return_phase,
                filter_edges,
                pad,
                n_times_out,
                valid_starts_per_freq,
                time_indices_out,
                wl_per_freq,
                fft_plan_padded_batch,
                ifft_plan_padded_batch,
                wavelet_ffts,
                df_lock,
                power_df,
                phase_df,
                gpu_buffers,
            )
        end

    if gpu_active
        for channel in channel_names
            process_channel(channel)
        end
    else
        Threads.@threads for channel in channel_names
            process_channel(channel)
        end
    end

    # Restore original column order which might be shuffled by multithreading
    expected_cols = vcat([:time, :freq], channel_names)
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
        :wavelet,
        nothing,  # baseline
        copy(dat.analysis_info),
    )
end

function _process_tf_channel!(
    channel::Symbol,
    dat::EpochData,
    n_trials::Int,
    n_samples_per_epoch::Int,
    n_conv_pow2::Int,
    num_frex::Int,
    return_trials::Bool,
    return_phase::Bool,
    filter_edges::Bool,
    pad::Union{Nothing,Symbol},
    n_times_out::Int,
    valid_starts_per_freq::Vector{Int},
    time_indices_out,
    wl_per_freq::Vector{Int},
    fft_plan_padded_batch,
    ifft_plan_padded_batch,
    wavelet_ffts::Vector{Vector{ComplexF64}},
    df_lock::ReentrantLock,
    power_df,
    phase_df,
    gpu_buffers,
)
    # Thread-local reusable buffers (CPU)
    local_data_padded = Matrix{Float64}(undef, n_conv_pow2, n_trials)
    local_data_fft = Matrix{ComplexF64}(undef, n_conv_pow2, n_trials)
    local_conv_result = Matrix{ComplexF64}(undef, n_conv_pow2, n_trials)
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
        adjusted_indices_cpu = Vector{Int}(undef, n_times_out)

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
        # Pre-allocate reusable output buffers
        eegpower_trials = return_trials ? zeros(Float64, n_trials, num_frex, n_times_out) : Array{Float64,3}(undef, 0, 0, 0)
        eegconv_trials =
            (return_trials && return_phase) ? zeros(ComplexF64, n_trials, num_frex, n_times_out) : Array{ComplexF64,3}(undef, 0, 0, 0)
        eegpower_avg = return_trials ? Matrix{Float64}(undef, 0, 0) : zeros(Float64, num_frex, n_times_out)
        eegconv_avg = (!return_trials && return_phase) ? zeros(ComplexF64, num_frex, n_times_out) : Matrix{ComplexF64}(undef, 0, 0)
    end

    # Pre-extract all trial data for this channel
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

    # Pad data to n_conv_pow2 (zero-padding at the end)
    fill!(local_data_padded, 0.0)

    n_pre_pad = (!isnothing(pad) && (pad == :both || pad == :pre)) ? n_samples_per_epoch - 1 : 0
    n_post_pad = (!isnothing(pad) && (pad == :both || pad == :post)) ? n_samples_per_epoch - 1 : 0
    n_padded_samples = n_pre_pad + n_samples_per_epoch + n_post_pad

    # Explicit loop for padding to avoid copy allocations (Virtual Padding)
    @inbounds for trial_idx = 1:n_trials
        for j = 1:n_padded_samples
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
    sqrt2 = sqrt(2)
    for fi = 1:num_frex

        # Frequency-domain convolution: multiply data FFT by wavelet FFT (broadcast across all trials)
        # We use explicit loops to avoid .* broadcast allocations on huge matrices
        curr_wavelet = wavelet_ffts[fi]
        if gpu_active
            # GPU inner loop broadcast
            curr_wavelet_gpu = gpu_buffers.curr_wavelet_gpu
            copyto!(curr_wavelet_gpu, curr_wavelet)
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
        else
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
    end

    if gpu_active
        if !return_trials
            eegpower_avg_gpu ./= n_trials
            if return_phase
                eegconv_avg_gpu ./= n_trials
            end
        end

        # Download from GPU exactly once
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
    else
        if !return_trials
            eegpower_avg ./= n_trials
            if return_phase
                eegconv_avg ./= n_trials
            end
        end
    end

    if filter_edges
        # Use padded indices for edge filtering - padding extends valid region
        _filter_edges!(
            return_trials ? eegpower_trials : eegpower_avg,
            return_trials ? (return_phase ? eegconv_trials : nothing) : (return_phase ? eegconv_avg : nothing),
            num_frex,
            time_indices_out,
            wl_per_freq,
            n_padded_samples,
        )
    end

    lock(df_lock) do
        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = copy(vec(@view eegpower_trials[trial_idx, :, :]))

                if return_phase
                    phase_df[trial_idx][!, channel] = copy(vec(angle.(@view eegconv_trials[trial_idx, :, :])))
                else
                    phase_df[trial_idx][!, channel] = fill(NaN, num_frex * n_times_out)
                end
            end
        else
            power_df[!, channel] = copy(vec(eegpower_avg))

            if return_phase
                phase_vec = Vector{Float64}(undef, num_frex * n_times_out)
                @inbounds @simd for i in eachindex(eegconv_avg)
                    phase_vec[i] = angle(eegconv_avg[i])
                end
                phase_df[!, channel] = phase_vec
            else
                phase_df[!, channel] = fill(NaN, num_frex * n_times_out)
            end
        end
    end
    return nothing
end

# Vector version: process each EpochData element
function tf_morlet(data_vec::Vector{<:EpochData}; kwargs...)
    return [tf_morlet(dat; kwargs...) for dat in data_vec]
end

# =============================================================================
#     BATCH PROCESSING FUNCTIONS
# =============================================================================

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
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    kwargs...,
)
    log_file = "tf_morlet.log"
    setup_global_logging(log_file)

    try
        @info "Batch tf_morlet started at \$(now())"
        @log_call "tf_morlet"

        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_morlet_\$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '\$file_pattern'"
            return nothing
        end

        @info "Found \$(length(files)) files for tf_morlet analysis"

        process_fn =
            (input_path, output_path) -> begin
                filename = basename(input_path)
                data = read_data(input_path)
                if isnothing(data) || !(data isa Vector{<:EpochData})
                    return BatchResult(false, filename, "Invalid data type")
                end
                data = _condition_select(data, condition_selection)
                tf_results = tf_morlet(data; kwargs...)
                jldsave(output_path; data = tf_results)
                return BatchResult(true, filename, "TF morlet analysis complete (\$(length(tf_results)) conditions)")
            end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "TF morlet")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
