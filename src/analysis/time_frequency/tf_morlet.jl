function tf_morlet(
    dat::EpochData;
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    frequencies::Union{AbstractRange,AbstractVector{<:Real}} = range(1, 40, length = 40),
    cycles::Union{Real,Tuple{Real,Real}} = 7,
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

    # Get original data time range (before padding) - these are the time points we want in output
    n_samples_original_unpadded = n_samples(dat)  # Store original unpadded length for edge filtering

    # Apply padding if requested 
    if !isnothing(pad)
        mirror!(dat, pad)
    end

    # Get number of trials/epochs
    n_trials = n_epochs(dat)
    n_samples_per_epoch = n_samples(dat) # Use full signal for convolution (may be padded)

    # Get all time points from processed (possibly padded) data
    times_all = time_vector(dat)

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
        time_indices_out = (n_pre_pad+1):(n_pre_pad+n_samples_original_unpadded)
    else
        time_indices_out = 1:n_samples_original_unpadded
    end
    n_times_out = length(time_indices_out)

    # Initialize output structures with only original time points as we unpad!
    times_out = times_all[time_indices_out]
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
    n_conv_pow2 = nextpow(2, max_wl + n_samples_per_epoch - 1)

    # Pre-allocate reusable buffers for single trial processing
    signal_padded = zeros(ComplexF64, n_conv_pow2)
    eegfft = zeros(ComplexF64, n_conv_pow2)
    wavelet_padded = zeros(ComplexF64, n_conv_pow2)
    conv_buffer = zeros(ComplexF64, n_conv_pow2)
    eegconv_buffer = zeros(ComplexF64, n_conv_pow2)

    # Create FFT plans using the actual buffers 
    p_fft_signal = plan_fft(signal_padded, flags = FFTW.MEASURE)
    p_fft_wavelet = plan_fft(wavelet_padded, flags = FFTW.MEASURE)
    p_ifft_conv = plan_ifft(eegconv_buffer, flags = FFTW.MEASURE)

    # Pre-compute some constants 
    inv_sr = 1.0 / dat.sample_rate
    two_pi = 2 * pi
    sqrt_pi = sqrt(pi)

    # Pre-compute wavelets and their FFTs once (same for all channels and trials)
    wavelet_ffts = Vector{Vector{ComplexF64}}(undef, num_frex)
    wl_per_freq = Vector{Int}(undef, num_frex)
    conv_indices_per_freq = Vector{Vector{Int}}(undef, num_frex)

    for fi = 1:num_frex
        freq_val = frequencies[fi]
        sigma = cycles_vec[fi] / (two_pi * freq_val)
        hw = ceil(Int, 6 * sigma * dat.sample_rate) ÷ 2
        wl = hw * 2 + 1
        wl_per_freq[fi] = wl
        valid_start = hw + 1

        # Pre-compute convolution indices for all samples in padded data
        conv_indices_per_freq[fi] = [valid_start + sample_idx - 1 for sample_idx = 1:n_samples_per_epoch]

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

    # Pre-allocate reusable output buffers (reused across all channels)
    # Only allocate for original unpadded samples
    if return_trials
        eegpower = zeros(Float64, num_frex, n_times_out, n_trials)
        eegconv = zeros(ComplexF64, num_frex, n_times_out, n_trials)
    else
        eegpower = zeros(Float64, num_frex, n_times_out)
        eegconv = zeros(ComplexF64, num_frex, n_times_out)
    end

    # TODO: initial testing shows this is just as fast as concatenating 
    # all trials and/or channels into a single matrix and then processing that
    # but it still seems too slow!!!

    # Process each selected channel and each trial separately
    for channel in channel_labels(dat)

        # Reset output arrays for this channel
        fill!(eegpower, 0)
        fill!(eegconv, 0)

        for trial_idx = 1:n_trials

            # Copy single trial data to padded buffer
            fill!(signal_padded, 0)
            signal_padded[1:n_samples_per_epoch] .= dat.data[trial_idx][!, channel]

            # FFT for this trial
            mul!(eegfft, p_fft_signal, signal_padded)

            # Loop through frequencies - reuse pre-computed wavelet FFTs
            for fi = 1:num_frex

                # TODO: here is the problem!!! But need to test again concat option(s)
                @inbounds @simd for i = 1:n_conv_pow2
                    conv_buffer[i] = wavelet_ffts[fi][i] * eegfft[i]
                end
                mul!(eegconv_buffer, p_ifft_conv, conv_buffer)

                # Extract only original unpadded samples from convolution
                conv_indices = conv_indices_per_freq[fi]
                @inbounds for (ti_out, ti_padded) in enumerate(time_indices_out)
                    val = eegconv_buffer[conv_indices[ti_padded]]

                    if return_trials
                        eegpower[fi, ti_out, trial_idx] = abs2(val)
                        eegconv[fi, ti_out, trial_idx] = val
                    else
                        eegpower[fi, ti_out] += abs2(val)
                        eegconv[fi, ti_out] += val
                    end

                end
            end
        end

        if !return_trials
            eegpower ./= n_trials
            eegconv ./= n_trials
        end

        if filter_edges
            # Use padded indices for edge filtering - padding extends valid region
            _filter_edges!(eegpower, eegconv, num_frex, time_indices_out, wl_per_freq, n_samples_per_epoch)
        end

        if return_trials # Store each trial separately
            for trial_idx = 1:n_trials
                power_df[trial_idx][!, channel] = copy(vec(@view eegpower[:, :, trial_idx]))
                phase_df[trial_idx][!, channel] = vec(angle.(@view eegconv[:, :, trial_idx]))
            end
        else
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
        :wavelet,
        nothing,  # baseline
        dat.analysis_info,
    )
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
- `file_pattern::String`: Pattern to match files (e.g., "epochs", "epochs_cleaned")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- All other keyword arguments are forwarded to `tf_morlet` (e.g., `frequencies`, `cycles`, `pad`, `filter_edges`)

# Examples
```julia
# Run Morlet wavelet TF on all epoch files
tf_morlet("epochs_cleaned"; cycles=7)

# Specific participants, custom frequencies
tf_morlet("epochs_cleaned";
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
        @info "Batch tf_morlet started at $(now())"
        @log_call "tf_morlet"

        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_morlet_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files for tf_morlet analysis"

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
                return BatchResult(true, filename, "TF morlet analysis complete ($(length(tf_results)) conditions)")
            end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "TF morlet")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
