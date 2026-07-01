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
    template_padded_batch = zeros(ComplexF64, n_trials, n_conv_pow2)
    fft_plan_padded_batch = plan_fft!(template_padded_batch, 2, flags = FFTW.MEASURE)
    ifft_plan_padded_batch = plan_ifft!(template_padded_batch, 2, flags = FFTW.MEASURE)

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
    Threads.@threads for channel in channel_names
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
        )
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
)
    # Pre-allocate reusable output buffers (reused across trials for this channel)
    eegpower_trials = return_trials ? zeros(Float64, n_trials, num_frex, n_times_out) : Array{Float64,3}(undef, 0, 0, 0)
    eegconv_trials =
        (return_trials && return_phase) ? zeros(ComplexF64, n_trials, num_frex, n_times_out) : Array{ComplexF64,3}(undef, 0, 0, 0)
    eegpower_avg = return_trials ? Matrix{Float64}(undef, 0, 0) : zeros(Float64, num_frex, n_times_out)
    eegconv_avg = (!return_trials && return_phase) ? zeros(ComplexF64, num_frex, n_times_out) : Matrix{ComplexF64}(undef, 0, 0)

    # Thread-local reusable buffers
    local_data_padded = Matrix{Float64}(undef, n_trials, n_conv_pow2)
    local_data_fft = Matrix{ComplexF64}(undef, n_trials, n_conv_pow2)
    local_conv_result = Matrix{ComplexF64}(undef, n_trials, n_conv_pow2)
    local_trial_signals = Matrix{Float64}(undef, n_trials, n_samples_per_epoch)

    # Pre-extract all trial data for this channel
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

    # Pad data to n_conv_pow2 (zero-padding at the end)
    fill!(local_data_padded, 0.0)

    # Explicit loop for padding to avoid copy allocations
    @inbounds for j = 1:n_samples_per_epoch
        for i = 1:n_trials
            local_data_padded[i, j] = local_trial_signals[i, j]
        end
    end

    # FFT entire padded data (batch process all trials at once)
    @inbounds @simd for i in eachindex(local_data_padded)
        local_data_fft[i] = ComplexF64(local_data_padded[i])
    end
    fft_plan_padded_batch * local_data_fft

    # Process each frequency
    sqrt2 = sqrt(2)
    for fi = 1:num_frex

        # Frequency-domain convolution: multiply data FFT by wavelet FFT (broadcast across all trials)
        # We use explicit loops to avoid .* broadcast allocations on huge matrices
        curr_wavelet = wavelet_ffts[fi]
        @inbounds for i = 1:n_conv_pow2
            wv = curr_wavelet[i]
            @simd for trial_idx = 1:n_trials
                local_conv_result[trial_idx, i] = local_data_fft[trial_idx, i] * wv
            end
        end

        # IFFT to get time-domain result (in-place along dim 2)
        ifft_plan_padded_batch * local_conv_result

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
                    val = local_conv_result[trial_idx, idx]
                    eegpower_trials[trial_idx, fi, ti_out] = abs2(val)
                    if return_phase
                        eegconv_trials[trial_idx, fi, ti_out] = val
                    end
                end
            else
                power_sum = 0.0
                conv_sum = 0.0im
                @inbounds @simd for trial_idx = 1:n_trials
                    val = local_conv_result[trial_idx, idx]
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

    if filter_edges
        # Use padded indices for edge filtering - padding extends valid region
        _filter_edges!(
            return_trials ? eegpower_trials : eegpower_avg,
            return_trials ? (return_phase ? eegconv_trials : nothing) : (return_phase ? eegconv_avg : nothing),
            num_frex,
            time_indices_out,
            wl_per_freq,
            n_samples_per_epoch,
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
