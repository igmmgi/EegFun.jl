"""
    _tf_build_row_index(df) -> Dict{Tuple{Float64,Float64}, Int}

Build a lookup from (rounded_freq, rounded_time) → row index for a TF DataFrame.
Uses `round(x, digits=6)` for robust floating-point matching.
"""
function _tf_build_row_index(df::DataFrame)
    idx = Dict{Tuple{Float64,Float64},Int}()
    n = nrow(df)
    sizehint!(idx, n)
    freq_col = df.freq
    time_col = df.time
    for i = 1:n
        key = (round(freq_col[i], digits = 6), round(time_col[i], digits = 6))
        if !haskey(idx, key)
            idx[key] = i
        end
    end
    return idx
end

"""
    _tf_df_to_matrix!(mat, df, channel, rounded_freqs, rounded_times, row_index)

In-place extraction of a `[n_freqs × n_time]` power matrix for one channel from a TF DataFrame.
Writes into pre-allocated `mat` which should be initialized (e.g. with `fill!(mat, NaN)`).
"""
function _tf_df_to_matrix!(
    mat::AbstractMatrix{Float64},
    df::DataFrame,
    channel::Symbol,
    rounded_freqs::Vector{Float64},
    rounded_times::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int},
)
    col_data = df[!, channel]::Vector{Float64}
    @inbounds for (ti, rt) in enumerate(rounded_times)
        for (fi, rf) in enumerate(rounded_freqs)
            row = get(row_index, (rf, rt), nothing)
            if !isnothing(row)
                mat[fi, ti] = col_data[row]
            end
        end
    end
    return mat
end

"""
    _tf_df_to_matrix(df, channel, frequencies, time_points[, row_index]) -> Matrix{Float64}

Extract a `[n_freqs × n_time]` power matrix for one channel from a TF DataFrame.

Uses a Dict-based row index for O(1) lookups. Optionally accepts a precomputed `row_index`
(from `_tf_build_row_index`) to avoid rebuilding it when extracting multiple channels
from the same DataFrame.

Missing cells are filled with `NaN`.
"""
function _tf_df_to_matrix(
    df::DataFrame,
    channel::Symbol,
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int} = _tf_build_row_index(df);
    rounded_freqs::Vector{Float64} = round.(frequencies, digits = 6),
    rounded_times::Vector{Float64} = round.(time_points, digits = 6),
)
    mat = fill(NaN, length(frequencies), length(time_points))
    return _tf_df_to_matrix!(mat, df, channel, rounded_freqs, rounded_times, row_index)
end

"""
    _tf_matrix_to_df!(df, channel, mat, frequencies, time_points[, row_index])

Write a `[n_freqs × n_time]` matrix back into the DataFrame column for `channel`.

Uses the same Dict-based row index as `_tf_df_to_matrix` for consistency.
Only writes cells that exist in the DataFrame (skips missing rows).
"""
function _tf_matrix_to_df!(
    df::DataFrame,
    channel::Symbol,
    mat::Matrix{Float64},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int} = _tf_build_row_index(df),
)
    col = df[!, channel]::Vector{Float64}
    for (ti, t) in enumerate(time_points)
        rt = round(t, digits = 6)
        for (fi, f) in enumerate(frequencies)
            row = get(row_index, (round(f, digits = 6), rt), nothing)
            if !isnothing(row)
                col[row] = mat[fi, ti]
            end
        end
    end

    return nothing
end

"""
    _filter_edges!(eegpower, eegconv, num_frex, time_indices, window_lengths_per_freq, n_samples_per_epoch)

Set edge regions to NaN where the analysis window extends beyond the data boundaries.
Applied to both power and (optionally) complex convolution results.
"""
function _filter_edges!(
    eegpower::AbstractArray,
    eegconv::Union{AbstractArray,Nothing},
    num_frex::Int,
    time_indices::AbstractVector{Int},
    window_lengths_per_freq::Union{Vector{Int},Vector{Float64}},
    n_samples_per_epoch::Int,
)
    for fi = 1:num_frex

        half_nsamplefreqoi = window_lengths_per_freq[fi] / 2.0
        min_valid_threshold = half_nsamplefreqoi
        max_valid_threshold = n_samples_per_epoch - half_nsamplefreqoi

        for ti in eachindex(time_indices)
            sample_idx = time_indices[ti]
            if !(sample_idx >= min_valid_threshold && sample_idx < max_valid_threshold)
                if ndims(eegpower) == 3  # return_trials = true, layout is [trial, freq, time]
                    @views eegpower[:, fi, ti] .= NaN
                    if !isnothing(eegconv)
                        @views eegconv[:, fi, ti] .= NaN * im
                    end
                else  # return_trials = false
                    eegpower[fi, ti] = NaN
                    if !isnothing(eegconv)
                        eegconv[fi, ti] = NaN * im
                    end
                end
            end
        end
    end
end

"""
    _validate_tf_pad(pad)

Throw an error unless `pad` is `nothing`, `:pre`, `:post` or `:both`.
"""
function _validate_tf_pad(pad)
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        @minimal_error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end
    return nothing
end

"""
    _tf_gpu_active(use_gpu::Bool) -> Bool

Return whether the GPU should be used. Warns and returns `false` if a GPU was
requested but CUDA is not available or functional.
"""
function _tf_gpu_active(use_gpu::Bool)
    gpu_active = use_gpu && is_gpu_available()
    if use_gpu && !gpu_active
        @minimal_warning "GPU requested but CUDA is not available or functional. Falling back to CPU."
    end
    return gpu_active
end

"""
    _subset_tf_data(dat::EpochData, channel_selection, interval_selection) -> EpochData

Subset `dat` by channel and interval selection, erroring if no data remain.
"""
function _subset_tf_data(dat::EpochData, channel_selection, interval_selection::Interval)
    dat = subset(dat; channel_selection = channel_selection, interval_selection = interval_selection)
    isempty(dat.data) && @minimal_error("No data remaining after subsetting")
    return dat
end

"""
    _validate_tf_frequencies(frequencies)

Throw an error unless `frequencies` is non-empty and strictly positive.
"""
function _validate_tf_frequencies(frequencies)
    if length(frequencies) == 0
        @minimal_error("`frequencies` must contain at least one frequency")
    end
    if any(f -> f <= 0, frequencies)
        @minimal_error("All frequencies in `frequencies` must be positive")
    end
    return nothing
end

"""
    _validate_tf_window_samples(n_window_samples)

Throw an error if any window (given in samples, scalar or vector) is shorter than 2 samples.
"""
function _validate_tf_window_samples(n_window_samples)
    min_samples = minimum(n_window_samples)
    if min_samples < 2
        @minimal_error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end
    return nothing
end

"""
    _tf_time_grid(dat::EpochData, time_steps::Real) -> (time_indices, times_out)

Determine the output time grid (sample indices and time values) for a step size of
`time_steps` seconds. Errors if no valid time points are found.
"""
function _tf_time_grid(dat::EpochData, time_steps::Real)
    times = time_vector(dat)
    time_min = minimum(times)
    time_max = maximum(times)
    time_steps_range = time_min:time_steps:time_max
    time_indices, times_out = find_times(times, time_steps_range)
    if isempty(time_indices)
        @minimal_error("No valid time points found with step size $time_steps in range ($time_min to $time_max seconds)")
    end
    return time_indices, times_out
end

# =============================================================================
# SHARED DECOMPOSITION SCAFFOLDING (tf_morlet / tf_stft / tf_multitaper)
#
# All three methods follow the same pipeline per channel:
#   load trials -> mirror-pad -> batch FFT -> per-frequency kernel (method-specific)
#   -> average over trials -> filter edges -> store into power/phase DataFrames
# Everything except the per-frequency kernel lives here.
# =============================================================================

"""
    _tf_prepare(dat, pad, use_gpu, channel_selection, interval_selection, frequencies) -> (dat, gpu_active)

Validate common TF arguments and subset `dat` to the requested channels and interval.
"""
function _tf_prepare(dat::EpochData, pad, use_gpu::Bool, channel_selection, interval_selection::Interval, frequencies)
    _validate_tf_pad(pad)
    gpu_active = _tf_gpu_active(use_gpu)
    dat = _subset_tf_data(dat, channel_selection, interval_selection)
    _validate_tf_frequencies(frequencies)
    return dat, gpu_active
end

"""
    _tf_padding(pad, n_samples) -> (n_pre_pad, n_post_pad, n_padded)

Number of mirrored samples added before/after each epoch for `pad`, and the padded length.
"""
function _tf_padding(pad::Union{Nothing,Symbol}, n_samples::Int)
    n_pre_pad = (pad === :pre || pad === :both) ? n_samples - 1 : 0
    n_post_pad = (pad === :post || pad === :both) ? n_samples - 1 : 0
    return n_pre_pad, n_post_pad, n_pre_pad + n_samples + n_post_pad
end

"""
    TFGeometry

Sizes and options shared by the per-channel TF kernels.
"""
struct TFGeometry
    n_trials::Int
    n_samples::Int       # samples per epoch (unpadded)
    n_pre_pad::Int       # mirrored samples prepended to each epoch
    n_padded::Int        # n_pre_pad + n_samples + n_post_pad
    n_fft::Int           # FFT length (>= n_padded)
    num_frex::Int
    n_times::Int         # number of output time points
    return_trials::Bool
    return_phase::Bool
    filter_edges::Bool
end

"""
    _tf_batch_fft_plans(n_fft, n_trials, gpu_active) -> (fft_plan, ifft_plan)

In-place FFT plans along dimension 1 of an `n_fft × n_trials` complex matrix.
"""
function _tf_batch_fft_plans(n_fft::Int, n_trials::Int, gpu_active::Bool)
    if gpu_active
        template = gpu_array(zeros(ComplexF64, n_fft, n_trials))
        return plan_fft!(template, 1), plan_ifft!(template, 1)
    end
    template = zeros(ComplexF64, n_fft, n_trials)
    return plan_fft!(template, 1, flags = FFTW.MEASURE), plan_ifft!(template, 1, flags = FFTW.MEASURE)
end

"""
    _tf_padded_fft(wavelet, n_fft, plan, scratch) -> Vector{ComplexF64}

FFT of `wavelet` zero-padded to `n_fft` (uses `scratch` as the padded input buffer).
"""
function _tf_padded_fft(wavelet::AbstractVector, plan, scratch::Vector{ComplexF64})
    fill!(scratch, 0)
    scratch[1:length(wavelet)] = wavelet
    out = zeros(ComplexF64, length(scratch))
    mul!(out, plan, scratch)
    return out
end

"""
    TFBuffers

Reusable per-thread host buffers for one channel. Output arrays have size zero when not needed
(depending on `return_trials`/`return_phase`, or on the GPU where outputs live on the device).
"""
struct TFBuffers
    trial_signals::Matrix{Float64}   # n_samples × n_trials
    data_padded::Matrix{Float64}     # n_fft × n_trials (mirror-padded, zero-filled tail)
    data_fft::Matrix{ComplexF64}     # n_fft × n_trials
    conv_result::Matrix{ComplexF64}  # n_fft × n_trials
    power_trials::Array{Float64,3}   # n_trials × num_frex × n_times
    conv_trials::Array{ComplexF64,3}
    power_avg::Matrix{Float64}       # num_frex × n_times
    conv_avg::Matrix{ComplexF64}
end

function TFBuffers(g::TFGeometry; host_outputs::Bool = true)
    rt, rp = g.return_trials, g.return_phase
    return TFBuffers(
        Matrix{Float64}(undef, g.n_samples, g.n_trials),
        Matrix{Float64}(undef, g.n_fft, g.n_trials),
        Matrix{ComplexF64}(undef, g.n_fft, g.n_trials),
        host_outputs ? Matrix{ComplexF64}(undef, g.n_fft, g.n_trials) : Matrix{ComplexF64}(undef, 0, 0),
        (host_outputs && rt) ? zeros(Float64, g.n_trials, g.num_frex, g.n_times) : Array{Float64,3}(undef, 0, 0, 0),
        (host_outputs && rt && rp) ? zeros(ComplexF64, g.n_trials, g.num_frex, g.n_times) : Array{ComplexF64,3}(undef, 0, 0, 0),
        (host_outputs && !rt) ? zeros(Float64, g.num_frex, g.n_times) : Matrix{Float64}(undef, 0, 0),
        (host_outputs && !rt && rp) ? zeros(ComplexF64, g.num_frex, g.n_times) : Matrix{ComplexF64}(undef, 0, 0),
    )
end

"""
    _tf_gpu_buffers(g::TFGeometry; extra...) -> NamedTuple

Device buffers shared by all GPU kernels, merged with method-specific `extra` buffers.
"""
function _tf_gpu_buffers(g::TFGeometry; extra...)
    rt, rp = g.return_trials, g.return_phase
    base = (
        local_data_fft_gpu = gpu_array(zeros(ComplexF64, g.n_fft, g.n_trials)),
        local_conv_result_gpu = gpu_array(zeros(ComplexF64, g.n_fft, g.n_trials)),
        curr_wavelet_gpu = gpu_array(zeros(ComplexF64, g.n_fft)),
        eegpower_trials_gpu = rt ? gpu_array(zeros(Float64, g.n_trials, g.num_frex, g.n_times)) : nothing,
        eegconv_trials_gpu = (rt && rp) ? gpu_array(zeros(ComplexF64, g.n_trials, g.num_frex, g.n_times)) : nothing,
        eegpower_avg_gpu = !rt ? gpu_array(zeros(Float64, g.num_frex, g.n_times)) : nothing,
        eegconv_avg_gpu = (!rt && rp) ? gpu_array(zeros(ComplexF64, g.num_frex, g.n_times)) : nothing,
        adjusted_indices_gpu = gpu_array(zeros(Int, g.n_times)),
    )
    return merge(base, values(extra))
end

"""Zero the output arrays (host `TFBuffers` or device NamedTuple) before processing a channel."""
function _tf_reset_outputs!(buf::TFBuffers, g::TFGeometry)
    if g.return_trials
        fill!(buf.power_trials, 0.0)
        g.return_phase && fill!(buf.conv_trials, 0.0im)
    else
        fill!(buf.power_avg, 0.0)
        g.return_phase && fill!(buf.conv_avg, 0.0im)
    end
    return nothing
end
function _tf_reset_outputs!(gb::NamedTuple, g::TFGeometry)
    if g.return_trials
        fill!(gb.eegpower_trials_gpu, 0.0)
        g.return_phase && fill!(gb.eegconv_trials_gpu, 0.0im)
    else
        fill!(gb.eegpower_avg_gpu, 0.0)
        g.return_phase && fill!(gb.eegconv_avg_gpu, 0.0im)
    end
    return nothing
end

"""
    _tf_load_channel_fft!(buf, dat, channel, g, fft_plan, gpu_buffers)

Copy `channel` from every trial, mirror-pad according to `g`, and compute the batch FFT
(into `buf.data_fft`, or into `gpu_buffers.local_data_fft_gpu` on the GPU).
"""
function _tf_load_channel_fft!(buf::TFBuffers, dat::EpochData, channel::Symbol, g::TFGeometry, fft_plan, gpu_buffers)
    n_samples, n_trials, n_pre_pad, n_padded = g.n_samples, g.n_trials, g.n_pre_pad, g.n_padded
    trial_signals, data_padded, data_fft = buf.trial_signals, buf.data_padded, buf.data_fft

    # Pre-extract all trial data for this channel
    for trial_idx = 1:n_trials
        col = dat.data[trial_idx][!, channel]
        if col isa Vector{Float64}
            col_f64 = col::Vector{Float64}
            @inbounds @simd for i = 1:n_samples
                trial_signals[i, trial_idx] = col_f64[i]
            end
        else
            @inbounds @simd for i = 1:n_samples
                trial_signals[i, trial_idx] = Float64(col[i])
            end
        end
    end

    # Virtual mirror padding (no copies): [reversed pre | signal | reversed post | zeros]
    fill!(data_padded, 0.0)
    @inbounds for trial_idx = 1:n_trials
        for j = 1:n_padded
            local src_j
            if j <= n_pre_pad
                src_j = n_samples - j + 1
            elseif j > n_pre_pad + n_samples
                idx_in_post = j - (n_pre_pad + n_samples)
                src_j = n_samples - idx_in_post
            else
                src_j = j - n_pre_pad
            end
            data_padded[j, trial_idx] = trial_signals[src_j, trial_idx]
        end
    end

    @inbounds @simd for i in eachindex(data_padded)
        data_fft[i] = ComplexF64(data_padded[i])
    end

    # Batch FFT of all trials at once
    if isnothing(gpu_buffers)
        mul!(data_fft, fft_plan, data_fft)
    else
        copyto!(gpu_buffers.local_data_fft_gpu, data_fft)
        mul!(gpu_buffers.local_data_fft_gpu, fft_plan, gpu_buffers.local_data_fft_gpu)
    end
    return nothing
end

"""
    _tf_channel_results(buf, gpu_buffers, g) -> (power, conv)

The finished result arrays for the current channel (downloaded from the GPU if needed):
`[trial × freq × time]` when `return_trials`, otherwise `[freq × time]`. `conv` is `nothing`
unless `return_phase`.
"""
function _tf_channel_results(buf::TFBuffers, gpu_buffers, g::TFGeometry)
    rt, rp = g.return_trials, g.return_phase
    if isnothing(gpu_buffers)
        return rt ? (buf.power_trials, rp ? buf.conv_trials : nothing) : (buf.power_avg, rp ? buf.conv_avg : nothing)
    end
    gb = gpu_buffers
    return rt ? (Array(gb.eegpower_trials_gpu), rp ? Array(gb.eegconv_trials_gpu) : nothing) :
           (Array(gb.eegpower_avg_gpu), rp ? Array(gb.eegconv_avg_gpu) : nothing)
end

"""
    TFOutput

Output DataFrames (power and phase; one per trial when `return_trials`) plus a lock for
concurrent column writes from channel threads.
"""
struct TFOutput{T<:Union{DataFrame,Vector{DataFrame}}}
    power::T
    phase::T
    lock::ReentrantLock
end

function TFOutput(times_out::AbstractVector, freqs::AbstractVector, g::TFGeometry)
    time_col = repeat(times_out, inner = length(freqs))
    freq_col = repeat(freqs, outer = length(times_out))
    new_df() = DataFrame(time = time_col, freq = freq_col, copycols = false)
    if g.return_trials
        return TFOutput([new_df() for _ = 1:g.n_trials], [new_df() for _ = 1:g.n_trials], ReentrantLock())
    end
    return TFOutput(new_df(), new_df(), ReentrantLock())
end

"""
    _tf_store_channel!(out::TFOutput, channel, power, conv, g)

Write one channel's power (and phase, or `NaN` if not requested) into the output DataFrames.
"""
function _tf_store_channel!(out::TFOutput, channel::Symbol, power::AbstractArray, conv, g::TFGeometry)
    n_cells = g.num_frex * g.n_times
    lock(out.lock)
    try
        if g.return_trials
            for trial_idx = 1:g.n_trials
                out.power[trial_idx][!, channel] = copy(vec(@view power[trial_idx, :, :]))
                out.phase[trial_idx][!, channel] = g.return_phase ? vec(angle.(@view conv[trial_idx, :, :])) : fill(NaN, n_cells)
            end
        else
            out.power[!, channel] = copy(vec(power))
            if g.return_phase
                phase_vec = Vector{Float64}(undef, n_cells)
                @inbounds @simd for i in eachindex(conv)
                    phase_vec[i] = angle(conv[i])
                end
                out.phase[!, channel] = phase_vec
            else
                out.phase[!, channel] = fill(NaN, n_cells)
            end
        end
    finally
        unlock(out.lock)
    end
    return nothing
end

"""
    _tf_foreach_channel(process!, channels, gpu_active, make_buffers)

Call `process!(channel, buffers)` for every channel. On the CPU, channels are split across
threads and each thread allocates its own buffers once via `make_buffers(true)`; on the GPU,
channels run sequentially with `make_buffers(false)` (outputs live on the device).
"""
function _tf_foreach_channel(process!::F, channels::Vector{Symbol}, gpu_active::Bool, make_buffers::B) where {F,B}
    gpu_active && return _tf_foreach_channel_gpu(process!, channels, make_buffers)
    n_threads = min(Threads.nthreads(), length(channels))
    channel_chunks = [Symbol[] for _ = 1:n_threads]
    for (i, ch) in enumerate(channels)
        push!(channel_chunks[mod1(i, n_threads)], ch)
    end
    # NOTE: no variable assigned inside this closure may also be assigned in the enclosing
    # function, otherwise Julia boxes it into a single binding shared by all threads.
    Threads.@threads for tid = 1:n_threads
        local thread_buffers = make_buffers(true)
        for channel in channel_chunks[tid]
            process!(channel, thread_buffers)
        end
    end
    return nothing
end

function _tf_foreach_channel_gpu(process!::F, channels::Vector{Symbol}, make_buffers::B) where {F,B}
    device_buffers = make_buffers(false)
    for channel in channels
        process!(channel, device_buffers)
    end
    return nothing
end

"""
    _tf_result(dat, out::TFOutput, channels, method, g) -> TimeFreqData or TimeFreqEpochData

Restore channel column order (threads may write columns in any order) and wrap the output.
"""
function _tf_result(dat::EpochData, out::TFOutput, channels::Vector{Symbol}, method::Symbol, g::TFGeometry)
    expected_cols = vcat([:time, :freq], channels)
    _foreach_df(df -> select!(df, expected_cols), out.power)
    _foreach_df(df -> select!(df, expected_cols), out.phase)
    return_type = g.return_trials ? TimeFreqEpochData : TimeFreqData
    return return_type(
        dat.file,
        dat.condition,
        dat.condition_name,
        out.power,
        out.phase,
        copy(dat.layout),
        dat.sample_rate,
        method,
        nothing,  # baseline
        copy(dat.analysis_info),
    )
end

"""
    _tf_batch_worker(tf_fn, label, condition_selection, kwargs) -> Function

Per-file worker for the TF batch methods: load `Vector{EpochData}`, select conditions,
run `tf_fn(data; kwargs...)` and save the results.
"""
function _tf_batch_worker(tf_fn::Function, label::String, condition_selection, kwargs)
    return (input_path, output_path) -> begin
        filename = basename(input_path)
        data = read_data(input_path)
        if isnothing(data) || !(data isa Vector{<:EpochData})
            return BatchResult(false, filename, "Invalid data type")
        end
        data = _condition_select(data, condition_selection)
        tf_results = tf_fn(data; kwargs...)
        jldsave(output_path; data = tf_results)
        return BatchResult(true, filename, "TF $label analysis complete ($(length(tf_results)) conditions)")
    end
end
