function _tf_build_row_index(df::DataFrame)
    idx = Dict{Tuple{Float64,Float64},Int}()
    n = nrow(df)
    sizehint!(idx, n)
    freq_col, time_col = df.freq, df.time
    @inbounds for i = 1:n
        key = (round(freq_col[i], digits = 6), round(time_col[i], digits = 6))
        haskey(idx, key) || (idx[key] = i)
    end
    return idx
end

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
            !isnothing(row) && (mat[fi, ti] = col_data[row])
        end
    end
    return mat
end

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

function _tf_matrix_to_df!(
    df::DataFrame,
    channel::Symbol,
    mat::Matrix{Float64},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int} = _tf_build_row_index(df),
)
    col = df[!, channel]::Vector{Float64}
    @inbounds for (ti, t) in enumerate(time_points)
        rt = round(t, digits = 6)
        for (fi, f) in enumerate(frequencies)
            row = get(row_index, (round(f, digits = 6), rt), nothing)
            !isnothing(row) && (col[row] = mat[fi, ti])
        end
    end
    return nothing
end

function _filter_edges!(
    eegpower::AbstractArray,
    eegconv::Union{AbstractArray,Nothing},
    num_frex::Int,
    time_indices::AbstractVector{Int},
    window_lengths_per_freq::Union{Vector{Int},Vector{Float64}},
    n_samples_per_epoch::Int,
)
    nan_c = NaN * im
    is_3d = ndims(eegpower) == 3
    n_tr = is_3d ? size(eegpower, 1) : 0

    @inbounds for fi = 1:num_frex
        half_win = (window_lengths_per_freq[fi] - 1.0) * 0.5
        min_valid = 1.0 + half_win
        max_valid = n_samples_per_epoch - half_win

        for ti in eachindex(time_indices)
            s = time_indices[ti]
            if !(min_valid <= s <= max_valid)
                if is_3d
                    @simd for tr = 1:n_tr
                        eegpower[tr, fi, ti] = NaN
                    end
                    if !isnothing(eegconv)
                        @simd for tr = 1:n_tr
                            eegconv[tr, fi, ti] = nan_c
                        end
                    end
                else
                    eegpower[fi, ti] = NaN
                    !isnothing(eegconv) && (eegconv[fi, ti] = nan_c)
                end
            end
        end
    end
    return nothing
end

function _validate_tf_pad(pad)
    (!isnothing(pad) && pad ∉ (:pre, :post, :both)) && @minimal_error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    return nothing
end

function _tf_gpu_active(use_gpu::Bool)
    gpu_active = use_gpu && is_gpu_available()
    use_gpu && !gpu_active && @minimal_warning "GPU requested but CUDA is not available or functional. Falling back to CPU."
    return gpu_active
end

function _subset_tf_data(dat::EpochData, channel_selection, interval_selection::Interval)
    sdat = subset(dat; channel_selection = channel_selection, interval_selection = interval_selection)
    isempty(sdat.data) && @minimal_error("No data remaining after subsetting")
    return sdat
end

function _validate_tf_frequencies(frequencies)
    isempty(frequencies) && @minimal_error("`frequencies` must contain at least one frequency")
    any(<=(0), frequencies) && @minimal_error("All frequencies in `frequencies` must be positive")
    return nothing
end

function _validate_tf_window_samples(n_window_samples)
    min_samples = minimum(n_window_samples)
    min_samples < 2 && @minimal_error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    return nothing
end

function _tf_time_grid(dat::EpochData, time_steps::Real)
    times = time_vector(dat)
    t_min, t_max = extrema(times)
    time_indices, times_out = find_times(times, t_min:time_steps:t_max)
    isempty(time_indices) && @minimal_error("No valid time points found with step size $time_steps in range ($t_min to $t_max seconds)")
    return time_indices, times_out
end

# =============================================================================
# SHARED DECOMPOSITION SCAFFOLDING (tf_morlet / tf_stft / tf_multitaper)
# =============================================================================

function _tf_prepare(dat::EpochData, pad, use_gpu::Bool, channel_selection, interval_selection::Interval, frequencies)
    _validate_tf_pad(pad)
    gpu_active = _tf_gpu_active(use_gpu)
    sdat = _subset_tf_data(dat, channel_selection, interval_selection)
    _validate_tf_frequencies(frequencies)
    return sdat, gpu_active
end

function _tf_padding(pad::Union{Nothing,Symbol}, n_samples::Int)
    n_pre_pad = (pad === :pre || pad === :both) ? n_samples - 1 : 0
    n_post_pad = (pad === :post || pad === :both) ? n_samples - 1 : 0
    return n_pre_pad, n_post_pad, n_pre_pad + n_samples + n_post_pad
end

struct TFGeometry
    n_trials::Int
    n_samples::Int
    n_pre_pad::Int
    n_padded::Int
    n_fft::Int
    num_frex::Int
    n_times::Int
    return_trials::Bool
    return_phase::Bool
    filter_edges::Bool
end

function _tf_batch_fft_plans(n_fft::Int, n_trials::Int, gpu_active::Bool)
    if gpu_active
        template = gpu_array(zeros(ComplexF64, n_fft, n_trials))
        return plan_fft!(template, 1), plan_ifft!(template, 1)
    end
    template = Matrix{ComplexF64}(undef, n_fft, n_trials)
    return plan_fft!(template, 1; flags = FFTW.MEASURE), plan_ifft!(template, 1; flags = FFTW.MEASURE)
end

function _tf_padded_fft(wavelet::AbstractVector, plan, scratch::Vector{ComplexF64})
    n_w = length(wavelet)
    @inbounds for i = 1:n_w
        scratch[i] = wavelet[i]
    end
    @inbounds @simd for i = (n_w+1):length(scratch)
        scratch[i] = zero(ComplexF64)
    end
    out = Vector{ComplexF64}(undef, length(scratch))
    mul!(out, plan, scratch)
    return out
end

struct TFBuffers
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
        Matrix{ComplexF64}(undef, g.n_fft, g.n_trials),
        host_outputs ? Matrix{ComplexF64}(undef, g.n_fft, g.n_trials) : Matrix{ComplexF64}(undef, 0, 0),
        (host_outputs && rt) ? Array{Float64,3}(undef, g.n_trials, g.num_frex, g.n_times) : Array{Float64,3}(undef, 0, 0, 0),
        (host_outputs && rt && rp) ? Array{ComplexF64,3}(undef, g.n_trials, g.num_frex, g.n_times) : Array{ComplexF64,3}(undef, 0, 0, 0),
        (host_outputs && !rt) ? Matrix{Float64}(undef, g.num_frex, g.n_times) : Matrix{Float64}(undef, 0, 0),
        (host_outputs && !rt && rp) ? Matrix{ComplexF64}(undef, g.num_frex, g.n_times) : Matrix{ComplexF64}(undef, 0, 0),
    )
end

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
        extracted_power_gpu = !rt ? gpu_array(zeros(Float64, g.n_times, g.n_trials)) : nothing,
        sum_power_gpu = !rt ? gpu_array(zeros(Float64, g.n_times, 1)) : nothing,
        sum_complex_gpu = (!rt && rp) ? gpu_array(zeros(ComplexF64, g.n_times, 1)) : nothing,
    )
    return merge(base, values(extra))
end

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

# Concrete function barrier for mirror-padded direct ComplexF64 extraction
@inline function _load_mirror_padded_trial!(
    data_fft::Matrix{ComplexF64},
    trial_idx::Int,
    col::AbstractVector,
    n_samples::Int,
    n_pre_pad::Int,
    n_padded::Int,
    n_fft::Int,
)
    if n_pre_pad == 0 && n_padded == n_samples
        @inbounds @simd for j = 1:n_samples
            data_fft[j, trial_idx] = ComplexF64(Float64(col[j]), 0.0)
        end
    else
        @inbounds @simd for j = 1:n_pre_pad
            data_fft[j, trial_idx] = ComplexF64(Float64(col[n_samples-j+1]), 0.0)
        end
        @inbounds @simd for j = 1:n_samples
            data_fft[n_pre_pad+j, trial_idx] = ComplexF64(Float64(col[j]), 0.0)
        end
        n_post = n_padded - n_pre_pad - n_samples
        @inbounds @simd for j = 1:n_post
            data_fft[n_pre_pad+n_samples+j, trial_idx] = ComplexF64(Float64(col[n_samples-j]), 0.0)
        end
    end
    @inbounds @simd for j = (n_padded+1):n_fft
        data_fft[j, trial_idx] = zero(ComplexF64)
    end
    return nothing
end

function _tf_load_channel_fft!(buf::TFBuffers, dat::EpochData, channel::Symbol, g::TFGeometry, fft_plan, gpu_buffers)
    data_fft = buf.data_fft
    for trial_idx = 1:g.n_trials
        _load_mirror_padded_trial!(data_fft, trial_idx, dat.data[trial_idx][!, channel], g.n_samples, g.n_pre_pad, g.n_padded, g.n_fft)
    end
    if isnothing(gpu_buffers)
        mul!(data_fft, fft_plan, data_fft)
    else
        copyto!(gpu_buffers.local_data_fft_gpu, data_fft)
        mul!(gpu_buffers.local_data_fft_gpu, fft_plan, gpu_buffers.local_data_fft_gpu)
    end
    return nothing
end

function _tf_channel_results(buf::TFBuffers, gpu_buffers, g::TFGeometry)
    rt, rp = g.return_trials, g.return_phase
    if isnothing(gpu_buffers)
        return rt ? (buf.power_trials, rp ? buf.conv_trials : nothing) : (buf.power_avg, rp ? buf.conv_avg : nothing)
    end
    gb = gpu_buffers
    return rt ? (Array(gb.eegpower_trials_gpu), rp ? Array(gb.eegconv_trials_gpu) : nothing) :
           (Array(gb.eegpower_avg_gpu), rp ? Array(gb.eegconv_avg_gpu) : nothing)
end

struct TFOutput{TFreq}
    time_col::Vector{Float64}
    freq_col::Vector{TFreq}
    power_avg::Vector{Vector{Float64}}
    phase_avg::Vector{Vector{Float64}}
    power_trials::Matrix{Vector{Float64}}
    phase_trials::Matrix{Vector{Float64}}
end

function TFOutput(times_out::AbstractVector, freqs::AbstractVector, n_ch::Int, g::TFGeometry)
    time_col = repeat(Float64.(times_out), inner = length(freqs))
    freq_col = repeat(freqs, outer = length(times_out))
    if g.return_trials
        return TFOutput(
            time_col,
            freq_col,
            Vector{Float64}[],
            Vector{Float64}[],
            Matrix{Vector{Float64}}(undef, n_ch, g.n_trials),
            Matrix{Vector{Float64}}(undef, n_ch, g.n_trials),
        )
    end
    return TFOutput(
        time_col,
        freq_col,
        Vector{Vector{Float64}}(undef, n_ch),
        Vector{Vector{Float64}}(undef, n_ch),
        Matrix{Vector{Float64}}(undef, 0, 0),
        Matrix{Vector{Float64}}(undef, 0, 0),
    )
end

function _tf_store_channel!(out::TFOutput, ch_idx::Int, power::AbstractArray, conv, g::TFGeometry)
    n_cells = g.num_frex * g.n_times
    if g.return_trials
        for tr = 1:g.n_trials
            p_col = Vector{Float64}(undef, n_cells)
            idx = 0
            @inbounds for ti = 1:g.n_times, fi = 1:g.num_frex
                idx += 1
                p_col[idx] = power[tr, fi, ti]
            end
            out.power_trials[ch_idx, tr] = p_col
            if g.return_phase
                ph_col = Vector{Float64}(undef, n_cells)
                idx = 0
                @inbounds for ti = 1:g.n_times, fi = 1:g.num_frex
                    idx += 1
                    ph_col[idx] = angle(conv[tr, fi, ti])
                end
                out.phase_trials[ch_idx, tr] = ph_col
            else
                out.phase_trials[ch_idx, tr] = fill(NaN, n_cells)
            end
        end
    else
        out.power_avg[ch_idx] = copy(vec(power))
        if g.return_phase
            ph_col = Vector{Float64}(undef, n_cells)
            @inbounds @simd for i in eachindex(conv)
                ph_col[i] = angle(conv[i])
            end
            out.phase_avg[ch_idx] = ph_col
        else
            out.phase_avg[ch_idx] = fill(NaN, n_cells)
        end
    end
    return nothing
end

function _tf_run_chunks!(process!::F, channels::Vector{Symbol}, chunks::Vector{UnitRange{Int}}, thread_buffers::Vector{TB}) where {F,TB}
    Threads.@threads :static for tid in eachindex(chunks)
        buf = thread_buffers[tid]
        for ci in chunks[tid]
            process!(ci, channels[ci], buf)
        end
    end
    return nothing
end

function _tf_foreach_channel(process!::F, channels::Vector{Symbol}, gpu_active::Bool, make_buffers::B) where {F,B}
    n_ch = length(channels)
    n_ch == 0 && return nothing
    if gpu_active
        device_buffers = make_buffers(false)
        for ci = 1:n_ch
            process!(ci, channels[ci], device_buffers)
        end
        return nothing
    end
    n_chunks = min(Threads.nthreads(), n_ch)
    if n_chunks == 1
        buf = make_buffers(true)
        for ci = 1:n_ch
            process!(ci, channels[ci], buf)
        end
        return nothing
    end
    chunk_size = cld(n_ch, n_chunks)
    chunks = [((k-1)*chunk_size+1):min(k*chunk_size, n_ch) for k = 1:n_chunks if (k - 1) * chunk_size < n_ch]
    thread_buffers = [make_buffers(true) for _ in eachindex(chunks)]
    _tf_run_chunks!(process!, channels, chunks, thread_buffers)
    return nothing
end

function _build_tf_df(
    time_col::Vector{Float64},
    freq_col::AbstractVector,
    ch_cols::AbstractVector{Vector{Float64}},
    col_names::Vector{Symbol},
)
    n_ch = length(ch_cols)
    cols = Vector{Any}(undef, n_ch + 2)
    cols[1] = time_col
    cols[2] = freq_col
    @inbounds for ci = 1:n_ch
        cols[ci+2] = ch_cols[ci]
    end
    return DataFrame(cols, col_names; copycols = false)
end

function _tf_result(dat::EpochData, out::TFOutput, channels::Vector{Symbol}, method::Symbol, g::TFGeometry)
    col_names = vcat([:time, :freq], channels)
    if g.return_trials
        pow_dfs = [_build_tf_df(out.time_col, out.freq_col, @view(out.power_trials[:, tr]), col_names) for tr = 1:g.n_trials]
        phs_dfs = [_build_tf_df(out.time_col, out.freq_col, @view(out.phase_trials[:, tr]), col_names) for tr = 1:g.n_trials]
        return TimeFreqEpochData(
            dat.file,
            dat.condition,
            dat.condition_name,
            pow_dfs,
            phs_dfs,
            copy(dat.layout),
            dat.sample_rate,
            method,
            nothing,
            copy(dat.analysis_info),
        )
    end
    pow_df = _build_tf_df(out.time_col, out.freq_col, out.power_avg, col_names)
    phs_df = _build_tf_df(out.time_col, out.freq_col, out.phase_avg, col_names)
    return TimeFreqData(
        dat.file,
        dat.condition,
        dat.condition_name,
        pow_df,
        phs_df,
        copy(dat.layout),
        dat.sample_rate,
        method,
        nothing,
        copy(dat.analysis_info),
    )
end

# Shared single-taper frequency convolution + time extraction kernel (used by tf_morlet and tf_stft)
function _tf_single_taper_channel!(
    ch_idx::Int,
    channel::Symbol,
    buf::TFBuffers,
    dat::EpochData,
    g::TFGeometry,
    time_indices::AbstractVector{Int},
    sample_offsets::Vector{Int},
    edge_time_indices::AbstractVector{Int},
    edge_window_lengths::Union{Vector{Int},Vector{Float64}},
    fft_plan::P1,
    ifft_plan::P2,
    wavelet_ffts::Vector{Vector{ComplexF64}},
    out::TFOutput,
    gpu_buffers,
) where {P1,P2}
    n_trials, n_fft, num_frex, n_times = g.n_trials, g.n_fft, g.num_frex, g.n_times
    return_trials, return_phase = g.return_trials, g.return_phase
    inv_n_trials = 1.0 / n_trials

    _tf_load_channel_fft!(buf, dat, channel, g, fft_plan, gpu_buffers)

    if !isnothing(gpu_buffers)
        _tf_reset_outputs!(gpu_buffers, g)
        gb = gpu_buffers
        adjusted_indices_cpu = Vector{Int}(undef, n_times)

        for fi = 1:num_frex
            copyto!(gb.curr_wavelet_gpu, wavelet_ffts[fi])
            gb.local_conv_result_gpu .= gb.local_data_fft_gpu .* gb.curr_wavelet_gpu
            mul!(gb.local_conv_result_gpu, ifft_plan, gb.local_conv_result_gpu)

            off = sample_offsets[fi]
            @inbounds for ti = 1:n_times
                adjusted_indices_cpu[ti] = time_indices[ti] + off
            end
            copyto!(gb.adjusted_indices_gpu, adjusted_indices_cpu)
            extracted_gpu = @view gb.local_conv_result_gpu[gb.adjusted_indices_gpu, :]

            if return_trials
                @views gb.eegpower_trials_gpu[:, fi, :] .= transpose(abs2.(extracted_gpu))
                return_phase && (@views gb.eegconv_trials_gpu[:, fi, :] .= transpose(extracted_gpu))
            else
                gb.extracted_power_gpu .= abs2.(extracted_gpu)
                sum!(gb.sum_power_gpu, gb.extracted_power_gpu)
                @views gb.eegpower_avg_gpu[fi, :] .= vec(gb.sum_power_gpu)
                if return_phase
                    sum!(gb.sum_complex_gpu, extracted_gpu)
                    @views gb.eegconv_avg_gpu[fi, :] .= vec(gb.sum_complex_gpu)
                end
            end
        end

        if !return_trials
            gb.eegpower_avg_gpu .*= inv_n_trials
            return_phase && (gb.eegconv_avg_gpu .*= inv_n_trials)
        end
    else
        local_data_fft = buf.data_fft
        local_conv_result = buf.conv_result
        eegpower_trials, eegconv_trials = buf.power_trials, buf.conv_trials
        eegpower_avg, eegconv_avg = buf.power_avg, buf.conv_avg

        for fi = 1:num_frex
            curr_wavelet = wavelet_ffts[fi]
            @inbounds for tr = 1:n_trials
                @simd for i = 1:n_fft
                    local_conv_result[i, tr] = local_data_fft[i, tr] * curr_wavelet[i]
                end
            end
            mul!(local_conv_result, ifft_plan, local_conv_result)

            off = sample_offsets[fi]
            if return_trials
                if return_phase
                    @inbounds for ti = 1:n_times
                        idx = time_indices[ti] + off
                        @simd for tr = 1:n_trials
                            val = local_conv_result[idx, tr]
                            eegpower_trials[tr, fi, ti] = abs2(val)
                            eegconv_trials[tr, fi, ti] = val
                        end
                    end
                else
                    @inbounds for ti = 1:n_times
                        idx = time_indices[ti] + off
                        @simd for tr = 1:n_trials
                            eegpower_trials[tr, fi, ti] = abs2(local_conv_result[idx, tr])
                        end
                    end
                end
            elseif return_phase
                @inbounds for ti = 1:n_times
                    idx = time_indices[ti] + off
                    p_sum = 0.0
                    c_sum = 0.0im
                    @simd for tr = 1:n_trials
                        val = local_conv_result[idx, tr]
                        p_sum += abs2(val)
                        c_sum += val
                    end
                    eegpower_avg[fi, ti] = p_sum * inv_n_trials
                    eegconv_avg[fi, ti] = c_sum * inv_n_trials
                end
            else
                @inbounds for ti = 1:n_times
                    idx = time_indices[ti] + off
                    p_sum = 0.0
                    @simd for tr = 1:n_trials
                        p_sum += abs2(local_conv_result[idx, tr])
                    end
                    eegpower_avg[fi, ti] = p_sum * inv_n_trials
                end
            end
        end
    end

    power, conv = _tf_channel_results(buf, gpu_buffers, g)
    g.filter_edges && _filter_edges!(power, conv, num_frex, edge_time_indices, edge_window_lengths, g.n_padded)
    _tf_store_channel!(out, ch_idx, power, conv, g)
    return nothing
end

function _tf_batch_worker(tf_fn::F, label::String, condition_selection, kwargs) where {F}
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
