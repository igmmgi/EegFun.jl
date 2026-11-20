struct CleanLineThreadBuffer
    data_tapered::Matrix{Float64}
    J::Matrix{ComplexF64}
    datafit::Vector{Float64}
    datafitwin::Vector{Float64}
    datafitwin0::Vector{Float64}
    fsig::Vector{Float64}
    amps::Vector{ComplexF64}
end

"""
    cleanline!(dat::ContinuousData; line_frequencies = [50.0, 100.0, 150.0, 200.0],
               bandwidth::Real = 2.0, sliding_win_length::Real = 4.0, sliding_win_step::Real = 2.0,
               time_bandwidth::Real = 3.0, k_tapers::Int = 5, p_value::Real = 0.05,
               pad::Int = 2, channel_selection = channels())

Remove sinusoidal line noise in-place via adaptive sliding-window multi-taper F-test regression
(Chronux / EEGLAB CleanLine algorithm). Harmonics `>= sample_rate / 2` are automatically skipped.
"""
function cleanline!(
    dat::ContinuousData;
    line_frequencies::AbstractVector{<:Real} = [50.0, 100.0, 150.0, 200.0],
    bandwidth::Real = 2.0,
    sliding_win_length::Real = 4.0,
    sliding_win_step::Real = 2.0,
    time_bandwidth::Real = 3.0,
    k_tapers::Int = 5,
    p_value::Real = 0.05,
    pad::Int = 2,
    channel_selection = channels(),
)
    k_tapers <= 1 && @minimal_error("k_tapers must be >= 2 for CleanLine F-test regression, got $k_tapers")
    sliding_win_length > 0 || @minimal_error("sliding_win_length must be positive, got $sliding_win_length")
    sliding_win_step > 0 || @minimal_error("sliding_win_step must be positive, got $sliding_win_step")

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    isempty(selected_channels) && return dat

    Fs = Float64(dat.sample_rate)
    nyquist = Fs / 2.0
    valid_frequencies = Float64[Float64(f) for f in line_frequencies if 0.0 < f < nyquist]
    if isempty(valid_frequencies)
        @minimal_warning "No line_frequencies below Nyquist ($nyquist Hz); CleanLine skipped."
        return dat
    end

    N = nrow(dat.data)
    Nwin = round(Int, sliding_win_length * Fs)
    if N < Nwin
        @minimal_warning "Data length ($N samples) is shorter than CleanLine window length ($Nwin samples). CleanLine skipped."
        return dat
    end
    Nstep = max(1, round(Int, sliding_win_step * Fs))
    Noverlap = max(0, Nwin - Nstep)
    nfft = max(nextpow(2, Nwin) * (2^pad), Nwin)
    freqs = rfftfreq(nfft, Fs)
    Nf = length(freqs)

    @info "Starting CleanLine Multi-Taper removal..."
    @info "  Window: $(sliding_win_length)s, Step: $(sliding_win_step)s"
    @info "  Tapers: TW=$(time_bandwidth), K=$(k_tapers), NFFT=$(nfft)"

    tapers = Matrix{Float64}(DSP.dpss(Nwin, time_bandwidth, k_tapers))
    Kodd = 1:2:k_tapers
    Keven = 2:2:k_tapers
    H0 = Float64[sum(@view tapers[:, k]) for k in Kodd]
    H0sq = sum(abs2, H0)

    smooth = if Noverlap > 0
        tau = 100.0
        Float64[1.0 / (1.0 + exp(-tau * (x - Noverlap / 2.0) / Noverlap)) for x = 1:Noverlap]
    else
        Float64[]
    end

    sig_threshold = quantile(FDist(2, 2 * k_tapers - 2), 1.0 - p_value / Nwin)

    # Precompute candidate frequency bin ranges once across all windows and channels
    fscanbw = Float64(bandwidth)
    bin_ranges = UnitRange{Int}[]
    for f_target in valid_frequencies
        if fscanbw > 0.0
            f_min = max(0.0, f_target - fscanbw / 2.0)
            f_max = min(nyquist - eps(nyquist), f_target + fscanbw / 2.0)
            idx_start = searchsortedfirst(freqs, f_min)
            idx_end = searchsortedlast(freqs, f_max)
            idx_start <= idx_end && push!(bin_ranges, idx_start:idx_end)
        else
            _, idx = findmin(f -> abs(f - f_target), freqs)
            freqs[idx] < nyquist && push!(bin_ranges, idx:idx)
        end
    end

    _foreach_chunk(length(selected_channels)) do ch_range
        buf = CleanLineThreadBuffer(
            zeros(Float64, nfft, k_tapers),
            zeros(ComplexF64, Nf, k_tapers),
            Vector{Float64}(undef, N),
            zeros(Float64, Nwin),
            zeros(Float64, Nwin),
            sizehint!(Float64[], length(bin_ranges)),
            sizehint!(ComplexF64[], length(bin_ranges)),
        )
        plan = plan_rfft(buf.data_tapered, 1)
        for idx in ch_range
            ch = selected_channels[idx]
            _cleanline_channel!(
                dat.data[!, ch],
                tapers,
                Fs,
                Nwin,
                Nstep,
                Noverlap,
                smooth,
                bin_ranges,
                sig_threshold,
                freqs,
                H0,
                H0sq,
                Kodd,
                Keven,
                plan,
                buf,
            )
        end
    end

    return dat
end

function cleanline!(dat::AbstractVector{<:ContinuousData}; kwargs...)
    for d in dat
        cleanline!(d; kwargs...)
    end
    return dat
end

"""
    cleanline(dat::ContinuousData, args...; kwargs...) -> ContinuousData
    cleanline(dat::AbstractVector{<:ContinuousData}, args...; kwargs...) -> Vector{ContinuousData}

Non-mutating version of [`cleanline!`](@ref).
"""
cleanline(dat::ContinuousData, args...; kwargs...) = cleanline!(copy(dat), args...; kwargs...)
cleanline(dat::AbstractVector{<:ContinuousData}, args...; kwargs...) = [cleanline(d, args...; kwargs...) for d in dat]

function _cleanline_channel!(
    data::AbstractVector{<:Real},
    tapers::Matrix{Float64},
    Fs::Float64,
    Nwin::Int,
    Nstep::Int,
    Noverlap::Int,
    smooth::Vector{Float64},
    bin_ranges::Vector{UnitRange{Int}},
    sig_threshold::Float64,
    freqs::AbstractVector{<:Real},
    H0::Vector{Float64},
    H0sq::Float64,
    Kodd::StepRange{Int,Int},
    Keven::StepRange{Int,Int},
    plan,
    buf::CleanLineThreadBuffer,
)
    N = length(data)
    winstart = 1:Nstep:(N-Nwin+1)
    nw = length(winstart)
    datafit = buf.datafit
    fill!(datafit, 0.0)

    for n = 1:nw
        s_idx = winstart[n]
        idx_range = s_idx:(s_idx+Nwin-1)
        _cleanline_ftestc_fit!(buf, @view(data[idx_range]), tapers, Fs, bin_ranges, sig_threshold, freqs, H0, H0sq, Kodd, Keven, plan)

        if n > 1 && Noverlap > 0
            offset = Nwin - Noverlap
            @inbounds @simd for i = 1:Noverlap
                w = smooth[i]
                buf.datafitwin[i] = w * buf.datafitwin[i] + (1.0 - w) * buf.datafitwin0[offset+i]
            end
        end

        @inbounds copyto!(datafit, s_idx, buf.datafitwin, 1, Nwin)
        @inbounds copyto!(buf.datafitwin0, 1, buf.datafitwin, 1, Nwin)
    end

    prev_end = winstart[nw] + Nwin - 1
    if prev_end < N
        tail_start = N - Nwin + 1
        tail_overlap = prev_end - tail_start + 1
        _cleanline_ftestc_fit!(buf, @view(data[tail_start:N]), tapers, Fs, bin_ranges, sig_threshold, freqs, H0, H0sq, Kodd, Keven, plan)

        if tail_overlap > 0
            tau = 100.0
            @inbounds for i = 1:tail_overlap
                w = 1.0 / (1.0 + exp(-tau * (i - tail_overlap / 2.0) / tail_overlap))
                buf.datafitwin[i] = w * buf.datafitwin[i] + (1.0 - w) * datafit[tail_start+i-1]
            end
        end

        @inbounds copyto!(datafit, tail_start, buf.datafitwin, 1, Nwin)
    end

    @inbounds @simd for i = 1:N
        data[i] -= datafit[i]
    end
    return nothing
end

function _cleanline_ftestc_fit!(
    buf::CleanLineThreadBuffer,
    datawin::AbstractVector{<:Real},
    tapers::Matrix{Float64},
    Fs::Float64,
    bin_ranges::Vector{UnitRange{Int}},
    sig_threshold::Float64,
    freqs::AbstractVector{<:Real},
    H0::Vector{Float64},
    H0sq::Float64,
    Kodd::StepRange{Int,Int},
    Keven::StepRange{Int,Int},
    plan,
)
    Nwin = length(datawin)
    K = size(tapers, 2)
    nfft = size(buf.data_tapered, 1)

    win_sum = 0.0
    @inbounds @simd for i = 1:Nwin
        win_sum += Float64(datawin[i])
    end
    win_mean = win_sum / Nwin

    @inbounds for k = 1:K
        @simd for i = 1:Nwin
            buf.data_tapered[i, k] = (Float64(datawin[i]) - win_mean) * tapers[i, k]
        end
        if Nwin < nfft
            fill!(@view(buf.data_tapered[(Nwin+1):nfft, k]), 0.0)
        end
    end

    mul!(buf.J, plan, buf.data_tapered)

    empty!(buf.fsig)
    empty!(buf.amps)
    inv_H0sq = 1.0 / H0sq
    num_scale = (K - 1) * H0sq
    n_odd = length(Kodd)

    # Evaluate Thomson F-statistic only at candidate frequency bins around line frequencies
    @inbounds for rng in bin_ranges
        max_fval = -1.0
        max_idx = -1
        best_A = zero(ComplexF64)
        for f in rng
            jp_h0 = zero(ComplexF64)
            for idx_odd = 1:n_odd
                jp_h0 += buf.J[f, Kodd[idx_odd]] * H0[idx_odd]
            end
            A_f = jp_h0 * inv_H0sq

            den_f = 0.0
            for idx_odd = 1:n_odd
                den_f += abs2(buf.J[f, Kodd[idx_odd]] - A_f * H0[idx_odd])
            end
            for k_even in Keven
                den_f += abs2(buf.J[f, k_even])
            end

            fval = den_f > 0.0 ? (num_scale * abs2(A_f)) / den_f : 0.0
            if fval > sig_threshold && fval > max_fval
                max_fval = fval
                max_idx = f
                best_A = A_f
            end
        end
        if max_idx != -1
            push!(buf.fsig, Float64(freqs[max_idx]))
            push!(buf.amps, best_A)
        end
    end

    fill!(buf.datafitwin, 0.0)
    @inbounds for i in eachindex(buf.fsig)
        phase_arg = 2π * buf.fsig[i] / Fs
        step = cis(phase_arg)
        z = 2.0 * buf.amps[i]
        for j = 1:Nwin
            buf.datafitwin[j] += real(z)
            z *= step
        end
    end
    return nothing
end
