# === ADAPTIVE MULTI-TAPER LINE NOISE REMOVAL (CLEANLINE) ===

struct CleanLineThreadBuffer
    data_tapered::Matrix{Float64}
    J::Matrix{ComplexF64}
    JpH0::Vector{ComplexF64}
    A::Vector{ComplexF64}
    Jhat::Matrix{ComplexF64}
    num::Vector{Float64}
    den::Vector{Float64}
    Fval::Vector{Float64}
    datafitwin::Vector{Float64}
    datafitwin0::Vector{Float64}
end

"""
    cleanline!(dat::ContinuousData; kwargs...)

Remove line noise using adaptive sliding-window multi-taper regression.
This is a Julia port of the Chronux / EEGLAB CleanLine plugin, utilizing 
`DSP.jl` for DPSS (Slepian) tapers and natively threaded sliding windows.

# Keyword Arguments
- `line_frequencies::AbstractVector{<:Real}`: Base frequencies to remove (default: `[50.0, 100.0, 150.0, 200.0]`). Harmonics `>= sample_rate / 2` are automatically skipped.
- `bandwidth::Real`: Bandwidth to scan around each line frequency (default: `2.0`).
- `sliding_win_length::Real`: Length of the sliding window in seconds (default: `4.0`).
- `sliding_win_step::Real`: Step size of the sliding window in seconds (default: `2.0`).
- `time_bandwidth::Real`: Time-bandwidth product (TW) for DPSS tapers (default: `3.0`).
- `k_tapers::Int`: Number of tapers to use (`k_tapers >= 2`, typically `<= 2*TW - 1`; default: `5`).
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
    if k_tapers <= 1
        @minimal_error("k_tapers must be >= 2 for CleanLine F-test regression, got $k_tapers")
    end
    sliding_win_length > 0 || @minimal_error("sliding_win_length must be positive, got $sliding_win_length")
    sliding_win_step > 0 || @minimal_error("sliding_win_step must be positive, got $sliding_win_step")

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        return dat
    end

    Fs = dat.sample_rate
    nyquist = Fs / 2.0
    valid_frequencies = Float64[Float64(f) for f in line_frequencies if 0.0 < f < nyquist]
    if isempty(valid_frequencies)
        @minimal_warning "No line_frequencies below Nyquist ($nyquist Hz); CleanLine skipped."
        return dat
    end

    Nwin = round(Int, sliding_win_length * Fs)
    if nrow(dat.data) < Nwin
        @minimal_warning "Data length ($(nrow(dat.data)) samples) is shorter than CleanLine window length ($Nwin samples). CleanLine skipped."
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

    # 1. Generate DPSS (Slepian) tapers
    tapers = DSP.dpss(Nwin, time_bandwidth, k_tapers)
    Kodd = 1:2:k_tapers
    Keven = 2:2:k_tapers
    H0 = vec(sum(tapers[:, Kodd], dims = 1))
    H0sq = sum(abs2, H0)

    # 2. Pre-calculate overlap smoothing function (sigmoidal)
    tau = 100.0
    smooth = if Noverlap > 0
        x = 1:Noverlap
        @. 1.0 / (1.0 + exp(-tau * (x - Noverlap / 2.0) / Noverlap))
    else
        Float64[]
    end

    # 4. Process each channel in parallel
    p_corrected = p_value / Nwin
    sig_threshold = quantile(FDist(2, 2 * k_tapers - 2), 1.0 - p_corrected)

    Threads.@threads for ch in selected_channels
        # Create FFT plan per-thread to ensure thread safety
        plan = plan_rfft(zeros(Float64, nfft, k_tapers), 1)

        buf = CleanLineThreadBuffer(
            zeros(Float64, nfft, k_tapers),
            zeros(ComplexF64, Nf, k_tapers),
            zeros(ComplexF64, Nf),
            zeros(ComplexF64, Nf),
            zeros(ComplexF64, Nf, length(Kodd)),
            zeros(Float64, Nf),
            zeros(Float64, Nf),
            zeros(Float64, Nf),
            zeros(Float64, Nwin),
            zeros(Float64, Nwin),
        )
        _cleanline_channel!(
            dat.data[!, ch],
            tapers,
            Fs,
            Nwin,
            Nstep,
            Noverlap,
            smooth,
            valid_frequencies,
            bandwidth,
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

Non-mutating version of [`cleanline!`](@ref). Returns a copy of `dat` with line noise removed.
"""
function cleanline(dat::ContinuousData, args...; kwargs...)
    return cleanline!(copy(dat), args...; kwargs...)
end

function cleanline(dat::AbstractVector{<:ContinuousData}, args...; kwargs...)
    return [cleanline(d, args...; kwargs...) for d in dat]
end

function _cleanline_channel!(
    data::AbstractVector{<:Real},
    tapers::AbstractMatrix{<:Real},
    Fs::Real,
    Nwin::Int,
    Nstep::Int,
    Noverlap::Int,
    smooth::AbstractVector{<:Real},
    f0::AbstractVector{<:Real},
    fscanbw::Real,
    sig_threshold::Real,
    freqs::AbstractVector{<:Real},
    H0::AbstractVector{<:Real},
    H0sq::Real,
    Kodd::AbstractVector{Int},
    Keven::AbstractVector{Int},
    plan,
    buf::CleanLineThreadBuffer,
)
    N = length(data)
    winstart = 1:Nstep:(N-Nwin+1)
    nw = length(winstart)

    datafit = zeros(Float64, N)

    for n = 1:nw
        idx_range = winstart[n]:(winstart[n]+Nwin-1)
        datawin = @view data[idx_range]

        _cleanline_ftestc_fit!(
            buf,
            datawin,
            tapers,
            Float64(Fs),
            f0,
            Float64(fscanbw),
            Float64(sig_threshold),
            freqs,
            H0,
            Float64(H0sq),
            Kodd,
            Keven,
            plan,
        )

        # Apply sigmoidal smoothing in the overlap region
        if n > 1 && Noverlap > 0
            for i = 1:Noverlap
                buf.datafitwin[i] = smooth[i] * buf.datafitwin[i] + (1.0 - smooth[i]) * buf.datafitwin0[Nwin-Noverlap+i]
            end
        end

        datafit[idx_range] .= buf.datafitwin
        buf.datafitwin0 .= buf.datafitwin
    end

    # Clean remaining tail samples if (N - Nwin) is not an exact multiple of Nstep
    prev_end = winstart[nw] + Nwin - 1
    if prev_end < N
        tail_start = N - Nwin + 1
        tail_range = tail_start:N
        tail_overlap = prev_end - tail_start + 1
        datawin_tail = @view data[tail_range]

        _cleanline_ftestc_fit!(
            buf,
            datawin_tail,
            tapers,
            Float64(Fs),
            f0,
            Float64(fscanbw),
            Float64(sig_threshold),
            freqs,
            H0,
            Float64(H0sq),
            Kodd,
            Keven,
            plan,
        )

        if tail_overlap > 0
            tau = 100.0
            for i = 1:tail_overlap
                w = 1.0 / (1.0 + exp(-tau * (i - tail_overlap / 2.0) / tail_overlap))
                buf.datafitwin[i] = w * buf.datafitwin[i] + (1.0 - w) * datafit[tail_start+i-1]
            end
        end

        datafit[tail_range] .= buf.datafitwin
    end

    # Subtract the fitted line noise from the original data
    for i = 1:N
        data[i] -= datafit[i]
    end
end

function _cleanline_ftestc_fit!(
    buf::CleanLineThreadBuffer,
    datawin::AbstractVector{<:Real},
    tapers::AbstractMatrix{<:Real},
    Fs::Real,
    f0::AbstractVector{<:Real},
    fscanbw::Real,
    sig_threshold::Real,
    freqs::AbstractVector{<:Real},
    H0::AbstractVector{<:Real},
    H0sq::Real,
    Kodd::AbstractVector{Int},
    Keven::AbstractVector{Int},
    plan,
)
    Nwin = length(datawin)
    K = size(tapers, 2)
    Nf = length(freqs)
    nyquist = Float64(Fs) / 2.0

    # Calculate window mean to prevent DC spectral leakage into the 50Hz band
    win_mean = sum(datawin) / Nwin

    fill!(buf.data_tapered, 0.0)
    for k = 1:K
        @inbounds @simd for i = 1:Nwin
            buf.data_tapered[i, k] = (Float64(datawin[i]) - win_mean) * Float64(tapers[i, k])
        end
    end

    # 1. Tapered FFT
    mul!(buf.J, plan, buf.data_tapered)

    Jp = @view buf.J[:, Kodd]
    J_even = @view buf.J[:, Keven]

    # 4. F-statistic
    fill!(buf.JpH0, zero(ComplexF64))
    for k = 1:length(Kodd)
        h = Float64(H0[k])
        @inbounds @simd for f = 1:Nf
            buf.JpH0[f] += Jp[f, k] * h
        end
    end

    h0sq_f = Float64(H0sq)
    @inbounds @simd for f = 1:Nf
        buf.A[f] = buf.JpH0[f] / h0sq_f
    end

    for k = 1:length(Kodd)
        h = Float64(H0[k])
        @inbounds @simd for f = 1:Nf
            buf.Jhat[f, k] = buf.A[f] * h
        end
    end

    fill!(buf.den, 0.0)
    for k = 1:length(Kodd)
        @inbounds @simd for f = 1:Nf
            buf.den[f] += abs2(Jp[f, k] - buf.Jhat[f, k])
        end
    end
    for k = 1:length(Keven)
        @inbounds @simd for f = 1:Nf
            buf.den[f] += abs2(J_even[f, k])
        end
    end

    @inbounds @simd for f = 1:Nf
        buf.num[f] = (K - 1) * abs2(buf.A[f]) * h0sq_f
        buf.Fval[f] = buf.den[f] > 0.0 ? (buf.num[f] / buf.den[f]) : 0.0
    end

    # 5. Significance and Peak Finding
    fsig = Float64[]
    amps = ComplexF64[]

    for f_target_raw in f0
        f_target = Float64(f_target_raw)
        (f_target <= 0.0 || f_target >= nyquist) && continue
        if fscanbw > 0.0
            f_min = max(0.0, f_target - fscanbw / 2.0)
            f_max = min(nyquist - eps(nyquist), f_target + fscanbw / 2.0)

            idx_start = searchsortedfirst(freqs, f_min)
            idx_end = searchsortedlast(freqs, f_max)

            if idx_start <= idx_end
                max_fval = -1.0
                max_idx = -1
                for i = idx_start:idx_end
                    if buf.Fval[i] > sig_threshold && buf.Fval[i] > max_fval
                        max_fval = buf.Fval[i]
                        max_idx = i
                    end
                end

                if max_idx != -1
                    push!(fsig, Float64(freqs[max_idx]))
                    push!(amps, buf.A[max_idx])
                end
            end
        else
            _, idx = findmin(abs.(freqs .- f_target))
            if freqs[idx] < nyquist && buf.Fval[idx] > sig_threshold
                push!(fsig, Float64(freqs[idx]))
                push!(amps, buf.A[idx])
            end
        end
    end

    # 6. Reconstruct time-domain signal
    fill!(buf.datafitwin, 0.0)
    if !isempty(fsig)
        for i = 1:length(fsig)
            f = fsig[i]
            amp = amps[i]
            phase_arg = 2π * f / Fs
            @inbounds @simd for j = 1:Nwin
                buf.datafitwin[j] += 2.0 * real(exp(im * phase_arg * (j - 1)) * amp)
            end
        end
    end
end
