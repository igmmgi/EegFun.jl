# =============================================================================
# ADAPTIVE MULTI-TAPER LINE NOISE REMOVAL (CLEANLINE)
# =============================================================================

"""
    CleanLineParams

Parameters for adaptive multi-taper line noise removal, based on the Chronux 
toolbox / EEGLAB cleanline plugin algorithms.
"""
Base.@kwdef struct CleanLineParams
    line_frequencies::Vector{Float64} = [50.0, 100.0, 150.0, 200.0]
    scan_for_lines::Bool = true
    bandwidth::Float64 = 2.0
    sliding_win_length::Float64 = 4.0
    sliding_win_step::Float64 = 2.0
    smoothing_tau::Float64 = 100.0
    pad::Int = 0
    p_value::Float64 = 0.05
    time_bandwidth::Float64 = 3.0
    k_tapers::Int = 5
end

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
- `line_frequencies::Vector{Float64}`: Base frequencies to remove (e.g., `[50, 100, 150]`).
- `bandwidth::Float64`: Bandwidth to scan around each line frequency.
- `sliding_win_length::Float64`: Length of the sliding window in seconds.
- `sliding_win_step::Float64`: Step size of the sliding window in seconds.
- `time_bandwidth::Float64`: Time-bandwidth product (TW) for DPSS tapers.
- `k_tapers::Int`: Number of tapers to use (K <= 2TW-1).
"""
function cleanline!(
    dat::ContinuousData;
    line_frequencies::Vector{Float64} = [50.0, 100.0, 150.0, 200.0],
    bandwidth::Float64 = 2.0,
    sliding_win_length::Float64 = 4.0,
    sliding_win_step::Float64 = 2.0,
    time_bandwidth::Float64 = 3.0,
    k_tapers::Int = 5,
    p_value::Float64 = 0.05,
    pad::Int = 2,
    channel_selection::Function = channels()
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta=false, include_extra=false)
    if isempty(selected_channels)
        return nothing
    end

    Fs = dat.sample_rate
    Nwin = round(Int, sliding_win_length * Fs)
    Nstep = round(Int, sliding_win_step * Fs)
    Noverlap = Nwin - Nstep
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
    H0 = vec(sum(tapers[:, Kodd], dims=1))
    H0sq = sum(abs2, H0)

    # 2. Pre-calculate overlap smoothing function (sigmoidal)
    tau = 100.0
    x = 1:Noverlap
    smooth = @. 1.0 / (1.0 + exp(-tau * (x - Noverlap / 2.0) / Noverlap))

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
            zeros(Float64, Nwin)
        )
        _cleanline_channel!(dat.data[!, ch], tapers, Fs, Nwin, Nstep, Noverlap, smooth, line_frequencies, bandwidth, sig_threshold, freqs, H0, H0sq, Kodd, Keven, plan, buf)
    end

    return nothing
end
@add_nonmutating cleanline!

function _cleanline_channel!(data::Vector{Float64}, tapers::Matrix{Float64}, Fs::Real, Nwin::Int, Nstep::Int, Noverlap::Int, smooth::Vector{Float64}, f0::Vector{Float64}, fscanbw::Float64, sig_threshold::Float64, freqs::AbstractVector{Float64}, H0::Vector{Float64}, H0sq::Float64, Kodd::StepRange{Int,Int}, Keven::StepRange{Int,Int}, plan, buf::CleanLineThreadBuffer)
    N = length(data)
    winstart = 1:Nstep:(N - Nwin + 1)
    nw = length(winstart)
    
    datafit = zeros(winstart[nw] + Nwin - 1)
    
    for n in 1:nw
        idx_range = winstart[n]:(winstart[n] + Nwin - 1)
        datawin = @view data[idx_range]

        _cleanline_ftestc_fit!(buf, datawin, tapers, Float64(Fs), f0, fscanbw, sig_threshold, freqs, H0, H0sq, Kodd, Keven, plan)
        
        # Apply sigmoidal smoothing in the overlap region
        if n > 1
            for i in 1:Noverlap
                buf.datafitwin[i] = smooth[i] * buf.datafitwin[i] + (1.0 - smooth[i]) * buf.datafitwin0[Nwin - Noverlap + i]
            end
        end
        
        datafit[idx_range] .= buf.datafitwin
        buf.datafitwin0 .= buf.datafitwin
    end
    
    # Subtract the fitted line noise from the original data
    for i in 1:length(datafit)
        data[i] -= datafit[i]
    end
end

function _cleanline_ftestc_fit!(buf::CleanLineThreadBuffer, datawin::AbstractVector{Float64}, tapers::Matrix{Float64}, Fs::Float64, f0::Vector{Float64}, fscanbw::Float64, sig_threshold::Float64, freqs::AbstractVector{Float64}, H0::Vector{Float64}, H0sq::Float64, Kodd::StepRange{Int,Int}, Keven::StepRange{Int,Int}, plan)
    Nwin = length(datawin)
    K = size(tapers, 2)
    Nf = length(freqs)
    
    # Calculate window mean to prevent DC spectral leakage into the 50Hz band
    win_mean = sum(datawin) / Nwin
    
    fill!(buf.data_tapered, 0.0)
    for k in 1:K
        @inbounds @simd for i in 1:Nwin
            buf.data_tapered[i, k] = (datawin[i] - win_mean) * tapers[i, k]
        end
    end
    
    # 1. Tapered FFT
    mul!(buf.J, plan, buf.data_tapered)
    
    Jp = @view buf.J[:, Kodd]
    J_even = @view buf.J[:, Keven]
    
    # 4. F-statistic
    fill!(buf.JpH0, zero(ComplexF64))
    for k in 1:length(Kodd)
        h = H0[k]
        @inbounds @simd for f in 1:Nf
            buf.JpH0[f] += Jp[f, k] * h
        end
    end
    
    @inbounds @simd for f in 1:Nf
        buf.A[f] = buf.JpH0[f] / H0sq
    end
    
    for k in 1:length(Kodd)
        h = H0[k]
        @inbounds @simd for f in 1:Nf
            buf.Jhat[f, k] = buf.A[f] * h
        end
    end
    
    fill!(buf.den, 0.0)
    for k in 1:length(Kodd)
        @inbounds @simd for f in 1:Nf
            buf.den[f] += abs2(Jp[f, k] - buf.Jhat[f, k])
        end
    end
    for k in 1:length(Keven)
        @inbounds @simd for f in 1:Nf
            buf.den[f] += abs2(J_even[f, k])
        end
    end
    
    @inbounds @simd for f in 1:Nf
        buf.num[f] = (K - 1) * abs2(buf.A[f]) * H0sq
        buf.Fval[f] = buf.num[f] / buf.den[f]
    end
    
    # 5. Significance and Peak Finding
    fsig = Float64[]
    amps = ComplexF64[]
    
    for f_target in f0
        if fscanbw > 0.0
            f_min = f_target - fscanbw / 2.0
            f_max = f_target + fscanbw / 2.0
            
            idx_start = searchsortedfirst(freqs, f_min)
            idx_end = searchsortedlast(freqs, f_max)
            
            if idx_start <= idx_end
                max_fval = -1.0
                max_idx = -1
                for i in idx_start:idx_end
                    if buf.Fval[i] > sig_threshold && buf.Fval[i] > max_fval
                        max_fval = buf.Fval[i]
                        max_idx = i
                    end
                end
                
                if max_idx != -1
                    push!(fsig, freqs[max_idx])
                    push!(amps, buf.A[max_idx])
                end
            end
        else
            _, idx = findmin(abs.(freqs .- f_target))
            if buf.Fval[idx] > sig_threshold
                push!(fsig, freqs[idx])
                push!(amps, buf.A[idx])
            end
        end
    end
    
    # 6. Reconstruct time-domain signal
    fill!(buf.datafitwin, 0.0)
    if !isempty(fsig)
        for i in 1:length(fsig)
            f = fsig[i]
            amp = amps[i]
            phase_arg = 2π * f / Fs
            @inbounds @simd for j in 1:Nwin
                buf.datafitwin[j] += 2.0 * real(exp(im * phase_arg * (j - 1)) * amp)
            end
        end
    end
end
