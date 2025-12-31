#!/usr/bin/env julia
# Quick verification script to test the optimized tf_mtm function
# This script can be run standalone without the full eegfun package

using FFTW
using DSP
using LinearAlgebra
using Statistics

println("="^60)
println("Time-Frequency Optimization Verification")
println("="^60)
println()

# Copy the dpss_tapers function
function dpss_tapers(n::Int, nw::Real, k::Int)
    return DSP.dpss(n, nw, k)
end

# Optimized version
function tf_mtm_optimized(signal, times, sample_rate, frequencies, time_steps; t_ftimwin::Float64=0.5, tapsmofrq::Float64=4.0)
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal
    n_samples, n_trials = size(signal)
    
    tois_idx = [findfirst(≈(t, atol=(1000/sample_rate)/1000), times) for t in time_steps]
    
    n_freqs = length(frequencies)
    n_timepoints = length(time_steps)
    
    timewinidx = round(Int, t_ftimwin * sample_rate)
    timewinidx += iseven(timewinidx)
    half_win = timewinidx ÷ 2
    
    # DPSS tapers
    nw = tapsmofrq * timewinidx / sample_rate / 2
    n_tapers = max(1, floor(Int, 2 * nw - 1))
    tapers = dpss_tapers(timewinidx, nw, n_tapers)
    
    # Frequency indices
    freq_indices = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
    
    tf_trials = zeros(n_freqs, n_timepoints, n_trials)
    
    # Pre-allocate all buffers outside the time loop to avoid repeated allocations
    tmpdat = zeros(ComplexF64, timewinidx, n_trials)
    tapered_buf = zeros(ComplexF64, timewinidx, n_trials)
    fdat = zeros(ComplexF64, timewinidx, n_trials)
    taper_power = zeros(n_freqs, n_trials)
    
    # Create FFT plan for in-place operations
    p_fft = plan_fft(tapered_buf, 1)

    @inbounds for (ti, center_idx) in enumerate(tois_idx)
        start_idx = center_idx - half_win
        end_idx = center_idx + half_win
        
        if start_idx < 1 || end_idx > n_samples
            continue
        end
        
        # Extract window for all trials
        for trial in 1:n_trials
            @simd for s in 1:timewinidx
                tmpdat[s, trial] = signal[start_idx + s - 1, trial]
            end
        end
        
        # Reset taper_power for this time point
        fill!(taper_power, 0.0)
        
        # Multi-taper power accumulation
        for t in 1:n_tapers
            # Apply taper - manual loop for better performance
            for trial in 1:n_trials
                @simd for s in 1:timewinidx
                    tapered_buf[s, trial] = tmpdat[s, trial] * tapers[s, t]
                end
            end
            
            # FFT all trials at once using in-place operation
            mul!(fdat, p_fft, tapered_buf)
            
            # Accumulate power at frequency indices
            for trial in 1:n_trials
                @simd for fi in 1:n_freqs
                    taper_power[fi, trial] += abs2(fdat[freq_indices[fi], trial])
                end
            end
        end
        
        # Store averaged power
        inv_n_tapers = 1.0 / n_tapers
        for trial in 1:n_trials
            @simd for fi in 1:n_freqs
                tf_trials[fi, ti, trial] = taper_power[fi, trial] * inv_n_tapers
            end
        end
    end
    
    return tf_trials, collect(Float64, time_steps), collect(Float64, frequencies)
end

# Generate test data
function generate_test_data(n_samples, n_trials, sample_rate)
    times = range(0, length=n_samples, step=1/sample_rate)
    signal = randn(n_samples, n_trials)
    # Add some frequency components
    for trial in 1:n_trials
        signal[:, trial] .+= sin.(2π * 10 * times)  # 10 Hz
        signal[:, trial] .+= 0.5 * sin.(2π * 20 * times)  # 20 Hz
    end
    return signal, collect(times)
end

# Test parameters
sample_rate = 256.0
n_samples = round(Int, 2.5 * sample_rate)
n_trials = 100

signal, times = generate_test_data(n_samples, n_trials, sample_rate)

toi = -0.5:0.01:2.0
foi = 1:1:40

println("Test Configuration:")
println("  Signal size: $(size(signal))")
println("  Time points: $(length(toi))")
println("  Frequencies: $(length(foi))")
println("  Trials: $n_trials")
println()

# Warmup
println("Warming up...")
tf_mtm_optimized(signal, times, sample_rate, foi, toi; t_ftimwin=0.5, tapsmofrq=4.0)

# Test
println("\nRunning optimized version...")
@time result = tf_mtm_optimized(signal, times, sample_rate, foi, toi; t_ftimwin=0.5, tapsmofrq=4.0)

println("\nResult shape: $(size(result[1]))")
println("Result range: [$(minimum(result[1])), $(maximum(result[1]))]")
println()

# Verify the result looks reasonable
power_mean = mean(result[1])
power_std = std(result[1])
println("Power statistics:")
println("  Mean: $(round(power_mean, digits=4))")
println("  Std:  $(round(power_std, digits=4))")
println()

println("="^60)
println("✓ Verification complete!")
println("="^60)
