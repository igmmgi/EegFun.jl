# Simplified benchmark that doesn't require full package
using FFTW
using DSP
using BenchmarkTools
using Statistics

# Copy the dpss_tapers function
function dpss_tapers(n::Int, nw::Real, k::Int)
    return DSP.dpss(n, nw, k)
end

# Copy the current tf_mtm implementation
function tf_mtm_current(signal, times, sample_rate, frequencies, time_steps; t_ftimwin::Float64=0.5, tapsmofrq::Float64=4.0)
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
    
    # Pre-allocate buffer and FFT plan
    tmpdat = zeros(ComplexF64, timewinidx, n_trials)
    p_fft = plan_fft(tmpdat, 1)
    
    # Temporary buffer for tapered signals
    tapered_buf = zeros(ComplexF64, timewinidx, n_trials)

    @inbounds for (ti, center_idx) in enumerate(tois_idx)
        start_idx = center_idx - half_win
        end_idx = center_idx + half_win
        
        if start_idx < 1 || end_idx > n_samples
            continue
        end
        
        # Extract window for all trials
        @views tmpdat .= signal[start_idx:end_idx, :]
        
        # Multi-taper power accumulation
        taper_power = zeros(n_freqs, n_trials)
        for t in 1:n_tapers
            # Apply taper
            @views tapered_buf .= tmpdat .* tapers[:, t]
            
            # FFT all trials at once
            fdat = p_fft * tapered_buf
            
            # Accumulate power
            @views taper_power .+= abs2.(fdat[freq_indices, :])
        end
        
        tf_trials[:, ti, :] .= taper_power ./ n_tapers
    end
    
    return tf_trials, collect(Float64, time_steps), collect(Float64, frequencies)
end

# Generate test data similar to what would come from EEG
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

# Test parameters matching user's request
sample_rate = 256.0
n_samples = round(Int, 2.5 * sample_rate)  # 2.5 seconds of data
n_trials = 100  # typical number of trials

signal, times = generate_test_data(n_samples, n_trials, sample_rate)

toi = -0.5:0.01:2.0
foi = 1:1:40

println("Signal size: $(size(signal))")
println("Number of time points: $(length(toi))")
println("Number of frequencies: $(length(foi))")
println("Number of trials: $n_trials")
println()

# Warmup
println("Warming up...")
tf_mtm_current(signal, times, sample_rate, foi, toi; t_ftimwin=0.5, tapsmofrq=4.0)

println("\nBenchmarking current implementation...")
@time result = tf_mtm_current(signal, times, sample_rate, foi, toi; t_ftimwin=0.5, tapsmofrq=4.0)

println("\nDetailed benchmark:")
@btime tf_mtm_current($signal, $times, $sample_rate, $foi, $toi; t_ftimwin=0.5, tapsmofrq=4.0)

println("\nResult size: $(size(result[1]))")
