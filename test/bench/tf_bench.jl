using eegfun
using BenchmarkTools

# Standard parameters
sample_rate = 256.0
n_trials = 100
time_window = [-1.0, 2.0]
foi = 2:2:40
toi = -0.5:0.02:1.5

times, signal = eegfun.generate_signal(
    n_trials, time_window, sample_rate, 
    [10.0, 25.0], [1.0, 2.0], [[0.0, 1.0], [1.0, 2.0]], 0.5
)

println("--- BENCHMARK RESULTS (BASELINE) ---")
println("n_trials: $n_trials, n_freqs: $(length(foi)), n_timepoints: $(length(toi))")

println("\n1. tf_wavelet (7 cycles):")
@btime eegfun.tf_wavelet($signal, $times, $sample_rate, $foi, $toi; width=7)

println("\n2. tf_superlet (order 5):")
@btime eegfun.tf_superlet($signal, $times, $sample_rate, $foi, $toi; order=5)

println("\n3. tf_multitaper (Hanning, 0.3s window):")
@btime eegfun.tf_multitaper($signal, $times, $sample_rate, $foi, $toi; taper=:hanning, t_ftimwin=0.3)

println("\n4. tf_multitaper (DPSS, 0.3s window, 4Hz smoothing):")
@btime eegfun.tf_multitaper($signal, $times, $sample_rate, $foi, $toi; taper=:dpss, t_ftimwin=0.3, tapsmofrq=4.0)
