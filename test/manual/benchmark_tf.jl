using eegfun
using BenchmarkTools
using DataFrames

println("="^60)
println("Time-Frequency Performance Test")
println("="^60)
println()

# Generate synthetic test data
println("Generating synthetic test data...")
sample_rate = 256.0
times, signal = eegfun.generate_signal(
    100,                                    # n_trials (typical for EEG)
    [-0.5, 2.5],                           # time_window
    sample_rate,                            # sample_rate
    [10.0, 20.0],                          # frequencies
    [1.0, 1.0],                            # amplitudes
    [[0.0, 1.0], [1.0, 2.0]],             # time windows for each freq
    0.5                                     # noise amplitude
)

# Convert to EpochData
epochs = eegfun.signal_to_epochs(times, signal, :TestChannel, Int(sample_rate))

# Time and frequency parameters
toi = -0.5:0.01:2.0
foi = 1:1:40

println("Test configuration:")
println("  Trials: $(length(epochs.data))")
println("  Samples per trial: $(nrow(epochs.data[1]))")
println("  Time points: $(length(toi))")
println("  Frequencies: $(length(foi))")
println()

# Warmup
println("Warming up...")
eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, taper=:dpss, t_ftimwin=0.5, tapsmofrq=4.0)

println("\n" * "="^60)
println("BENCHMARK: Multitaper (DPSS) Method")
println("="^60)
println()

# Run benchmark
println("Running benchmark (this may take a few seconds)...")
@time tf_data = eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, taper=:dpss, t_ftimwin=0.5, tapsmofrq=4.0)

println("\nDetailed benchmark:")
@btime eegfun.tf_analysis($epochs, $foi, $toi; method=:multitaper, taper=:dpss, t_ftimwin=0.5, tapsmofrq=4.0)

println("\n" * "="^60)
println("Results:")
println("="^60)
println("Output size: $(size(tf_data.data))")
println("Channels: $(names(tf_data.data)[3:end])")
println()
println("Expected performance: ~2-3 seconds (down from ~7 seconds)")
println("="^60)
