using eegfun
using BenchmarkTools
using JLD2
using DataFrames
using Statistics

println("="^60)
println("Benchmarking All Time-Frequency Methods")
println("="^60)

# 1. Setup Data Directory
data_dir = "/home/ian/Documents/Julia/output_data"
filename = "simulation_epochs_cleaned.jld2"
filepath = joinpath(data_dir, filename)

# Ensure data exists (reuse reproduction script logic if needed, but file should complicate exist)
if !isfile(filepath)
    println("Data file not found at $filepath. Generating it now...")
    # Generate Realistic Multi-Channel Data
    n_channels = 32
    n_trials = 100
    sample_rate = 256
    duration = 3.0   # seconds
    
    # Generate base signal times
    time_vec = range(-1.0, duration-1.0, length=round(Int, duration*sample_rate))
    
    # Create multiple channels
    println("Creating $n_channels channels, $n_trials trials each...")
    times, signal = eegfun.generate_signal(
        n_trials, [-1.0, 2.0], sample_rate, 
        [10.0], [1.0], [[0.0, 1.0]], 0.5
    )
    base_epochs = eegfun.signal_to_epochs(times, signal, :Ch1, Int(sample_rate))
    
    channel_names = [Symbol("Ch$i") for i in 2:n_channels]
    for df in base_epochs.data
        for col in channel_names
            df[!, col] = randn(nrow(df))
        end
    end
    mkpath(data_dir)
    jldsave(filepath; data=[base_epochs])
end

# Load Data
println("Loading data from $filepath...")
epochs = eegfun.load_data(filepath)[1]

# Configuration
toi = -0.5:0.02:2.0
foi = 1:1:40
println("\nConfiguration:")
println("  Trials: $(length(epochs.data))")
println("  Channels: $(length(names(epochs.data[1])) - 2)")
println("  Frequencies: $(length(foi))")
println("  Time points: $(length(toi))")

println("\n" * "-"^60)
println("1. Method: :multitaper (Hanning / Default)")
println("-"^60)
println("Warming up...")
eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, width=7) # width ignored
println("Benchmarking...")
@time eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, width=7)

println("\n" * "-"^60)
println("2. Method: :wavelet")
println("-"^60)
println("Warming up...")
eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)
println("Benchmarking...")
@time eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)

println("\n" * "-"^60)
println("3. Method: :superlet")
println("-"^60)
# Superlet can be slow, let's see. 
# Usage: order_min=1, order_max=30, cycles=...
# Check kwarg defaults or provide reasonable ones
println("Warming up...")
eegfun.tf_analysis(epochs, foi, toi; method=:superlet, order_min=1, order_max=3)
println("Benchmarking...")
@time eegfun.tf_analysis(epochs, foi, toi; method=:superlet, order_min=1, order_max=3)

println("\n" * "="^60)
println("Done.")
