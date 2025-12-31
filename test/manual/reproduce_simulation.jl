using eegfun
using BenchmarkTools
using JLD2
using DataFrames

println("="^60)
println("Reproducing User Scenario with Realistic Data")
println("="^60)

# 1. Setup Data Directory
data_dir = "/home/ian/Documents/Julia/output_data"
mkpath(data_dir)

# 2. Generate Realistic Multi-Channel Data
println("\nGenerating realistic multi-channel dataset...")
n_channels = 32 # Typical EEG setup
n_trials = 100
sample_rate = 256
duration = 3.0   # seconds

# Generate base signal times
time_vec = range(-1.0, duration-1.0, length=round(Int, duration*sample_rate))

# Create multiple channels
println("Creating $n_channels channels, $n_trials trials each...")
# We use eegfun.signal_to_epochs for 1 channel, then we'll add more manually
# to simulate a real loaded file
times, signal = eegfun.generate_signal(
    n_trials, [-1.0, 2.0], sample_rate, 
    [10.0], [1.0], [[0.0, 1.0]], 0.5
)
base_epochs = eegfun.signal_to_epochs(times, signal, :Ch1, Int(sample_rate))

# Add more channels to the DataFrames
channel_names = [Symbol("Ch$i") for i in 2:n_channels]
for df in base_epochs.data
    for col in channel_names
        df[!, col] = randn(nrow(df)) # Add random noise channels
    end
end
# Update layout (hacky but works for this test)
# In real EpochData, layout matches channels. eegfun might rely on layout or just columns.
# tf_analysis selects channels. By default it might select all?
# The user's code: `epochs = eegfun.load_data(epoch_file)[1]` then `tf_analysis(epochs, ...)`
# tf_analysis default channel_selection is `channels()`, which usually implies all data channels.

# Save to file to mimic user workflow
filename = "simulation_epochs_cleaned.jld2"
filepath = joinpath(data_dir, filename)
println("Saving to $filepath...")
jldsave(filepath; data=[base_epochs]) # Save as Vector{EpochData}

println("\n" * "="^60)
println("Running User's Code")
println("="^60)

# 3. User's Exact Code Block
# ---------------------------------------------------------
# Try to load epoched data from output_data
# data_dir defined above
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
epoch_file = joinpath(data_dir, epoch_files[1])
epochs = eegfun.load_data(epoch_file)[1] # take single epoch

# Time and frequency parameters
toi = -0.5:0.01:2
foi = 1:1:40

# Precision check: method=:multitaper with width=7
# Note: width=7 is ignored by multitaper, so this runs default Hanning taper
println("Configuration:")
println("  Method: :multitaper (width=7 ignored -> uses defaults)")
println("  Channels: $(length(names(epochs.data[1])) - 2) (excluding time/epoch)")
println("  Trials: $(length(epochs.data))")
println("  Frequencies: $(length(foi))")
println("  Time points: $(length(toi))")
println()

println("Warming up...")
eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, width=7)

println("\nBenchmarking...")
@time tf_data1 = eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, width=7)

println("\nDetailed benchmark:")
@btime eegfun.tf_analysis($epochs, $foi, $toi; method=:multitaper, width=7)

println("\n" * "="^60)
println("Done.")
