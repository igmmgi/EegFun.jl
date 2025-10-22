#!/usr/bin/env julia

using Pkg
Pkg.activate(".")

using eegfun
using DataFrames

# Create some test data
println("Creating test data...")
n_samples = 1000
n_channels = 4
time = collect(range(0, 10, length=n_samples))
channel_names = [:Fz, :Cz, :Pz, :Oz]

# Create random EEG-like data
data_matrix = randn(n_samples, n_channels) .* 50  # 50 microvolts amplitude
data_df = DataFrame(time=time)
for (i, ch) in enumerate(channel_names)
    data_df[!, ch] = data_matrix[:, i]
end

# Create layout
layout = create_layout(channel_names)

# Create ContinuousData
dat = ContinuousData(data_df, layout)

println("Original data stats:")
println("  Mean Fz: $(mean(dat.data.Fz))")
println("  Std Fz: $(std(dat.data.Fz))")

# Test the mutating version
println("\nTesting plot_databrowser! (mutating version)...")
fig, ax, modified_dat = plot_databrowser!(dat)

println("After plot_databrowser!:")
println("  Mean Fz: $(mean(modified_dat.data.Fz))")
println("  Std Fz: $(std(modified_dat.data.Fz))")
println("  Same object? $(dat === modified_dat)")

println("\nTest completed!")
