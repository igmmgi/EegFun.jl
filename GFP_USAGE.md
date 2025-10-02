# Global Field Power (GFP) - Usage Guide

## Overview

The `gfp()` function calculates Global Field Power, a reference-independent measure of response strength computed as the standard deviation across all channels at each time point. This package also provides Global Dissimilarity calculations, which measure the rate of topographic change.

## Basic Usage

### Single ERP dataset

```julia
using eegfun
using JLD2

# Load your ERP data
erp_data = load("participant_1_erps.jld2", "erps")[1]

# Calculate GFP using all channels
gfp_result = gfp(erp_data)

# The result is a DataFrame with columns: time, gfp, and metadata
println("Time points: ", size(gfp_result, 1))
println("GFP at t=0: ", gfp_result.gfp[gfp_result.time .== 0.0][1])
```

### Normalized GFP (0-100%)

```julia
# Calculate GFP normalized to 0-100% range
# (useful for comparing across datasets or conditions)
gfp_result = gfp(erp_data, normalize = true)
```

### With Channel Selection

```julia
# Calculate GFP using specific channels
gfp_result = gfp(erp_data, 
                 channel_selection = channels([:C3, :C4, :Cz, :CP3, :CP4, :CPz]))

# Use pattern matching for channel selection
gfp_result = gfp(erp_data,
                 channel_selection = channels(x -> startswith.(string.(x), "C")))
```

### Multiple Conditions

```julia
# Load ERP data for multiple conditions
erps = load("participant_1_erps.jld2", "erps")

# Calculate GFP for all conditions
gfp_results = gfp(erps, normalize = true)

# gfp_results is a Vector{DataFrame}, one for each condition
for (i, gfp_data) in enumerate(gfp_results)
    println("Condition $i peak GFP: ", maximum(gfp_data.gfp))
end
```

## Global Dissimilarity

Global Dissimilarity measures how quickly the topographic distribution changes over time.

```julia
# Calculate Global Dissimilarity
gd_result = global_dissimilarity(erp_data)

# With normalization
gd_result = global_dissimilarity(erp_data, normalize = true)

# The result has columns: time, dissimilarity, and metadata
```

## Combined GFP and Dissimilarity

```julia
# Calculate both metrics at once (more efficient)
result = gfp_and_dissimilarity(erp_data, normalize = true)

# Result has columns: time, gfp, dissimilarity, and metadata
time_vector = result.time
gfp_values = result.gfp
dissimilarity_values = result.dissimilarity
```

## Visualization

### Basic GFP Plot

```julia
using CairoMakie

# Plot GFP directly from ERP data
plot_gfp(erp_data)

# Plot normalized GFP
plot_gfp(erp_data, normalize = true)
```

### Plot with Individual Channel Traces

```julia
# Show ERP traces above the GFP
plot_gfp(erp_data, show_erp_traces = true)
```

### Plot with Global Dissimilarity

```julia
# Include Global Dissimilarity panel
plot_gfp(erp_data, show_dissimilarity = true)
```

### Plot All Panels (Similar to MATLAB calculateMyGFP)

```julia
# Three-panel plot: ERP traces, GFP, and dissimilarity
plot_gfp(erp_data, 
         show_erp_traces = true, 
         show_dissimilarity = true,
         normalize = true)
```

### Plot Multiple Conditions

```julia
# Load multiple conditions
erps = load("participant_1_erps.jld2", "erps")

# Plot all conditions on same axes
plot_gfp(erps, normalize = true)
```

### Custom Styling

```julia
# Customize plot appearance
plot_gfp(erp_data,
         color = :blue,
         linewidth = 3,
         xlim = (-0.2, 0.8),
         title = "Global Field Power - Participant 1")
```

### Plot Pre-computed Results

```julia
# Calculate GFP separately
gfp_result = gfp(erp_data, normalize = true)

# Plot later (useful for saving/loading results)
plot_gfp(gfp_result)

# Or with both metrics
result = gfp_and_dissimilarity(erp_data, normalize = true)
plot_gfp(result, show_dissimilarity = true)
```

## Complete Example Workflow

### Single Participant Analysis

```julia
using eegfun
using JLD2
using CairoMakie

# Load ERP data
erp_data = load("participant_05_erps.jld2", "erps")[1]

# Calculate GFP and dissimilarity
result = gfp_and_dissimilarity(erp_data, normalize = true)

# Plot
fig = plot_gfp(erp_data, 
               show_erp_traces = true,
               show_dissimilarity = true,
               normalize = true)

# Save plot
save("participant_05_gfp.png", fig)

# Save GFP data
save("participant_05_gfp_data.jld2", "gfp_data", result)
```

### Multiple Participants

```julia
using eegfun
using JLD2
using Statistics

participants = 1:20

# Calculate GFP for all participants
all_gfp = []
for participant in participants
    filename = "participant_$(participant)_erps.jld2"
    
    if !isfile(filename)
        @warn "Skipping participant $participant (file not found)"
        continue
    end
    
    erps = load(filename, "erps")
    gfp_result = gfp(erps[1], normalize = true)  # Condition 1
    push!(all_gfp, gfp_result)
    
    println("✓ Processed participant $participant")
end

# Compute grand average GFP
# Extract GFP values from all participants
gfp_matrix = hcat([df.gfp for df in all_gfp]...)
grand_avg_gfp = vec(mean(gfp_matrix, dims = 2))

# Plot grand average
using DataFrames
grand_avg_df = DataFrame(
    time = all_gfp[1].time,
    gfp = grand_avg_gfp,
    condition = fill(1, length(grand_avg_gfp)),
    condition_name = fill("grand_average", length(grand_avg_gfp))
)

plot_gfp(grand_avg_df)
```

### Comparing Conditions

```julia
using eegfun
using JLD2
using CairoMakie

# Load data with multiple conditions
erps = load("participant_05_erps.jld2", "erps")

# Calculate GFP for each condition
gfp_cond1 = gfp(erps[1], normalize = true)
gfp_cond2 = gfp(erps[2], normalize = true)

# Plot both conditions
fig = Figure()
ax = Axis(fig[1, 1],
          xlabel = "Time (s)",
          ylabel = "GFP (%)",
          title = "GFP Comparison")

lines!(ax, gfp_cond1.time, gfp_cond1.gfp, label = "Condition 1", color = :blue)
lines!(ax, gfp_cond2.time, gfp_cond2.gfp, label = "Condition 2", color = :red)
vlines!(ax, [0.0], color = :black, linestyle = :dash)
axislegend(ax)

display(fig)
```

## Interpreting Results

### Global Field Power (GFP)

- **High GFP**: Indicates strong, synchronized activity across channels
- **Low GFP**: Indicates weak or desynchronized activity
- **GFP Peaks**: Often correspond to major ERP components (e.g., P1, N1, P3)
- **Reference-independent**: Unlike individual channels, GFP is not affected by reference choice

### Global Dissimilarity

- **High Dissimilarity**: Indicates rapid changes in topographic distribution
- **Dissimilarity Peaks**: May indicate transitions between different brain states or ERP components
- **Use for segmentation**: Can help identify stable vs. transitional periods in the ERP

### Normalization

Normalizing to 0-100% is useful for:
- Comparing across participants with different overall amplitudes
- Comparing across different electrode montages
- Statistical comparisons where absolute amplitude differences are not of interest

## References

- Lehmann, D. & Skrandies, W. (1980). Reference-free identification of components of checkerboard-evoked multichannel potential fields. *Electroencephalography and Clinical Neurophysiology*, 48, 609-621.

- Lehmann, D. & Skrandies, W. (1984). Spatial analysis of evoked potentials in man--a review. *Progress in Neurobiology*, 23, 227-250.

- Skrandies, W. (1990). Global Field Power and Topographic Similarity. *Brain Topography*, 3(1), 137-141.

## Common Patterns

### Finding GFP Peaks

```julia
using eegfun

# Calculate GFP
gfp_result = gfp(erp_data, normalize = true)

# Find peaks in specific time window
time_window = (gfp_result.time .>= 0.1) .& (gfp_result.time .<= 0.3)
peak_gfp = maximum(gfp_result.gfp[time_window])
peak_time = gfp_result.time[time_window][argmax(gfp_result.gfp[time_window])]

println("Peak GFP: $peak_gfp at time $peak_time s")
```

### Exporting for Statistical Analysis

```julia
using CSV

# Calculate GFP for all participants and conditions
results = DataFrame(
    participant = Int[],
    condition = Int[],
    mean_gfp = Float64[],
    peak_gfp = Float64[],
    peak_latency = Float64[]
)

for participant in 1:20
    erps = load("participant_$(participant)_erps.jld2", "erps")
    
    for (cond_idx, erp) in enumerate(erps)
        gfp_result = gfp(erp, normalize = true)
        
        # Extract metrics in time window of interest
        time_window = (gfp_result.time .>= 0.0) .& (gfp_result.time .<= 0.5)
        
        push!(results, (
            participant = participant,
            condition = cond_idx,
            mean_gfp = mean(gfp_result.gfp[time_window]),
            peak_gfp = maximum(gfp_result.gfp[time_window]),
            peak_latency = gfp_result.time[time_window][argmax(gfp_result.gfp[time_window])]
        ))
    end
end

# Save for statistical analysis
CSV.write("gfp_metrics.csv", results)
```

## Troubleshooting

### Error: "No channels selected"

Make sure your channel selection actually matches channels in your data:

```julia
# Check available channels
println(channel_labels(erp_data))

# Use correct channel selection
gfp_result = gfp(erp_data, channel_selection = channels([:Fz, :Cz, :Pz]))
```

### GFP is constant (flat line)

This can happen if all channels have identical values. Check your data:

```julia
# Verify channels have different values
using Statistics
channel_vars = [var(erp_data.data[!, ch]) for ch in channel_labels(erp_data)]
println("Channel variances: ", channel_vars)
```

### Normalization produces NaN

This occurs when GFP has zero range (constant value). The function will skip normalization with a warning.

