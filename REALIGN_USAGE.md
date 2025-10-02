# Epoch Realignment - Usage Guide

## Overview

The `realign()` function takes stimulus-locked (or any reference-locked) epoched EEG data and realigns it to a different time point specified in a DataFrame column. This is particularly useful for creating response-locked waveforms from stimulus-locked epochs, which is essential for response-locked LRP analysis.

## Why Realignment?

When analyzing event-related potentials, you typically epoch your data around stimulus onset (stimulus-locked). However, you may want to analyze activity time-locked to the response instead (response-locked). Realignment allows you to:

1. **Maintain trial identity**: Use the exact same trials for both stimulus-locked and response-locked analyses
2. **Preserve data**: Keep all preprocessing (artifact rejection, filtering, etc.) that was done on stimulus-locked epochs
3. **Ensure comparability**: Guarantee that S-locked and R-locked LRPs come from identical trial sets

## Basic Usage

### Single Dataset Realignment

```julia
using eegfun
using JLD2

# Load stimulus-locked epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Realign to response times (stored in :rt column)
# Non-mutating version (original data unchanged)
realigned_epochs = realign(epochs, :rt)

# Now time=0 corresponds to the response for each trial
```

### In-Place Realignment

```julia
# Mutating version (modifies original data)
realign!(epochs, :rt)

# epochs is now response-locked
```

### Different Realignment Columns

The realignment column can be named anything, as long as it exists in your epoch DataFrames:

```julia
# Realign to saccade onset
realigned_epochs = realign(epochs, :saccade_time)

# Realign to any custom event
realigned_epochs = realign(epochs, :button_press_time)

# Realign to stimulus offset
realigned_epochs = realign(epochs, :stimulus_offset)
```

## Requirements

### Data Structure

1. **EpochData format**: Input must be `EpochData` (not averaged ERP data)
2. **Realignment column**: Must exist in all epoch DataFrames
3. **Constant values**: Realignment time must be constant within each epoch
4. **Finite values**: No NaN or Inf values in realignment column
5. **Sufficient data**: Original epochs must be long enough to accommodate all realignment times

### Realignment Column Values

The realignment column should contain times **relative to the current time zero**:

```julia
# Example: If your epochs are from -0.5s to 2.0s (stimulus-locked)
# And RT is 0.6s after stimulus
# Then the :rt column should contain 0.6

# After realignment:
# Time 0 will correspond to the response
# The epoch will span approximately -0.6s to 1.4s (response-locked)
```

## How It Works

The realignment process involves three steps:

### Step 1: Individual Epoch Realignment

Each epoch's time vector is shifted so that the realignment event becomes time zero:

```julia
# For each epoch:
# new_time = old_time - realignment_time

# Example:
# Original time: [-0.5, -0.4, ..., 1.9, 2.0]
# RT: 0.6s
# New time: [-1.1, -1.0, ..., 1.3, 1.4]
```

### Step 2: Find Common Time Window

After realignment, different trials have different time ranges (because RTs vary). The function finds the time window that is valid for **all** trials:

```julia
# Trial 1: RT = 0.4s → After realignment: -0.9 to 1.6s
# Trial 2: RT = 0.8s → After realignment: -1.3 to 1.2s
# Trial 3: RT = 0.6s → After realignment: -1.1 to 1.4s

# Common window: max(-0.9, -1.3, -1.1) to min(1.6, 1.2, 1.4)
#              = -0.9 to 1.2s
```

### Step 3: Crop to Common Window

All epochs are cropped to the common time window, ensuring all epochs have identical time vectors.

## Complete Example Workflows

### Response-Locked LRP Analysis

```julia
using eegfun
using JLD2

# 1. Load stimulus-locked epochs
epochs = load("participant_05_epochs_cleaned.jld2", "epochs")

println("Original epochs:")
println("  Time range: $(minimum(epochs.data[1].time)) to $(maximum(epochs.data[1].time)) s")
println("  N samples: $(nrow(epochs.data[1]))")
println("  N epochs: $(length(epochs.data))")

# 2. Realign to response
realigned_epochs = realign(epochs, :rt)

println("\nRealigned epochs:")
println("  Time range: $(minimum(realigned_epochs.data[1].time)) to $(maximum(realigned_epochs.data[1].time)) s")
println("  N samples: $(nrow(realigned_epochs.data[1]))")
println("  N epochs: $(length(realigned_epochs.data))")

# 3. Average to create response-locked ERPs
erps_response_locked = average_epochs(realigned_epochs)

# 4. Calculate response-locked LRP
lrp_r = lrp(erps_response_locked[1], erps_response_locked[2])

# 5. Save
save("participant_05_epochs_response_locked.jld2", "epochs", realigned_epochs)
save("participant_05_erps_response_locked.jld2", "erps", erps_response_locked)
save("participant_05_lrp_response_locked.jld2", "lrp", lrp_r)
```

### Batch Processing Multiple Participants

```julia
using eegfun

# Realign all participants' epochs to response times
realign("epochs_cleaned", :rt,
        input_dir = "/data/study1",
        output_dir = "/data/study1/realigned")

# This processes all files matching "epochs_cleaned" pattern
# and saves realigned data to a new directory
```

### Comparing Stimulus-Locked vs Response-Locked

```julia
using eegfun
using JLD2
using CairoMakie

# Load both versions
erps_s = load("participant_05_erps_stimulus_locked.jld2", "erps")
erps_r = load("participant_05_erps_response_locked.jld2", "erps")

# Plot comparison
fig = Figure(size = (1000, 400))

# Stimulus-locked
ax1 = Axis(fig[1, 1], 
          xlabel = "Time from stimulus (s)",
          ylabel = "Amplitude (μV)",
          title = "Stimulus-Locked")
plot_erp!(ax1, erps_s[1], channels = [:C3, :C4])
vlines!(ax1, [0], color = :black, linestyle = :dash, label = "Stimulus")

# Response-locked
ax2 = Axis(fig[1, 2],
          xlabel = "Time from response (s)", 
          ylabel = "Amplitude (μV)",
          title = "Response-Locked")
plot_erp!(ax2, erps_r[1], channels = [:C3, :C4])
vlines!(ax2, [0], color = :red, linestyle = :dash, label = "Response")

display(fig)
```

## Advanced Usage

### Checking Data Before Realignment

```julia
using Statistics

# Check RT distribution
rts = [epoch.rt[1] for epoch in epochs.data]
println("RT range: $(minimum(rts)) to $(maximum(rts)) s")
println("RT mean: $(mean(rts)) s")
println("RT std: $(std(rts)) s")

# Check if original epochs are long enough
original_duration = maximum(epochs.data[1].time) - minimum(epochs.data[1].time)
println("Original epoch duration: $(original_duration) s")

# Estimate final epoch length after realignment
# This is a rough estimate
min_rt = minimum(rts)
max_rt = maximum(rts)
time_before_response = -(minimum(epochs.data[1].time) - max_rt)
time_after_response = maximum(epochs.data[1].time) - min_rt

println("Estimated response-locked epoch: $(time_before_response) to $(time_after_response) s")
```

### Handling Multiple Realignment Columns

If you have multiple possible realignment targets, you can process them separately:

```julia
# Realign to first button press
realigned_button1 = realign(epochs, :button1_time)

# Realign to second button press
realigned_button2 = realign(epochs, :button2_time)

# Save both versions
save("epochs_button1_locked.jld2", "epochs", realigned_button1)
save("epochs_button2_locked.jld2", "epochs", realigned_button2)
```

### Quality Control After Realignment

```julia
# Verify that all epochs have the same time vector
time_vectors_match = all([
    all(epochs.data[i].time .≈ epochs.data[1].time) 
    for i in 2:length(epochs.data)
])
println("All time vectors match: $time_vectors_match")

# Check realignment column is now zero
rt_is_zero = all([
    all(abs.(epoch.rt) .< 1e-10) 
    for epoch in epochs.data
])
println("RT column is zero: $rt_is_zero")

# Check final epoch length
final_length = nrow(epochs.data[1])
final_duration = maximum(epochs.data[1].time) - minimum(epochs.data[1].time)
println("Final epoch: $final_length samples, $(final_duration) s")
```

## Common Issues and Solutions

### Issue: "No common time window found"

**Cause**: Original epochs are too short to accommodate all realignment times.

**Solution**: Use longer epochs when initially creating your stimulus-locked data.

```julia
# When epoching, use a longer window:
# Instead of: epoch_window = (-0.5, 1.5)  # 2.0s total
# Use:        epoch_window = (-1.0, 2.5)  # 3.5s total

# This gives more flexibility for realignment
```

### Issue: "Realignment column has varying values"

**Cause**: The realignment column has different values within a single epoch.

**Solution**: Ensure the realignment column is constant within each epoch. This usually means the column was created correctly during epoching.

```julia
# Check if values vary within epochs
for (i, epoch) in enumerate(epochs.data)
    if !all(epoch.rt .≈ epoch.rt[1])
        println("Epoch $i has varying RT values!")
        println("Values: ", unique(epoch.rt))
    end
end
```

### Issue: Very short final epochs

**Cause**: Wide range of realignment times + insufficient original epoch length.

**Solution**: 
1. Use longer original epochs
2. Consider removing trials with very early or very late RTs
3. Accept shorter epochs if necessary

```julia
using Statistics

# Remove outlier RTs before realignment
rts = [epoch.rt[1] for epoch in epochs.data]
rt_mean = mean(rts)
rt_std = std(rts)

# Keep only trials within 2 SD
epochs_filtered = EpochData(
    [epoch for epoch in epochs.data if abs(epoch.rt[1] - rt_mean) < 2 * rt_std],
    epochs.layout,
    epochs.sample_rate,
    epochs.analysis_info
)

# Now realign
realigned = realign(epochs_filtered, :rt)
```

## Batch Processing

### Basic Batch Realignment

```julia
# Process all participants
realign("epochs_cleaned", :rt)

# Creates output directory: realigned_epochs_cleaned_rt/
```

### With Options

```julia
# Specific participants
realign("epochs_cleaned", :rt,
        participants = 1:20)

# Custom directories
realign("epochs_cleaned", :rt,
        input_dir = "/data/study1/epochs",
        output_dir = "/data/study1/response_locked")

# Full example
realign("epochs_cleaned", :rt,
        input_dir = "/data/study1",
        participants = 1:30,
        output_dir = "/data/study1/response_locked_epochs")
```

### Batch Processing Output

```
your_data_directory/
├── 1_epochs_cleaned.jld2
├── 2_epochs_cleaned.jld2
├── 3_epochs_cleaned.jld2
└── realigned_epochs_cleaned_rt/     ← NEW!
    ├── 1_epochs_cleaned.jld2        ← Realigned epochs
    ├── 2_epochs_cleaned.jld2        ← Realigned epochs
    ├── 3_epochs_cleaned.jld2        ← Realigned epochs
    └── realign.log                  ← Processing log
```

## Full Workflow: S-LRP and LRP-R

Here's a complete workflow for calculating both stimulus-locked LRP (S-LRP) and response-locked LRP (LRP-R) using the same trials:

```julia
using eegfun
using JLD2

participant = 5

# 1. Load stimulus-locked epochs
epochs_s = load("participant_$(participant)_epochs_cleaned.jld2", "epochs")

# 2. Calculate stimulus-locked ERPs and LRP
erps_s = average_epochs(epochs_s)
lrp_s = lrp(erps_s[1], erps_s[2])  # Assuming conditions 1=left, 2=right
save("participant_$(participant)_lrp_stimulus.jld2", "lrp", lrp_s)

# 3. Realign to response
epochs_r = realign(epochs_s, :rt)

# 4. Calculate response-locked ERPs and LRP
erps_r = average_epochs(epochs_r)
lrp_r = lrp(erps_r[1], erps_r[2])
save("participant_$(participant)_lrp_response.jld2", "lrp", lrp_r)

# 5. Plot both
using CairoMakie

fig = Figure(size = (1200, 400))

ax1 = Axis(fig[1, 1], 
          xlabel = "Time from stimulus (s)",
          ylabel = "LRP (μV)",
          title = "S-LRP")
lines!(ax1, lrp_s.data.time, lrp_s.data.C3, label = "C3")
vlines!(ax1, [0], color = :black, linestyle = :dash)

ax2 = Axis(fig[1, 2],
          xlabel = "Time from response (s)",
          ylabel = "LRP (μV)", 
          title = "LRP-R")
lines!(ax2, lrp_r.data.time, lrp_r.data.C3, label = "C3")
vlines!(ax2, [0], color = :red, linestyle = :dash)

display(fig)
```

## Performance Considerations

- Realignment is relatively fast (< 1 second for typical datasets)
- The cropping step ensures all epochs have identical length, which is required for averaging
- Memory usage is minimal as the function works with DataFrames efficiently

## References

This implementation is based on the approach used in FieldTrip's `ft_redefinetrial` function with response-locking functionality.

## See Also

- LRP calculation: See `LRP_USAGE.md`
- Epoching: See documentation for `epochs()` function
- Averaging: See documentation for `average_epochs()` function

