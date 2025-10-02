# Data Mirroring - Usage Guide

## Overview

The `mirror()` function adds time-reversed (mirrored) copies of data before and/or after epochs. This is primarily used to reduce edge artifacts when filtering EEG data, as the mirrored sections create smooth transitions at epoch boundaries.

## Why Mirror Data?

When filtering epoched data, edge artifacts can occur at the beginning and end of epochs due to the abrupt transitions. Mirroring solves this by:

1. **Extending epochs smoothly**: Reversed data maintains signal continuity
2. **Reducing edge effects**: Filter transients occur in the mirrored sections, not the original data
3. **Preserving original data**: After filtering, simply remove the mirrored sections

## Basic Usage

### With Epoched Data

```julia
using eegfun, JLD2

# Load epochs
epochs = load("participant_1_epochs.jld2", "epochs")

# Mirror on both sides (recommended)
mirror!(epochs, :both)

# Apply processing (e.g., filtering)
filter!(epochs, 0.1, 30.0)

# Remove mirrored sections
unmirror!(epochs, :both)

# Continue with analysis
save("participant_1_epochs_filtered.jld2", "epochs", epochs)
```

### With ERP Data

```julia
# Load averaged ERP
erp = load("participant_1_erp.jld2", "erp")

# Mirror, process, unmirror
mirror!(erp, :both)
filter!(erp, 0.1, 30.0)
unmirror!(erp, :both)
```

## How It Works

### Original Data
```
Time:  [-0.2, -0.1, 0.0, 0.1, 0.2]
Data:  [  1,    2,   3,   4,   5 ]
```

### :pre (Mirror Before)
```
Time:  [-0.6, -0.5, -0.4, -0.3, -0.1, 0.0, 0.1, 0.2]
Data:  [  5,    4,    3,    2,    2,   3,   4,   5 ]
        └────── mirrored ──────┘    └─── original ──┘
```

### :post (Mirror After)
```
Time:  [-0.2, -0.1, 0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
Data:  [  1,    2,   3,   4,   4,   3,   2,   1 ]
        └─── original ──┘    └──── mirrored ─────┘
```

### :both (Mirror Both Sides)
```
Time:  [-0.6, -0.5, -0.4, -0.3, -0.1, 0.0, 0.1, 0.3, 0.4, 0.5, 0.6]
Data:  [  5,    4,    3,    2,    2,   3,   4,   4,   3,   2,   1 ]
        └────── pre ───────┘    └ original ┘  └──── post ────┘
```

## Mirroring Sides

### :both (Default - Recommended)
```julia
mirror!(epochs, :both)
# Mirrors on both sides
# Use when: You need maximum protection from edge effects
# Results in: ~3× epoch length
```

### :pre (Mirror Before Only)
```julia
mirror!(epochs, :pre)
# Mirrors only at the beginning
# Use when: Edge effects mainly at epoch start
# Results in: ~2× epoch length
```

### :post (Mirror After Only)
```julia
mirror!(epochs, :post)
# Mirrors only at the end
# Use when: Edge effects mainly at epoch end
# Results in: ~2× epoch length
```

## Complete Workflows

### Standard Filtering Workflow

```julia
using eegfun, JLD2

# Load epochs
epochs = load("participant_1_epochs.jld2", "epochs")

println("Original epoch length: $(nrow(epochs.data[1])) samples")

# Step 1: Mirror data
mirror!(epochs, :both)
println("After mirroring: $(nrow(epochs.data[1])) samples")

# Step 2: Apply filter
filter!(epochs, 0.1, 30.0)  # High-pass 0.1 Hz, low-pass 30 Hz

# Step 3: Remove mirrors
unmirror!(epochs, :both)
println("After unmirroring: $(nrow(epochs.data[1])) samples")

# Save filtered data
save("participant_1_epochs_filtered.jld2", "epochs", epochs)
```

### Non-Mutating Version

```julia
# If you want to keep original data
epochs_original = load("participant_1_epochs.jld2", "epochs")

# Create mirrored copy
epochs_mirrored = mirror(epochs_original, :both)

# Filter the mirrored copy
filter!(epochs_mirrored, 0.1, 30.0)

# Unmirror
epochs_filtered = unmirror(epochs_mirrored, :both)

# Now you have both original and filtered
```

### Batch Processing

```julia
using eegfun, JLD2

participants = 1:20

for participant in participants
    filename = "participant_$(participant)_epochs.jld2"
    
    # Load
    epochs = load(filename, "epochs")
    
    # Mirror → Filter → Unmirror
    mirror
!(epochs, :both)
    filter!(epochs, 0.1, 30.0)
    unmirror
!(epochs, :both)
    
    # Save
    output_file = "participant_$(participant)_epochs_filtered.jld2"
    save(output_file, "epochs", epochs)
    
    println("✓ Processed participant $participant")
end
```

### With Baseline Correction

```julia
# Recommended order: Mirror → Baseline → Filter → Unmirror

epochs = load("participant_1_epochs.jld2", "epochs")

# Mirror first
mirror!(epochs, :both)

# Baseline correction (works on mirrored data)
baseline!(epochs, IntervalTime(-0.2, 0.0))

# Filter
filter!(epochs, 0.1, 30.0)

# Unmirror
unmirror!(epochs, :both)

# Save
save("participant_1_epochs_processed.jld2", "epochs", epochs)
```

### ERP Processing

```julia
# For averaged ERP data
erp = load("grand_average.jld2", "grand_avg")

# Mirror, filter, unmirror
mirror!(erp, :both)
filter!(erp, 0.5, 30.0)
unmirror!(erp, :both)

# Save
save("grand_average_filtered.jld2", "grand_avg", erp)
```

## Important Notes

### 1. Always Unmirror!

```julia
# ✓ CORRECT
mirror!(epochs, :both)
filter!(epochs, 0.1, 30.0)
unmirror!(epochs, :both)  # Don't forget!

# ✗ WRONG - leaves mirrored sections
mirror!(epochs, :both)
filter!(epochs, 0.1, 30.0)
# Oops! Data is now 3× longer than it should be
```

### 2. Match Mirror and Unmirror Sides

```julia
# ✓ CORRECT
mirror!(epochs, :both)
# ... processing ...
unmirror!(epochs, :both)  # Same side

# ✗ WRONG - sides don't match
mirror!(epochs, :both)
# ... processing ...
unmirror!(epochs, :pre)  # Wrong! Should be :both
```

### 3. Memory Considerations

Mirroring increases memory usage:
- `:pre` or `:post`: ~2× memory
- `:both`: ~3× memory

For large datasets:
```julia
# Process one condition at a time
for cond in conditions
    epochs_cond = load("participant_1_cond_$(cond)_epochs.jld2", "epochs")
    mirror
!(epochs_cond, :both)
    filter!(epochs_cond, 0.1, 30.0)
    unmirror
!(epochs_cond, :both)
    save("participant_1_cond_$(cond)_filtered.jld2", "epochs", epochs_cond)
end
```

## When to Use Mirroring

### Use Mirroring When:
- ✓ Applying filters (especially high-pass filters)
- ✓ Performing FFT or wavelet analysis
- ✓ Any processing sensitive to edge effects
- ✓ Short epochs where edge effects are significant

### Skip Mirroring When:
- ✗ Epochs are very long (edge effects negligible)
- ✗ Only performing amplitude-based operations
- ✗ Memory is severely constrained
- ✗ Processing speed is critical

## Advanced Usage

### Custom Filtering Pipeline

```julia
function custom_filter_with_mirror!(epochs, hp_freq, lp_freq)
    # Mirror
    mirror
!(epochs, :both)
    
    # Apply high-pass
    if hp_freq > 0
        filter!(epochs, hp_freq, nothing)
    end
    
    # Apply low-pass
    if lp_freq < Inf
        filter!(epochs, nothing, lp_freq)
    end
    
    # Unmirror
    unmirror
!(epochs, :both)
    
    return epochs
end

# Use it
epochs = load("participant_1_epochs.jld2", "epochs")
custom_filter_with_mirror!(epochs, 0.1, 30.0)
```

### Comparing With and Without Mirroring

```julia
using CairoMakie

epochs_original = load("participant_1_epochs.jld2", "epochs")

# Filter without mirroring
epochs_no_mirror = copy(epochs_original)
filter!(epochs_no_mirror, 0.1, 30.0)

# Filter with mirroring
epochs_with_mirror = copy(epochs_original)
mirror!(epochs_with_mirror, :both)
filter!(epochs_with_mirror, 0.1, 30.0)
unmirror!(epochs_with_mirror, :both)

# Compare edge effects
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "μV", 
         title = "Edge Effects: With vs Without Mirroring")

# Plot first 50ms to see edge effects
erp_no_mirror = average_epochs(epochs_no_mirror)
erp_with_mirror = average_epochs(epochs_with_mirror)

time_mask = erp_no_mirror.data.time .< -0.15  # Look at epoch start
lines!(ax, erp_no_mirror.data.time[time_mask], erp_no_mirror.data.Cz[time_mask],
      label = "Without mirroring", color = :red)
lines!(ax, erp_with_mirror.data.time[time_mask], erp_with_mirror.data.Cz[time_mask],
      label = "With mirroring", color = :green)

axislegend(ax)
display(fig)
```

### Conditional Mirroring

```julia
function smart_filter!(epochs, hp_freq, lp_freq; use_mirror = :auto)
    # Decide whether to use mirroring based on filter parameters
    if use_mirror == :auto
        # Use mirroring for aggressive high-pass filters
        use_mirror = (hp_freq !== nothing && hp_freq > 0.05)
    end
    
    if use_mirror
        mirror
    !(epochs, :both)
    end
    
    filter!(epochs, hp_freq, lp_freq)
    
    if use_mirror
        unmirror
    !(epochs, :both)
    end
    
    return epochs
end

# Use automatic decision
smart_filter!(epochs, 0.1, 30.0)  # Will use mirroring
smart_filter!(epochs, nothing, 30.0)  # Won't use mirroring
```

## Verification

### Check Mirror Worked Correctly

```julia
# Before mirroring
original_length = nrow(epochs.data[1])
original_time_start = epochs.data[1].time[1]
original_time_end = epochs.data[1].time[end]

# Mirror
mirror!(epochs, :both)

# Check lengths
mirrored_length = nrow(epochs.data[1])
@assert mirrored_length ≈ 3 * original_length - 4  # Approximate

# Unmirror
unmirror!(epochs, :both)

# Should be back to original
@assert nrow(epochs.data[1]) == original_length
@assert epochs.data[1].time[1] ≈ original_time_start
@assert epochs.data[1].time[end] ≈ original_time_end

println("✓ Mirroring/unmirroring worked correctly")
```

## Troubleshooting

### Issue: Epoch length wrong after unmirroring

**Cause**: Mirror and unmirror sides don't match

**Solution**: 
```julia
# Always use matching sides
mirror!(epochs, :both)
# ... processing ...
unmirror!(epochs, :both)  # Must match!
```

### Issue: Time vector doesn't make sense

**Cause**: Forgot to unmirror

**Solution**:
```julia
# Check if data is still mirrored
println("Current epoch length: $(nrow(epochs.data[1]))")
println("Expected original length: ???")

# If too long, unmirror with appropriate side
unmirror!(epochs, :both)
```

### Issue: Edge artifacts still present

**Possible causes**:
1. Mirror wasn't applied before filtering
2. Epochs too short for effective mirroring
3. Filter is too aggressive

**Solutions**:
```julia
# Ensure mirroring is applied
mirror!(epochs, :both)
filter!(epochs, hp_freq, lp_freq)
unmirror!(epochs, :both)

# Or use longer epochs if possible
```

## Performance Considerations

Mirroring operations are fast:
- Mirroring: ~1-2 ms per epoch
- Unmirroring: ~1-2 ms per epoch
- Memory: Increases by 2-3× during mirroring

The processing time increase is negligible compared to filtering benefits.

## References

This implementation is based on standard practice in EEG analysis for reducing edge artifacts during filtering. The approach is similar to what's used in FieldTrip and EEGLAB.

## See Also

- Filtering: `filter()` function
- Baseline correction: `baseline()` function  
- Epoch extraction: `extract_epochs()` function
- ERP averaging: `average_epochs()` function

