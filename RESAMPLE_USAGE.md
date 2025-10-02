# Resample Function Usage Guide

## Overview

The `resample` function downsamples EEG data by a specified integer factor, reducing the sampling rate while preserving all metadata including triggers, trial information, and other columns. It works with all three data types: `ContinuousData`, `EpochData`, and `ErpData`.

## Key Features

- **Works with all data types**: Continuous, Epoched, and ERP data
- **Factor-based downsampling**: Specify integer factors (e.g., 2, 4, 8)
- **Metadata preservation**: All columns preserved including triggers
- **Batch processing**: Process multiple participant files at once
- **Mutating and non-mutating**: Choose between in-place or copy operations

## Basic Usage

### Single File Processing

#### Continuous Data
```julia
using eegfun

# Load continuous data at 512 Hz
data = load("participant_1_continuous.jld2", "data")

# Downsample to 256 Hz (factor of 2)
data_256hz = resample(data, 2)

# Original data unchanged
@assert data.sample_rate == 512
@assert data_256hz.sample_rate == 256
```

#### Epoched Data
```julia
# Load epoched data at 512 Hz
epochs = load("participant_1_epochs.jld2", "epochs")

# Downsample to 128 Hz (factor of 4)
epochs_128hz = resample(epochs, 4)

# All epochs are resampled
@assert epochs_128hz.sample_rate == 128
```

#### ERP Data
```julia
# Load ERP data at 1024 Hz
erp = load("participant_1_erp.jld2", "erp")

# Downsample to 256 Hz (factor of 4)
erp_256hz = resample(erp, 4)

# n_epochs preserved
@assert erp_256hz.n_epochs == erp.n_epochs
```

### In-Place Resampling

```julia
# Modify data in-place to save memory
data = load("continuous.jld2", "data")

resample!(data, 2)  # Modifies data directly

# data.sample_rate is now 256 Hz
```

## Batch Processing

### Basic Batch Resampling
```julia
# Resample all continuous data files by factor of 2
resample("continuous", 2)

# Output directory: ./resampled_by_2_continuous/
```

### Specific Participants
```julia
# Resample only participants 1, 2, and 3
resample("epochs", 2, participants = [1, 2, 3])
```

### Custom Input/Output Directories
```julia
resample("erp", 4,
        input_dir = "/data/study1/preprocessed",
        output_dir = "/data/study1/resampled_256hz")
```

### Full Batch Example
```julia
# Resample all participants in a study
resample("continuous", 2,
        input_dir = "/data/flanker_task/filtered",
        participants = 1:30,
        output_dir = "/data/flanker_task/resampled")
```

## Choosing a Downsampling Factor

### Common Factors

| Original Rate | Factor | New Rate | Use Case |
|--------------|--------|----------|----------|
| 512 Hz | 2 | 256 Hz | Standard ERP analysis |
| 1024 Hz | 4 | 256 Hz | High-resolution to standard |
| 256 Hz | 2 | 128 Hz | Time-frequency analysis |
| 500 Hz | 5 | 100 Hz | Slow wave analysis |

### Guidelines

1. **Powers of 2 work best**: Factors like 2, 4, 8 (512→256→128→64)
2. **Must result in integer sample rates**: 500 Hz ÷ 2 = 250 Hz ✓, 500 Hz ÷ 3 = 166.67 Hz ✗
3. **Consider Nyquist frequency**: New rate should be > 2× highest frequency of interest
4. **Match target analysis**: Some algorithms require specific sampling rates

## Metadata Preservation

### All Columns Are Preserved

The resample function keeps **all** DataFrame columns:

```julia
# Data with triggers and metadata
data = DataFrame(
    time = ...,
    C3 = ...,
    C4 = ...,
    trigger = [0, 0, 1, 0, 2, 0, ...],  # Trigger codes
    trial = [1, 1, 1, 1, 1, 1, ...],    # Trial numbers
    condition = ["A", "A", "A", ...],    # Conditions
    rt = [0.5, 0.5, 0.5, ...]           # Response times
)

resampled = resample(continuous_data, 2)

# All columns preserved
@assert "trigger" in names(resampled.data)
@assert "trial" in names(resampled.data)
@assert "condition" in names(resampled.data)
@assert "rt" in names(resampled.data)
```

### Trigger Preservation

Triggers are preserved in the downsampled data:

```julia
# Original: triggers at samples 100 and 500
# Factor 2: triggers may appear at samples 50 and 250 (if they fall on kept samples)

# Check trigger preservation
original_triggers = sum(data.data.trigger .> 0)
resampled_triggers = sum(resampled.data.trigger .> 0)

@info "Original triggers: $original_triggers"
@info "Resampled triggers: $resampled_triggers"
```

**Note**: If a trigger falls on a sample that's not kept during downsampling, it won't appear in the resampled data. For critical triggers, ensure they're positioned on samples that will be kept, or use factor 1 (no resampling).

## Avoiding Aliasing

### The Problem

Downsampling without filtering can cause **aliasing** - high-frequency signals appearing as low-frequency artifacts.

### Solution: Filter Before Resampling

```julia
using eegfun

# Load data at 512 Hz
data = load("participant_1_continuous.jld2", "data")

# Filter to avoid aliasing
# For downsampling to 256 Hz, low-pass filter at ~120 Hz (< Nyquist/2)
filter!(data, hp = 0.1, lp = 120.0)

# Now safely downsample
resample!(data, 2)

# Result: Clean 256 Hz data without aliasing
```

### Recommended Workflow

```julia
# 1. Load data
data = load("continuous.jld2", "data")  # 512 Hz

# 2. High-pass filter
filter!(data, hp = 0.1, lp = 0.0)

# 3. Low-pass filter (below new Nyquist)
# New rate will be 256 Hz, so Nyquist = 128 Hz
# Filter at ~120 Hz to be safe
filter!(data, hp = 0.0, lp = 120.0)

# 4. Downsample
resample!(data, 2)

# 5. Save
save("continuous_256hz.jld2", "data", data)
```

## Complete Processing Pipeline

### Example: Preprocessing with Resampling

```julia
using eegfun

# Load raw data at 1024 Hz
data = load("participant_1_raw.jld2", "data")

# 1. High-pass filter
filter!(data, hp = 0.1, lp = 0.0)

# 2. Rereference
rereference!(data, :avg)

# 3. Low-pass filter before downsampling
filter!(data, hp = 0.0, lp = 120.0)

# 4. Downsample to 256 Hz
resample!(data, 4)  # 1024 / 4 = 256

# 5. Epoch
conditions = [
    EpochCondition(
        name = "target",
        trigger_sequences = [[1]],
    ),
    EpochCondition(
        name = "standard",
        trigger_sequences = [[2]],
    )
]

epochs = epoch_data(data, conditions, window = (-0.2, 0.8))

# 6. Baseline correction
baseline!(epochs, window = (-0.2, 0.0))

# 7. Artifact rejection
clean_epochs = reject_epochs_automatic(epochs, 2.5)

# 8. Average to ERP
erp = analyse(clean_epochs)

# Save final ERP
save("participant_1_erp_256hz.jld2", "erp", erp)
```

## Batch Processing Workflow

### Full Study Resampling

```julia
using eegfun

# Define study parameters
input_dir = "/data/study1/filtered"
output_dir = "/data/study1/resampled_256hz"
participants = 1:30

# Batch resample all continuous data
resample("continuous", 2,
        input_dir = input_dir,
        output_dir = output_dir,
        participants = participants)

# Log file saved to: /data/study1/resampled_256hz/resample.log
```

### Processing Multiple Stages

```julia
# Stage 1: Resample continuous data
resample("continuous", 2, input_dir = "data/filtered")

# Stage 2: Resample epochs (from different directory)
resample("epochs", 2, input_dir = "data/epoched")

# Stage 3: Resample ERPs
resample("erp", 2, input_dir = "data/erps")
```

## Performance Considerations

### Memory Usage

- **Non-mutating (`resample`)**: Creates a copy, uses ~2× memory
- **Mutating (`resample!`)**: Modifies in-place, memory-efficient

```julia
# Memory-efficient for large datasets
data = load("large_file.jld2", "data")
resample!(data, 2)  # Modifies in-place
save("large_file_resampled.jld2", "data", data)
```

### Speed

Downsampling is fast - it's simple decimation (taking every nth sample):

```julia
# Typical speed: ~0.1-1 second per file
# 1 million samples at 1024 Hz → 250k samples at 256 Hz: ~0.5 seconds
```

## Common Patterns

### 1. Standard ERP Analysis
```julia
# 512 Hz → 256 Hz for typical ERP analysis
resample("continuous", 2)
```

### 2. High-Resolution to Standard
```julia
# 1024 Hz → 256 Hz for standard processing
resample("continuous", 4)
```

### 3. Time-Frequency Analysis
```julia
# 512 Hz → 128 Hz for faster TF decomposition
resample("epochs", 4)
```

### 4. Matching Across Studies
```julia
# Study 1: 512 Hz data
resample("continuous", 2, input_dir = "study1")  # → 256 Hz

# Study 2: 1024 Hz data
resample("continuous", 4, input_dir = "study2")  # → 256 Hz

# Now both studies at 256 Hz for combined analysis
```

## Troubleshooting

### Error: "Sample rate not evenly divisible"
```julia
# Problem: 500 Hz ÷ 3 = 166.67 Hz (not integer)
resample(data, 3)  # ✗ Error

# Solution: Use factor that results in integer
resample(data, 2)  # ✓ 500 Hz → 250 Hz
resample(data, 5)  # ✓ 500 Hz → 100 Hz
```

### Missing Triggers After Resampling
```julia
# If triggers fall on non-kept samples, they disappear

# Check trigger alignment before resampling
trigger_indices = findall(data.data.trigger .> 0)
@info "Trigger indices: $trigger_indices"

# Ensure important triggers fall on kept samples (multiples of factor)
```

### Aliasing Artifacts
```julia
# Always filter before downsampling!

# Before
filter!(data, hp = 0.1, lp = 120.0)  # Low-pass at < Nyquist/2

# Then
resample!(data, 2)  # Safe downsampling
```

## Summary

1. **Choose appropriate factor**: Must result in integer sample rate
2. **Filter first**: Avoid aliasing by low-pass filtering
3. **Check triggers**: Ensure critical triggers are preserved
4. **Use batch processing**: Efficient for multiple participants
5. **Consider memory**: Use `resample!` for large datasets

The resample function makes it easy to adjust sampling rates while preserving all your data's metadata and structure!

