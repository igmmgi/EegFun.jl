# Automatic Epoch Rejection - Usage Guide

## Overview

The `reject_epochs_automatic()` function automatically rejects epochs containing artifacts based on statistical criteria. It calculates six metrics (variance, max, min, absolute max, range, kurtosis) for each epoch and rejects those that exceed a z-score threshold for any metric.

## Why Automatic Rejection?

Manual artifact rejection can be time-consuming and subjective. Automatic rejection:

1. **Saves time**: Process hundreds of epochs in seconds
2. **Objective**: Uses statistical criteria rather than visual inspection
3. **Consistent**: Same criteria applied across all participants
4. **Documented**: Provides detailed reports of what was rejected and why
5. **Flexible**: Adjustable thresholds for different data quality needs

## Basic Usage

### Single Dataset Rejection

```julia
using eegfun
using JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Reject epochs with z-score > 2.0
# Non-mutating version (preserves original)
clean_epochs, rejection_info = reject_epochs_automatic(epochs, 2.0)

# Check results
println("Original epochs: $(rejection_info.n_original)")
println("Remaining epochs: $(rejection_info.n_remaining)")
println("Rejected epochs: $(rejection_info.rejected_epochs)")

# Save cleaned data
save("participant_1_epochs_cleaned.jld2", "epochs", clean_epochs)
```

### In-Place Rejection

```julia
# Mutating version (modifies original data)
rejection_info = reject_epochs_automatic!(epochs, 2.0)

# epochs now contains only clean epochs
```

## Z-Score Criterion

The z-criterion determines how aggressive the rejection is:

| Z-Criterion | Interpretation | Use Case |
|-------------|----------------|----------|
| 1.5 | Very aggressive | Very noisy data, need maximum cleaning |
| 2.0 | Aggressive | Standard rejection, good balance |
| 2.5 | Moderate | High-quality data, remove only obvious artifacts |
| 3.0 | Conservative | Excellent data quality, minimal rejection |
| 3.5+ | Very conservative | Keep maximum number of trials |

**Rule of thumb**: Start with 2.0 and adjust based on results.

```julia
# More aggressive rejection
clean_epochs, info = reject_epochs_automatic(epochs, 1.5)

# More conservative rejection  
clean_epochs, info = reject_epochs_automatic(epochs, 3.0)
```

## How It Works

### Step 1: Calculate Metrics

For each epoch, six metrics are calculated across all channels:

1. **Variance**: Measures signal variability
2. **Maximum**: Highest voltage value
3. **Minimum**: Lowest voltage value (as absolute minimum)
4. **Absolute Maximum**: Largest absolute voltage
5. **Range**: Difference between max and min
6. **Kurtosis**: Measures "spikiness" of distribution

For each metric, the **maximum across channels** is taken to detect global artifacts.

### Step 2: Z-Score Calculation

```julia
# For each metric:
# 1. Calculate metric for all epochs
# 2. Take maximum across channels for each epoch
# 3. Calculate z-score across epochs
# 4. Reject epochs where z-score > criterion
```

### Step 3: Combine Results

An epoch is rejected if it exceeds the z-criterion for **ANY** of the six metrics (OR logic).

```julia
# Example:
# Epoch 5: variance z-score = 1.5, max z-score = 3.2, ...
# Result: REJECTED (because max exceeds 2.0)

# Epoch 10: all z-scores < 2.0
# Result: KEPT
```

## Understanding Rejection Info

The `EpochRejectionInfo` structure tells you exactly what happened:

```julia
clean_epochs, info = reject_epochs_automatic(epochs, 2.0)

# Display summary
println(info)

# Output:
# EpochRejectionInfo:
#   Z-criterion: 2.0
#   Original epochs: 50
#   Remaining epochs: 43
#   Rejected epochs: 7
#
#   Rejection breakdown:
#     Variance:  2 epochs
#     Maximum:   3 epochs
#     Minimum:   1 epochs
#     Absolute:  3 epochs
#     Range:     2 epochs
#     Kurtosis:  4 epochs

# Access specific information
println("Rejected by variance: ", info.rejected_by_variance)
println("Rejected by kurtosis: ", info.rejected_by_kurtosis)
println("All rejected: ", info.rejected_epochs)
```

## Saving Rejection Reports

Generate detailed text reports for documentation:

```julia
clean_epochs, info = reject_epochs_automatic(epochs, 2.0)

# Save detailed report
save_rejection_report(info, "participant_1_rejection_report.txt")
```

The report includes:
- Summary statistics
- Breakdown by rejection criterion
- List of all rejected epoch numbers
- Percentage of epochs rejected

## Channel Selection

By default, all channels are used for artifact detection. You can specify which channels to analyze:

```julia
# Use only frontal channels for detection
clean_epochs, info = reject_epochs_automatic(epochs, 2.0,
    channel_selection = channels(x -> startswith.(string.(x), "F")))

# Use specific channels
clean_epochs, info = reject_epochs_automatic(epochs, 2.0,
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4, :Fz]))

# Exclude EOG channels from analysis
eeg_channels = setdiff(channel_labels(epochs), [:HEOG, :VEOG])
clean_epochs, info = reject_epochs_automatic(epochs, 2.0,
    channel_selection = channels(eeg_channels))
```

## Complete Example Workflows

### Single Participant Processing

```julia
using eegfun
using JLD2

# Load epoched data
epochs = load("participant_05_epochs.jld2", "epochs")

println("Starting with $(length(epochs.data)) epochs")

# Apply automatic rejection
clean_epochs, rejection_info = reject_epochs_automatic(epochs, 2.0)

println("After rejection: $(length(clean_epochs.data)) epochs")
println("Rejection rate: $(round(100 * length(rejection_info.rejected_epochs) / rejection_info.n_original, digits=1))%")

# Save cleaned epochs
save("participant_05_epochs_cleaned.jld2", "epochs", clean_epochs)

# Save rejection report
save_rejection_report(rejection_info, "participant_05_rejection_report.txt")

# Continue with ERP analysis
erps = average_epochs(clean_epochs)
save("participant_05_erps.jld2", "erps", erps)
```

### Multiple Participants (Manual Loop)

```julia
using eegfun
using JLD2

participants = 1:20
all_rejection_info = Dict{Int, EpochRejectionInfo}()

for participant in participants
    filename = "participant_$(participant)_epochs.jld2"
    
    if !isfile(filename)
        @warn "Skipping participant $participant (file not found)"
        continue
    end
    
    # Load and clean
    epochs = load(filename, "epochs")
    clean_epochs, info = reject_epochs_automatic(epochs, 2.0)
    
    # Save cleaned data
    save("participant_$(participant)_epochs_cleaned.jld2", "epochs", clean_epochs)
    save_rejection_report(info, "participant_$(participant)_rejection_report.txt")
    
    # Store info for summary
    all_rejection_info[participant] = info
    
    println("✓ Participant $participant: $(info.n_remaining)/$(info.n_original) epochs kept")
end

# Summary across all participants
total_original = sum(info.n_original for info in values(all_rejection_info))
total_remaining = sum(info.n_remaining for info in values(all_rejection_info))
println("\nOverall: $(total_remaining)/$(total_original) epochs kept ($(round(100 * total_remaining / total_original, digits=1))%)")
```

### Batch Processing (Easiest Method)

```julia
using eegfun

# Process all participants at once
reject_epochs_automatic("epochs", 2.0)

# This automatically:
# - Finds all files matching "epochs" pattern
# - Applies rejection with z=2.0
# - Saves to new directory "rejected_z2p0_epochs/"
# - Creates rejection reports for each participant
# - Creates a log file
```

### Batch Processing with Options

```julia
# Specific participants
reject_epochs_automatic("epochs", 2.5,
                       participants = 1:20)

# Custom directories
reject_epochs_automatic("epochs", 2.0,
                       input_dir = "/data/study1/epochs",
                       output_dir = "/data/study1/clean_epochs")

# Only use specific channels for detection
reject_epochs_automatic("epochs", 2.0,
                       channel_selection = channels(x -> startswith.(string.(x), "F")))

# Full example
reject_epochs_automatic("epochs", 2.5,
                       input_dir = "/data/study1",
                       participants = 1:30,
                       channel_selection = channels(),
                       output_dir = "/data/study1/epochs_cleaned")
```

## Interpreting Results

### Good Rejection Pattern

```
Original epochs: 100
Remaining epochs: 92
Rejected epochs: 8 (8%)

Rejection breakdown:
  Variance:  3 epochs
  Maximum:   2 epochs
  Minimum:   1 epochs
  Absolute:  2 epochs
  Range:     2 epochs
  Kurtosis:  3 epochs
```

**Interpretation**: Reasonable rejection rate with artifacts caught by multiple criteria. This is typical for good data quality.

### High Rejection Rate

```
Original epochs: 100
Remaining epochs: 65
Rejected epochs: 35 (35%)

Rejection breakdown:
  Variance:  15 epochs
  Maximum:   20 epochs
  ...
```

**Possible causes**:
- Very noisy data
- Z-criterion too aggressive (try 2.5 or 3.0)
- Recording issues during session
- Movement artifacts

**Actions**:
1. Check data quality visually
2. Consider more conservative z-criterion
3. Check if specific channels are problematic
4. Review recording notes for issues

### Very Low Rejection Rate

```
Original epochs: 100
Remaining epochs: 99
Rejected epochs: 1 (1%)
```

**Possible interpretations**:
- Excellent data quality
- Z-criterion too conservative (try 2.0 or 1.5)
- Artifacts not captured by these metrics

**Actions**:
1. Visually inspect some epochs
2. Consider additional manual review
3. May be fine if data quality is genuinely excellent

## Advanced Usage

### Comparing Different Criteria

```julia
# Test multiple z-criteria to find optimal setting
z_values = [1.5, 2.0, 2.5, 3.0]
results = Dict()

for z in z_values
    clean, info = reject_epochs_automatic(epochs, z)
    results[z] = info
    println("z=$z: $(info.n_remaining) epochs remaining ($(round(100 * info.n_remaining / info.n_original, digits=1))%)")
end

# Choose criterion that keeps 80-90% of epochs
```

### Identifying Problematic Channels

```julia
# Test rejection channel by channel to identify problematic electrodes
channel_results = Dict()

for ch in channel_labels(epochs)
    clean, info = reject_epochs_automatic(epochs, 2.0,
        channel_selection = channels([ch]))
    channel_results[ch] = length(info.rejected_epochs)
end

# Sort by number of rejections
sorted_channels = sort(collect(channel_results), by = x -> x[2], rev = true)
println("Most problematic channels:")
for (ch, n_rejected) in sorted_channels[1:5]
    println("  $ch: $n_rejected epochs")
end
```

### Combining with Other Preprocessing

```julia
using eegfun, JLD2

# Complete preprocessing pipeline
epochs = load("participant_1_epochs_raw.jld2", "epochs")

# 1. Automatic rejection
epochs_cleaned, rejection_info = reject_epochs_automatic(epochs, 2.0)
println("After automatic rejection: $(length(epochs_cleaned.data)) epochs")

# 2. Manual review (if needed)
# ... visual inspection ...

# 3. Average to ERPs
erps = average_epochs(epochs_cleaned)

# 4. Baseline correction
baseline!(erps, IntervalTime(-0.2, 0.0))

# 5. Further analysis
save("participant_1_erps_final.jld2", "erps", erps)
```

## Quality Control

### Visual Verification

After automatic rejection, spot-check some epochs:

```julia
# Check a few kept epochs
kept_indices = setdiff(1:rejection_info.n_original, rejection_info.rejected_epochs)
sample_kept = sample(kept_indices, min(5, length(kept_indices)))

# Check a few rejected epochs (from original data)
sample_rejected = sample(rejection_info.rejected_epochs, min(5, length(rejection_info.rejected_epochs)))

println("Kept epoch examples: ", sample_kept)
println("Rejected epoch examples: ", sample_rejected)

# Visually inspect these in your original data
```

### Track Rejection Across Participants

```julia
# Create summary table
using DataFrames, CSV

summary = DataFrame(
    participant = Int[],
    n_original = Int[],
    n_remaining = Int[],
    percent_kept = Float64[],
    n_rejected_variance = Int[],
    n_rejected_max = Int[],
    n_rejected_kurtosis = Int[]
)

for participant in 1:20
    info = load("participant_$(participant)_rejection_info.jld2", "info")
    
    push!(summary, (
        participant = participant,
        n_original = info.n_original,
        n_remaining = info.n_remaining,
        percent_kept = 100 * info.n_remaining / info.n_original,
        n_rejected_variance = length(info.rejected_by_variance),
        n_rejected_max = length(info.rejected_by_max),
        n_rejected_kurtosis = length(info.rejected_by_kurtosis)
    ))
end

# Save summary
CSV.write("rejection_summary.csv", summary)

# Identify outliers
mean_kept = mean(summary.percent_kept)
std_kept = std(summary.percent_kept)
outliers = summary[abs.(summary.percent_kept .- mean_kept) .> 2 * std_kept, :]
println("Participants with unusual rejection rates:")
println(outliers)
```

## Common Issues and Solutions

### Issue: Too many epochs rejected

**Solution 1**: Use more conservative z-criterion
```julia
# Instead of z=2.0, try z=2.5 or 3.0
clean_epochs, info = reject_epochs_automatic(epochs, 2.5)
```

**Solution 2**: Check for problematic channels
```julia
# Remove or interpolate bad channels before rejection
```

**Solution 3**: Review data quality
```julia
# May indicate fundamental issues with recording quality
```

### Issue: Too few epochs rejected

**Solution 1**: Use more aggressive z-criterion
```julia
clean_epochs, info = reject_epochs_automatic(epochs, 1.5)
```

**Solution 2**: Check if artifacts are being caught
```julia
# Visually inspect some epochs to verify they're clean
```

### Issue: Specific metric rejecting many epochs

**Example**: Kurtosis rejecting 50% of epochs

**Possible causes**:
- Muscle artifacts (high kurtosis)
- Electrode pops
- Specific channel issues

**Solution**: Check which channels are contributing
```julia
# Test individual channels to identify culprits
```

## References

This implementation is based on the automatic rejection approach used in FieldTrip and EEGLAB, using z-scored statistical metrics across epochs.

## See Also

- Manual epoch rejection: `remove_bad_epochs()` function
- Artifact detection: `artifact_detection.jl`
- ICA for artifact removal: `ica.jl`
- Autoreject algorithm: `autoreject.jl` (if implemented)

