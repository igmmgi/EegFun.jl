# Automated Artifact Handling

Low-quality data is an unavoidable part of EEG research. `EegFun.jl` provides a structured workflow for identifying artifacts, repairing what can be salvaged, and rejecting what cannot.

This tutorial follows the standard pipeline for automated artifact handling.

## 1. Automatic Detection

The first step is identifying problematic segments. `EegFun.jl` uses statistical z-scores and absolute voltage thresholds to flag "bad" epochs.

```julia
# Detect bad epochs using default criteria:
# - Z-scores > 3.0 for variance, range, kurtosis, etc.
# - Absolute voltage > 100 μV
artifacts = detect_bad_epochs_automatic(epochs)

# You can customize these thresholds:
artifacts = detect_bad_epochs_automatic(
    epochs, 
    z_criterion = 2.5,     # More aggressive z-score
    abs_criterion = 80.0   # 80 μV threshold
)
```

The returned `EpochRejectionInfo` object contains a breakdown of exactly why each epoch was flagged.

## 2. Channel Repair

If only one or two electrodes are noisy in an otherwise clean trial, you can "repair" them using information from their neighbors.

### Analyze Repairability

Before repairing, use `channel_repairable!` to check if bad electrodes have enough clean neighbors to be salvaged.

```julia
# Identify which artifacts can be fixed via neighbors
channel_repairable!(artifacts, epochs.layout)
```

### Apply Repair

You can use simple neighbor interpolation or more advanced spherical splines.

```julia
# Method 1: Weighted Neighbor Interpolation (Fast)
repair_artifacts!(epochs, artifacts, method = :neighbor_interpolation)

# Method 2: Spherical Spline Interpolation (High Quality)
repair_artifacts!(epochs, artifacts, method = :spherical_spline)
```

> [!NOTE]
> Spherical Spline interpolation requires a 3D layout (Cartesian coordinates).

## 3. Epoch Rejection

After repairing what is possible, you must remove the epochs that are still flagged as "bad" (e.g., those where too many neighbors were also noisy).

```julia
# Remove the remaining bad epochs in-place
reject_epochs!(epochs, artifacts)
```

## 4. Participant-Level Subsetting

Finally, if a participant has lost too much data across a condition, it may be better to exclude them from the study entirely.

```julia
# Move participants with < 70% data retention to an "excluded" folder
subset_bad_data("preprocessed_files", 70.0)
```

This updates the `epoch_summary_subset.jld2` and `file_summary_subset.jld2` files, ensuring your group-level analysis only includes high-quality participants.

## Summary Workflow

A typical artifact handling script looks like this:

```julia
# 1. Detect
artifacts = detect_bad_epochs_automatic(epochs, z_criterion=3, abs_criterion=100)

# 2. Repair (Salvage single bad channels)
channel_repairable!(artifacts, epochs.layout)
repair_artifacts!(epochs, artifacts, method=:neighbor_interpolation)

# 3. Reject (Remove non-repairable trials)
reject_epochs!(epochs, artifacts)
```
