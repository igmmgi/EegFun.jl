# Automated Artifact Handling

Artifacts are an unavoidable aspect of EEG research. `EegFun.jl` provides a structured workflow for identifying artifacts, repairing what can be salvaged, and rejecting what cannot.

## Continuous Data: Sample-Level Detection

Before epoching, you can flag artifact-contaminated samples in continuous data. This is useful for marking segments to exclude from downstream analysis (e.g., ICA training) or for data quality reporting.

### Extreme Values

`is_extreme_value!` flags samples where the absolute voltage exceeds a threshold. This catches amplifier saturation, large movement artifacts, and disconnected electrodes.

```julia
# Flag samples exceeding ±100 μV (adds column :is_extreme_value_100)
EegFun.is_extreme_value!(dat, 100)

# Check how many samples were flagged
EegFun.n_values(dat, :is_extreme_value_100)

# Only check specific channels
EegFun.is_extreme_value!(dat, 100, channel_selection = EegFun.channels([:Fp1, :Fp2]))
```

### Step Values

`is_step_value!` flags samples where the voltage jump between consecutive samples exceeds a threshold. This detects sudden discontinuities (cable pulls, electrode pops) that extreme value detection might miss if the signal returns quickly.

```julia
# Flag jumps > 50 μV between consecutive samples (adds column :is_step_value_50.0)
EegFun.is_step_value!(dat, 50.0)

# Check how many step artifacts were detected
EegFun.n_values(dat, :is_step_value_50.0)
```

> [!TIP]
> Both functions also have non-mutating versions (`is_extreme_value`, `is_step_value`) that return a boolean mask without modifying the data.

Once flagged, these artifact columns can be visualised in `EegFun.plot_databrowser` — select the relevant column from the **Overlay** dropdown menu to highlight flagged samples on the EEG traces.

## Epoch-Level: Automatic Detection

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
