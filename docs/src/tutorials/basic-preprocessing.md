# Basic EEG Preprocessing

This tutorial covers a complete preprocessing pipeline for EEG data, from raw data to clean, analysis-ready signals.

## Overview

A typical preprocessing workflow includes:

1. **High-pass filtering** - Remove slow drifts
2. **Rereferencing** - Choose reference scheme
3. **Artifact detection** - Mark bad channels/timepoints
4. **Low-pass filtering** - Remove high-frequency noise
5. **Visual inspection** - Verify data quality

## Step 1: High-Pass Filtering

Remove slow drifts and DC offsets:

```julia
using EegFun

# Apply 1 Hz high-pass filter
EegFun.highpass_filter!(dat, 1.0)
```

> **TODO**: Explain why 1 Hz is common, discuss alternatives

## Step 2: Rereferencing

Choose an appropriate reference scheme:

```julia
# Average reference (common choice)
EegFun.rereference!(dat, :avg)

# Alternative: specific electrode(s)
# EegFun.rereference!(dat, [:M1, :M2])  # Mastoid reference
```

> **TODO**: Add section explaining reference choices

## Step 3: Artifact Detection

Mark extreme values and bad channels:

```julia
# Mark samples exceeding ±100 µV
EegFun.is_extreme_value!(dat, 100)

# TODO: Add bad channel detection example
```

## Step 4: Low-Pass Filtering

Remove high-frequency noise:

```julia
# Apply 40 Hz low-pass filter
EegFun.lowpass_filter!(dat, 40.0)
```

## Step 5: Visual Inspection

> **TODO**: Add plotting examples for quality control
> - Plot continuous data with artifacts marked
> - Show before/after filtering comparison
> - Channel-by-channel visualization

## Complete Example

```julia
using EegFun

# Load data
dat = EegFun.read_raw_data("data.bdf")
layout = EegFun.read_layout("biosemi64.csv")
dat = EegFun.create_eeg_dataframe(dat, layout)

# Preprocessing pipeline
EegFun.highpass_filter!(dat, 1.0)
EegFun.rereference!(dat, :avg)
EegFun.is_extreme_value!(dat, 100)
EegFun.lowpass_filter!(dat, 40.0)

# Save preprocessed data
# TODO: Add save example
```

## Next Steps

- [Extract epochs](erp-analysis.md) from preprocessed data
- [Apply ICA](ica-workflow.md) for advanced artifact removal
- Learn about [filtering concepts](../explanations/filtering.md)
