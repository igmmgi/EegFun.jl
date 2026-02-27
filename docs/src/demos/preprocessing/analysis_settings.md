# Analysis Settings

This demo shows how to save and replay preprocessing settings using `AnalysisSettings`.

### What are Analysis Settings?

`AnalysisSettings` stores a preprocessing recipe that can be applied to data in one step. This is useful for:

- Reproducing the same preprocessing across multiple datasets
- Storing settings from interactive exploration (e.g., databrowser GUI)
- Standardising processing across participants in a study

### Tracked Settings

An `AnalysisSettings` object can contain:

- **`hp_filter`**: High-pass filter cutoff frequency
- **`lp_filter`**: Low-pass filter cutoff frequency
- **`reference`**: Re-referencing scheme (e.g., `:avg`, `:Cz`)
- **`repaired_channels`**: Channels to interpolate
- **`repair_method`**: Interpolation method
- **`removed_ica_components`**: ICA components to remove
- **`selected_regions`**: ROI definitions

### Applying Settings

- **`apply_analysis_settings!(dat, settings)`**: Apply in-place (mutating)
- **`apply_analysis_settings(dat, settings)`**: Return a processed copy (non-mutating)
- **`apply_analysis_settings!(dat, ica, settings)`**: Apply with ICA component removal

### Analysis Info Tracking

Every EegFun data object automatically tracks its preprocessing history via `analysis_info`. This records what filters, references, and other operations have been applied.

## Workflow Summary

This demo covers:

- Inspecting analysis info on data objects
- Creating AnalysisSettings with filter and reference parameters
- Applying settings to data (both mutating and non-mutating)


## Code Examples

::: details Show Code

```julia
# Demo: Analysis Settings
# Shows how to save and replay preprocessing settings
# (filter, rereference, ICA) using AnalysisSettings.

using EegFun

#######################################################################
# SETUP
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)


#######################################################################
# INSPECTING ANALYSIS INFO
#######################################################################

# Every EegFun data object tracks its preprocessing state
dat.analysis_info  # shows reference, filter info

# Apply some preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.rereference!(dat, :avg)

# Check that analysis info is updated
dat.analysis_info  # now shows hp_filter = 0.1, reference = :avg


#######################################################################
# USING ANALYSIS SETTINGS
#######################################################################

# AnalysisSettings stores a recipe for preprocessing
settings = EegFun.AnalysisSettings(1.0, 30.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[])

# Apply settings to fresh data (non-mutating — returns a processed copy)
dat_fresh = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
dat_fresh = EegFun.create_eegfun_data(dat_fresh, layout)

dat_processed = EegFun.apply_analysis_settings(dat_fresh, settings)
dat_processed.analysis_info

# Or apply in-place (mutating)
# EegFun.apply_analysis_settings!(dat_fresh, settings)
```

:::

## See Also

- [API Reference](../../reference/index.md)
