# Plot ERP Measurements

This demo shows how to plot ERP waveforms with measurement overlays using `EegFun.plot_erp_measurements`.

## What Does It Do?

`plot_erp_measurements` combines ERP waveform visualization with quantified measurements in a single plot. It highlights the analysis interval and displays the computed measurement directly on the waveform, making it easy to see exactly what is being measured and where.

## Supported Measurements

- **Mean amplitude** — average voltage within the analysis interval
- **Peak amplitude** — maximum or minimum voltage within the interval
- **Peak latency** — time point of the peak within the interval

## Key Options

| Parameter | Purpose |
| --- | --- |
| `analysis_interval` | Time window for measurement (seconds) |
| `baseline_interval` | Baseline correction window |
| `channel_selection` | Which channels to plot |
| `condition_selection` | Which conditions to include |
| `layout` | `:overlay` (default) or `:grid` |

## Workflow Summary

This demo covers:

## Mean Amplitude Overlay

- Plot ERP with shaded analysis window and computed mean

## Peak Latency

- Show peak markers on specific channels

## Grid Layout

- Display each channel in its own panel with measurement annotations

## File Path Input

- Load data directly from JLD2 file path without pre-loading


## Code Examples

::: details Show Code

```julia
# Demo: Plotting ERPs with Measurement Overlays
# Computes and visualizes ERP measurements in a single call.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_final.jld2"))

# Mean amplitude
# EegFun.plot_erp_measurements(dat, "mean_amplitude", analysis_interval = (0.3, 0.5), baseline_interval = (-0.2, 0.0)) # not really useful as too crowded!

# Peak latency with specific channels
EegFun.plot_erp_measurements(
    dat,
    "max_peak_latency",
    analysis_interval = (0.0, 1.0),
    baseline_interval = (-0.2, 0.0),
    channel_selection = EegFun.channels([:Cz, :Pz]),
)

# Grid layout
EegFun.plot_erp_measurements(dat, "max_peak_amplitude", analysis_interval = (0.0, 1.0), baseline_interval = (-0.2, 0.0), layout = :grid)

# Load from file path
EegFun.plot_erp_measurements(
    EegFun.example_path("data/julia/erps/example1_erps_final.jld2"),
    "max_peak_amplitude",
    analysis_interval = (0.0, 1.0),
    baseline_interval = (-0.2, 0.0),
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Cz, :Pz]),
    average_channels = true,
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
