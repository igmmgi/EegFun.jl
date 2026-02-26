This demo shows how to plot ERP waveforms with measurement overlays using `EegFun.plot_erp_measurements`.

### What Does It Do?

`plot_erp_measurements` combines ERP waveform visualization with quantified measurements in a single plot. It highlights the analysis interval and displays the computed measurement directly on the waveform, making it easy to see exactly what is being measured and where.

### Supported Measurements

- **Mean amplitude** — average voltage within the analysis interval
- **Peak amplitude** — maximum or minimum voltage within the interval
- **Peak latency** — time point of the peak within the interval

### Key Options

| Parameter | Purpose |
| --- | --- |
| `analysis_interval` | Time window for measurement (seconds) |
| `baseline_interval` | Baseline correction window |
| `channel_selection` | Which channels to plot |
| `condition_selection` | Which conditions to include |
| `layout` | `:overlay` (default) or `:grid` |

## Workflow Summary

This demo covers:

### Mean Amplitude Overlay

- Plot ERP with shaded analysis window and computed mean

### Peak Latency

- Show peak markers on specific channels

### Grid Layout

- Display each channel in its own panel with measurement annotations

### File Path Input

- Load data directly from JLD2 file path without pre-loading
