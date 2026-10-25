# Plot ERP Filter GUI

This demo shows how to use the interactive filter GUI to explore the effect of different filter settings on ERP waveforms.

## Why an Interactive Filter Tool?

- **Teaching** — show students how filters affect ERP signals in real time
- **Parameter selection** — experiment with cutoff frequency and filter order before committing
- **Method comparison** — toggle between IIR (Butterworth) and FIR filters, `filtfilt` vs `filt`

## Key Features

| Feature | Description |
| --- | --- |
| Lowpass / Highpass toggles | Enable each filter independently |
| Cutoff sliders | Adjust frequency in real time |
| Order sliders | Change filter steepness |
| Method menu | Switch between Butterworth and FIR |
| Mirror toggle | Enable mirror padding to reduce edge artifacts |
| Multi-condition | Compare filter effects across conditions side-by-side |

## What You'll Learn

1. Launching the GUI for a single ERP
2. Focusing on a specific channel
3. Comparing filter effects across multiple conditions


## Code Examples

::: details Show Code

```julia
# Demo: ERP Filter GUI
# Shows how to launch the interactive filter exploration tool
# to compare the effect of different filter settings on ERP waveforms.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# Launch with default demo data (bundled example ERPs, channel :Cz):
gui = EegFun.plot_erp_filter_gui()

# Alternatively, specify a channel or custom data / file:
# gui = EegFun.plot_erp_filter_gui(channel = :Pz)
# dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
# gui = EegFun.plot_erp_filter_gui(dat, channel = :Cz)

# Keep window open when run from terminal
if !isinteractive() && !isnothing(gui.screen)
    wait(gui.screen)
end
```

:::

## See Also

- [API Reference](../../reference/index.md)
