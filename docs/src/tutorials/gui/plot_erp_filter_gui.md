# Plot ERP Filter GUI

This demo shows how to use the interactive filter GUI to explore the effect of different filter settings on ERP waveforms.

## Why an Interactive Filter Tool?

- **Teaching** — show students how filters affect ERP signals in real time
- **Parameter selection** — experiment with cutoff frequency and filter order before committing
- **Method comparison** — toggle between IIR (Butterworth) and FIR filters, and compare two-pass `filtfilt`, MNE/EEGLAB-style one-pass zero-phase FIR (`filt_zerophase`), and causal one-pass `filt`

## Key Features

| Feature | Description |
| --- | --- |
| Lowpass / Highpass toggles | Enable each filter independently |
| Cutoff sliders | Adjust frequency in real time |
| Order / FIR tap readout | Adjust IIR filter steepness or inspect live FIR tap count and group delay |
| Method & Function menus | Compare Butterworth IIR vs FIR Hamming across `filtfilt`, `filt_zerophase`, and `filt` |
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
# Shows how to launch the interactive filter exploration GUI to compare
# IIR (Butterworth) and FIR filter settings on ERP and continuous waveforms.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# gui = EegFun.plot_erp_filter_gui("/path/to/participant_1_erps.jld2")

using EegFun
using GLMakie

#######################################################################
# 1. QUICK DEMO (BUNDLED ERP DATA)
#######################################################################

# Controls to test:
# - Channel dropdown, Filter Method (Butterworth IIR vs FIR Hamming)
# - Filter Function: filtfilt (2-pass zero-phase), filt_zerophase (1-pass FIR, MNE/EEGLAB style), filt (1-pass causal)
# - Mirror padding toggle, Low-pass & High-pass enable toggles, Cutoff & Order sliders (plus live FIR tap/delay readout)
gui = EegFun.plot_erp_filter_gui()

#######################################################################
# 2. MULTI-CONDITION ERP & SPECIFIC CHANNEL
#######################################################################

erps = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
# gui = EegFun.plot_erp_filter_gui(erps, channel = :Pz)

#######################################################################
# 3. CONTINUOUS DATA FILTER EXPLORATION
#######################################################################

# raw = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
# dat = EegFun.create_eegfun_data(raw)
# EegFun.rereference!(dat, :avg)
# gui = EegFun.plot_erp_filter_gui(dat, channel = :Cz)

# Keep window open when run as a script from terminal
if !isinteractive() && !isnothing(gui.screen)
    wait(gui.screen)
end
```

:::

## See Also

- [API Reference](../../reference/index.md)
