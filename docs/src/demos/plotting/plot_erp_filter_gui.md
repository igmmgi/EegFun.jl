# Plot ERP Filter GUI

This demo shows how to use the interactive filter GUI to explore the effect of different filter settings on ERP waveforms.

### Why an Interactive Filter Tool?

- **Teaching** — show students how filters affect ERP signals in real time
- **Parameter selection** — experiment with cutoff frequency and filter order before committing
- **Method comparison** — toggle between IIR (Butterworth) and FIR filters, `filtfilt` vs `filt`

### Key Features

| Feature | Description |
| --- | --- |
| Lowpass / Highpass toggles | Enable each filter independently |
| Cutoff sliders | Adjust frequency in real time |
| Order sliders | Change filter steepness |
| Method menu | Switch between Butterworth and FIR |
| Mirror toggle | Enable mirror padding to reduce edge artifacts |
| Multi-condition | Compare filter effects across conditions side-by-side |

### What You'll Learn

1. Launching the GUI for a single ERP
2. Focusing on a specific channel
3. Comparing filter effects across multiple conditions


## Code Examples

::: details Show Code

```julia
# Demo: ERP Filter GUI
# Shows how to launch the interactive filter exploration tool
# to compare the effect of different filter settings on ERP waveforms.

using EegFun

#######################################################################
# 1. LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# epoch and average
epochs = EegFun.epoch_data(dat, [:trigger1], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))
erp = EegFun.average(epochs)

#######################################################################
# 2. LAUNCH FILTER GUI — SINGLE CONDITION
#######################################################################

# interactive GUI with sliders for cutoff, order, method
EegFun.plot_erp_filter_gui(erp)

#######################################################################
# 3. SPECIFY A CHANNEL
#######################################################################

# focus on a specific channel
EegFun.plot_erp_filter_gui(erp, channel = :Cz)

#######################################################################
# 4. COMPARE MULTIPLE CONDITIONS
#######################################################################

# pass a vector of ERPs to see filter effects side-by-side
epochs2 = EegFun.epoch_data(dat, [:trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs2, (-0.2, 0.0))
erp2 = EegFun.average(epochs2)

EegFun.plot_erp_filter_gui([erp, erp2], channel = :Cz)
```

:::

## See Also

- [API Reference](../../reference/index.md)
