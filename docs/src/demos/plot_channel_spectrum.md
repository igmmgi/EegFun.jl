# Plot Channel Spectrum

## Overview

## Overview

This demo demonstrates frequency spectrum visualization for individual channels.

### Frequency Analysis

Visualize the power spectral density of EEG channels:
- **Frequency bands**: Delta, theta, alpha, beta, gamma
- **Peak detection**: Identify dominant frequencies
- **Channel comparison**: Compare spectral profiles across electrodes

### Use Cases

- **Quality control**: Detect line noise (50/60 Hz) or other artifacts
- **Band-specific analysis**: Quantify alpha peaks, theta power, etc.
- **Pre/post filtering validation**: Verify filter effects


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Plots
EegFun.plot_channel_spectrum(dat)
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1]), title = "Fp1 Power Spectrum")
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :F4]), title = "Frontal Channels Power Spectrum")
```

:::

## See Also

- [API Reference](../reference/index.md)
