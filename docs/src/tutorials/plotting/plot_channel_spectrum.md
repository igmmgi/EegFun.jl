# Plot Channel Spectrum

This demo shows how to plot power spectra from `SpectrumData` objects.

## `plot_channel_spectrum` Multiple Dispatch

| Function | Input | Best For |
| --- | --- | --- |
| `plot_channel_spectrum(dat)` | Raw/preprocessed data | Quick interactive inspection with Welch's method |
| `plot_channel_spectrum(spectrum)` | `SpectrumData` | Plotting pre-computed spectra with full control |

## Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `channel_selection` | `channels()` | Which channels to plot |
| `x_scale` | `:linear` | `:linear` or `:log10` |
| `y_scale` | `:linear` | `:linear` or `:log10` |
| `unit` | `:linear` | `:linear` (μV²/Hz) or `:dB` |
| `max_freq` | `nothing` | Maximum frequency to display |
| `show_legend` | `true` | Show channel name legend |

## What You'll Learn

1. Plotting spectra for all or selected channels
2. Switching between linear and log scales
3. Using decibel units
4. Limiting the frequency range


## Code Examples

::: details Show Code

```julia
# Demo: Channel Power Spectrum Plotting
# Shows how to visualise power spectra from both raw continuous data 
# and pre-computed SpectrumData, featuring interactive controls for 
# axis scaling, unit conversions, and frequency band overlays.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
# LOAD AND PREPARE DATA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.5)

#######################################################################
# WORKFLOW 1: COMPUTING ON THE FLY FROM RAW DATA
#######################################################################
# plot_channel_spectrum can compute the spectrum dynamically from raw data

# Plot all channels
EegFun.plot_channel_spectrum(dat)

# Plot a single channel with a specific title
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1]), title = "Fp1 Power Spectrum")

# Plot multiple frontal channels
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :F4]), title = "Frontal Channels Power Spectrum")


#######################################################################
# WORKFLOW 2: PLOTTING PRE-COMPUTED SPECTRUM DATA
#######################################################################
# If you have already run freq_spectrum, you can pass the result directly

# Compute spectrum first
spectrum = EegFun.freq_spectrum(dat)

# Plot all channels from pre-computed data
EegFun.plot_channel_spectrum(spectrum)

# Single channel from pre-computed data
EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels(:Cz))

#######################################################################
# INTERACTIVE/INITIAL PARAMETER OVERRIDES
#######################################################################
# Note: You can toggle these visually in the GUI, but you can also 
# set them as starting defaults via keyword arguments.

# Start with log-scale axes
EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels(:Cz), x_scale = :log10, y_scale = :log10)

# Start with Decibel (dB) units instead of linear power
EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels([:Cz, :Pz]), unit = :dB, max_freq = 100.0)

# Custom styling for lines
EegFun.plot_channel_spectrum(
    spectrum,
    channel_selection = EegFun.channels([:Cz, :Oz]),
    linewidth = 3,
    line_alpha = 0.6,
    title = "Custom Styled Power Spectrum",
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
