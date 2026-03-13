# Plot Frequency Spectrum

This demo shows how to plot power spectra from `SpectrumData` objects.

### When to Use `plot_frequency_spectrum` vs `plot_channel_spectrum`?

| Function | Input | Best For |
| --- | --- | --- |
| `plot_channel_spectrum(dat)` | Raw/preprocessed data | Quick interactive inspection with Welch's method |
| `plot_frequency_spectrum(spectrum)` | `SpectrumData` | Plotting pre-computed spectra with full control |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `channel_selection` | `channels()` | Which channels to plot |
| `x_scale` | `:linear` | `:linear` or `:log10` |
| `y_scale` | `:linear` | `:linear` or `:log10` |
| `unit` | `:linear` | `:linear` (μV²/Hz) or `:dB` |
| `max_freq` | `nothing` | Maximum frequency to display |
| `show_legend` | `true` | Show channel name legend |

### What You'll Learn

1. Plotting spectra for all or selected channels
2. Switching between linear and log scales
3. Using decibel units
4. Limiting the frequency range


## Code Examples

::: details Show Code

```julia
# Demo: Plotting Power Spectrum from SpectrumData
# Shows how to visualise frequency spectra computed by freq_spectrum,
# with options for channel selection, axis scaling, and dB units.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

#######################################################################
# LOAD DATA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.5)

#######################################################################
# COMPUTE SPECTRUM
#######################################################################

spectrum = EegFun.freq_spectrum(dat)

#######################################################################
# BASIC SPECTRUM PLOT — ALL CHANNELS
#######################################################################

EegFun.plot_frequency_spectrum(spectrum)

#######################################################################
# SINGLE CHANNEL
#######################################################################

EegFun.plot_frequency_spectrum(spectrum, channel_selection = EegFun.channels(:Cz))

#######################################################################
# LOG SCALE AXES
#######################################################################

EegFun.plot_frequency_spectrum(spectrum, channel_selection = EegFun.channels(:Cz), x_scale = :log10, y_scale = :log10)

#######################################################################
# DECIBEL UNITS
#######################################################################

EegFun.plot_frequency_spectrum(spectrum, channel_selection = EegFun.channels([:Cz, :Pz]), unit = :dB, max_freq = 100.0)

#######################################################################
# CUSTOM STYLING
#######################################################################

EegFun.plot_frequency_spectrum(
    spectrum,
    channel_selection = EegFun.channels([:Cz, :Oz]),
    linewidth = 3,
    line_alpha = 0.6,
    title = "Power Spectrum",
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
