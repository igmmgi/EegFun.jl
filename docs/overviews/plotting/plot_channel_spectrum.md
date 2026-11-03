This demo shows how to plot power spectra from `SpectrumData` objects.

### `plot_channel_spectrum` Multiple Dispatch

| Function | Input | Best For |
| --- | --- | --- |
| `plot_channel_spectrum(dat)` | Raw/preprocessed data | Quick interactive inspection with Welch's method |
| `plot_channel_spectrum(spectrum)` | `SpectrumData` | Plotting pre-computed spectra with full control |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `channel_selection` | `EegFun.channels()` | Which channels to plot |
| `xscale` | `:linear` | `:linear` or `:log10` |
| `yscale` | `:linear` | `:linear` or `:log10` |
| `unit` | `:linear` | `:linear` (μV²/Hz) or `:dB` |
| `max_freq` | `nothing` | Maximum frequency to display |
| `legend` | `true` | Show channel name legend |

### What You'll Learn

1. Plotting spectra for all or selected channels
2. Switching between linear and log scales
3. Using decibel units
4. Limiting the frequency range
