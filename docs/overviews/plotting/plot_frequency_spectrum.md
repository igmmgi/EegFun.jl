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
