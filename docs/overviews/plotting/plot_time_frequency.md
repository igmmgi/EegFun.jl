This demo shows how to plot time-frequency decomposition results as heatmaps.

### When to Use

Use `plot_time_frequency` to visualise the output of `tf_morlet`, `tf_multitaper`, or `tf_stft` — the three time-frequency decomposition methods in EegFun.

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `channel` | First channel | Which channel to plot |
| `baseline_interval` | `nothing` | Pre-stimulus interval for on-the-fly baseline correction |
| `baseline_method` | `:db` | `:db`, `:percent`, or `:relchange` |
| `colormap` | `:viridis` | Colour map for the heatmap |
| `colorrange` | auto | Explicit `(min, max)` colour range |
| `ylogscale` | `false` | Log-scale the frequency axis |
| `colorbar` | `true` | Show colour bar |

### What You'll Learn

1. Plotting time-frequency data for a specific channel
2. Applying baseline correction on the fly (dB, percent, relative change)
3. Using log-scaled frequency axes
4. Customising colour maps and colour range for publication figures
