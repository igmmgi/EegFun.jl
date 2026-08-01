# Plot Time-Frequency

This demo shows how to plot time-frequency decomposition results as heatmaps.

## When to Use

Use `plot_tf` to visualise the output of `tf_morlet`, `tf_multitaper`, or `tf_stft` — the three time-frequency decomposition methods in EegFun.

## Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `channel` | First channel | Which channel to plot |
| `baseline_interval` | `nothing` | Pre-stimulus interval for on-the-fly baseline correction |
| `baseline_method` | `:db` | `:db`, `:percent`, `:relchange`, `:absolute`, `:relative`, or `:zscore` |
| `colormap` | `:viridis` | Colour map for the heatmap |
| `colorrange` | auto | Explicit `(min, max)` colour range |
| `ylogscale` | `false` | Log-scale the frequency axis |
| `colorbar` | `true` | Show colour bar |

## What You'll Learn

1. Plotting time-frequency data for a specific channel
2. Applying baseline correction on the fly (dB, percent, relative change)
3. Using log-scaled frequency axes
4. Customising colour maps and colour range for publication figures


## Code Examples

::: details Show Code

```julia
# Demo: Plotting Time-Frequency Data
# Shows how to plot time-frequency decomposition results as heatmaps,
# with options for baseline correction, colour maps, and log-scaled axes.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

#######################################################################
# LOAD DATA AND RUN TF DECOMPOSITION
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.5)
EegFun.lowpass_filter!(dat, 100.0)

epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))  # -200 to 1000 ms

EegFun.baseline!(epochs, (-0.2, 0.0))

# Morlet wavelet decomposition
tf_data = EegFun.tf_morlet(epochs, frequencies = 2:2:60)

#######################################################################
# BASIC TIME-FREQUENCY PLOT
#######################################################################

# plot first available channel
EegFun.plot_tf(tf_data)

#######################################################################
# SPECIFIC CHANNEL
#######################################################################

EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz))

#######################################################################
# WITH BASELINE CORRECTION
#######################################################################

# apply dB baseline 
EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz), baseline_interval = (-0.3, 0.0), baseline_method = :db)

# percentage change baseline
EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz), baseline_interval = (-0.3, 0.0), baseline_method = :percent)

#######################################################################
# LOG-SCALED FREQUENCY AXIS
#######################################################################

EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz), baseline_interval = (-0.3, 0.0), ylogscale = true)

#######################################################################
# CUSTOM COLOUR MAP AND RANGE
#######################################################################

EegFun.plot_tf(
    tf_data,
    channel_selection = EegFun.channels(:Cz),
    baseline_interval = (-0.3, 0.0),
    colormap = :RdBu,
    colorrange = (-15.0, 15.0),
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
