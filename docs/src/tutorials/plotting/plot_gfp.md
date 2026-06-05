# Plot GFP

This demo shows how to visualise Global Field Power (GFP) and Global Dissimilarity for ERP data.

### What is GFP?

- **GFP** is the spatial standard deviation across all electrodes at each time point
- It summarises the overall strength of the scalp field, independent of polarity
- **Global Dissimilarity** measures how much the scalp topography changes between time points

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `plot_gfp(erp)` | Plot GFP for one condition | Quick inspection |
| `plot_gfp([erp1, erp2])` | Compare GFP across conditions | Condition effects |
| `plot_gfp(gfp_dataframe)` | Plot pre-computed GFP | Custom workflows |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `normalize` | `true` | Percentage (true) or raw μV (false) |
| `show_erp_traces` | `false` | Add overlaid ERP waveforms |
| `show_dissimilarity` | `false` | Add Global Dissimilarity panel |
| `channel_selection` | `channels()` | Subset of channels for GFP |

### What You'll Learn

1. Plotting GFP for single and multiple conditions
2. Showing ERP traces and dissimilarity in a multi-panel layout
3. Using pre-computed GFP results
4. Selecting specific channels for GFP calculation


## Code Examples

::: details Show Code

```julia
# Demo: Plotting Global Field Power
# Shows how to visualise GFP (Global Field Power) and Global Dissimilarity
# across conditions, with optional ERP trace overlay.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

#######################################################################
# LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# extract epochs
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

EegFun.baseline!(epochs, (-0.2, 0.0))
erps = EegFun.average_epochs(epochs)

# plot GFP for a single condition
EegFun.plot_gfp(erps)

```

:::

## See Also

- [API Reference](../../reference/index.md)
