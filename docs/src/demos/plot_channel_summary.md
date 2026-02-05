# Plot Channel Summary

## Overview

## Overview

This demo shows how to visualize channel summary statistics across the montage.

### Summary Visualizations

Plot aggregate statistics for all channels:
- **Topographic maps**: Spatial distribution of variance, kurtosis, etc.
- **Bar plots**: Channel-by-channel comparison
- **Outlier detection**: Highlight problematic electrodes

### Applications

- Quality control visualization
- Identify spatial patterns in noise
- Document data characteristics


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

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# basic channel summary statistics
cs = EegFun.channel_summary(dat)
cs = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
cs = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]), sample_selection = x -> x.sample .< 2000)
cs = EegFun.channel_summary(dat, channel_selection = x -> endswith.(string.(x), "z")) # all midline channels 
cs = EegFun.channel_summary(dat, channel_selection = x -> .!(endswith.(string.(x), "z"))) # all non-midline channels 

# Plotting Channel Summaries
EegFun.plot_channel_summary(cs, :range)
EegFun.plot_channel_summary(cs, :min)
EegFun.plot_channel_summary(cs, :min, bar_color = :red)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar])

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

cs = EegFun.channel_summary(epochs[1])a

EegFun.plot_channel_summary(cs, :range, average_over = :epoch)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar], average_over = :epoch)
```

:::

## See Also

- [API Reference](../reference/index.md)
