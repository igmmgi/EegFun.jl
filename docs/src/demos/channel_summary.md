# Channel Summary

## Overview

## Overview

This demo shows how to generate summary statistics across all channels.

### Channel Summary Statistics

Summarize data characteristics across channels:
- **Amplitude statistics**: Mean, variance, min, max per channel
- **Temporal statistics**: Time-domain characteristics
- **Quality metrics**: Identify outlier channels

### Use Cases

- **Data quality overview**: Quick assessment of recording quality
- **Channel comparison**: Identify systematic differences across the montage
- **Documentation**: Generate summary statistics for methods sections


## Code Examples

::: details Show Code

```julia
using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file)a

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

summary = EegFun.channel_summary(dat)
EegFun.log_pretty_table(summary; title = "Initial Channel Summary")

EegFun.is_extreme_value!(dat, 100);
summary = EegFun.channel_summary(dat, sample_selection = EegFun.samples_not(:is_extreme_value_100))
```

:::

## See Also

- [API Reference](../reference/index.md)
