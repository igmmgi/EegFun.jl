# Channel Summary Demo

Demonstrates generating summary statistics for all channels.

## Overview

Creates channel-wise summary statistics with optional artifact masking to exclude bad samples from the analysis.

## Source Code

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

## See Also

- [Preprocessing Reference](../reference/preprocessing.md)
