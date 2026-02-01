# Resample

Resample data to different sampling rates.

## Overview

Demonstrates Resample data to different sampling rates.

## Source Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

EegFun.sample_rate(dat)   # current sample rate
EegFun.trigger_count(dat) # current triggers in file

dat_new = EegFun.resample(dat, 2) # downsample by a factor of 2
EegFun.sample_rate(dat_new)       # should = original ÷ 2
EegFun.trigger_count(dat_new)     # triggers should be preserved

dat_new = EegFun.resample(dat, 4) # downsample by a factor of 4
EegFun.sample_rate(dat_new)       # should = original ÷ 4
EegFun.trigger_count(dat_new)     # triggers should be preserved
```

## See Also

- [API Reference](../reference/index.md)
