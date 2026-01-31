# Trigger Marker Visualization

Visualizes event markers and triggers in continuous data.

## Overview

Demonstrates Visualizes event markers and triggers in continuous data.

## Source Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# count from raw file
count = EegFun.trigger_count(dat)

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

count = EegFun.trigger_count(dat)

EegFun.plot_trigger_overview(dat)
EegFun.plot_trigger_overview(dat; ignore_triggers = [3])

EegFun.plot_trigger_timing(dat)
EegFun.plot_trigger_timing(dat; ignore_triggers = [3])
```

## See Also

- [API Reference](../reference/index.md)
