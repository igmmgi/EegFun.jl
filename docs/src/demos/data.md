# Data Loading and Access Demo

Demonstrates basic EEG data loading, structure creation, and data access functions.

## Overview

Shows how to load raw EEG data from BDF files, create EegFun data structures with electrode layouts, and use data access and subsetting functions.

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

# data/channel access functions
EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.channel_data(dat)
EegFun.extra_data(dat)

# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

EegFun.all_data(epochs)
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.meta_data(epochs)
EegFun.meta_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.channel_data(epochs)
EegFun.channel_data(epochs, epoch_selection = EegFun.epochs(1:2))

# test subsetting
dat_subset = EegFun.subset(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
dat_subset = EegFun.subset(dat, sample_selection = x -> x.sample .<= 10_000) # first 10000 samples
dat_subset = EegFun.subset(dat, sample_selection = x -> x.time .<= 10) # first 10 seconds
```

## See Also

- [Data Loading Reference](../reference/data-loading.md)
- [Data Structures](../explanations/data-structures.md)
