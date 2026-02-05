# Data

## Overview

## Overview

This demo demonstrates loading EEG data and understanding the basic data structures.

### Supported File Formats

- **BioSemi (.bdf)**: High-resolution EEG recordings
- **BrainVision (.vhdr)**: BrainProducts format
- **EEGLAB (.set)**: Matlab-based EEG format
- **FieldTrip (.mat)**: Matlab-based format

### Core Data Structures

**ContinuousData:**
- Raw EEG time series
- Contains electrode data, triggers, sampling rate, metadata

**EpochData:**
- Segmented trials around events
- Organized by experimental conditions

**ErpData:**
- Averaged event-related potentials
- One waveform per condition


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

# data/channel access functions
EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.channel_data(dat)
EegFun.extra_data(dat)

# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

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

:::

## See Also

- [API Reference](../reference/index.md)
