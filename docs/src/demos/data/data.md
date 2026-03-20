# Data

This demo demonstrates the core data structures in EegFun.jl and how to access and subset data.

### Core Data Structures

**ContinuousData**:

- Raw EEG time series from continuous recording
- Contains electrode data, triggers, sampling rate, metadata
- Single DataFrame

**EpochData**:

- Segmented trials around experimental events
- Organized by experimental conditions
- Vector of DataFrames (one per trial)

**ErpData**:

- Averaged event-related potentials
- Single DataFrame

### Data Access Functions

**all_data**: Complete data including channels and metadata

**channel_data**: Only EEG channel columns

**meta_data**: Only metadata columns (time, sample, triggers, etc.)

**extra_data**: Derived/computed channels (EOG, artifacts, boolean flags, etc.)

## Workflow Summary

This demo shows basic data operations:

### Load and Structure Data

- Load raw BioSemi data
- Load electrode layout
- Create ContinuousData structure

### Access Data Components

- Extract all data, channel data, and metadata
- Understand different data access patterns

### Create Epochs

- Define epoch conditions around triggers
- Extract segmented trials
- Access epoch data with selection functions

### Subset Data

- Subset by channel selection (specific electrodes)
- Subset by sample selection (time intervals)
- Create smaller datasets for analysis


## Code Examples

::: details Show Code

```julia
# Demo: Data Manipulation
# Shows basic data access, subsetting, and manipulation functions.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# data/channel access functions
EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.channel_data(dat)
EegFun.extra_data(dat) # empty

# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

EegFun.all_data(epochs)
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.meta_data(epochs)
EegFun.meta_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.channel_data(epochs)
EegFun.channel_data(epochs, epoch_selection = EegFun.epochs(1:2))

# We can subset out EegFun datatypes
dat_subset = EegFun.subset(dat, channel_selection = EegFun.channels([:Fp1, :Fp2])) # only Fp1 and Fp2
EegFun.all_data(dat_subset)

dat_subset = EegFun.subset(dat, sample_selection = x -> x.sample .<= 10_000)       # first 10000 samples
EegFun.all_data(dat_subset)

dat_subset = EegFun.subset(dat, sample_selection = x -> x.time .<= 10)             # first 10 seconds
EegFun.all_data(dat_subset)
```

:::

## See Also

- [API Reference](../../reference/index.md)
