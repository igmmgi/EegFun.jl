# Data

This demo demonstrates understanding the core data structures in EegFun.jl and how to access and subset data.

This demo demonstrates understanding the core data structures in EegFun.jl and how to access and subset data.

### Core Data Structures

**ContinuousData**:

- Raw EEG time series from continuous recording
- Contains electrode data, triggers, sampling rate, metadata
- Used for preprocessing, filtering, artifact detection

**EpochData**:

- Segmented trials around experimental events
- Organized by experimental conditions
- Vector of DataFrames (one per trial)

**ErpData**:

- Averaged event-related potentials
- One waveform per condition
- Used for statistical analysis and plotting

### Data Access Functions

**all_data**: Complete data including channels and metadata

**channel_data**: Only EEG channel columns

**meta_data**: Only metadata columns (time, sample, triggers, etc.)

**extra_data**: Derived/computed channels (EOG, artifacts, boolean flags, etc.)

## Workflow Summary

This demo shows basic data operations:

### 1. Load and Structure Data

- Load raw BioSemi data
- Load electrode layout
- Create ContinuousData structure

### 2. Access Data Components

- Extract all data, channel data, and metadata
- Understand different data access patterns

### 3. Create Epochs

- Define epoch conditions around triggers
- Extract segmented trials
- Access epoch data with selection functions

### 4. Subset Data

- Subset by channel selection (specific electrodes)
- Subset by sample selection (time windows)
- Create smaller datasets for analysis


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

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

- [API Reference](../reference/index.md)
