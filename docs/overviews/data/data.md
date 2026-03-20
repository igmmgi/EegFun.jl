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
