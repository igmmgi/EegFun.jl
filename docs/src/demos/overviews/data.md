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
