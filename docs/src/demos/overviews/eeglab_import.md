This demo demonstrates importing EEGLAB `.set` files into EegFun.jl, including both continuous and epoched data with ICA decompositions.

### What is EEGLAB Format?

EEGLAB is a widely-used MATLAB toolbox for EEG analysis. The `.set` format stores:

- EEG data (continuous or epoched)
- Channel locations and metadata
- Event markers and triggers
- ICA decompositions (if computed)
- Preprocessing history
- Lots more ???

### Import Capabilities

**Continuous data**:

- Raw EEG time series
- Triggers and event markers
- Channel locations
- Sampling rate and metadata

**Epoched data**:

- Segmented trials
- Trial-specific event information
- Baseline windows
- Rejection markers

**ICA information**:

- Independent component activations
- Unmixing and mixing matrices
- Component labels and classifications

### Data Mapping

**EegFun.read_eeglab** converts EEGLAB structures to native EegFun types:

- EEGLAB continuous → `ContinuousData`
- EEGLAB epoched → `EpochData`
- EEGLAB ICA → `InfoIca`

All EegFun functions should work with imported EEGLAB data.

## Important Notes

EEGLAB import is functional but considered work-in-progress. The implementation has only been tested with sample datasets from `eeglab/sample_data`, but should handle common use cases.

## Workflow Summary

This demo shows EEGLAB import workflows:

### 1. Import Continuous Data

- Load raw continuous EEGLAB file
- Visualize in databrowser
- Check trigger counts

### 2. Import Epoched Data with ICA

- Load epoched EEGLAB file
- Extract both epoch data and ICA decomposition
- Visualize epochs and ICA components
