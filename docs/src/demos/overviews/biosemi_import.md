# BioSemi Import

This demo demonstrates importing BioSemi BDF files into EegFun.jl.

### What is BioSemi BDF Format?

BioSemi is a popular manufacturer of active electrode EEG systems. The BDF (BioSemi Data Format) is a 24-bit variant of the European Data Format (EDF) designed specifically for BioSemi's ActiveTwo and ActiveOne systems.

**Key features**:

- 24-bit resolution for high-precision recordings
- Continuous raw data storage
- Trigger channel for event markers
- Standard format widely supported across analysis tools

### Import Capabilities

**Data loading**:

- Raw continuous EEG time series (24-bit precision)
- Sampling rate and channel metadata
- Trigger/event markers (encoded in Status channel)
- Automatic channel configuration

**What you need**:

- `.bdf` file (contains all data)
- Electrode layout file matching your cap configuration
- BioSemi systems typically use standard cap layouts (16, 32, 64, 72, or 128 channels)

### Data Mapping

**EegFun.read_raw_data** (or **EegFun.read_bdf**) loads BioSemi data, then **create_eegfun_data** converts to native EegFun types:

- BioSemi BDF → `BiosemiData` (intermediate) → `ContinuousData` (EegFun)
- Triggers extracted from Status channel
- Layout coordinates added separately

All EegFun functions work seamlessly with imported BioSemi data.

## Important Notes

### Layout Files

BioSemi doesn't embed electrode positions in the data file, so you must provide a matching layout. EegFun includes standard BioSemi layouts in `resources/layouts/biosemi/`:

- `biosemi16.csv`
- `biosemi32.csv`
- `biosemi64.csv`
- `biosemi72.csv` (64 + 8 external channels)
- `biosemi128.csv`

Choose the layout matching your recording setup.

### Trigger Encoding

BioSemi encodes triggers in a dedicated Status channel using binary encoding. EegFun automatically extracts these during import, making them available via `trigger_count()` and epoch extraction functions.

## Workflow Summary

This demo shows the complete BioSemi import workflow:

### 1. Load Raw BDF Data

- Use `read_raw_data()` with path to `.bdf` file
- Automatic format detection based on file extension

### 2. Prepare Layout

- Load appropriate BioSemi layout file
- Convert polar coordinates to cartesian for visualization

### 3. Create EegFun Data Structure

- Combine raw data with layout using `create_eegfun_data()`
- Results in `ContinuousData` ready for analysis
