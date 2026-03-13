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

**What you need**:

- `.bdf` file (contains all data)
- Electrode layout file matching your cap configuration
- BioSemi systems typically use standard cap layouts (16, 32, 64, 72, or 128 channels)

### Data Mapping

**EegFun.read_raw_data** loads BioSemi data, then **create_eegfun_data** converts to native EegFun types:

- BioSemi BDF → `BiosemiData` (intermediate) → `ContinuousData` (EegFun)
- Triggers extracted from Status channel
- Layout coordinates added separately

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

### Load Raw BDF Data

- Use `read_raw_data()` with path to `.bdf` file
- Automatic format detection based on file extension

### Prepare Layout

- Load appropriate BioSemi layout file
- Convert polar coordinates to cartesian for visualization

### Create EegFun Data Structure

- Combine raw data with layout using `create_eegfun_data()`
- Results in `ContinuousData` ready for analysis


## Code Examples

::: details Show Code

```julia
# Demo: Loading and Processing BioSemi BDF Files
# Shows how to load BioSemi .bdf files, create EegFun data structures
# with layouts, and get the triggers/events.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

# Load raw BDF data
# read_raw_data automatically detects the .bdf extension
raw_data = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))

# EegFun uses the Julia package: BiosemiDataFormat.jl
# https://github.com/igmmgi/BiosemiDataFormat.jl


raw_data.header.channel_labels  # electrode labels
raw_data.header.sample_rate     # sampling rate
raw_data.header.num_channels    # number of channels
raw_data.data
raw_data.triggers

# Load and prepare electrode layout
# BioSemi systems typically use their standard cap configurations
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)

layout.data # DataFrame with labels and positions

# Create EegFun data structure
dat = EegFun.create_eegfun_data(raw_data, layout)

EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.all_labels(dat)
EegFun.channel_labels(dat)
EegFun.meta_labels(dat)
EegFun.extra_labels(dat) # empty

# Check trigger information
EegFun.trigger_count(dat) # BioSemi data uses 8 bit trigger lines
```

:::

## See Also

- [API Reference](../../reference/index.md)
