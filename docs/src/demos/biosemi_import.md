# BioSemi Import

# BioSemi Import

# BioSemi Import

This demo demonstrates importing BioSemi BDF files into EegFun.jl, covering the complete workflow from raw data loading through preprocessing and visualization.

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

### 4. Basic Preprocessing

- Apply common average reference
- High-pass filter to remove slow drifts
- Additional preprocessing as needed

### 5. Visualization and Analysis

- Browse data with `plot_databrowser()`
- Extract epochs around triggers
- Create ERPs and compare conditions


## Code Examples

::: details Show Code

```julia
"""
Demo: Loading and Processing BioSemi BDF Files

This demo shows how to:
- Load BioSemi .bdf files (raw continuous data format)
- Create EegFun data structures with layouts
- Apply basic preprocessing
- Visualize the data
- Work with triggers/events

BioSemi is a popular EEG system manufacturer. The .bdf format is the 
BioSemi Data Format, a 24-bit variant of the European Data Format (EDF).

Once loaded, all EegFun functions work seamlessly with BioSemi data.
"""

using EegFun

# Step 1: Load raw BDF data
# read_raw_data automatically detects the .bdf extension
println("Loading BioSemi .bdf file...")
raw_data = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")

# Step 2: Load and prepare electrode layout
# BioSemi systems typically use their standard cap configurations
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)

# Step 3: Create EegFun data structure
println("Creating EegFun data structure...")
dat = EegFun.create_eegfun_data(raw_data, layout_file)

# Step 4: Check trigger information
println("\nTrigger summary:")
EegFun.trigger_count(dat)

# Step 5: Basic preprocessing
println("\nApplying preprocessing...")
EegFun.rereference!(dat, :avg)  # Common average reference
EegFun.highpass_filter!(dat, 0.1)  # Remove slow drifts

# Step 6: Visualize the data
println("\nOpening data browser...")
EegFun.plot_databrowser(dat)

# Optional: Extract and visualize epochs
println("\nExtracting epochs...")
epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]

epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
EegFun.baseline!(epochs, (-0.2, 0.0))

# Visualize epochs
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Cz]))

# Compare conditions
erps = EegFun.average_epochs(epochs)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz]))
```

:::

## See Also

- [API Reference](../reference/index.md)
