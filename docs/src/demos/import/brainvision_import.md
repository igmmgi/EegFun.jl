# BrainVision Import

This demo demonstrates importing BrainVision format files into EegFun.jl.

### What is BrainVision Format?


BrainVision is the data format used by Brain Products GmbH EEG systems (e.g., BrainAmp, actiCHamp, LiveAmp). It's a flexible, text-header based format widely used in research.

**The format consists of three files**:

- `.vhdr` - Header file (ASCII text) containing metadata, channel information, and data format specifications
- `.eeg` - Binary data file containing the actual EEG samples
- `.vmrk` - Marker file (ASCII text) with trigger/event information

**Key features**:

- Human-readable header and marker files
- Flexible data encoding (16-bit, 32-bit integer or float)
- Separate marker file for triggers and annotations
- Wide support across analysis platforms

### Import Capabilities

**Data loading**:

- Raw continuous EEG time series
- Sampling rate and channel metadata
- Trigger/event markers with descriptions

**What you need**:

- BrainVision triplet (`.vhdr`, `.eeg`, `.vmrk` files in same directory)
- Electrode layout file matching your cap configuration
- Brain Products systems often use actiCap or EasyCap layouts

### Data Mapping

**EegFun.read_raw_data** loads BrainVision data, then **create_eegfun_data** converts to native EegFun types:

- BrainVision → `BrainVisionData` (intermediate) → `ContinuousData` (EegFun)
- Markers extracted from `.vmrk` file
- Layout coordinates added separately


## Important Notes

### File Organization

All three files must be in the same directory with matching base names:

```text
my_recording.vhdr
my_recording.eeg
my_recording.vmrk
```

You can specify any of the three files to `read_raw_data()` - it will automatically find the others.

### Layout Files

BrainVision doesn't embed 3D electrode positions, so you must provide a matching layout. Depending on your cap manufacturer:

**Brain Products actiCap**: Check `resources/layouts/brainproducts_acticap/`  
**EasyCap**: Check `resources/layouts/easycap/`  

Choose the layout matching your recording configuration.

### Marker Types

BrainVision supports rich marker information including:

- Stimulus markers (trigger codes)
- Response markers (button presses)
- Comments and annotations
- Segment boundaries

EegFun imports all markers and makes them available for epoch extraction and analysis.

## Workflow Summary

This demo shows the complete BrainVision import workflow:

### 1. Load Raw BrainVision Data

- Use `read_raw_data()` with path to any file in the triplet (`.vhdr`, `.eeg`, or `.vmrk`)
- Automatic format detection based on file extension
- All three files loaded automatically

### 2. Prepare Layout

- Load appropriate layout file for your cap configuration
- Convert polar coordinates to cartesian for visualization

### 3. Create EegFun Data Structure

- Combine raw data with layout using `create_eegfun_data()`
- Results in `ContinuousData` ready for analysis


## Code Examples

::: details Show Code

```julia
"""
Demo: Loading and Processing BrainVision Files

This demo shows how to:
- Load BrainVision format files (.vhdr, .eeg, .vmrk)
- Create EegFun data structures 
- Apply basic preprocessing
- Get the triggers/events

The format consists of three files:
- .vhdr (header file with metadata)
- .eeg (binary data file)
- .vmrk (marker/trigger file)
"""
# Demo: BrainVision Import
# Shows how to load BrainVision .vhdr files and create EegFun data structures.

using EegFun

# Load raw BrainVision data
# read_raw_data automatically detects BrainVision files (.vhdr, .eeg, or .vmrk)
# You can specify any of the three files; it will find the others, assuming 
# they are in the same directory
dat = EegFun.read_raw_data("./resources/data/brainvision/example1.vhdr")

# EegFun uses the Julia package: BrainVisionDataFormat.jl
# https://github.com/igmmgi/BrainVisionDataFormat.jl
raw_data.data
raw_data.header
raw_data.markers
raw_data.markers[1]
raw_data.markers[2]
raw_data.markers[3]
# and so on

# Create EegFun data structure
dat = EegFun.create_eegfun_data(raw_data)

EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.all_labels(dat)
EegFun.channel_labels(dat)
EegFun.meta_labels(dat)
EegFun.extra_labels(dat)

# Check trigger information
# EegFun creates a hash of the marker names for the numeric trigger column
# triggers_info contains the mapping of trigger names to numeric values
println("\nTrigger summary:")
EegFun.trigger_count(dat)


```

:::

## See Also

- [API Reference](../../reference/index.md)
