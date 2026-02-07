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
- Channel information and units
- Reference channel configuration

**What you need**:

- BrainVision triplet (`.vhdr`, `.eeg`, `.vmrk` files in same directory)
- Electrode layout file matching your cap configuration
- Brain Products systems often use actiCap or EasyCap layouts

### Data Mapping

**EegFun.read_raw_data** (or **EegFun.read_brainvision**) loads BrainVision data, then **create_eegfun_data** converts to native EegFun types:

- BrainVision → `BrainVisionData` (intermediate) → `ContinuousData` (EegFun)
- Markers extracted from `.vmrk` file
- Layout coordinates added separately

All EegFun functions work seamlessly with imported BrainVision data.

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
**BioSemi**: Use if you have a custom setup with BioSemi caps

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
