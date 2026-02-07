# BrainVision Import

# BrainVision Import

# BrainVision Import

This demo demonstrates importing BrainVision format files into EegFun.jl, covering the complete workflow from raw data loading through preprocessing and visualization.

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
Demo: Loading and Processing BrainVision Files

This demo shows how to:
- Load BrainVision format files (.vhdr, .eeg, .vmrk)
- Create EegFun data structures with layouts
- Apply basic preprocessing
- Visualize the data
- Work with triggers/events

BrainVision is the data format used by Brain Products GmbH EEG systems.
The format consists of three files:
- .vhdr (header file with metadata)
- .eeg (binary data file)
- .vmrk (marker/trigger file)

Once loaded, all EegFun functions work seamlessly with BrainVision data.
"""

using EegFun

# Step 1: Load raw BrainVision data
# read_raw_data automatically detects BrainVision files (.vhdr, .eeg, or .vmrk)
# You can specify any of the three files; it will find the others
println("Loading BrainVision files...")
raw_data = EegFun.read_raw_data("./resources/data/brainvision/example1.vhdr")

# Alternative: You can also use .eeg or .vmrk extension
# raw_data = EegFun.read_raw_data("./resources/data/brainvision/example1.eeg")

# Step 2: Load and prepare electrode layout
# Brain Products systems often use actiCap or EasyCap configurations
# Choose the layout that matches your recording setup
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")  # Adjust as needed
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
