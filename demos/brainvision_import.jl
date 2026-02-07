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

# Load raw BrainVision data
# read_raw_data automatically detects BrainVision files (.vhdr, .eeg, or .vmrk)
# You can specify any of the three files; it will find the others
raw_data = EegFun.read_raw_data("./resources/data/brainvision/example1.vhdr")

# Create EegFun data structure
dat = EegFun.create_eegfun_data(raw_data)

# Check trigger information
println("\nTrigger summary:")
EegFun.trigger_count(dat)


