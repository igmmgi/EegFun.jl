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

using EegFun

# Load raw BrainVision data
# read_raw_data automatically detects BrainVision files (.vhdr, .eeg, or .vmrk)
# You can specify any of the three files; it will find the others, assuming 
# they are in the same directory
raw_data = EegFun.read_raw_data("./resources/data/brainvision/example1.vhdr")

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


