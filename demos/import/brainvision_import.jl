# Demo: Loading and Processing BrainVision Files
# Shows how to load BrainVision format files (.vhdr, .eeg, .vmrk),
# create EegFun data structures, apply basic preprocessing,
# and get the triggers/events.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.vhdr")

using EegFun

# Load raw BrainVision data
# read_raw_data automatically detects BrainVision files (.vhdr, .eeg, or .vmrk)
# You can specify any of the three files; it will find the others, assuming 
# they are in the same directory
dat = EegFun.read_raw_data(EegFun.example_path("data/brainvision/example1.vhdr"))

# EegFun uses the Julia package: BrainVisionDataFormat.jl
# https://github.com/igmmgi/BrainVisionDataFormat.jl
dat.data
dat.header

# Create EegFun data structure
eeg = EegFun.create_eegfun_data(dat)

EegFun.all_data(eeg)
EegFun.meta_data(eeg)
EegFun.all_labels(eeg)
EegFun.channel_labels(eeg)
EegFun.meta_labels(eeg)
EegFun.extra_labels(eeg)

# Check trigger information
# EegFun creates a hash of the marker names for the numeric trigger column
# triggers_info contains the mapping of trigger names to numeric values
println("\nTrigger summary:")
EegFun.trigger_count(eeg)


