# Demo: Loading and Processing BioSemi BDF Files
# Shows how to load BioSemi .bdf files, create EegFun data structures
# with layouts, and get the triggers/events.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

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
println("\nTrigger summary:")
EegFun.trigger_count(dat) # BioSemi data uses 8 bit trigger lines

# Basic preprocessing
EegFun.highpass_filter!(dat, 1.0) # 1 Hz high-pass
EegFun.rereference!(dat, :avg)    # Average reference

# Open interactive databrowser
EegFun.plot_databrowser(dat);

