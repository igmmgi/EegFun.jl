"""
Demo: Loading and Processing BioSemi BDF Files

This demo shows how to:
- Load BioSemi .bdf files (raw continuous data format)
- Create EegFun data structures with layouts
- Get the triggers/events
"""

using EegFun

# Load raw BDF data
# read_raw_data automatically detects the .bdf extension
raw_data = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")

# EegFun uses the Julia package: BiosemiDataFormat.jl
# https://github.com/igmmgi/BiosemiDataFormat.jl

raw_data.header.xxx # tab-autocomplete in REPPL
raw_data.data
raw_data.triggers

# Load and prepare electrode layout
# BioSemi systems typically use their standard cap configurations
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
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
