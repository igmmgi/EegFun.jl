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

# Load raw BDF data
# read_raw_data automatically detects the .bdf extension
raw_data = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")

# Load and prepare electrode layout
# BioSemi systems typically use their standard cap configurations
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)

# Create EegFun data structure
dat = EegFun.create_eegfun_data(raw_data, layout_file)

# Check trigger information
EegFun.trigger_count(dat)
