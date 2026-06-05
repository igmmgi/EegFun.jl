# Comprehensive Data Import

## Overview

Demonstrates Comprehensive Data Import functionality.


## Code Examples

::: details Show Code

```julia
# Demo: Importing and Preprocessing Various Data Formats
# Shows how to load different file formats, create EegFun data structures,
# apply basic preprocessing (filter, rereference), and visualize.

using EegFun

# ==============================================================================
# 1. Biosemi (.bdf)
# ==============================================================================
# Load raw BDF data
raw_bdf = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))

# Load a standard Biosemi 72-channel layout for 3D/2D visualization
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)

# Create EegFun data structure
dat_bdf = EegFun.create_eegfun_data(raw_bdf, layout)

# Basic preprocessing
EegFun.highpass_filter!(dat_bdf, 1.0) # 1 Hz high-pass
EegFun.rereference!(dat_bdf, :avg)    # average reference

# Check triggers and display
println("BDF Triggers: ", EegFun.trigger_count(dat_bdf))
# EegFun.plot_databrowser(dat_bdf)

# ==============================================================================
# 2. European Data Format (.edf)
# ==============================================================================
raw_edf = EegFun.read_raw_data(EegFun.example_path("data/edf/test.edf"))

# Using the fallback auto-generated layout (if you don't have a specific .csv layout)
dat_edf = EegFun.create_eegfun_data(raw_edf)

EegFun.highpass_filter!(dat_edf, 1.0)
EegFun.rereference!(dat_edf, :avg)

println("EDF Triggers: ", EegFun.trigger_count(dat_edf))
# EegFun.plot_databrowser(dat_edf)

# ==============================================================================
# 3. BrainVision (.vhdr)
# ==============================================================================
raw_bv = EegFun.read_raw_data(EegFun.example_path("data/brainvision/example1.vhdr"))
dat_bv = EegFun.create_eegfun_data(raw_bv)

EegFun.highpass_filter!(dat_bv, 1.0)
EegFun.rereference!(dat_bv, :avg)

println("BrainVision Triggers: ", EegFun.trigger_count(dat_bv))
# EegFun.plot_databrowser(dat_bv)

# ==============================================================================
# 4. Functional Image Format (.fif)
# ==============================================================================
raw_fif = EegFun.read_raw_data(EegFun.example_path("data/fif/test_raw.fif"))
dat_fif = EegFun.create_eegfun_data(raw_fif)

EegFun.highpass_filter!(dat_fif, 1.0)
EegFun.rereference!(dat_fif, :avg)

println("FIF Triggers: ", EegFun.trigger_count(dat_fif))
# EegFun.plot_databrowser(dat_fif)

# ==============================================================================
# 5. Extensible Data Format (.xdf)
# ==============================================================================
raw_xdf = EegFun.read_raw_data(EegFun.example_path("data/xdf/test1.xdf"))
dat_xdf = EegFun.create_eegfun_data(raw_xdf)

EegFun.highpass_filter!(dat_xdf, 1.0)
EegFun.rereference!(dat_xdf, :avg)

println("XDF Triggers: ", EegFun.trigger_count(dat_xdf))
# EegFun.plot_databrowser(dat_xdf)

# ==============================================================================
# 6. EEGLAB (.set) & FieldTrip (.mat)
# ==============================================================================
# EEGLAB and FieldTrip data is read directly into the EegFunData structure 
# (bypassing the `create_eegfun_data` step) because they are typically
# preprocessed data saved natively in MATLAB.
# Therefore, we can often skip the basic preprocessing steps.

dat_eeglab = EegFun.read_raw_data(EegFun.example_path("data/eeglab/eeglab_data.set"))
println("EEGLAB Triggers: ", EegFun.trigger_count(dat_eeglab))
# EegFun.plot_databrowser(dat_eeglab)

dat_fieldtrip = EegFun.read_raw_data(EegFun.example_path("data/fieldtrip/continuous.mat"))
println("FieldTrip Triggers: ", EegFun.trigger_count(dat_fieldtrip))
# EegFun.plot_databrowser(dat_fieldtrip)
```

:::

## See Also

- [API Reference](../../reference/index.md)
