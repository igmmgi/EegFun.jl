# Demo: Analysis Settings
# Shows how to save and replay preprocessing settings
# (filter, rereference, ICA) using AnalysisSettings.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
# SETUP
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)


#######################################################################
# INSPECTING ANALYSIS INFO
#######################################################################

# Every EegFun data object tracks its preprocessing state
dat.analysis_info  # shows reference, filter info

# Apply some preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.rereference!(dat, :avg)

# Check that analysis info is updated
dat.analysis_info  # now shows hp_filter = 0.1, reference = :avg


#######################################################################
# USING ANALYSIS SETTINGS
#######################################################################

# AnalysisSettings stores a recipe for preprocessing
settings = EegFun.AnalysisSettings(1.0, 30.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[])

# Apply settings to fresh data (non-mutating — returns a processed copy)
dat_fresh = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
dat_fresh = EegFun.create_eegfun_data(dat_fresh, layout)

dat_processed = EegFun.apply_analysis_settings(dat_fresh, settings)
dat_processed.analysis_info

# Or apply in-place (mutating)
# EegFun.apply_analysis_settings!(dat_fresh, settings)
