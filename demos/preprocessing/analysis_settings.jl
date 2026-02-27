# Demo: Analysis Settings
# Shows how to save and replay preprocessing settings
# (filter, rereference, ICA) using AnalysisSettings.

using EegFun

#######################################################################
# SETUP
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
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
settings = EegFun.AnalysisSettings(hp_filter = 1.0, lp_filter = 30.0, reference = :avg)

# Apply settings to fresh data (non-mutating — returns a processed copy)
dat_fresh = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
dat_fresh = EegFun.create_eegfun_data(dat_fresh, layout)

dat_processed = EegFun.apply_analysis_settings(dat_fresh, settings)
dat_processed.analysis_info

# Or apply in-place (mutating)
# EegFun.apply_analysis_settings!(dat_fresh, settings)
