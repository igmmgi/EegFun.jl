# Demo: Plotting Global Field Power
# Shows how to visualise GFP (Global Field Power) and Global Dissimilarity
# across conditions, with optional ERP trace overlay.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
# LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# extract epochs
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

EegFun.baseline!(epochs, (-0.2, 0.0))
erps = EegFun.average_epochs(epochs)

# plot GFP for a single condition
EegFun.plot_gfp(erps)

