# Demo: Plotting Decoding Results
# Shows how to visualise MVPA decoding accuracy over time,
# compare individual subjects, and overlay statistical significance.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# LOAD DATA AND PREPARE FOR DECODING
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# basic preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# extract epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-1, 2))
EegFun.baseline!(epochs, (-0.2, 0.0))

# decode condition from EEG patterns
decoded = EegFun.decode_libsvm(epochs)

# plot accuracy over time with error shading
EegFun.plot_decoding(decoded)
