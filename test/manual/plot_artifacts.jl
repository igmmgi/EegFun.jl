using EegFun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

# Artifacts
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
EegFun.plot_artifact_detection(epochs[1], artifacts[1])
