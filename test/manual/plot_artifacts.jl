using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

# Artifacts
artifacts = eegfun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
eegfun.plot_artifact_detection(epochs[1], artifacts[1])
