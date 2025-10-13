using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

# Calculate EOG signals
eegfun.channel_difference!(dat, channel_selection1 = eegfun.channels([:Fp1, :Fp2]), channel_selection2 = eegfun.channels([:IO1, :IO2]), channel_out = :vEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.channel_difference!(dat, channel_selection1 = eegfun.channels([:F9]),        channel_selection2 = eegfun.channels([:F10]),       channel_out = :hEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

dat.data[!, :is_vEOG]
sum(dat.data[!, :is_vEOG])
dat.data[!, :is_hEOG]
sum(dat.data[!, :is_hEOG])

eegfun.is_extreme_value!(dat, 100);
sum(dat.data[!, :is_extreme_value_100])

eegfun.n_extreme_value(dat, 100)
eegfun.n_extreme_value(dat, 100, mode = :separate)


# artifact detection in epochs
# some epoched data
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[3]]),
]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)


bad_epochs = eegfun.detect_bad_epochs(epochs)
eegfun.get_rejected_epochs(bad_epochs[1])
eegfun.get_rejected_epochs(bad_epochs)


bad_epochs = eegfun.detect_bad_epochs(epochs[1])
eegfun.get_rejected_epochs(bad_epochs)
bad_epochs = eegfun.detect_bad_epochs(epochs)

eegfun.plot_artifact_detection(epochs[1], bad_epochs[1])

# repair
epochs_repaired = eegfun.repair_artifacts(epochs[1], bad_epochs[1], method = :neighbor_interpolation)
epochs_repaired = eegfun.repair_artifacts(epochs[1], bad_epochs[1], method = :spherical_spline)
epochs_repaired = eegfun.repair_artifacts(epochs[1], bad_epochs[1], method = :reject)

eegfun.plot_artifact_repair(epochs[1], epochs_repaired, bad_epochs[1])






