using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

# Calculate EOG signals
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:Fp1, :Fp2]),
    channel_selection2 = eegfun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:F9]),
    channel_selection2 = eegfun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

dat.data[!, :is_vEOG]
sum(dat.data[!, :is_vEOG])
dat.data[!, :is_hEOG]
sum(dat.data[!, :is_hEOG])

eegfun.is_extreme_value!(dat, 100);
sum(dat.data[!, :is_extreme_value_100])

eegfun.n_extreme_value(dat, 100)
eegfun.n_extreme_value(dat, 50, mode = :separate)


# add bool columns to the data frame
eegfun.is_extreme_value!(dat, 100);
eegfun.is_extreme_value!(dat, 100; channel_selection = eegfun.channels_not([:Fp1, :Fp2]));
eegfun.is_extreme_value!(dat, 100; channel_selection = x -> endswith.(string.(x), "z"));
eegfun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")));
eegfun.is_extreme_value!(
    dat,
    100;
    channel_selection = x -> .!(endswith.(string.(x), "z")),
    sample_selection = x -> x.sample .> 1000,
);

# retrun count of extreme values at specific electrodes at different thresholds
eegfun.n_extreme_value(dat, 100)
eegfun.n_extreme_value(dat, 100, channel_selection = eegfun.channels([:Fp1, :Fp2])) # count extreme values at Fp1 at 100 uV threshold
eegfun.n_extreme_value(dat, 100, channel_selection = x -> endswith.(string.(x), "z"))
eegfun.n_extreme_value(
    dat,
    100,
    channel_selection = x -> .!(endswith.(string.(x), "z")),
    sample_selection = x -> x.sample .< 10,
)





# artifact detection in epochs
# some epoched data
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[3]]),
]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)


bad_epochs = eegfun.detect_bad_epochs_automatic(epochs)

eegfun.plot_artifact_detection(epochs[1], bad_epochs[1])
eegfun.plot_artifact_detection(epochs[1], bad_epochs[1], ylimits = (-100, 100))
eegfun.plot_artifact_detection(epochs[1], bad_epochs[1], xlimits = (-5, 15), ylimits = (-100, 100))


eegfun.get_rejected_epochs(bad_epochs[1])
eegfun.get_rejected_epochs(bad_epochs)

eegfun.unique_rejections(bad_epochs[1].rejected_epochs)
eegfun.unique_channels(bad_epochs[1].rejected_epochs)
eegfun.unique_epochs(bad_epochs[1].rejected_epochs)


bad_epochs_automatic = eegfun.detect_bad_epochs_automatic(epochs[1])
eegfun.get_rejected_epochs(bad_epochs_automatic)
eegfun.plot_artifact_detection(epochs[1], bad_epochs)

bad_epochs_manual = eegfun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))
bad_epochs_manual = eegfun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs, dims = (4, 4))


eegfun.unique_rejections(bad_epochs.rejected_epochs)
eegfun.unique_channels(bad_epochs.rejected_epochs)
eegfun.unique_epochs(bad_epochs.rejected_epochs)

# repair
bad_epochs_automatic = eegfun.detect_bad_epochs_automatic(epochs)

epochs_repaired = eegfun.repair_artifacts(epochs, bad_epochs_automatic, method = :neighbor_interpolation)
epochs_repaired = eegfun.repair_artifacts(epochs, bad_epochs_automatic, method = :spherical_spline)
epochs_repaired = eegfun.repair_artifacts(epochs, bad_epochs_automatic, method = :reject)

eegfun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1])
eegfun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1], ylimits = (-100, 100))











bad_epochs_manual = eegfun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))
bad_epochs_manual = eegfun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4), artifact_info = bad_epochs_automatic)
bad_epochs_manual = eegfun.detect_bad_epochs_interactive(
    epochs[1],
    dims = (4, 4),
    artifact_info = bad_epochs_automatic,
    colormap = :seaborn_colorblind,
)
