using EegFun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.read_bdf(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);

EegFun.rereference!(dat, :avg)
EegFun.filter_data!(dat, "hp", 1)

# Calculate EOG signals
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

dat.data[!, :is_vEOG]
sum(dat.data[!, :is_vEOG])
dat.data[!, :is_hEOG]
sum(dat.data[!, :is_hEOG])

EegFun.is_extreme_value!(dat, 100);
sum(dat.data[!, :is_extreme_value_100])

EegFun.n_extreme_value(dat, 100)
EegFun.n_extreme_value(dat, 50, mode = :separate)


# add bool columns to the data frame
EegFun.is_extreme_value!(dat, 100);
EegFun.is_extreme_value!(dat, 100; channel_selection = EegFun.channels_not([:Fp1, :Fp2]));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> endswith.(string.(x), "z"));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")));
EegFun.is_extreme_value!(
    dat,
    100;
    channel_selection = x -> .!(endswith.(string.(x), "z")),
    sample_selection = x -> x.sample .> 1000,
);

# retrun count of extreme values at specific electrodes at different thresholds
EegFun.n_extreme_value(dat, 100)
EegFun.n_extreme_value(dat, 100, channel_selection = EegFun.channels([:Fp1, :Fp2])) # count extreme values at Fp1 at 100 uV threshold
EegFun.n_extreme_value(dat, 100, channel_selection = x -> endswith.(string.(x), "z"))
EegFun.n_extreme_value(
    dat,
    100,
    channel_selection = x -> .!(endswith.(string.(x), "z")),
    sample_selection = x -> x.sample .< 10,
)

# artifact detection in epochs
# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[3]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

bad_epochs = EegFun.detect_bad_epochs_automatic(epochs)


# determine which channels can be repaired
test = EegFun.channel_repairable(bad_epochs, dat.layout)
EegFun.channel_repairable!(bad_epochs, dat.layout)

epochs_repaired = EegFun.repair_artifacts(epochs, test, method = :neighbor_interpolation)




bad_epochs
bad_epochs[1]

EegFun.unique_rejections(bad_epochs)
EegFun.unique_rejections(bad_epochs[1])
EegFun.unique_channels(bad_epochs)
EegFun.unique_channels(bad_epochs[1])
EegFun.unique_epochs(bad_epochs)
EegFun.unique_epochs(bad_epochs[1])


bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 0)
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])

bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 150)
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])

bad_epochs =
    EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 1, abs_criterion = 0, z_measures = [:variance, :range])
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1], ylim = (-100, 100))



EegFun.get_rejected(bad_epochs[1])
EegFun.get_rejected(bad_epochs)

EegFun.unique_rejections(bad_epochs[1].rejected)
EegFun.unique_channels(bad_epochs[1].rejected)
EegFun.unique_epochs(bad_epochs[1].rejected)

# automatic vs. and/or visual/manual
EegFun.get_rejected(bad_epochs)
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])

bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 200, z_criterion = 0)
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs[1], dims = (4, 4))
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(
    epochs[1],
    artifact_info = bad_epochs[1],
    dims = (4, 4),
    ylim = (-100, 100),
    xlim = (-1, 2),
)


EegFun.unique_rejections(bad_epochs[1].rejected)
EegFun.unique_rejections(bad_epochs)
EegFun.unique_rejections(bad_epochs[1])

EegFun.unique_channels(bad_epochs[1].rejected)
EegFun.unique_channels(bad_epochs)
EegFun.unique_channels(bad_epochs[1])

EegFun.unique_epochs(bad_epochs[1].rejected)
EegFun.unique_epochs(bad_epochs)
EegFun.unique_epochs(bad_epochs[1])

# repair
bad_epochs_automatic = EegFun.detect_bad_epochs_automatic(epochs)

epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs_automatic, method = :neighbor_interpolation)
epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs_automatic, method = :spherical_spline)
epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs_automatic, method = :reject)

EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1])
EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1], ylim = (-100, 100))






bad_epochs_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4), artifact_info = bad_epochs_automatic)
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(
    epochs[1],
    dims = (4, 4),
    artifact_info = bad_epochs_automatic,
    colormap = :seaborn_colorblind,
)

# reject 
epochs_rejected = EegFun.reject_epochs(epochs, bad_epochs_automatic)
EegFun.plot_artifact_rejection(epochs[1], epochs_rejected[1], bad_epochs_automatic[1])
