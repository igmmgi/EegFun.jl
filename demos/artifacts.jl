using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/julia/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# Plot data in data browser
EegFun.plot_databrowser(dat)

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

# How many vEOG/hEOG onsets were detected?
EegFun.n_values(dat, :is_vEOG)
EegFun.n_values(dat, :is_hEOG)

# We should now be able to see newly calculated EOG signals in 
# the databrowser under "Extra Channels"
EegFun.plot_databrowser(dat)

# Look for extreme values (adds Bool column to dataframe :is_extreme_value_100)
EegFun.is_extreme_value!(dat, 100); # 100 mV criterion

# How many extreme values
EegFun.n_extreme_value(dat, 100)                    # across all electrodes
EegFun.n_extreme_value(dat, 100, mode = :separate)  # separately for each electrode

# Additional examples with channel selection predicate
EegFun.is_extreme_value!(dat, 100; channel_selection = EegFun.channels_not([:Fp1, :Fp2]));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> endswith.(string.(x), "z"));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")), sample_selection = x -> x.sample .> 1000);

# Artifact Detection in Epoched Data
EegFun.trigger_count(dat) # our markers or triggers
EegFun.mark_epoch_windows!(dat, [1, 2, 3, 5, 6], [-0.2, 1.0])

# Can use marked epoch window as a basis for automatic artifact detection
EegFun.is_extreme_value!(dat, 20; sample_selection = EegFun.samples(:epoch_window), channel_out = :is_extreme_value_epoch)

# We can now see the epoch windows and artifacts within them in the databrowser under "Epoch Windows"
EegFun.plot_databrowser(dat)


# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

# Detect bad epochs
# Default is combination of absolute and z-score criteria with z-score based on all epochs using the 
# following metrics: [:variance, :max, :min, :abs, :range, :kurtosis],
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs)
bad_epochs[1]

# But we can change the metrics (here only an absolute criterion is used)
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

# We can inspect the automatic detection visually
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])

# We can also perform visual artifact detection
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))


# Our a combination of both
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)
bad_epochs_manual_visual = EegFun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs[1], dims = (4, 4))


# Artifact Info
# First condition
EegFun.unique_rejections(bad_epochs[1])
EegFun.unique_channels(bad_epochs[1])
EegFun.unique_epochs(bad_epochs[1])
EegFun.get_rejected(bad_epochs[1]) # TODO: this is just the same as unique_rejections?

# repair
# TODO: repair seems broken here?
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :neighbor_interpolation)
# epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :spherical_spline)
# epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :reject)

EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1])
EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1], ylim = (-100, 100))

# reject 
epochs_rejected = EegFun.reject_epochs(epochs, bad_epochs)
EegFun.plot_artifact_rejection(epochs[1], epochs_rejected[1], bad_epochs[1])
