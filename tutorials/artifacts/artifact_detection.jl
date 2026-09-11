# Demo: Artifact Detection, Repair, and Rejection
# Comprehensive walkthrough of artifact handling in EegFun:
# continuous-data artifact detection, EOG channels, epoched artifact detection,
# interactive review, artifact repair, and epoch rejection.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
# LOAD AND PREPROCESS DATA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Browse raw data (press "i" for interaction help)
EegFun.plot_databrowser(dat)


#######################################################################
# EOG CHANNEL CREATION AND ONSET DETECTION
#######################################################################

# Calculate vertical and horizontal EOG channels
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
)

# Automatically detect EOG onsets and mark intervals
EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 50, :hEOG, :is_hEOG)
EegFun.mark_epoch_intervals!(dat, :is_vEOG, [-0.05, 0.4])

# How many EOG onsets were detected?
EegFun.n_values(dat, :is_vEOG)
EegFun.n_values(dat, :is_hEOG)


#######################################################################
# CONTINUOUS-DATA ARTIFACT DETECTION
#######################################################################

# Extreme value detection (adds Bool column :is_extreme_value_100)
EegFun.is_extreme_value!(dat, 100)
EegFun.n_extreme_value(dat, 100)                    # count across all channels
EegFun.n_extreme_value(dat, 100, mode = :separate)  # per channel

# Step artifact detection — catches sudden voltage jumps (e.g., cable disconnections)
EegFun.is_step_value!(dat, 50.0)
EegFun.n_step_value(dat, 50.0)                      # count across all channels
EegFun.n_step_value(dat, 50.0, mode = :separate)    # per channel

# Channel-specific detection examples
EegFun.is_extreme_value!(dat, 100; channel_selection = EegFun.channels_not([:Fp1, :Fp2]))
EegFun.is_step_value!(dat, 50.0; channel_selection = EegFun.channels([:Fp1, :Fp2]), channel_out = :is_step_frontal)

# Epoch interval marking for continuous data
EegFun.mark_epoch_intervals!(dat, [1, 2, 3, 5, 6], [-0.2, 1.0])
EegFun.is_extreme_value!(dat, 50; sample_selection = EegFun.samples(:epoch_interval), channel_out = :is_extreme_value_epoch)

# Browse data with artifact annotations visible
# EegFun.plot_databrowser(dat)


#######################################################################
# EXTRACT AND BASELINE EPOCHS
#######################################################################

epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
EegFun.baseline!(epochs)


#######################################################################
# AUTOMATIC EPOCH ARTIFACT DETECTION
#######################################################################

# Default: combined z-score (z > 3) and absolute threshold (> 100 μV)
# Uses all z-score measures: [:variance, :max, :min, :abs, :range, :kurtosis]
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs)

# Print rejection summary
bad_epochs[1]  # first condition
bad_epochs[2]  # second condition

# Only absolute threshold (no z-score)
bad_abs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

# Only z-score (no absolute threshold)
bad_z = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2, abs_criterion = 0)

# Only variance-based z-score (common for ERP studies)
bad_var = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2.5, z_measures = [:variance])

# Stricter thresholds
bad_strict = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2, abs_criterion = 50)

# Channel-specific detection
bad_frontal = EegFun.detect_bad_epochs_automatic(epochs, channel_selection = EegFun.channels([:Fp1, :Fp2, :Fz]))


#######################################################################
# INSPECT REJECTION RESULTS
#######################################################################

# Get summary information
EegFun.unique_rejections(bad_epochs[1])  # all unique (channel, epoch) pairs
EegFun.unique_channels(bad_epochs[1])    # which channels had artifacts
EegFun.unique_epochs(bad_epochs[1])      # which epoch indices were flagged
EegFun.get_rejected(bad_epochs[1])       # indices of rejected epochs


#######################################################################
# VISUALISE DETECTION RESULTS
#######################################################################

# Plot detection summary for each condition
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])
EegFun.plot_artifact_detection(epochs[2], bad_epochs[2])


#######################################################################
# INTERACTIVE REVIEW (MANUAL + AUTOMATIC)
#######################################################################

# Interactive grid review — start with automatic detections
# bad_manual = EegFun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs[1], dims = (4, 4))

# Purely manual review (no pre-filled detections)
# bad_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))


#######################################################################
# ARTIFACT REPAIR
#######################################################################

bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

# Repair via neighbor interpolation
epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :neighbor_interpolation)
# epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :spherical_spline)
# epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :reject)

# Visualise repair results
EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1])
EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1], ylim = (-100, 100))


#######################################################################
# EPOCH REJECTION
#######################################################################

epochs_rejected = EegFun.reject_epochs(epochs, bad_epochs)
EegFun.plot_artifact_rejection(epochs[1], epochs_rejected[1], bad_epochs[1])
