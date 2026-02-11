# Demo: Artifact Detection
# Shows automatic epoch rejection using statistical criteria (z-score and absolute threshold),
# how to inspect rejection results, and the available detection metrics.

using EegFun

#######################################################################
# 1. LOAD AND PREPROCESS DATA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)


#######################################################################
# 2. EXTRACT EPOCHS
#######################################################################

epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))

# Apply baseline correction
EegFun.baseline!(epochs)


#######################################################################
# 3. DETECT ARTIFACTS — DEFAULT SETTINGS
#######################################################################

# Default: combined z-score (z > 3) and absolute threshold (> 100 μV)
# Uses all z-score measures: [:variance, :max, :min, :abs, :range, :kurtosis]
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs)

# Print rejection summary
bad_epochs[1]  # display summary for first condition
bad_epochs[2]  # display summary for second condition


#######################################################################
# 4. CUSTOMISE DETECTION CRITERIA
#######################################################################

# Only absolute threshold (no z-score)
bad_abs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

# Only z-score (no absolute threshold)
bad_z = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2, abs_criterion = 0)

# Only variance-based z-score (common for ERP studies)
bad_var = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2.5, z_measures = [:variance])

# Stricter thresholds
bad_strict = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2, abs_criterion = 50)


#######################################################################
# 5. CHANNEL-SPECIFIC DETECTION
#######################################################################

# Detect artifacts only in specific channels
bad_frontal = EegFun.detect_bad_epochs_automatic(
    epochs,
    channel_selection = EegFun.channels([:Fp1, :Fp2, :Fz]),
)


#######################################################################
# 6. INSPECT REJECTION RESULTS
#######################################################################

# Get summary information
EegFun.unique_rejections(bad_epochs[1])  # all unique (channel, epoch) pairs
EegFun.unique_channels(bad_epochs[1])    # which channels had artifacts
EegFun.unique_epochs(bad_epochs[1])      # which epoch indices were flagged
EegFun.get_rejected(bad_epochs[1])       # indices of rejected epochs


#######################################################################
# 7. VISUALISE DETECTION RESULTS
#######################################################################

# Plot detection summary for a condition
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])
EegFun.plot_artifact_detection(epochs[2], bad_epochs[2])


#######################################################################
# 8. INTERACTIVE REVIEW (MANUAL + AUTOMATIC)
#######################################################################

# Interactive grid review — start with automatic detections
# bad_manual = EegFun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs[1], dims = (4, 4))

# Purely manual review (no pre-filled detections)
# bad_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))
