"""
Demo: Artifact Detection and Management

This demo showcases EegFun.jl's tools for identifying and handling artifacts, including:
- EOG (eye movement) onset detection
- Extreme value detection (peak-to-peak/rejection thresholds)
- Automatic epoch rejection based on statistical measures (Z-scores)

Artifact management is a critical step in EEG preprocessing to ensure 
that eye blinks, muscle activity, and technical noise do not bias your results.
"""

using EegFun

# 1. Load data
# We'll use the sample BioSemi data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Basic preprocessing before artifact detection
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# 2. EOG Detection (on continuous data)
# Useful for identifying blinks or saccades if you have EOG channels
# Here we demonstrate the principle using a frontal channel as a proxy
println("Detecting EOG-like onsets...")
eog_cfg = Dict(
    "vEOG_criterion" => 50.0,
    "hEOG_criterion" => 50.0,
    "vEOG_channels" => [["Fp1"], ["Fp1"], ["vEOG_proxy"]],
    "hEOG_channels" => [["F7"], ["F8"], ["hEOG_proxy"]],
)
# Note: In a real scenario, you'd have dedicated EOG channels
# detect_eog_signals!(dat, eog_cfg)

# 3. Extreme Value Detection
# Flag points exceeding a certain voltage (e.g., 100 μV)
println("Flagging extreme values (> 100 μV)...")
EegFun.is_extreme_value!(dat, 100, channel_out = :is_extreme)
n_artifacts = EegFun.n_values(dat, :is_extreme)
println("Found $n_artifacts samples exceeding threshold.")

# 4. Automatic Epoch Rejection
# Typically done after segmenting the data
println("\nExtracting epochs for automatic rejection demo...")
epoch_cfg = [EegFun.EpochCondition(name = "Stimulus", trigger_sequences = [[1], [2]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))

# Automatically detect bad epochs using z-scores across multiple measures
# Measures include: variance, max, min, abs, range, and kurtosis
println("Running automatic rejection...")
rejection_info = EegFun.detect_bad_epochs_automatic(
    epochs[1],
    z_criterion = 3.0,     # Reject if > 3 standard deviations
    abs_criterion = 120.0,  # Also reject if any point > 120 μV
)

println("Rejection Summary for $(epochs[1].condition_name):")
println("- Total Epochs: $(rejection_info.n_epochs)")
println("- Rejected: $(rejection_info.n_artifacts)")
println("- Rejected Epoch IDs: $(rejection_info.rejected)")

# 5. Visualization
# The databrowser can often highlight these flagged artifacts if configured
EegFun.plot_databrowser(dat)
