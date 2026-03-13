# Demo: Time-Frequency Operations
# Shows baseline correction, channel operations, condition operations,
# and grand averaging for time-frequency data.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

#######################################################################
# CREATE EXAMPLE TF DATA
#######################################################################

# Load and preprocess
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Extract epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))

# Compute time-frequency decomposition (Morlet wavelets)
tf_data = EegFun.tf_morlet(epochs, frequencies = 4:1:30, cycles = 3)


#######################################################################
# TF BASELINE CORRECTION
#######################################################################

# Decibel baseline (default, most common)
tf_bl = EegFun.tf_baseline(tf_data[1], (-0.2, 0.0); method = :db)

# Other methods: :absolute, :relative, :relchange, :normchange, :zscore, :percent
tf_pct = EegFun.tf_baseline(tf_data[1], (-0.2, 0.0); method = :percent)

# In-place version
# EegFun.tf_baseline!(tf_data[1], (-0.2, 0.0); method = :db)

# Plot the result
EegFun.plot_tf(tf_bl, channel_selection = EegFun.channels([:Cz]))


#######################################################################
# TF CHANNEL AVERAGE
#######################################################################

# Average specific channel groups (creates new column, preserves originals)
tf_avg = EegFun.channel_average(tf_data[1], channel_selections = [EegFun.channels([:Fz, :Cz, :Pz])], output_labels = [:midline])

# Reduce to only the averaged channel (drops originals)
tf_reduced =
    EegFun.channel_average(tf_data[1], channel_selections = [EegFun.channels([:Fz, :Cz, :Pz])], output_labels = [:midline], reduce = true)


#######################################################################
# TF CHANNEL DIFFERENCE
#######################################################################

# Laterality: left minus right hemisphere
tf_lat = EegFun.channel_difference(
    tf_data[1],
    channel_selection1 = EegFun.channels([:C3]),
    channel_selection2 = EegFun.channels([:C4]),
    channel_out = :laterality,
)


#######################################################################
# TF CONDITION DIFFERENCE
#######################################################################

# Subtract condition 2 from condition 1
tf_diff = EegFun.condition_difference(tf_data, [(1, 2)])

# Multiple differences
# tf_diff = EegFun.condition_difference(tf_data, [(1, 2), (3, 4)])


#######################################################################
# TF CONDITION AVERAGE
#######################################################################

# Average conditions into one
tf_avg_cond = EegFun.condition_average(tf_data, [[1, 2]])


#######################################################################
# TF GRAND AVERAGE
#######################################################################

# Grand average across participants (requires multi-participant data)
# all_tf = EegFun.read_all_data(EegFun.TimeFreqData, "derivatives/tf/tf_morlet")
# grand_avg = EegFun.grand_average(all_tf)

# With condition selection
# grand_avg = EegFun.grand_average(all_tf, condition_selection = EegFun.conditions([1, 2]))
