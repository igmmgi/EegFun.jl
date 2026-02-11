# Demo: Condition Operations
# Shows how to combine epoch conditions and compute ERP difference waves.
# Covers condition_combine and condition_difference.

using EegFun

#######################################################################
# 1. CREATE EXAMPLE DATA
#######################################################################

# Load raw data and do minimal preprocessing
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)

dat = EegFun.create_eegfun_data(dat, layout)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Extract epochs for four conditions
epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
    EegFun.EpochCondition(name = "Trigger3", trigger_sequences = [[3]]),
    EegFun.EpochCondition(name = "Trigger4", trigger_sequences = [[5]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
EegFun.baseline!(epochs, (-0.2, 0.0))

# Combine conditions 1+2 into one group and 3+4 into another
epochs_combined = EegFun.condition_combine(epochs, [[1, 2], [3, 4]])


#######################################################################
# 3. CONDITION DIFFERENCE WAVES (IN MEMORY)
#######################################################################

# Create ERPs from the original four conditions
erps = EegFun.average_epochs(epochs)

# Difference wave: condition 1 minus condition 2
diff_waves = EegFun.condition_difference(erps, [(1, 2)])
EegFun.plot_erp(diff_waves, channel_selection = EegFun.channels([:Cz]))

# Multiple difference pairs
diff_waves = EegFun.condition_difference(erps, [(1, 2), (3, 4)])


