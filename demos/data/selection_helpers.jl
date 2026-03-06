# Demo: Selection Helpers
# Shows how to use EegFun's selection helper functions (channels, conditions,
# participants, epochs, times, samples) for filtering and subsetting data.

using EegFun

#######################################################################
# SETUP
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Create epochs and ERPs
epoch_cfg =
    [EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]]), EegFun.EpochCondition(name = "Cond2", trigger_sequences = [[2]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
erps = EegFun.average_epochs(epochs)


#######################################################################
# CHANNELS() — Channel Selection
#######################################################################

# Select all channels (default)
EegFun.subset(erps[1], channel_selection = EegFun.channels())

# Select specific channels by name
EegFun.subset(erps[1], channel_selection = EegFun.channels(:Cz, :Pz, :Fz))

# Select channels by vector
EegFun.subset(erps[1], channel_selection = EegFun.channels([:Fp1, :Fp2]))



#######################################################################
# TIMES() — Time Window Selection
#######################################################################

# Select all time points (default)
EegFun.subset(erps[1], interval_selection = EegFun.times())

# Select a time window (in seconds)
EegFun.subset(erps[1], interval_selection = EegFun.times(0.0, 0.5))

# Combine channel and time selection
sub = EegFun.subset(erps[1], channel_selection = EegFun.channels(:Cz, :Pz), interval_selection = EegFun.times(0.1, 0.4))


#######################################################################
# EPOCHS() — Epoch Selection
#######################################################################

# Select all epochs (default)
EegFun.all_data(epochs, epoch_selection = EegFun.epochs())

# Select specific epoch indices
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:5))

# Subset epochs from an EpochData object
EegFun.subset(epochs[1], epoch_selection = EegFun.epochs(1:3))


#######################################################################
# SAMPLES() — Sample-level Predicates
#######################################################################

# Custom predicate on metadata columns
EegFun.subset(dat, sample_selection = x -> x.sample .<= 5000)

# Time-based predicate
EegFun.subset(dat, sample_selection = x -> x.time .<= 5.0)


#######################################################################
# COMBINING SELECTIONS
#######################################################################

# Multiple selections at once
sub = EegFun.subset(
    epochs[1],
    channel_selection = EegFun.channels(:Cz),
    epoch_selection = EegFun.epochs(1:3),
    interval_selection = EegFun.times(0.0, 0.5),
)

# Use with plotting
EegFun.plot_erp(erps, channel_selection = EegFun.channels(:Cz, :Pz), interval_selection = EegFun.times(-0.1, 0.8))
