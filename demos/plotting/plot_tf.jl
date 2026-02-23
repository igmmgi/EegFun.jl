# Demo: Plotting Time-Frequency Data
# Shows how to plot time-frequency decomposition results as heatmaps,
# with options for baseline correction, colour maps, and log-scaled axes.

using EegFun

#######################################################################
# 1. LOAD DATA AND RUN TF DECOMPOSITION
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.5)
EegFun.lowpass_filter!(dat, 100.0)

epochs = EegFun.epoch_data(dat, [:trigger1], (-0.5, 1.0))
EegFun.baseline!(epochs, (-0.2, 0.0))

# Morlet wavelet decomposition
tf_data = EegFun.tf_morlet(epochs, 2:2:60)

#######################################################################
# 2. BASIC TIME-FREQUENCY PLOT
#######################################################################

# plot first available channel
EegFun.plot_tf(tf_data)

#######################################################################
# 3. SPECIFIC CHANNEL
#######################################################################

EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz))

#######################################################################
# 4. WITH BASELINE CORRECTION
#######################################################################

# apply dB baseline on the fly
EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz), baseline_interval = (-0.3, 0.0), baseline_method = :db)

# percentage change baseline
EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz), baseline_interval = (-0.3, 0.0), baseline_method = :percent)

#######################################################################
# 5. LOG-SCALED FREQUENCY AXIS
#######################################################################

EegFun.plot_tf(tf_data, channel_selection = EegFun.channels(:Cz), baseline_interval = (-0.3, 0.0), ylogscale = true)

#######################################################################
# 6. CUSTOM COLOUR MAP AND RANGE
#######################################################################

EegFun.plot_tf(
    tf_data,
    channel_selection = EegFun.channels(:Cz),
    baseline_interval = (-0.3, 0.0),
    colormap = :RdBu,
    colorrange = (-3.0, 3.0),
    title = "Alpha/Beta ERD",
)
