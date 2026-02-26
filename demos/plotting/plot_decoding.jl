# Demo: Plotting Decoding Results
# Shows how to visualise MVPA decoding accuracy over time,
# compare individual subjects, and overlay statistical significance.

using EegFun

#######################################################################
# LOAD DATA AND PREPARE FOR DECODING
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# basic preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# epoch
epochs = EegFun.epoch_data(dat, [:trigger1, :trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))

#######################################################################
# RUN DECODING
#######################################################################

# decode condition from EEG patterns
decoded = EegFun.decode(epochs)

#######################################################################
# BASIC DECODING PLOT
#######################################################################

# plot accuracy over time with error shading
EegFun.plot_decoding(decoded)

#######################################################################
# CUSTOM STYLING
#######################################################################

# change colour, line width, title
EegFun.plot_decoding(decoded,
    color = :red,
    linewidth = 3,
    title = "Face vs. Object Decoding"
)

# hide error shading
EegFun.plot_decoding(decoded, show_error = false)

#######################################################################
# MULTI-SUBJECT SUBPLOT GRID
#######################################################################

# when you have a list of decoded results (one per subject),
# each is plotted in its own subplot
# EegFun.plot_decoding(all_decoded, title = "Individual Subjects")

#######################################################################
# DECODING WITH SIGNIFICANCE
#######################################################################

# after running statistics on the decoded data:
# stats = EegFun.test_against_chance(decoded_list, alpha = 0.05,
#              correction_method = :bonferroni)
# EegFun.plot_decoding(grand_avg, stats, show_significance = true)
