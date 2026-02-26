# Demo: Plotting Statistical Test Results
# Shows how to visualise analytic and permutation test results with
# ERP waveforms, difference waves, t-values, and significance markers.

using EegFun

#######################################################################
# LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

epochs = EegFun.epoch_data(dat, [:trigger1, :trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))

#######################################################################
# RUN STATISTICAL TEST
#######################################################################

# prepare for comparison and run analytic test
prepared = EegFun.prepare_statistics(epochs)
result = EegFun.analytic_test(prepared, correction_method = :bonferroni)

#######################################################################
# BASIC PLOT — ERP WAVEFORMS (single channel)
#######################################################################

# plot ERP waveforms for both conditions at a channel
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true)

#######################################################################
# ADD DIFFERENCE WAVE
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true, plot_difference = true)

#######################################################################
# ADD SIGNIFICANCE MARKERS
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true, plot_difference = true, plot_significance = true)

#######################################################################
# SHOW T-VALUES
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Pz), plot_tvalues = true, plot_critical_t = true)

#######################################################################
# SHIFT DIFFERENCE WAVE
#######################################################################

# offset the difference wave for visibility (similar to MATLAB)
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true, plot_difference = true, difference_offset = 3.0)

#######################################################################
# SIGNIFICANCE BAR POSITIONING
#######################################################################

# position significance bars at the bottom, at zero, or at a custom y
EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels(:Cz),
    plot_erp = true,
    plot_significance = true,
    significance_position = :bottom,
    significance_color = (:red, 0.5),
)

#######################################################################
# GRID LAYOUT — multiple channels
#######################################################################

EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels([:Cz, :Pz, :Oz, :Fz]),
    layout = :grid,
    plot_erp = true,
    plot_significance = true,
)
