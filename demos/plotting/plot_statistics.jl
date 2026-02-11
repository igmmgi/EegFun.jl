# Demo: Plotting Statistical Test Results
# Shows how to visualise analytic and permutation test results with
# ERP waveforms, difference waves, t-values, and significance markers.

using EegFun

#######################################################################
# 1. LOAD DATA AND CREATE ERPS
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
# 2. RUN STATISTICAL TEST
#######################################################################

# prepare for comparison and run analytic test
prepared = EegFun.prepare_statistics(epochs)
result = EegFun.analytic_test(prepared, correction_method = :bonferroni)

#######################################################################
# 3. BASIC PLOT — ERP WAVEFORMS
#######################################################################

# plot ERP waveforms for both conditions at a channel
EegFun.plot_analytic_test(result,
    channel = :Cz,
    plot_erp = true
)

#######################################################################
# 4. ADD DIFFERENCE WAVE
#######################################################################

EegFun.plot_analytic_test(result,
    channel = :Cz,
    plot_erp = true,
    plot_difference = true
)

#######################################################################
# 5. ADD SIGNIFICANCE MARKERS
#######################################################################

EegFun.plot_analytic_test(result,
    channel = :Cz,
    plot_erp = true,
    plot_difference = true,
    show_significance = true
)

#######################################################################
# 6. SHOW T-VALUES
#######################################################################

EegFun.plot_analytic_test(result,
    channel = :Pz,
    plot_tvalues = true,
    show_critical_t = true
)

#######################################################################
# 7. SHIFT DIFFERENCE WAVE
#######################################################################

# offset the difference wave for visibility (similar to MATLAB)
EegFun.plot_analytic_test(result,
    channel = :Cz,
    plot_erp = true,
    plot_difference = true,
    shift_difference = true
)

#######################################################################
# 8. SIGNIFICANCE BAR POSITIONING
#######################################################################

# position significance bars at the bottom, at zero, or at a custom y
EegFun.plot_analytic_test(result,
    channel = :Cz,
    plot_erp = true,
    show_significance = true,
    sig_bar_position = :bottom,
    sig_bar_color = (:red, 0.5)
)
