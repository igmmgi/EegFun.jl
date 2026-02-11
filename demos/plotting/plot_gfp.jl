# Demo: Plotting Global Field Power
# Shows how to visualise GFP (Global Field Power) and Global Dissimilarity
# across conditions, with optional ERP trace overlay.

using EegFun

#######################################################################
# 1. LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# epoch and average
epochs = EegFun.epoch_data(dat, [:trigger1, :trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))
erp1 = EegFun.average(epochs, condition_selection = conditions(:trigger1))
erp2 = EegFun.average(epochs, condition_selection = conditions(:trigger2))

#######################################################################
# 2. BASIC GFP PLOT
#######################################################################

# plot GFP for a single condition
EegFun.plot_gfp(erp1)

#######################################################################
# 3. COMPARE CONDITIONS
#######################################################################

# overlay GFP for multiple conditions
EegFun.plot_gfp([erp1, erp2])

#######################################################################
# 4. WITH ERP TRACES AND DISSIMILARITY
#######################################################################

# show all three panels: ERP traces, GFP, and Global Dissimilarity
EegFun.plot_gfp([erp1, erp2],
    show_erp_traces = true,
    show_dissimilarity = true
)

#######################################################################
# 5. RAW VS NORMALISED
#######################################################################

# plot in raw microvolts instead of percentage
EegFun.plot_gfp(erp1, normalize = false)

#######################################################################
# 6. PRE-COMPUTED GFP
#######################################################################

# compute GFP separately, then plot the result
gfp_result = EegFun.gfp(erp1, normalize = true)
EegFun.plot_gfp(gfp_result, color = :red, linewidth = 3)

#######################################################################
# 7. CHANNEL SELECTION
#######################################################################

# compute GFP over a subset of channels
EegFun.plot_gfp(erp1, channel_selection = channels([:Cz, :Pz, :Oz]))
