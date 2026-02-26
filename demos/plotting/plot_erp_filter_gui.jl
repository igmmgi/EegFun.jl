# Demo: ERP Filter GUI
# Shows how to launch the interactive filter exploration tool
# to compare the effect of different filter settings on ERP waveforms.

using EegFun

#######################################################################
# LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# Create ERP
erp = EegFun.average_epochs(epochs)

#######################################################################
# LAUNCH FILTER GUI — SINGLE CONDITION
#######################################################################

# interactive GUI with sliders for cutoff, order, method
EegFun.plot_erp_filter_gui(erp)

#######################################################################
# SPECIFY A CHANNEL
#######################################################################

# focus on a specific channel
EegFun.plot_erp_filter_gui(erp, channel = :Cz)

#######################################################################
# COMPARE MULTIPLE CONDITIONS
#######################################################################

# pass a vector of ERPs to see filter effects side-by-side
epochs2 = EegFun.epoch_data(dat, [:trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs2, (-0.2, 0.0))
erp2 = EegFun.average(epochs2)

EegFun.plot_erp_filter_gui([erp, erp2], channel = :Cz)
