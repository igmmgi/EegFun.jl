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
dat = EegFun.create_test_eegfun_data(dat, layout)

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)

# EPOCHS
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# Create ERP
erp = EegFun.average_epochs(epochs)

# interactive GUI with sliders for cutoff, order, method
EegFun.plot_erp_filter_gui(erp)

