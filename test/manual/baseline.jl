using EegFun

# read raw data
dat = EegFun.read_raw_data("./data/files/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Plot 
EegFun.plot_databrowser(dat) # DC shift visible

# Baseline stuff and replot 
EegFun.baseline!(dat)
EegFun.plot_databrowser(dat) # now zero mean over all samples

# baseline to first sample 
EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 1));
EegFun.plot_databrowser(dat)

# EpochData
# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.2, 1.0)  # -200 to 1000 ms

EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(epochs, (0.0, 0.0)) # baseline to t=0
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(epochs, (0.5, 0.5)) # baseline to t=0.5
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1]))

# Create some ERP data
eprs = EegFun.average_epochs(epochs)

EegFun.plot_erp(eprs, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(eprs, (0.0, 0.0)) # baseline to t=0
EegFun.plot_erp(eprs, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(eprs, (0.5, 0.5)) # baseline to t=0.5
EegFun.plot_erp(eprs, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(eprs, (-0.2, 0)) # baseline to t=-0.2 to 0.0
EegFun.plot_erp(eprs, channel_selection = EegFun.channels([:Fp1]))

