using EegFun
using JLD2

# read raw data
dat = EegFun.read_raw_data("./data/raw_files/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.5, 1.0)

# Basic plots for Epochs
EegFun.plot_epochs(epochs[1])
EegFun.plot_epochs(epochs[1], layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels([:Fp1]))

EegFun.plot_epochs(epochs[1], layout = :grid)
EegFun.plot_epochs(epochs[1], layout = :grid, channel_selection = EegFun.channels([:Fp1, :Fp2]))

EegFun.plot_epochs(epochs[1], layout = :topo)
EegFun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = false, theme_fontsize = 10, layout_topo_show_scale = true)

EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]); layout = :topo)

