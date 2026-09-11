# Demo: Epoch Plotting
# Shows epoch visualization with channel selection and overlay options.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# Basic plots for Epochs
EegFun.plot_epochs(epochs) # not that useful as crowded!
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:PO7, :PO8]))

EegFun.plot_epochs(epochs, layout = :grid)
EegFun.plot_epochs(epochs, layout = :grid, channel_selection = EegFun.channels([:Fp1, :Fp2, :PO7, :PO8]))

EegFun.plot_epochs(epochs, layout = :topo)
EegFun.plot_epochs(epochs, layout = :topo, add_xy_origin = false, theme_fontsize = 10, layout_topo_show_scale = true)

EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]); layout = :topo)






