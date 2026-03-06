# Demo: ERP Image Plots
# Shows trial-by-trial ERP activity as image/heatmap plots.

using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_test_eegfun_data(dat, layout_file)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

EegFun.plot_erp_image(epochs, layout = :single)  # average of all channels
EegFun.plot_erp_image(epochs, layout = :single, channel_selection = EegFun.channels([:PO7, :PO8]))
EegFun.plot_erp_image(epochs, layout = :single, channel_selection = EegFun.channels([:PO7, :PO8]), interpolate = true)

EegFun.plot_erp_image(epochs, layout = :grid, colorbar_plot = false)
EegFun.plot_erp_image(epochs, layout = :grid, colorbar_plot = true)

EegFun.plot_erp_image(epochs, layout = :topo)
EegFun.plot_erp_image(epochs, channel_selection = EegFun.channels([:PO7]))
EegFun.plot_erp_image(epochs, channel_selection = EegFun.channels([:Fp1]), plot_erp = false)


EegFun.plot_erp_image(epochs, layout = :single)

fig, axes = EegFun.plot_erp_image(
    epochs,
    channel_selection = EegFun.channels([:PO7, :PO8]),
    layout = :single,
    boxcar_average = 2,
    colorrange = (-50, 50),
)

fig, axes = EegFun.plot_erp_image(epochs, layout = :topo, boxcar_average = 20, colorrange = (-50, 50))
