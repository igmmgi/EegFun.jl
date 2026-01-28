using EegFun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_17.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_bdf(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

EegFun.plot_erp_image(epochs[1], layout = :single)
EegFun.plot_erp_image(epochs[1], layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))

EegFun.plot_erp_image(epochs[1], layout = :grid, colorbar_plot = false)
EegFun.plot_erp_image(epochs[1], layout = :grid, colorbar_plot = true)

EegFun.plot_erp_image(epochs[1], layout = :topo)
EegFun.plot_erp_image(epochs[1], channel_selection = EegFun.channels([:PO7]))
EegFun.plot_erp_image(epochs[1], channel_selection = EegFun.channels([:Fp1]), plot_erp = false)


EegFun.plot_erp_image(epochs[1], layout = :single)

fig, axes = EegFun.plot_erp_image(
    epochs[1],
    # channel_selection = EegFun.channels([:Fp1, :Fp2]),
    # channel_selection = EegFun.channels([:Fp1, :Fp2]),
    layout = :topo,
    boxcar_average = 20,
    colorrange = (-50, 50),
)


# TODO: no electrode labels in title
EegFun.plot_erp_image(epochs[1], layout = :topo)
EegFun.plot_erp_image(epochs[1], layout = :grid)
