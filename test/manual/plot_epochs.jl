using EegFun
using GLMakie
# using CairoMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_6.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)

EegFun.plot_layout_2d(layout_file)

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)
EegFun.is_extreme_value!(dat, 500);
EegFun.mark_epoch_windows!(dat, [1, 2, 3, 4], [-0.5, 2.0]) # simple epoch marking with trigger 1 and 3

# EPOCHS
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1], [2]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[3], [4]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.5, 2)

EegFun.plot_epochs(epochs[1], layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels([:Fp1]))
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))


EegFun.plot_epochs(epochs[1], layout = :single)
EegFun.plot_epochs(epochs, layout = :single)
EegFun.plot_epochs(epochs[1], layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))


EegFun.plot_epochs(epochs[1], layout = :grid)
EegFun.plot_epochs(epochs, layout = :grid)
EegFun.plot_epochs(epochs[1], layout = :topo)
EegFun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = false)
EegFun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = true, theme_fontsize = 10, layout_topo_show_scale = true)
EegFun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = false, theme_fontsize = 30)

EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]); layout = :topo)


