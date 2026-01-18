using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_17.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

# EPOCHS
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

eegfun.plot_erp_image(epochs[1], layout = :single)
eegfun.plot_erp_image(epochs[1], layout = :single, channel_selection = eegfun.channels([:Fp1, :Fp2]))

eegfun.plot_erp_image(epochs[1], layout = :grid, colorbar_plot = false)
eegfun.plot_erp_image(epochs[1], layout = :grid, colorbar_plot = true)

eegfun.plot_erp_image(epochs[1], layout = :topo)
eegfun.plot_erp_image(epochs[1], channel_selection = eegfun.channels([:Fp1]))
eegfun.plot_erp_image(epochs[1], channel_selection = eegfun.channels([:Fp1]), plot_erp = false)


eegfun.plot_erp_image(epochs[1], layout = :single)

fig, axes = eegfun.plot_erp_image(
    epochs[1],
    channel_selection = eegfun.channels([:Fp1, :Fp2]),
    layout = :single,
    boxcar_average = 20,
    colorrange = (-20, 20),
)


# TODO: no electrode labels in title
eegfun.plot_erp_image(epochs[1], layout = :topo)
eegfun.plot_erp_image(epochs[1], layout = :grid)
