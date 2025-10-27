using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)
eegfun.is_extreme_value!(dat, 500);
eegfun.mark_epoch_windows!(dat, [1, 3, 4, 5], [-1, 1.0]) # simple epoch marking with trigger 1 and 3

# EPOCHS
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

# TODO: fix epoch slider 
eegfun.plot_databrowser(epochs[1])

eegfun.plot_epochs(epochs[1], layout = :single)
eegfun.plot_epochs(epochs[1], layout = :single, channel_selection = eegfun.channels([:Fp1, :Fp2]))

eegfun.plot_epochs(epochs[1], layout = :grid)
eegfun.plot_epochs(epochs[1], layout = :topo)
eegfun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = false)
eegfun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = false, theme_fontsize = 10)

# TODO: ERP plotting bring to front
eegfun.plot_epochs(epochs[1], channel_selection = eegfun.channels([:Fp1, :Fp2, :Cz]))
eegfun.plot_epochs(epochs[1], channel_selection = eegfun.channels([:Fp1, :Fp2, :Cz]); layout = :topo)
