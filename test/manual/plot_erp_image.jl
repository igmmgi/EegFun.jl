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

# EPOCHS
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = []
for (idx, epoch) in enumerate(epoch_cfg)
    push!(epochs, eegfun.extract_epochs(dat, idx, epoch, -2, 4))
end

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
