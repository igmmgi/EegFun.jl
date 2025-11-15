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

epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)
erps = eegfun.average_epochs(epochs)

# ERP Plot


fig, ax = eegfun.plot_erp(erps, average_channels = false)
fig, ax = eegfun.plot_erp(erps, average_channels = true)
fig, ax = eegfun.plot_erp(erps, average_channels = false, layout = :grid)

# ERP Plot
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8]), average_channels = false)
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8]), average_channels = true)


# ERP Plot
fig, ax = eegfun.plot_erp(erps[1], channel_selection = eegfun.channels([:Cz, :PO7, :PO8, :C1]), layout = [2, 2])
fig, ax = eegfun.plot_erp([erps[1], copy(erps[1])])
fig, ax = eegfun.plot_erp(erps[1], layout = :grid)
fig, ax = eegfun.plot_erp(erps, layout = :grid, title = "Custom Title")
fig, ax = eegfun.plot_erp(erps[1], layout = :topo)
fig, ax = eegfun.plot_erp(erps, layout = :topo, channel_selection = eegfun.channels([:Fp1, :Fp2, :PO8]))
fig, ax = eegfun.plot_erp(erps, layout = :grid, channel_selection = x -> startswith.(string.(x), "F"))
fig, ax = eegfun.plot_erp([erps[1], copy(erps[1])], channel_selection = eegfun.channels([:Fp1, :Fp2]), layout = :grid)
fig, ax = eegfun.plot_erp(erps[1], linewidth = 2)
fig, ax = eegfun.plot_erp(erps[1], layout = :grid)
fig, ax = eegfun.plot_erp(erps[1], layout = :topo)

# Combined plots
fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 1], width = Relative(0.2), height = Relative(0.2), halign = 0, valign = 0)
eegfun.plot_erp!(fig, ax1, erps, average_channels = true)
eegfun.plot_topography!(
    fig,
    ax2,
    erps[1];
    point_plot = false,
    label_plot = false,
    colorbar_plot = true,
    colorbar_width = Relative(0.03),
    colorbar_height = Relative(0.2),
    colorbar_tellheight = false,
    colorbar_tellwidth = false,
    colorbar_position = (1, 1),
    colorbar_halign = 0.25,
    colorbar_valign = 0,
    colorbar_flipaxis = true,
    radius = 1.5,
)
# eegfun.plot_topography!(fig, ax2, erps[1]; point_plot = false, label_plot = false)
fig


GLMakie.closeall()
