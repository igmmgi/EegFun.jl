using eegfun
using GLMakie
# using CairoMakie
# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)
# EPOCHS -> ERPs
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]), 
             eegfun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
             eegfun.EpochCondition(name = "ExampleEpoch3", trigger_sequences = [[3]]), 
             eegfun.EpochCondition(name = "ExampleEpoch4", trigger_sequences = [[4]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)
erps = eegfun.average_epochs(epochs)
# ERP Plots (:single)
fig, ax = eegfun.plot_erp(erps, average_channels = false, colormap = :viridis, legend_nbanks = 12)
fig, ax = eegfun.plot_erp(erps[1], average_channels = false, colormap = :viridis, legend_nbanks = 12)

fig, ax = eegfun.plot_erp(erps, average_channels = true)
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8]), average_channels = true)
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8]), average_channels = false)
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8]), average_channels = false, legend_nbanks = 3)
fig, ax = eegfun.plot_erp([erps[1], erps[1]], channel_selection = eegfun.channels([:PO8]))

# ERP Plots (:grid)
fig, ax = eegfun.plot_erp(erps, layout = :grid)
fig, ax = eegfun.plot_erp(erps, channel_selection = eegfun.channels([:Cz, :PO7, :PO8, :Fp1, :Fp2, :F3]), layout = :grid)

# EPR Plots (:topo)
fig, ax = eegfun.plot_erp(erps, layout = :topo)
fig, ax = eegfun.plot_erp([erps[1], erps[1]], layout = :topo)
fig, ax = eegfun.plot_erp(erps, layout = :topo, channel_selection = eegfun.channels([:Fp1, :Fp2, :PO8]))

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
)
fig

GLMakie.closeall()
