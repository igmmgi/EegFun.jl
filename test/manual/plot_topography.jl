using eegfun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

#################################
# Single DataFrameEeg
#################################
eegfun.plot_topography(dat)
eegfun.plot_topography(dat, ylim = (-0.01, 0.01))
eegfun.plot_topography(dat, colorbar_plot = false)

eegfun.plot_topography(dat, method = :spherical_spline)
eegfun.plot_topography(dat, method = :multiquadratic)
eegfun.plot_topography(dat, method = :multiquadratic, head_radius = 0.9)
eegfun.plot_topography(dat, gridscale = 50)
eegfun.plot_topography(dat, gridscale = 1000)
eegfun.plot_topography(dat, colormap = :inferno)
eegfun.plot_topography(dat, title = "Custom Title", title_fontsize = 30)
eegfun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)
eegfun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, method = :spherical_spline)
eegfun.plot_topography(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :Cz]))
eegfun.plot_topography(dat, label_fontsize = 130)
eegfun.plot_topography(dat, point_markersize = 30, point_marker = :x)
eegfun.plot_topography(dat, colorbar_ticksize = 130)
eegfun.plot_topography(dat, colorbar_labelcolor = :red)
eegfun.plot_topography(dat, colorbar_size = 50)
eegfun.plot_topography(dat, colorbar_size = 20, colorbar_position = (2, 1), colorbar_vertical=false)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    eegfun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
    eegfun.EpochCondition(name = "ExampleEpoch3", trigger_sequences = [[3]]),
    eegfun.EpochCondition(name = "ExampleEpoch4", trigger_sequences = [[4]]),
]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

eegfun.plot_topography(epochs, 2)
eegfun.plot_topography(epochs[1], 2)



eegfun.plot_topography(epochs, 1)
eegfun.plot_topography(epochs[2], 1)
eegfun.plot_topography(epochs[1], 1, gridscale = 50)
eegfun.plot_topography(epochs[2], 1, gridscale = 1000)
eegfun.plot_topography(epochs[1], 1, colormap = :inferno)
eegfun.plot_topography(epochs[2], 1, title = "Custom Title", title_fontsize = 30)
eegfun.plot_topography(epochs[1], 1, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)
eegfun.plot_topography(epochs[2], 1, channel_selection = eegfun.channels([:Fp1, :Fp2, :Cz]))
eegfun.plot_topography(epochs[1], 1, label_fontsize = 30)
eegfun.plot_topography(epochs[2], 1, point_markersize = 30, point_marker = :x)
eegfun.plot_topography(epochs[1], 1, ylim = (-10, 10))
eegfun.plot_topography(epochs[2], 1, ylim = (-1, 1))


#################################
# ERP like data
#################################
erps = eegfun.average_epochs(epochs)

eegfun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2))

eegfun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2), 
colorbar_datasets = [4], dims = (1, 4))

eegfun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2), 
colorbar_plot = false)

eegfun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2), 
colorbar_plot = true, colorbar_position = (2, 1), colorbar_vertical = false)


eegfun.plot_topography(erps[1])

eegfun.plot_topography(erps[1:3])
eegfun.plot_topography(erps[1:2])


eegfun.plot_topography(erps, ylim = (-2, 2))
eegfun.plot_topography(erps[2])

eegfun.plot_topography(erps[1], gridscale = 50)
eegfun.plot_topography(erps[2], gridscale = 1000)
eegfun.plot_topography(erps[1], colormap = :inferno)
eegfun.plot_topography(erps[2], title = "Custom Title", title_fontsize = 30)



##################################################################
# Combine Plot Examples

fig = Figure(size = (700, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
eegfun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 2))
eegfun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_position = (2, 2))
fig

fig = Figure(size = (1400, 400))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
eegfun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 3))
eegfun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_plot = false)
fig

fig = Figure(size = (1400, 400))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 3])
eegfun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 2))
eegfun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_position = (1, 4))
fig

fig = Figure(size = (1000, 1000))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 3])
ax3 = Axis(fig[2, 1:4])
eegfun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 2))
eegfun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_position = (1, 4))
eegfun.plot_topography!(fig, ax3, epochs[2], 1, display_plot = false, colorbar_plot = false)
fig

GLMakie.closeall()
