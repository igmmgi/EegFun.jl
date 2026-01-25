using EegFun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file, normalization_radius = 1.0)
dat = EegFun.read_bdf(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);
EegFun.rereference!(dat, :avg)
EegFun.filter_data!(dat, "hp", 1)


EegFun.plot_topography(
    dat,
    method = :spherical_spline,
    sample_selection = x -> x.time .>= 7.984 .&& x.time .<= 8.168,
    gridscale = 100,
)

EegFun.plot_topography(
    dat,
    sample_selection = x -> x.time .>= 7.984 .&& x.time .<= 8.168,
    ylim = (-30, 30),
    gridscale = 100,
)

EegFun.plot_databrowser(dat)


#################################
# Single DataFrameEeg
#################################
EegFun.plot_topography(dat)
EegFun.plot_topography(dat, ylim = (-0.01, 0.01))
EegFun.plot_topography(dat, colorbar_plot = false)

EegFun.plot_topography(dat, method = :nearest)
EegFun.plot_topography(dat, method = :shepard)
EegFun.plot_topography(dat, method = :spherical_spline)
EegFun.plot_topography(dat, method = :thin_plate)
EegFun.plot_topography(dat, method = :multiquadratic, head_radius = 0.9)
EegFun.plot_topography(dat, gridscale = 50)
EegFun.plot_topography(dat, gridscale = 1000)
EegFun.plot_topography(dat, colormap = :inferno)
EegFun.plot_topography(dat, title = "Custom Title", title_fontsize = 30)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, method = :spherical_spline)
EegFun.plot_topography(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_topography(dat, label_fontsize = 130)
EegFun.plot_topography(dat, point_markersize = 30, point_marker = :x)
EegFun.plot_topography(dat, colorbar_ticksize = 130)
EegFun.plot_topography(dat, colorbar_labelcolor = :red)
EegFun.plot_topography(dat, colorbar_size = 50)
EegFun.plot_topography(dat, colorbar_size = 20, colorbar_position = (2, 1), colorbar_vertical = false)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
    EegFun.EpochCondition(name = "ExampleEpoch3", trigger_sequences = [[3]]),
    EegFun.EpochCondition(name = "ExampleEpoch4", trigger_sequences = [[4]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

EegFun.plot_topography(epochs, 2)
EegFun.plot_topography(epochs[1], 2)



EegFun.plot_topography(epochs, 1)
EegFun.plot_topography(epochs[2], 1)
EegFun.plot_topography(epochs[1], 1, gridscale = 50)
EegFun.plot_topography(epochs[2], 1, gridscale = 1000)
EegFun.plot_topography(epochs[1], 1, colormap = :inferno)
EegFun.plot_topography(epochs[2], 1, title = "Custom Title", title_fontsize = 30)
EegFun.plot_topography(epochs[1], 1, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)
EegFun.plot_topography(epochs[2], 1, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_topography(epochs[1], 1, label_fontsize = 30)
EegFun.plot_topography(epochs[2], 1, point_markersize = 30, point_marker = :x)
EegFun.plot_topography(epochs[1], 1, ylim = (-10, 10))
EegFun.plot_topography(epochs[2], 1, ylim = (-1, 1))


#################################
# ERP like data
#################################
erps = EegFun.average_epochs(epochs)

EegFun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2))

EegFun.plot_topography(
    erps,
    sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6,
    ylim = (-2, 2),
    colorbar_plot_numbers = [4],
    dims = (1, 4),
)

EegFun.plot_topography(
    erps,
    sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6,
    ylim = (-2, 2),
    colorbar_plot = false,
)

EegFun.plot_topography(
    erps,
    sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6,
    ylim = (-2, 2),
    colorbar_plot = true,
    colorbar_position = (2, 1),
    colorbar_vertical = false,
)


EegFun.plot_topography(erps[1])

EegFun.plot_topography(erps[1:3])
EegFun.plot_topography(erps[1:2])


EegFun.plot_topography(erps, ylim = (-2, 2))
EegFun.plot_topography(erps[2])

EegFun.plot_topography(erps[1], gridscale = 50)
EegFun.plot_topography(erps[2], gridscale = 1000)
EegFun.plot_topography(erps[1], colormap = :inferno)
EegFun.plot_topography(erps[2], title = "Custom Title", title_fontsize = 30)



##################################################################
# Combine Plot Examples

fig = Figure(size = (700, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
EegFun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 2))
EegFun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_position = (2, 2))
fig

fig = Figure(size = (1400, 400))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
EegFun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 3))
EegFun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_plot = false)
fig

fig = Figure(size = (1400, 400))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 3])
EegFun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 2))
EegFun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_position = (1, 4))
fig

fig = Figure(size = (1000, 1000))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 3])
ax3 = Axis(fig[2, 1:4])
EegFun.plot_topography!(fig, ax1, epochs[1], 1, display_plot = false, colorbar_position = (1, 2))
EegFun.plot_topography!(fig, ax2, epochs[2], 1, display_plot = false, colorbar_position = (1, 4))
EegFun.plot_topography!(fig, ax3, epochs[2], 1, display_plot = false, colorbar_plot = false)
fig

GLMakie.closeall()
