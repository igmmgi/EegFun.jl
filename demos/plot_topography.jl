using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
# EegFun.polar_to_cartesian_xy!(layout_file, preserve_radial_distance = true)
EegFun.polar_to_cartesian_xy!(layout_file, preserve_radial_distance = true)


# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

EegFun.plot_topography(dat)
EegFun.plot_topography(dat, method = :nearest)

EegFun.plot_topography(
    dat,
    sample_selection = x -> x.time .>= 5.973 .&& x.time .<= 6.02,
    gridscale = 75,
    ylim = (-200, 200),
    head_radius = 1.0,
    method = :thin_plate,
)

EegFun.plot_topography(dat, interval_selection = (5.973, 6.02), gridscale = 100)
EegFun.plot_topography(dat, interval_selection = (5, 6), gridscale = 100)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 5.973 .&& x.time .<= 6.02, gridscale = 100, ylim = (-100, 100))

# Various combinations
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
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

EegFun.plot_topography(epochs, 1) # epoch 1
EegFun.plot_topography(epochs, 2) # epoch 2

# EegFun.plot_topography(epochs) # TODO: would this be useful?

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
    colorbar_plot_numbers = [2],
    dims = (1, 2),
)

EegFun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2), colorbar_plot = false)

EegFun.plot_topography(
    erps,
    sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6,
    ylim = (-2, 2),
    colorbar_plot = true,
    colorbar_position = (2, 1),
    colorbar_vertical = false,
)

EegFun.plot_topography(erps)
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
