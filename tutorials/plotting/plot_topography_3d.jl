# Demo: 3D Topographic Maps
# Shows interactive 3D scalp topography visualization for ERP data and components.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# Important: To interact with 3D plots (rotate, pan, zoom), you should use GLMakie!
using GLMakie

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
# Note: plot_topography_3d uses the spherical coordinates (inc, azi) from the layout,
# rather than the 2D projected (x2, y2) coordinates. 

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# visually selected blink like artifact via interval selection times
EegFun.plot_topography_3d(dat, interval_selection = EegFun.times(5.973, 6.02), ylim = (-200, 200))
EegFun.plot_topography_3d(dat, interval_selection = EegFun.times(6), ylim = (-200, 200))

# blink like artifact via sample selection predicate
EegFun.plot_topography_3d(dat, sample_selection = x -> x.time .>= 5.973 .&& x.time .<= 6.02, ylim = (-200, 200))

# Various stylistic combinations (sharing aesthetics with plot_topography)
EegFun.plot_topography_3d(dat, colorbar_plot = false)
EegFun.plot_topography_3d(dat, colormap = :inferno)
EegFun.plot_topography_3d(dat, colormap = :redblue)
EegFun.plot_topography_3d(dat, title = "3D Custom Title", title_fontsize = 30)
EegFun.plot_topography_3d(dat, label_plot = false, point_plot = false)

# You can also change the default starting camera angle (in degrees). 
# For example, to look at the top of the head from the side:
EegFun.plot_topography_3d(dat, camera_azimuth = 90, camera_elevation = 45)

# Customizing the 3D markers and text
EegFun.plot_topography_3d(
    dat,
    sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6,
    point_plot = true,
    point_markersize = 10,
    point_color = :white,
    label_plot = true,
    label_fontsize = 14,
    label_color = :black,
)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Separate 3D plots
EegFun.plot_topography_3d(epochs[1], 1, interval_selection = (0.4, 0.6)) # Condition 1, Epoch 1
EegFun.plot_topography_3d(epochs[2], 1, interval_selection = (0.4, 0.6)) # Condition 2, Epoch 1

#################################
# ERP like data
#################################
erps = EegFun.average_epochs(epochs)

EegFun.plot_topography_3d(erps[1], sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2))
EegFun.plot_topography_3d(erps[1], interval_selection = (0.4, 0.6), ylim = (-2, 2))

EegFun.plot_topography_3d(
    erps[2],
    interval_selection = (0.4, 0.6),
    ylim = (-2, 2),
    colorbar_plot = true,
    colorbar_position = (1, 2), # Note: GLMakie LScene placement
)
