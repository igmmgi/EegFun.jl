# Demo: Topographic Maps
# Shows scalp topography visualization for ERP data and components.
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and prepare layout file
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# visually selected blink like artifact via interval selection times
EegFun.plot_topography(dat, interval_selection = EegFun.times(5.973, 6.02), ylim = (-200, 200))
EegFun.plot_topography(dat, interval_selection = EegFun.times(6), ylim = (-200, 200))

# blink like artifact via sample selection predicate
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 5.973 .&& x.time .<= 6.02, ylim = (-200, 200))

# blink like artifact across multiple time points
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 5.5 .&& x.time .<= 6.5, ylim = (-75, 75), n_topo = 5, dims = (1, 5))
EegFun.plot_topography(
    dat,
    sample_selection = x -> x.time .>= 5.5 .&& x.time .<= 6.5,
    ylim = (-75, 75),
    n_topo = 5,
    dims = (1, 5),
    colorbar_plot_numbers = 5,
)

# Various combinations
EegFun.plot_topography(dat, colorbar_plot = false, head_radius = 1.25)
EegFun.plot_topography(dat, gridscale = 250)
EegFun.plot_topography(dat, colormap = :inferno)
EegFun.plot_topography(dat, colormap = :redblue)
EegFun.plot_topography(dat, title = "Custom Title", title_fontsize = 30)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, method = :spherical_spline)
EegFun.plot_topography(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
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

# Separate plots
EegFun.plot_topography(epochs[1], 1) # Condition 1, Epoch 1
EegFun.plot_topography(epochs[2], 1) # Condition 2, Epoch 1
EegFun.plot_topography(epochs, ylim = (-0.1, 0.1)) # TODO: aspect ration?; global scale?
EegFun.plot_topography(epochs, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)



#################################
# ERP like data
#################################
erps = EegFun.average_epochs(epochs)

EegFun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2), time_unit = :ms)
EegFun.plot_topography(erps, interval_selection = (0.4, 0.6), ylim = (-2, 2))

EegFun.plot_topography(erps, interval_selection = (0.4, 0.6), ylim = (-2, 2), colorbar_plot_numbers = [1, 2])

EegFun.plot_topography(
    erps,
    interval_selection = (0.4, 0.6),
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

