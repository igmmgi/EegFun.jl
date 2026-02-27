# Demo: Plotting ICA Component Topographies
# Shows how to visualise ICA component scalp maps as topographic plots,
# select specific components, and adjust display options.

using EegFun

#######################################################################
# LOAD DATA AND RUN ICA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# highpass at 1 Hz for ICA (critical for good decomposition)
EegFun.highpass_filter!(dat, 1.0)

# average reference
EegFun.rereference!(dat, :avg)

# create additional Bool colum for v. extreme values
EegFun.is_extreme_value!(dat, 200)

# run ICA (excludinng the v. extreme samples)
ica = EegFun.run_ica(dat, sample_selection = EegFun.samples_not(:is_extreme_value_200))

#######################################################################
# BASIC TOPOGRAPHY GRID
#######################################################################

# plot all components in a grid
EegFun.plot_topography(ica, label_plot = false, point_plot = false)

#######################################################################
# SELECT SPECIFIC COMPONENTS
#######################################################################

# plot only the first 10 components
EegFun.plot_topography(ica, component_selection = EegFun.components(1:10))

# plot even-numbered components
EegFun.plot_topography(ica, component_selection = EegFun.components(2:2:20))

#######################################################################
# CUSTOM GRID SIZE
#######################################################################

# control rows × columns layout
EegFun.plot_topography(ica, component_selection = EegFun.components(1:12), dims = (3, 4))

#######################################################################
# DISPLAY OPTIONS
#######################################################################

# shared colour scale across all components
EegFun.plot_topography(ica, component_selection = EegFun.components(1:10), use_global_scale = true)

# adjust colour map
EegFun.plot_topography(ica, component_selection = EegFun.components(1:6), colormap = :RdBu)
