using EegFun

layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout, preserve_radial_distance = true);

# basic plot
EegFun.plot_layout_2d(layout)

# some customisations
EegFun.plot_layout_2d(layout; head_color = :grey, head_linewidth = 5, head_radius = 1.5)
EegFun.plot_layout_2d(layout; point_plot = false, label_plot = false)
EegFun.plot_layout_2d(layout; point_marker = :x, point_markersize = 20)
EegFun.plot_layout_2d(layout; point_color = :red)

EegFun.plot_layout_2d(layout; label_fontsize = 20, label_color = :grey)

# Annotating plots with ROIs
fig, ax = EegFun.plot_layout_2d(layout)
EegFun.add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], roi_border_size = 0.05)
EegFun.add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], roi_border_size = 0.1)
EegFun.add_topo_rois!(
    ax,
    layout,
    [[:Fp1]],
    roi_border_size = 0.10,
    roi_fill = [true],
    roi_linecolor = [:black],
    roi_fillcolor = [:red],
    roi_fillalpha = [0.2],
)
EegFun.add_topo_rois!(
    ax,
    layout,
    [[:CPz, :C2, :FCz, :C1]],
    roi_border_size = 0.15,
    roi_linewidth = [5],
    roi_fill = [true],
    roi_linecolor = [:blue],
    roi_fillcolor = [:green],
    roi_fillalpha = [0.5],
)

# Plots with neighbours (hover mouse over channel to show neighbour connections)
EegFun.get_neighbours_xy!(layout, 0.5);
EegFun.plot_layout_2d(layout, neighbours = true)

EegFun.get_neighbours_xyz!(layout, 0.5);
EegFun.plot_layout_3d(layout, neighbours = true)

# how to print neighbours to a file
EegFun.print_layout_neighbours(layout, "electrode_neighbours_1.toml")
EegFun.print_layout_neighbours(layout.neighbours, "electrode_neighbours_2.toml")

# save a basic figure
# NB. for vector graphics, use CairoMakie
# save("topo_roi.png", fig)
