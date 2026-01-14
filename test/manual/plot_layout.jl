using eegfun
using GLMakie
layout = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout);

eegfun.plot_layout_2d(layout; head_color = :blue)
eegfun.plot_layout_2d(layout; head_linewidth = 5)
eegfun.plot_layout_2d(layout; head_radius = 1.8)

eegfun.plot_layout_2d(layout; point_plot = false)
eegfun.plot_layout_2d(layout; point_marker = :x)
eegfun.plot_layout_2d(layout; point_markersize = 15)
eegfun.plot_layout_2d(layout; point_color = :red)

eegfun.plot_layout_2d(layout; label_plot = false)
eegfun.plot_layout_2d(layout; label_fontsize = 30)
eegfun.plot_layout_2d(layout; label_color = :green)
eegfun.plot_layout_2d(layout; label_xoffset = 0.1)
eegfun.plot_layout_2d(layout; label_yoffset = 0.1)
eegfun.plot_layout_2d(layout; label_zoffset = 2)


fig, ax = eegfun.plot_layout_2d(layout)
eegfun.add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], roi_border_size = 0.05)
eegfun.add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], roi_border_size = 0.1)
eegfun.add_topo_rois!(
    ax,
    layout,
    [[:Fp1]],
    roi_border_size = 0.10,
    roi_fill = [true],
    roi_linecolor = [:black],
    roi_fillcolor = [:red],
    roi_fillalpha = [0.2],
)
eegfun.add_topo_rois!(
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

# Plots with neighbours
eegfun.get_layout_neighbours_xy!(layout, 0.5);
eegfun.plot_layout_2d(layout, neighbours = true)

# Plots with neighbours
eegfun.get_layout_neighbours_xyz!(layout, 0.5);
eegfun.plot_layout_3d(layout, neighbours = true)

# how to print neighbours to a file
eegfun.print_layout_neighbours(layout, "electrode_neighbours_1.toml")
eegfun.print_layout_neighbours(layout.neighbours, "electrode_neighbours_2.toml")

GLMakie.closeall()

# save a basic figure
# NB. for vector graphics, use CairoMakie
# save("topo_roi.png", fig)
