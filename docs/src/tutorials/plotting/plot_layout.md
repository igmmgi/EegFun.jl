# Plot Layout

This demo demonstrates visualizing electrode layout configurations in 2D and 3D space.

## What is Electrode Layout?

Electrode layout defines the spatial positions of EEG sensors on the scalp using:

- **2D coordinates** (x, y): Projected positions for topographic maps
- **3D coordinates** (x, y, z): Actual positions on spherical head model
- **Spherical coordinates** (incidence, azimuth): Angular positions

## Layout Visualization Functions

**plot_layout_2d**:

- Standard topographic view from above
- Shows electrode positions and labels
- Customizable head outline, markers, and labels
- Supports region-of-interest highlighting

**plot_layout_3d**:

- Full 3D spherical visualization
- Interactive rotation and zoom
- Shows true spatial arrangement

**Neighbor visualization**:

- Display spatial relationships between electrodes
- Helps verify interpolation neighborhoods
- Interactive hover to show connections

## Use Cases

**Verify electrode montage**:

- Confirm correct channel positions
- Check for mislabeled electrodes
- Validate custom layouts

**Publication figures**:

- Document electrode positions
- Highlight regions of interest
- Create publication-ready diagrams

**Analysis planning**:

- Understand spatial sampling density
- Plan interpolation strategies
- Define ROIs for analysis

## Customization Options

**Head outline**:

- Color, line width, radius

**Electrode markers**:

- Marker style, size, color
- Show/hide markers

**Labels**:

- Font size, color
- Position offsets
- Show/hide labels

**Regions of interest (ROIs)**:

- Highlight electrode groups
- Customizable borders and fills
- Multiple ROIs with different styles

## Neighbor Detection

Calculate spatial neighbors based on distance thresholds:

- **2D neighbors**: Using projected coordinates
- **3D neighbors**: Using spherical coordinates
- Export neighbor definitions for channel repair

## Workflow Summary

This demo shows electrode layout visualization:

## Basic 2D Visualization

- Load electrode layout
- Convert to 2D coordinates
- Plot with default settings

## Customize Appearance

- Adjust head outline (color, width, radius)
- Modify markers (style, size, color)
- Customize labels (size, color, offsets)

## Add Regions of Interest

- Highlight electrode groups
- Customize ROI borders and fills
- Create multiple ROIs with different styles

## Neighbor Visualization

- Calculate neighbor relationships
- Plot with neighbor connections
- Export neighbor definitions

## Save Figures

- Export publication-ready figures
- Use CairoMakie for vector graphics


## Code Examples

::: details Show Code

```julia
# Demo: Layout Visualization
# Shows electrode layout plotting and manipulation.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
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
result = EegFun.plot_layout_2d(layout)
EegFun.add_topo_rois!(result.axes[1], layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], roi_border_size = 0.05)
EegFun.add_topo_rois!(result.axes[1], layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], roi_border_size = 0.1)
EegFun.add_topo_rois!(
    result.axes[1],
    layout,
    [[:Fp1]],
    roi_border_size = 0.10,
    roi_fill = [true],
    roi_linecolor = [:black],
    roi_fillcolor = [:red],
    roi_fillalpha = [0.2],
)
EegFun.add_topo_rois!(
    result.axes[1],
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

# Print neighbours while 2D neighbours are still computed
# (the 3D coordinate conversion below clears cached neighbours)
EegFun.print_layout_neighbours(layout, joinpath(DEMO_OUTPUT, "electrode_neighbours_1.toml"))
EegFun.print_layout_neighbours(layout.neighbours, joinpath(DEMO_OUTPUT, "electrode_neighbours_2.toml"))

# 3D layout with neighbours
EegFun.get_neighbours_xyz!(layout, 0.5);
EegFun.plot_layout_3d(layout, neighbours = true)

# save a basic figure
# NB. for vector graphics, use CairoMakie
# save("topo_roi.png", result.fig)
```

:::

## See Also

- [API Reference](../../reference/index.md)
