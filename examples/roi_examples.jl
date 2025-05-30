using EEGfun
using GLMakie

# Read and prepare layout
layout = read_layout("./layouts/biosemi64.csv")
polar_to_cartesian_xy!(layout)

# Create figure and plot basic layout
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1, 1], title="ROI Examples")
plot_layout_2d!(fig, ax, layout)

# Add different types of ROIs

# 1. Simple outline ROI using LibGEOS
add_topo_rois!(ax, layout,
    [[:Fp1, :Fp2, :AF3]],
    border_size = 10,
    roi_kwargs = Dict(
        :color => [:black],
        :linewidth => [2]
    )
)

# 2. Filled ROIs using Graham's Scan
add_topo_rois!(ax, layout,
    [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]],
    border_size = 10,
    use_libgeos = false,
    roi_kwargs = Dict(
        :fill => [true, true],
        :fillcolor => [:red, :blue],
        :fillalpha => [0.2, 0.2],
        :color => [:darkred, :darkblue],
        :linewidth => [2, 2]
    )
)

# 3. Central ROI with different style
add_topo_rois!(ax, layout,
    [[:Cz, :C1, :C2, :FCz]],
    border_size = 8,
    roi_kwargs = Dict(
        :fill => [true],
        :fillcolor => [:green],
        :fillalpha => [0.15],
        :color => [:darkgreen],
        :linewidth => [3]
    )
)

# Add a title to explain the ROIs
Label(fig[0, 1],
    """ROI Examples:
    Black: Simple outline (LibGEOS)
    Red/Blue: Filled regions (Graham's Scan)
    Green: Central region""",
    tellheight = true)

# Display the figure
display(fig) 