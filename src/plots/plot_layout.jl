########################################################
# 2D layout
########################################################
"""
    plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame;
                  head_kwargs::Dict=Dict(), point_kwargs::Dict=Dict(),
                  label_kwargs::Dict=Dict())

Plot a 2D EEG electrode layout with customizable head shape, electrode points, and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: DataFrame containing electrode positions with columns x2, y2, and label
- `head_kwargs`: Keyword arguments for head shape rendering (color, linewidth, etc.)
- `point_kwargs`: Keyword arguments for electrode points (plot_points, marker, size, color etc.)
- `label_kwargs`: Keyword arguments for electrode labels (plot_labels, fontsize, color, xoffset, yoffset, etc.)

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xy!(layout)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout)
"""
function plot_layout_2d!(
    fig::Figure,
    ax::Axis,
    layout::DataFrame;
    head_kwargs::Dict = Dict(),
    point_kwargs::Dict = Dict(),
    label_kwargs::Dict = Dict(),
)

    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    # Default kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    merged_head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    merged_point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = merged_point_kwargs[:plot_points]
    point_plot_kwargs = filter(p -> p.first != :plot_points, merged_point_kwargs)

    label_default_kwargs = Dict(:plot_labels => true, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0)
    merged_label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = merged_label_kwargs[:plot_labels]
    xoffset = merged_label_kwargs[:xoffset]
    yoffset = merged_label_kwargs[:yoffset]
    label_plot_kwargs = filter(p -> p.first ∉ (:plot_labels, :xoffset, :yoffset), merged_label_kwargs)

    # Head shape - Use hardcoded radius
    radius = 88 # mm
    arc!(ax, Point2f(0), radius * 2, -π, π; merged_head_kwargs...) # head
    arc!(ax, Point2f(radius * 2, 0), radius * 2 / 7, -π / 2, π / 2; merged_head_kwargs...) # ear right
    arc!(ax, Point2f(-radius * 2, 0), -radius * 2 / 7, π / 2, -π / 2; merged_head_kwargs...) # ear left
    lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * 4; merged_head_kwargs...) # nose

    # Regular points
    if plot_points
        scatter!(ax, layout[!, :x2], layout[!, :y2]; point_plot_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(ax, position = (label.x2 + xoffset, label.y2 + yoffset), String(label.label); label_plot_kwargs...)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return fig, ax
end


"""
    plot_layout_2d(layout::DataFrame; kwargs...)

Create a new figure and plot a 2D EEG electrode layout.

# Arguments
- `layout`: DataFrame containing electrode positions
- `kwargs...`: Keyword arguments passed to the base plot_layout_2d function

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xy!(layout)
    plot_layout_2d(layout)
"""
function plot_layout_2d(layout::DataFrame; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_layout_2d(layout, neighbours; kwargs...)

Create a new figure and plot a 2D EEG electrode layout with interactive points showing electrode connections.

# Arguments
- `layout`: DataFrame containing electrode positions
- `neighbours`: OrderedDict mapping electrode symbols to their neighboring electrodes
- `kwargs...`: Additional keyword arguments passed to the plot_layout_2d! function

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    neighbours, nneighbours = get_electrode_neighbours_xy(layout, 80)
    fig, ax = plot_layout_2d(layout, neighbours)
"""
function plot_layout_2d(layout, neighbours; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout, neighbours; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame, 
                   neighbours::OrderedDict; kwargs...)

Create a 2D EEG electrode layout with interactive points showing electrode connections.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: DataFrame containing electrode positions
- `neighbours`: OrderedDict mapping electrode symbols to their neighboring electrodes
- `kwargs...`: Additional keyword arguments passed to the base plot_layout_2d function

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    neighbours, nneighbours = get_electrode_neighbours_xy(layout, 80)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout, neighbours)
"""
function plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame, neighbours::OrderedDict; kwargs...)
    plot_layout_2d!(fig, ax, layout; point_kwargs = Dict(:plot_points => false), kwargs...)
    positions = Observable(Point2f.(layout.x2, layout.y2))
    add_interactive_points!(fig, ax, layout, neighbours, positions)
    return fig, ax
end

# using LibGEOS package for convexhull
# TODO: concave hull?
"""
    create_convex_hull(xpos::Vector{<:Real}, ypos::Vector{<:Real}, border_size::Real)

Create a convex hull around a set of 2D points with a specified border size.

# Arguments
- `xpos`: Array of x-coordinates
- `ypos`: Array of y-coordinates
- `border_size`: Size of the border around points

# Returns
- A LibGEOS convex hull polygon

"""
function create_convex_hull(xpos::Vector{<:Real}, ypos::Vector{<:Real}, border_size::Real)
    circle_points = 0:2*pi/361:2*pi
    xs = (border_size.*sin.(circle_points).+transpose(xpos))[:]
    ys = (border_size.*cos.(circle_points).+transpose(ypos))[:]
    xys = [[xs[i], ys[i]] for i in eachindex(xs)]
    push!(xys, xys[1])
    poly = LibGEOS.Polygon([xys])
    hull = LibGEOS.convexhull(poly)
    return hull
end

"""
    add_topo_rois!(ax::Axis, layout::DataFrame, rois::Vector{<:Vector{Symbol}};
                  border_size::Real=10, roi_kwargs::Dict=Dict())

Add regions of interest (ROIs) to a topographic plot based on groups of electrodes.

# Arguments
- `ax`: The axis to add ROIs to
- `layout`: DataFrame containing electrode positions (needs x2, y2)
- `rois`: Array of arrays, where each inner array contains electrode labels for a ROI
- `border_size`: Size of the border around ROI points (default: 10)
- `roi_kwargs`: Keyword arguments for ROI styling (color, linewidth, fill, fillcolor, fillalpha)

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xy!(layout)
    fig, ax = plot_layout_2d(layout)
    # Add ROIs for two electrode groups
    add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 10)
    # Add a filled ROI
    add_topo_rois!(ax, layout, [[:Fp1]], border_size = 10,
                  roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))

"""
function add_topo_rois!(
    ax::Axis,
    layout::DataFrame,
    rois::Vector{<:Vector{Symbol}};
    border_size::Real = 10,
    roi_kwargs::Dict = Dict(),
)

    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    roi_default_kwargs = Dict(
        :color => repeat([:black], length(rois)),
        :linewidth => repeat([2], length(rois)),
        :fill => repeat([false], length(rois)),
        :fillcolor => repeat([:black], length(rois)),
        :fillalpha => repeat([0.2], length(rois)),
    )

    merged_roi_kwargs = merge(roi_default_kwargs, roi_kwargs)
    for (idx, roi) in enumerate(rois)
        xpos = filter(row -> row.label ∈ roi, layout).x2
        ypos = filter(row -> row.label ∈ roi, layout).y2
        border = create_convex_hull(xpos, ypos, border_size)
        lines!(ax, border, linewidth = merged_roi_kwargs[:linewidth][idx], color = merged_roi_kwargs[:color][idx])
        if merged_roi_kwargs[:fill][idx]
            poly!(ax, border, color = merged_roi_kwargs[:fillcolor][idx], alpha = merged_roi_kwargs[:fillalpha][idx])
        end
    end
end

########################################################
# 3D layout
########################################################
"""
    plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame;
                  point_kwargs::Dict=Dict(), label_kwargs::Dict=Dict())

Plot a 3D EEG electrode layout with customizable electrode points and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: DataFrame containing electrode positions with columns x3, y3, z3, and label
- `point_kwargs`: Keyword arguments for electrode points (plot_points, marker, size, color etc.)
- `label_kwargs`: Keyword arguments for electrode labels (plot_labels, fontsize, color, xoffset, yoffset, zoffset etc.)

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xyz!(layout)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout)
"""
function plot_layout_3d!(
    fig::Figure,
    ax::Axis3,
    layout::DataFrame;
    point_kwargs::Dict = Dict(),
    label_kwargs::Dict = Dict(),
)

    if (:x3 ∉ propertynames(layout) || :y3 ∉ propertynames(layout) || :z3 ∉ propertynames(layout))
        polar_to_cartesian_xyz!(layout)
    end

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    merged_point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = merged_point_kwargs[:plot_points]
    point_plot_kwargs = filter(p -> p.first != :plot_points, merged_point_kwargs)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0, :zoffset => 0)
    merged_label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = merged_label_kwargs[:plot_labels]

    xoffset = merged_label_kwargs[:xoffset]
    yoffset = merged_label_kwargs[:yoffset]
    zoffset = merged_label_kwargs[:zoffset]
    label_plot_kwargs = filter(p -> p.first ∉ (:plot_labels, :xoffset, :yoffset, :zoffset), merged_label_kwargs)

    if plot_points
        scatter!(ax, layout[!, :x3], layout[!, :y3], layout[!, :z3]; point_plot_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(
                ax,
                position = (label.x3 + xoffset, label.y3 + yoffset, label.z3 + zoffset),
                String(label.label);
                label_plot_kwargs...,
            )
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return fig, ax
end

"""
    plot_layout_3d(layout; kwargs...)

Create a new figure and plot a 3D EEG electrode layout.

# Arguments
- `layout`: DataFrame containing electrode positions
- `kwargs...`: Additional keyword arguments passed to the plot_layout_3d function

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xyz!(layout)
    fig, ax = plot_layout_3d(layout)
"""
function plot_layout_3d(layout; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout; kwargs...)
    display(fig)
    return fig, ax
end

"""
    plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame, 
                   neighbours::OrderedDict; kwargs...)

Create a 3D EEG electrode layout with interactive points showing electrode connections.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: DataFrame containing electrode positions
- `neighbours`: OrderedDict mapping electrode symbols to their neighboring electrodes
- `kwargs...`: Additional keyword arguments passed to the base plot_layout_3d function

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xyz!(layout)
    neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 40)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout, neighbours)
"""
function plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame, neighbours::OrderedDict; kwargs...)
    plot_layout_3d!(fig, ax, layout; point_kwargs = Dict(:plot_points => false), kwargs...)
    positions = Observable(Point3f.(layout.x3, layout.y3, layout.z3))
    add_interactive_points!(fig, ax, layout, neighbours, positions, true)
    return fig, ax
end

function plot_layout_3d(layout::DataFrame, neighbours::OrderedDict; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout, neighbours; kwargs...)
    display(fig)
    return fig, ax
end


"""
    add_interactive_points!(fig::Figure, ax::Union{Axis, Axis3}, layout::DataFrame,
                          neighbours::OrderedDict, positions::Observable, is_3d::Bool=false)

Add interactive electrode points that highlight and show connections to neighboring electrodes on hover.

# Arguments
- `fig`: The figure to add interactivity to
- `ax`: The axis to add interactivity to (can be 2D or 3D)
- `layout`: DataFrame containing electrode information
- `neighbours`: OrderedDict mapping electrode symbols to their neighboring electrodes
- `positions`: Observable containing point positions (Point2f or Point3f)
- `is_3d`: Boolean indicating if the plot is 3D (default: false)

# Returns
- The figure and axis objects

"""
function add_interactive_points!(
    fig::Figure,
    ax::Union{Axis,Axis3},
    layout::DataFrame,
    neighbours::OrderedDict,
    positions::Observable,
    is_3d::Bool = false,
)
    base_size = 15
    hover_size = 25
    sizes = Observable(fill(base_size, length(layout.label)))

    # Add interactive scatter points
    p = scatter!(ax, positions; color = :black, markersize = sizes, inspectable = true, markerspace = :pixel)

    # Initialize line segments
    linesegments = Observable(is_3d ? Point3f[] : Point2f[])
    lines!(ax, linesegments, color = :gray, linewidth = 3)

    # Add hover interaction
    on(events(fig).mouseposition) do mp
        plt, i = pick(fig)
        if plt == p
            # Reset all sizes to base size
            new_sizes = fill(base_size, length(layout.label))
            new_sizes[i] = hover_size
            sizes[] = new_sizes

            # Create lines to neighboring electrodes
            hovered_pos = positions[][i]
            new_lines = is_3d ? Point3f[] : Point2f[]
            for neighbor in neighbours[Symbol(layout.label[i])].electrodes
                neighbor_idx = findfirst(==(neighbor), layout.label)
                push!(new_lines, hovered_pos, positions[][neighbor_idx])
            end
            linesegments[] = new_lines
        end
    end

    return fig, ax
end

# Basic Tests
# TODO: Implement proper tests
#
# basic layouts
# layout = read_layout("./layouts/biosemi72.csv");
# layout = read_layout("./layouts/biosemi64.csv");
#  
# # 2D layout
# polar_to_cartesian_xy!(layout)
# plot_layout_2d(layout);
# 
# # 2D layout with ROIs
# fig, ax = plot_layout_2d(layout)
# add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 10)
# add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 5)
# add_topo_rois!(ax, layout, [[:Fp1]], border_size = 10, roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))
# add_topo_rois!(ax, layout, [[:CPz, :C2, :FCz,  :C1]], border_size = 5, roi_kwargs = Dict(:fill => [true], :fillcolor => [:blue], :fillalpha => [0.2]))
# 
# # 2D layout with neighbours
# neighbours, nneighbours = get_electrode_neighbours_xy(layout, 80);
# plot_layout_2d(layout, neighbours)
# 
# # 3D layout
# polar_to_cartesian_xyz!(layout)
# plot_layout_3d(layout)
# 
# # 3D layout with neighbours
# neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 40);
# plot_layout_3d(layout, neighbours)