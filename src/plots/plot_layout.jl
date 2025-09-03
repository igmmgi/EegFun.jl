########################################################
# Layout plotting functions
########################################################

# Default parameters for layout plots with descriptions
# Dict is used for documentation and for defaults
# Using prefixed keys to avoid conflicts between different component types
const PLOT_LAYOUT_HEAD_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :head_color => (DEFAULT_HEAD_COLOR, "Color of the head shape outline."),
    :head_linewidth => (DEFAULT_HEAD_LINEWIDTH, "Line width of the head shape outline."),
)

const PLOT_LAYOUT_POINT_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :point_plot => (true, "Whether to plot electrode points."),
    :point_marker => (DEFAULT_POINT_MARKER, "Marker style for electrode points."),
    :point_markersize => (DEFAULT_POINT_SIZE, "Size of electrode point markers."),
    :point_color => (DEFAULT_POINT_COLOR, "Color of electrode points."),
)

const PLOT_LAYOUT_LABEL_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :label_plot => (true, "Whether to plot electrode labels."),
    :label_fontsize => (DEFAULT_LABEL_FONTSIZE, "Font size for electrode labels."),
    :label_color => (DEFAULT_LABEL_COLOR, "Color of electrode labels."),
    :label_xoffset => (0, "X-axis offset for electrode labels."),
    :label_yoffset => (0, "Y-axis offset for electrode labels."),
    :label_zoffset => (0, "Z-axis offset for electrode labels (3D only)."),
)

const PLOT_LAYOUT_ROI_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :roi_color => (:black, "Color of ROI outline."),
    :roi_linewidth => (2, "Line width of ROI outline."),
    :roi_fill => (false, "Whether to fill the ROI area."),
    :roi_fillcolor => (:gray, "Color of ROI fill."),
    :roi_fillalpha => (0.2, "Transparency of ROI fill."),
)




"""
    plot_layout_2d!(fig::Figure, ax::Axis, layout::Layout; neighbours::Bool=false, kwargs...)

Plot a 2D EEG electrode layout with customizable head shape, electrode points, and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: Layout containing electrode positions with columns x2, y2, and label
- `neighbours`: Boolean to show interactive neighbour connections (default: false)

$(generate_kwargs_doc(PLOT_LAYOUT_HEAD_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_POINT_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_LABEL_KWARGS))

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    layout = Layout("biosemi64.csv")
    polar_to_cartesian_xy!(layout)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout)
"""
function plot_layout_2d!(
    fig::Figure,
    ax::Axis,
    layout::Layout;
    neighbours::Bool = false,
    kwargs...
)

    _ensure_coordinates_2d!(layout)

    # Handle each component's kwargs directly using prefixes (no cross-component validation)
    head_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_HEAD_KWARGS, kwargs; validate=false)
    point_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_POINT_KWARGS, kwargs; validate=false)
    label_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_LABEL_KWARGS, kwargs; validate=false)
    
    # Extract control variables
    plot_points = point_kwargs[:point_plot]
    plot_labels = label_kwargs[:label_plot]
    xoffset = label_kwargs[:label_xoffset]
    yoffset = label_kwargs[:label_yoffset]

    # Head shape - Use constants
    radius = DEFAULT_HEAD_RADIUS
    arc!(ax, Point2f(0), radius * 2, -π, π; color=head_kwargs[:head_color], linewidth=head_kwargs[:head_linewidth]) # head
    arc!(ax, Point2f(radius * 2, 0), radius * 2 * HEAD_EAR_RATIO, -π / 2, π / 2; color=head_kwargs[:head_color], linewidth=head_kwargs[:head_linewidth]) # ear right
    arc!(ax, Point2f(-radius * 2, 0), -radius * 2 * HEAD_EAR_RATIO, π / 2, -π / 2; color=head_kwargs[:head_color], linewidth=head_kwargs[:head_linewidth]) # ear left
    lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * HEAD_NOSE_SCALE; color=head_kwargs[:head_color], linewidth=head_kwargs[:head_linewidth]) # nose

    # Regular points
    if plot_points
        scatter!(ax, layout.data[!, :x2], layout.data[!, :y2]; 
                marker=point_kwargs[:point_marker], 
                markersize=point_kwargs[:point_markersize], 
                color=point_kwargs[:point_color])
    end

    if plot_labels
        x_coords = layout.data[!, :x2] .+ xoffset
        y_coords = layout.data[!, :y2] .+ yoffset
        labels = String.(layout.data[!, :label])
        for i in eachindex(labels)
            text!(ax, position = (x_coords[i], y_coords[i]), labels[i]; 
                  fontsize=label_kwargs[:label_fontsize], 
                  color=label_kwargs[:label_color])
        end
    end

    if neighbours
        if isnothing(layout.neighbours)
            @minimal_warning "Layout has no neighbours data. Set neighbours=false or calculate neighbours first."
        else
            positions = Point2f.(layout.data.x2, layout.data.y2)
            _add_interactive_points!(fig, ax, layout.data, layout.neighbours, positions)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return nothing
end


"""
    plot_layout_2d(layout::Layout; neighbours::Bool=false, display_plot::Bool=true, kwargs...)

Create a new figure and plot a 2D EEG electrode layout.

# Arguments
- `layout`: Layout containing electrode positions
- `neighbours`: Boolean to show interactive neighbour connections (default: false)
- `display_plot`: Boolean to display the plot (default: true)

$(generate_kwargs_doc(PLOT_LAYOUT_HEAD_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_POINT_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_LABEL_KWARGS))

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xy!(layout)
    plot_layout_2d(layout)
    # With neighbour interactivity
    plot_layout_2d(layout, neighbours=true)
"""
function plot_layout_2d(layout::Layout; neighbours::Bool = false, display_plot::Bool = true, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])

    plot_layout_2d!(fig, ax, layout; neighbours = neighbours, kwargs...)

    if display_plot
        display_figure(fig)
    end
    return fig, ax
end



"""
    _create_convex_hull_graham(xpos::Vector{<:Real}, ypos::Vector{<:Real}, border_size::Real)

Create a convex hull around a set of 2D points with a specified border size using Graham's Scan algorithm.

# Arguments
- `xpos`: Array of x-coordinates
- `ypos`: Array of y-coordinates
- `border_size`: Size of the border around points

# Returns
- A Vector of 2D points forming the convex hull
"""
function _create_convex_hull_graham(xpos::Vector{<:Real}, ypos::Vector{<:Real}, border_size::Real)
    # Generate points around each electrode with the border
    circle_points = 0:(2π/361):2π
    xs = (border_size .* sin.(circle_points) .+ transpose(xpos))[:]
    ys = (border_size .* cos.(circle_points) .+ transpose(ypos))[:]

    # Convert to array of points
    points = [[xs[i], ys[i]] for i in eachindex(xs)]
    n = length(points)

    # Find the bottommost point (and leftmost if tied)
    ymin = minimum(p -> p[2], points)
    p0 = points[findfirst(p -> p[2] == ymin, points)]

    # Sort points by polar angle with respect to p0
    sort!(points, by = p -> begin
        if p == p0
            return -Inf
        end
        return atan(p[2] - p0[2], p[1] - p0[1])
    end)

    # Initialize stack for Graham's scan
    stack = Vector{Vector{Float64}}()
    push!(stack, points[1])
    push!(stack, points[2])

    # Process remaining points
    for i = 3:n
        while length(stack) > 1 && _orientation(stack[end-1], stack[end], points[i]) != 2
            pop!(stack)
        end
        push!(stack, points[i])
    end

    # Close the hull by connecting back to the first point
    push!(stack, stack[1])

    return stack
end


"""
    add_topo_rois!(ax::Axis, layout::DataFrame, rois::Vector{<:Vector{Symbol}}; border_size::Real=10, kwargs...)

Add regions of interest (ROIs) to a topographic plot based on groups of electrodes.
Uses Graham's Scan algorithm to create convex hulls around electrode groups.

# Arguments
- `ax`: The axis to add ROIs to
- `layout`: DataFrame containing electrode positions (needs x2, y2)
- `rois`: Array of arrays, where each inner array contains electrode labels for a ROI
- `border_size`: Size of the border around ROI points (default: 10)

$(generate_kwargs_doc(PLOT_LAYOUT_ROI_KWARGS))

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xy!(layout)
    fig, ax = plot_layout_2d(layout)
    # Add simple ROI
    add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1]], border_size=10)
    # Add filled ROI
    add_topo_rois!(ax, layout, [[:Fp1]], border_size=5, fill=true, fillcolor=:red, fillalpha=0.2)
"""
function add_topo_rois!(
    ax::Axis,
    layout::DataFrame,
    rois::Vector{<:Vector{Symbol}};
    border_size::Real = 10,
    kwargs...
)
    # Merge user kwargs with defaults for ROI component only
    roi_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_ROI_KWARGS, kwargs; validate=false)
    
    # Handle multiple ROIs with different styles
    n_rois = length(rois)

    # Expand single values to arrays if needed
    merged_kwargs = Dict{Symbol,Vector}()
    for (key, value) in roi_kwargs
        # Expand single values to arrays if needed
        if !isa(value, Vector)
            merged_kwargs[key] = fill(value, n_rois)
        elseif length(value) != n_rois
            throw(ArgumentError("Length of $(key) ($(length(value))) must match number of ROIs ($n_rois)"))
        else
            merged_kwargs[key] = value
        end
    end

    # Create ROIs
    for (i, roi) in enumerate(rois)
        # Get coordinates for ROI electrodes
        roi_idx = findall(in(roi), layout.label)
        if isempty(roi_idx)
            @minimal_warning "No electrodes found for ROI $i"
            continue
        end

        # Create convex hull
        hull_points = _create_convex_hull_graham(layout.x2[roi_idx], layout.y2[roi_idx], border_size)
        xs = [p[1] for p in hull_points]
        ys = [p[2] for p in hull_points]

        # Draw the ROI
        if merged_kwargs[:roi_fill][i]
            poly!(
                ax,
                xs,
                ys;
                color = (merged_kwargs[:roi_fillcolor][i], merged_kwargs[:roi_fillalpha][i]),
                strokecolor = merged_kwargs[:roi_color][i],
                strokewidth = merged_kwargs[:roi_linewidth][i],
            )
        else
            lines!(ax, xs, ys; color = merged_kwargs[:roi_color][i], linewidth = merged_kwargs[:roi_linewidth][i])
        end
    end
end



"""
    plot_layout_3d!(fig::Figure, ax::Axis3, layout::Layout; neighbours::Bool=false, kwargs...)

Plot a 3D EEG electrode layout with customizable head shape, electrode points, and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: Layout containing electrode positions with columns x3, y3, z3, and label
- `neighbours`: Boolean to show interactive neighbour connections (default: false)

$(generate_kwargs_doc(PLOT_LAYOUT_HEAD_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_POINT_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_LABEL_KWARGS))

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    layout = Layout("biosemi64.csv")
    polar_to_cartesian_xyz!(layout)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout)
"""
function plot_layout_3d!(
    fig::Figure,
    ax::Axis3,
    layout::Layout;
    neighbours::Bool = false,
    kwargs...
)

    _ensure_coordinates_3d!(layout)

    # Handle each component's kwargs directly using prefixes (no cross-component validation)
    head_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_HEAD_KWARGS, kwargs; validate=false)
    point_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_POINT_KWARGS, kwargs; validate=false)
    label_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_LABEL_KWARGS, kwargs; validate=false)
    
    # Extract control variables
    plot_points = point_kwargs[:point_plot]
    plot_labels = label_kwargs[:label_plot]
    xoffset = label_kwargs[:label_xoffset]
    yoffset = label_kwargs[:label_yoffset]
    zoffset = label_kwargs[:label_zoffset]

    # Regular points
    if plot_points
        scatter!(ax, layout.data[!, :x3], layout.data[!, :y3], layout.data[!, :z3]; 
                marker=point_kwargs[:point_marker], 
                markersize=point_kwargs[:point_markersize], 
                color=point_kwargs[:point_color])
    end

    if plot_labels
        x_coords = layout.data[!, :x3] .+ xoffset
        y_coords = layout.data[!, :y3] .+ yoffset
        z_coords = layout.data[!, :z3] .+ zoffset
        labels = String.(layout.data[!, :label])
        for i in eachindex(labels)
            text!(ax, position = (x_coords[i], y_coords[i], z_coords[i]), labels[i]; 
                  fontsize=label_kwargs[:label_fontsize], 
                  color=label_kwargs[:label_color])
        end
    end

    if neighbours
        if isnothing(layout.neighbours)
            @minimal_warning "Layout has no neighbours data. Set neighbours=false or calculate neighbours first."
        else
            positions = Point3f.(layout.data.x3, layout.data.y3, layout.data.z3)
            _add_interactive_points!(fig, ax, layout.data, layout.neighbours, positions, true)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return fig, ax
end


"""
    plot_layout_3d(layout::Layout; neighbours::Bool=false, display_plot::Bool=true, kwargs...)

Create a new figure and plot a 3D EEG electrode layout.

# Arguments
- `layout`: Layout containing electrode positions
- `neighbours`: Boolean to show interactive neighbour connections (default: false)
- `display_plot`: Boolean to display the plot (default: true)

$(generate_kwargs_doc(PLOT_LAYOUT_HEAD_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_POINT_KWARGS))
$(generate_kwargs_doc(PLOT_LAYOUT_LABEL_KWARGS))

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    polar_to_cartesian_xyz!(layout)
    plot_layout_3d(layout)
    # With neighbour interactivity
    plot_layout_3d(layout, neighbours=true)
"""
function plot_layout_3d(layout::Layout; neighbours::Bool = false, display_plot::Bool = true, kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])

    plot_layout_3d!(fig, ax, layout; neighbours = neighbours, kwargs...)

    if display_plot
        display_figure(fig)
    end

    return fig, ax

end






"""
    _add_interactive_points!(fig::Figure, ax::Union{Axis, Axis3}, layout::DataFrame,
                          neighbours::OrderedDict, positions::Union{Vector{Point2f}, Vector{Point3f}}, is_3d::Bool=false)

Add interactive electrode points that highlight and show connections to neighboring electrodes on hover.

# Arguments
- `fig`: The figure to add interactivity to
- `ax`: The axis to add interactivity to (can be 2D or 3D)
- `layout`: DataFrame containing electrode information
- `neighbours`: OrderedDict mapping electrode symbols to their neighboring electrodes
- `positions`: Vector containing point positions (Point2f or Point3f) - no longer Observable
- `is_3d`: Boolean indicating if the plot is 3D (default: false)

# Returns
- The figure and axis objects
"""
function _add_interactive_points!(
    fig::Figure,
    ax::Union{Axis,Axis3},
    layout::DataFrame,
    neighbours::OrderedDict,
    positions::Union{Vector{Point2f},Vector{Point3f}},  # Keep as regular vector
    is_3d::Bool = false,
)
    base_size = 15
    hover_size = 25
    sizes = Observable(fill(base_size, length(layout.label)))  # Observable for reactivity

    # Pre-compute neighbor indices for better performance
    neighbor_indices = Vector{Vector{Int}}(undef, length(layout.label))
    for (i, label) in enumerate(layout.label)
        neighbor_indices[i] = Int[]
        if haskey(neighbours, Symbol(label))
            for neighbor in neighbours[Symbol(label)].electrodes
                neighbor_idx = findfirst(==(neighbor), layout.label)
                if neighbor_idx !== nothing
                    push!(neighbor_indices[i], neighbor_idx)
                end
            end
        end
    end

    # Pre-allocate vectors for better performance
    new_sizes = fill(base_size, length(layout.label))
    max_neighbors = maximum(length.(neighbor_indices), init = 0)
    new_lines = is_3d ? Vector{Point3f}(undef, max_neighbors * 2) : Vector{Point2f}(undef, max_neighbors * 2)

    # Add interactive scatter points
    p = scatter!(ax, positions; color = :black, markersize = sizes, inspectable = true, markerspace = :pixel)

    # Initialize line segments
    linesegments = Observable(is_3d ? Point3f[] : Point2f[])  # Observable for reactivity
    lines!(ax, linesegments, color = :gray, linewidth = 3)

    # Add hover interaction
    on(events(fig).mouseposition) do mp
        plt, i = pick(fig)
        if plt == p
            # Reset all sizes to base size and set hover size
            fill!(new_sizes, base_size)
            new_sizes[i] = hover_size
            sizes[] = new_sizes

            # Create lines to neighboring electrodes
            hovered_pos = positions[i]
            line_count = 0

            # Use pre-computed neighbor indices
            for neighbor_idx in neighbor_indices[i]
                new_lines[line_count+1] = hovered_pos
                new_lines[line_count+2] = positions[neighbor_idx]
                line_count              += 2
            end

            # Only update with the actual lines we created
            linesegments[] = new_lines[1:line_count]
        end
    end
end
