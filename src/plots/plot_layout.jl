########################################################
# Layout plotting functions
########################################################

########################################################
# Helper functions for kwargs handling
########################################################

"""
    _merge_kwargs(defaults::Dict, user_kwargs::Dict)

Helper function to merge default keyword arguments with user-provided arguments.

# Arguments
- `defaults::Dict`: Dictionary of default keyword arguments
- `user_kwargs::Dict`: Dictionary of user-provided keyword arguments

# Returns
- `Dict`: Merged dictionary with user kwargs taking precedence
"""
function _merge_kwargs(defaults::Dict, user_kwargs::Dict)
    return merge(defaults, user_kwargs)
end

"""
    _extract_plot_kwargs(merged_kwargs::Dict, exclude_keys::Vector{Symbol})

Helper function to extract plotting keyword arguments by excluding specific keys.

# Arguments
- `merged_kwargs::Dict`: Dictionary of merged keyword arguments
- `exclude_keys::Vector{Symbol}`: Keys to exclude from the result

# Returns
- `Dict`: Filtered dictionary without the excluded keys
"""
function _extract_plot_kwargs(merged_kwargs::Dict, exclude_keys::Vector{Symbol})
    return filter(p -> p.first ∉ exclude_keys, merged_kwargs)
end

"""
    _validate_layout_columns(layout::DataFrame, required_cols::Vector{Symbol})

Helper function to validate that a layout DataFrame contains required columns.

# Arguments
- `layout::DataFrame`: The layout DataFrame to validate
- `required_cols::Vector{Symbol}`: Vector of required column names

# Throws
- `ArgumentError`: If any required columns are missing
"""
function _validate_layout_columns(layout::Layout, required_cols::Vector{Symbol})
    missing_cols = setdiff(required_cols, propertynames(layout.data))
    if !isempty(missing_cols)
        throw(ArgumentError("Missing required columns: $missing_cols"))
    end
end

########################################################
# 2D layout
########################################################
"""
    _ensure_coordinates_2d!(layout::Layout)

Helper function to ensure 2D coordinates exist in the layout DataFrame.
Converts polar coordinates to Cartesian if needed.

# Arguments
- `layout::DataFrame`: The layout DataFrame to check and potentially convert
"""
function _ensure_coordinates_2d!(layout::Layout)
    if !has_2d_coords(layout)
        @info "Converting polar coordinates to 2D Cartesian coordinates"
        polar_to_cartesian_xy!(layout)
    end
end

"""
    plot_layout_2d!(fig::Figure, ax::Axis, layout::Layout;
                  head_kwargs::Dict=Dict(), point_kwargs::Dict=Dict(),
                  label_kwargs::Dict=Dict())

Plot a 2D EEG electrode layout with customizable head shape, electrode points, and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: Layout containing electrode positions with columns x2, y2, and label
- `head_kwargs`: Keyword arguments for head shape rendering (color, linewidth, etc.)
- `point_kwargs`: Keyword arguments for electrode points (plot_points, marker, size, color etc.)
- `label_kwargs`: Keyword arguments for electrode labels (plot_labels, fontsize, color, xoffset, yoffset, etc.)

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
    head_kwargs::Dict = Dict(),
    point_kwargs::Dict = Dict(),
    label_kwargs::Dict = Dict(),
)

    _validate_layout_columns(layout, [:label])
    _ensure_coordinates_2d!(layout)

    # Get the DataFrame from Layout
    df = layout.data

    # Use helper functions for kwargs handling
    head_default_kwargs = Dict(:color => DEFAULT_HEAD_COLOR, :linewidth => DEFAULT_HEAD_LINEWIDTH)
    merged_head_kwargs = _merge_kwargs(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => DEFAULT_POINT_MARKER, :markersize => DEFAULT_POINT_SIZE, :color => DEFAULT_POINT_COLOR)
    merged_point_kwargs = _merge_kwargs(point_default_kwargs, point_kwargs)
    plot_points = merged_point_kwargs[:plot_points]
    point_plot_kwargs = _extract_plot_kwargs(merged_point_kwargs, [:plot_points])

    label_default_kwargs = Dict(:plot_labels => true, :fontsize => DEFAULT_LABEL_FONTSIZE, :color => DEFAULT_LABEL_COLOR, :xoffset => 0, :yoffset => 0)
    merged_label_kwargs = _merge_kwargs(label_default_kwargs, label_kwargs)
    plot_labels = merged_label_kwargs[:plot_labels]
    xoffset = merged_label_kwargs[:xoffset]
    yoffset = merged_label_kwargs[:yoffset]
    label_plot_kwargs = _extract_plot_kwargs(merged_label_kwargs, [:plot_labels, :xoffset, :yoffset])

    # Head shape - Use constants
    radius = DEFAULT_HEAD_RADIUS
    arc!(ax, Point2f(0), radius * 2, -π, π; merged_head_kwargs...) # head
    arc!(ax, Point2f(radius * 2, 0), radius * 2 * HEAD_EAR_RATIO, -π / 2, π / 2; merged_head_kwargs...) # ear right
    arc!(ax, Point2f(-radius * 2, 0), -radius * 2 * HEAD_EAR_RATIO, π / 2, -π / 2; merged_head_kwargs...) # ear left
    lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * HEAD_NOSE_SCALE; merged_head_kwargs...) # nose

    # Regular points
    if plot_points
        scatter!(ax, df[!, :x2], df[!, :y2]; point_plot_kwargs...)
    end

    if plot_labels
        # More efficient: access columns directly instead of iterating rows
        x_coords = df[!, :x2] .+ xoffset
        y_coords = df[!, :y2] .+ yoffset
        labels = String.(df[!, :label])
        
        for i in eachindex(labels)
            text!(ax, position = (x_coords[i], y_coords[i]), labels[i]; label_plot_kwargs...)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return nothing
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
function plot_layout_2d(layout::Layout; 
    display_plot::Bool = true,
    kwargs...
)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout; kwargs...)
    if display_plot
        display_figure(fig)
    end
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
- `nothing` (modifies the provided figure and axis in-place)

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    neighbours, nneighbours = get_electrode_neighbours_xy(layout, 80)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout, neighbours)
"""
function plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame, neighbours::OrderedDict; kwargs...)
    plot_layout_2d!(fig, ax, layout; point_kwargs = Dict(:plot_points => false), kwargs...)
    positions = Point2f.(layout.x2, layout.y2)
    _add_interactive_points!(fig, ax, layout, neighbours, positions)
    return nothing
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
    circle_points = 0:2π/361:2π
    xs = (border_size.*sin.(circle_points).+transpose(xpos))[:]
    ys = (border_size.*cos.(circle_points).+transpose(ypos))[:]
    
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
    for i in 3:n
        while length(stack) > 1 && orientation(stack[end-1], stack[end], points[i]) != 2
            pop!(stack)
        end
        push!(stack, points[i])
    end
    
    # Close the hull by connecting back to the first point
    push!(stack, stack[1])
    
    return stack
end


"""
    add_topo_rois!(ax::Axis, layout::DataFrame, rois::Vector{<:Vector{Symbol}};
                  border_size::Real=10, roi_kwargs::Dict=Dict())

Add regions of interest (ROIs) to a topographic plot based on groups of electrodes.
Uses Graham's Scan algorithm to create convex hulls around electrode groups.

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
    # Add simple ROI
    add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1]], border_size=10)
    # Add filled ROI
    add_topo_rois!(ax, layout, [[:Fp1]], border_size=5,
                  roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))
"""
function add_topo_rois!(
    ax::Axis,
    layout::DataFrame,
    rois::Vector{<:Vector{Symbol}};
    border_size::Real = 10,
    roi_kwargs::Dict = Dict(),
)
    # Default kwargs
    default_kwargs = Dict(
        :color => :black,
        :linewidth => 2,
        :fill => false,
        :fillcolor => :gray,
        :fillalpha => 0.2,
    )
    
    # Handle multiple ROIs with different styles
    n_rois = length(rois)
    
    # Expand single values to arrays if needed
    for (key, value) in roi_kwargs
        if !isa(value, Vector)
            roi_kwargs[key] = fill(value, n_rois)
        elseif length(value) != n_rois
            throw(ArgumentError("Length of $(key) ($(length(value))) must match number of ROIs ($n_rois)"))
        end
    end
    
    # Merge with defaults, ensuring all values are vectors
    merged_kwargs = Dict{Symbol,Vector}()
    for (key, default_value) in default_kwargs
        if haskey(roi_kwargs, key)
            merged_kwargs[key] = roi_kwargs[key]
        else
            merged_kwargs[key] = fill(default_value, n_rois)
        end
    end
    
    # Create ROIs
    for (i, roi) in enumerate(rois)
        # Get coordinates for ROI electrodes
        roi_idx = findall(in(roi), layout.label)
        if isempty(roi_idx)
            @warn "No electrodes found for ROI $i"
            continue
        end
        
        # Create convex hull
        hull_points = _create_convex_hull_graham(
            layout.x2[roi_idx],
            layout.y2[roi_idx],
            border_size
        )
        xs = [p[1] for p in hull_points]
        ys = [p[2] for p in hull_points]
        
        # Draw the ROI
        if merged_kwargs[:fill][i]
            poly!(
                ax,
                xs,
                ys;
                color = (merged_kwargs[:fillcolor][i], merged_kwargs[:fillalpha][i]),
                strokecolor = merged_kwargs[:color][i],
                strokewidth = merged_kwargs[:linewidth][i]
            )
        else
            lines!(
                ax,
                xs,
                ys;
                color = merged_kwargs[:color][i],
                linewidth = merged_kwargs[:linewidth][i]
            )
        end
    end
end

########################################################
# 3D layout
########################################################
"""
    _ensure_coordinates_3d!(layout::Layout)

Helper function to ensure 3D coordinates exist in the layout DataFrame.
Converts polar coordinates to Cartesian if needed.

# Arguments
- `layout::DataFrame`: The layout DataFrame to check and potentially convert
"""
function _ensure_coordinates_3d!(layout::Layout)
    if !has_3d_coords(layout)
        @info "Converting polar coordinates to 3D Cartesian coordinates"
        polar_to_cartesian_xyz!(layout)
    end
end

"""
    plot_layout_3d!(fig::Figure, ax::Axis3, layout::Layout;
                  head_kwargs::Dict = Dict(),
                  point_kwargs::Dict = Dict(),
                  label_kwargs::Dict = Dict(),
                  display_plot = true,
)

Plot a 3D EEG electrode layout with customizable head shape, electrode points, and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: Layout containing electrode positions with columns x3, y3, z3, and label
- `head_kwargs`: Keyword arguments for head shape rendering (color, linewidth, etc.)
- `point_kwargs`: Keyword arguments for electrode points (plot_points, marker, size, color etc.)
- `label_kwargs`: Keyword arguments for electrode labels (plot_labels, fontsize, color, xoffset, yoffset, zoffset etc.)
- `display_plot`: Boolean indicating whether to display the plot (default: true)

# Returns
- The figure and axis objects

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
    head_kwargs::Dict = Dict(),
    point_kwargs::Dict = Dict(),
    label_kwargs::Dict = Dict(),
    display_plot::Bool = true,
)
    # Validate required columns
    _validate_layout_columns(layout, [:label])

    _ensure_coordinates_3d!(layout)

    # Get the DataFrame from Layout
    df = layout.data

    # Use helper functions for kwargs handling
    head_default_kwargs = Dict(:color => DEFAULT_HEAD_COLOR, :linewidth => DEFAULT_HEAD_LINEWIDTH)
    merged_head_kwargs = _merge_kwargs(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => DEFAULT_POINT_MARKER, :markersize => DEFAULT_POINT_SIZE, :color => DEFAULT_POINT_COLOR)
    merged_point_kwargs = _merge_kwargs(point_default_kwargs, point_kwargs)
    plot_points = merged_point_kwargs[:plot_points]
    point_plot_kwargs = _extract_plot_kwargs(merged_point_kwargs, [:plot_points])

    label_default_kwargs = Dict(:plot_labels => true, :fontsize => DEFAULT_LABEL_FONTSIZE, :color => DEFAULT_LABEL_COLOR, :xoffset => 0, :yoffset => 0, :zoffset => 0)
    merged_label_kwargs = _merge_kwargs(label_default_kwargs, label_kwargs)
    plot_labels = merged_label_kwargs[:plot_labels]
    xoffset = merged_label_kwargs[:xoffset]
    yoffset = merged_label_kwargs[:yoffset]
    zoffset = merged_label_kwargs[:zoffset]
    label_plot_kwargs = _extract_plot_kwargs(merged_label_kwargs, [:plot_labels, :xoffset, :yoffset, :zoffset])

    # Head shape - Use constants
    radius = DEFAULT_HEAD_RADIUS
    arc!(ax, Point3f(0, 0, 0), radius * 2, -π, π; merged_head_kwargs...) # head
    arc!(ax, Point3f(radius * 2, 0, 0), radius * 2 * HEAD_EAR_RATIO, -π / 2, π / 2; merged_head_kwargs...) # ear right
    arc!(ax, Point3f(-radius * 2, 0, 0), -radius * 2 * HEAD_EAR_RATIO, π / 2, -π / 2; merged_head_kwargs...) # ear left
    lines!(ax, Point3f[(-0.05, 0.5, 0), (0.0, 0.6, 0), (0.05, 0.5, 0)] .* radius * HEAD_NOSE_SCALE; merged_head_kwargs...) # nose

    # Regular points
    if plot_points
        scatter!(ax, df[!, :x3], df[!, :y3], df[!, :z3]; point_plot_kwargs...)
    end

    if plot_labels
        # More efficient: access columns directly instead of iterating rows
        x_coords = df[!, :x3] .+ xoffset
        y_coords = df[!, :y3] .+ yoffset
        z_coords = df[!, :z3] .+ zoffset
        labels = String.(df[!, :label])
        
        for i in eachindex(labels)
            text!(ax, position = (x_coords[i], y_coords[i], z_coords[i]), labels[i]; label_plot_kwargs...)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    if display_plot
        display_figure(fig)
    end

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
    if get(kwargs, :display_plot, true) # force a new figure window
        display_figure(fig)
    end
    return fig, ax
end

"""
    plot_layout_3d(layout, neighbours; kwargs...)

Create a new figure and plot a 3D EEG electrode layout with interactive points showing electrode connections.

# Arguments
- `layout`: DataFrame containing electrode positions
- `neighbours`: OrderedDict mapping electrode symbols to their neighboring electrodes
- `kwargs...`: Additional keyword arguments passed to the plot_layout_3d! function

# Returns
- The figure and axis objects

# Example
    layout = read_layout("./layouts/biosemi64.csv")
    neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 80)
    fig, ax = plot_layout_3d(layout, neighbours)
"""
function plot_layout_3d(layout, neighbours; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout, neighbours; kwargs...)
    if get(kwargs, :display_plot, true)
        display_figure(fig)
    end
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
    neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 80)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout, neighbours)
"""
function plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame, neighbours::OrderedDict; kwargs...)
    plot_layout_3d!(fig, ax, layout; point_kwargs = Dict(:plot_points => false), kwargs...)
    positions = Point3f.(layout.x3, layout.y3, layout.z3)
    _add_interactive_points!(fig, ax, layout, neighbours, positions, true)
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
    positions::Union{Vector{Point2f}, Vector{Point3f}},  # Keep as regular vector
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
    max_neighbors = maximum(length.(neighbor_indices), init=0)
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
                new_lines[line_count + 1] = hovered_pos
                new_lines[line_count + 2] = positions[neighbor_idx]
                line_count += 2
            end

            # Only update with the actual lines we created
            linesegments[] = new_lines[1:line_count]
        end
    end
end