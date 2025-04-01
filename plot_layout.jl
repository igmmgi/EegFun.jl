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
- `point_kwargs`: Keyword arguments for electrode points (marker, size, etc.)
- `label_kwargs`: Keyword arguments for electrode labels (fontsize, color, etc.)

# Returns
- The figure and axis objects
"""
function plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame; 
                        head_kwargs::Dict = Dict(), 
                        point_kwargs::Dict = Dict(), 
                        label_kwargs::Dict = Dict())
    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    # Default kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = pop!(point_kwargs, :plot_points)

    label_default_kwargs = Dict(:plot_labels => true, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = pop!(label_kwargs, :plot_labels)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    # Head shape
    radius = 88 # mm
    arc!(ax, Point2f(0), radius * 2, -π, π; head_kwargs...) # head
    arc!(Point2f(radius * 2, 0), radius * 2 / 7, -π / 2, π / 2; head_kwargs...) # ear right
    arc!(Point2f(-radius * 2, 0), -radius * 2 / 7, π / 2, -π / 2; head_kwargs...) # ear left
    lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * 4; head_kwargs...) # nose

    # Regular points
    if plot_points
        scatter!(ax, layout[!, :x2], layout[!, :y2]; point_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(ax, position = (label.x2 + xoffset, label.y2 + yoffset), String(label.label); label_kwargs...)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return fig, ax
end

"""
    plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame; 
                  point_kwargs::Dict=Dict(), label_kwargs::Dict=Dict())

Plot a 3D EEG electrode layout with customizable electrode points and labels.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: DataFrame containing electrode positions with columns x3, y3, z3, and label
- `point_kwargs`: Keyword arguments for electrode points (marker, size, etc.)
- `label_kwargs`: Keyword arguments for electrode labels (fontsize, color, etc.)

# Returns
- The figure and axis objects
"""
function plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame; 
                        point_kwargs::Dict = Dict(), 
                        label_kwargs::Dict = Dict())
    if (:x3 ∉ propertynames(layout) || :y3 ∉ propertynames(layout) || :z3 ∉ propertynames(layout))
        polar_to_cartesian_xyz!(layout)
    end

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = pop!(point_kwargs, :plot_points)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0, :zoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = pop!(label_kwargs, :plot_labels)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)
    zoffset = pop!(label_kwargs, :zoffset)

    if plot_points
        scatter!(ax, layout[!, :x3], layout[!, :y3], layout[!, :z3]; point_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(
                ax,
                position = (label.x3 + xoffset, label.y3 + yoffset, label.z3 + zoffset),
                String(label.label);
                label_kwargs...,
            )
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

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
function add_interactive_points!(fig::Figure, ax::Union{Axis, Axis3}, 
                               layout::DataFrame, 
                               neighbours::OrderedDict, 
                               positions::Observable, 
                               is_3d::Bool = false)
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

"""
    plot_layout_2d(layout; kwargs...)

Create a new figure and plot a 2D EEG electrode layout.

# Arguments
- `layout`: DataFrame containing electrode positions
- `kwargs...`: Keyword arguments passed to the base plot_layout_2d function

# Returns
- The figure and axis objects
"""
function plot_layout_2d(layout::DataFrame; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d!(fig, ax, layout; kwargs...)
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
"""
function plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame, 
                        neighbours::OrderedDict; kwargs...)
    plot_layout_2d!(fig, ax, layout; point_kwargs = Dict(:plot_points => false), kwargs...)
    positions = Observable(Point2f.(layout.x2, layout.y2))
    add_interactive_points!(fig, ax, layout, neighbours, positions)
    return fig, ax
end

function plot_layout_3d(layout::DataFrame; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d!(fig, ax, layout; kwargs...)
    display(fig)
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
"""
function plot_layout_3d!(fig::Figure, ax::Axis3, layout::DataFrame, 
                        neighbours::OrderedDict; kwargs...)
    plot_layout_3d!(fig, ax, layout; point_kwargs = Dict(:plot_points => false), kwargs...)
    positions = Observable(Point3f.(layout.x3, layout.y3, layout.z3))
    add_interactive_points!(fig, ax, layout, neighbours, positions, true)
    return fig, ax
end

# using LibGEOS package (convexhull)
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
- `layout`: DataFrame containing electrode positions
- `rois`: Array of arrays, where each inner array contains electrode labels for a ROI
- `border_size`: Size of the border around ROI points (default: 10)
- `roi_kwargs`: Keyword arguments for ROI styling (color, fill, etc.)
"""
function add_topo_rois!(ax::Axis, layout::DataFrame, rois::Vector{<:Vector{Symbol}}; 
                       border_size::Real = 10, roi_kwargs::Dict = Dict())
    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end
    roi_default_kwargs = Dict(
        :color => repeat([:black], length(rois)),
        :linewidth => repeat([2], length(rois)),
        :fill => repeat([false], length(rois)),
        :fillcolor => repeat([:black], length(rois)),
        :fillalpha => repeat([0.2], length(rois)),
    )
    roi_kwargs = merge(roi_default_kwargs, roi_kwargs)
    for (idx, roi) in enumerate(rois)
        xpos = filter(row -> row.label ∈ roi, layout).x2
        ypos = filter(row -> row.label ∈ roi, layout).y2
        border = create_convex_hull(xpos, ypos, border_size)
        lines!(ax, border, linewidth = roi_kwargs[:linewidth][idx], color = roi_kwargs[:color][idx])
        if roi_kwargs[:fill][idx]
            poly!(ax, border, color = roi_kwargs[:fillcolor][idx], alpha = roi_kwargs[:fillalpha][idx])
        end
    end
end

function plot_layout_2d(fig::Figure, ax::Axis, layout::DataFrame, 
                        neighbours::OrderedDict; kwargs...)
    plot_layout_2d!(fig, ax, layout, neighbours; kwargs...)
    return fig, ax
end

function plot_layout_3d(fig::Figure, ax::Axis3, layout::DataFrame, 
                        neighbours::OrderedDict; kwargs...)
    plot_layout_3d!(fig, ax, layout, neighbours; kwargs...)
    return fig, ax
end

"""
Structure representing head shape configuration for EEG visualizations.

# Fields
- `radius`: Base radius for head shape in mm
- `style`: Head shape style (:standard, :simple, :none)
- `ear_ratio`: Ratio of ear radius to head radius
- `nose_size`: Size of nose relative to head radius
- `nose_points`: Points defining nose shape 
- `kwargs`: Additional styling parameters
"""
struct HeadShape
    radius::Float64
    style::Symbol
    ear_ratio::Float64
    nose_size::Float64
    nose_points::Vector{Point2f}
    kwargs::Dict{Symbol,Any}
    
    # Constructor with defaults
    function HeadShape(;
        radius::Real = 88.0,
        style::Symbol = :standard,
        ear_ratio::Real = 1/7,
        nose_size::Real = 0.2,
        nose_points = nothing,
        kwargs...)
        
        # Default nose points if not provided
        if isnothing(nose_points)
            nose_points = [Point2f(-0.05, 0.5), Point2f(0.0, 0.6), Point2f(0.05, 0.5)]
        end
        
        new(Float64(radius), style, Float64(ear_ratio), Float64(nose_size), nose_points, Dict{Symbol,Any}(kwargs))
    end
end

"""
    draw_head!(ax, head::HeadShape)

Draw head shape on the given axis according to the head shape configuration.

# Arguments
- `ax`: The axis to draw on
- `head`: HeadShape configuration
"""
function draw_head!(ax, head::HeadShape)
    # Extract parameters
    radius = head.radius
    kwargs = head.kwargs
    
    if head.style == :none
        return  # Don't draw anything
    elseif head.style == :simple
        # Simple circle head
        arc!(ax, Point2f(0), radius * 2, -π, π; kwargs...)
    elseif head.style == :standard
        # Standard head with ears and nose
        arc!(ax, Point2f(0), radius * 2, -π, π; kwargs...)  # head
        arc!(Point2f(radius * 2, 0), radius * 2 * head.ear_ratio, -π/2, π/2; kwargs...)  # ear right
        arc!(Point2f(-radius * 2, 0), -radius * 2 * head.ear_ratio, π/2, -π/2; kwargs...)  # ear left
        lines!(ax, head.nose_points .* (radius * 4 * head.nose_size); kwargs...)  # nose
    else
        error("Unknown head style: $(head.style)")
    end
end

"""
    plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame; 
                  head::Union{HeadShape,Nothing}=nothing,
                  style::Union{Symbol,StylePreset}=:default,
                  head_kwargs::Dict=Dict(),
                  point_kwargs::Dict=Dict(), 
                  label_kwargs::Dict=Dict())

Plot a 2D EEG electrode layout with customizable styling and appearance.

# Arguments
- `fig`: The figure to plot on
- `ax`: The axis to plot on
- `layout`: DataFrame containing electrode positions with columns x2, y2, and label
- `head`: Optional HeadShape configuration
- `style`: Style preset name or StylePreset object
- `head_kwargs`: Keyword arguments for head shape rendering (color, linewidth, etc.)
- `point_kwargs`: Keyword arguments for electrode points (marker, size, etc.)
- `label_kwargs`: Keyword arguments for electrode labels (fontsize, color, etc.)

# Returns
- The figure and axis objects
"""
function plot_layout_2d!(fig::Figure, ax::Axis, layout::DataFrame; 
                       head::Union{HeadShape,Nothing}=nothing,
                       style::Union{Symbol,StylePreset}=:default,
                       head_kwargs::Dict=Dict(),
                       point_kwargs::Dict=Dict(), 
                       label_kwargs::Dict=Dict())
    # Process style preset
    if isa(style, Symbol)
        preset = get_style_preset(style)
    else
        preset = style
    end
    
    # Apply style to figure/axis
    apply_style_preset!(fig, ax, preset)
    
    # Merge style settings with user-provided kwargs
    head_kwargs = merge(preset.head, head_kwargs)
    point_kwargs = merge(preset.electrodes, point_kwargs)
    label_kwargs = merge(preset.labels, label_kwargs)
    
    # Coordinate conversion if needed
    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    # Handle head shape
    if isnothing(head)
        # Use kwargs for head parameters
        head = HeadShape(; head_kwargs...)
    end
    
    # Draw the head
    draw_head!(ax, head)

    # Regular points
    point_default_kwargs = Dict(:plot_points => true)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = pop!(point_kwargs, :plot_points)

    label_default_kwargs = Dict(:plot_labels => true, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = pop!(label_kwargs, :plot_labels)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    if plot_points
        scatter!(ax, layout[!, :x2], layout[!, :y2]; point_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(ax, position = (label.x2 + xoffset, label.y2 + yoffset), String(label.label); label_kwargs...)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return fig, ax
end

"""
A collection of preset styles for EEG visualizations.
Contains predefined combinations of colors, sizes, and formatting options
for consistent visual appearance across different plots.
"""
struct StylePreset
    name::Symbol
    head::Dict{Symbol, Any}
    electrodes::Dict{Symbol, Any}
    labels::Dict{Symbol, Any}
    lines::Dict{Symbol, Any}
    background::Dict{Symbol, Any}
    
    # Constructor with defaults
    function StylePreset(; 
                         name::Symbol = :default,
                         head::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                         electrodes::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                         labels::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                         lines::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                         background::Dict{Symbol, Any} = Dict{Symbol, Any}())
        new(name, head, electrodes, labels, lines, background)
    end
end

"""
    get_style_preset(name::Symbol)

Get a predefined style preset by name.

Available presets:
- `:default`: Standard black and white visualization style
- `:publication`: Clean style suitable for academic publications
- `:presentation`: High contrast style for presentations
- `:dark`: Dark mode with light elements
- `:colorblind`: Colorblind-friendly color scheme

# Returns
- StylePreset object
"""
function get_style_preset(name::Symbol)
    if name == :default
        # Standard black and white style
        return StylePreset(
            name = :default,
            head = Dict(:color => :black, :linewidth => 2),
            electrodes = Dict(:color => :black, :markersize => 12, :marker => :circle),
            labels = Dict(:color => :black, :fontsize => 20),
            lines = Dict(:color => :gray40, :linewidth => 1.5),
            background = Dict(:color => :white)
        )
    elseif name == :publication
        # Clean style for academic publications
        return StylePreset(
            name = :publication,
            head = Dict(:color => :black, :linewidth => 1.5),
            electrodes = Dict(:color => :black, :markersize => 8, :marker => :circle),
            labels = Dict(:color => :black, :fontsize => 14),
            lines = Dict(:color => :gray60, :linewidth => 1.0),
            background = Dict(:color => :white)
        )
    elseif name == :presentation
        # High contrast for presentations
        return StylePreset(
            name = :presentation,
            head = Dict(:color => :black, :linewidth => 3),
            electrodes = Dict(:color => :royalblue, :markersize => 14, :marker => :circle),
            labels = Dict(:color => :black, :fontsize => 22, :font => :bold),
            lines = Dict(:color => :gray30, :linewidth => 2.0),
            background = Dict(:color => :white)
        )
    elseif name == :dark
        # Dark mode with light elements
        return StylePreset(
            name = :dark,
            head = Dict(:color => :white, :linewidth => 2),
            electrodes = Dict(:color => :cyan, :markersize => 12, :marker => :circle),
            labels = Dict(:color => :white, :fontsize => 20),
            lines = Dict(:color => :gray70, :linewidth => 1.5),
            background = Dict(:color => :black)
        )
    elseif name == :colorblind
        # Colorblind-friendly scheme
        return StylePreset(
            name = :colorblind,
            head = Dict(:color => :black, :linewidth => 2),
            electrodes = Dict(:color => :dodgerblue, :markersize => 12, :marker => :circle),
            labels = Dict(:color => :black, :fontsize => 20),
            lines = Dict(:color => :orange, :linewidth => 1.5),
            background = Dict(:color => :white)
        )
    else
        error("Unknown style preset: $name")
    end
end

"""
    apply_style_preset!(fig, ax, preset::StylePreset)

Apply a style preset to a figure and axis.

# Arguments
- `fig`: The figure to apply styling to
- `ax`: The axis to apply styling to
- `preset`: StylePreset to apply
"""
function apply_style_preset!(fig, ax, preset::StylePreset)
    # Apply background
    if haskey(preset.background, :color)
        ax.backgroundcolor = preset.background[:color]
    end
    
    # Other global styling options could be applied here
    return fig, ax
end

