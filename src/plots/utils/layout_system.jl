# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const LAYOUT_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Topo layout parameters
    :plot_width => (0.10, "Width of individual plots in topo layout (as fraction of figure width)"),
    :plot_height => (0.10, "Height of individual plots in topo layout (as fraction of figure height)"),
    :margin => (0.02, "Margin between plots in topo layout"),
    :scale_offset_factor => (0.1, "Offset factor for scale plot position"),
    :fallback_scale_x => (0.8, "Fallback x position for scale plot"),
    :fallback_scale_y => (-0.8, "Fallback y position for scale plot"),
    
    # Grid layout parameters
    :grid_dims => (nothing, "Grid dimensions as (rows, cols) tuple. If nothing, automatically calculated"),
    :grid_rowgap => (10, "Gap between rows in grid layout (in pixels)"),
    :grid_colgap => (10, "Gap between columns in grid layout (in pixels)"),
    :grid_padding => ((10, 10, 10, 10), "Padding around grid layout as (left, right, top, bottom) tuple (in pixels)"),
    
    # General layout parameters
    :figure_padding => ((10, 10, 10, 10), "Padding around entire figure as (left, right, top, bottom) tuple (in pixels)"),
    :aspect_ratio => (nothing, "Aspect ratio for individual plots (width/height). If nothing, automatically determined"),
)


"""
    PlotLayout

A struct that defines how to arrange plots in a figure.
"""
struct PlotLayout
    type::Symbol  # :single, :grid, :topo, :custom
    rows::Int
    cols::Int
    positions::Vector{Tuple{Float64,Float64}}  # Custom positions for :custom type
    channels::Vector{Symbol}     # Channels to plot
    metadata::Dict{Symbol,Any}  # Layout configuration parameters (from LAYOUT_KWARGS)
end

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

"""
    _create_layout_metadata(kwargs)

Create metadata dictionary from layout kwargs by merging with defaults.
Accepts kwargs as a NamedTuple, Dict, or Base.Pairs.
"""
function _create_layout_metadata(kwargs)
    layout_kwargs = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)
    return copy(layout_kwargs)
end

"""
    _apply_aspect_ratio!(ax::Axis, plot_layout::PlotLayout)

Apply aspect ratio from plot_layout metadata to axis if specified.
"""
function _apply_aspect_ratio!(ax::Axis, plot_layout::PlotLayout)
    aspect_val = plot_layout.metadata[:aspect_ratio]
    if aspect_val !== nothing
        ax.aspect = aspect_val
    end
end

"""
    _create_axis_with_grid(fig, position; kwargs...)

Create an Axis with grid settings from kwargs.
"""
function _create_axis_with_grid(fig, position; kwargs...)
    return Axis(position,
        xgridvisible = kwargs[:xgrid],
        ygridvisible = kwargs[:ygrid],
        xminorgridvisible = kwargs[:xminorgrid],
        yminorgridvisible = kwargs[:yminorgrid],
    )
end

"""
    _create_topo_axis(fig, plot_width, plot_height, halign, valign; kwargs...)

Create an Axis for topo layout with relative positioning and grid settings.
"""
function _create_topo_axis(fig, plot_width, plot_height, halign, valign; kwargs...)
    return Axis(
        fig[1, 1],
        width = Relative(plot_width),
        height = Relative(plot_height),
        halign = halign,
        valign = valign,
        xgridvisible = kwargs[:xgrid],
        ygridvisible = kwargs[:ygrid],
        xminorgridvisible = kwargs[:xminorgrid],
        yminorgridvisible = kwargs[:yminorgrid],
    )
end

"""
    _apply_grid_spacing!(fig::Figure, plot_layout::PlotLayout)

Apply grid spacing (rowgap/colgap) from plot_layout metadata to figure layout.
"""
function _apply_grid_spacing!(fig::Figure, plot_layout::PlotLayout)
    rowgap_val = plot_layout.metadata[:grid_rowgap]
    if rowgap_val !== nothing
        rowgap!(fig.layout, rowgap_val)
    end
    colgap_val = plot_layout.metadata[:grid_colgap]
    if colgap_val !== nothing
        colgap!(fig.layout, colgap_val)
    end
end

"""
    _add_axis_and_channel!(axes::Vector{Axis}, channels::Vector{Symbol}, ax::Axis, channel::Symbol)

Helper to add axis and channel to their respective vectors.
"""
function _add_axis_and_channel!(axes::Vector{Axis}, channels::Vector{Symbol}, ax::Axis, channel::Symbol)
    push!(axes, ax)
    push!(channels, channel)
end

"""
    create_single_layout(channels::Vector{Symbol}; kwargs...)

Create a single plot layout for multiple channels (all shown on one plot).

# Keyword Arguments
- `kwargs...`: Additional layout parameters (e.g., `figure_padding`, `aspect_ratio`)
"""
function create_single_layout(channels::Vector{Symbol}; kwargs...)
    metadata = _create_layout_metadata(kwargs)
    return PlotLayout(:single, 1, 1, [(0.0, 0.0)], channels, metadata)
end


"""
    create_grid_layout(channels::Vector{Symbol}; rows::Int = 0, cols::Int = 0, kwargs...)

Create a grid layout for multiple channels.
If rows or cols is 0, it will be calculated automatically.

# Keyword Arguments
- `kwargs...`: Additional layout parameters (e.g., `grid_rowgap`, `grid_colgap`, `grid_padding`)
"""
function create_grid_layout(channels::Vector{Symbol}; rows::Int = 0, cols::Int = 0, kwargs...)
    n_channels = length(channels)
    
    # Create metadata first to get grid_dims if provided
    metadata = _create_layout_metadata(kwargs)
    
    # Parse grid_dims from metadata or use function parameters
    grid_dims = get(metadata, :grid_dims, nothing)
    rows, cols = _parse_dims_to_rows_cols(grid_dims, n_channels; default_rows = rows, default_cols = cols)

    if rows * cols < n_channels
        throw(ArgumentError("Grid ($(rows)×$(cols)=$(rows*cols)) is too small for $n_channels channels"))
    end

    println("create_grid_layout: creating $(rows)×$(cols) grid for $(n_channels) channels")
    
    return PlotLayout(:grid, rows, cols, [], channels, metadata)
end

"""
    create_topo_layout(layout::Layout, channels::Vector{Symbol}; kwargs...)

Create a topographic layout based on channel positions.

# Keyword Arguments
$(generate_kwargs_doc(LAYOUT_KWARGS))
"""
function create_topo_layout(
    layout::Layout,
    channels::Vector{Symbol};
    kwargs...,
)
    # Create metadata from kwargs
    topo_kwargs = _create_layout_metadata(kwargs)

    _ensure_coordinates_2d!(layout)

    # Get layout coordinates for bounds calculation (simplified)
    layout_x_coords = layout.data.x2
    layout_y_coords = layout.data.y2

    # Calculate bounds from original layout (fixed, won't change as we add positions)
    layout_minx = isempty(layout_x_coords) ? 0.0 : minimum(layout_x_coords)
    layout_maxx = isempty(layout_x_coords) ? 0.0 : maximum(layout_x_coords)
    layout_miny = isempty(layout_y_coords) ? 0.0 : minimum(layout_y_coords)
    layout_maxy = isempty(layout_y_coords) ? 0.0 : maximum(layout_y_coords)
    layout_xrange = layout_maxx - layout_minx
    layout_yrange = layout_maxy - layout_miny

    # Get channel positions from layout data
    positions = Tuple{Float64,Float64}[]
    new_plot_count = 0  # Count how many new plots we've placed
    last_new_plot_pos = nothing  # Track position of last new plot
    
    for channel in channels
        idx = findfirst(==(channel), layout.data.label)
        if idx !== nothing
            x = layout.data.x2[idx]
            y = layout.data.y2[idx]
            pos = (x, y)
            push!(positions, pos)
        else
            @minimal_warning "Channel $channel not found in layout, finding non-overlapping position"
            new_plot_count += 1
            
            # Calculate spacing based on layout range
            spacing = max(layout_xrange, layout_yrange) * 0.15  # Spacing between new plots
            
            if new_plot_count == 1 # First new plot: to the right, near the top
                new_pos = (layout_maxx + spacing, layout_maxy - spacing * 0.5)
            else # Subsequent new plots: directly below the previous new plot
                new_pos = (last_new_plot_pos[1], last_new_plot_pos[2] - spacing)
            end
            
            push!(positions, new_pos)
            last_new_plot_pos = new_pos  # Update for next new plot
        end
    end

    return PlotLayout(:topo, 0, 0, positions, channels, topo_kwargs)
end

"""
    _normalize_coordinates(positions::Vector{Tuple{Float64, Float64}}, margin::Float64)

Normalize coordinates to [0,1] range and handle edge cases.
Returns normalized coordinates and ranges for layout positioning.
"""
function _normalize_coordinates(positions::Vector{Tuple{Float64,Float64}}, margin::Float64)
    if isempty(positions)
        return Float64[], Float64[], Float64[], Float64[], 1.0, 1.0
    end

    x_coords = [pos[1] for pos in positions]
    y_coords = [pos[2] for pos in positions]
    minx, maxx = extrema(x_coords)
    miny, maxy = extrema(y_coords)

    # Handle case where all coordinates are the same
    xrange = maxx - minx == 0 ? 1.0 : maxx - minx
    yrange = maxy - miny == 0 ? 1.0 : maxy - miny

    return x_coords, y_coords, minx, maxx, miny, maxy, xrange, yrange
end

"""
    create_layout(layout_spec, channels, eeg_layout; kwargs...)

Create a PlotLayout object based on the layout specification.
This is a generic function that can be used by any plot type.

# Keyword Arguments
- `kwargs...`: Additional keyword arguments passed to the specific layout creation function
  (e.g., `plot_width`, `plot_height` for topo layouts)
"""
function create_layout(
    layout_spec::Union{Symbol,PlotLayout,Vector{Int}},
    channels::Vector{Symbol},
    eeg_layout::Union{Layout,Nothing};
    kwargs...,
)

    if layout_spec === :single
        return create_single_layout(channels; kwargs...)
    elseif layout_spec === :grid
        return create_grid_layout(channels; kwargs...)
    elseif layout_spec === :topo
        return create_topo_layout(eeg_layout, channels; kwargs...)
    elseif layout_spec isa Vector{Int}
        if length(layout_spec) != 2
            throw(ArgumentError("layout must be a 2-element vector [rows, cols]"))
        end
        
        # Use _parse_dims_to_rows_cols to handle Vector{Int} consistently
        n_channels = length(channels)
        rows, cols = _parse_dims_to_rows_cols(layout_spec, n_channels)

        # If both dimensions are specified but grid size doesn't fit, warn and use best_rect
        if rows * cols > n_channels
            @minimal_warning "Grid size [$(rows), $(cols)] is too large for $(n_channels) channels. Using optimal grid dimensions instead."
            return create_grid_layout(channels; kwargs...)  # Use best_rect
        elseif rows * cols < n_channels
            @minimal_warning "Grid size [$(rows), $(cols)] is too small for $(n_channels) channels. Using optimal grid dimensions instead."
            return create_grid_layout(channels; kwargs...)  # Use best_rect
        end

        return create_grid_layout(channels; rows = rows, cols = cols, kwargs...)
    elseif layout_spec isa PlotLayout
        return layout_spec
    else
        throw(ArgumentError("Invalid layout specification: $layout_spec"))
    end
end

"""
    create_custom_layout(positions::Vector{Tuple{Float64, Float64}}, channels::Vector{Symbol})

Create a custom layout with specific positions for each channel.
"""
function create_custom_layout(positions::Vector{Tuple{Float64,Float64}}, channels::Vector{Symbol})
    if length(positions) != length(channels)
        throw(
            ArgumentError(
                "Number of positions ($(length(positions))) must match number of channels ($(length(channels)))",
            ),
        )
    end
    return PlotLayout(:custom, 0, 0, positions, channels, Dict{Symbol,Any}())
end

"""
    _parse_dims_to_rows_cols(dims, n_items::Int; default_rows::Int = 0, default_cols::Int = 0)

Parse dims parameter (Tuple, Vector, or nothing) and convert to (rows, cols).
If dims is nothing or invalid, uses default_rows/default_cols or auto-calculates.

# Arguments
- `dims`: Tuple/Vector of (rows, cols), or nothing
- `n_items::Int`: Number of items to arrange (for auto-calculation)
- `default_rows::Int`: Default rows if dims is nothing (0 = auto-calculate)
- `default_cols::Int`: Default cols if dims is nothing (0 = auto-calculate)

# Returns
- `(rows::Int, cols::Int)`: Parsed dimensions
"""
function _parse_dims_to_rows_cols(dims, n_items::Int; default_rows::Int = 0, default_cols::Int = 0)
    rows = default_rows
    cols = default_cols
    
    if dims !== nothing
        if dims isa Tuple && length(dims) == 2
            rows = dims[1] <= 0 ? 0 : dims[1]
            cols = dims[2] <= 0 ? 0 : dims[2]
        elseif dims isa Vector && length(dims) == 2
            rows = dims[1] <= 0 ? 0 : dims[1]
            cols = dims[2] <= 0 ? 0 : dims[2]
        end
    end
    
    # Auto-calculate if both are 0 or not specified
    if rows == 0 && cols == 0
        rows, cols = best_rect(n_items)
    elseif rows == 0
        rows = ceil(Int, n_items / cols)
    elseif cols == 0
        cols = ceil(Int, n_items / rows)
    end
    
    return rows, cols
end

"""
    best_rect(n::Int)

Find the best rectangular layout for n items.
Returns (rows, cols) that minimizes the difference between rows and cols,
preferring arrangements that are as square-like as possible.
For numbers with poor factors (like primes or near-primes), uses approximate
square-like arrangements even if they have a few extra spaces.
"""
function best_rect(n::Int)
    if n <= 0
        throw(ArgumentError("n must be positive"))
    end

    if n == 1
        return (1, 1)
    end

    # Find exact factors of n
    exact_factors = Tuple{Int,Int}[]
    for i = 1:isqrt(n)
        if n % i == 0
            push!(exact_factors, (i, n ÷ i))
        end
    end

    # Filter out 1×n arrangements (unless it's the only option)
    good_factors = Base.filter(f -> f[1] > 1, exact_factors)
    
    if !isempty(good_factors)
        # Use the factor pair with the smallest difference (most square-like)
        best_exact = argmin(good_factors) do (r, c)
            abs(r - c)
        end
        
        # Check if the best exact factor is reasonably square-like
        # If the aspect ratio is too extreme (> 3:1), prefer approximate square
        aspect_ratio = max(best_exact[1], best_exact[2]) / min(best_exact[1], best_exact[2])
        if aspect_ratio <= 3.0
            return best_exact
        end
    end

    # No good exact factors, or aspect ratio too extreme - use sqrt approach
    # This gives a more square-like arrangement even if it has extra spaces
    rows = ceil(Int, sqrt(n))
    cols = ceil(Int, n / rows)
    
    # Ensure we have enough space
    if rows * cols < n
        cols = ceil(Int, n / rows)
    end
    
    return (rows, cols)
end

"""
    apply_layout!(fig::Figure, plot_layout::PlotLayout; kwargs...)

Apply a plot layout to a figure, creating and positioning axes.
Returns the created axes and their associated channels.
The actual plotting should be done separately by the calling function.
"""
function apply_layout!(fig::Figure, plot_layout::PlotLayout; kwargs...)

    axes = Axis[]
    channels = Symbol[]
    
    # Note: figure_padding must be set at Figure creation time, not here
    # It's stored in metadata for reference but should be applied by the caller

    if plot_layout.type == :single
        ax = _create_axis_with_grid(fig, fig[1, 1]; kwargs...)
        _apply_aspect_ratio!(ax, plot_layout)
        _add_axis_and_channel!(axes, channels, ax, plot_layout.channels[1])

    elseif plot_layout.type == :grid

        @info "apply_layout! grid: plot_layout.rows=$(plot_layout.rows), plot_layout.cols=$(plot_layout.cols), channels=$(length(plot_layout.channels))"
        
        # Apply grid spacing if specified in metadata
        _apply_grid_spacing!(fig, plot_layout)
        # Note: grid_padding is stored but not yet applied (future feature)
        
        for (idx, channel) in enumerate(plot_layout.channels)
            row = fld(idx-1, plot_layout.cols) + 1
            col = mod(idx-1, plot_layout.cols) + 1

            ax = _create_axis_with_grid(fig, fig[row, col]; kwargs...)
            _apply_aspect_ratio!(ax, plot_layout)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end
        
        # Apply grid padding if specified (this affects the content area)
        # Note: grid_padding is less commonly used, but we store it for potential future use

    elseif plot_layout.type == :topo
        plot_width = plot_layout.metadata[:plot_width]
        plot_height = plot_layout.metadata[:plot_height]
        margin = plot_layout.metadata[:margin]

        # Use the coordinate normalization helper
        x_coords, y_coords, minx, maxx, miny, maxy, xrange, yrange =
            _normalize_coordinates(plot_layout.positions, margin)

        # For topographic layout, create individual axes positioned using halign/valign with Relative sizing
        for (idx, (channel, pos)) in enumerate(zip(plot_layout.channels, plot_layout.positions))
            # Normalize position to [0,1]
            norm_x = (pos[1] - minx) / xrange
            norm_y = (pos[2] - miny) / yrange

            # Clamp positions to margin bounds
            halign = clamp(norm_x, margin, 1 - margin)
            valign = clamp(norm_y, margin, 1 - margin)

            # Create axis with relative positioning and grid settings
            ax = _create_topo_axis(fig, plot_width, plot_height, halign, valign; kwargs...)
            _apply_aspect_ratio!(ax, plot_layout)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end

    elseif plot_layout.type == :custom
        for (idx, (channel, pos)) in enumerate(zip(plot_layout.channels, plot_layout.positions))
            ax = _create_axis_with_grid(fig, fig[pos[1], pos[2]]; kwargs...)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end
    end

    return axes, channels
end

"""
    _set_grid_axis_properties!(ax::Axis, plot_layout::PlotLayout, channel::Symbol, 
                              row::Int, col::Int, total_rows::Int, total_cols::Int; 
                              kwargs...)

Set properties for axes in a grid layout.
"""
function _set_grid_axis_properties!(
    ax::Axis,
    plot_layout::PlotLayout,
    channel::Symbol,
    row::Int,
    col::Int,
    total_rows::Int,
    total_cols::Int;
    kwargs...,
)

    # Use channel name unless user explicitly specified a title
    # For grid layouts, this means each subplot shows its channel name
    # unless a global title is specified
    # Note: DEFAULT_ERP_KWARGS sets title="", so we need to check if user actually specified something
    if haskey(kwargs, :title) && kwargs[:title] != ""
        ax.title = kwargs[:title]
    else
        ax.title = string(channel)
    end

    # Only show y-axis labels on the leftmost column
    if col != 1
        ax.ylabel = ""
        ax.yticklabelsvisible = false
    end

    # Only show x-axis labels on the bottom row
    if row != total_rows
        ax.xlabel = ""
        ax.xticklabelsvisible = false
    end

end



"""
    _apply_layout_axis_properties!(axes::Vector{Axis}, plot_layout::PlotLayout; kwargs...)

Apply layout-specific axis properties to all axes after plotting.
This ensures layout properties override any auto-set properties from Makie.
"""
function _apply_layout_axis_properties!(axes::Vector{Axis}, plot_layout::PlotLayout; kwargs...)
    if plot_layout.type == :grid
        for (idx, ax) in enumerate(axes)
            row = fld(idx-1, plot_layout.cols) + 1
            col = mod(idx-1, plot_layout.cols) + 1
            _set_grid_axis_properties!(
                ax,
                plot_layout,
                plot_layout.channels[idx],
                row,
                col,
                plot_layout.rows,
                plot_layout.cols;
                kwargs...,
            )
        end
    elseif plot_layout.type == :topo
        # For topo layouts, remove axis labels and ticks, and add scale plot
        for (idx, ax) in enumerate(axes)
            # Set the title to show the channel name
            ax.title = string(plot_layout.channels[idx])
            hidedecorations!(ax, grid = false, ticks = true, ticklabels = true)
            hidespines!(ax)
        end
    end
end

"""
    _apply_axis_properties!(ax::Axis; kwargs...)

Apply common axis properties from keyword arguments.
"""
function _apply_axis_properties!(ax::Axis; kwargs...)

    ax.title = kwargs[:title]
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]

    # Axis limits and y +ve/ve inversions
    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    ax.yreversed = kwargs[:yreversed]

    # Apply grid settings
    ax.xgridvisible = kwargs[:xgrid]
    ax.ygridvisible = kwargs[:ygrid]
    ax.xminorgridvisible = kwargs[:xminorgrid]
    ax.yminorgridvisible = kwargs[:yminorgrid]

end


"""
    _add_scale_axis!(fig::Figure, plot_layout::PlotLayout; 
                    scale_axis_x::Float64 = 0.95,
                    scale_axis_y::Float64 = 0.05,
                    scale_axis_width::Float64 = 0.10,
                    scale_axis_height::Float64 = 0.10)

Add a scale axis to show the scale of the plots.
"""
function _add_scale_axis!(
    fig::Figure,
    plot_layout::PlotLayout;
    scale_axis_x::Float64 = 0.95,
    scale_axis_y::Float64 = 0.05,
    scale_axis_width::Float64 = 0.10,
    scale_axis_height::Float64 = 0.10,
)

    # Create scale axis
    scale_ax = Axis(fig[scale_axis_x, scale_axis_y], width = scale_axis_width, height = scale_axis_height)

    # Hide decorations
    scale_ax.xlabel = ""
    scale_ax.ylabel = ""
    scale_ax.xticklabelsvisible = false
    scale_ax.yticklabelsvisible = false

    # Add scale information
    scale_ax.title = "Scale"

    return scale_ax
end
