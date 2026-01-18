# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const LAYOUT_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Topo layout parameters
    :topo_plot_width => (0.10, "Width of individual plots in topo layout (as fraction of figure width)"),
    :topo_plot_height => (0.10, "Height of individual plots in topo layout (as fraction of figure height)"),
    :topo_margin => (0.02, "Margin between plots in topo layout"),
    :topo_scale_offset => (0.1, "Offset factor for scale plot position"),
    :topo_scale_pos => ((0.8, -0.8), "Fallback position for scale plot as (x, y) tuple"),

    # Grid layout parameters
    :grid_dims => (nothing, "Grid dimensions as (rows, cols) tuple. If nothing, automatically calculated"),
    :grid_rowgap => (10, "Gap between rows in grid layout (in pixels)"),
    :grid_colgap => (10, "Gap between columns in grid layout (in pixels)"),
    :grid_skip_positions =>
        (nothing, "Positions to skip in grid layout as vector of (row, col) tuples, e.g., [(2,1), (2,3)]"),
)


"""
    PlotLayout

A struct that defines how to arrange plots in a figure.
"""
struct PlotLayout
    type::Symbol                               # :single, :grid, :topo, :custom
    dims::AbstractVector{Int}                  # [rows, cols] or (rows, cols) - only meaningful for :grid and :single types
    positions::Vector{Tuple{Float64,Float64}}  # Custom positions for :custom type
    channels::Vector{Symbol}                   # Channels to plot
    metadata::Dict{Symbol,Any}                 # Layout configuration parameters 
end

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


"""
    _create_axis(fig, position; kwargs...)

Create an Axis at the given position with settings from kwargs.
"""
function _create_axis(fig, position; kwargs...)
    return Axis(
        position;
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
Sets gaps for all rows and columns in the grid.
"""
function _apply_grid_spacing!(fig::Figure, plot_layout::PlotLayout)
    rowgap_val = plot_layout.metadata[:grid_rowgap]
    colgap_val = plot_layout.metadata[:grid_colgap]
    rows, cols = plot_layout.dims

    # Set row gaps for all rows (gaps are between rows, so set for rows 1 to rows-1)
    if rowgap_val !== nothing && rows > 1
        for row = 1:(rows-1)
            rowgap!(fig.layout, row, rowgap_val)
        end
    end

    # Set column gaps for all columns (gaps are between columns, so set for cols 1 to cols-1)
    if colgap_val !== nothing && cols > 1
        for col = 1:(cols-1)
            colgap!(fig.layout, col, colgap_val)
        end
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
    _create_single_layout(channels::Vector{Symbol}; kwargs...)

Create a single plot layout for multiple channels (all shown on one plot).

# Keyword Arguments
- `kwargs...`: Additional layout parameters (e.g., `figure_padding`, `aspect_ratio`)
"""
function _create_single_layout(channels::Vector{Symbol}; kwargs...)
    metadata = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)
    return PlotLayout(:single, [1, 1], [(0.0, 0.0)], channels, metadata)
end


"""
    _create_grid_layout(channels::Vector{Symbol}; kwargs...)

Create a grid layout for multiple channels.
Grid dimensions can be specified via `grid_dims` in kwargs (e.g., `grid_dims = (2, 3)`).
If not specified, dimensions are calculated automatically.

# Keyword Arguments
- `kwargs...`: Layout parameters including `grid_dims` (vector/tuple of (rows, cols)), `grid_rowgap`, `grid_colgap`,
"""
function _create_grid_layout(channels::Vector{Symbol}; kwargs...)
    isempty(channels) && throw(ArgumentError("Cannot create grid layout with empty channel list"))
    metadata = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)

    # Get skip positions if provided
    skip_positions = metadata[:grid_skip_positions]
    skip_count = 0
    if skip_positions !== nothing
        # Validate skip positions are tuples
        for pos in skip_positions
            if !isa(pos, Tuple) || length(pos) != 2
                throw(ArgumentError("grid_skip_positions must be a vector of (row, col) tuples, got: $pos"))
            end
        end
        # Convert to Set for faster lookup
        skip_set = Set(skip_positions)
        skip_count = length(skip_set)
        metadata[:grid_skip_positions] = skip_set
    end

    # Calculate grid dimensions
    # If grid_dims is specified, validate it has enough space (accounting for skipped positions)
    # If not specified, calculate dimensions for n_channels + skip_count positions
    n_channels = length(channels)
    if metadata[:grid_dims] !== nothing
        # Validate and auto-correct grid dimensions if too small
        # Need enough positions for n_channels + skip_count (skip positions still need to be in grid)
        rows, cols = _validate_dims(metadata[:grid_dims], n_channels + skip_count)

        # Validate skip positions are within grid bounds
        if skip_positions !== nothing
            for (row, col) in skip_positions
                if row < 1 || row > rows || col < 1 || col > cols
                    throw(ArgumentError("Skip position ($row, $col) is outside grid bounds ($rows×$cols)"))
                end
            end
        end

        # Final check: ensure we have enough available positions after skipping
        available_positions = rows * cols - skip_count
        if available_positions < n_channels
            throw(
                ArgumentError(
                    "Grid size ($rows×$cols) with $skip_count skipped positions has only $available_positions available positions, but need $n_channels",
                ),
            )
        end
    else
        # Auto-calculate dimensions for n_channels + skip_count total positions
        rows, cols = best_rect(n_channels + skip_count)

        # Validate skip positions are within auto-calculated grid bounds
        if skip_positions !== nothing
            for (row, col) in skip_positions
                if row < 1 || row > rows || col < 1 || col > cols
                    throw(
                        ArgumentError(
                            "Skip position ($row, $col) is outside auto-calculated grid bounds ($rows×$cols)",
                        ),
                    )
                end
            end
        end
    end

    return PlotLayout(:grid, [rows, cols], [], channels, metadata)
end

"""
    _create_topo_layout(layout::Layout, channels::Vector{Symbol}; kwargs...)

Create a topographic layout based on channel positions.

# Keyword Arguments
$(generate_kwargs_doc(LAYOUT_KWARGS))
"""
function _create_topo_layout(layout::Layout, channels::Vector{Symbol}; kwargs...)
    # Create metadata from kwargs
    topo_kwargs = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)

    _ensure_coordinates_2d!(layout)

    # Get layout coordinates for bounds calculation 
    layout_x_coords = layout.data.x2
    layout_y_coords = layout.data.y2

    # Calculate bounds from original layout 
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
            push!(positions, (layout.data.x2[idx], layout.data.y2[idx]))
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

    return PlotLayout(:topo, [0, 0], positions, channels, topo_kwargs)
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
    layout_spec::Union{Symbol,PlotLayout},
    channels::Vector{Symbol},
    eeg_layout::Union{Layout,Nothing};
    kwargs...,
)

    if layout_spec === :single
        return _create_single_layout(channels; kwargs...)
    elseif layout_spec === :grid
        return _create_grid_layout(channels; kwargs...)
    elseif layout_spec === :topo
        return _create_topo_layout(eeg_layout, channels; kwargs...)
    elseif layout_spec isa PlotLayout
        return layout_spec
    else
        throw(
            ArgumentError(
                "Invalid layout specification: $layout_spec. Must be :single, :grid, :topo, or a PlotLayout object.",
            ),
        )
    end
end

"""
    create_custom_layout(positions::Vector{Tuple{Float64, Float64}}, channels::Vector{Symbol})

Create a custom layout with specific positions for each channel.
"""
function create_custom_layout(positions::Vector{Tuple{Float64,Float64}}, channels::Vector{Symbol})
    if length(positions) != length(channels)
        throw(ArgumentError("n positions ($(length(positions))) must match n channels ($(length(channels)))"))
    end
    return PlotLayout(:custom, [0, 0], positions, channels, Dict{Symbol,Any}())
end

"""
    _validate_dims(dims, n_items::Int)

Validate dimensions and ensure they fit n_items.
If dims is nothing, returns best_rect(n_items).
If dims is too small for n_items, warns and returns best_rect(n_items).
Otherwise returns dims as-is.

# Arguments
- `dims`: Tuple/Vector of (rows, cols), or nothing
- `n_items::Int`: Number of items to arrange (for validation)

# Returns
- `(rows::Int, cols::Int)`: Valid dimensions (guaranteed to fit n_items)
"""
function _validate_dims(dims, n_items::Int)

    dims === nothing && return best_rect(n_items)

    length(dims) !== 2 && throw(ArgumentError("dims must be a tuple of two positive integers (rows, cols), got: $dims"))

    # If grid is too small, use best_rect instead
    if dims[1] * dims[2] < n_items
        @minimal_warning "Grid size ($(dims[1])×$(dims[2])) is too small for $n_items items. Using optimal grid dimensions instead."
        return best_rect(n_items)
    else
        return dims
    end

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

    # Handle edge cases
    n <= 0 && throw(ArgumentError("n must be positive"))
    n == 1 && return (1, 1)

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

    return (rows, cols)
end

"""
    _apply_layout!(fig::Figure, plot_layout::PlotLayout; kwargs...)

Apply a plot layout to a figure, creating and positioning axes.
Returns the created axes and their associated channels.
The actual plotting should be done separately by the calling function.
"""
function _apply_layout!(fig::Figure, plot_layout::PlotLayout; kwargs...)

    axes = Axis[]
    channels = Symbol[]

    if plot_layout.type == :single
        ax = _create_axis(fig, fig[1, 1]; kwargs...)
        _add_axis_and_channel!(axes, channels, ax, plot_layout.channels[1])

    elseif plot_layout.type == :grid
        rows, cols = plot_layout.dims
        skip_positions = get(plot_layout.metadata, :grid_skip_positions, nothing)

        channel_idx = 1
        for row = 1:rows
            for col = 1:cols
                is_skipped = skip_positions !== nothing && (row, col) in skip_positions
                has_channel = channel_idx <= length(plot_layout.channels)
                if has_channel && !is_skipped
                    channel = plot_layout.channels[channel_idx]
                    ax = _create_axis(fig, fig[row, col]; kwargs...)
                    _add_axis_and_channel!(axes, channels, ax, channel)
                    channel_idx += 1
                else # empty axis
                    ax = _create_axis(fig, fig[row, col]; kwargs...)
                    hidespines!(ax)
                    hidedecorations!(ax)
                end
            end
        end

        # apply grid spacing after axes are created (grid structure is now complete)
        _apply_grid_spacing!(fig, plot_layout)

    elseif plot_layout.type == :topo
        plot_width = plot_layout.metadata[:topo_plot_width]
        plot_height = plot_layout.metadata[:topo_plot_height]

        # Compute bounds for normalization
        positions = plot_layout.positions
        isempty(positions) &&
            throw(ArgumentError("Cannot create topo layout with empty positions. Problem with layout creation!"))

        minx, maxx = extrema(pos[1] for pos in positions)
        miny, maxy = extrema(pos[2] for pos in positions)
        xrange = maxx - minx == 0 ? 1.0 : maxx - minx
        yrange = maxy - miny == 0 ? 1.0 : maxy - miny

        # For topographic layout, create individual axes positioned using halign/valign with Relative sizing
        for (channel, pos) in zip(plot_layout.channels, positions)
            # Normalize to [0,1] - no margin needed since halign/valign centers the plot at the position
            halign = (pos[1] - minx) / xrange
            valign = (pos[2] - miny) / yrange

            # Create axis with relative positioning and grid settings
            ax = _create_topo_axis(fig, plot_width, plot_height, halign, valign; kwargs...)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end

    elseif plot_layout.type == :custom
        for (channel, pos) in zip(plot_layout.channels, plot_layout.positions)
            ax = _create_axis(fig, fig[pos[1], pos[2]]; kwargs...)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end
    end

    return axes, channels
end

"""
    _set_grid_axis_properties!(ax::Axis, channel::Symbol, 
                              row::Int, col::Int, total_rows::Int, total_cols::Int; 
                              kwargs...)

Set properties for axes in a grid layout.
"""
function _set_grid_axis_properties!(
    ax::Axis,
    channel::Symbol,
    row::Int,
    col::Int,
    total_rows::Int,
    total_cols::Int;
    kwargs...,
)

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
        rows, cols = plot_layout.dims
        for (idx, ax) in enumerate(axes)
            row = fld(idx-1, cols) + 1
            col = mod(idx-1, cols) + 1
            _set_grid_axis_properties!(ax, plot_layout.channels[idx], row, col, rows, cols; kwargs...)
        end
    elseif plot_layout.type == :topo
        for (idx, ax) in enumerate(axes)
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
    _add_scale_axis!(fig::Figure; 
                    scale_axis_x::Float64 = 0.95,
                    scale_axis_y::Float64 = 0.05,
                    scale_axis_width::Float64 = 0.10,
                    scale_axis_height::Float64 = 0.10)

Add a scale axis to show the scale of the plots.
"""
function _add_scale_axis!(
    fig::Figure;
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
