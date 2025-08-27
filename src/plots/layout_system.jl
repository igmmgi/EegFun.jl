"""
    PlotLayout

A struct that defines how to arrange plots in a figure.
"""
struct PlotLayout
    type::Symbol  # :single, :grid, :topo, :custom
    rows::Int     
    cols::Int     
    positions::Vector{Tuple{Float64, Float64}}  # Custom positions for :custom type
    channels::Vector{Symbol}     # Channels to plot
    metadata::Dict{Symbol, Any}  # Additional layout-specific parameters
end



"""
    create_single_layout(channels::Vector{Symbol})

Create a single plot layout for multiple channels (all shown on one plot).
"""
function create_single_layout(channels::Vector{Symbol})
    return PlotLayout(:single, 1, 1, [(0.0, 0.0)], channels, Dict{Symbol, Any}())
end

"""
    create_grid_layout(channels::Vector{Symbol}; rows::Int = 0, cols::Int = 0)

Create a grid layout for multiple channels.
If rows or cols is 0, it will be calculated automatically.
"""
function create_grid_layout(channels::Vector{Symbol}; rows::Int = 0, cols::Int = 0)
    n_channels = length(channels)
    if rows == 0 && cols == 0
        rows, cols = best_rect(n_channels)
    elseif rows == 0
        rows = ceil(Int, n_channels / cols)
    elseif cols == 0
        cols = ceil(Int, n_channels / rows)
    end
    
    if rows * cols < n_channels
        throw(ArgumentError("Grid ($(rows)×$(cols)=$(rows*cols)) is too small for $n_channels channels"))
    end
    
    return PlotLayout(:grid, rows, cols, [], channels, Dict{Symbol, Any}())
end

"""
    create_topo_layout(layout::Layout, channels::Vector{Symbol}; 
                       plot_width::Float64 = 0.10, 
                       plot_height::Float64 = 0.10, 
                       margin::Float64 = 0.02,
                       scale_offset_factor::Float64 = 0.1,
                       fallback_scale_x::Float64 = 0.8,
                       fallback_scale_y::Float64 = -0.8)

Create a topographic layout based on channel positions.
"""
function create_topo_layout(layout::Layout, channels::Vector{Symbol}; 
                           plot_width::Float64 = 0.10, 
                           plot_height::Float64 = 0.10, 
                           margin::Float64 = 0.02,
                           scale_offset_factor::Float64 = 0.1,
                           fallback_scale_x::Float64 = 0.8,
                           fallback_scale_y::Float64 = -0.8)
   
    _ensure_coordinates_2d!(layout)
    
    # Get channel positions from layout data
    positions = []
    for channel in channels
        idx = findfirst(==(channel), layout.data.label)
        if idx !== nothing
            x = layout.data.x2[idx]
            y = layout.data.y2[idx]
            push!(positions, (x, y))
        else
            @minimal_warning "Channel $channel not found in layout, using default position"
            push!(positions, (0.0, 0.0))
        end
    end
    
    # Add scale plot position (bottom right area) using existing positions
    if !isempty(positions)
        x_coords = [pos[1] for pos in positions]
        y_coords = [pos[2] for pos in positions]
        minx, maxx = extrema(x_coords)
        miny, maxy = extrema(y_coords)
        
        # Place scale plot in bottom right area
        scale_x = maxx + (maxx - minx) * scale_offset_factor  # Slightly to the right
        scale_y = miny - (maxy - miny) * scale_offset_factor  # Slightly below
        scale_position = (scale_x, scale_y)
    else
        scale_position = (fallback_scale_x, fallback_scale_y)
    end
    
    # Store metadata for topographic layout
    metadata = Dict{Symbol, Any}(
        :plot_width => plot_width,
        :plot_height => plot_height,
        :margin => margin,
        :scale_offset_factor => scale_offset_factor,
        :fallback_scale_x => fallback_scale_x,
        :fallback_scale_y => fallback_scale_y
    )
    
    return PlotLayout(:topo, 0, 0, positions, channels, metadata)
end

"""
    _normalize_coordinates(positions::Vector{Tuple{Float64, Float64}}, margin::Float64)

Normalize coordinates to [0,1] range and handle edge cases.
Returns normalized coordinates and ranges for layout positioning.
"""
function _normalize_coordinates(positions::Vector{Tuple{Float64, Float64}}, margin::Float64)
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
    create_layout(layout_spec, channels, eeg_layout)

Create a PlotLayout object based on the layout specification.
This is a generic function that can be used by any plot type.
"""
function create_layout(layout_spec::Union{Symbol, PlotLayout, Vector{Int}}, 
                       channels::Vector{Symbol}, 
                       eeg_layout::Union{Layout, Nothing})
    if layout_spec === :single
        return create_single_layout(channels)
    elseif layout_spec === :grid
        return create_grid_layout(channels)
    elseif layout_spec === :topo
        return create_topo_layout(eeg_layout, channels)
    elseif layout_spec isa Vector{Int}
        if length(layout_spec) != 2
            throw(ArgumentError("layout must be a 2-element vector [rows, cols]"))
        end
        return create_grid_layout(channels, rows = layout_spec[1], cols = layout_spec[2])
    elseif layout_spec isa PlotLayout
        return layout_spec
    else
        throw(ArgumentError("Invalid layout specification: $layout_spec"))
    end
end

"""
    create_custom_layout(positions::Vector{Tuple{Float64, Float64}}, 
                        channels::Vector{Symbol})

Create a custom layout with specific positions for each channel.
"""
function create_custom_layout(positions::Vector{Tuple{Float64, Float64}}, channels::Vector{Symbol})
    if length(positions) != length(channels)
        throw(ArgumentError("Number of positions ($(length(positions))) must match number of channels ($(length(channels)))"))
    end
    return PlotLayout(:custom, 0, 0, positions, channels, Dict{Symbol, Any}())
end

"""
    best_rect(n::Int)

Find the best rectangular layout for n items.
Returns (rows, cols) that minimizes the difference between rows and cols.
"""
function best_rect(n::Int)
    if n <= 0
        throw(ArgumentError("n must be positive"))
    end
    
    if n == 1
        return (1, 1)
    end
    
    # Find factors of n
    factors = []
    for i in 1:isqrt(n)
        if n % i == 0
            push!(factors, (i, n ÷ i))
        end
    end
    
    if isempty(factors)
        # n is prime, find closest rectangular arrangement
        rows = ceil(Int, sqrt(n))
        cols = ceil(Int, n / rows)
        return (rows, cols)
    end
    
    # Return the factor pair with the smallest difference
    return argmin(factors) do (r, c)
        abs(r - c)
    end
end

"""
    apply_layout!(fig::Figure, plot_layout::PlotLayout, plot_function::Function, 
                 data, args...; kwargs...)

Apply a plot layout to a figure and call the plot function for each channel.
Returns the created axes.
"""
function apply_layout!(fig::Figure, plot_layout::PlotLayout, plot_function::Function, data, args...; kwargs...)
    
    axes = Axis[]
    
    if plot_layout.type == :single
        ax = Axis(fig[1, 1])
        push!(axes, ax)
        
        # Call plot function for the single channel
        plot_function(ax, data, plot_layout.channels, args...; kwargs...)
        
    elseif plot_layout.type == :grid
        for (idx, channel) in enumerate(plot_layout.channels)
            row = fld(idx-1, plot_layout.cols) + 1
            col = mod(idx-1, plot_layout.cols) + 1
            
            ax = Axis(fig[row, col])
            push!(axes, ax)
            
            # Set grid-specific axis properties (clean labels)
            _set_grid_axis_properties!(ax, plot_layout, channel, row, col, plot_layout.rows, plot_layout.cols; kwargs...)
            
            # Call plot function for this channel
            plot_function(ax, data, [channel], args...; kwargs...)
        end
        
    elseif plot_layout.type == :topo
        plot_width = get(plot_layout.metadata, :plot_width, 0.10)
        plot_height = get(plot_layout.metadata, :plot_height, 0.10)
        margin = get(plot_layout.metadata, :margin, 0.02)
        
        # Use the new coordinate normalization helper
        x_coords, y_coords, minx, maxx, miny, maxy, xrange, yrange = _normalize_coordinates(plot_layout.positions, margin)
        
        # For topographic layout, use the same approach as plot_epochs
        # Create individual axes positioned using halign/valign with Relative sizing
        for (idx, (channel, pos)) in enumerate(zip(plot_layout.channels, plot_layout.positions))
            # Normalize position to [0,1]
            norm_x = (pos[1] - minx) / xrange
            norm_y = (pos[2] - miny) / yrange
            
            # Clamp positions to margin bounds
            halign = clamp(norm_x, margin, 1 - margin)
            valign = clamp(norm_y, margin, 1 - margin)
            
            # Create axis with relative positioning
            ax = Axis(
                fig[1, 1],
                width = Relative(plot_width),
                height = Relative(plot_height),
                halign = halign,
                valign = valign
            )
            push!(axes, ax)
            plot_function(ax, data, [channel], args...; kwargs...)
        end
        
    elseif plot_layout.type == :custom
        for (channel, pos) in zip(plot_layout.channels, plot_layout.positions)
            ax = Axis(fig[pos[1], pos[2]])
            push!(axes, ax)
            plot_function(ax, data, [channel], args...; kwargs...)
        end
    end
    
    return axes
end

"""
    _set_grid_axis_properties!(ax::Axis, plot_layout::PlotLayout, channel::Symbol, 
                              row::Int, col::Int, total_rows::Int, total_cols::Int; 
                              kwargs...)

Set properties for axes in a grid layout.
"""
function _set_grid_axis_properties!(ax::Axis, plot_layout::PlotLayout, channel::Symbol, 
                                  row::Int, col::Int, total_rows::Int, total_cols::Int; 
                                  kwargs...)
    
    ax.title = string(channel)
    
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
    
    # Apply other properties
    if haskey(kwargs, :ylim)
        ylims!(ax, kwargs[:ylim])
    end
    
    if haskey(kwargs, :xlim)
        xlims!(ax, kwargs[:xlim])
    end
end



"""
    _add_scale_axis!(fig::Figure, plot_layout::PlotLayout; 
                    scale_axis_x::Float64 = 0.95,
                    scale_axis_y::Float64 = 0.05,
                    scale_axis_width::Float64 = 0.10,
                    scale_axis_height::Float64 = 0.10)

Add a scale axis to show the scale of the plots.
"""
function _add_scale_axis!(fig::Figure, plot_layout::PlotLayout; 
                         scale_axis_x::Float64 = 0.95,
                         scale_axis_y::Float64 = 0.05,
                         scale_axis_width::Float64 = 0.10,
                         scale_axis_height::Float64 = 0.10)
    
    # Create scale axis
    scale_ax = Axis(fig[scale_axis_x, scale_axis_y], 
                    width = scale_axis_width, height = scale_axis_height)
    
    # Hide decorations
    scale_ax.xlabel = ""
    scale_ax.ylabel = ""
    scale_ax.xticklabelsvisible = false
    scale_ax.yticklabelsvisible = false
    
    # Add scale information
    scale_ax.title = "Scale"
    
    return scale_ax
end

"""
    _apply_axis_properties!(ax::Axis; kwargs...)

Apply common axis properties from keyword arguments.
"""
function _apply_axis_properties!(ax::Axis; kwargs...)
    
    # Apply limits
    if haskey(kwargs, :xlim) && !isnothing(kwargs[:xlim])
        xlims!(ax, kwargs[:xlim])
    end
    
    if haskey(kwargs, :ylim) && !isnothing(kwargs[:ylim])
        ylims!(ax, kwargs[:ylim])
    end
    
    # Apply other properties
    if haskey(kwargs, :yreversed) && kwargs[:yreversed]
        ax.yreversed = true
    end
    
    # Apply fontsize directly to the axis (avoiding Makie theming issues)
    if haskey(kwargs, :theme_fontsize)
        try
            ax.fontsize = kwargs[:theme_fontsize]
        catch
            @minimal_warning "Could not set fontsize directly on axis, skipping fontsize setting"
        end
    end
end
