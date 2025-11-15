
function display_figure(fig)
    display(get_makie_screen(get_makie_backend()), fig)
end

function get_makie_backend()
    backend_str = string(Makie.current_backend())
    if occursin("GLMakie", backend_str)
        return :GLMakie
    elseif occursin("CairoMakie", backend_str)
        return :CairoMakie
    else
        @minimal_error "Makie backend not found"
    end
end

function get_makie_screen(makie_backend::Symbol)
    return getfield(Main, makie_backend).Screen()
end

"""
    _get_colorbar_defaults()

Get all default values for Colorbar attributes by creating a single Colorbar instance.
Returns a dictionary mapping attribute names to their default values.
"""
function _get_colorbar_defaults()
    # Create a minimal figure
    fig = Figure()
    cb = Colorbar(fig)

    # Get all attribute values at once
    defaults = Dict{Symbol,Any}()
    for attr in propertynames(Colorbar)
        defaults[attr] = getproperty(cb, attr)
    end

    return defaults
end

# Cache the colorbar defaults
const COLORBAR_DEFAULTS = _get_colorbar_defaults()

"""
    _get_legend_defaults()

Get all default values for Legend attributes by creating a temporary Legend instance.
Returns a dictionary mapping attribute names to their default values.
"""
function _get_legend_defaults()
    # Create a minimal figure and axis with a dummy plot that has labels
    # This is necessary because axislegend() requires plots with labels
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, [1, 2], [1, 2], label = "dummy")  # Create a plot with a label
    leg = axislegend(ax)

    # Get all attribute values at once
    defaults = Dict{Symbol,Any}()
    for attr in propertynames(Legend)
        defaults[attr] = getproperty(leg, attr)
    end

    return defaults
end

# Cache the legend defaults
const LEGEND_DEFAULTS = _get_legend_defaults()

"""
    _extract_colorbar_kwargs!(plot_kwargs::Dict{Symbol, Any})

Extract all colorbar-related parameters from plot_kwargs and return a clean dictionary
suitable for passing to Colorbar constructor.

# Arguments
- `plot_kwargs`: Dictionary of plot parameters (modified in-place)

# Returns
- `Dict{Symbol, Any}`: Cleaned colorbar parameters with invalid attributes removed
"""
function _extract_colorbar_kwargs!(plot_kwargs::Dict{Symbol,Any})
    colorbar_kwargs = Dict{Symbol,Any}()
    colorbar_attrs = propertynames(Colorbar)

    for attr in colorbar_attrs
        colorbar_key = Symbol("colorbar_$(attr)")
        if haskey(plot_kwargs, colorbar_key)
            value = pop!(plot_kwargs, colorbar_key)
            if value !== nothing  # Only add if not the default nothing
                colorbar_kwargs[attr] = value
            end
        end
    end

    # These cannot be passed to colorbar kwargs
    pop!(colorbar_kwargs, :colormap, nothing)
    pop!(colorbar_kwargs, :limits, nothing)
    pop!(colorbar_kwargs, :highclip, nothing)
    pop!(colorbar_kwargs, :lowclip, nothing)

    return colorbar_kwargs
end

"""
    _extract_legend_kwargs(plot_kwargs::Dict{Symbol, Any}; exclude_positioning::Bool=false)

Extract legend-related parameters from plot_kwargs and return a new dictionary
suitable for passing to axislegend().

Does not mutate plot_kwargs - only reads from it.

# Arguments
- `plot_kwargs`: Dictionary of plot parameters (read-only)
- `exclude_positioning`: If true, exclude positioning attributes (halign, valign, alignmode) 
  that conflict with the `position` parameter in axislegend()

# Returns
- `Dict{Symbol, Any}`: New dictionary containing extracted legend parameters (with legend_ prefix removed)
"""
function _extract_legend_kwargs(plot_kwargs::Dict{Symbol,Any}; exclude_positioning::Bool = false)
    legend_kwargs = Dict{Symbol,Any}()
    legend_attrs = propertynames(Legend)

    for attr in legend_attrs
        # attr in positioning_attrs && continue  # Skip positioning attributes if requested
        legend_key = Symbol("legend_$(attr)")
        if haskey(plot_kwargs, legend_key)
            value = plot_kwargs[legend_key]  # Read from plot_kwargs
            if value in [:legend_halign, :legend_valign, :legend_alignmode]
                println("Skipping positioning attribute: $attr")
                continue
            end
            if value !== nothing  # Only add if not the default nothing
                legend_kwargs[attr] = value  # Store in new dict without "legend_" prefix
            end
        end
    end
  
    # TODO: I do not know why this is needed? Bug here? Bug Makie axislegend?
    # Without it legend_position is ignored!
    # Remove positioning attributes that conflict with explicit position parameter
    [pop!(legend_kwargs, attr, nothing) for attr in [:halign, :valign, :alignmode]]

    return legend_kwargs
end

# =============================================================================
# AXIS STYLING FUNCTIONS
# =============================================================================

"""
    _set_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)

Apply grid settings to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `xgrid`: Whether to show x-axis grid
- `ygrid`: Whether to show y-axis grid  
- `xminorgrid`: Whether to show x-axis minor grid
- `yminorgrid`: Whether to show y-axis minor grid
"""
function _set_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)
    ax.xgridvisible = xgrid
    ax.ygridvisible = ygrid
    ax.xminorgridvisible = xminorgrid
    ax.yminorgridvisible = yminorgrid
end

"""
    _set_axis_properties!(ax; xlim = nothing, ylim = nothing, xlabel = "", ylabel = "", yreversed = false)

Apply axis limits, labels, and direction to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `xlim`: X-axis limits as (min, max) tuple or nothing for auto-scaling
- `ylim`: Y-axis limits as (min, max) tuple or nothing for auto-scaling
- `xlabel`: Label for x-axis (default: empty string)
- `ylabel`: Label for y-axis (default: empty string)
- `yreversed`: Whether to reverse the y-axis (default: false)
"""
function _set_axis_properties!(ax; xlim = nothing, ylim = nothing, xlabel = "", ylabel = "", yreversed = false)

    # Set axis labels
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.yreversed = yreversed
    
    # Set axis limits
    if xlim !== nothing
        xmin, xmax = xlim
        xlims!(ax, xmin, xmax)
    end
    
    if ylim !== nothing
        ymin, ymax = ylim
        ylims!(ax, ymin, ymax)
    end
end

"""
    _set_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 0.5, alpha = 0.7)

Add origin lines at x=0 and y=0 to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `add_xy_origin`: Whether to add origin lines at x=0 and y=0
- `color`: Color of the origin lines
- `linewidth`: Line width of the origin lines
- `alpha`: Transparency of the origin lines
"""
function _set_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 0.5, alpha = 0.7)
    if add_xy_origin
        hlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
        vlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
    end
end
