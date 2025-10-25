
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

# =============================================================================
# AXIS STYLING FUNCTIONS
# =============================================================================

"""
    _setup_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)

Apply grid settings to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `xgrid`: Whether to show x-axis grid
- `ygrid`: Whether to show y-axis grid  
- `xminorgrid`: Whether to show x-axis minor grid
- `yminorgrid`: Whether to show y-axis minor grid
"""
function _setup_axis_grid!(ax; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)
    ax.xgridvisible = xgrid
    ax.ygridvisible = ygrid
    ax.xminorgridvisible = xminorgrid
    ax.yminorgridvisible = yminorgrid
end

"""
    _setup_axis_limits!(ax; xlim = nothing, ylim = nothing)

Apply axis limits to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `xlim`: X-axis limits as (min, max) tuple or nothing for auto-scaling
- `ylim`: Y-axis limits as (min, max) tuple or nothing for auto-scaling
"""
function _setup_axis_limits!(ax; xlim = nothing, ylim = nothing)
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
    _setup_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 0.5, alpha = 0.7)

Add origin lines at x=0 and y=0 to the axis.

# Arguments
- `ax`: Makie Axis object

# Keyword Arguments
- `add_xy_origin`: Whether to add origin lines at x=0 and y=0
- `color`: Color of the origin lines
- `linewidth`: Line width of the origin lines
- `alpha`: Transparency of the origin lines
"""
function _setup_origin_lines!(ax; add_xy_origin = true, color = :gray, linewidth = 0.5, alpha = 0.7)
    if add_xy_origin
        hlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
        vlines!(ax, 0, color = color, linewidth = linewidth, alpha = alpha)
    end
end
