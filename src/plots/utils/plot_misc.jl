
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
