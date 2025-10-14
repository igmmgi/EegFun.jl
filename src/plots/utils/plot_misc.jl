
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
    defaults = Dict{Symbol, Any}()
    for attr in propertynames(Colorbar)
        defaults[attr] = getproperty(cb, attr)
    end
    
    return defaults
end

# Cache the colorbar defaults
const COLORBAR_DEFAULTS = _get_colorbar_defaults()
