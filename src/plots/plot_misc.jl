
function display_figure(fig)
    display(get_makie_screen(get_makie_backend()), fig)
end

function get_makie_backend()
    backend_str = string(Makie.current_backend())
    if occursin("WGLMakie", backend_str)
        return :WGLMakie
    elseif occursin("GLMakie", backend_str)
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
