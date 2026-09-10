module EegFunGLMakieExt

using EegFun
using GLMakie

function __init__()
    EegFun.MAKIE_EXT_STATE.glmakie_active = true
    EegFun.MAKIE_EXT_STATE.create_screen = (; size = nothing) -> begin
        return isnothing(size) ? GLMakie.Screen() : GLMakie.Screen(size = size)
    end
    EegFun.MAKIE_EXT_STATE.display_screen = (fig; size = nothing) -> begin
        screen = isnothing(size) ? GLMakie.Screen() : GLMakie.Screen(size = size)
        display(screen, fig)
        return screen
    end
    EegFun.MAKIE_EXT_STATE.set_title = (title::String) -> begin
        try
            GLMakie.activate!(title = title)
        catch
        end
    end
end

end # module
