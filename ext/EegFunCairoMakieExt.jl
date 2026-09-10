module EegFunCairoMakieExt

using EegFun
using CairoMakie

function __init__()
    EegFun.MAKIE_EXT_STATE.cairomakie_active = true
    EegFun.MAKIE_EXT_STATE.save_vector = (path, fig; pt_per_unit = 1.0) -> begin
        Makie.save(path, fig; backend = CairoMakie, pt_per_unit = pt_per_unit)
    end
end

end # module
