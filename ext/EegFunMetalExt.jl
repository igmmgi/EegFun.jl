module EegFunMetalExt

using EegFun
using Metal

function __init__()
    dev_name = try
        Metal.name(Metal.current_device())
    catch
        "Apple Silicon GPU"
    end
    EegFun.register_gpu_backend!(:Metal, dev_name, Metal.MetalBackend(), Metal.MtlArray)
end

end # module
