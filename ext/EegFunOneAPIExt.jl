module EegFunOneAPIExt

using EegFun
using oneAPI

function __init__()
    dev_name = try
        string(oneAPI.device())
    catch
        "Intel GPU"
    end
    EegFun.register_gpu_backend!(:oneAPI, dev_name, oneAPI.oneAPIBackend(), oneAPI.oneArray)
end

end # module
