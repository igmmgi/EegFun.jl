module EegFunCUDAExt

using EegFun
using CUDA

function __init__()
    dev_name = try
        CUDA.name(CUDA.device())
    catch
        "NVIDIA GPU"
    end
    EegFun.register_gpu_backend!(:CUDA, dev_name, CUDA.CUDABackend(), CUDA.CuArray)
end

end # module
