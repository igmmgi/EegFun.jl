# CUDA.jl backend extension for EegFun GPU acceleration.

module EegFunCUDAExt

using EegFun

using CUDA

######################################
######### Internal Functions #########
######################################

function __init__()
    if CUDA.functional()
        dev_name = try
            CUDA.name(CUDA.device())
        catch
            "NVIDIA GPU"
        end
        EegFun._register_gpu_backend!(:CUDA, dev_name, CUDA.CUDABackend(), CUDA.CuArray)
    end
end

end # module
