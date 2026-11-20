# AMDGPU.jl backend extension for EegFun GPU acceleration.

module EegFunAMDGPUExt

using EegFun

using AMDGPU

######################################
######### Internal Functions #########
######################################

function __init__()
    if AMDGPU.functional()
        dev_name = try
            string(AMDGPU.device())
        catch
            "AMD GPU"
        end
        EegFun._register_gpu_backend!(:AMDGPU, dev_name, AMDGPU.ROCBackend(), AMDGPU.ROCArray)
    end
end

end # module
