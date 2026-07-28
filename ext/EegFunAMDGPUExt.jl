module EegFunAMDGPUExt

using EegFun
using AMDGPU

function __init__()
    dev_name = try
        string(AMDGPU.device())
    catch
        "AMD GPU"
    end
    EegFun.register_gpu_backend!(:AMDGPU, dev_name, AMDGPU.ROCBackend(), AMDGPU.ROCArray)
end

end # module
