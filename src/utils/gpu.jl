# =============================================================================
# UNIFIED GPU INTERFACE AND BACKEND REGISTRY
# =============================================================================
# EegFun uses KernelAbstractions as its primary GPU kernel abstraction engine.
# Vendor-specific packages (AMDGPU.jl, CUDA.jl, Metal.jl) are handled
# via Julia 1.9+ package extensions. When a user loads a GPU package (e.g., `using CUDA`),
# the corresponding extension registers its array constructor and KA backend with EegFun.

mutable struct GpuState
    active::Bool
    backend_name::Symbol
    device_name::String
    backend_obj::Any
    converter::Any
end

const GLOBAL_GPU_STATE = GpuState(false, :CPU, "CPU", CPU(), identity)

"""
    register_gpu_backend!(backend_name::Symbol, device_name::String, backend_obj, converter)

Internal registration function called by vendor package extensions (`ext/EegFun...Ext.jl`).
"""
function register_gpu_backend!(backend_name::Symbol, device_name::String, backend_obj, converter)
    GLOBAL_GPU_STATE.active = true
    GLOBAL_GPU_STATE.backend_name = backend_name
    GLOBAL_GPU_STATE.device_name = device_name
    GLOBAL_GPU_STATE.backend_obj = backend_obj
    GLOBAL_GPU_STATE.converter = converter
    @info "[GPU BACKEND REGISTERED] $backend_name ($device_name)"
    return nothing
end

"""
    is_gpu_available() -> Bool

Check whether a functional GPU extension (AMDGPU, CUDA, Metal) has been loaded and registered.
"""
function is_gpu_available()
    return GLOBAL_GPU_STATE.active
end

"""
    gpu_backend() -> KernelAbstractions.Backend

Get the active registered KernelAbstractions backend object (or KernelAbstractions.CPU() by default).
"""
function gpu_backend()
    return GLOBAL_GPU_STATE.backend_obj
end

"""
    gpu_array(x::AbstractArray)

Convert a CPU Array to the active GPU array type (or return `x` unchanged if no GPU extension is registered).
"""
function gpu_array(x::AbstractArray)
    if GLOBAL_GPU_STATE.active && GLOBAL_GPU_STATE.converter !== identity
        return GLOBAL_GPU_STATE.converter(x)
    end
    return x
end

"""
    gpu_device_name() -> String

Return the device name of the active GPU backend.
"""
function gpu_device_name()
    return GLOBAL_GPU_STATE.device_name
end
