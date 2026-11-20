# GPU backend registration, status queries, and device array conversion hooks.

######################################
######### Types & Constants ##########
######################################

mutable struct GpuState
    active::Bool
    supports_fft::Bool
    backend_name::Symbol
    device_name::String
    backend_obj::Any
    converter::Any
end

const GLOBAL_GPU_STATE = GpuState(false, false, :CPU, "CPU", CPU(), identity)

######################################
########### User Functions ###########
######################################

"""
    is_gpu_available() -> Bool

Return `true` when an accelerated GPU backend (`CUDA`, `AMDGPU`, or `Metal`) is active.
"""
is_gpu_available() = GLOBAL_GPU_STATE.active

"""
    is_gpu_fft_available() -> Bool

Return `true` when an active GPU backend supports device FFT plans.
"""
is_gpu_fft_available() = GLOBAL_GPU_STATE.active && GLOBAL_GPU_STATE.supports_fft

"""
    gpu_backend()

Return the active `KernelAbstractions` backend (`CPU()` when no GPU extension is loaded).
"""
gpu_backend() = GLOBAL_GPU_STATE.backend_obj

"""
    gpu_device_name() -> String

Return the human-readable name of the active compute device (`"CPU"` by default).
"""
gpu_device_name() = GLOBAL_GPU_STATE.device_name

######################################
######### Internal Functions #########
######################################

function _register_gpu_backend!(backend_name::Symbol, device_name::String, backend_obj, converter; supports_fft::Bool = true)
    GLOBAL_GPU_STATE.active = true
    GLOBAL_GPU_STATE.supports_fft = supports_fft
    GLOBAL_GPU_STATE.backend_name = backend_name
    GLOBAL_GPU_STATE.device_name = device_name
    GLOBAL_GPU_STATE.backend_obj = backend_obj
    GLOBAL_GPU_STATE.converter = converter
    @info "[GPU BACKEND REGISTERED] $backend_name ($device_name)"
    return nothing
end

function _gpu_array(x::AbstractArray)
    if GLOBAL_GPU_STATE.active && GLOBAL_GPU_STATE.converter !== identity
        return GLOBAL_GPU_STATE.converter(x)
    end
    return x
end
