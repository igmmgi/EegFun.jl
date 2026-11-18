mutable struct GpuState
    active::Bool
    supports_fft::Bool
    backend_name::Symbol
    device_name::String
    backend_obj::Any
    converter::Any
end

const GLOBAL_GPU_STATE = GpuState(false, false, :CPU, "CPU", CPU(), identity)

function register_gpu_backend!(backend_name::Symbol, device_name::String, backend_obj, converter; supports_fft::Bool = true)
    GLOBAL_GPU_STATE.active = true
    GLOBAL_GPU_STATE.supports_fft = supports_fft
    GLOBAL_GPU_STATE.backend_name = backend_name
    GLOBAL_GPU_STATE.device_name = device_name
    GLOBAL_GPU_STATE.backend_obj = backend_obj
    GLOBAL_GPU_STATE.converter = converter
    @info "[GPU BACKEND REGISTERED] $backend_name ($device_name)"
    return nothing
end

is_gpu_available() = GLOBAL_GPU_STATE.active
is_gpu_fft_available() = GLOBAL_GPU_STATE.active && GLOBAL_GPU_STATE.supports_fft
gpu_backend() = GLOBAL_GPU_STATE.backend_obj
gpu_device_name() = GLOBAL_GPU_STATE.device_name

function gpu_array(x::AbstractArray)
    if GLOBAL_GPU_STATE.active && GLOBAL_GPU_STATE.converter !== identity
        return GLOBAL_GPU_STATE.converter(x)
    end
    return x
end
