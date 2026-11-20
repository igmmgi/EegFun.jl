# Metal.jl backend extension for Apple Silicon GPU acceleration.

module EegFunMetalExt

using EegFun

using Metal

######################################
######### Internal Functions #########
######################################

_to_metal(x::AbstractArray{Float64}) = Metal.MtlArray(Float32.(x))

_to_metal(x::AbstractArray{ComplexF64}) = Metal.MtlArray(ComplexF32.(x))

_to_metal(x::AbstractArray) = Metal.MtlArray(x)

function __init__()
    if Metal.functional()
        dev_name = try
            Metal.name(Metal.current_device())
        catch
            "Apple Silicon GPU"
        end
        EegFun._register_gpu_backend!(:Metal, dev_name, Metal.MetalBackend(), _to_metal; supports_fft = false)
    end
end

end # module
