"""
    RIDE (Residue Iteration DEcomposition) Utilities

Contains high-performance array shifting and padding functions
crucial for native Julia RIDE iteration.
"""

"""
    shift_signal!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T}

Shift a signal by `shift` samples, zero-padding the vacated elements.
Positive `shift` moves data forward in time (right).
Negative `shift` moves data backward in time (left).
"""
function shift_signal!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T}
    n = length(src)
    
    # If shift is larger than array, return all zeros
    if abs(shift) >= n
        fill!(dest, zero(T))
        return dest
    end
    
    if shift > 0
        # Shift right
        copyto!(dest, shift + 1, src, 1, n - shift)
        fill!(view(dest, 1:shift), zero(T))
    elseif shift < 0
        # Shift left
        s = abs(shift)
        copyto!(dest, 1, src, s + 1, n - s)
        fill!(view(dest, (n - s + 1):n), zero(T))
    else
        # No shift
        copyto!(dest, src)
    end
    return dest
end

"""
    shift_signal(src::AbstractVector{T}, shift::Int) where {T}

Out-of-place version of `shift_signal!`.
"""
function shift_signal(src::AbstractVector{T}, shift::Int) where {T}
    dest = similar(src)
    shift_signal!(dest, src, shift)
    return dest
end

"""
    shift_matrix!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}, shift::Int) where {T}

Shift each column of a matrix (samples x channels) by `shift` samples.
"""
function shift_matrix!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}, shift::Int) where {T}
    @inbounds for col in 1:size(src, 2)
        v_src = view(src, :, col)
        v_dest = view(dest, :, col)
        shift_signal!(v_dest, v_src, shift)
    end
    return dest
end

function shift_matrix(src::AbstractMatrix{T}, shift::Int) where {T}
    dest = similar(src)
    shift_matrix!(dest, src, shift)
    return dest
end
