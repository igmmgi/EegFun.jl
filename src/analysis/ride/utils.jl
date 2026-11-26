# Signal shifting, Tukey tapering, and cross-correlation helpers for RIDE.

######################################
########### User Functions ###########
######################################

"""
    shift_signal(src::AbstractVector{T}, shift::Int) where {T} -> Vector{T}

Shift `src` by `shift` samples into a newly allocated vector, zero-padding vacated elements (`shift > 0` shifts right, `shift < 0` shifts left).
"""
function shift_signal(src::AbstractVector{T}, shift::Int) where {T}
    dest = similar(src)
    return shift_signal!(dest, src, shift)
end

"""
    shift_signal!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T}

In-place version of [`shift_signal`](@ref) writing into `dest`.
"""
function shift_signal!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T}
    n = length(src)
    if abs(shift) >= n
        fill!(dest, zero(T))
    elseif shift > 0
        copyto!(dest, shift + 1, src, 1, n - shift)
        fill!(view(dest, 1:shift), zero(T))
    elseif shift < 0
        s = -shift
        copyto!(dest, 1, src, s + 1, n - s)
        fill!(view(dest, (n-s+1):n), zero(T))
    else
        copyto!(dest, src)
    end
    return dest
end

"""
    shift_matrix(src::AbstractMatrix{T}, shift::Int) where {T} -> Matrix{T}

Shift each column of `src` (`samples × channels`) by `shift` samples into a newly allocated matrix, zero-padding vacated elements.
"""
function shift_matrix(src::AbstractMatrix{T}, shift::Int) where {T}
    dest = similar(src)
    return shift_matrix!(dest, src, shift)
end

"""
    shift_matrix!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}, shift::Int) where {T}

In-place version of [`shift_matrix`](@ref) writing into `dest`.
"""
function shift_matrix!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}, shift::Int) where {T}
    @inbounds for col in axes(src, 2)
        shift_signal!(view(dest, :, col), view(src, :, col), shift)
    end
    return dest
end

######################################
######### Internal Functions #########
######################################

@inline function _accum_shifted!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int, ::Val{Sign}) where {T,Sign}
    n = length(dest)
    abs(shift) >= n && return dest
    s_sgn = T(Sign)
    if shift > 0
        @inbounds @simd for s = 1:(n-shift)
            dest[s+shift] += s_sgn * src[s]
        end
    elseif shift < 0
        s_val = -shift
        @inbounds @simd for s = 1:(n-s_val)
            dest[s] += s_sgn * src[s+s_val]
        end
    else
        @inbounds @simd for s = 1:n
            dest[s] += s_sgn * src[s]
        end
    end
    return dest
end

@inline _subtract_shifted!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T} =
    _accum_shifted!(dest, src, shift, Val(-1))

@inline _add_shifted!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T} = _accum_shifted!(dest, src, shift, Val(1))
