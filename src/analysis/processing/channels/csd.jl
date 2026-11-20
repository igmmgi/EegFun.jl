const _CsdData = Union{ContinuousData,ErpData,EpochData}

"""
    compute_csd!(data; m::Int=4, lambda::Real=1e-5, n_legendre_terms::Int=50)

Compute the Current Source Density (CSD / surface Laplacian) of `ContinuousData`, `EpochData`, or `ErpData`
in-place using spherical splines (Perrin et al., 1989; 1990 Corrigenda). Sets `data.analysis_info.reference = :csd`.
"""
function compute_csd!(data::_CsdData; m::Int = 4, lambda::Real = 1e-5, n_legendre_terms::Int = 50)
    channels = channel_labels(data)
    isempty(channels) && @minimal_error("No EEG channels found in data for CSD computation")
    n_ch = length(channels)

    coords = _spherical_coords(data.layout, channels)
    cosang = coords * coords'
    G = _calc_g(cosang, m, n_legendre_terms)
    H = _calc_h(cosang, m, n_legendre_terms)

    if lambda > 0
        lambda_scaled = lambda * tr(G) / n_ch
        @inbounds for i = 1:n_ch
            G[i, i] += lambda_scaled
        end
    end

    C = zeros(Float64, n_ch + 1, n_ch + 1)
    C[1:n_ch, 1:n_ch] .= G
    C[1:n_ch, n_ch+1] .= 1.0
    C[n_ch+1, 1:n_ch] .= 1.0

    C_inv = pinv(C)
    filter_matrix = H * view(C_inv, 1:n_ch, 1:n_ch)

    _apply_csd!(data, filter_matrix, channels)
    data.analysis_info.reference = :csd
    return data
end

function compute_csd!(dat::AbstractVector{<:_CsdData}, args...; kwargs...)
    for d in dat
        compute_csd!(d, args...; kwargs...)
    end
    return dat
end

function _spherical_coords(layout::Layout, channels::AbstractVector{Symbol} = Symbol.(layout.data.label); warn_zero::Bool = true)
    use_polar = hasproperty(layout.data, :inc) && hasproperty(layout.data, :azi)
    !use_polar && _ensure_coordinates_3d!(layout)
    df = layout.data
    label_to_idx = Dict{Symbol,Int}(Symbol(df.label[i]) => i for i = 1:nrow(df))
    coords = zeros(Float64, length(channels), 3)

    for (i, ch) in enumerate(channels)
        row_idx = get(label_to_idx, ch, nothing)
        if isnothing(row_idx)
            @minimal_warning "Channel $ch not found in layout; defaulting unit-sphere coordinate to (0, 0, 1)"
            coords[i, 3] = 1.0
            continue
        end
        if use_polar
            si, ci = sincos(deg2rad(Float64(df.inc[row_idx])))
            sa, ca = sincos(deg2rad(Float64(df.azi[row_idx])))
            coords[i, 1] = si * ca
            coords[i, 2] = si * sa
            coords[i, 3] = ci
        else
            x, y, z = Float64(df.x3[row_idx]), Float64(df.y3[row_idx]), Float64(df.z3[row_idx])
            norm_v = hypot(x, y, z)
            if norm_v > 1e-12
                inv_norm = 1.0 / norm_v
                coords[i, 1] = x * inv_norm
                coords[i, 2] = y * inv_norm
                coords[i, 3] = z * inv_norm
            else
                warn_zero && @minimal_warning "Channel $ch has zero-norm 3D coordinates (0, 0, 0); defaulting to (0, 0, 1)"
                coords[i, 3] = 1.0
            end
        end
    end
    return coords
end

_calc_h(cosang::AbstractArray, m::Int = 4, n_legendre_terms::Int = 50) = _calc_g(cosang, m - 1, n_legendre_terms)

const _CSD_BLOCK_SIZE = 2048

function _ensure_float64_cols(df::DataFrame, channels::AbstractVector{Symbol})
    cols = Vector{Vector{Float64}}(undef, length(channels))
    for (c, ch) in enumerate(channels)
        col = df[!, ch]
        if col isa Vector{Float64}
            cols[c] = col
        else
            fcol = Float64.(col)
            df[!, ch] = fcol
            cols[c] = fcol
        end
    end
    return cols
end

function _apply_csd_df!(df::DataFrame, FT::Matrix{Float64}, channels::AbstractVector{Symbol}, X::Matrix{Float64}, Y::Matrix{Float64})
    cols = _ensure_float64_cols(df, channels)
    n_samp = nrow(df)
    n_ch = length(cols)
    blk = size(X, 1)
    for r0 = 1:blk:n_samp
        len = min(blk, n_samp - r0 + 1)
        rows = r0:(r0+len-1)
        Xv = len == blk ? X : view(X, 1:len, :)
        Yv = len == blk ? Y : view(Y, 1:len, :)
        for c = 1:n_ch
            copyto!(view(Xv, :, c), view(cols[c], rows))
        end
        mul!(Yv, Xv, FT)
        for c = 1:n_ch
            copyto!(view(cols[c], rows), view(Yv, :, c))
        end
    end
    return nothing
end

function _apply_csd!(data::SingleDataFrameEeg, filter_matrix::Matrix{Float64}, channels::AbstractVector{Symbol})
    blk = max(1, min(nrow(data.data), _CSD_BLOCK_SIZE))
    n_ch = length(channels)
    FT = collect(transpose(filter_matrix))
    _apply_csd_df!(data.data, FT, channels, Matrix{Float64}(undef, blk, n_ch), Matrix{Float64}(undef, blk, n_ch))
    return nothing
end

function _apply_csd!(data::EpochData, filter_matrix::Matrix{Float64}, channels::AbstractVector{Symbol})
    isempty(data.data) && return nothing
    FT = collect(transpose(filter_matrix))
    n_ch = length(channels)
    blk = max(1, min(maximum(nrow, data.data), _CSD_BLOCK_SIZE))
    n_chunks = min(Threads.nthreads(), length(data.data))
    @sync for chunk in Iterators.partition(data.data, cld(length(data.data), n_chunks))
        Threads.@spawn begin
            X = Matrix{Float64}(undef, blk, n_ch)
            Y = Matrix{Float64}(undef, blk, n_ch)
            for epoch_df in chunk
                _apply_csd_df!(epoch_df, FT, channels, X, Y)
            end
        end
    end
    return nothing
end

"""
    compute_csd(dat, args...; kwargs...)

Non-mutating version of [`compute_csd!`](@ref).
"""
compute_csd(dat::_CsdData, args...; kwargs...) = compute_csd!(copy(dat), args...; kwargs...)
compute_csd(dat::AbstractVector{<:_CsdData}, args...; kwargs...) = [compute_csd(d, args...; kwargs...) for d in dat]
