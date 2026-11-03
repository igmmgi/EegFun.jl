const _CsdData = Union{ContinuousData,ErpData,EpochData}

"""
    compute_csd!(data; m=4, lambda=1e-5, n_legendre_terms=50)

Compute the Current Source Density (CSD) / Surface Laplacian of the data in-place.
This acts as a spatial high-pass filter that removes volume conduction and 
sharpens the topography.

Uses spherical spline interpolation according to Perrin et al. (1989; 1990 Corrigenda).

# Arguments
- `data`: EEG data object (`ContinuousData`, `EpochData`, or `ErpData`)
- `m::Int`: Stiffness parameter for the spherical splines (default: 4)
- `lambda::Real`: Regularization parameter (default: 1e-5)
- `n_legendre_terms::Int`: Number of Legendre polynomial terms (default: 50)
"""
function compute_csd!(data::_CsdData; m::Int = 4, lambda::Real = 1e-5, n_legendre_terms::Int = 50)
    # 1. Ensure 3D coordinates are available and extract unit-sphere coordinates
    channels = channel_labels(data)
    isempty(channels) && @minimal_error("No EEG channels found in data for CSD computation")
    n_channels = length(channels)

    coords = _spherical_coords(data.layout, channels)

    # 2. Calculate cosine angles between all electrodes
    cosang = coords * coords'

    # 3. Calculate G and H matrices
    G = _calc_g(cosang, m, n_legendre_terms)
    H = _calc_h(cosang, m, n_legendre_terms)

    # 4. Add regularization to diagonal of G (MNE/FieldTrip standard scaling)
    if lambda > 0
        trace_G = sum(G[i, i] for i = 1:n_channels)
        lambda_scaled = lambda * trace_G / n_channels
        for i = 1:n_channels
            G[i, i] += lambda_scaled
        end
    end

    # 5. Construct system matrix C
    # C = [G     ones(n)]
    #     [ones' 0      ]
    C = zeros(n_channels + 1, n_channels + 1)
    C[1:n_channels, 1:n_channels] = G
    C[1:n_channels, n_channels+1] .= 1.0
    C[n_channels+1, 1:n_channels] .= 1.0
    C[n_channels+1, n_channels+1] = 0.0

    # Compute inverse
    C_inv = pinv(C)

    # 6. Compute the spatial filter matrix
    # filter_matrix = H * C_inv[1:n_channels, 1:n_channels]
    filter_matrix = H * C_inv[1:n_channels, 1:n_channels]

    # 7. Apply filter to data and record reference
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

# Extract unit-sphere 3D coordinates for channels, guarding zero-norm vectors
function _spherical_coords(layout::Layout, channels::AbstractVector{Symbol} = Symbol.(layout.data.label))
    _ensure_coordinates_3d!(layout)
    df = layout.data
    label_to_idx = Dict(Symbol(df.label[i]) => i for i = 1:nrow(df))
    n_channels = length(channels)
    coords = zeros(Float64, n_channels, 3)

    for (i, ch) in enumerate(channels)
        row_idx = get(label_to_idx, ch, nothing)
        if isnothing(row_idx)
            @minimal_warning "Channel $ch not found in layout; defaulting unit-sphere coordinate to (0, 0, 1)"
            coords[i, 3] = 1.0
            continue
        end
        x = Float64(df.x3[row_idx])
        y = Float64(df.y3[row_idx])
        z = Float64(df.z3[row_idx])
        norm_v = sqrt(x^2 + y^2 + z^2)
        if norm_v > 1e-12
            inv_norm = 1.0 / norm_v
            coords[i, 1] = x * inv_norm
            coords[i, 2] = y * inv_norm
            coords[i, 3] = z * inv_norm
        else
            @minimal_warning "Channel $ch has zero-norm 3D coordinates (0, 0, 0); defaulting to (0, 0, 1)"
            coords[i, 1] = 0.0
            coords[i, 2] = 0.0
            coords[i, 3] = 1.0
        end
    end
    return coords
end

# Internal function for calculating H matrix (negative surface Laplacian of G; Perrin et al., 1990 Corrigenda)
function _calc_h(cosang::AbstractArray, m::Int = 4, n_legendre_terms::Int = 50)
    # H(cos θ) = (1 / 4π) * Σ[ (2n+1) / (n^(m-1) * (n+1)^(m-1)) * P_n(cos θ) ]
    factors = [(2n + 1) / (float(n)^(m - 1) * float(n + 1)^(m - 1) * 4π) for n = 1:n_legendre_terms]
    C = Matrix{Float64}(cosang)
    return _g_from_cos!(C, factors)
end

# --- Type-specific application methods ---

function _apply_csd!(data::SingleDataFrameEeg, filter_matrix::Matrix{Float64}, channels::AbstractVector{Symbol})
    n_samples = nrow(data.data)
    n_chans = length(channels)
    X = Matrix{Float64}(undef, n_samples, n_chans)
    for (c, ch) in enumerate(channels)
        X[:, c] = data.data[!, ch]
    end
    FT = collect(transpose(filter_matrix))
    Y = Matrix{Float64}(undef, n_samples, n_chans)
    mul!(Y, X, FT)
    for (c, ch) in enumerate(channels)
        data.data[!, ch] = Y[:, c]
    end
end

function _apply_csd!(data::EpochData, filter_matrix::Matrix{Float64}, channels::AbstractVector{Symbol})
    FT = collect(transpose(filter_matrix))
    n_chans = length(channels)
    Threads.@threads for epoch_df in data.data
        n_samples = nrow(epoch_df)
        X = Matrix{Float64}(undef, n_samples, n_chans)
        for (c, ch) in enumerate(channels)
            X[:, c] = epoch_df[!, ch]
        end
        Y = Matrix{Float64}(undef, n_samples, n_chans)
        mul!(Y, X, FT)
        for (c, ch) in enumerate(channels)
            epoch_df[!, ch] = Y[:, c]
        end
    end
end

"""
    compute_csd(dat, args...; kwargs...)

Non-mutating version of `compute_csd!`.
"""
function compute_csd(dat::_CsdData, args...; kwargs...)
    return compute_csd!(copy(dat), args...; kwargs...)
end

function compute_csd(dat::AbstractVector{<:_CsdData}, args...; kwargs...)
    return [compute_csd(d, args...; kwargs...) for d in dat]
end
