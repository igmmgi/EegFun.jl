"""
    compute_csd!(data; m=4, lambda=1e-5, n_legendre_terms=50)

Compute the Current Source Density (CSD) / Surface Laplacian of the data in-place.
This acts as a spatial high-pass filter that removes volume conduction and 
sharpens the topography.

Uses spherical spline interpolation according to Perrin et al. (1989).

# Arguments
- `data`: EEG data object (ContinuousData, EpochData, or ErpData)
- `m::Int`: Stiffness parameter for the spherical splines (default: 4)
- `lambda::Real`: Regularization parameter (default: 1e-5)
- `n_legendre_terms::Int`: Number of Legendre polynomial terms (default: 50)
"""
function compute_csd!(data::EegData; m::Int = 4, lambda::Real = 1e-5, n_legendre_terms::Int = 50)
    # 1. Ensure 3D coordinates are available
    _ensure_coordinates_3d!(data.layout)

    channels = Symbol.(data.layout.data.label)
    n_channels = length(channels)

    # 2. Extract coordinates and normalize to unit sphere
    coords = zeros(Float64, n_channels, 3)
    x3 = data.layout.data.x3
    y3 = data.layout.data.y3
    z3 = data.layout.data.z3

    for i = 1:n_channels
        x = Float64(x3[i])
        y = Float64(y3[i])
        z = Float64(z3[i])
        inv_norm = 1.0 / sqrt(x^2 + y^2 + z^2)
        coords[i, 1] = x * inv_norm
        coords[i, 2] = y * inv_norm
        coords[i, 3] = z * inv_norm
    end

    # 3. Calculate cosine angles between all electrodes
    cosang = coords * coords'

    # 4. Calculate G and H matrices
    G = _calc_g(cosang, m, n_legendre_terms)
    H = _calc_h(cosang, m, n_legendre_terms)

    # 5. Add regularization to diagonal of G (MNE/FieldTrip standard scaling)
    if lambda > 0
        trace_G = sum(G[i, i] for i = 1:n_channels)
        lambda_scaled = lambda * trace_G / n_channels
        for i = 1:n_channels
            G[i, i] += lambda_scaled
        end
    end

    # 6. Construct system matrix C
    # C = [G     ones(n)]
    #     [ones' 0      ]
    C = zeros(n_channels + 1, n_channels + 1)
    C[1:n_channels, 1:n_channels] = G
    C[1:n_channels, n_channels+1] .= 1.0
    C[n_channels+1, 1:n_channels] .= 1.0
    C[n_channels+1, n_channels+1] = 0.0

    # Compute inverse
    C_inv = pinv(C)

    # 7. Compute the spatial filter matrix
    # filter_matrix = H * C_inv[1:n_channels, 1:n_channels]
    filter_matrix = H * C_inv[1:n_channels, 1:n_channels]

    # 8. Apply filter to data
    _apply_csd!(data, filter_matrix, channels)

    return data
end

function compute_csd!(dat::AbstractVector{<:EegData}, args...; kwargs...)
    for d in dat
        compute_csd!(d, args...; kwargs...)
    end
    return dat
end

# Internal function for calculating H matrix (second spatial derivative of G)
function _calc_h(cosang::AbstractArray, m::Int = 4, n_legendre_terms::Int = 50)
    # H(cos θ) = - (1 / 4π) * Σ[ (2n+1) / (n^(m-1) * (n+1)^(m-1)) * P_n(cos θ) ]
    factors = [-(2n + 1) / (float(n)^(m - 1) * float(n + 1)^(m - 1) * 4π) for n = 1:n_legendre_terms]
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
function compute_csd(dat::EegData, args...; kwargs...)
    return compute_csd!(copy(dat), args...; kwargs...)
end

function compute_csd(dat::AbstractVector{<:EegData}, args...; kwargs...)
    return [compute_csd(d, args...; kwargs...) for d in dat]
end
