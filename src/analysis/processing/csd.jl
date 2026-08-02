"""
    compute_csd!(data; m=4, lambda=1e-5, n_legendre_terms=50)

Compute the Current Source Density (CSD) / Surface Laplacian of the data in-place.
This acts as a spatial high-pass filter that removes volume conduction and 
sharpens the topography.

Uses spherical spline interpolation according to Perrin et al. (1989).

# Arguments
- `data`: EEG data object (ContinuousData, EpochData, or ErpData)
- `m::Int`: Stiffness parameter for the spherical splines (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)
- `n_legendre_terms::Int`: Number of Legendre polynomial terms (default: 50)
"""
function compute_csd!(data; m::Int = 4, lambda::Float64 = 1e-5, n_legendre_terms::Int = 50)
    # 1. Ensure 3D coordinates are available
    _ensure_coordinates_3d!(data.layout)

    # Performance Note:
    # The G and H matrices depend purely on the physical coordinates of the electrodes.
    # Currently, they are re-calculated on every function call. For optimization, we could 
    # cache `filter_matrix` directly inside the `Layout` object in the future.

    channels = Symbol.(data.layout.data.label)
    n_channels = length(channels)

    # 2. Extract coordinates
    coords = zeros(Float64, n_channels, 3)
    for (i, ch) in enumerate(channels)
        idx = findfirst(x -> x == ch, data.layout.data.label)
        coords[i, :] = [data.layout.data.x3[idx], data.layout.data.y3[idx], data.layout.data.z3[idx]]
    end

    # Normalize to unit sphere
    for i = 1:n_channels
        norm_val = sqrt(sum(coords[i, :] .^ 2))
        coords[i, :] ./= norm_val
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

# Internal function for calculating H matrix (second spatial derivative of G)
function _calc_h(cosang::AbstractArray, m::Int = 4, n_legendre_terms::Int = 50)
    # H(cos θ) = - (1 / 4π) * Σ[ (2n+1) / (n^(m-1) * (n+1)^(m-1)) * P_n(cos θ) ]
    factors = zeros(n_legendre_terms + 1)
    factors[1] = 0.0

    for n = 1:n_legendre_terms
        factors[n+1] = -(2 * n + 1) / (n^(m-1) * (n + 1)^(m-1) * 4 * π)
    end

    return _eval_legendre(cosang, factors)
end

# --- Type-specific application methods ---

function _apply_csd!(data::ContinuousData, filter_matrix::Matrix{Float64}, channels::Vector{Symbol})
    # Extract data matrix (channels x timepoints)
    data_matrix = Matrix(data.data[:, channels])'

    # Apply spatial filter
    filtered_matrix = filter_matrix * data_matrix

    # Update data
    data.data[:, channels] = filtered_matrix'
end

function _apply_csd!(data::EpochData, filter_matrix::Matrix{Float64}, channels::Vector{Symbol})
    for epoch_df in data.data
        data_matrix = Matrix(epoch_df[:, channels])'
        filtered_matrix = filter_matrix * data_matrix
        epoch_df[:, channels] = filtered_matrix'
    end
end

function _apply_csd!(data::ErpData, filter_matrix::Matrix{Float64}, channels::Vector{Symbol})
    data_matrix = Matrix(data.data[:, channels])'
    filtered_matrix = filter_matrix * data_matrix
    data.data[:, channels] = filtered_matrix'
end

@add_nonmutating compute_csd!
