##########################################
# 2D topographic plot
##########################################
function _plot_topography!(
    fig::Figure,
    ax::Axis,
    dat::DataFrame,
    layout::Layout;
    ylim = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    method = :multiquadratic,  # :multiquadratic or :spherical_spline
)

    # ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    # deal with kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 200)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    # actual data interpolation
    channel_data = mean.(eachcol(dat[!, layout.data.label]))
    if method == :spherical_spline
        data = _data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
    elseif method == :multiquadratic
        data = _data_interpolation_topo_multiquadratic(channel_data, layout, gridscale)
    end

    if isnothing(ylim)
        ylim = minimum(data[.!isnan.(data)]), maximum(data[.!isnan.(data)])
    end

    # Clear the axis 
    empty!(ax)

    # Use different ranges based on interpolation method
    contour_range = method == :spherical_spline ? DEFAULT_HEAD_RADIUS * 4 : DEFAULT_HEAD_RADIUS * 2

    co = contourf!(
        range(-contour_range, contour_range, length = gridscale),
        range(-contour_range, contour_range, length = gridscale),
        data,
        levels = range(ylim[1], ylim[2], div(gridscale, 2));
        topo_kwargs...,
    )

    if plot_colorbar
        Colorbar(fig[1, 2], co; colorbar_kwargs...)
    end

    # head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    return fig, ax

end


"""
    plot_topography!(fig, ax, dat::EpochData, epoch::Int; kwargs...)

Add a topographic plot to existing figure/axis from epoched EEG data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: EpochData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topography!(fig, ax, dat::EpochData, epoch::Int; kwargs...)
    plot_topography!(fig, ax, dat.data[epoch], dat.layout; kwargs...)
end

function plot_topography(dat::EpochData, epoch::Int; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topography!(fig, ax, dat.data[epoch], dat.layout; kwargs...)
    return fig, ax
end



"""
    plot_topography!(fig, ax, dat::SingleDataFrameEeg; kwargs...)

Add a topographic plot to existing figure/axis from single DataFrame EEG data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: SingleDataFrameEeg object (ContinuousData or ErpData)
- `kwargs...`: Additional keyword arguments
"""
function plot_topography!(
    fig,
    ax,
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    kwargs...,
)
    dat_subset = subset(dat, channel_selection = channel_selection, sample_selection = sample_selection)
    plot_topography!(fig, ax, dat_subset.data, dat_subset.layout; kwargs...)
end

function plot_topography(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot = true,
    kwargs...,
)
    fig = Figure()
    ax = Axis(fig[1, 1])
    dat_subset = subset(dat, channel_selection = channel_selection, sample_selection = sample_selection)
    _plot_topography!(fig, ax, dat_subset.data, dat_subset.layout; kwargs...)
    if display_plot
        display_figure(fig)
    end
    return fig, ax
end


"""
    _circle_mask!(dat::Matrix{<:AbstractFloat}, grid_scale::Int)

Apply a circular mask to a topographic data matrix, setting values outside
the circle to NaN.

# Arguments
- `dat::Matrix{Float64}`: Data matrix to mask (modified in-place)
- `grid_scale::Int`: Size of the grid (must be square matrix)

# Throws
- `ArgumentError`: If grid_scale is not positive
- `DimensionMismatch`: If matrix is not square
"""
function _circle_mask!(dat::Matrix{<:AbstractFloat}, grid_scale::Int)
    if grid_scale <= 0
        throw(ArgumentError("grid_scale must be positive"))
    end
    if size(dat) != (grid_scale, grid_scale)
        throw(DimensionMismatch("Data matrix must be $(grid_scale)×$(grid_scale)"))
    end

    center = grid_scale / 2
    @inbounds for col = 1:grid_scale
        for row = 1:grid_scale
            x_dist = center - col
            y_dist = center - row
            if sqrt(x_dist^2 + y_dist^2) > center
                dat[col, row] = NaN
            end
        end
    end
    return dat
end



"""
    _data_interpolation_topo(dat::Vector{Float64}, points::Matrix{Float64}, grid_scale::Int)

Interpolate EEG data using scattered interpolation 

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `points::Matrix{<:AbstractFloat}`: 2×N matrix of electrode coordinates
- `grid_scale::Int`: Size of the output grid
"""
function _data_interpolation_topo_multiquadratic(
    dat::Vector{<:AbstractFloat},
    layout::Layout,
    grid_scale::Int,
)
    # Check input data
    if any(isnan, dat) || any(isinf, dat)
        throw(ArgumentError("Input data contains NaN or Inf values"))
    end

    points = permutedims(Matrix(layout.data[!, [:x2, :y2]]))
    
    # Create grid more efficiently - avoid collect() and Iterators.product
    x_range = range(-DEFAULT_HEAD_RADIUS * 2, DEFAULT_HEAD_RADIUS * 2, length = grid_scale)
    y_range = range(-DEFAULT_HEAD_RADIUS * 2, DEFAULT_HEAD_RADIUS * 2, length = grid_scale)

    # Create regular grid more efficiently
    grid_points = zeros(2, grid_scale^2)
    idx = 1
    @inbounds for y in y_range
        for x in x_range
            grid_points[1, idx] = x
            grid_points[2, idx] = y
            idx += 1
        end
    end

    # Perform scattered interpolation (essential for proper EEG topography)
    try
        itp = ScatteredInterpolation.interpolate(Multiquadratic(), points, dat)
        result = reshape(ScatteredInterpolation.evaluate(itp, grid_points), grid_scale, grid_scale)
        _circle_mask!(result, grid_scale)
        return result
    catch e
        throw(ErrorException("Interpolation failed: $(e)"))
    end

end




"""
    _data_interpolation_topo_spherical_spline(dat::Vector{Float64}, layout::DataFrame, grid_scale::Int; m::Int=4, lambda::Float64=1e-5)

Interpolate EEG data using spherical spline interpolation for topographic plotting.
Implementation follows MNE-Python exactly.

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `layout::DataFrame`: Layout information with x3, y3, z3 coordinates
- `grid_scale::Int`: Size of the output grid
- `m::Int=4`: Order of Legendre polynomials (stiffness parameter)
- `lambda::Float64=1e-5`: Regularization parameter

# Returns
- `Matrix{Float64}`: Interpolated data on a regular grid
"""
function _data_interpolation_topo_spherical_spline(dat::Vector{<:AbstractFloat}, layout::Layout, grid_scale::Int;)

    # Extract 3D coordinates as they are (no normalization needed)
    n_channels = length(dat)
    coords = zeros(Float64, n_channels, 3)
    for i = 1:n_channels
        coords[i, :] = [layout.data.x3[i], layout.data.y3[i], layout.data.z3[i]]
    end

    # Find the actual radius of the electrode positions
    electrode_radius = mean([sqrt(sum(coords[i, :] .^ 2)) for i = 1:n_channels])

    # Create a 2D grid for plotting
    x_range = range(-DEFAULT_HEAD_RADIUS * 2, DEFAULT_HEAD_RADIUS * 2, length = grid_scale)
    y_range = range(-DEFAULT_HEAD_RADIUS * 2, DEFAULT_HEAD_RADIUS * 2, length = grid_scale)

    # For spherical spline calculation, we need unit sphere coordinates
    # So we normalize for the G matrix calculation only
    coords_unit = copy(coords)
    for i = 1:n_channels
        norm_factor = sqrt(sum(coords_unit[i, :] .^ 2))
        if norm_factor > 0
            coords_unit[i, :] ./= norm_factor
        end
    end

    # Compute cosine angles between all electrode pairs (exactly like MNE)
    cosang = coords_unit * coords_unit'  # This is the dot product matrix

    # Compute G matrix using MNE's exact g-function
    G = _calc_g_matrix(cosang)

    # Add regularization to diagonal (exactly like MNE)
    for i = 1:n_channels
        G[i, i] += 1e-5
    end

    # Add constraint rows/columns (exactly like MNE)
    G_extended = zeros(Float64, n_channels+1, n_channels+1)
    G_extended[1:n_channels, 1:n_channels] = G
    G_extended[1:n_channels, n_channels+1] .= 1.0
    G_extended[n_channels+1, 1:n_channels] .= 1.0
    G_extended[n_channels+1, n_channels+1] = 0.0

    # Solve the system to get weights (exactly like MNE)
    data_vector = vcat(dat, 0.0)
    weights = G_extended \ data_vector

    # Pre-compute all grid points for vectorized operations
    # Use the same order as the original nested loops: for x in x_range, y in y_range
    grid_x = vec(repeat(x_range', length(y_range), 1))
    grid_y = repeat(y_range, length(x_range))

    # Calculate distances from center for all grid points at once
    r_2d = sqrt.(grid_x .^ 2 .+ grid_y .^ 2)

    # Initialize result array
    interpolated_values = fill(NaN, length(grid_x))

    # Find valid grid points (within head and plotting area)
    valid_mask = (r_2d .<= DEFAULT_HEAD_RADIUS) .& (r_2d .<= DEFAULT_HEAD_RADIUS * 2.0)
    valid_indices = findall(valid_mask)

    if !isempty(valid_indices)
        # Extract valid grid points
        valid_x = grid_x[valid_indices]
        valid_y = grid_y[valid_indices]
        valid_r = r_2d[valid_indices]

        # Pre-compute stereographic projection for all valid points at once
        r_norm = valid_r ./ DEFAULT_HEAD_RADIUS

        # Vectorized stereographic projection
        # Fix coordinate conversion to prevent 90-degree offset
        z3 = electrode_radius .* (1.0 .- r_norm .^ 2) ./ (1.0 .+ r_norm .^ 2)
        # Swap x and y to fix the 90-degree rotation
        x3 = valid_y .* (1.0 .+ z3 ./ electrode_radius)
        y3 = valid_x .* (1.0 .+ z3 ./ electrode_radius)

        # Stack into 3D coordinates - ensure proper coordinate order
        grid_points_3d = hcat(x3, y3, z3)

        # Normalize to unit sphere for all points at once
        norms = sqrt.(sum(grid_points_3d .^ 2, dims = 2))
        grid_points_unit = grid_points_3d ./ norms

        # Compute cosine angles to all electrodes for all grid points at once
        # Ensure proper matrix multiplication order
        cosang_grid = grid_points_unit * coords_unit'

        # Clamp cosine angles
        cosang_grid = clamp.(cosang_grid, -1.0, 1.0)

        # Compute g-function values for all grid points and electrodes at once
        g_values = _calc_g_function.(cosang_grid)

        # Add the constant term (1.0) for each grid point
        g_values_extended = hcat(g_values, ones(size(g_values, 1)))

        # Compute interpolation for all valid points at once using matrix multiplication
        interpolated_valid = g_values_extended * weights

        # Store results back
        interpolated_values[valid_indices] = interpolated_valid
    end

    # Reshape to grid and apply circular mask
    result = reshape(interpolated_values, grid_scale, grid_scale)
    _circle_mask!(result, grid_scale)

    return result
end

# Legendre polynomial calculation using iterative method
function _legendre_polynomial(n::Int, x::Float64)
    n == 0 && return 1.0
    n == 1 && return x

    p_prev = 1.0  # P_0
    p_curr = x    # P_1

    for i = 2:n
        p_next = ((2i - 1) * x * p_curr - (i - 1) * p_prev) / i
        p_prev = p_curr
        p_curr = p_next
    end

    return p_curr

end

# Legendre polynomial evaluation for scalar input
function _legendre_val(x::Float64, factors::Vector{Float64})
    """Evaluate Legendre polynomial series for scalar input."""
    result = 0.0
    for (i, factor) in enumerate(factors)
        if i == 1  # Skip the first factor (0.0)
            continue
        end
        n = i - 1  # Legendre polynomial order
        result += factor * _legendre_polynomial(n, x)
    end
    return result
end

# MNE-Python's exact g-function calculation for EEG topography (m=4)
function _calc_g_function(cosang::Float64, n_legendre_terms::Int = 15)
    """Calculate spherical spline g function between points on a sphere.

    This is the exact implementation from MNE-Python, optimized for EEG topography (m=4).
    """
    cosang ≈ 1.0 && return 0.0

    # Use m=4 (standard for EEG) and fewer terms for speed
    factors = [(2 * n + 1) / (n^4 * (n + 1)^4 * 4 * π) for n = 1:n_legendre_terms]

    # Use Legendre polynomial evaluation
    return _legendre_val(cosang, [0.0; factors])
end

# MNE-Python's exact G matrix calculation for EEG topography (m=4)
function _calc_g_matrix(cosang::Matrix{Float64}, n_legendre_terms::Int = 15)
    """Calculate spherical spline G matrix between points on a sphere.

    This is the exact implementation from MNE-Python, optimized for EEG topography (m=4).
    """
    factors = [(2 * n + 1) / (n^4 * (n + 1)^4 * 4 * π) for n = 1:n_legendre_terms]

    # Use Legendre polynomial evaluation for the entire matrix
    return _legendre_val.(cosang, Ref([0.0; factors]))
end
