##########################################
# 2D topographic plot
##########################################
function plot_topoplot!(
    fig, 
    ax,
    dat,
    layout;
    xlim = nothing,
    ylim = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    method = :multiquadratic,  # :multiquadratic or :spherical_spline
    m = 4,  # for spherical spline (MNE-Python default)
    lambda = 1e-7,  # for spherical spline (more conservative regularization)
)


    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end
    
    # Ensure we have 3D coordinates for spherical spline
    if method == :spherical_spline && (:x3 ∉ propertynames(layout) || :y3 ∉ propertynames(layout) || :z3 ∉ propertynames(layout))
        polar_to_cartesian_xyz!(layout)
    end

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 400)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    if isnothing(xlim)
        xlim = [dat.time[1], dat.time[end]]
        xlim_idx = 1:nrow(dat)
    end

    # convert xlim to index
    xlim_idx = find_idx_range(dat.time, xlim[1], xlim[2])

    # interpolate data using chosen method
    if method == :spherical_spline
        # Ensure we have both 2D and 3D coordinates
        if !all(col -> col in propertynames(layout), [:x3, :y3, :z3])
            throw(ArgumentError("Layout must contain x3, y3, z3 columns"))
        end
        if !all(col -> col in propertynames(layout), [:x2, :y2])
            polar_to_cartesian_xy!(layout)
        end
        
        data = data_interpolation_topo_spherical_spline(
            mean.(eachcol(dat[xlim_idx, layout.label])),
            layout,
            gridscale,
            lambda=lambda
        )
    else  # default to multiquadratic
        # Ensure we have 2D coordinates
        if !all(col -> col in propertynames(layout), [:x2, :y2])
            polar_to_cartesian_xy!(layout)
        end
        
        data = data_interpolation_topo(
            mean.(eachcol(dat[xlim_idx, layout.label])),
            permutedims(Matrix(layout[!, [:x2, :y2]])),
            gridscale,
        )
    end

    if isnothing(ylim)
        ylim = minimum(data[.!isnan.(data)]), maximum(data[.!isnan.(data)])
    end

    radius = 88 # mm
    
    # Clear the axis to prevent double plotting
    empty!(ax)
    
    # Use different ranges based on interpolation method
    contour_range = method == :spherical_spline ? radius * 4 : radius * 2
    
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
    plot_layout_2d!(fig, ax, layout, head_kwargs = head_kwargs, point_kwargs = point_kwargs, label_kwargs = label_kwargs)

    return fig, ax

end

function plot_topoplot(
    dat,
    layout;
    kwargs...
)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat, layout; kwargs...)
    return fig, ax
end

"""
    plot_topoplot(dat::ContinuousData; kwargs...)

Create a topographic plot from continuous EEG data.

# Arguments
- `dat`: ContinuousData object
- `kwargs...`: Additional keyword arguments passed to plot_topoplot!

# Returns
- Figure and Axis objects
"""
function plot_topoplot(dat::ContinuousData; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
    return fig, ax
end

"""
    plot_topoplot!(fig, ax, dat::ContinuousData; kwargs...)

Add a topographic plot to existing figure/axis from continuous EEG data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: ContinuousData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topoplot!(fig, ax, dat::ContinuousData; kwargs...)
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
end

"""
    plot_topoplot(dat::EpochData, epoch::Int; kwargs...)

Create a topographic plot from epoched EEG data.

# Arguments
- `dat`: EpochData object
- `kwargs...`: Additional keyword arguments passed to plot_topoplot!

# Returns
- Figure and Axis objects
"""
function plot_topoplot(dat::EpochData, epoch::Int; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat.data[epoch], dat.layout; kwargs...)
    return fig, ax
end

"""
    plot_topoplot!(fig, ax, dat::EpochData, epoch::Int; kwargs...)

Add a topographic plot to existing figure/axis from epoched EEG data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: EpochData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topoplot!(fig, ax, dat::EpochData, epoch::Int; kwargs...)
    plot_topoplot!(fig, ax, dat.data[epoch], dat.layout; kwargs...)
end

"""
    plot_topoplot(dat::ErpData; kwargs...)

Create a topographic plot from ERP data.

# Arguments
- `dat`: ErpData object
- `kwargs...`: Additional keyword arguments passed to plot_topoplot!

# Returns
- Figure and Axis objects
"""
function plot_topoplot(dat::ErpData; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
    return fig, ax
end

"""
    plot_topoplot!(fig, ax, dat::ErpData; kwargs...)

Add a topographic plot to existing figure/axis from ERP data.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `dat`: ErpData object
- `kwargs...`: Additional keyword arguments
"""
function plot_topoplot!(fig, ax, dat::ErpData; kwargs...)
    plot_topoplot!(fig, ax, dat.data, dat.layout; kwargs...)
end


"""
    circle_mask!(dat::Matrix{<:AbstractFloat}, grid_scale::Int)

Apply a circular mask to a topographic data matrix, setting values outside
the circle to NaN.

# Arguments
- `dat::Matrix{Float64}`: Data matrix to mask (modified in-place)
- `grid_scale::Int`: Size of the grid (must be square matrix)

# Throws
- `ArgumentError`: If grid_scale is not positive
- `DimensionMismatch`: If matrix is not square
"""
function circle_mask!(dat::Matrix{<:AbstractFloat}, grid_scale::Int)
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
    data_interpolation_topo(dat::Vector{Float64}, points::Matrix{Float64}, grid_scale::Int)

Interpolate EEG data using scattered interpolation 

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `points::Matrix{<:AbstractFloat}`: 2×N matrix of electrode coordinates
- `grid_scale::Int`: Size of the output grid
"""
function data_interpolation_topo(dat::Vector{<:AbstractFloat}, points::Matrix{<:AbstractFloat}, grid_scale::Int)
    # Check input data
    if any(isnan, dat) || any(isinf, dat)
        throw(ArgumentError("Input data contains NaN or Inf values"))
    end
    if any(isnan, points) || any(isinf, points)
        throw(ArgumentError("Input points contain NaN or Inf values"))
    end

    # Create grid using the same coordinate range as the plotting function
    radius = 88.0  # mm - same as used in plotting function
    x_range = range(-radius * 2, radius * 2, length = grid_scale)
    y_range = range(-radius * 2, radius * 2, length = grid_scale)

    # Create regular grid
    grid_points = zeros(2, grid_scale^2)
    @inbounds for (idx, (i, j)) in enumerate(Iterators.product(x_range, y_range))
        grid_points[1, idx] = i
        grid_points[2, idx] = j
    end

    # Perform interpolation
    try
        itp = ScatteredInterpolation.interpolate(Multiquadratic(), points, dat)
        result = reshape(ScatteredInterpolation.evaluate(itp, grid_points), grid_scale, grid_scale)
        circle_mask!(result, grid_scale)
        return result
    catch e
        throw(ErrorException("Interpolation failed: $(e)"))
    end

end



"""
    data_interpolation_topo_spherical_spline(dat::Vector{Float64}, layout::DataFrame, grid_scale::Int; 
                                           m::Int=4, lambda::Float64=1e-5)

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
function data_interpolation_topo_spherical_spline(
    dat::Vector{<:AbstractFloat}, 
    layout::DataFrame, 
    grid_scale::Int;
    lambda::Float64=1e-5
)
    # Ensure we have 3D coordinates
    if !all(col -> col in propertynames(layout), [:x3, :y3, :z3])
        throw(ArgumentError("Layout must contain x3, y3, z3 columns"))
    end
    
    n_channels = length(dat)
    radius = 88.0  # mm - same as used in plotting function
    
    # Extract 3D coordinates as they are (no normalization needed)
    coords = zeros(Float64, n_channels, 3)
    for i in 1:n_channels
        coords[i, :] = [layout.x3[i], layout.y3[i], layout.z3[i]]
    end
    
    # Find the actual radius of the electrode positions
    electrode_radius = mean([sqrt(sum(coords[i,:].^2)) for i in 1:n_channels])
    
    # Create a 2D grid for plotting
    x_range = range(-radius * 2, radius * 2, length=grid_scale)
    y_range = range(-radius * 2, radius * 2, length=grid_scale)
    
    # For spherical spline calculation, we need unit sphere coordinates
    # So we normalize for the G matrix calculation only
    coords_unit = copy(coords)
    for i in 1:n_channels
        norm_factor = sqrt(sum(coords_unit[i,:].^2))
        if norm_factor > 0
            coords_unit[i,:] ./= norm_factor
        end
    end
    
    # Compute cosine angles between all electrode pairs (exactly like MNE)
    cosang = coords_unit * coords_unit'  # This is the dot product matrix
    
    # Compute G matrix using MNE's exact g-function
    G = calc_g_matrix(cosang)
    
    # Add regularization to diagonal (exactly like MNE)
    for i in 1:n_channels
        G[i,i] += lambda
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
    r_2d = sqrt.(grid_x.^2 .+ grid_y.^2)
    
    # Initialize result array
    interpolated_values = fill(NaN, length(grid_x))
    
    # Find valid grid points (within head and plotting area)
    valid_mask = (r_2d .<= radius) .& (r_2d .<= radius * 2.0)
    valid_indices = findall(valid_mask)
    
    if !isempty(valid_indices)
        # Extract valid grid points
        valid_x = grid_x[valid_indices]
        valid_y = grid_y[valid_indices]
        valid_r = r_2d[valid_indices]
        
        # Pre-compute stereographic projection for all valid points at once
        r_norm = valid_r ./ radius
        
        # Vectorized stereographic projection
        # Fix coordinate conversion to prevent 90-degree offset
        z3 = electrode_radius .* (1.0 .- r_norm.^2) ./ (1.0 .+ r_norm.^2)
        x3 = valid_x .* (1.0 .+ z3./electrode_radius)
        y3 = valid_y .* (1.0 .+ z3./electrode_radius)
        
        # Stack into 3D coordinates - ensure proper coordinate order
        grid_points_3d = hcat(x3, y3, z3)
        
        # Normalize to unit sphere for all points at once
        norms = sqrt.(sum(grid_points_3d.^2, dims=2))
        grid_points_unit = grid_points_3d ./ norms
        
        # Compute cosine angles to all electrodes for all grid points at once
        # Ensure proper matrix multiplication order
        cosang_grid = grid_points_unit * coords_unit'
        
        # Clamp cosine angles
        cosang_grid = clamp.(cosang_grid, -1.0, 1.0)
        
        # Compute g-function values for all grid points and electrodes at once
        g_values = calc_g_function.(cosang_grid)
        
        # Add the constant term (1.0) for each grid point
        g_values_extended = hcat(g_values, ones(size(g_values, 1)))
        
        # Compute interpolation for all valid points at once using matrix multiplication
        interpolated_valid = g_values_extended * weights
        
        # Store results back
        interpolated_values[valid_indices] = interpolated_valid
    end
    
    # Reshape to grid and apply circular mask
    result = reshape(interpolated_values, grid_scale, grid_scale)
    circle_mask!(result, grid_scale)
    
    return result
end

# Legendre polynomial calculation using iterative method
function legendre_polynomial(n::Int, x::Float64)
    n == 0 && return 1.0
    n == 1 && return x

    p_prev = 1.0  # P_0
    p_curr = x    # P_1
        
    for i in 2:n
        p_next = ((2i - 1) * x * p_curr - (i - 1) * p_prev) / i
        p_prev = p_curr
        p_curr = p_next
    end
        
    return p_curr

end

# Legendre polynomial evaluation for scalar input
function legendre_val(x::Float64, factors::Vector{Float64})
    """Evaluate Legendre polynomial series for scalar input."""
    result = 0.0
    for (i, factor) in enumerate(factors)
        if i == 1  # Skip the first factor (0.0)
            continue
        end
        n = i - 1  # Legendre polynomial order
        result += factor * legendre_polynomial(n, x)
    end
    return result
end

# MNE-Python's exact g-function calculation for EEG topography (m=4)
function calc_g_function(cosang::Float64, n_legendre_terms::Int=15)
    """Calculate spherical spline g function between points on a sphere.
    
    This is the exact implementation from MNE-Python, optimized for EEG topography (m=4).
    """
    cosang ≈ 1.0 && return 0.0
    
    # Use m=4 (standard for EEG) and fewer terms for speed
    factors = [
        (2 * n + 1) / (n^4 * (n + 1)^4 * 4 * π)
        for n in 1:n_legendre_terms
    ]
    
    # Use Legendre polynomial evaluation
    return legendre_val(cosang, [0.0; factors])
end

# MNE-Python's exact G matrix calculation for EEG topography (m=4)
function calc_g_matrix(cosang::Matrix{Float64}, n_legendre_terms::Int=15)
    """Calculate spherical spline G matrix between points on a sphere.
    
    This is the exact implementation from MNE-Python, optimized for EEG topography (m=4).
    """
    factors = [
        (2 * n + 1) / (n^4 * (n + 1)^4 * 4 * π)
        for n in 1:n_legendre_terms
    ]
    
    # Use Legendre polynomial evaluation for the entire matrix
    return legendre_val.(cosang, Ref([0.0; factors]))
end