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
    m = 4,  # for spherical spline
    lambda = 1e-5,  # for spherical spline
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

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
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
            m=m,
            lambda=lambda
        )
        println(size(data))
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
    
    co = contourf!(
       range(-radius * 2, radius * 2, length = gridscale),
       range(-radius * 2, radius * 2, length = gridscale),
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



# epochs = []
# for (idx, epoch) in enumerate([1, 4, 5, 3])
#      push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
# end
# 
# # Continuous Data
# plot_topoplot(dat)
# plot_topoplot(dat.data, dat.layout)
# 
# # Epoch Data
# epoch = extract_epochs(dat, 1, 1, -2, 4)
# plot_topoplot(epoch, 1) # 1st epoch
# plot_topoplot(epoch, 2) # 2nd epoch
# 
# # ERP Data
# erp = average_epochs(epochs)
# plot_topoplot(erp)
# 
# # Try some keyword arguments
# plot_topoplot(erp, xlim = (-0.1, 0.2), ylim = (-10, 10), 
#     head_kwargs = Dict(:linewidth => 5),
#     point_kwargs = Dict(:markersize => 20),
#     label_kwargs = Dict(:fontsize => 12),
#     topo_kwargs = Dict(:colormap => :viridis),
#     colorbar_kwargs = Dict(:width => 20),
# )


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

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `layout::DataFrame`: Layout information with x3, y3, z3 coordinates
- `grid_scale::Int`: Size of the output grid
- `m::Int=4`: Order of Legendre polynomials
- `lambda::Float64=1e-5`: Regularization parameter

# Returns
- `Matrix{Float64}`: Interpolated data on a regular grid
"""
function data_interpolation_topo_spherical_spline(
    dat::Vector{<:AbstractFloat}, 
    layout::DataFrame, 
    grid_scale::Int;
    m::Int=4, 
    lambda::Float64=1e-5
)
    # Ensure we have both 2D and 3D coordinates
    if !all(col -> col in propertynames(layout), [:x3, :y3, :z3])
        throw(ArgumentError("Layout must contain x3, y3, z3 columns"))
    end
    if !all(col -> col in propertynames(layout), [:x2, :y2])
        polar_to_cartesian_xy!(layout)
    end
    
    # Extract and normalize coordinates
    n_channels = length(dat)
    coords = zeros(Float64, n_channels, 3)
    for i in 1:n_channels
        coords[i, :] = [layout.x3[i], layout.y3[i], layout.z3[i]]
        # Normalize to unit sphere
        norm_factor = sqrt(sum(coords[i,:].^2))
        if norm_factor > 0
            coords[i,:] ./= norm_factor
        end
    end
    
    # Create a 2D grid that matches the plotting function exactly
    # Use the same coordinate range as the plotting function expects
    radius = 88.0  # mm - same as used in plotting function
    x_range = range(-radius * 2, radius * 2, length=grid_scale)
    y_range = range(-radius * 2, radius * 2, length=grid_scale)
    
    # Create 2D grid points
    grid_points_2d = zeros(2, grid_scale^2)
    grid_idx = 1
    for x in x_range, y in y_range
        grid_points_2d[1, grid_idx] = x
        grid_points_2d[2, grid_idx] = y
        grid_idx += 1
    end
    
    # Pre-compute g-function matrix for all electrode pairs
    G = zeros(Float64, n_channels+1, n_channels+1)
    for i in 1:n_channels, j in 1:n_channels
        if i != j
            cos_angle = dot(coords[i,:], coords[j,:])
            cos_angle = clamp(cos_angle, -1.0, 1.0)
            G[i,j] = fast_g_function(cos_angle, m)
        end
    end
    
    # Add regularization to diagonal
    for i in 1:n_channels
        G[i,i] = fast_g_function(1.0, m) + lambda
    end
    
    # Add constraint rows/columns
    G[1:n_channels, n_channels+1] .= 1.0
    G[n_channels+1, 1:n_channels] .= 1.0
    G[n_channels+1, n_channels+1] = 0.0
    
    # Pre-compute factorization
    G_fact = lu(G)
    
    # Create RHS for the data
    data_vector = vcat(dat, 0.0)
    
    # Solve the system to get weights
    weights = G_fact \ data_vector
    
    # Interpolate to grid points
    interpolated_values = zeros(Float64, grid_scale^2)
    for i in 1:grid_scale^2
        x2, y2 = grid_points_2d[1, i], grid_points_2d[2, i]
        
        # Calculate distance from center
        r_2d = sqrt(x2^2 + y2^2)
        
        if r_2d <= radius * 2.0  # Within full plotting area
            # Use spherical spline with proper 3D projection for accuracy
            if r_2d <= radius
                # Within head - use proper 3D spherical projection
                # Swap x and y to match the corrected 3D coordinate system
                z3 = sqrt(radius^2 - r_2d^2)
                grid_point_3d = [y2, x2, z3]  # Swap x2 and y2 to match 3D coords
            else
                # Beyond head - project to sphere surface by scaling
                scale_factor = radius / r_2d
                x_scaled = x2 * scale_factor
                y_scaled = y2 * scale_factor
                # Ensure we don't get negative values under sqrt due to floating point precision
                z_squared = radius^2 - x_scaled^2 - y_scaled^2
                z3 = z_squared > 0 ? sqrt(z_squared) : 0.0
                # Swap x and y to match the corrected 3D coordinate system
                grid_point_3d = [y_scaled, x_scaled, z3]  # Swap x_scaled and y_scaled
            end
            
            # Normalize to unit sphere for interpolation
            norm_factor = sqrt(sum(grid_point_3d.^2))
            if norm_factor > 0
                grid_point_3d ./= norm_factor
            end
            
            # Compute g-function values using proper 3D angles
            g_values = zeros(Float64, n_channels+1)
            for j in 1:n_channels
                cos_angle = dot(grid_point_3d, coords[j,:])
                cos_angle = clamp(cos_angle, -1.0, 1.0)
                g_values[j] = fast_g_function(cos_angle, m)
            end
            g_values[n_channels+1] = 1.0
            
            # Use the pre-computed spherical spline weights
            interpolated_values[i] = dot(g_values, weights)
        else
            # Outside plotting area - set to NaN
            interpolated_values[i] = NaN
        end
    end
    
    # Reshape to grid and apply circular mask
    result = reshape(interpolated_values, grid_scale, grid_scale)
    circle_mask!(result, grid_scale)
    
    return result
end

# Spherical spline interpolation based on Perrin et al. (1989)
# Fast approximations derived from MNE-Python/EEGLAB implementations

# Much faster g-function with pre-computed coefficients
function fast_g_function(cos_angle::Float64, m::Int)
    if cos_angle ≈ 1.0
        return 0.0
    end
    
    # Use a table-based approach with common m values
    if m == 4  # Most common case
        # These coefficients work well for m=4 and avoid computing Legendre polynomials
        if cos_angle > 0.9999
            return 0.0
        elseif cos_angle > 0.99
            return 0.05 * (1.0 - cos_angle)
        else
            # Approximation that works well for m=4
            return 0.25 * (1.0 - cos_angle) * (log(1.0 - cos_angle) + 0.5)
        end
    else
        # For other m values, use a faster approximation
        # This is based on a simplified version of the sum used in typical implementations
        if cos_angle ≈ 1.0
            return 0.0
        end
        
        # Use a fixed number of terms - first 7-10 terms usually sufficient
        sum_value = 0.0
        for n in 1:10
            # Use a simpler polynomial approximation instead of full Legendre
            # For common angles, this is accurate enough
            P_n = legendre_approx(n, cos_angle)
            sum_value += ((2*n + 1) / (n^m * (n+1)^m)) * P_n
        end
        
        return sum_value / (4 * π)
    end
end

# Approximation of Legendre polynomials for common cases
function legendre_approx(n::Int, x::Float64)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    elseif n == 2
        return (3x^2 - 1)/2
    elseif n == 3
        return (5x^3 - 3x)/2
    elseif n == 4
        return (35x^4 - 30x^2 + 3)/8
    else
        # For higher n, use the iterative method
        p_prev = (35x^4 - 30x^2 + 3)/8  # P_4
        p_curr = (63x^5 - 70x^3 + 15x)/8  # P_5
        
        for i in 6:n
            p_next = ((2i - 1) * x * p_curr - (i - 1) * p_prev) / i
            p_prev = p_curr
            p_curr = p_next
        end
        
        return p_curr
    end
end