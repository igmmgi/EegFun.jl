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
function data_interpolation_topo(dat::Vector{<:AbstractFloat}, points::Matrix{<:AbstractFloat}, grid_scale::Int;)
    # Check input data
    if any(isnan, dat) || any(isinf, dat)
        throw(ArgumentError("Input data contains NaN or Inf values"))
    end
    if any(isnan, points) || any(isinf, points)
        throw(ArgumentError("Input points contain NaN or Inf values"))
    end

    # Create grid based on the range of input points
    x_range = range(minimum(points[1, :]), maximum(points[1, :]), length = grid_scale)
    y_range = range(minimum(points[2, :]), maximum(points[2, :]), length = grid_scale)

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