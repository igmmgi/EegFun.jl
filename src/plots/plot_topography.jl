# Note: PLOT_TOPOGRAPHY_KWARGS is defined in plot_ica.jl

##########################################
# 2D topographic plot
##########################################
function _plot_topography!(fig::Figure, ax::Axis, dat::DataFrame, layout::Layout; kwargs...)

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    # actual data interpolation
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colorbar_position = pop!(plot_kwargs, :colorbar_position)
    ylim = pop!(plot_kwargs, :ylim)

    # Set title based on user preferences and data
    if plot_kwargs[:show_title]
        if plot_kwargs[:title] != ""
            ax.title = plot_kwargs[:title] # Use user-provided title
        else # get time from data 
            time_min, time_max = extrema(dat.time)
            ax.title = @sprintf("%.3f to %.3f s", time_min, time_max)
        end
        ax.titlesize = plot_kwargs[:title_fontsize]
    end

    # Compute interpolated data
    channel_data = mean.(eachcol(dat[!, layout.data.label]))
    if method == :multiquadratic || method == :inverse_multiquadratic || method == :gaussian || 
       method == :inverse_quadratic || method == :thin_plate || method == :polyharmonic ||
       method == :shepard || method == :nearest
        data = _data_interpolation_topo_multiquadratic(channel_data, layout, gridscale; rbf_type = method)
    elseif method == :spherical_spline
        data = _data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
    else
        throw(
            ArgumentError(
                "Unknown interpolation method: $method. Supported: :multiquadratic, :inverse_multiquadratic, :gaussian, :inverse_quadratic, :thin_plate, :polyharmonic, :shepard, :nearest, :spherical_spline",
            ),
        )
    end

    # Calculate ylim if not provided (must be after data is computed)
    if isnothing(ylim)
        # Make ylim symmetric around 0 for balanced topographic visualization
        data_min, data_max = extrema(data[.!isnan.(data)])
        max_abs = max(abs(data_min), abs(data_max))
        ylim = (-max_abs, max_abs)
    end

    co = contourf!(
        ax,
        range(-1.0, 1.0, length = gridscale),
        range(-1.0, 1.0, length = gridscale),
        data,
        levels = range(ylim[1], ylim[2], div(gridscale, 2));
        extendlow = :auto,
        extendhigh = :auto,
        colormap = pop!(plot_kwargs, :colormap),
        nan_color = :transparent,
    )
    co.colorrange = ylim

    # colorbar
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    if pop!(plot_kwargs, :colorbar_plot)
        Colorbar(fig[colorbar_position...], co; colorbar_kwargs...)
    end

    # head shape
    plot_layout_2d!(fig, ax, layout; plot_kwargs...)

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
function plot_topography(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot = true,
    interactive = true,
    kwargs...,
)

    set_window_title(_generate_window_title(dat))
    fig = Figure()
    ax = Axis(fig[1, 1])

    plot_topography!(
        fig,
        ax,
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        kwargs...,
    )

    # only enable interactivity for ErpData (context menu requires ErpData)
    if interactive && dat isa ErpData
        _setup_topo_interactivity!(fig, ax, dat)
    end
    display_plot && display_figure(fig)

    set_window_title("Makie")
    return fig, ax
end

function plot_topography!(
    fig,
    ax,
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    kwargs...,
)
    # Merge user kwargs with defaults to get all parameters
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    dat_subset = subset(dat, channel_selection = channel_selection, sample_selection = sample_selection)
    _plot_topography!(fig, ax, dat_subset.data, dat_subset.layout; plot_kwargs...)

end

"""
    plot_topography(dat::Vector{ErpData}; kwargs...)

Create topographic plots from a vector of ERP datasets by broadcasting across conditions.

# Arguments
- `dat`: Vector of ErpData objects (e.g., different conditions)
- `kwargs...`: Additional keyword arguments passed to plot_topography
"""
function plot_topography(
    dat::Vector{<:SingleDataFrameEeg};
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot = true,
    interactive = true,
    kwargs...,
)

    n_datasets = length(dat)
    n_datasets == 0 && @minimal_error_throw "Cannot plot empty vector of datasets"

    # Check if colorbars are enabled and get colorbar position
    kwargs_dict = Dict{Symbol,Any}(kwargs)
    colorbar_enabled = get(kwargs_dict, :colorbar_plot, true)
    user_colorbar_position = get(kwargs_dict, :colorbar_position, nothing)
    colorbar_plot_numbers = get(kwargs_dict, :colorbar_plot_numbers, [])
    dims = get(kwargs_dict, :dims, nothing)

    # Calculate grid dimensions for plots
    if isnothing(dims)
        plot_rows, plot_cols = best_rect(n_datasets)
    else
        if length(dims) != 2 || any(dims .<= 0)
            throw(ArgumentError("Invalid dimensions: $dims. Expected [rows, cols] with positive values."))
        end
        plot_rows, plot_cols = dims
        # Validate dimensions
        total_cells = plot_rows * plot_cols
        if total_cells < n_datasets
            throw(ArgumentError("Grid dimensions $dims provide $total_cells cells but need $n_datasets."))
        end
    end

    # Determine layout based on colorbar position
    if colorbar_enabled && user_colorbar_position !== nothing
        # User provided custom colorbar position - use it
        cb_row_offset, cb_col_offset = user_colorbar_position
        # Calculate total grid size needed
        # If colorbars are below (cb_row_offset > 1), we need to double the rows
        if cb_row_offset > 1
            total_rows = plot_rows * 2
            total_cols = plot_cols
        else
            # Colorbars to the right or other positions
            total_rows = plot_rows
            total_cols = plot_cols * 2
        end
    elseif colorbar_enabled
        # Default: colorbars to the right
        total_rows = plot_rows
        total_cols = plot_cols * 2
    else
        # No colorbars
        total_rows = plot_rows
        total_cols = plot_cols
    end

    # Create single figure with subplots
    set_window_title(_generate_window_title(dat))
    fig = Figure()
    axes = Axis[]

    # Plot each dataset in its own subplot
    for (idx, dataset) in enumerate(dat)
        base_row = div(idx - 1, plot_cols) + 1
        base_col = mod1(idx, plot_cols)

        if colorbar_enabled && user_colorbar_position !== nothing
            # User provided custom colorbar position
            cb_row_offset, cb_col_offset = user_colorbar_position
            if cb_row_offset > 1
                # Colorbars below: each plot needs 2 rows
                plot_row = (base_row - 1) * 2 + 1
                plot_col = base_col
                colorbar_row = plot_row + (cb_row_offset - 1)
                colorbar_col = plot_col + (cb_col_offset - 1)
            else
                # Colorbars to the right or other positions
                plot_row = base_row
                plot_col = (base_col - 1) * 2 + 1
                colorbar_row = plot_row + (cb_row_offset - 1)
                colorbar_col = plot_col + (cb_col_offset - 1)
            end
        elseif colorbar_enabled
            # Default: colorbar to the right - each subplot gets 2 columns
            plot_row = base_row
            plot_col = (base_col - 1) * 2 + 1
            colorbar_row = plot_row
            colorbar_col = plot_col + 1
        else
            # No colorbars
            plot_row = base_row
            plot_col = base_col
        end

        ax = Axis(fig[plot_row, plot_col])
        push!(axes, ax)

        # Set subplot title to condition name if available
        if hasproperty(dataset, :condition_name) && dataset.condition_name !== ""
            ax.title = dataset.condition_name
        end

        # Prepare kwargs for this subplot
        subplot_kwargs = copy(kwargs_dict)
        # Determine if this dataset should have a colorbar
        # If colorbar_plot_numbers is empty, show colorbar for all datasets
        # Otherwise, only show for datasets whose index (1-based) is in the list
        should_show_colorbar = colorbar_enabled && (isempty(colorbar_plot_numbers) || idx in colorbar_plot_numbers)
        if should_show_colorbar
            subplot_kwargs[:colorbar_position] = (colorbar_row, colorbar_col)
        else
            # Disable colorbar for this specific dataset
            subplot_kwargs[:colorbar_plot] = false
        end

        # Plot the topography in this subplot
        plot_topography!(
            fig,
            ax,
            dataset;
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            subplot_kwargs...,
        )
    end

    # Set column sizes only if colorbars are enabled and to the right (default)
    if colorbar_enabled && (user_colorbar_position === nothing || user_colorbar_position[1] <= 1)
        # Make colorbar columns narrower than plot columns
        # This ensures colorbars don't take up too much space while keeping plots visible
        for col = 1:total_cols
            if col % 2 == 1
                colsize!(fig.layout, col, Auto())
            else
                colsize!(fig.layout, col, Fixed(50))
            end
        end
    end

    # Only enable interactivity if all datasets are ErpData (context menu requires ErpData)
    if interactive && all(d isa ErpData for d in dat)
        shared_selection_state = TopoSelectionState(axes)
        _setup_shared_topo_interactivity!(fig, axes, dat, shared_selection_state)
    end

    display_plot && display_figure(fig)

    set_window_title("Makie")
    return fig, axes
end

function plot_topography(dat::Vector{EpochData}, epoch::Int; kwargs...)
    @info "Plotting epoch $(epoch) for each dataset in the vector"
    plot_topography.(dat, Ref(epoch))
end


function plot_topography!(
    fig,
    ax,
    dat::MultiDataFrameEeg,
    epoch::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    kwargs...,
)
    plot_topography!(
        fig,
        ax,
        convert(dat, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        kwargs...,
    )
end

function plot_topography(
    dat::MultiDataFrameEeg,
    epoch::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot = true,
    kwargs...,
)

    set_window_title(_generate_window_title(dat) * " - Epoch $epoch")
    fig = Figure()
    ax = Axis(fig[1, 1])

    plot_topography!(
        fig,
        ax,
        convert(dat, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        kwargs...,
    )

    display_plot && display_figure(fig)

    set_window_title("Makie")
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

    center = (grid_scale + 1) / 2
    radius_squared = (grid_scale / 2)^2

    @inbounds for col = 1:grid_scale
        x_dist_squared = (center - col)^2
        for row = 1:grid_scale
            y_dist_squared = (center - row)^2
            if x_dist_squared + y_dist_squared > radius_squared
                dat[col, row] = NaN
            end
        end
    end
    return dat
end



"""
    _data_interpolation_topo_multiquadratic(dat::Vector{<:AbstractFloat}, layout::Layout, grid_scale::Int; rbf_type::Symbol=:multiquadratic)

Interpolate EEG data using scattered interpolation with various radial basis functions.

For detailed information about these interpolation methods, see:
https://eljungsk.github.io/ScatteredInterpolation.jl/dev/methods/

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `layout::Layout`: Layout containing electrode 2D coordinates
- `grid_scale::Int`: Size of the output grid
- `rbf_type::Symbol`: Type of radial basis function or interpolation method. Options:
  - `:multiquadratic` - Multiquadratic (default): φ(r) = √(1 + (εr)²)
  - `:inverse_multiquadratic` - Inverse Multiquadratic: φ(r) = 1/√(1 + (εr)²)
  - `:gaussian` - Gaussian: φ(r) = exp(-(εr)²)
  - `:inverse_quadratic` - Inverse Quadratic: φ(r) = 1/(1 + (εr)²)
  - `:thin_plate` - Thin Plate Spline: φ(r) = r²ln(r)
  - `:polyharmonic` - Polyharmonic Spline (k=3): φ(r) = r³
  - `:shepard` - Inverse Distance Weighting (Shepard interpolation): wᵢ(x) = 1/||x-xᵢ||ᵖ
  - `:nearest` - Nearest Neighbor: returns value of nearest sample point

# Returns
- `Matrix{Float64}`: Interpolated data on a regular grid
"""
function _data_interpolation_topo_multiquadratic(dat::Vector{<:AbstractFloat}, layout::Layout, grid_scale::Int; rbf_type::Symbol=:multiquadratic)

    if any(isnan, dat) || any(isinf, dat)
        throw(ArgumentError("Input data contains NaN or Inf values"))
    end

    points = permutedims(Matrix(layout.data[!, [:x2, :y2]]))
    # Use normalized coordinate ranges since layout coordinates are normalized to [-1, 1]
    x_range = range(-1.0, 1.0, length = grid_scale)
    y_range = range(-1.0, 1.0, length = grid_scale)

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

    # Select radial basis function type or interpolation method
    rbf = if rbf_type == :multiquadratic
        ScatteredInterpolation.Multiquadratic()
    elseif rbf_type == :inverse_multiquadratic
        ScatteredInterpolation.InverseMultiquadratic()
    elseif rbf_type == :gaussian
        ScatteredInterpolation.Gaussian()
    elseif rbf_type == :inverse_quadratic
        ScatteredInterpolation.InverseQuadratic()
    elseif rbf_type == :thin_plate
        ScatteredInterpolation.ThinPlate()
    elseif rbf_type == :polyharmonic
        ScatteredInterpolation.Polyharmonic(3)  # k=3 for r³
    elseif rbf_type == :shepard
        ScatteredInterpolation.Shepard()  # Inverse Distance Weighting
    elseif rbf_type == :nearest
        ScatteredInterpolation.NearestNeighbor()
    else
        throw(ArgumentError("Unknown RBF type: $rbf_type. Supported: :multiquadratic, :inverse_multiquadratic, :gaussian, :inverse_quadratic, :thin_plate, :polyharmonic, :shepard, :nearest"))
    end

    itp = ScatteredInterpolation.interpolate(rbf, points, dat)
    result = reshape(ScatteredInterpolation.evaluate(itp, grid_points), grid_scale, grid_scale)
    _circle_mask!(result, grid_scale)
    return result

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

    # Ensure 3D coordinates exist
    _ensure_coordinates_3d!(layout)

    # Extract 3D coordinates efficiently from DataFrame
    n_channels = length(dat)
    x3_col = layout.data.x3::Vector{Float64}
    y3_col = layout.data.y3::Vector{Float64}
    z3_col = layout.data.z3::Vector{Float64}

    coords = Matrix{Float64}(undef, n_channels, 3)
    @inbounds for i = 1:n_channels
        coords[i, 1] = x3_col[i]
        coords[i, 2] = y3_col[i]
        coords[i, 3] = z3_col[i]
    end

    # Find the actual radius of the electrode positions
    electrode_radius = 0.0
    @inbounds for i = 1:n_channels
        electrode_radius += sqrt(coords[i, 1]^2 + coords[i, 2]^2 + coords[i, 3]^2)
    end
    electrode_radius /= n_channels

    # Pre-normalize electrode coordinates to unit sphere
    coords_unit = Matrix{Float64}(undef, n_channels, 3)
    @inbounds for i = 1:n_channels
        norm_factor = sqrt(coords[i, 1]^2 + coords[i, 2]^2 + coords[i, 3]^2)
        if norm_factor > 0
            coords_unit[i, 1] = coords[i, 1] / norm_factor
            coords_unit[i, 2] = coords[i, 2] / norm_factor
            coords_unit[i, 3] = coords[i, 3] / norm_factor
        else
            coords_unit[i, 1], coords_unit[i, 2], coords_unit[i, 3] = coords[i, 1], coords[i, 2], coords[i, 3]
        end
    end

    # Step 1: Solve for weights
    cosang = coords_unit * coords_unit'
    G = _calc_g_matrix(cosang)

    G_extended = Matrix{Float64}(undef, n_channels + 1, n_channels + 1)
    @inbounds for i = 1:n_channels
        for j = 1:n_channels
            G_extended[i, j] = G[i, j]
        end
        G_extended[i, i] += 1e-5 # λ regularization
        G_extended[i, n_channels+1] = 1.0
        G_extended[n_channels+1, i] = 1.0
    end
    G_extended[n_channels+1, n_channels+1] = 0.0

    data_vector = Vector{Float64}(undef, n_channels + 1)
    @inbounds for i = 1:n_channels
        data_vector[i] = dat[i]
    end
    data_vector[n_channels+1] = 0.0

    weights = G_extended \ data_vector

    # Step 2: Prepare grid points for interpolation
    x_range = range(-1.0, 1.0, length = grid_scale)
    y_range = range(-1.0, 1.0, length = grid_scale)

    # Find valid indices first to minimize work
    valid_indices = Tuple{Int,Int}[]
    sizehint!(valid_indices, grid_scale^2)
    for (iy, y) in enumerate(y_range), (ix, x) in enumerate(x_range)
        if x^2 + y^2 <= 1.0
            push!(valid_indices, (ix, iy))
        end
    end

    n_valid = length(valid_indices)
    grid_points_unit = Matrix{Float64}(undef, n_valid, 3)

    @inbounds for i = 1:n_valid
        ix, iy = valid_indices[i]
        x, y = x_range[ix], y_range[iy]
        r2 = x^2 + y^2
        r = sqrt(r2)

        # Stereographic projection
        z3 = electrode_radius * (1.0 - r^2) / (1.0 + r^2)
        px = x * (1.0 + z3 / electrode_radius)
        py = y * (1.0 + z3 / electrode_radius)
        pz = z3

        # Normalize to unit sphere
        p_norm = sqrt(px^2 + py^2 + pz^2)
        if p_norm > 0
            grid_points_unit[i, 1] = px / p_norm
            grid_points_unit[i, 2] = py / p_norm
            grid_points_unit[i, 3] = pz / p_norm
        else
            grid_points_unit[i, 1], grid_points_unit[i, 2], grid_points_unit[i, 3] = px, py, pz
        end
    end

    # Step 3: Compute G matrix between valid grid points and channels using BLAS
    cosang_grid = grid_points_unit * coords_unit'

    # Step 4: Apply g-function in-place with pre-calculated factors
    factors = _get_g_factors(15)
    @inbounds for i in eachindex(cosang_grid)
        cosang_grid[i] = _legendre_val(clamp(cosang_grid[i], -1.0, 1.0), factors)
    end

    # Step 5: Compute final interpolation values using BLAS
    w_channels = @view weights[1:n_channels]
    interpolated_valid = cosang_grid * w_channels
    interpolated_valid .+= weights[n_channels+1]

    # Step 6: Map back to 2D grid
    result = fill(NaN, grid_scale, grid_scale)
    @inbounds for i = 1:n_valid
        ix, iy = valid_indices[i]
        result[ix, iy] = interpolated_valid[i]
    end

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

# Optimized Legendre polynomial series evaluation using recurrence
# Pre-calculate recurrence coefficients to avoid divisions in the inner loop
const _LEGENDRE_COEFF_A = Float64[(2n - 1) / n for n = 0:17]
const _LEGENDRE_COEFF_B = Float64[(n - 1) / n for n = 0:17]

@inline function _legendre_val(x::Float64, factors::Vector{Float64})
    """Evaluate Legendre polynomial series efficiently using recurrence relation.

    factors[1] = 0.0 (for n=0, skipped)
    factors[2] = factor for n=1
    factors[3] = factor for n=2
    etc.
    """
    max_order = length(factors) - 1
    if max_order < 1
        return 0.0
    end

    # Initialize recurrence: P_0 = 1, P_1 = x
    p_prev = 1.0  # P_0
    p_curr = x     # P_1

    # Start with n=1 term: factors[2] * P_1(x) = factors[2] * x
    @inbounds result = factors[2] * p_curr

    # Compute remaining polynomials using recurrence (n=2 to max_order)
    # Recurrence: P_n(x) = ((2n-1)/n * x * P_{n-1}(x) - (n-1)/n * P_{n-2}(x))
    @inbounds for n = 2:max_order
        p_next = _LEGENDRE_COEFF_A[n+1] * x * p_curr - _LEGENDRE_COEFF_B[n+1] * p_prev
        result += factors[n+1] * p_next
        p_prev = p_curr
        p_curr = p_next
    end

    return result
end

# Pre-computed factors for m=4, cached to avoid recomputation
const _G_FACTORS_CACHE = Dict{Int,Vector{Float64}}()

function _get_g_factors(n_legendre_terms::Int = 15)
    """Get pre-computed factors for g-function calculation (cached)."""
    if !haskey(_G_FACTORS_CACHE, n_legendre_terms)
        factors = [(2 * n + 1) / (n^4 * (n + 1)^4 * 4 * π) for n = 1:n_legendre_terms]
        _G_FACTORS_CACHE[n_legendre_terms] = [0.0; factors]  # Prepend 0.0 for n=0
    end
    return _G_FACTORS_CACHE[n_legendre_terms]::Vector{Float64}
end

# MNE-Python's exact g-function calculation for EEG topography (m=4)
function _calc_g_function(cosang::Float64, n_legendre_terms::Int = 15)
    """Calculate spherical spline g function between points on a sphere.

    This is the exact implementation from MNE-Python, optimized for EEG topography (m=4).
    """
    cosang ≈ 1.0 && return 0.0

    # Use cached factors and Legendre evaluation
    factors = _get_g_factors(n_legendre_terms)
    return _legendre_val(cosang, factors)
end

# MNE-Python's exact G matrix calculation for EEG topography (m=4)
function _calc_g_matrix(cosang::Matrix{Float64}, n_legendre_terms::Int = 15)
    """Calculate spherical spline G matrix between points on a sphere.

    This is the exact implementation from MNE-Python, optimized for EEG topography (m=4).
    """
    factors = _get_g_factors(n_legendre_terms)

    # Use Legendre polynomial evaluation for the entire matrix
    G = similar(cosang)
    @inbounds for i in eachindex(cosang)
        G[i] = _legendre_val(cosang[i], factors)
    end
    return G
end

# =============================================================================
# INTERACTIVITY FUNCTIONS
# =============================================================================

"""
    _setup_topo_keyboard_handlers!(fig::Figure, axes::Union{Axis, Vector{Axis}})

Set up keyboard event handlers for topographic plots.
Handles both single axis and multiple axes.
"""
function _setup_topo_keyboard_handlers!(fig::Figure, axes::Union{Axis,Vector{Axis}})
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.i
                show_plot_help(:topography)
            elseif event.key == Keyboard.up
                axes isa Vector ? _topo_scale_up!.(axes) : _topo_scale_up!(axes)
            elseif event.key == Keyboard.down
                axes isa Vector ? _topo_scale_down!.(axes) : _topo_scale_down!(axes)
            end
        end
    end
end

"""
    _setup_topo_interactivity!(fig::Figure, ax::Axis)

Set up keyboard interactivity for topographic plots.
"""
function _setup_topo_interactivity!(fig::Figure, ax::Axis, original_data = nothing)
    deregister_interaction!(ax, :rectanglezoom)
    _setup_topo_keyboard_handlers!(fig, ax)
    _setup_topo_selection!(fig, ax, original_data)
end

# =============================================================================
# TOPO SELECTION STATE
# =============================================================================

"""
    TopoSelectionState

Simple state for spatial region selection in topographic plots.
"""
mutable struct TopoSelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64,Float64,Float64}}  # x_min, y_min, x_max, y_max
    visible::Observable{Bool}
    rectangles::Vector{Vector{Makie.Poly}}  # Store rectangles for each axis: rectangles[axis_idx][rect_idx]
    bounds_list::Observable{Vector{Tuple{Float64,Float64,Float64,Float64}}}  # Store all selection bounds
    temp_rectangles::Vector{Union{Makie.Poly,Nothing}}  # Temporary rectangles for each axis
    selected_channels::Vector{Symbol}  # Store selected channel names
    axes::Vector{Axis}  # All axes that share this selection state

    function TopoSelectionState(axes::Vector{Axis})
        n_axes = length(axes)
        # Initialize with empty lists for multiple selections
        new(
            Observable(false),
            Observable((0.0, 0.0, 0.0, 0.0)),
            Observable(false),
            [Makie.Poly[] for _ = 1:n_axes],  # Empty vectors for rectangles per axis
            Observable{Tuple{Float64,Float64,Float64,Float64}}[],  # Empty vector for bounds
            [nothing for _ = 1:n_axes],  # No temporary rectangles initially
            Symbol[],  # Empty vector for selected channels
            axes,  # Store all axes
        )
    end
end

"""
    _setup_shared_topo_interactivity!(fig::Figure, axes::Vector{Axis}, datasets::Vector, shared_selection_state)

Set up shared interactivity for multiple topographic plots.
"""
function _setup_shared_topo_interactivity!(
    fig::Figure,
    axes::Vector{Axis},
    datasets::Vector,
    shared_selection_state::TopoSelectionState,
)
    deregister_interaction!.(axes, :rectanglezoom)
    _setup_topo_keyboard_handlers!(fig, axes)
    _setup_shared_topo_selection!(fig, datasets, shared_selection_state)
end

"""
    _scale_topo_levels!(ax::Axis, scale_factor::Float64)

Scale the topographic plot levels by the given factor.
- scale_factor < 1.0: zoom in (compress range)
- scale_factor > 1.0: zoom out (expand range)
"""
function _scale_topo_levels!(ax::Axis, scale_factor::Float64)
    # Find the Contourf plot in the axis
    for plot in ax.scene.plots
        current_levels = plot.levels[]
        if !isnothing(current_levels)
            # Calculate new range while keeping it centered
            level_min, level_max = extrema(current_levels)
            center = (level_min + level_max) / 2
            range_size = level_max - level_min

            # Scale the range by the factor but keep it centered
            new_range = range_size * scale_factor
            new_min = center - new_range / 2
            new_max = center + new_range / 2

            # Create new levels with the same density but scaled range
            plot.levels[] = range(new_min, new_max, length = length(current_levels))
            break
        end
    end
end

"""
    _topo_scale_up!(ax::Axis)

Increase the scale of the topographic plot (zoom in on color range).
"""
_topo_scale_up!(ax::Axis) = _scale_topo_levels!(ax, 0.8)  # Compress range by 20%

"""
    _topo_scale_down!(ax::Axis)

Decrease the scale of the topographic plot (zoom out from color range).
"""
_topo_scale_down!(ax::Axis) = _scale_topo_levels!(ax, 1.25)  # Expand range by 25%

# =============================================================================
# REGION SELECTION FOR TOPO PLOTS
# =============================================================================


"""
    _create_rectangles_on_all_axes!(selection_state::TopoSelectionState, rect_points::Vector{Point2f}, 
                                     is_temporary::Bool)

Create rectangles on all axes with the given points.
- If `is_temporary`, stores in `temp_rectangles`
- Otherwise, creates permanent rectangles and stores in `rectangles`
"""
function _create_rectangles_on_all_axes!(
    selection_state::TopoSelectionState,
    rect_points::Vector{Point2f},
    is_temporary::Bool,
)
    for (idx, other_ax) in enumerate(selection_state.axes)
        rect = poly!(
            other_ax,
            rect_points,
            color = (:blue, 0.3),
            strokecolor = :black,
            strokewidth = 1,
            visible = true,
            overdraw = true,
        )

        if is_temporary
            selection_state.temp_rectangles[idx] = rect
        else
            push!(selection_state.rectangles[idx], rect)
        end
    end
end

"""
    _setup_topo_selection!(fig::Figure, ax::Axis, original_data)

Set up simple region selection for topographic plots.
This is a convenience wrapper for single-axis plots that calls the shared version.
"""
function _setup_topo_selection!(fig::Figure, ax::Axis, original_data)
    selection_state = TopoSelectionState([ax])
    _setup_shared_topo_selection!(fig, [original_data], selection_state)
end

"""
    _setup_shared_topo_selection!(fig::Figure, datasets::Vector, shared_selection_state)

Set up shared region selection for multiple topographic plots.
Event handlers are set up only once for all axes.
"""
function _setup_shared_topo_selection!(fig::Figure, datasets::Vector, shared_selection_state::TopoSelectionState)
    # Track Shift key state
    shift_pressed = Observable(false)

    # Handle keyboard events to track Shift key
    on(events(fig).keyboardbutton) do event
        if event.key == Keyboard.left_shift
            if event.action == Keyboard.press
                shift_pressed[] = true
            elseif event.action == Keyboard.release
                shift_pressed[] = false
            end
        end
    end

    # Handle mouse events at the figure level
    on(events(fig).mousebutton) do event
        # Check if mouse is over any of the axes in the shared state
        mouse_pos = events(fig).mouseposition[]
        active_ax, active_dataset = _find_active_axis_with_dataset(shared_selection_state.axes, mouse_pos, datasets)

        # Only process if mouse is over one of the shared axes
        active_ax === nothing && return

        if event.button == Mouse.left
            if event.action == Mouse.press
                if shift_pressed[]
                    _start_topo_selection!(active_ax, shared_selection_state)
                else
                    _clear_all_topo_selections!(shared_selection_state)
                end
            elseif event.action == Mouse.release
                if shared_selection_state.active[]
                    _finish_topo_selection!(active_ax, shared_selection_state, active_dataset)
                end
            end
        end
        if event.button == Mouse.right && event.action == Mouse.press
            selected_channels = shared_selection_state.selected_channels
            if !isempty(selected_channels)
                # Interactivity is only enabled for ErpData, so we can safely call the context menu
                _show_topo_context_menu!(datasets, selected_channels)
            else
                @minimal_warning "No channels selected. Select channels with Shift+Left Click+Drag, then right-click to plot ERP"
            end
        end
    end

    # Handle mouse movement for selection dragging
    on(events(fig).mouseposition) do pos
        # Only update if selection is active and mouse is over any shared axis
        if shared_selection_state.active[]
            active_ax = _find_active_axis(shared_selection_state.axes, pos)
            if active_ax !== nothing
                _update_topo_selection!(active_ax, shared_selection_state)
            end
        end
    end
end

"""
    _start_topo_selection!(ax::Axis, selection_state::TopoSelectionState)

Start spatial region selection in topographic plot.
"""
function _start_topo_selection!(ax::Axis, selection_state::TopoSelectionState)
    selection_state.active[] = true
    selection_state.visible[] = true

    # Get mouse position in axis coordinates
    mouse_pos = mouseposition(ax)
    mouse_x, mouse_y = mouse_pos[1], mouse_pos[2]

    # Store axis coordinates for spatial selection
    selection_state.bounds[] = (mouse_x, mouse_y, mouse_x, mouse_y)

    # Create temporary rectangles for all axes (all points same initially)
    initial_points = [Point2f(mouse_x, mouse_y) for _ = 1:4]
    _create_rectangles_on_all_axes!(selection_state, initial_points, true)

    _update_topo_selection!(ax, selection_state)
end

"""
    _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, original_data)

Finish spatial region selection in topographic plot.
"""
function _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, original_data = nothing)
    # Get final mouse position in axis coordinates
    mouse_pos = mouseposition(ax)
    mouse_x, mouse_y = mouse_pos[1], mouse_pos[2]

    # Get start position
    start_x, start_y = selection_state.bounds[][1], selection_state.bounds[][2]

    # Update bounds with final position (x_min, y_min, x_max, y_max) in axis coords
    x_min, x_max = minmax(start_x, mouse_x)
    y_min, y_max = minmax(start_y, mouse_y)
    final_bounds = (x_min, y_min, x_max, y_max)
    selection_state.bounds[] = final_bounds

    # Store this selection in the bounds list
    current_bounds = selection_state.bounds_list[]
    push!(current_bounds, final_bounds)
    selection_state.bounds_list[] = current_bounds

    # Create permanent rectangles on ALL axes with the same bounds
    start_x, start_y, end_x, end_y = final_bounds

    rect_points = Point2f[
        Point2f(Float64(start_x), Float64(start_y)),
        Point2f(Float64(end_x), Float64(start_y)),
        Point2f(Float64(end_x), Float64(end_y)),
        Point2f(Float64(start_x), Float64(end_y)),
    ]

    # Create permanent rectangles on all axes
    _create_rectangles_on_all_axes!(selection_state, rect_points, false)

    # Remove the temporary rectangles from all axes
    for (idx, temp_rect) in enumerate(selection_state.temp_rectangles)
        if !isnothing(temp_rect)
            delete!(temp_rect.parent, temp_rect)
            selection_state.temp_rectangles[idx] = nothing
        end
    end

    # Find electrodes within ALL selected spatial regions
    all_selected_electrodes = Symbol[]
    for bounds in selection_state.bounds_list[]
        x_min, y_min, x_max, y_max = bounds
        region_electrodes = _find_electrodes_in_region(x_min, y_min, x_max, y_max, original_data)
        append!(all_selected_electrodes, region_electrodes)
    end

    unique_electrodes = unique(all_selected_electrodes)
    @info "$(length(selection_state.bounds_list[])) regions; Channels: $unique_electrodes"

    # Store selected channels in the state
    selection_state.selected_channels = unique_electrodes

    # Reset active state
    selection_state.active[] = false
end

"""
    _update_topo_selection!(ax::Axis, selection_state::TopoSelectionState)

Update the visual selection rectangle for spatial selection.
"""
function _update_topo_selection!(ax::Axis, selection_state::TopoSelectionState)
    if selection_state.active[]
        # Get current mouse position in axis coordinates
        axis_pos = mouseposition(ax)
        start_x, start_y = selection_state.bounds[][1], selection_state.bounds[][2]

        # Update bounds with the axis coordinates
        end_x, end_y = axis_pos[1], axis_pos[2]
        selection_state.bounds[] = (start_x, start_y, end_x, end_y)

        # Update rectangle points for the temporary rectangles
        rect_points = Point2f[
            Point2f(Float64(start_x), Float64(start_y)),
            Point2f(Float64(end_x), Float64(start_y)),
            Point2f(Float64(end_x), Float64(end_y)),
            Point2f(Float64(start_x), Float64(end_y)),
        ]

        # Update temporary rectangles on all axes
        for temp_rect in selection_state.temp_rectangles
            if !isnothing(temp_rect)
                temp_rect[1] = rect_points
            end
        end

    end
end

"""
    _clear_all_topo_selections!(selection_state::TopoSelectionState)

Clear all topographic selections and remove all rectangles.
"""
function _clear_all_topo_selections!(selection_state::TopoSelectionState)
    # remove all rectangles from all axes
    for axis_rects in selection_state.rectangles
        for rect in axis_rects
            delete!(rect.parent, rect)
        end
        empty!(axis_rects)  # Clear the list after deleting
    end

    # clear the state data
    selection_state.bounds_list[] = Tuple{Float64,Float64,Float64,Float64}[]
    empty!(selection_state.selected_channels)

    # reset state observables
    selection_state.active[] = false
    selection_state.visible[] = false
end


"""
    _find_electrodes_in_region(x_min::Float64, y_min::Float64, x_max::Float64, y_max::Float64, original_data)

Find electrodes within the selected spatial region using actual electrode coordinates.
This approach uses the real layout data from the topographic plot.
"""
function _find_electrodes_in_region(x_min::Float64, y_min::Float64, x_max::Float64, y_max::Float64, original_data)
    # Filter electrodes that are inside the selection rectangle
    selected_rows =
        Base.filter(row -> x_min <= row.x2 <= x_max && y_min <= row.y2 <= y_max, eachrow(original_data.layout.data))
    return [Symbol(row.label) for row in selected_rows]
end

"""
    _show_topo_context_menu!(datasets::Union{ErpData, Vector{ErpData}}, selected_channels::Vector{Symbol})

Show a context menu for plotting selected channels from topography plot.
Supports both single dataset and multiple datasets (conditions).
"""
function _show_topo_context_menu!(datasets::Union{ErpData,Vector{ErpData}}, selected_channels::Vector{Symbol})

    datasets_vec = datasets isa Vector ? datasets : [datasets]
    has_multiple_conditions = length(datasets_vec) > 1
    has_multiple_channels = length(selected_channels) > 1

    plot_types = String[]
    plot_configs = Tuple{Bool,Bool,String}[]

    if has_multiple_conditions && has_multiple_channels
        # Multiple conditions and multiple channels: show all 4 options
        push!(plot_types, "Separate channels, separate conditions")
        push!(plot_types, "Average channels, separate conditions")
        push!(plot_types, "Separate channels, average conditions")
        push!(plot_types, "Average channels, average conditions")
        plot_configs = [
            (false, false, "separate channels, separate conditions"),
            (false, true, "average channels, separate conditions"),
            (true, false, "separate channels, average conditions"),
            (true, true, "average channels, average conditions"),
        ]
    elseif has_multiple_conditions
        # Multiple conditions but single channel: only condition averaging options
        push!(plot_types, "Separate conditions")
        push!(plot_types, "Average conditions")
        plot_configs = [(false, false, "separate conditions"), (true, false, "average conditions")]
    elseif has_multiple_channels
        # Single condition but multiple channels: only channel averaging options
        push!(plot_types, "Plot Individual Channels")
        push!(plot_types, "Plot Averaged Channels")
        plot_configs = [(false, false, "individual channels"), (false, true, "averaged channels")]
    else
        # Single condition and single channel: just plot it
        push!(plot_types, "Plot Channel")
        plot_configs = [(false, false, "single channel")]
    end

    menu_fig = Figure(size = (400, 200))
    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for (idx, btn) in enumerate(menu_buttons)
        on(btn.clicks) do n
            avg_conditions, avg_channels, msg = plot_configs[idx]

            # Prepare data: average conditions if needed
            data_to_plot = avg_conditions ? _average_conditions(datasets_vec) : datasets_vec
            @info "Plotting ERP: $msg: $selected_channels"
            plot_erp(data_to_plot; channel_selection = channels(selected_channels), average_channels = avg_channels)
        end
    end

    new_screen = GLMakie.Screen()
    display(new_screen, menu_fig)
end

##########################################
# 3D topographic plot
##########################################

"""
    _read_obj_file(filepath::String) -> (vertices::Matrix{Float64}, faces::Matrix{Int})

Read a 3D mesh from an OBJ file.

# Arguments
- `filepath::String`: Path to the OBJ file

# Returns
- `vertices::Matrix{Float64}`: [n_vertices × 3] matrix of vertex coordinates (x, y, z)
- `faces::Matrix{Int}`: [n_faces × 3] matrix of triangle face indices (1-indexed)

# Notes
- Only supports triangular faces (f v1 v2 v3)
- Ignores texture coordinates and normals
- Vertex indices in OBJ are 1-indexed, returned as 1-indexed
"""
function _read_obj_file(filepath::String)
    vertices = Vector{Vector{Float64}}()
    faces = Vector{Vector{Int}}()
    
    open(filepath, "r") do file
        for line in eachline(file)
            # Remove comments
            line = strip(split(line, '#')[1])
            isempty(line) && continue
            
            parts = split(line)
            if length(parts) == 0
                continue
            end
            
            if parts[1] == "v" && length(parts) >= 4
                # Vertex: v x y z [w]
                x = parse(Float64, parts[2])
                y = parse(Float64, parts[3])
                z = parse(Float64, parts[4])
                push!(vertices, [x, y, z])
            elseif parts[1] == "f" && length(parts) >= 4
                # Face: f v1 v2 v3 [v4 ...]
                # Handle formats: f v1 v2 v3 or f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
                face_verts = Int[]
                for i = 2:length(parts)  # Read all vertices in the face
                    vert_str = split(parts[i], '/')[1]  # Get vertex index, ignore texture/normal
                    vert_idx = parse(Int, vert_str)
                    # OBJ uses 1-indexing, but negative indices are relative to end
                    if vert_idx < 0
                        vert_idx = length(vertices) + vert_idx + 1
                    end
                    push!(face_verts, vert_idx)
                end
                # Triangulate faces: quads become 2 triangles, n-gons become n-2 triangles
                if length(face_verts) == 3
                    # Already a triangle
                    push!(faces, face_verts)
                elseif length(face_verts) == 4
                    # Quad: split into two triangles (fan triangulation)
                    push!(faces, [face_verts[1], face_verts[2], face_verts[3]])
                    push!(faces, [face_verts[1], face_verts[3], face_verts[4]])
                elseif length(face_verts) > 4
                    # N-gon: fan triangulation from first vertex
                    for i = 2:(length(face_verts) - 1)
                        push!(faces, [face_verts[1], face_verts[i], face_verts[i + 1]])
                    end
                end
            end
        end
    end
    
    # Convert to matrices
    n_vertices = length(vertices)
    n_faces = length(faces)
    
    if n_vertices == 0
        throw(ArgumentError("OBJ file contains no vertices"))
    end
    
    vertices_matrix = Matrix{Float64}(undef, n_vertices, 3)
    @inbounds for i = 1:n_vertices
        vertices_matrix[i, :] = vertices[i]
    end
    
    # Normalize vertices to fit within unit sphere (scale to [-1, 1] range)
    # Find bounding box
    min_coords = [minimum(vertices_matrix[:, i]) for i = 1:3]
    max_coords = [maximum(vertices_matrix[:, i]) for i = 1:3]
    center = [(min_coords[i] + max_coords[i]) / 2.0 for i = 1:3]
    ranges = [max_coords[i] - min_coords[i] for i = 1:3]
    max_range = maximum(ranges)
    
    if max_range > 0
        # Center and scale to fit in unit sphere
        @inbounds for i = 1:n_vertices
            for j = 1:3
                vertices_matrix[i, j] = (vertices_matrix[i, j] - center[j]) / (max_range / 2.0)
            end
        end
    end
    
    faces_matrix = Matrix{Int}(undef, n_faces, 3)
    @inbounds for i = 1:n_faces
        faces_matrix[i, :] = faces[i]
    end
    
    return vertices_matrix, faces_matrix
end

"""
    _generate_parametric_head_mesh(radius::Float64=1.0, n_points::Int=50)

Generate a parametric 3D head mesh (ellipsoid with simple nose and ears).

# Arguments
- `radius::Float64`: Base radius of the head (default: 1.0)
- `n_points::Int`: Number of points along each angular dimension (default: 50)

# Returns
- `vertices::Matrix{Float64}`: [n_points × 3] matrix of vertex coordinates
- `faces::Matrix{Int}`: [n_faces × 3] matrix of triangle face indices
- `vertex_normals::Matrix{Float64}`: [n_points × 3] matrix of vertex normals
"""
function _generate_parametric_head_mesh(radius::Float64 = 1.0, n_points::Int = 50)
    # Generate spherical coordinates
    # For EEG: z is up (top of head), x is forward (nose), y is left-right (ears)
    # theta: azimuth angle in x-y plane (0 = +x/nose, π/2 = +y/right ear, π = -x/back, 3π/2 = -y/left ear)
    # phi: elevation from z-axis (0 = top of head, π/2 = equator, π = bottom)
    theta = range(0, 2π, length = n_points)
    phi = range(0, π, length = n_points)

    n_vertices = n_points * n_points
    vertices = Matrix{Float64}(undef, n_vertices, 3)
    vertex_normals = Matrix{Float64}(undef, n_vertices, 3)

    # Generate vertices with clearly head-like shape
    # Store as grid: phi (rows) × theta (cols) - phi varies first (outer loop), then theta
    idx = 1
    for (i, p) in enumerate(phi), (j, t) in enumerate(theta)
        # Base ellipsoid - head is wider than it is tall
        # Standard spherical coordinates: x = r*sin(phi)*cos(theta), y = r*sin(phi)*sin(theta), z = r*cos(phi)
        x_base = radius * sin(p) * cos(t)
        y_base = radius * sin(p) * sin(t)
        z_base = radius * cos(p)
        
        # Make head shape: slightly wider horizontally, taller vertically
        x_base *= 0.92  # Narrower front-to-back
        y_base *= 0.92  # Narrower left-to-right  
        z_base *= 1.2   # Taller (top to bottom)

        # Add nose protrusion (front = +x direction, theta ≈ 0)
        nose_factor = 0.0
        # Front is where theta ≈ 0 (cos(theta) ≈ 1, sin(theta) ≈ 0)
        if abs(t) < 0.5 || abs(t - 2π) < 0.5  # Front region
            # Nose should be in upper-middle face
            if p > 0.4 && p < 1.0
                nose_strength = exp(-(min(abs(t), abs(t - 2π)) / 0.3)^2)  # Gaussian falloff
                vertical_strength = sin((p - 0.4) / 0.6 * π)
                nose_factor = 0.35 * radius * nose_strength * vertical_strength
            end
        end
        # Nose extends forward in +x direction
        x = x_base + nose_factor
        y = y_base
        z = z_base

        # Add ear protrusions (sides = ±y direction, theta ≈ π/2 and 3π/2)
        ear_factor = 0.0
        # Right ear: theta ≈ π/2 (y positive)
        if abs(t - π/2) < 0.6
            if p > 0.5 && p < 1.3
                ear_strength = exp(-((t - π/2) / 0.4)^2)
                vertical_strength = sin((p - 0.5) / 0.8 * π)
                ear_factor = 0.4 * radius * ear_strength * vertical_strength
            end
            y = y + ear_factor  # Extend in +y (right)
        # Left ear: theta ≈ 3π/2 (y negative)
        elseif abs(t - 3π/2) < 0.6
            if p > 0.5 && p < 1.3
                ear_strength = exp(-((t - 3π/2) / 0.4)^2)
                vertical_strength = sin((p - 0.5) / 0.8 * π)
                ear_factor = 0.4 * radius * ear_strength * vertical_strength
            end
            y = y - ear_factor  # Extend in -y (left)
        end

        vertices[idx, 1] = x
        vertices[idx, 2] = y
        vertices[idx, 3] = z

        # Compute normal (pointing outward)
        norm = sqrt(x^2 + y^2 + z^2)
        if norm > 0
            vertex_normals[idx, 1] = x / norm
            vertex_normals[idx, 2] = y / norm
            vertex_normals[idx, 3] = z / norm
        else
            vertex_normals[idx, :] = [0.0, 0.0, 1.0]
        end

        idx += 1
    end

    # Generate triangular faces
    faces = Vector{Int}[]
    for i = 1:(n_points - 1)
        for j = 1:(n_points - 1)
            # Two triangles per quad
            v1 = (i - 1) * n_points + j
            v2 = (i - 1) * n_points + j + 1
            v3 = i * n_points + j
            v4 = i * n_points + j + 1

            push!(faces, [v1, v2, v3])
            push!(faces, [v2, v4, v3])
        end
    end

    # Convert faces to matrix
    n_faces = length(faces)
    faces_matrix = Matrix{Int}(undef, n_faces, 3)
    for (i, face) in enumerate(faces)
        faces_matrix[i, :] = face
    end

    return vertices, faces_matrix, vertex_normals
end

"""
    _interpolate_data_on_3d_head(
        channel_data::Vector{Float64},
        layout::Layout,
        vertices::Matrix{Float64},
    )

Interpolate EEG channel data onto 3D head mesh vertices using spherical spline.

# Arguments
- `channel_data::Vector{Float64}`: EEG values at electrode positions
- `layout::Layout`: Layout containing electrode 3D coordinates
- `vertices::Matrix{Float64}`: [n_vertices × 3] head mesh vertex coordinates

# Returns
- `interpolated_values::Vector{Float64}`: Interpolated values at each vertex
"""
function _interpolate_data_on_3d_head(
    channel_data::Vector{Float64},
    layout::Layout,
    vertices::Matrix{Float64},
)
    # Ensure 3D coordinates exist
    _ensure_coordinates_3d!(layout)

    # Extract 3D coordinates efficiently from DataFrame
    n_channels = length(channel_data)
    x3_col = layout.data.x3::Vector{Float64}
    y3_col = layout.data.y3::Vector{Float64}
    z3_col = layout.data.z3::Vector{Float64}

    coords = Matrix{Float64}(undef, n_channels, 3)
    @inbounds for i = 1:n_channels
        coords[i, 1] = x3_col[i]
        coords[i, 2] = y3_col[i]
        coords[i, 3] = z3_col[i]
    end

    # Find the actual radius of the electrode positions
    electrode_radius = 0.0
    @inbounds for i = 1:n_channels
        electrode_radius += sqrt(coords[i, 1]^2 + coords[i, 2]^2 + coords[i, 3]^2)
    end
    electrode_radius /= n_channels

    # Pre-normalize electrode coordinates to unit sphere
    coords_unit = Matrix{Float64}(undef, n_channels, 3)
    @inbounds for i = 1:n_channels
        norm_factor = sqrt(coords[i, 1]^2 + coords[i, 2]^2 + coords[i, 3]^2)
        if norm_factor > 0
            coords_unit[i, 1] = coords[i, 1] / norm_factor
            coords_unit[i, 2] = coords[i, 2] / norm_factor
            coords_unit[i, 3] = coords[i, 3] / norm_factor
        else
            coords_unit[i, 1], coords_unit[i, 2], coords_unit[i, 3] = coords[i, 1], coords[i, 2], coords[i, 3]
        end
    end

    # Step 1: Solve for weights (same as spherical spline)
    cosang = coords_unit * coords_unit'
    G = _calc_g_matrix(cosang)

    G_extended = Matrix{Float64}(undef, n_channels + 1, n_channels + 1)
    @inbounds for i = 1:n_channels
        for j = 1:n_channels
            G_extended[i, j] = G[i, j]
        end
        G_extended[i, i] += 1e-5 # λ regularization
        G_extended[i, n_channels+1] = 1.0
        G_extended[n_channels+1, i] = 1.0
    end
    G_extended[n_channels+1, n_channels+1] = 0.0

    data_vector = Vector{Float64}(undef, n_channels + 1)
    @inbounds for i = 1:n_channels
        data_vector[i] = channel_data[i]
    end
    data_vector[n_channels+1] = 0.0

    weights = G_extended \ data_vector

    # Step 2: Normalize vertex coordinates to unit sphere for interpolation
    n_vertices = size(vertices, 1)
    vertices_unit = Matrix{Float64}(undef, n_vertices, 3)
    @inbounds for i = 1:n_vertices
        norm_factor = sqrt(vertices[i, 1]^2 + vertices[i, 2]^2 + vertices[i, 3]^2)
        if norm_factor > 0
            vertices_unit[i, 1] = vertices[i, 1] / norm_factor
            vertices_unit[i, 2] = vertices[i, 2] / norm_factor
            vertices_unit[i, 3] = vertices[i, 3] / norm_factor
        else
            vertices_unit[i, :] = vertices[i, :]
        end
    end

    # Step 3: Compute G matrix between vertices and channels
    cosang_vertices = vertices_unit * coords_unit'

    # Step 4: Apply g-function
    factors = _get_g_factors(15)
    @inbounds for i in eachindex(cosang_vertices)
        cosang_vertices[i] = _legendre_val(clamp(cosang_vertices[i], -1.0, 1.0), factors)
    end

    # Step 5: Compute interpolated values
    w_channels = @view weights[1:n_channels]
    interpolated_values = cosang_vertices * w_channels
    interpolated_values .+= weights[n_channels+1]

    return interpolated_values
end

"""
    _project_electrodes_to_head_surface(
        layout::Layout,
        vertices::Matrix{Float64},
        faces::Matrix{Int},
        mesh_type::Union{Symbol, Nothing} = nothing,
    )

Project electrode positions onto the head mesh surface.

For OBJ meshes, uses ray-triangle intersection to find where each electrode's
direction vector intersects the mesh surface. For parametric meshes, uses
the mesh vertices directly.

# Arguments
- `layout::Layout`: Layout containing electrode 3D coordinates
- `vertices::Matrix{Float64}`: [n_vertices × 3] head mesh vertex coordinates
- `faces::Matrix{Int}`: [n_faces × 3] face indices (1-indexed)
- `mesh_type::Union{Symbol, Nothing}`: `:obj` for OBJ mesh, `nothing` for parametric

# Returns
- `electrode_positions::Matrix{Float64}`: [n_electrodes × 3] projected electrode positions
"""
function _project_electrodes_to_head_surface(
    layout::Layout,
    vertices::Matrix{Float64},
    faces::Matrix{Int},
    mesh_type::Union{Symbol, Nothing} = nothing,
)
    _ensure_coordinates_3d!(layout)
    
    n_electrodes = nrow(layout.data)
    electrode_positions = Matrix{Float64}(undef, n_electrodes, 3)
    
    if mesh_type == :obj
        # For OBJ mesh: find closest vertex in each electrode's direction
        # Find Cz index and its direction to determine top of head
        cz_idx = findfirst(==(Symbol("Cz")), layout.data.label)
        
        # Get Cz direction vector (should point to top of head)
        if cz_idx !== nothing
            cz_x = layout.data.x3[cz_idx]
            cz_y = layout.data.y3[cz_idx]
            cz_z = layout.data.z3[cz_idx]
            cz_norm = sqrt(cz_x^2 + cz_y^2 + cz_z^2)
            if cz_norm > 0
                cz_dir = [cz_x / cz_norm, cz_y / cz_norm, cz_z / cz_norm]
            else
                cz_dir = [0.0, 0.0, 1.0]  # Default to +z if Cz has no direction
            end
        else
            cz_dir = [0.0, 0.0, 1.0]  # Default to +z if Cz not found
        end
        
        # Find vertex that best matches Cz direction (top of head)
        best_dot = -Inf
        top_vertex_idx = 1
        n_vertices = size(vertices, 1)
        @inbounds for v = 1:n_vertices
            vertex = @view vertices[v, :]
            vertex_norm = sqrt(vertex[1]^2 + vertex[2]^2 + vertex[3]^2)
            if vertex_norm > 0
                vertex_dir = [vertex[1] / vertex_norm, vertex[2] / vertex_norm, vertex[3] / vertex_norm]
                dot_product = cz_dir[1] * vertex_dir[1] + cz_dir[2] * vertex_dir[2] + cz_dir[3] * vertex_dir[3]
                if dot_product > best_dot
                    best_dot = dot_product
                    top_vertex_idx = v
                end
            end
        end
        top_vertex = @view vertices[top_vertex_idx, :]
        top_radius = sqrt(top_vertex[1]^2 + top_vertex[2]^2 + top_vertex[3]^2)
        
        @inbounds for i = 1:n_electrodes
            # Get electrode direction (normalized)
            x_e = layout.data.x3[i]
            y_e = layout.data.y3[i]
            z_e = layout.data.z3[i]
            norm_e = sqrt(x_e^2 + y_e^2 + z_e^2)
            
            if norm_e > 0
                dir = [x_e / norm_e, y_e / norm_e, z_e / norm_e]
            else
                dir = [x_e, y_e, z_e]
            end
            
            # Find closest vertex in this direction
            best_dot = -Inf
            best_vertex = [0.0, 0.0, 0.0]
            @inbounds for v = 1:n_vertices
                vertex = @view vertices[v, :]
                vertex_norm = sqrt(vertex[1]^2 + vertex[2]^2 + vertex[3]^2)
                if vertex_norm > 0
                    vertex_dir = [vertex[1] / vertex_norm, vertex[2] / vertex_norm, vertex[3] / vertex_norm]
                    dot_product = dir[1] * vertex_dir[1] + dir[2] * vertex_dir[2] + dir[3] * vertex_dir[3]
                    if dot_product > best_dot
                        best_dot = dot_product
                        best_vertex = vertex
                    end
                end
            end
            
            # Scale to match head size (use top radius as reference)
            vertex_radius = sqrt(best_vertex[1]^2 + best_vertex[2]^2 + best_vertex[3]^2)
            if vertex_radius > 0 && top_radius > 0
                scale = top_radius / vertex_radius
                # For Cz, ensure it's exactly at the top
                if cz_idx !== nothing && i == cz_idx
                    electrode_positions[i, :] = top_vertex
                else
                    electrode_positions[i, :] = best_vertex .* scale
                end
            else
                electrode_positions[i, :] = best_vertex
            end
        end
    else
        # For parametric mesh: electrodes are already on the sphere, use as-is
        @inbounds for i = 1:n_electrodes
            electrode_positions[i, 1] = layout.data.x3[i]
            electrode_positions[i, 2] = layout.data.y3[i]
            electrode_positions[i, 3] = layout.data.z3[i]
        end
    end
    
    return electrode_positions
end

"""
    plot_topography_3d(
        dat::SingleDataFrameEeg;
        channel_selection::Function = channels(),
        sample_selection::Function = samples(),
        display_plot::Bool = true,
        kwargs...,
    )

Create a 3D topographic plot on a parametric head model or custom OBJ mesh.

# Arguments
- `dat::SingleDataFrameEeg`: EEG data (ContinuousData or ErpData)
- `channel_selection::Function`: Channel selection predicate (default: all channels)
- `sample_selection::Function`: Sample selection predicate (default: all samples)
- `display_plot::Bool`: Whether to display the plot (default: true)
- `sphere_resolution::Int`: Resolution of parametric head mesh (default: 50). Ignored if `head_mesh_file` is provided.
- `head_mesh_file::Union{String, Nothing}`: Path to OBJ file containing 3D head mesh (default: nothing). If provided, uses this mesh instead of generating a parametric one.
- `colormap`: Colormap for data visualization (default: :RdBu)
- `colorrange`: Color range for data (default: automatic)
- `show_electrodes::Bool`: Whether to show electrode positions (default: true)
- `kwargs...`: Additional keyword arguments

# Returns
- `fig::Figure`: The figure object
- `ax::Axis3`: The 3D axis object

# Examples
```julia
# Plot 3D topography with parametric head model
fig, ax = plot_topography_3d(erp_data; sample_selection=samples((0.1, 0.2)))

# Plot 3D topography with custom OBJ head model
fig, ax = plot_topography_3d(erp_data; head_mesh_file="path/to/head_model.obj")
```
"""
function _generate_sphere_mesh_from_electrodes(layout::Layout, resolution::Int = 50)
    _ensure_coordinates_3d!(layout)
    
    # Get electrode coordinates (same as plot_layout_3d uses)
    x3_electrodes = layout.data.x3::Vector{Float64}
    y3_electrodes = layout.data.y3::Vector{Float64}
    z3_electrodes = layout.data.z3::Vector{Float64}
    
    # Convert electrode positions to spherical coordinates to understand their distribution
    n_channels = length(x3_electrodes)
    inc_electrodes = Vector{Float64}(undef, n_channels)
    azi_electrodes = Vector{Float64}(undef, n_channels)
    radius_electrodes = Vector{Float64}(undef, n_channels)
    
    @inbounds for i = 1:n_channels
        x, y, z = x3_electrodes[i], y3_electrodes[i], z3_electrodes[i]
        r = sqrt(x^2 + y^2 + z^2)
        radius_electrodes[i] = r
        if r > 0
            inc_electrodes[i] = acos(clamp(z / r, -1.0, 1.0))
            azi_electrodes[i] = atan(y, x)
            if azi_electrodes[i] < 0
                azi_electrodes[i] += 2π
            end
        else
            inc_electrodes[i] = 0.0
            azi_electrodes[i] = 0.0
        end
    end
    
    # Generate dense grid in spherical coordinates
    inc_range = range(0, π, length = resolution)
    azi_range = range(0, 2π, length = resolution)
    
    n_vertices = resolution * resolution
    vertices = Matrix{Float64}(undef, n_vertices, 3)
    
    idx = 1
    for (i, inc) in enumerate(inc_range), (j, azi) in enumerate(azi_range)
        # For each grid point, find the radius by interpolating from nearby electrodes
        # Use inverse distance weighting to get radius at this (inc, azi)
        total_weight = 0.0
        weighted_radius = 0.0
        
        @inbounds for k = 1:n_channels
            # Angular distance between grid point and electrode
            inc_diff = inc - inc_electrodes[k]
            azi_diff = azi - azi_electrodes[k]
            # Handle azimuth wrap-around
            if azi_diff > π
                azi_diff -= 2π
            elseif azi_diff < -π
                azi_diff += 2π
            end
            
            # Angular distance (simplified - uses great circle distance approximation)
            angular_dist = sqrt(inc_diff^2 + (azi_diff * sin(inc))^2)
            
            # Inverse distance weighting (with small epsilon to avoid division by zero)
            if angular_dist < 0.01
                # Very close to an electrode, use its radius directly
                weighted_radius = radius_electrodes[k]
                total_weight = 1.0
                break
            else
                weight = 1.0 / (angular_dist^2 + 0.01)  # Small epsilon for stability
                weighted_radius += weight * radius_electrodes[k]
                total_weight += weight
            end
        end
        
        radius = total_weight > 0 ? weighted_radius / total_weight : 1.0
        
        # Convert back to Cartesian using the interpolated radius
        # Use the same coordinate system as the electrodes (no rotation needed)
        vertices[idx, 1] = radius * sin(inc) * cos(azi)
        vertices[idx, 2] = radius * sin(inc) * sin(azi)
        vertices[idx, 3] = radius * cos(inc)
        idx += 1
    end
    
    # Generate triangular faces
    faces = Vector{Int}[]
    for i = 1:(resolution - 1)
        for j = 1:(resolution - 1)
            v1 = (i - 1) * resolution + j
            v2 = (i - 1) * resolution + j + 1
            v3 = i * resolution + j
            v4 = i * resolution + j + 1
            
            push!(faces, [v1, v2, v3])
            push!(faces, [v2, v4, v3])
        end
    end
    
    # Convert faces to matrix
    n_faces = length(faces)
    faces_matrix = Matrix{Int}(undef, n_faces, 3)
    for (i, face) in enumerate(faces)
        faces_matrix[i, :] = face
    end
    
    return vertices, faces_matrix
end

function plot_topography_3d(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot::Bool = true,
    sphere_resolution::Int = 50,
    head_mesh_file::Union{String, Nothing} = nothing,
    colormap = :RdBu,
    colorrange = nothing,
    show_electrodes::Bool = true,
    kwargs...,
)
    # Subset data
    dat_subset = subset(dat; channel_selection = channel_selection, sample_selection = sample_selection)

    # Ensure 3D coordinates exist
    _ensure_coordinates_3d!(dat_subset.layout)

    # Compute average channel data
    channel_data = mean.(eachcol(dat_subset.data[!, dat_subset.layout.data.label]))

    # Load head mesh from OBJ file or generate parametric mesh
    if !isnothing(head_mesh_file)
        vertices, faces = _read_obj_file(head_mesh_file)
    else
        vertices, faces = _generate_sphere_mesh_from_electrodes(dat_subset.layout, sphere_resolution)
    end

    # Interpolate data onto sphere mesh using existing function
    interpolated_values = _interpolate_data_on_3d_head(channel_data, dat_subset.layout, vertices)

    # Determine color range
    if isnothing(colorrange)
        data_min, data_max = extrema(interpolated_values[.!isnan.(interpolated_values)])
        max_abs = max(abs(data_min), abs(data_max))
        colorrange = (-max_abs, max_abs)
    end

    # Create figure and axis
    set_window_title(_generate_window_title(dat))
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)

    # For parametric mesh, reshape to grid format for surface plotting
    # For OBJ mesh, use mesh plotting instead
    if isnothing(head_mesh_file)
        # Parametric mesh: reshape vertices and values to grid for surface plotting
        # The vertices are generated in order: for each inc (i), for each azi (j)
        # This means: vertices[1:resolution] = all azi at inc=0 (top), 
        #            vertices[resolution+1:2*resolution] = all azi at inc=π/resolution, etc.
        # When reshaping, this creates a matrix where:
        # - Rows correspond to inclination (inc) - from top (inc=0) to bottom (inc=π)
        # - Columns correspond to azimuth (azi) - from 0 to 2π
        n_points = sphere_resolution
        
        # Rotate vertices to match 2D topography orientation
        # In 2D, frontal electrodes are at the top. Rotate 45 degrees around y-axis:
        # x' = x*cos(45°) + z*sin(45°), y' = y, z' = -x*sin(45°) + z*cos(45°)
        # cos(45°) = sin(45°) = 1/√2
        cos_45 = 1.0 / sqrt(2.0)
        sin_45 = 1.0 / sqrt(2.0)
        vertices_rotated = similar(vertices)
        vertices_rotated[:, 1] = vertices[:, 1] * cos_45 + vertices[:, 3] * sin_45  # x' = x*cos + z*sin
        vertices_rotated[:, 2] = vertices[:, 2]                                     # y' = y
        vertices_rotated[:, 3] = -vertices[:, 1] * sin_45 + vertices[:, 3] * cos_45 # z' = -x*sin + z*cos
        
        x_grid = reshape(vertices_rotated[:, 1], n_points, n_points)
        y_grid = reshape(vertices_rotated[:, 2], n_points, n_points)
        z_grid = reshape(vertices_rotated[:, 3], n_points, n_points)
        values_grid = reshape(interpolated_values, n_points, n_points)
        
        # Plot surface with interpolated colors
        surface_plot = surface!(
            ax,
            x_grid,
            y_grid,
            z_grid;
            color = values_grid,
            colormap = colormap,
            colorrange = colorrange,
            shading = true,
            alpha = 1.0,
        )
    else
        # OBJ mesh: use mesh plotting with vertex colors
        # Convert faces to Makie-compatible format (0-indexed, flat vector)
        # OBJ files use 1-indexing, so subtract 1 to get 0-indexed
        # Reverse winding order to fix backface culling (OBJ typically uses CCW, Makie may expect CW)
        n_faces = size(faces, 1)
        faces_makie = Vector{UInt32}(undef, n_faces * 3)
        @inbounds for i = 1:n_faces
            idx = (i - 1) * 3 + 1
            # OBJ indices are 1-indexed, convert to 0-indexed for Makie
            # Ensure indices are valid (>= 1) before converting
            v1 = faces[i, 1]
            v2 = faces[i, 2]
            v3 = faces[i, 3]
            if v1 < 1 || v2 < 1 || v3 < 1
                throw(ArgumentError("Invalid face indices in OBJ file: face $i has indices ($v1, $v2, $v3). OBJ files must use 1-indexed vertex indices."))
            end
            # Keep original winding order (OBJ files typically use counter-clockwise)
            faces_makie[idx] = UInt32(v1 - 1)
            faces_makie[idx + 1] = UInt32(v2 - 1)
            faces_makie[idx + 2] = UInt32(v3 - 1)
        end
        
        # Plot mesh with interpolated vertex colors
        # Filter out NaN values and replace with a default color
        valid_mask = .!isnan.(interpolated_values)
        if !all(valid_mask)
            # Replace NaN values with the mean of valid values for solid rendering
            valid_mean = mean(interpolated_values[valid_mask])
            interpolated_values_clean = copy(interpolated_values)
            interpolated_values_clean[.!valid_mask] .= valid_mean
        else
            interpolated_values_clean = interpolated_values
        end
        
        # Plot mesh with interpolated vertex colors
        mesh_plot = mesh!(
            ax,
            vertices,
            faces_makie;
            color = interpolated_values_clean,
            colormap = colormap,
            colorrange = colorrange,
            shading = true,
            alpha = 1.0,  # Fully opaque
            transparency = false,
            overdraw = false,  # Don't draw over other objects
        )
    end


    # Optionally show electrode positions - plot AFTER surface so they're on top
    if show_electrodes
        # Project electrodes onto head mesh surface
        electrode_positions = _project_electrodes_to_head_surface(
            dat_subset.layout,
            vertices,
            faces,
            isnothing(head_mesh_file) ? nothing : :obj  # Different projection for OBJ vs parametric
        )
        
        # Offset electrodes slightly outward from surface to ensure they're always visible
        scale_factor = 1.02
        x_electrodes = electrode_positions[:, 1] .* scale_factor
        y_electrodes = electrode_positions[:, 2] .* scale_factor
        z_electrodes = electrode_positions[:, 3] .* scale_factor
        
        scatter!(
            ax,
            x_electrodes,
            y_electrodes,
            z_electrodes;
            color = :black,
            markersize = 10,
            marker = :circle,
            strokewidth = 1,
            strokecolor = :black,
            overdraw = false,  # Disable depth testing so electrodes always appear on top
        )
    end

    # Add colorbar
    Colorbar(fig[1, 2], colormap = colormap, limits = colorrange, label = "Amplitude")

    # Set title
    time_min, time_max = extrema(dat_subset.data.time)
    ax.title = @sprintf("%.3f to %.3f s", time_min, time_max)

    hidedecorations!(ax)
    hidespines!(ax)

    display_plot && display_figure(fig)
    set_window_title("Makie")

    return fig, ax
end
