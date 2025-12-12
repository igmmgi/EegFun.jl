# Note: PLOT_TOPOGRAPHY_KWARGS is defined in plot_ica.jl

##########################################
# 2D topographic plot
##########################################
function _plot_topography!(fig::Figure, ax::Axis, dat::DataFrame, layout::Layout; kwargs...)

    # ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

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
    if method == :multiquadratic
        data = _data_interpolation_topo_multiquadratic(channel_data, layout, gridscale)
    elseif method == :spherical_spline
        data = _data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
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
        nan_color = :transparent
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
    kwargs_dict = Dict{Symbol, Any}(kwargs)
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
        for col in 1:total_cols
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
    _data_interpolation_topo(dat::Vector{Float64}, points::Matrix{Float64}, grid_scale::Int)

Interpolate EEG data using scattered interpolation 

# Arguments
- `dat::Vector{<:AbstractFloat}`: EEG values at electrode positions
- `points::Matrix{<:AbstractFloat}`: 2×N matrix of electrode coordinates
- `grid_scale::Int`: Size of the output grid
"""
function _data_interpolation_topo_multiquadratic(dat::Vector{<:AbstractFloat}, layout::Layout, grid_scale::Int)

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

    itp = ScatteredInterpolation.interpolate(Multiquadratic(), points, dat)
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

    # Extract 3D coordinates as they are (no normalization needed)
    n_channels = length(dat)
    coords = zeros(Float64, n_channels, 3)
    for i = 1:n_channels
        coords[i, :] = [layout.data.x3[i], layout.data.y3[i], layout.data.z3[i]]
    end

    # Find the actual radius of the electrode positions
    electrode_radius = mean([sqrt(sum(coords[i, :] .^ 2)) for i = 1:n_channels])

    # Create a 2D grid for plotting using normalized coordinates
    x_range = range(-1.0, 1.0, length = grid_scale)
    y_range = range(-1.0, 1.0, length = grid_scale)

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
    # Since coordinates are normalized to [-1, 1], the head radius is 1.0
    valid_mask = r_2d .<= 1.0
    valid_indices = findall(valid_mask)

    if !isempty(valid_indices)
        # Extract valid grid points
        valid_x = grid_x[valid_indices]
        valid_y = grid_y[valid_indices]
        valid_r = r_2d[valid_indices]

        # Pre-compute stereographic projection for all valid points at once
        # Since coordinates are normalized, use 1.0 as the head radius
        r_norm = valid_r ./ 1.0

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

# =============================================================================
# INTERACTIVITY FUNCTIONS
# =============================================================================

"""
    _setup_topo_keyboard_handlers!(fig::Figure, axes::Union{Axis, Vector{Axis}})

Set up keyboard event handlers for topographic plots.
Handles both single axis and multiple axes.
"""
function _setup_topo_keyboard_handlers!(fig::Figure, axes::Union{Axis, Vector{Axis}})
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
            [Makie.Poly[] for _ in 1:n_axes],  # Empty vectors for rectangles per axis
            Observable{Tuple{Float64,Float64,Float64,Float64}}[],  # Empty vector for bounds
            [nothing for _ in 1:n_axes],  # No temporary rectangles initially
            Symbol[],  # Empty vector for selected channels
            axes,  # Store all axes
        )
    end
end

"""
    _setup_shared_topo_interactivity!(fig::Figure, axes::Vector{Axis}, datasets::Vector, shared_selection_state)

Set up shared interactivity for multiple topographic plots.
"""
function _setup_shared_topo_interactivity!(fig::Figure, axes::Vector{Axis}, datasets::Vector, shared_selection_state::TopoSelectionState)
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
function _create_rectangles_on_all_axes!(selection_state::TopoSelectionState, rect_points::Vector{Point2f}, is_temporary::Bool)
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
        active_ax, active_dataset = _find_active_axis_with_dataset(
            shared_selection_state.axes, mouse_pos, datasets
        )
        
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
    initial_points = [Point2f(mouse_x, mouse_y) for _ in 1:4]
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
    selected_rows = Base.filter(row -> x_min <= row.x2 <= x_max && y_min <= row.y2 <= y_max, eachrow(original_data.layout.data))
    return [Symbol(row.label) for row in selected_rows]
end

"""
    _show_topo_context_menu!(datasets::Union{ErpData, Vector{ErpData}}, selected_channels::Vector{Symbol})

Show a context menu for plotting selected channels from topography plot.
Supports both single dataset and multiple datasets (conditions).
"""
function _show_topo_context_menu!(datasets::Union{ErpData, Vector{ErpData}}, selected_channels::Vector{Symbol})

    datasets_vec = datasets isa Vector ? datasets : [datasets]
    has_multiple_conditions = length(datasets_vec) > 1
    has_multiple_channels = length(selected_channels) > 1
    
    plot_types = String[]
    plot_configs = Tuple{Bool, Bool, String}[]
    
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
        plot_configs = [
            (false, false, "separate conditions"),
            (true, false, "average conditions"),
        ]
    elseif has_multiple_channels
        # Single condition but multiple channels: only channel averaging options
        push!(plot_types, "Plot Individual Channels")
        push!(plot_types, "Plot Averaged Channels")
        plot_configs = [
            (false, false, "individual channels"),
            (false, true, "averaged channels"),
        ]
    else
        # Single condition and single channel: just plot it
        push!(plot_types, "Plot Channel")
        plot_configs = [
            (false, false, "single channel"),
        ]
    end
    
    menu_fig = Figure(size = (400, 200))
    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]
    
    for (idx, btn) in enumerate(menu_buttons)
        on(btn.clicks) do n
            avg_conditions, avg_channels, msg = plot_configs[idx]
            
            # Prepare data: average conditions if needed
            data_to_plot = avg_conditions ? _average_conditions(datasets_vec) : datasets_vec
            @info "Plotting ERP: $msg: $selected_channels"
            plot_erp(
                data_to_plot;
                channel_selection = channels(selected_channels),
                average_channels = avg_channels,
            )
        end
    end
    
    new_screen = GLMakie.Screen()
    display(new_screen, menu_fig)
end
