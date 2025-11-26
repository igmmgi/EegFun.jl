# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_TOPOGRAPHY_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    :figure_title => ("EEG Topography Plot", "Title for the plot window"),
    :interactive => (true, "Whether to enable interactive features"),

    # Topography-specific parameters
    :method => (:multiquadratic, "Interpolation method: :multiquadratic or :spherical_spline"),
    :gridscale => (200, "Grid resolution for interpolation"),
    :colormap => (:jet, "Colormap for the topography"),
    :ylim => (nothing, "Y-axis limits (nothing for auto)"),

    # Title parameters
    :title => ("", "Plot title"),
    :title_fontsize => (16, "Font size for the title"),
    :show_title => (true, "Whether to show the title"),

    # Colorbar parameters - get all Colorbar attributes with their actual defaults
    # This allows users to control any Colorbar parameter
    [
        Symbol("colorbar_$(attr)") => (get(COLORBAR_DEFAULTS, attr, nothing), "Colorbar $(attr) parameter") for
        attr in propertynames(Colorbar)
    ]...,

    # Override specific colorbar parameters with custom defaults
    :colorbar_plot => (true, "Whether to display the colorbar"),
    :colorbar_position => ((1, 2), "Position of the colorbar as (row, col) tuple"),
    :colorbar_label => ("μV", "Label for the colorbar"),

    # Head shape parameters (reusing layout kwargs)
    :head_color => (:black, "Color of the head shape outline."),
    :head_linewidth => (2, "Line width of the head shape outline."),
    :head_radius => (1.0, "Radius of the head shape in mm."),

    # Electrode point parameters
    :point_plot => (true, "Whether to plot electrode points."),
    :point_marker => (:circle, "Marker style for electrode points."),
    :point_markersize => (12, "Size of electrode point markers."),
    :point_color => (:black, "Color of electrode points."),

    # Electrode label parameters
    :label_plot => (true, "Whether to plot electrode labels."),
    :label_fontsize => (20, "Font size for electrode labels."),
    :label_color => (:black, "Color of electrode labels."),
    :label_xoffset => (0, "X-axis offset for electrode labels."),
    :label_yoffset => (0, "Y-axis offset for electrode labels."),
)


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
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    channel_data = mean.(eachcol(dat[!, layout.data.label]))
    if method == :multiquadratic
        data = _data_interpolation_topo_multiquadratic(channel_data, layout, gridscale)
    elseif method == :spherical_spline
        data = _data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
    end

    ylim = pop!(plot_kwargs, :ylim)
    if isnothing(ylim)
        # Make ylim symmetric around 0 for balanced topographic visualization
        data_min, data_max = extrema(data[.!isnan.(data)])
        max_abs = max(abs(data_min), abs(data_max))
        ylim = (-max_abs, max_abs)
    end

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

    interactive && _setup_topo_interactivity!(fig, ax, dat)
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
    
    # Calculate optimal grid dimensions
    rows, cols = best_rect(n_datasets)
    
    # Create single figure with subplots
    set_window_title(_generate_window_title(dat))
    fig = Figure()
    axes = Axis[]
    
    # Plot each dataset in its own subplot
    for (idx, dataset) in enumerate(dat)

        base_row = div(idx - 1, cols) + 1
        base_col = mod1(idx, cols)
        
        # Each subplot gets 2 columns: plot in first column, colorbar in second
        plot_row = base_row
        plot_col = (base_col - 1) * 2 + 1
        colorbar_row = plot_row
        colorbar_col = plot_col + 1
        
        ax = Axis(fig[plot_row, plot_col])
        push!(axes, ax)
        
        # Set subplot title to condition name if available
        if hasproperty(dataset, :condition_name) && dataset.condition_name !== ""
            ax.title = dataset.condition_name
        end
        
        # Calculate colorbar position relative to this subplot
        kwargs_dict = Dict{Symbol, Any}(kwargs)
        kwargs_dict[:colorbar_position] = (colorbar_row, colorbar_col)
        
        # Plot the topography in this subplot
        plot_topography!(
            fig,
            ax,
            dataset;
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            kwargs_dict...,
        )
    end
    
    # Set column sizes: make colorbar columns narrower than plot columns
    # This ensures colorbars don't take up too much space while keeping plots visible
    total_cols = cols * 2
    for col in 1:total_cols
        if col % 2 == 1
            colsize!(fig.layout, col, Auto())
        else
            colsize!(fig.layout, col, Fixed(50))
        end
    end
    
    # Set up interactivity for all axes if requested
    if interactive
        shared_selection_state = TopoSelectionState(axes)
        _setup_shared_topo_interactivity!(fig, axes, dat, shared_selection_state)
    end
    
    display_plot && display_figure(fig)
    
    set_window_title("Makie")
    return fig, axes
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
    interactive = true,
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

    interactive && _setup_topo_interactivity!(fig, ax, dat)
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
    valid_mask = (r_2d .<= 1.0) .& (r_2d .<= 1.0)
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
    _setup_topo_interactivity!(fig::Figure, ax::Axis)

Set up keyboard interactivity for topographic plots.
"""
function _setup_topo_interactivity!(fig::Figure, ax::Axis, original_data = nothing)

    deregister_interaction!(ax, :rectanglezoom)

    # Handle keyboard events at the figure level
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.i
                show_plot_help(:topography)
            elseif event.key == Keyboard.up
                _topo_scale_up!(ax)
            elseif event.key == Keyboard.down
                _topo_scale_down!(ax)
            end
        end
    end

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
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.i
                show_plot_help(:topography)
            elseif event.key == Keyboard.up
                _topo_scale_up!.(axes)
            elseif event.key == Keyboard.down
                _topo_scale_down!.(axes)
            end
        end
    end
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
    _setup_topo_selection!(fig::Figure, ax::Axis, original_data)

Set up simple region selection for topographic plots.
"""
function _setup_topo_selection!(fig::Figure, ax::Axis, original_data)
    selection_state = TopoSelectionState([ax])

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
        # Check if mouse is over this specific axis
        mouse_pos = events(fig).mouseposition[]
        if !_is_mouse_in_axis(ax, mouse_pos)
            return
        end
        
        if event.button == Mouse.left
            if event.action == Mouse.press
                # Check if Shift is held down
                if shift_pressed[]
                    _start_topo_selection!(ax, selection_state)
                end
            elseif event.action == Mouse.release
                if selection_state.active[]
                    _finish_topo_selection!(ax, selection_state, original_data)
                end
            end
        end
        if event.button == Mouse.left && event.action == Mouse.press && !shift_pressed[]
            _clear_all_topo_selections!(selection_state)
        end
        if event.button == Mouse.right && event.action == Mouse.press
            selected_channels = selection_state.selected_channels
            if !isempty(selected_channels) && original_data isa ErpData
                _show_topo_context_menu!(original_data, selected_channels)
            elseif !isempty(selected_channels)
                @warn "Cannot plot ERP: data is not ErpData (got $(typeof(original_data)))"
            else
                @info "No channels selected. Select channels with Shift+Left Click+Drag, then right-click to plot ERP"
            end
        end
    end

    # Handle mouse movement for selection dragging
    on(events(fig).mouseposition) do pos
        # Only update if mouse is over this axis and selection is active
        if selection_state.active[] && _is_mouse_in_axis(ax, pos)
            _update_topo_selection!(ax, selection_state)
        end
    end
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
        active_ax = nothing
        active_ax_idx = nothing
        active_dataset = nothing
        
        for (idx, check_ax) in enumerate(shared_selection_state.axes)
            if _is_mouse_in_axis(check_ax, mouse_pos)
                active_ax = check_ax
                active_ax_idx = idx
                if idx <= length(datasets)
                    active_dataset = datasets[idx]
                end
                break
            end
        end
        
        # Only process if mouse is over one of the shared axes
        if active_ax === nothing
            return
        end
        
        if event.button == Mouse.left
            if event.action == Mouse.press
                # Check if Shift is held down
                if shift_pressed[]
                    _start_topo_selection!(active_ax, shared_selection_state)
                end
            elseif event.action == Mouse.release
                if shared_selection_state.active[]
                    _finish_topo_selection!(active_ax, shared_selection_state, active_dataset)
                end
            end
        end
        if event.button == Mouse.left && event.action == Mouse.press && !shift_pressed[]
            _clear_all_topo_selections!(shared_selection_state)
        end
        if event.button == Mouse.right && event.action == Mouse.press
            selected_channels = shared_selection_state.selected_channels
            # Check if all datasets are ErpData
            all_erp_data = all(d isa ErpData for d in datasets)
            if !isempty(selected_channels) && all_erp_data
                _show_topo_context_menu!(datasets, selected_channels)
            elseif !isempty(selected_channels)
                @warn "Cannot plot ERP: data is not ErpData"
            else
                @info "No channels selected. Select channels with Shift+Left Click+Drag, then right-click to plot ERP"
            end
        end
    end

    # Handle mouse movement for selection dragging
    on(events(fig).mouseposition) do pos
        # Only update if selection is active and mouse is over any shared axis
        if shared_selection_state.active[]
            for check_ax in shared_selection_state.axes
                if _is_mouse_in_axis(check_ax, pos)
                    _update_topo_selection!(check_ax, shared_selection_state)
                    break
                end
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

    # Find the index of this axis
    ax_idx = findfirst(isequal(ax), selection_state.axes)
    if ax_idx === nothing
        ax_idx = 1
    end
    
    # Create temporary rectangles for all axes
    initial_points =
        [Point2f(mouse_x, mouse_y), Point2f(mouse_x, mouse_y), Point2f(mouse_x, mouse_y), Point2f(mouse_x, mouse_y)]
    
    # Create temp rectangle on the active axis
    selection_state.temp_rectangles[ax_idx] = poly!(
        ax,
        initial_points,
        color = (:blue, 0.3),    
        strokecolor = :black,    
        strokewidth = 1,         
        visible = true,
        overdraw = true,          
    )
    
    # Create temp rectangles on all other axes with same bounds
    for (other_idx, other_ax) in enumerate(selection_state.axes)
        if other_idx != ax_idx
            selection_state.temp_rectangles[other_idx] = poly!(
                other_ax,
                initial_points,
                color = (:blue, 0.3),
                strokecolor = :black,
                strokewidth = 1,
                visible = true,
                overdraw = true,
            )
        end
    end

    _update_topo_selection!(ax, selection_state)
end

"""
    _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, original_data)

Finish spatial region selection in topographic plot.
"""
function _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, original_data = nothing)
    selection_state.active[] = false

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
    bounds = selection_state.bounds[]
    start_x, start_y = bounds[1], bounds[2]
    end_x, end_y = bounds[3], bounds[4]

    rect_points = Point2f[
        Point2f(Float64(start_x), Float64(start_y)),
        Point2f(Float64(end_x), Float64(start_y)),
        Point2f(Float64(end_x), Float64(end_y)),
        Point2f(Float64(start_x), Float64(end_y)),
    ]

    # Find the index of this axis
    ax_idx = findfirst(isequal(ax), selection_state.axes)
    if ax_idx === nothing
        ax_idx = 1
    end

    # Create permanent rectangles on all axes
    for (other_idx, other_ax) in enumerate(selection_state.axes)
        permanent_rectangle = poly!(
            other_ax,
            rect_points,
            color = (:blue, 0.3),    
            strokecolor = :black,     
            strokewidth = 1,          
            visible = true,
            overdraw = true,           
        )
        push!(selection_state.rectangles[other_idx], permanent_rectangle)
    end

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
        selection_state.bounds[] = (start_x, start_y, axis_pos[1], axis_pos[2])

        # Update the temporary rectangles on all axes during dragging
        bounds = selection_state.bounds[]
        start_x, start_y = bounds[1], bounds[2]
        end_x, end_y = bounds[3], bounds[4]

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
    selected_electrodes = Symbol[]
    for row in eachrow(original_data.layout.data)
        # Check if this electrode is inside the selection rectangle
        if x_min <= row.x2 <= x_max && y_min <= row.y2 <= y_max
            push!(selected_electrodes, Symbol(row.label))
        end
    end
    return selected_electrodes
end

"""
    _show_topo_context_menu!(datasets::Union{ErpData, Vector{ErpData}}, selected_channels::Vector{Symbol})

Show a context menu for plotting selected channels from topography plot.
Supports both single dataset and multiple datasets (conditions).
"""
function _show_topo_context_menu!(datasets::Union{ErpData, Vector{ErpData}}, selected_channels::Vector{Symbol})

    datasets_vec = datasets isa Vector ? datasets : [datasets]
    has_multiple_conditions = length(datasets_vec) > 1
    
    plot_types = String[]
    if has_multiple_conditions
        push!(plot_types, "Separate electrodes, separate conditions")
        push!(plot_types, "Average electrodes, separate conditions")
        push!(plot_types, "Separate electrodes, average conditions")
        push!(plot_types, "Average electrodes, average conditions")
    else
        push!(plot_types, "Plot Individual Channels")
        push!(plot_types, "Plot Averaged Channels")
    end
    
    menu_fig = Figure(size = (400, 200))
    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]
    
    # Define plot configurations: (average_conditions, average_channels, info_msg)
    plot_configs = has_multiple_conditions ? [
        (false, false, "separate electrodes, separate conditions"),
        (false, true, "average electrodes, separate conditions"),
        (true, false, "separate electrodes, average conditions"),
        (true, true, "average electrodes, average conditions"),
    ] : [
        (false, false, "individual channels"),
        (false, true, "averaged channels"),
    ]
    
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
    
    new_screen = getfield(Main, :GLMakie).Screen()
    display(new_screen, menu_fig)
end
