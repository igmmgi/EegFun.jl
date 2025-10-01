# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_TOPOGRAPHY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    :interactive => (true, "Whether to enable interactive features"),

    # Topography-specific parameters
    :method => (:multiquadratic, "Interpolation method: :multiquadratic or :spherical_spline"),
    :gridscale => (200, "Grid resolution for interpolation"),
    :colormap => (:jet, "Colormap for the topography"),
    :ylim => (nothing, "Y-axis limits (nothing for auto)"),
    :colorrange => (nothing, "Color range for the topography. If nothing, automatically determined"),
    :nan_color => (:transparent, "Color for NaN values"),

    # Title parameters
    :title => ("", "Plot title"),
    :title_fontsize => (16, "Font size for the title"),
    :show_title => (true, "Whether to show the title"),

    # Axis labels and styling
    :xlabel => ("", "Label for x-axis"),
    :ylabel => ("", "Label for y-axis"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),

    # Colorbar parameters
    :plot_colorbar => (true, "Whether to display the colorbar"),
    :colorbar_width => (30, "Width of the colorbar"),
    :colorbar_label => ("μV", "Label for the colorbar"),
    :colorbar_fontsize => (12, "Font size for colorbar labels"),

    # Grid parameters
    :grid_visible => (false, "Whether to show grid"),
    :grid_alpha => (0.3, "Transparency of grid"),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),

    # Head shape parameters (reusing layout kwargs)
    :head_color => (:black, "Color of the head shape outline."),
    :head_linewidth => (2, "Line width of the head shape outline."),
    :head_radius => (1.0, "Radius of the head shape in mm."),
    :head_ear_ratio => (1/7, "Ratio of ear size to head radius."),
    :head_nose_scale => (4.0, "Scale factor for nose size."),

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

    # Extract specific values
    method = plot_kwargs[:method]
    gridscale = plot_kwargs[:gridscale]
    ylim = plot_kwargs[:ylim]
    plot_colorbar = plot_kwargs[:plot_colorbar]
    colorrange = plot_kwargs[:colorrange]
    nan_color = plot_kwargs[:nan_color]

    # actual data interpolation
    channel_data = mean.(eachcol(dat[!, layout.data.label]))
    if method == :spherical_spline
        data = _data_interpolation_topo_spherical_spline(channel_data, layout, gridscale)
    elseif method == :multiquadratic
        data = _data_interpolation_topo_multiquadratic(channel_data, layout, gridscale)
    end

    if isnothing(ylim)
        # Make ylim symmetric around 0 for balanced topographic visualization
        data_min, data_max = extrema(data[.!isnan.(data)])
        max_abs = max(abs(data_min), abs(data_max))
        ylim = (-max_abs, max_abs)
    end

    # Set title based on user preferences and data
    if plot_kwargs[:show_title]
        if plot_kwargs[:title] != ""
            # Use user-provided title
            ax.title = plot_kwargs[:title]
        elseif hasproperty(dat, :time) && !isempty(dat.time)
            # Set default title showing time range if data has time column
            time_min, time_max = extrema(dat.time)
            time_title = @sprintf("%.3f to %.3f s", time_min, time_max)
            ax.title = time_title
        end
        ax.titlesize = plot_kwargs[:title_fontsize]
    end

    # Clear the axis 
    empty!(ax)

    # Use normalized coordinate ranges since layout coordinates are normalized to [-1, 1]
    contour_range = 1.0  # Since coordinates are normalized to [-1, 1]

    # Set contour parameters
    contour_kwargs = Dict{Symbol,Any}(:extendlow => :auto, :extendhigh => :auto, :colormap => plot_kwargs[:colormap])

    if !isnothing(nan_color)
        contour_kwargs[:nan_color] = nan_color
    end

    co = contourf!(
        range(-contour_range, contour_range, length = gridscale),
        range(-contour_range, contour_range, length = gridscale),
        data,
        levels = range(ylim[1], ylim[2], div(gridscale, 2));
        contour_kwargs...,
    )

    # Set color range on the contour plot if specified
    if !isnothing(colorrange)
        co.colorrange = colorrange
    end

    if plot_colorbar
        Colorbar(
            fig[1, 2],
            co;
            width = plot_kwargs[:colorbar_width],
            label = plot_kwargs[:colorbar_label],
            labelsize = plot_kwargs[:colorbar_fontsize],
        )
    end

    # head shape
    plot_layout_2d!(fig, ax, layout; plot_kwargs...)

    # Apply axis styling
    if plot_kwargs[:xlabel] != ""
        ax.xlabel = plot_kwargs[:xlabel]
    end
    if plot_kwargs[:ylabel] != ""
        ax.ylabel = plot_kwargs[:ylabel]
    end

    # Set font sizes
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    # Set grid properties
    ax.xgridvisible = plot_kwargs[:xgrid]
    ax.ygridvisible = plot_kwargs[:ygrid]
    if plot_kwargs[:grid_visible]
        ax.xgridvisible = true
        ax.ygridvisible = true
    end
    ax.xgridcolor = (:gray, plot_kwargs[:grid_alpha])
    ax.ygridcolor = (:gray, plot_kwargs[:grid_alpha])

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
    interactive = true,
    kwargs...,
)
    # Merge user kwargs with defaults to get all parameters
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # Override with function parameters if provided
    plot_kwargs[:display_plot] = display_plot
    plot_kwargs[:interactive] = interactive

    fig = Figure()
    ax = Axis(fig[1, 1])
    dat_subset = subset(dat, channel_selection = channel_selection, sample_selection = sample_selection)
    _plot_topography!(fig, ax, dat_subset.data, dat_subset.layout; plot_kwargs...)

    if plot_kwargs[:interactive]
        _setup_topo_interactivity!(fig, ax, dat)
    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    return fig, ax
end

"""
    plot_topography(dat::Vector{ErpData}; kwargs...)

Create topographic plots from a vector of ERP datasets by broadcasting across conditions.

# Arguments
- `dat`: Vector of ErpData objects (e.g., different conditions)
- `kwargs...`: Additional keyword arguments passed to plot_topography
"""
function plot_topography(
    dat::Vector{ErpData};
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot = true,
    interactive = true,
    kwargs...,
)
    return plot_topography.(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        display_plot = display_plot,
        interactive = interactive,
        kwargs...,
    )
end



function plot_topography!(
    fig,
    ax,
    dat::EpochData,
    epoch::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    kwargs...,
)
    plot_topography!(
        fig,
        ax,
        convert(epoch_data, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        kwargs...,
    )
end

function plot_topography(
    dat::EpochData,
    epoch::Int;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    display_plot = true,
    interactive = true,
    kwargs...,
)
    # Merge user kwargs with defaults to get all parameters
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # Override with function parameters if provided
    plot_kwargs[:display_plot] = display_plot
    plot_kwargs[:interactive] = interactive

    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_topography!(
        fig,
        ax,
        convert(epoch_data, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        plot_kwargs...,
    )

    if plot_kwargs[:interactive]
        _setup_topo_interactivity!(fig, ax, dat)
    end

    if plot_kwargs[:display_plot]
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
    radius = grid_scale / 2  # For normalized coordinates, the radius is half the grid size
    @inbounds for col = 1:grid_scale
        for row = 1:grid_scale
            x_dist = center - col
            y_dist = center - row
            if sqrt(x_dist^2 + y_dist^2) > radius
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
            if event.key == Keyboard.up
                _topo_scale_up!(ax)
            elseif event.key == Keyboard.down
                _topo_scale_down!(ax)
            end
        end
    end

    _setup_topo_selection!(fig, ax, original_data)
end

"""
    _topo_scale_up!(ax::Axis)

Increase the scale of the topographic plot (zoom in on color range).
"""
function _topo_scale_up!(ax::Axis)

    # Find the Contourf plot in the axis
    for plot in ax.scene.plots
        if plot isa Makie.Contourf

            # Get current levels
            current_levels = plot.levels[]
            if !isnothing(current_levels)

                # For zoom in: compress the range around 0 for better contrast
                level_min, level_max = extrema(current_levels)
                center = (level_min + level_max) / 2
                range_size = level_max - level_min

                # Compress the range by 20% but keep it centered
                new_range = range_size * 0.8
                new_min = center - new_range / 2
                new_max = center + new_range / 2

                # Create new levels with the same density but compressed range
                new_levels = range(new_min, new_max, length = length(current_levels))
                plot.levels[] = new_levels

                break
            else
            end
        end
    end
end

"""
    _topo_scale_down!(ax::Axis)

Decrease the scale of the topographic plot (zoom out from color range).
"""
function _topo_scale_down!(ax::Axis)

    # Find the Contourf plot in the axis
    for plot in ax.scene.plots
        if plot isa Makie.Contourf

            # Get current levels
            current_levels = plot.levels[]
            if !isnothing(current_levels)

                # For zoom out: expand the range around 0 for less contrast
                level_min, level_max = extrema(current_levels)
                center = (level_min + level_max) / 2
                range_size = level_max - level_min

                # Expand the range by 25% but keep it centered
                new_range = range_size * 1.25
                new_min = center - new_range / 2
                new_max = center + new_range / 2

                # Create new levels with the same density but expanded range
                new_levels = range(new_min, new_max, length = length(current_levels))
                plot.levels[] = new_levels

                break
            end
        end
    end
end

# =============================================================================
# REGION SELECTION FOR TOPO PLOTS
# =============================================================================

"""
    TopoSelectionState

Simple state for spatial region selection in topographic plots.
"""
mutable struct TopoSelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64,Float64,Float64}}  # x_min, y_min, x_max, y_max
    visible::Observable{Bool}
    rectangles::Vector{Makie.Poly}  # Store multiple selection rectangles
    bounds_list::Observable{Vector{Tuple{Float64,Float64,Float64,Float64}}}  # Store all selection bounds
    temp_rectangle::Union{Makie.Poly,Nothing}  # Temporary rectangle for dragging

    function TopoSelectionState(ax::Axis)
        # Initialize with empty lists for multiple selections
        new(
            Observable(false),
            Observable((0.0, 0.0, 0.0, 0.0)),
            Observable(false),
            Makie.Poly[],  # Empty vector for rectangles
            Observable{Tuple{Float64,Float64,Float64,Float64}}[],  # Empty vector for bounds
            nothing,  # No temporary rectangle initially
        )
    end
end

"""
    _setup_topo_selection!(fig::Figure, ax::Axis, original_data)

Set up simple region selection for topographic plots.
"""
function _setup_topo_selection!(fig::Figure, ax::Axis, original_data)
    selection_state = TopoSelectionState(ax)

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
        if event.button == Mouse.left
            if event.action == Mouse.press
                # Check if Shift is held down
                if shift_pressed[]
                    _start_topo_selection!(ax, selection_state, event)
                end
            elseif event.action == Mouse.release
                if selection_state.active[]
                    _finish_topo_selection!(ax, selection_state, event, original_data)
                end
            end
        end
        if event.button == Mouse.left && event.action == Mouse.press && !shift_pressed[]
            _clear_all_topo_selections!(selection_state)
        end
        if event.button == Mouse.right
            # TODO: implement right click functionality
            @info "TODO! :-)"
        end
    end

    # Handle mouse movement for selection dragging
    on(events(fig).mouseposition) do pos
        if selection_state.active[]
            _update_topo_selection!(ax, selection_state, pos)
        end
    end
end

"""
    _start_topo_selection!(ax::Axis, selection_state::TopoSelectionState, event)

Start spatial region selection in topographic plot.
"""
function _start_topo_selection!(ax::Axis, selection_state::TopoSelectionState, event)
    selection_state.active[] = true
    selection_state.visible[] = true

    # Get mouse position in screen coordinates and convert to axis coordinates
    screen_pos = events(ax.scene).mouseposition[]
    mouse_pos = mouseposition(ax)
    mouse_x, mouse_y = mouse_pos[1], mouse_pos[2]

    # Store axis coordinates for spatial selection
    selection_state.bounds[] = (mouse_x, mouse_y, mouse_x, mouse_y)

    # Create a temporary rectangle for dragging
    initial_points =
        [Point2f(mouse_x, mouse_y), Point2f(mouse_x, mouse_y), Point2f(mouse_x, mouse_y), Point2f(mouse_x, mouse_y)]
    selection_state.temp_rectangle = poly!(
        ax,
        initial_points,
        color = (:blue, 0.3),    # Blue with transparency
        strokecolor = :black,    # Black border
        strokewidth = 1,         # Thin border
        visible = true,
        overdraw = true,          # Ensure it's drawn on top
    )

    _update_topo_selection!(ax, selection_state, mouse_pos)
end

"""
    _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, event)

Finish spatial region selection in topographic plot.
"""
function _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, event, original_data = nothing)
    selection_state.active[] = false

    # Get final mouse position in screen coordinates and convert to axis coordinates
    screen_pos = events(ax.scene).mouseposition[]
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

    # Create the permanent rectangle for this selection
    bounds = selection_state.bounds[]
    start_x, start_y = bounds[1], bounds[2]
    end_x, end_y = bounds[3], bounds[4]

    rect_points = Point2f[
        Point2f(Float64(start_x), Float64(start_y)),
        Point2f(Float64(end_x), Float64(start_y)),
        Point2f(Float64(end_x), Float64(end_y)),
        Point2f(Float64(start_x), Float64(end_y)),
    ]

    # Create permanent rectangle
    permanent_rectangle = poly!(
        ax,
        rect_points,
        color = (:blue, 0.3),    # Blue with transparency
        strokecolor = :black,     # Black border
        strokewidth = 1,          # Thin border
        visible = true,
        overdraw = true,           # Ensure it's drawn on top
    )

    # Store the permanent rectangle
    push!(selection_state.rectangles, permanent_rectangle)

    # Remove the temporary rectangle
    if !isnothing(selection_state.temp_rectangle)
        delete!(selection_state.temp_rectangle.parent, selection_state.temp_rectangle)
        selection_state.temp_rectangle = nothing
    end

    # Find electrodes within ALL selected spatial regions
    all_selected_electrodes = Symbol[]
    for bounds in selection_state.bounds_list[]
        x_min, y_min, x_max, y_max = bounds
        region_electrodes = _find_electrodes_in_region(ax, x_min, y_min, x_max, y_max, original_data)
        append!(all_selected_electrodes, region_electrodes)
    end

    unique_electrodes = unique(all_selected_electrodes)
    @info "N selections: $(length(selection_state.bounds_list[])); Electrodes found: $unique_electrodes"
end

"""
    _update_topo_selection!(ax::Axis, selection_state::TopoSelectionState, pos)

Update the visual selection rectangle for spatial selection.
"""
function _update_topo_selection!(ax::Axis, selection_state::TopoSelectionState, pos)
    if selection_state.active[]
        # pos might be screen coordinates, convert to axis coordinates
        # Use mouseposition(ax) to get consistent axis coordinates
        axis_pos = mouseposition(ax)
        mouse_x, mouse_y = axis_pos[1], axis_pos[2]
        start_x, start_y = selection_state.bounds[][1], selection_state.bounds[][2]

        # Update bounds with the axis coordinates
        selection_state.bounds[] = (start_x, start_y, mouse_x, mouse_y)

        # Update the temporary rectangle during dragging
        if !isnothing(selection_state.temp_rectangle)
            bounds = selection_state.bounds[]
            start_x, start_y = bounds[1], bounds[2]
            end_x, end_y = bounds[3], bounds[4]

            # Update rectangle points for the temporary rectangle
            rect_points = Point2f[
                Point2f(Float64(start_x), Float64(start_y)),
                Point2f(Float64(end_x), Float64(start_y)),
                Point2f(Float64(end_x), Float64(end_y)),
                Point2f(Float64(start_x), Float64(end_y)),
            ]

            # Update the temporary rectangle
            selection_state.temp_rectangle[1] = rect_points
        end


    end
end

"""
    _clear_all_topo_selections!(selection_state::TopoSelectionState)

Clear all topographic selections and remove all rectangles.
"""
function _clear_all_topo_selections!(selection_state::TopoSelectionState)
    # Remove all rectangles from the scene
    for rect in selection_state.rectangles
        delete!(rect.parent, rect)
    end

    # Remove temporary rectangle if it exists
    if !isnothing(selection_state.temp_rectangle)
        delete!(selection_state.temp_rectangle.parent, selection_state.temp_rectangle)
        selection_state.temp_rectangle = nothing
    end

    # Clear the lists
    empty!(selection_state.rectangles)
    selection_state.bounds_list[] = Tuple{Float64,Float64,Float64,Float64}[]

    # Reset state
    selection_state.active[] = false
    selection_state.visible[] = false


end


"""
    _create_position_channel_map(ax::Axis, original_data=nothing)

Create a mapping from electrode coordinates to electrode labels for the topographic plot.
This uses the actual layout data from the original plot data.
"""
function _create_position_channel_map(ax::Axis, original_data = nothing)

    layout = original_data.layout
    _ensure_coordinates_2d!(layout)

    # Create mapping from electrode coordinates to labels
    position_channel_map = Dict()
    for (i, label) in enumerate(layout.data.label)
        x = layout.data.x2[i]
        y = layout.data.y2[i]
        position_channel_map[(x, y)] = Symbol(label)
    end
    return position_channel_map
end

"""
    _find_electrodes_in_region(ax::Axis, x_min::Float64, y_min::Float64, x_max::Float64, y_max::Float64, original_data=nothing)

Find electrodes within the selected spatial region using actual electrode coordinates.
This approach uses the real layout data from the topographic plot.
"""
function _find_electrodes_in_region(
    ax::Axis,
    x_min::Float64,
    y_min::Float64,
    x_max::Float64,
    y_max::Float64,
    original_data = nothing,
)
    position_channel_map = _create_position_channel_map(ax, original_data)

    if isempty(position_channel_map)
        @minimal_warning "No electrode mapping available, returning placeholder"
        return nothing
    end

    # Find electrodes within the selected region
    selected_electrodes = Symbol[]
    for ((x, y), channel) in position_channel_map
        # Check if this electrode is inside the selection rectangle
        if x_min <= x <= x_max && y_min <= y <= y_max
            push!(selected_electrodes, channel)
        end
    end

    return selected_electrodes
end
