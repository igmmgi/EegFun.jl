# Currently I have settings with :right, :below, :same for colorbar_position
# or offset positions (e.g., 1, 2). Actually, the whole colorbar stuff feels a bit off!
# MUST BE BETTER WAY TO HANDLE THIS!!!

# Single shared kwargs for all topography plots (both ICA and standard)
const PLOT_TOPOGRAPHY_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    :figure_title => ("Topography Plot", "Title for the plot window"),
    :interactive => (true, "Whether to enable interactive features"),
    :theme_fontsize => (24, "Font size for theme"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),

    # Topography parameters
    :method => (
        :thin_plate,
        "Interpolation method: :multiquadratic, :inverse_multiquadratic, :gaussian, :inverse_quadratic, :thin_plate, :polyharmonic, :shepard, :nearest, :spherical_spline. See https://eljungsk.github.io/ScatteredInterpolation.jl/dev/methods/ for details on methods.",
    ),
    :colormap => (:jet, "Colormap for the topography"),
    :gridscale => (75, "Grid resolution for interpolation"),
    :dims => (nothing, "Grid dimensions (rows, cols). If nothing, calculates best square-ish grid"),
    :ylim => (nothing, "Y-axis limits (nothing for auto). For ICA plots, use num_levels instead."),
    :num_levels => (20, "Number of contour levels (for ICA plots). For standard plots, use ylim instead."),

    # Title parameters
    :title => ("", "Plot title"),
    :title_fontsize => (16, "Font size for the title"),
    :show_title => (true, "Whether to show the title"),
    :time_unit =>
        (:s, "Time unit for display labels (:s or :ms). Only affects title strings — all intervals and selections remain in seconds."),

    # Head shape parameters
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

    # Colorbar parameters - get all Colorbar attributes with their actual defaults
    [
        Symbol("colorbar_$(attr)") => (get(COLORBAR_DEFAULTS, attr, nothing), "Colorbar $(attr) parameter") for
        attr in propertynames(Colorbar)
    ]...,

    # Override specific colorbar parameters with custom defaults
    :colorbar_plot => (true, "Whether to display the colorbar"),
    :colorbar_position => ((1, 2), "Colorbar position as (row, col) tuple, or :right, :below"),
    :colorbar_label => ("μV", "Label for the colorbar"),
    :colorbar_plot_numbers => ([], "Plot indices for which to show colorbars. Empty list shows colorbars for all plots."),

    # Channel highlighting
    :highlight_channels => (
        nothing,
        "Highlight channel groups as a NamedTuple or Vector of NamedTuples. Each group: (channels=[:Cz, :Pz], color=:white, size=8, marker=:circle)",
    ),

    # ICA-specific parameters (ignored for standard topography plots)
    :use_global_scale => (false, "Do multiple ICA topoplots share the same color scale based on min/max across all components?"),
    :component_selection => (components(), "Function that returns boolean vector for component filtering"),
)


"""
    _plot_topography!(fig::Figure, ax::Axis, ica::InfoIca, component::Int; kwargs...)

Internal function to plot a single ICA component topography on existing figure/axis.

# Arguments
- `fig::Figure`: The Figure object
- `ax::Axis`: The Axis object  
- `ica::InfoIca`: The ICA result object (contains layout information).
- `component::Int`: The component index to plot (1-based).

# Keyword Arguments
$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- `co`: The contour plot object returned by `contourf!`.
"""
function _plot_topography!(fig::Figure, ax::Axis, ica::InfoIca, component::Int; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # Extract commonly used kwargs
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)

    # Ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(ica.layout)
    _ensure_coordinates_3d!(ica.layout)

    # Prepare data for this component (now returns bounds too)
    data, x_bounds, y_bounds = _prepare_ica_topo_data(ica, component, method, gridscale)

    # Calculate levels
    levels = _calculate_topo_levels(data; num_levels = pop!(plot_kwargs, :num_levels))

    # Create contour plot with adaptive bounds
    co = contourf!(
        ax,
        range(x_bounds[1], x_bounds[2], length = size(data, 1)),
        range(y_bounds[1], y_bounds[2], length = size(data, 2)),
        data,
        levels = levels;
        colormap = pop!(plot_kwargs, :colormap),
        nan_color = :transparent,
    )

    # Add colorbar if requested
    # If colorbar_plot_numbers is empty, show colorbar for all components
    # Otherwise, only show for components in the list
    colorbar_plot_numbers = plot_kwargs[:colorbar_plot_numbers]
    should_show_colorbar = plot_kwargs[:colorbar_plot] && (isempty(colorbar_plot_numbers) || component in colorbar_plot_numbers)
    if should_show_colorbar
        colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
        colorbar_position = pop!(plot_kwargs, :colorbar_position, (1, 2))
        Colorbar(fig[colorbar_position...], co; colorbar_kwargs..., tellwidth = true, tellheight = false)
    end

    # Draw smooth circle to hide jagged interpolation edge
    _draw_smooth_circle_mask!(ax, x_bounds, y_bounds)

    # Add head shape and electrode markers
    plot_layout_2d!(fig, ax, ica.layout; plot_kwargs...)

    return co
end

"""
    plot_topography(ica::InfoIca; ...)

Plot multiple ICA component topographies in a grid layout within a new Figure.

# Arguments
- `ica::InfoIca`: The ICA result object (contains layout information).

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- `fig::Figure`: The generated Makie Figure containing the grid of topoplots.

# Examples

## Basic Usage
```julia
# Plot first 10 components (default)
fig = plot_topography(ica)

# Plot specific range of components
fig = plot_topography(ica, component_selection = components(5:15))

# Plot specific components
fig = plot_topography(ica, component_selection = components([1, 3, 5, 7]))

# Plot all components (if screen can handle it)
fig = plot_topography(ica, component_selection = components())
```

## Advanced Selection
```julia
# Plot components with custom selection
fig = plot_topography(ica, 
    component_selection = components(1:10)  # First 10 components
)

# Plot even-numbered components
fig = plot_topography(ica, 
    component_selection = components(2:2:20)  # Even components 2, 4, 6, ..., 20
)
```
"""
function plot_topography(ica::InfoIca; component_selection = components(), kwargs...)

    # ICA weights are dimensionless — override the default "μV" label
    merged_kwargs = Dict{Symbol,Any}(:colorbar_label => "a.u.", pairs(kwargs)...)

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, merged_kwargs)
    dims = pop!(plot_kwargs, :dims)
    display_plot = pop!(plot_kwargs, :display_plot)

    # ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(ica.layout)
    _ensure_coordinates_3d!(ica.layout)

    # Get selected components using the helper function
    comps = get_selected_components(ica, component_selection)
    if isempty(comps)
        @minimal_warning "No components selected for topography plot"
        return (fig = Figure(figure_padding = plot_kwargs[:figure_padding]),)
    end

    # Get colorbar settings to adjust grid if needed
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot)

    # Create figure
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])
    fig = Figure(figure_padding = plot_kwargs[:figure_padding])

    # Deal with plot dimensions
    isnothing(dims) && (dims = _best_rect(length(comps)))

    # Validate dimensions
    if length(dims) != 2 || any(dims .<= 0)
        throw(ArgumentError("Invalid dimensions: $dims. Expected [rows, cols] with positive values."))
    end

    # Ensure we have enough grid cells
    total_cells = dims[1] * dims[2]
    if total_cells < length(comps)
        throw(ArgumentError("Grid dimensions $dims provide $total_cells cells but need $(length(comps))."))
    end

    # Create axes and plot each component
    for i in eachindex(comps)
        # Calculate base row/col indices
        base_row, base_col = divrem(i - 1, dims[2]) .+ (1, 1)

        # Get colorbar position for this component
        colorbar_position = get(plot_kwargs, :colorbar_position, :right)

        # Convert symbol to tuple or use tuple directly
        if colorbar_position isa Symbol
            if colorbar_position == :right
                colorbar_offset = (1, 2)
            elseif colorbar_position == :below
                colorbar_offset = (2, 1)
            elseif colorbar_position == :same
                colorbar_offset = (1, 1)
            else
                throw(ArgumentError("colorbar_position must be :right, :below, :same, or a tuple (row, col), got: $colorbar_position"))
            end
        elseif colorbar_position isa Tuple
            # User provided tuple directly (row_offset, col_offset)
            colorbar_offset = colorbar_position
        else
            throw(ArgumentError("colorbar_position must be :right, :below, :same, or a tuple (row, col), got: $colorbar_position"))
        end

        # Calculate plot and colorbar positions
        if colorbar_plot
            if colorbar_offset[1] < colorbar_offset[2]
                # Colorbars to the right: each component needs 2 columns
                plot_row = base_row
                plot_col = (base_col - 1) * 2 + 1
                colorbar_row = plot_row + colorbar_offset[1] - 1
                colorbar_col = plot_col + colorbar_offset[2] - 1
            else
                # Colorbars below: each component needs 2 rows
                plot_row = (base_row - 1) * 2 + 1
                plot_col = base_col
                colorbar_row = plot_row + colorbar_offset[1] - 1
                colorbar_col = plot_col + colorbar_offset[2] - 1
            end
        else
            plot_col = base_col
            colorbar_row = base_row
            colorbar_col = base_col
        end

        # Create axis with title
        if colorbar_plot
            ax = Axis(fig[plot_row, plot_col], title = @sprintf("IC %d (%.1f%%)", comps[i], ica.variance[comps[i]] * 100))
        else
            ax = Axis(fig[base_row, base_col], title = @sprintf("IC %d (%.1f%%)", comps[i], ica.variance[comps[i]] * 100))
        end

        # Use the internal plotting function with colorbar position
        _plot_topography!(
            fig,
            ax,
            ica,
            comps[i];
            plot_kwargs...,
            colorbar_plot = colorbar_plot,
            colorbar_position = (colorbar_row, colorbar_col),
        )
    end

    # Add keyboard event handling for scaling
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key in (Keyboard.up, Keyboard.down)
            axes = filter(ax -> ax isa Axis, fig.content)
            if event.key == Keyboard.up
                _topo_scale_up!.(axes, plot_kwargs[:zoom_step])
            else
                _topo_scale_down!.(axes, plot_kwargs[:zoom_step])
            end
        end
    end

    if display_plot
        _display_figure(fig)
    end

    _set_window_title("Makie")
    return (fig = fig,)

end


"""Mutable state for the interactive ICA component activation viewer."""
mutable struct IcaOverlayState
    names::Vector{Symbol}
    plots::Vector{Vector{Any}}
end

mutable struct IcaArtifactState
    data::Union{Nothing,ArtifactComponents}
    toggles::Union{Nothing,Dict{Symbol,Observable{Bool}}}
    textbox::Union{Nothing,Textbox}
end

mutable struct IcaComponentState
    # Data
    dat::ContinuousData
    ica::InfoIca
    component_data::Matrix{Float64}

    # View settings
    n_visible_components::Int
    window_size::Int

    # Observables
    comp_start::Observable{Int}
    xrange::Observable{UnitRange{Int}}
    ylims::Observable{Tuple{Float64,Float64}}
    channel_data::Observable{Vector{Float64}}
    show_channel::Observable{Bool}
    channel_yscale::Observable{Float64}
    use_global_scale::Observable{Bool}
    invert_scale::Observable{Bool}

    # Components to display
    components::Vector{Int}

    # Plot parameters
    plot_kwargs::Dict{Symbol,Any}

    # Plot elements
    axs::Vector{Axis}
    channel_axs::Vector{Union{Axis,Nothing}}  # Allow for nothing values
    topo_axs::Vector{Axis}
    lines_obs::Vector{Observable{Vector{Float64}}}

    # Sub-states
    overlay::IcaOverlayState
    artifacts::IcaArtifactState

    # Constructor updated to accept and store topo_kwargs
    function IcaComponentState(
        dat,
        ica,
        component_selection,
        n_visible_components,
        window_size;
        method = :thin_plate,
        gridscale = 100,  # Grid resolution for interpolation
        kwargs...,
    )
        # Prepare data matrix
        dat_matrix = _prepare_ica_data_matrix(dat, ica)
        component_data = ica.unmixing * dat_matrix
        total_components = size(component_data, 1)

        # Create observables
        comp_start = Observable(1)
        use_global_scale = Observable(false)
        invert_scale = Observable(false)

        # Find index closest to time 0 to center the initial view 
        time_zero_idx = searchsortedfirst(dat.data.time, 0.0)
        time_zero_idx = clamp(time_zero_idx, 1, length(dat.data.time))
        half_window = div(window_size, 2)
        start_idx = max(1, time_zero_idx - half_window)
        end_idx = min(size(component_data, 2), start_idx + window_size - 1)

        # Adjust start_idx if end_idx reached the boundary
        if end_idx == size(component_data, 2)
            start_idx = max(1, end_idx - window_size + 1)
        end

        xrange = Observable(start_idx:end_idx)

        # Use the component selection predicate to get the actual component indices
        component_mask = component_selection(1:total_components)
        comps_to_use = findall(component_mask)  # Convert boolean mask to actual indices

        # Calculate initial y-range based on components we'll show
        if !isempty(comps_to_use) && all(idx -> idx <= total_components, comps_to_use)
            initial_range_data = component_data[comps_to_use, start_idx:end_idx]
        else
            initial_range_data = zeros(0, length(start_idx:end_idx)) # Empty if no valid components 
        end
        if !isempty(initial_range_data)
            initial_range = maximum(abs.(extrema(initial_range_data)))
        else
            initial_range = 1.0 # Default range if no data
        end
        ylims = Observable((-initial_range, initial_range))
        channel_data = Observable(zeros(size(dat.data, 1)))
        show_channel = Observable(false)
        channel_yscale = Observable(1.0)

        # Store plot kwargs directly, including the extracted parameters
        plot_kwargs = copy(kwargs)
        plot_kwargs[:method] = method
        plot_kwargs[:gridscale] = gridscale

        # Initialize empty plot element arrays
        axs = Vector{Axis}()
        channel_axs = Vector{Union{Axis,Nothing}}()
        topo_axs = Vector{Axis}()
        lines_obs = Vector{Observable{Vector{Float64}}}()
        overlay = IcaOverlayState(Symbol[], Vector{Any}[])
        artifacts = IcaArtifactState(nothing, nothing, nothing)

        new(
            dat,
            ica,
            component_data,
            n_visible_components,
            window_size,
            comp_start,
            xrange,
            ylims,
            channel_data,
            show_channel,
            channel_yscale,
            use_global_scale,
            invert_scale,
            comps_to_use,
            plot_kwargs,
            axs,
            channel_axs,
            topo_axs,
            lines_obs,
            overlay,
            artifacts,
        )
    end
end


"""
    plot_ica_component_activation(dat::ContinuousData, ica::InfoIca; ...)

# Interactive Viewer for ICA Component Activation

# Arguments
- `dat::ContinuousData`: Continuous EEG data object.
- `ica::InfoIca`: ICA results object.

# Keyword Arguments
- `artifacts::Union{Nothing,ArtifactComponents}=nothing`: Optional artifact detection results. When provided, adds category checkboxes for filtering components.
- `component_selection::Function=components(:all)`: Predicate to select which components to display.
- `n_visible_components::Int=10`: Number of components visible at once (auto-adjusted for selected components).
- `window_size::Int=2000`: Initial time window size in samples.
- `topo_kwargs::Dict=Dict()`: Keyword arguments passed down for topography plots (see `_plot_topo_on_axis!`).
- `head_kwargs::Dict=Dict()`: Keyword arguments passed down for head outlines.
- `point_kwargs::Dict=Dict()`: Keyword arguments passed down for channel markers.
- `label_kwargs::Dict=Dict()`: Keyword arguments passed down for channel labels.

# Returns
- `fig::Figure`: The Makie Figure object containing the interactive plot.
"""
function plot_ica_component_activation(dat::ContinuousData, ica::InfoIca; artifacts::Union{Nothing,ArtifactComponents} = nothing, kwargs...)
    # Generate window title from dataset
    title_str = _generate_window_title(dat)
    _set_window_title(title_str)

    # ICA weights are dimensionless — override the default "μV" label
    merged_kwargs = Dict{Symbol,Any}(:colorbar_label => "a.u.", pairs(kwargs)...)

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, merged_kwargs)

    # ensure coordinates are 2d
    _ensure_coordinates_2d!(ica.layout)
    _ensure_coordinates_3d!(ica.layout)

    # Extract commonly used kwargs that are NOT plot call specific from plot_kwargs
    component_selection = pop!(plot_kwargs, :component_selection)

    # Get selected components to determine layout
    component_mask = component_selection(1:length(ica.ica_label))
    selected_components = findall(component_mask)

    # Adapt n_visible_components to number of selected components (max 10 for scrollable view)
    n_visible_components = min(length(selected_components), 10)

    window_size = min(2000, n_samples(dat))
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    display_plot = pop!(plot_kwargs, :display_plot)

    # change some defaults
    plot_kwargs[:point_plot] = false
    plot_kwargs[:label_plot] = false

    # Pass kwargs to constructor
    state = IcaComponentState(
        dat,
        ica,
        component_selection,
        n_visible_components,
        window_size;
        method = method,
        gridscale = gridscale,
        plot_kwargs...,
    )

    # Create figure and UI
    fig = Figure(figure_padding = plot_kwargs[:figure_padding])
    _create_component_activation_plots!(fig, state)
    _add_navigation_controls!(fig, state)
    _add_navigation_sliders!(fig, state)
    _add_channel_menu!(fig, state)

    # Add artifact category checkboxes if artifacts provided
    if !isnothing(artifacts)
        _add_artifact_category_checkboxes!(fig, state, artifacts)
    end

    _setup_keyboard_interactions!(fig, state)

    colsize!(fig.layout, 1, Relative(0.15))
    colsize!(fig.layout, 2, Relative(0.85))
    rowgap!(fig.layout, 5)  # Tight spacing between component rows
    rowgap!(fig.layout, state.n_visible_components, 80)  # Larger gap before controls

    if display_plot
        _display_figure(fig)
    end
    _set_window_title("Makie")
    return (fig = fig,)
end



"""
    _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca)

Selects, centers, scales, and transposes data for ICA unmixing.

# Arguments
- `dat::ContinuousData`: Input EEG data.
- `ica::InfoIca`: Corresponding ICA result containing labels and scaling factors.

# Returns
- `Matrix{Float64}`: Data matrix ready for `ica.unmixing * dat_matrix`. (channels x samples)
"""
function _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca)
    dat_matrix = permutedims(Matrix(dat.data[!, ica.layout.data.label]))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica.scale
    return dat_matrix
end


# Internal helper to calculate contour levels
"""
    _calculate_topo_levels(data::AbstractMatrix{<:Real}; ...)

Calculate appropriate levels for a contour plot based on data values.

Handles local vs. global scaling, ignores NaNs, and ensures the range is non-zero.

# Arguments
- `data::AbstractMatrix{<:Real}`: The 2D data matrix for the contour plot.

# Keyword Arguments
- `use_global_scale::Bool=false`: If true, use `global_min` and `global_max` for scaling.
- `global_min=nothing`: The global minimum value (used if `use_global_scale` is true).
- `global_max=nothing`: The global maximum value (used if `use_global_scale` is true).
- `num_levels::Int=20`: The desired number of contour levels.

# Returns
- `AbstractRange`: A range object suitable for the `levels` argument of `contourf!`.
"""
function _calculate_topo_levels(
    data::AbstractMatrix{<:Real};
    use_global_scale = false,
    global_min = nothing,
    global_max = nothing,
    num_levels = 20,
)
    # Input validation
    num_levels <= 0 && throw(ArgumentError("num_levels must be positive, got $num_levels"))

    # Find local min/max, ignoring NaNs
    local_min, local_max = -1.0, 1.0
    valid_values = filter(!isnan, data)
    if !isempty(valid_values)
        local_min = minimum(valid_values)
        local_max = maximum(valid_values)
    end

    # Determine final min/max for levels
    if use_global_scale && !isnothing(global_min) && !isnothing(global_max)
        final_min, final_max = global_min, global_max
    else
        final_min, final_max = local_min, local_max
    end

    # Ensure min != max to avoid range error
    if final_min == final_max
        final_min -= 0.1
        final_max += 0.1
    end

    return range(final_min, final_max, length = num_levels)
end

# Internal helper to plot topo data and head shape on a given axis
"""
    _plot_topo_on_axis!(ax::Axis, fig::Figure, data::AbstractMatrix{<:Real}, layout::DataFrame, levels; ...)

Draws the topographic contour plot and head shape onto a specified `Axis`.

# Arguments
- `ax::Axis`: The Makie Axis to plot on.
- `fig::Figure`: The parent Figure (needed for `plot_layout_2d!`).
- `data::AbstractMatrix{<:Real}`: The 2D interpolated topographic data.
- `layout::DataFrame`: The channel layout DataFrame containing positions.
- `levels`: The contour levels (typically calculated by `_calculate_topo_levels`).

# Keyword Arguments
- `gridscale::Int=300`: Resolution of the interpolation grid.
- `colormap=:jet`: Colormap for the contour plot.
- `head_kwargs=Dict()`: Keyword arguments passed to `plot_layout_2d!` for the head outline.
- `point_kwargs=Dict()`: Keyword arguments passed to `plot_layout_2d!` for channel points.
- `label_kwargs=Dict()`: Keyword arguments passed to `plot_layout_2d!` for channel labels.

# Returns
- `co`: The contour plot object returned by `contourf!`.
"""
function _plot_topo_on_axis!(
    ax::Axis,
    fig::Figure,
    data::AbstractMatrix{<:Real},
    layout::Layout,
    levels;
    gridscale = 75,
    colormap = :jet,
    head_color = :black,
    head_linewidth = 2,
    head_radius = 1.0,
    point_plot = false,
    point_marker = :circle,
    point_markersize = 12,
    point_color = :black,
    label_plot = false,
    label_fontsize = 20,
    label_color = :black,
    label_xoffset = 0,
    label_yoffset = 0,
    kwargs...,
)

    # Validate gridscale
    gridscale <= 0 && throw(ArgumentError("gridscale must be positive, got $gridscale"))

    # Use adaptive bounds from prepare_ica_topo_data
    x_bounds = kwargs[:x_bounds]
    y_bounds = kwargs[:y_bounds]
    x_range = range(x_bounds[1], x_bounds[2], length = size(data, 1))
    y_range = range(y_bounds[1], y_bounds[2], length = size(data, 2))
    co = contourf!(ax, x_range, y_range, data; levels = levels, colormap = colormap, nan_color = :transparent)

    # Draw smooth circle to hide jagged interpolation edge
    _draw_smooth_circle_mask!(ax, x_bounds, y_bounds)

    plot_layout_2d!(
        fig,
        ax,
        layout;
        head_color = head_color,
        head_linewidth = head_linewidth,
        head_radius = head_radius,
        point_plot = point_plot,
        point_marker = point_marker,
        point_markersize = point_markersize,
        point_color = point_color,
        label_plot = label_plot,
        label_fontsize = label_fontsize,
        label_color = label_color,
        label_xoffset = label_xoffset,
        label_yoffset = label_yoffset,
    )

    hidedecorations!(ax, grid = false, label = false)

    return co

end

"""Plot a single ICA topoplot inside the component activation viewer."""
function _plot_ica_topo_in_viewer!(
    fig,
    topo_ax,
    ica,
    comp_idx;
    use_global_scale = false,
    gridscale = 100,
    colormap = :jet,
    num_levels = 20,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    method = :spherical_spline,
    pre_calculated_levels = nothing,
    kwargs...,
)
    # Prepare data using the new internal function (returns tuple: data, x_bounds, y_bounds)
    topo_data, x_bounds, y_bounds = _prepare_ica_topo_data(ica, comp_idx, method, gridscale)

    # Use pre-calculated levels if provided, otherwise calculate levels
    if isnothing(pre_calculated_levels)
        levels = _calculate_topo_levels(topo_data; use_global_scale = use_global_scale, num_levels = num_levels)
    else
        levels = pre_calculated_levels
    end

    # Convert dictionary kwargs to individual parameters
    head_color = get(head_kwargs, :head_color, :black)
    head_linewidth = get(head_kwargs, :head_linewidth, 2)
    head_radius = get(head_kwargs, :head_radius, 1.0)

    point_plot = get(point_kwargs, :point_plot, false)
    point_marker = get(point_kwargs, :point_marker, :circle)
    point_markersize = get(point_kwargs, :point_markersize, 12)
    point_color = get(point_kwargs, :point_color, :black)

    label_plot = get(label_kwargs, :label_plot, false)
    label_fontsize = get(label_kwargs, :label_fontsize, 20)
    label_color = get(label_kwargs, :label_color, :black)
    label_xoffset = get(label_kwargs, :label_xoffset, 0)
    label_yoffset = get(label_kwargs, :label_yoffset, 0)

    # Plot using the generic topo function
    co = _plot_topo_on_axis!(
        topo_ax,
        fig,
        topo_data,
        ica.layout,
        levels;
        x_bounds = x_bounds,
        y_bounds = y_bounds,
        gridscale = gridscale,
        colormap = colormap,
        head_color = head_color,
        head_linewidth = head_linewidth,
        head_radius = head_radius,
        point_plot = point_plot,
        point_marker = point_marker,
        point_markersize = point_markersize,
        point_color = point_color,
        label_plot = label_plot,
        label_fontsize = label_fontsize,
        label_color = label_color,
        label_xoffset = label_xoffset,
        label_yoffset = label_yoffset,
        kwargs...,
    )
    hidedecorations!(topo_ax, grid = false, label = false)

    return co
end

"""Create time-series and topoplot axes for each visible component."""
function _create_component_activation_plots!(fig, state)

    # Number of plots
    if length(state.components) == size(state.component_data, 1)
        num_plots = state.n_visible_components
    else
        num_plots = length(state.components)  # Show all components in the subset
    end

    for i = 1:num_plots

        topo_ax = Axis(fig[i, 1], aspect = DataAspect(), titlevisible = false)
        push!(state.topo_axs, topo_ax)

        # Get the actual component number
        comp_idx = _get_component_index(state, i)

        # Time series axis creation (now on the right)
        is_last_component = (i == state.n_visible_components)

        ax = Axis(
            fig[i, 2],
            ylabel = @sprintf("IC %d", comp_idx),
            yaxisposition = :left,

            # Tick visibility
            yticklabelsvisible = false,
            yticksvisible = true,
            xticklabelsvisible = is_last_component,
            xticksvisible = is_last_component,

            # Grid settings
            xgridvisible = false,
            ygridvisible = false,
            xminorgridvisible = false,
            yminorgridvisible = false,

            # Spacing
            ylabelpadding = 0.0,
            yticklabelpad = 0.0,
            yticklabelspace = 0.0,
            yautolimitmargin = (0, 0),
        )
        push!(state.axs, ax)

        # Always create channel overlay axis (keep spines hidden)
        ax_channel = Axis(
            fig[i, 2],
            yaxisposition = :right,
            xaxisposition = :top,

            # Tick visibility - all hidden
            yticklabelsvisible = false,
            yticksvisible = false,
            xticklabelsvisible = false,
            xticksvisible = false,

            # Grid settings - all disabled
            xgridvisible = false,
            ygridvisible = false,
            xminorgridvisible = false,
            yminorgridvisible = false,

            # Spine visibility - all hidden for overlay effect
            bottomspinevisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            leftspinevisible = false,
        )
        push!(state.channel_axs, ax_channel)
        push!(state.overlay.plots, Any[])  # empty overlay list for this axis

        linkaxes!(ax, ax_channel)

        # Set initial x-axis limits (overlay lines are added dynamically via _update_channel_overlays!)
        xlims!(ax_channel, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))

        # Observable creation for component plot
        # Handle potential invalid comp_idx robustly
        if comp_idx <= size(state.component_data, 1)
            comp_data = state.component_data[comp_idx, :]
        else
            comp_data = zeros(Float64, size(state.component_data, 2)) # Placeholder data
        end
        lines_obs = Observable(comp_data)
        push!(state.lines_obs, lines_obs)

        # Component line plot
        lines!(ax, @lift(state.dat.data.time[$(state.xrange)]), @lift($(lines_obs)[$(state.xrange)]), color = :black)

        # Set initial x-axis limits for component plot
        xlims!(ax, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))

        # Create the topo plot using the dedicated viewer function
        if comp_idx <= size(state.component_data, 1)
            _plot_ica_topo_in_viewer!(fig, topo_ax, state.ica, comp_idx; use_global_scale = state.use_global_scale[], state.plot_kwargs...)
            topo_ax.ylabel = @sprintf("IC %d\n(%.1f%%)", comp_idx, state.ica.variance[comp_idx] * 100)
            topo_ax.ylabelvisible = true
        else
            empty!(topo_ax) # Clear axis if comp_idx is invalid
            topo_ax.ylabel = @sprintf("Invalid IC %d", comp_idx)
        end

        # Add left-click handler to plot component spectrum (only when clicking on the topoplot)
        # Store the plot index i to dynamically get the current component index (handles page navigation)
        plot_index = i
        on(events(topo_ax).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                # Get mouse position in axis coordinates
                mouse_pos = mouseposition(topo_ax)
                mouse_x, mouse_y = mouse_pos[1], mouse_pos[2]

                # Get head radius from plot kwargs (default is 1.0)
                head_radius = get(state.plot_kwargs, :head_radius, 1.0)

                # Check if click is within the head circle (centered at 0,0 with radius head_radius)
                distance_from_center = sqrt(mouse_x^2 + mouse_y^2)
                if distance_from_center <= head_radius
                    # Dynamically get the current component index for this plot position
                    # This ensures correct component is selected even after page navigation
                    current_comp_idx = _get_component_index(state, plot_index)
                    plot_ica_component_spectrum(state.dat, state.ica, component_selection = components(current_comp_idx))
                end
            end
        end
    end
end

"""Add navigation buttons, component textbox, and scale checkboxes."""
function _add_navigation_controls!(fig, state)

    # Add navigation buttons below topo plots in column 1
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1], tellheight = false)

    # Column 1 spacing 
    colsize!(topo_nav, 1, 40) # Set width of the empty first column

    # Navigation buttons now in row 1, columns 2 & 3
    prev_topo = Button(topo_nav[1, 2], label = "◄ Previous", tellheight = false)
    next_topo = Button(topo_nav[1, 3], label = "Next ►", tellheight = false)

    # Component selection now in row 2, columns 2, 3, 4
    text_label = Label(topo_nav[2, 2], "Components:", tellheight = false, halign = :right)
    text_input = Textbox(topo_nav[2, 3], placeholder = "e.g. 1,3-5,8", tellheight = false)
    apply_button = Button(topo_nav[2, 4], label = "Apply", tellheight = false)

    # Store textbox in state for access by other UI elements
    state.artifacts.textbox = text_input

    # Global scale checkbox now in row 3, columns 2 & 3
    global_scale_check = Checkbox(topo_nav[3, 2], checked = state.use_global_scale[], tellheight = false)
    Label(topo_nav[3, 3], "Use Global Scale", tellwidth = false, tellheight = false)

    # Invert scale checkbox now in row 4, columns 2 & 3
    invert_scale_check = Checkbox(topo_nav[4, 2], checked = state.invert_scale[], tellheight = false)
    Label(topo_nav[4, 3], "Invert Scale", tellwidth = false, tellheight = false)

    # Add column/row gaps for better spacing
    colgap!(topo_nav, 2, 10)
    colgap!(topo_nav, 3, 5)
    rowgap!(topo_nav, 1, 35)
    rowgap!(topo_nav, 2, 45)
    rowgap!(topo_nav, 3, 35)

    # Connect checkboxes to state
    on(global_scale_check.checked) do checked
        state.use_global_scale[] = checked
        _update_components!(state)
    end

    on(invert_scale_check.checked) do checked
        state.invert_scale[] = checked
        _update_components!(state)
    end

    # Connect navigation buttons
    on(prev_topo.clicks) do _
        new_start = max(1, state.comp_start[] - state.n_visible_components)
        state.comp_start[] = new_start
        _update_components!(state)
    end

    on(next_topo.clicks) do _
        new_start = min(size(state.component_data, 1) - state.n_visible_components + 1, state.comp_start[] + state.n_visible_components)
        state.comp_start[] = new_start
        _update_components!(state)
    end

    # Connect apply button
    on(apply_button.clicks) do _
        text_value = text_input.displayed_string[]
        # Skip if empty or if it's the placeholder text
        if !isempty(text_value) && !startswith(text_value, "e.g.")
            comps = _parse_string_to_ints(text_value, size(state.component_data, 1))
            if !isempty(comps)
                @info "Creating new plot with components: $comps"
                plot_ica_component_activation(state.dat, state.ica, component_selection = components(comps))
            end
        end
    end
end

"""Add a position slider for x-axis scrolling."""
function _add_navigation_sliders!(fig, state)
    # Create new row for position slider below the navigation buttons
    slider_row = state.n_visible_components + 3  # Move slider to a new row

    # Calculate the step size for the slider (1% of the data length)
    step_size = max(1, div(length(state.dat.data.time), 100))

    # Use a more compact style
    x_slider = Slider(
        fig[slider_row, 2],
        range = 1:step_size:length(state.dat.data.time),
        startvalue = first(state.xrange[]),
        tellwidth = false,
        tellheight = false,
        width = Auto(),
    )

    on(x_slider.value) do x
        _update_view_range!(state, Int(round(x)))
    end
end

"""Update `xrange` and all axis x-limits from a new start position."""
function _update_view_range!(state, start_pos)
    # Ensure we stay within data bounds
    if start_pos + state.window_size > length(state.dat.data.time)
        start_pos = length(state.dat.data.time) - state.window_size + 1
    end
    end_pos = start_pos + state.window_size - 1

    # Update range
    state.xrange[] = start_pos:end_pos

    # Update axis limits based on the current view range
    first_idx = clamp(first(state.xrange[]), 1, length(state.dat.data.time))
    last_idx = clamp(last(state.xrange[]), 1, length(state.dat.data.time))
    new_xlims = (state.dat.data.time[first_idx], state.dat.data.time[last_idx])

    # Update all axes including channel axes
    for ax in state.axs
        xlims!(ax, new_xlims)
    end
    for ax in state.channel_axs
        if !isnothing(ax)  # Only update if the axis exists
            xlims!(ax, new_xlims)
        end
    end

end

# Colour palette for channel overlays — cycles when more channels are selected.
const _CHANNEL_OVERLAY_COLORS = [:grey, :dodgerblue, :darkorange, :mediumpurple, :forestgreen, :crimson, :goldenrod, :deeppink]

"""Open a popup window with a checkbox grid for selecting additional channel overlays."""
function _show_channel_popup(state)
    all_channels = Symbol.(names(state.dat.data))

    menu_fig = Figure(size = (900, 700))
    Label(menu_fig[1, 1], "Select Additional Channels", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_channels, 9, (ch, _) -> ch in state.overlay.names)  # pre-select active channels

    _add_group_buttons!(menu_fig, 3, cbs, all_channels, [("None (Clear)", _ -> false)])

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 150)
    on(btn_apply.clicks) do _
        selected = all_channels[findall(cb -> cb.checked[], cbs)]
        _update_channel_overlays!(state, selected)
    end

    _display_popup(menu_fig)
end

"""Add a button that opens the channel-overlay popup window."""
function _add_channel_menu!(fig, state)
    menu_row = state.n_visible_components + 3
    btn = Button(fig[menu_row, 1], label = "Additional Channels", fontsize = 18, tellheight = false)
    on(btn.clicks) do _
        _show_channel_popup(state)
    end
    return btn
end

"""Remove all channel overlay plots from every component axis."""
function _clear_channel_overlays!(state)
    for i in eachindex(state.overlay.plots)
        ax_ch = state.channel_axs[i]
        if !isnothing(ax_ch)
            for p in state.overlay.plots[i]
                delete!(ax_ch, p)
            end
        end
        empty!(state.overlay.plots[i])
    end
    empty!(state.overlay.names)
end

"""Apply a list of channel overlays — each channel gets a distinct colour."""
function _update_channel_overlays!(state, selected_channels::Vector{Symbol})
    _clear_channel_overlays!(state)

    # Store names for re-opening popup with pre-selected state
    append!(state.overlay.names, selected_channels)

    for (ch_idx, ch_name) in enumerate(selected_channels)
        color = _CHANNEL_OVERLAY_COLORS[mod1(ch_idx, length(_CHANNEL_OVERLAY_COLORS))]
        col_data = state.dat.data[!, ch_name]
        is_bool = eltype(col_data) == Bool

        for i = 1:min(state.n_visible_components, length(state.channel_axs))
            ax_ch = state.channel_axs[i]
            isnothing(ax_ch) && continue

            if is_bool
                # Boolean channel → vertical lines at all true positions
                true_times = state.dat.data.time[findall(col_data)]
                if !isempty(true_times)
                    p = vlines!(ax_ch, true_times, color = (color, 0.5), linewidth = 1)
                    push!(state.overlay.plots[i], p)
                end
            else
                # Numeric channel → reactive line that follows xrange scrolling
                ch_vec = Vector{Float64}(col_data)
                p = lines!(
                    ax_ch,
                    @lift(state.dat.data.time[$(state.xrange)]),
                    @lift(ch_vec[$(state.xrange)] .* $(state.channel_yscale)),
                    color = color,
                )
                push!(state.overlay.plots[i], p)
            end
        end
    end
end


"""
    _add_artifact_category_checkboxes!(fig, state, artifacts)

Add checkboxes to filter components by artifact category.
Sets up artifact state and creates UI checkboxes in the navigation panel.
function _get_artifact_categories(artifacts::ArtifactComponents)
    categories = []
    for prop in propertynames(artifacts)
        val = getproperty(artifacts, prop)
        if val isa Dict
            for (k, v) in val
                push!(categories, (k, string(k), v))
            end
        else
            lbl = join(titlecase.(split(string(prop), "_")), " ")
            push!(categories, (prop, lbl, val))
        end
    end
    return categories
end

"""
function _add_artifact_category_checkboxes!(fig, state, artifacts)
    # Store artifacts in state
    state.artifacts.data = artifacts
    state.artifacts.toggles = Dict{Symbol,Observable{Bool}}()

    # Access the topo_nav GridLayout
    topo_nav = content(fig[state.n_visible_components+1, 1])

    categories = _get_artifact_categories(artifacts)
    start_row = 5  # Start after row 4 (Invert Scale)

    for (idx, (key, label, comps)) in enumerate(categories)
        row = start_row + idx - 1
        comp_str = isempty(comps) ? "none" : join(string.(comps), ", ")

        # Create toggle observable and checkbox
        state.artifacts.toggles[key] = Observable(false)
        checkbox = Checkbox(topo_nav[row, 2], checked = false, tellheight = false)
        Label(topo_nav[row, 3], "$label ($comp_str)", tellwidth = false, tellheight = false)

        # Wire checkbox -> toggle -> textbox update
        on(checkbox.checked) do checked
            state.artifacts.toggles[key][] = checked
            _update_artifact_textbox!(state)
        end
    end

    # Add extra space between Invert Scale (row 4) and artifact checkboxes (row 5+)
    rowgap!(topo_nav, 4, 45)

    return nothing
end


"""
    _update_artifact_textbox!(state)

Populate component textbox with component numbers from active artifact toggles.
"""
function _update_artifact_textbox!(state)
    isnothing(state.artifacts.textbox) && return

    comps = Int[]
    categories = _get_artifact_categories(state.artifacts.data)
    for (key, _, source) in categories
        state.artifacts.toggles[key][] && append!(comps, source)
    end

    state.artifacts.textbox.displayed_string[] = join(string.(sort(unique(comps))), ",")
    return nothing
end


"""Register keyboard shortcuts for scrolling, zooming, and paging."""
function _setup_keyboard_interactions!(fig, state)

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            _show_plot_help(:ica)
        elseif event.action in (Keyboard.press, Keyboard.repeat)
            if event.key == Keyboard.left || event.key == Keyboard.right

                # Handle x-axis scrolling
                current_range = state.xrange[]
                data_length = size(state.component_data, 2)
                window_size = state.window_size
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - window_size)
                elseif event.key == Keyboard.right
                    new_start = min(data_length - window_size + 1, first(current_range) + window_size)
                end

                # Ensure we don't go beyond data bounds
                new_start = max(1, min(new_start, data_length - window_size + 1))
                state.xrange[] = new_start:(new_start+window_size-1)

                # Update x-axis limits for all axes
                # Ensure the indices are within bounds
                first_idx = clamp(first(state.xrange[]), 1, length(state.dat.data.time))
                last_idx = clamp(last(state.xrange[]), 1, length(state.dat.data.time))
                new_xlims = (state.dat.data.time[first_idx], state.dat.data.time[last_idx])
                xlims!.(state.axs, Ref(new_xlims))

            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed = (Keyboard.left_shift in events(fig).keyboardstate) || (Keyboard.right_shift in events(fig).keyboardstate)

                if !shift_pressed
                    zoom_step = get(state.plot_kwargs, :zoom_step, 0.2)
                    if event.key == Keyboard.up
                        _ymore!.(state.axs, zoom_step)
                    else
                        _yless!.(state.axs, zoom_step)
                    end
                    state.ylims[] = state.axs[1].yaxis.attributes.limits[]

                    # Force a redraw of the plots with scale inversion
                    # For subsets, update all components in the subset
                    if length(state.components) == size(state.component_data, 1)
                        num_plots = state.n_visible_components
                    else
                        num_plots = length(state.components)  # Update all components in the subset
                    end

                    for i = 1:num_plots
                        comp_idx = _get_component_index(state, i)
                        if comp_idx <= size(state.component_data, 1)
                            component_data = state.component_data[comp_idx, :]
                            if state.invert_scale[]
                                component_data = -component_data
                            end
                            state.lines_obs[i][] = component_data
                        end
                    end
                else
                    # HACK! With shift - adjust ONLY channel scale without changing axis limits
                    if event.key == Keyboard.up && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] * 1.1
                    elseif event.key == Keyboard.down && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] / 1.1
                    end
                end

            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                # Only handle page up/down if we're not using fixed components
                if length(state.components) == size(state.component_data, 1)
                    # Handle component scrolling
                    current_start = state.comp_start[]
                    if event.key == Keyboard.page_up
                        new_start = max(1, current_start - state.n_visible_components)
                    else  # page_down
                        new_start =
                            min(size(state.component_data, 1) - state.n_visible_components + 1, current_start + state.n_visible_components)
                    end

                    if new_start != current_start
                        state.comp_start[] = new_start
                        _update_components!(state)
                    end
                end
            end
        end
    end
end

"""Refresh time-series and topoplot data after component navigation or scale change."""
function _update_components!(state)
    if length(state.components) == size(state.component_data, 1)
        num_plots = state.n_visible_components
    else
        num_plots = length(state.components)  # Update all components in the subset
    end

    # Pre-allocate data vector
    all_data = Vector{Matrix{Float64}}(undef, num_plots)
    data_count = 0

    for i = 1:num_plots
        comp_idx = _get_component_index(state, i)
        if comp_idx <= size(state.component_data, 1)
            data_count += 1
            all_data[data_count], _, _ =
                _prepare_ica_topo_data(state.ica, comp_idx, state.plot_kwargs[:method], state.plot_kwargs[:gridscale])
        end
    end

    # Resize to actual count
    resize!(all_data, data_count)

    # Calculate levels using the simplified function
    levels_result = _calculate_ica_topo_levels(all_data, state.use_global_scale[], state.plot_kwargs[:num_levels])

    for i = 1:num_plots

        comp_idx = _get_component_index(state, i)

        if comp_idx <= size(state.component_data, 1)
            # Update component time series data with possible inversion
            component_data = state.component_data[comp_idx, :]
            if state.invert_scale[]
                component_data = -component_data
            end
            # Check if lines_obs exists for this index before updating
            if i <= length(state.lines_obs)
                state.lines_obs[i][] = component_data
            else
                @minimal_warning "Trying to update non-existent lines_obs at index $i"
            end

            # Clear and redraw topography using the viewer function
            if i <= length(state.topo_axs)
                topo_ax = state.topo_axs[i]
                # Determine which level to use for this component
                if state.use_global_scale[]
                    level_to_use = levels_result  # Use global levels for all components
                else
                    level_to_use = levels_result[i]  # Use local levels for this specific component
                end

                _plot_ica_topo_in_viewer!(
                    topo_ax.parent, # Pass the figure associated with the axis
                    topo_ax,
                    state.ica,
                    comp_idx;
                    use_global_scale = state.use_global_scale[],
                    pre_calculated_levels = level_to_use,
                    state.plot_kwargs...,
                )
                topo_ax.ylabel = @sprintf("IC %d\n(%.1f%%)", comp_idx, state.ica.variance[comp_idx] * 100)
                topo_ax.ylabelvisible = true
            else
                @minimal_warning "Trying to update non-existent topo_axs at index $i"
            end
            # Update y-axis label
            if i <= length(state.axs)
                state.axs[i].ylabel = @sprintf("IC %d", comp_idx)
            end
        else
            # Handle invalid comp_idx: clear plots and labels
            if i <= length(state.lines_obs)
                state.lines_obs[i][] = zeros(Float64, size(state.component_data, 2)) # Clear line
            end
            if i <= length(state.topo_axs)
                empty!(state.topo_axs[i])
                state.topo_axs[i].ylabel = @sprintf("Invalid IC %d", comp_idx)
                state.topo_axs[i].ylabelvisible = true
            end
            if i <= length(state.axs)
                state.axs[i].ylabel = ""
            end
        end
    end
end




"""Interpolate ICA mixing weights onto a 2-D grid for topographic display."""
function _prepare_ica_topo_data(ica::InfoIca, comp_idx::Int, method::Symbol, gridscale::Int)
    supported_methods =
        [:multiquadratic, :inverse_multiquadratic, :gaussian, :inverse_quadratic, :thin_plate, :polyharmonic, :shepard, :nearest]
    if method ∈ supported_methods
        data, x_bounds, y_bounds = _data_interpolation_topo(ica.mixing[:, comp_idx], ica.layout, gridscale, method = method)
        return data, x_bounds, y_bounds
    elseif method == :spherical_spline
        data = _data_interpolation_topo_spherical_spline(ica.mixing[:, comp_idx], ica.layout, gridscale)
        # Calculate circular bounds for spherical spline
        x_coords = ica.layout.data.x2
        y_coords = ica.layout.data.y2
        max_radius = maximum(sqrt.(x_coords .^ 2 .+ y_coords .^ 2))
        margin = max_radius * 0.05
        plot_radius = max_radius + margin
        return data, (-plot_radius, plot_radius), (-plot_radius, plot_radius)
    else
        throw(ArgumentError("Unknown interpolation method: $method. Supported: $supported_methods"))
    end
end

"""Calculate contour levels — global across all components or local per component."""
function _calculate_ica_topo_levels(
    all_data::Vector{Matrix{Float64}},
    use_global_scale::Bool,
    num_levels::Int;
    global_min = nothing,
    global_max = nothing,
)
    if use_global_scale # Find global min/max across all data 
        valid_data = filter(!isnan, vcat(all_data...))
        if !isempty(valid_data)
            actual_min, actual_max = extrema(valid_data)
            # Use provided min/max if given, otherwise use actual min/max
            min_val = isnothing(global_min) ? actual_min : global_min
            max_val = isnothing(global_max) ? actual_max : global_max
            # Calculate global levels using base function
            return _calculate_topo_levels(
                reshape([min_val, max_val], 1, 2);
                use_global_scale = true,
                global_min = min_val,
                global_max = max_val,
                num_levels = num_levels,
            )
        else
            throw(ArgumentError("No valid (non-NaN) data. Cannot calculate topographic levels."))
        end
    else # Calculate local levels for each component
        return [_calculate_topo_levels(data; num_levels = num_levels) for data in all_data]
    end
end

"""Map plot-row index `i` to the actual component number."""
function _get_component_index(state, i::Int)
    if length(state.components) == size(state.component_data, 1)
        return state.comp_start[] + i - 1
    else
        return state.components[i]
    end
end


################################################################################
# Default parameters for ICA quality assessment plots with descriptions
const PLOT_ICA_QUALITY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Statistical thresholds
    :z_threshold => (3.0, "Z-score threshold for identifying components"),
    :min_harmonic_power => (1.5, "Minimum harmonic power threshold for line noise detection"),
    :max_ibi_std_s => (0.2, "Maximum inter-beat interval standard deviation for ECG detection"),
    :min_peak_ratio => (0.7, "Minimum peak ratio for ECG detection"),

    # Plot styling
    :title_fontsize => (16, "Font size for plot titles"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :legend_fontsize => (12, "Font size for legend"),

    # Line styling
    :linewidth => (2, "Line width for plot lines"),
    :line_alpha => (0.8, "Transparency for plot lines"),

    # Marker styling
    :marker_size => (8, "Size of markers"),
    :marker_alpha => (0.7, "Transparency of markers"),

    # Color schemes
    :correlation_color => (:steelblue, "Color for correlation plots"),
    :threshold_color => (:red, "Color for threshold lines"),
    :identified_color => (:green, "Color for identified components"),
    :rejected_color => (:red, "Color for rejected components"),

    # Layout parameters
)
