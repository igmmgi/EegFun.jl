# Add imports at the top of the file
"""
    plot_ica_topoplot_single(fig, position, ica, comp_idx, layout; ...)

Plot a single ICA component topography on a given figure and position.

# Arguments
- `fig::Figure`: The target Makie Figure.
- `position`: The position within the figure's layout (e.g., `fig[1, 1]`).
- `ica::InfoIca`: The ICA result object.
- `comp_idx::Int`: The index of the component to plot.
- `layout::DataFrame`: DataFrame containing channel layout information (needs `x`, `y` or `x2`, `y2`).

# Keyword Arguments
- `colorbar_kwargs::Dict`: Controls the colorbar (e.g., `plot_colorbar=true`, `width=10`).
- `head_kwargs::Dict`: Controls the head outline (e.g., `color=:black`, `linewidth=2`).
- `point_kwargs::Dict`: Controls channel markers (e.g., `plot_points=false`).
- `label_kwargs::Dict`: Controls channel labels (e.g., `plot_labels=false`).
- `topo_kwargs::Dict`: Controls the topography plot (e.g., `colormap=:jet`, `gridscale=300`, `radius=88`, `levels=...`, `nan_color=:transparent`). If `levels` is not provided, local scaling is used.

# Returns
- `ax::Axis`: The Makie Axis containing the plot.
- `co`: The contour plot object.
"""
function plot_ica_topoplot_single(
    fig,
    position,
    ica,
    comp_idx,
    layout;
    colorbar_kwargs = Dict(),
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
)

    # default kwargs and merged kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    merged_head_kwargs = merge(head_default_kwargs, head_kwargs) # Keep full merged dict

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    merged_point_kwargs = merge(point_default_kwargs, point_kwargs) # Keep full merged dict

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0) # Removed duplicate :color
    merged_label_kwargs = merge(label_default_kwargs, label_kwargs) # Keep full merged dict

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300, :radius => 88) # Added radius default
    merged_topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = merged_topo_kwargs[:gridscale]     # Access gridscale
    radius = merged_topo_kwargs[:radius]           # Access radius
    merged_topo_kwargs = filter(p -> p.first ∉ (:gridscale, :radius, :levels), merged_topo_kwargs)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    merged_colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = merged_colorbar_kwargs[:plot_colorbar] # Access flag
    colorbar_default_kwargs = filter(p -> p.first != :plot_colorbar, merged_colorbar_kwargs)

    # Create a consistent GridLayout at this position - always the same structure
    gl = fig[position...] = GridLayout()

    # Create main axis in the first column
    ax = Axis(gl[1, 1], title = @sprintf("IC %d (%.1f%%)", comp_idx, ica.variance[comp_idx] * 100))

    # Extract layout data
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :] # Corrected bracket

    # Create the topo data
    data = data_interpolation_topo(ica.mixing[:, comp_idx], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)

    # Calculate levels for the plot
    # Use get on merged_topo_kwargs to check if user provided levels
    levels = get(merged_topo_kwargs, :levels, _calculate_topo_levels(data))

    # Create the plot using filtered kwargs
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data;
        levels = levels, # Pass levels explicitly
        merged_topo_kwargs..., # Pass filtered styling args
    )

    # Add head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = merged_head_kwargs,
        point_kwargs = merged_point_kwargs,
        label_kwargs = merged_label_kwargs,
    )

    # Hide decorations for cleaner look
    hidedecorations!(ax, grid = false)

    # Add colorbar if requested
    if plot_colorbar # Use accessed flag
        width = get(colorbar_plot_kwargs, :width, 10)
        height = get(colorbar_plot_kwargs, :height, Relative(0.8))
        ticksize = get(colorbar_plot_kwargs, :ticklabelsize, 10)

        Colorbar(gl[1, 2], co; ticklabelsize = ticksize, merged_colorbar_kwargs...)

        # Set column sizes
        colsize!(gl, 1, Relative(0.85))  # Plot column
        colsize!(gl, 2, Relative(0.15))  # Colorbar column
    else
        # If no colorbar, make the plot take the full width
        colsize!(gl, 1, Relative(1.0))
    end

    return ax, coa

end

"""
    plot_ica_topoplot(ica, layout; ...)

Plot multiple ICA component topographies in a grid layout within a new Figure.

# Arguments
- `ica::InfoIca`: The ICA result object.
- `layout::DataFrame`: DataFrame containing channel layout information (needs `x`, `y` or `x2`, `y2`).

# Keyword Arguments
- `comps=nothing`: Vector of component indices to plot. If `nothing`, plots all components.
- `dims=nothing`: Tuple specifying grid dimensions (rows, cols). If `nothing`, calculates best square-ish grid.
- `head_kwargs::Dict`: Controls the head outline for all plots.
- `point_kwargs::Dict`: Controls channel markers for all plots.
- `label_kwargs::Dict`: Controls channel labels for all plots.
- `topo_kwargs::Dict`: Controls the topography plots (e.g., `colormap`, `gridscale`, `radius`, `num_levels`, `nan_color`).
- `colorbar_kwargs::Dict`: Controls colorbars (e.g., `plot_colorbar`, `width`, `colorbar_plot_numbers`). `colorbar_plot_numbers` specifies which plot indices (1-based) should have a visible colorbar.
- `use_global_scale::Bool=false`: If `true`, all topoplots share the same color scale based on the min/max across all plotted components. If `false`, each plot is scaled independently.

# Returns
- `fig::Figure`: The generated Makie Figure containing the grid of topoplots.
"""
function plot_ica_topoplot(
    ica,
    layout;
    comps = nothing,
    dims = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    use_global_scale = false,
)

    # default kwargs and merged kwargs
    head_default_kwargs = Dict()
    merged_head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    merged_point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0)
    merged_label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300, :radius => 88, :num_levels => 20, :nan_color => :transparent)
    merged_topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = merged_topo_kwargs[:gridscale]
    radius = merged_topo_kwargs[:radius]
    num_levels = merged_topo_kwargs[:num_levels]
    merged_topo_kwargs = filter(p -> p.first ∉ (:gridscale, :radius, :num_levels, :levels), merged_topo_kwargs)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30, :colorbar_plot_numbers => [])
    merged_colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = merged_colorbar_kwargs[:plot_colorbar]
    colorbar_plot_numbers = merged_colorbar_kwargs[:colorbar_plot_numbers]
    merged_colorbar_kwargs = filter(p -> p.first ∉ (:plot_colorbar, :colorbar_plot_numbers), merged_colorbar_kwargs)

    # Process inputs (keep propertynames check for now)
    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    if isnothing(comps)
        comps = 1:size(ica.mixing)[2]
    end

    # Create figure with reduced margins
    fig = Figure(
        figure_padding = (0, 60, 0, 0), # Keep padding for right margin
        backgroundcolor = :white,
    )

    # Calculate layout dimensions
    if isnothing(dims)
        dims = best_rect(length(comps))
    end

    # Extract layout data once for all plots
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]

    # Calculate all topo data first
    all_data = []
    for i in eachindex(comps)
        # Use gridscale accessed earlier
        data =
            data_interpolation_topo(ica.mixing[:, comps[i]], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)
        push!(all_data, data)
    end

    # Calculate levels for all plots (global or local)
    global_levels = nothing
    all_local_levels = []

    if use_global_scale
        # Find global min/max once
        all_values = [v for data in all_data for v in data if !isnan(v)]
        global_min, global_max = -1.0, 1.0
        if !isempty(all_values)
            global_min = minimum(all_values)
            global_max = maximum(all_values)
        end
        # Use helper to calculate global levels
        global_levels = _calculate_topo_levels(
            reshape(all_values, 1, :); # Pass dummy matrix
            use_global_scale = true,
            global_min = global_min,
            global_max = global_max,
            num_levels = num_levels,
        )
    else
        # Pre-calculate all local levels using helper
        for data in all_data
            # Use the helper which handles NaN safely
            push!(all_local_levels, _calculate_topo_levels(data; num_levels = num_levels)) # Use accessed num_levels
        end
    end

    # First create all the GridLayouts for consistent sizing
    grids = []
    for i in eachindex(comps)
        row = ceil(Int, i / dims[2])
        col = ((i - 1) % dims[2]) + 1
        gl = fig[row, col] = GridLayout()
        push!(grids, gl)
    end

    # Now add the plots and colorbars
    for i in eachindex(comps)
        gl = grids[i]

        # Create main axis in the first column
        ax = Axis(gl[1, 1], title = @sprintf("IC %d (%.1f%%)", comps[i], ica.variance[comps[i]] * 100))

        # Get pre-calculated data and levels
        data = all_data[i]
        levels_to_use = use_global_scale ? global_levels : all_local_levels[i]

        co = _plot_topo_on_axis!(
            ax,
            fig,
            data,
            layout,
            levels_to_use;
            gridscale = gridscale,
            radius = radius,
            head_kwargs = merged_head_kwargs,
            point_kwargs = merged_point_kwargs,
            label_kwargs = merged_label_kwargs,
            merged_topo_kwargs...,
        )

        # Set up colorbar space for ALL plots for consistent sizing
        if plot_colorbar # Use accessed flag
            # Default values from filtered kwargs
            width = get(merged_colorbar_kwargs, :width, 10)
            height = get(merged_colorbar_kwargs, :height, Relative(0.8))
            ticksize = get(merged_colorbar_kwargs, :ticklabelsize, 10)

            # Create colorbar or placeholder in second column
            if i in colorbar_plot_numbers # Use accessed list
                # Pass filtered styling args
                Colorbar(
                    gl[1, 2],
                    co;
                    width = width,
                    height = height,
                    ticklabelsize = ticksize,
                    merged_colorbar_kwargs...,
                )
            else
                # Pass filtered styling args (though Box might not use them all)
                Box(
                    gl[1, 2];
                    width = width,
                    height = height,
                    color = :transparent,
                    strokewidth = 0,
                    merged_colorbar_kwargs...,
                )
            end

            # Set column sizes after creating the colorbar column
            colsize!(gl, 1, Relative(0.85))  # Plot column
            colsize!(gl, 2, Relative(0.15))  # Colorbar column
        else
            # If no colorbar, make the plot take the full width
            colsize!(gl, 1, Relative(1.0))
        end
    end

    return fig

end

# Version for the component viewer to use
function plot_ica_topoplot(
    fig,
    ax,
    ica,
    comp_idx,
    layout;
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)

    # default kwargs and merged kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    merged_head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    merged_point_kwargs = merge(point_defaults, point_kwargs)

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0) # Removed duplicate :color
    merged_label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300, :radius => 88, :nan_color => :transparent) # Added radius, nan_color defaults
    merged_topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = merged_topo_kwargs[:gridscale]     # Access gridscale
    radius = merged_topo_kwargs[:radius]           # Access radius
    merged_topo_kwargs = filter(p -> p.first ∉ (:gridscale, :radius), merged_topo_kwargs)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30, :colorbar_plot_numbers => [])
    merged_colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)

    # Clear the axis before drawing
    empty!(ax)

    # Extract layout data
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]

    # Create the topo data
    data = data_interpolation_topo(ica.mixing[:, comp_idx], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)

    # Create the plot directly on the existing axis using filtered kwargs
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data;
        merged_topo_kwargs...,
    )

    # Add head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = merged_head_kwargs,
        point_kwargs = merged_point_kwargs,
        label_kwargs = merged_label_kwargs,
    )

    ax.title = string(ica.ica_label[comp_idx]) 

    # Hide decorations for cleaner look
    hidedecorations!(ax, grid = false)

    return co

end


# Create a state structure to hold the visualization state
mutable struct IcaComponentState
    # Data
    dat::ContinuousData
    ica_result::InfoIca
    components::Matrix{Float64}
    total_components::Int

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

    # New field for specific components
    specific_components::Union{Nothing,Vector{Int}}

    #  Topoplot Kwargs 
    topo_gridscale::Int
    topo_radius::Real
    topo_colormap::Symbol
    topo_num_levels::Int
    topo_nan_color::Union{Symbol,Makie.Colorant}
    topo_head_kwargs::Dict
    topo_point_kwargs::Dict
    topo_label_kwargs::Dict

    # Plot elements
    axs::Vector{Axis}
    channel_axs::Vector{Union{Axis,Nothing}}  # Allow for nothing values
    topo_axs::Vector{Axis}
    lines_obs::Vector{Observable{Vector{Float64}}}

    # New field for boolean indicators
    channel_bool_indicators::Dict{Int,Any}

    # Constructor updated to accept and store topo_kwargs
    function IcaComponentState(
        dat,
        ica_result,
        n_visible_components,
        window_size,
        specific_components = nothing;
        topo_kwargs = Dict(),
        head_kwargs = Dict(),
        point_kwargs = Dict(),
        label_kwargs = Dict(),
    ) 
        # Prepare data matrix
        dat_matrix = prepare_ica_data_matrix(dat, ica_result)
        components = ica_result.unmixing * dat_matrix
        total_components = size(components, 1)

        # Create observables
        comp_start = Observable(1)
        use_global_scale = Observable(false)
        invert_scale = Observable(false)

        # Find index closest to time 0 to center the initial view
        time_zero_idx = findmin(abs.(dat.data.time))[2]
        half_window = div(window_size, 2)
        start_idx = max(1, time_zero_idx - half_window)
        end_idx = min(size(components, 2), start_idx + window_size - 1)

        # Adjust start_idx if end_idx reached the boundary
        if end_idx == size(components, 2)
            start_idx = max(1, end_idx - window_size + 1)
        end

        xrange = Observable(start_idx:end_idx)

        # Set initial range based on specific components if provided
        comps_to_use = if !isnothing(specific_components)
            specific_components
        else
            1:min(n_visible_components, total_components) # Ensure range is valid
        end

        # Calculate initial y-range based on components we'll show
        initial_range_data = if !isempty(comps_to_use) && all(idx -> idx <= total_components, comps_to_use)
            components[comps_to_use, start_idx:end_idx]
        else
            zeros(0, length(start_idx:end_idx)) # Empty if no valid components
        end
        initial_range = if !isempty(initial_range_data)
            maximum(abs.(extrema(initial_range_data)))
        else
            1.0 # Default range if no data
        end
        ylims = Observable((-initial_range, initial_range))
        channel_data = Observable(zeros(size(dat.data, 1)))
        show_channel = Observable(false)
        channel_yscale = Observable(1.0)

        # --- Process and store Topo Kwargs ---
        topo_defaults =
            Dict(:gridscale => 300, :radius => 88, :colormap => :jet, :num_levels => 20, :nan_color => :transparent)
        processed_topo_kwargs = merge(topo_defaults, topo_kwargs)

        head_defaults = Dict(:color => :black, :linewidth => 2)
        processed_head_kwargs = merge(head_defaults, head_kwargs)

        point_defaults = Dict(:plot_points => false)
        processed_point_kwargs = merge(point_defaults, point_kwargs)

        label_defaults = Dict(:plot_labels => false)
        processed_label_kwargs = merge(label_defaults, label_kwargs)

        # Store processed values
        topo_gridscale = processed_topo_kwargs[:gridscale]
        topo_radius = processed_topo_kwargs[:radius]
        topo_colormap = processed_topo_kwargs[:colormap]
        topo_num_levels = processed_topo_kwargs[:num_levels]
        topo_nan_color = processed_topo_kwargs[:nan_color]
        # --- End Process and store ---

        # Initialize empty plot element arrays
        axs = Vector{Axis}()
        channel_axs = Vector{Union{Axis,Nothing}}()
        topo_axs = Vector{Axis}()
        lines_obs = Vector{Observable{Vector{Float64}}}()
        channel_bool_indicators = Dict{Int,Any}()

        new(
            dat,
            ica_result,
            components,
            total_components,
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
            specific_components,
            # Pass stored kwargs
            topo_gridscale,
            topo_radius,
            topo_colormap,
            topo_num_levels,
            topo_nan_color,
            processed_head_kwargs,
            processed_point_kwargs,
            processed_label_kwargs,
            # Plot elements
            axs,
            channel_axs,
            topo_axs,
            lines_obs,
            channel_bool_indicators,
        )
    end
end

"""
    plot_ica_component_activation(dat::ContinuousData, ica_result::InfoIca; ...)

Create an interactive visualization of ICA components with topographic maps and time series plots.

Allows scrolling through components and time, adjusting scales, and overlaying raw channel data.

# Arguments
- `dat::ContinuousData`: Continuous EEG data (must contain a `layout`).
- `ica_result::InfoIca`: ICA result object (must match `dat`).

# Keyword Arguments
- `n_visible_components::Int=10`: Number of components visible vertically at once.
- `window_size::Int=2000`: Initial time window size in samples.
- `topo_kwargs::Dict=Dict()`: Keyword arguments passed down for topography plots (see `_plot_topo_on_axis!`).
- `head_kwargs::Dict=Dict()`: Keyword arguments passed down for head outlines.
- `point_kwargs::Dict=Dict()`: Keyword arguments passed down for channel markers.
- `label_kwargs::Dict=Dict()`: Keyword arguments passed down for channel labels.
- `specific_components=nothing`: An optional vector of specific component indices to display initially. If provided, `n_visible_components` is ignored, and only these components are shown without scrolling capability (PageUp/PageDown disabled).

# Returns
- `fig::Figure`: The Makie Figure object containing the interactive plot.
"""
function plot_ica_component_activation(
    dat::ContinuousData,
    ica_result::InfoIca;
    n_visible_components::Int = 10,
    window_size::Int = 2000,
    topo_kwargs = Dict(),
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    specific_components = nothing,
)
    # Create state, using specific components if provided
    if !isnothing(specific_components)
        n_visible_components = length(specific_components)
    end

    # Pass kwargs to constructor
    state = IcaComponentState(
        dat,
        ica_result,
        n_visible_components,
        window_size,
        specific_components;
        topo_kwargs = topo_kwargs,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    # Create figure with padding on the right for margin
    fig = Figure(
        figure_padding = (0, 60, 0, 0), # Keep right padding
        backgroundcolor = :white,
    )

    # Setup plots (no need to pass topo_kwargs here anymore)
    create_component_plots!(fig, state) # Removed topo_kwargs argument

    # Add controls
    add_navigation_controls!(fig, state)

    # Add navigation sliders
    add_navigation_sliders!(fig, state)

    # Add channel menu
    add_channel_menu!(fig, state)

    # Add keyboard interactions
    setup_keyboard_interactions!(fig, state)

    # --- Layout Adjustments ---
    colsize!(fig.layout, 1, Relative(0.15))
    colsize!(fig.layout, 2, Relative(0.85))
    rowgap!(fig.layout, state.n_visible_components, 30)
    # --- End Layout Adjustments ---

    display(fig)
    return fig
end

"""
    prepare_ica_data_matrix(dat::ContinuousData, ica_result::InfoIca)

Selects, centers, scales, and transposes data for ICA unmixing.

# Arguments
- `dat::ContinuousData`: Input EEG data.
- `ica_result::InfoIca`: Corresponding ICA result containing labels and scaling factors.

# Returns
- `Matrix{Float64}`: Data matrix ready for `ica_result.unmixing * dat_matrix`. (channels x samples)
"""
function prepare_ica_data_matrix(dat::ContinuousData, ica_result::InfoIca)
    dat_matrix = permutedims(Matrix(dat.data[!, ica_result.data_label]))
    dat_matrix .-= mean(dat_matrix, dims = 2) # TODO: check if this is correct
    dat_matrix ./= ica_result.scale
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
    local_min, local_max = NaN, NaN

    # Find local min/max, ignoring NaNs
    valid_values = [v for v in data if !isnan(v)]
    if !isempty(valid_values)
        local_min = minimum(valid_values)
        local_max = maximum(valid_values)
    else
        # Default if data is all NaN
        local_min = -1.0
        local_max = 1.0
    end

    # Determine final min/max for levels
    final_min, final_max = if use_global_scale && !isnothing(global_min) && !isnothing(global_max)
        global_min, global_max
    else
        local_min, local_max
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
- `radius::Real=88`: Radius used for coordinate scaling.
- `colormap=:jet`: Colormap for the contour plot.
- `nan_color=:transparent`: Color for NaN values in the data.
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
    layout::DataFrame,
    levels;
    gridscale = 300,
    radius = 88,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(), # These are NOW the merged versions from the caller
    # Remaining kwargs are for contourf!
    kwargs...,
) # Use kwargs... to capture filtered contourf args

    # Plot contour using the passed kwargs...
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data;
        levels = levels,
        kwargs..., # Pass the filtered contourf args directly
    )

    # Add head shape using plot_layout_2d!
    # Pass the ALREADY MERGED kwargs it received directly to plot_layout_2d!
    plot_layout_2d!(
        fig,
        ax,
        layout;
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    # Hide decorations
    hidexdecorations!(ax, grid = false)
    hideydecorations!(ax, grid = false)

    return co
end

# Create a simpler override specifically for the component viewer
# Now uses the internal helpers
function plot_topoplot_in_viewer!(
    fig,
    topo_ax,
    ica_result,
    comp_idx,
    layout;
    use_global_scale = false,
    global_min = nothing,
    global_max = nothing,
    gridscale = 300,
    radius = 88,
    colormap = :jet,
    num_levels = 20,
    nan_color = :transparent,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
)
    # Clear the axis
    empty!(topo_ax)

    # Extract layout data
    tmp_layout = layout[(in.(layout.label, Ref(ica_result.data_label))), :] # Use 'in.'

    # Create the topo data
    data = data_interpolation_topo(
        ica_result.mixing[:, comp_idx],
        permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
        gridscale,
    )

    # Calculate levels using helper
    levels = _calculate_topo_levels(
        data;
        use_global_scale = use_global_scale,
        global_min = global_min,
        global_max = global_max,
        num_levels = num_levels,
    )

    # Plot using helper
    co = _plot_topo_on_axis!(
        topo_ax,
        fig,
        data,
        layout,
        levels;
        gridscale = gridscale,
        radius = radius,
        colormap = colormap,
        nan_color = nan_color,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    # No title setting here, handled by calling function (create_component_plots!)

    return co
end

# Update create_component_plots! - Remove topo_kwargs argument
function create_component_plots!(fig, state) # Removed topo_kwargs argument
    # --- No local merging needed anymore ---

    # Create axes for the time series plots
    for i = 1:state.n_visible_components
        # Create topo plot axis first (now on the left)
        topo_ax = Axis(fig[i, 1])
        push!(state.topo_axs, topo_ax)

        # Get the actual component number
        comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
            state.specific_components[i]
        else
            state.comp_start[] + i - 1
        end

        # Time series axis creation (now on the right)
        ax = Axis(
            fig[i, 2],
            ylabel = @sprintf("IC %d", comp_idx),
            yaxisposition = :left,
            yticklabelsvisible = false,
            yticksvisible = true,
            xticklabelsvisible = (i == state.n_visible_components),
            xticksvisible = (i == state.n_visible_components),
            xgridvisible = false,
            ygridvisible = false,
            xminorgridvisible = false,
            yminorgridvisible = false,
            ylabelpadding = 0.0,
            yticklabelpad = 0.0,
            yticklabelspace = 0.0,
        )
        push!(state.axs, ax)

        # Always create channel overlay axis (keep spines hidden)
        ax_channel = Axis(
            fig[i, 2],
            yticklabelsvisible = false,
            yticksvisible = false,
            yaxisposition = :right,
            xaxisposition = :top,
            xticklabelsvisible = false,
            xticksvisible = false,
            xgridvisible = false,
            ygridvisible = false,
            xminorgridvisible = false,
            yminorgridvisible = false,
            bottomspinevisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            leftspinevisible = false,
        )
        push!(state.channel_axs, ax_channel)

        # Link axes
        linkyaxes!(ax, ax_channel)
        linkxaxes!(ax, ax_channel)

        # Channel overlay plot
        lines!(
            ax_channel,
            @lift(state.dat.data.time[$(state.xrange)]),
            @lift(
                $(state.show_channel) ? $(state.channel_data)[$(state.xrange)] .* $(state.channel_yscale) :
                zeros(Float64, length($(state.xrange)))
            ),
            color = :grey,
        )

        # Set initial x-axis limits
        xlims!(ax_channel, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))

        # Observable creation for component plot
        # Handle potential invalid comp_idx robustly
        comp_data = if comp_idx <= state.total_components
            state.components[comp_idx, :]
        else
            zeros(Float64, size(state.components, 2)) # Placeholder data
        end
        lines_obs = Observable(comp_data)
        push!(state.lines_obs, lines_obs)

        # Component line plot
        lines!(ax, @lift(state.dat.data.time[$(state.xrange)]), @lift($(lines_obs)[$(state.xrange)]), color = :black)

        # Set initial x-axis limits for component plot
        xlims!(ax, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))

        # Create the topo plot using the dedicated viewer function
        if comp_idx <= state.total_components
            # --- Pass parameters from state ---
            plot_topoplot_in_viewer!(
                fig,
                topo_ax,
                state.ica_result,
                comp_idx,
                state.dat.layout;
                use_global_scale = state.use_global_scale[],
                gridscale = state.topo_gridscale,
                radius = state.topo_radius,
                colormap = state.topo_colormap,
                num_levels = state.topo_num_levels,
                nan_color = state.topo_nan_color,
                head_kwargs = state.topo_head_kwargs,
                point_kwargs = state.topo_point_kwargs,
                label_kwargs = state.topo_label_kwargs,
            )
            # --- End Pass parameters ---
            # Set title here
            topo_ax.title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
        else
            empty!(topo_ax) # Clear axis if comp_idx is invalid
            topo_ax.title = @sprintf("Invalid IC %d", comp_idx)
        end
    end
end

# Update add_navigation_controls! to include the global scale checkbox
function add_navigation_controls!(fig, state)
    # Add navigation buttons below topo plots in column 1
    # Remove padding here
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1], tellheight = false)

    # --- Add an empty column 1 for spacing ---
    colsize!(topo_nav, 1, 40) # Set width of the empty first column

    # --- Shift all widgets one column to the right ---
    # Navigation buttons now in row 1, columns 2 & 3
    prev_topo = Button(topo_nav[1, 2], label = "◄ Previous", tellheight = false)
    next_topo = Button(topo_nav[1, 3], label = "Next ►", tellheight = false)

    # Component selection now in row 2, columns 2, 3, 4
    text_label = Label(topo_nav[2, 2], "Components:", tellheight = false, halign = :right)
    text_input = Textbox(topo_nav[2, 3], placeholder = "e.g. 1,3-5,8", tellheight = false)
    apply_button = Button(topo_nav[2, 4], label = "Apply", tellheight = false)

    # Global scale checkbox now in row 3, columns 2 & 3
    global_scale_check = Checkbox(topo_nav[3, 2], checked = state.use_global_scale[], tellheight = false)
    Label(topo_nav[3, 3], "Use Global Scale", tellwidth = false, tellheight = false)

    # Invert scale checkbox now in row 4, columns 2 & 3
    invert_scale_check = Checkbox(topo_nav[4, 2], checked = state.invert_scale[], tellheight = false)
    Label(topo_nav[4, 3], "Invert Scale", tellwidth = false, tellheight = false)

    # --- Gap settings (adjust column indices) ---
    # Add column gaps for better spacing
    colgap!(topo_nav, 2, 10) # Gap after the *new* column 2 (was 1)
    colgap!(topo_nav, 3, 5)  # Gap after the *new* column 3 (was 2)

    # Add row gaps for vertical spacing (unchanged values, indices okay)
    rowgap!(topo_nav, 1, 35)
    rowgap!(topo_nav, 2, 45)
    rowgap!(topo_nav, 3, 35)
    # --- End Gap settings ---


    # Connect checkboxes to state
    on(global_scale_check.checked) do checked
        state.use_global_scale[] = checked
        update_components!(state)
    end

    on(invert_scale_check.checked) do checked
        state.invert_scale[] = checked
        update_components!(state)
    end

    # Connect navigation buttons
    on(prev_topo.clicks) do _
        new_start = max(1, state.comp_start[] - state.n_visible_components)
        state.comp_start[] = new_start
        update_components!(state)
    end

    on(next_topo.clicks) do _
        new_start = min(
            state.total_components - state.n_visible_components + 1,
            state.comp_start[] + state.n_visible_components,
        )
        state.comp_start[] = new_start
        update_components!(state)
    end

    # Connect apply button
    on(apply_button.clicks) do _
        text_value = text_input.displayed_string[]
        if !isempty(text_value)
            comps = parse_component_input(text_value, state.total_components)
            if !isempty(comps)
                println("Creating new plot with components: $comps")
                current_channel = state.channel_data[]
                show_channel = state.show_channel[]
                use_global = state.use_global_scale[]
                invert = state.invert_scale[]
                new_fig = plot_ica_component_activation(
                    state.dat,
                    state.ica_result,
                    specific_components = comps,
                    n_visible_components = length(comps),
                    window_size = state.window_size,
                )
            end
        end
    end
end

# Function to parse component input text into a list of component indices
"""
    parse_component_input(text::String, total_components::Int)

Parses a comma-separated string potentially containing ranges (e.g., "1,3-5,8")
into a sorted, unique vector of valid component indices.

# Arguments
- `text::String`: The input string.
- `total_components::Int`: The maximum valid component index.

# Returns
- `Vector{Int}`: Sorted, unique vector of valid component indices found in the text.
                 Returns an empty vector if input is empty or invalid.
"""
function parse_component_input(text::String, total_components::Int)
    components = Int[]
    if isempty(text)
        return components
    end

    try
        # Split by comma
        parts = strip.(split(text, ','))
        for part in parts
            if occursin('-', part)
                # Handle ranges like "1-5"
                range_parts = strip.(split(part, '-'))
                if length(range_parts) == 2
                    start_num = parse(Int, range_parts[1])
                    end_num = parse(Int, range_parts[2])
                    if 1 <= start_num <= end_num <= total_components
                        append!(components, start_num:end_num)
                    end
                end
            else
                # Handle single numbers
                num = parse(Int, part)
                if 1 <= num <= total_components
                    push!(components, num)
                end
            end
        end
    catch e
        # Silently handle parsing errors
    end

    # Remove duplicates and sort
    unique!(sort!(components))
    return components
end

# Update add_navigation_sliders! to match new layout
function add_navigation_sliders!(fig, state)
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
        update_view_range!(state, Int(round(x)))
    end
end

# Helper function to update view range (similar to plot_databrowser style)
function update_view_range!(state, start_pos)
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

    # If we have a boolean channel selected, update its indicators
    if state.show_channel[] && !isempty(state.channel_bool_indicators)
        # Clear existing indicators
        for (i, indicator) in state.channel_bool_indicators
            if !isnothing(indicator)
                delete!(state.channel_axs[i], indicator)
            end
        end
        empty!(state.channel_bool_indicators)

        # Find the current channel
        current_channel = nothing
        for col in names(state.dat.data)
            if state.dat.data[!, col] == state.channel_data[]
                current_channel = Symbol(col)
                break
            end
        end

        # If we found a boolean channel, redraw its indicators
        if !isnothing(current_channel) && eltype(state.dat.data[!, current_channel]) == Bool
            add_boolean_indicators!(state, current_channel)
        end
    end
end

# Update add_channel_menu! to match new layout
function add_channel_menu!(fig, state)
    # Create a menu layout in column 2
    menu_row = state.n_visible_components + 2  # Keep menu in its own row

    # Create a simple grid layout
    menu_layout = GridLayout(fig[menu_row, 2], tellheight = false)

    # Create a simple label and menu
    Label(menu_layout[1, 1], "Additional Channel:", fontsize = 18, tellheight = false, width = 150)  # Fixed width for label
    channel_menu =
        Menu(menu_layout[1, 2], options = ["None"; names(state.dat.data)], default = "None", tellheight = false)

    on(channel_menu.selection) do selected
        update_channel_selection!(state, selected)
    end

    return menu_layout
end

# Helper function to update channel selection (matches databrowser pattern)
function update_channel_selection!(state, selected)
    # Clear previous channel visualizations from all axes
    for i = 1:state.n_visible_components
        if i <= length(state.channel_axs) &&
           haskey(state.channel_bool_indicators, i) &&
           !isnothing(state.channel_bool_indicators[i])
            delete!(state.channel_axs[i], state.channel_bool_indicators[i])
            state.channel_bool_indicators[i] = nothing
        end
    end
    empty!(state.channel_bool_indicators)

    if selected == "None"
        state.show_channel[] = false
        state.channel_data[] = zeros(Float64, size(state.dat.data, 1))
    else
        selected_sym = Symbol(selected)
        state.channel_data[] = state.dat.data[!, selected_sym]

        # Only show overlay plot for non-Boolean channels
        if eltype(state.dat.data[!, selected_sym]) == Bool
            state.show_channel[] = false
            add_boolean_indicators!(state, selected_sym)
        else
            state.show_channel[] = true
        end
    end
end

# Helper function to add boolean indicators (matching databrowser pattern)
function add_boolean_indicators!(state, channel_sym)
    # For each component axis, create a vertical line at each true position
    for i = 1:state.n_visible_components
        if i <= length(state.channel_axs)
            ax_channel = state.channel_axs[i]

            # Find all time points where the boolean is true
            true_indices = findall(state.dat.data[!, channel_sym])

            if !isempty(true_indices)
                # Get the time values for the true positions
                true_times = state.dat.data.time[true_indices]

                # Create vertical lines at each true position
                # Only create lines within the current view range
                current_range = state.xrange[]
                visible_times = true_times[true_times.>=state.dat.data.time[first(
                    current_range,
                )].&&true_times.<=state.dat.data.time[last(current_range)]]

                if !isempty(visible_times)
                    lines = vlines!(ax_channel, visible_times, color = :red, linewidth = 1)

                    # Store the reference to the lines
                    state.channel_bool_indicators[i] = lines
                end
            end
        end
    end
end

# Setup keyboard interactions
function setup_keyboard_interactions!(fig, state)
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            if event.key == Keyboard.left || event.key == Keyboard.right
                # Handle x-axis scrolling
                current_range = state.xrange[]
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - state.window_size)
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                else  # right
                    new_start =
                        min(size(state.components, 2) - state.window_size + 1, first(current_range) + state.window_size)
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                end

                # Update x-axis limits for all axes
                # Ensure the indices are within bounds
                first_idx = clamp(first(state.xrange[]), 1, length(state.dat.data.time))
                last_idx = clamp(last(state.xrange[]), 1, length(state.dat.data.time))
                new_xlims = (state.dat.data.time[first_idx], state.dat.data.time[last_idx])
                for ax in state.axs
                    xlims!(ax, new_xlims)
                end
                for ax in state.channel_axs
                    if !isnothing(ax)
                        xlims!(ax, new_xlims)
                    end
                end

            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed =
                    (Keyboard.left_shift in events(fig).keyboardstate) ||
                    (Keyboard.right_shift in events(fig).keyboardstate)

                if !shift_pressed
                    # Handle y-axis scaling
                    current_range = state.ylims[][2]  # Just take the positive limit
                    if event.key == Keyboard.up
                        # Zoom in - decrease range by 20%
                        new_range = current_range * 0.8
                    else  # down
                        # Zoom out - increase range by 20%
                        new_range = current_range * 1.2
                    end

                    # Keep centered on zero
                    state.ylims[] = (-new_range, new_range)

                    # Update y-axis limits for all axes
                    for ax in state.axs
                        ylims!(ax, state.ylims[])
                    end
                    # Also update channel axes to maintain alignment
                    for ax in state.channel_axs
                        if !isnothing(ax)
                            ylims!(ax, state.ylims[])
                        end
                    end

                    # Force a redraw of the plots with scale inversion
                    for i = 1:state.n_visible_components
                        # Get the correct component index based on whether we're using specific components
                        comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
                            state.specific_components[i]
                        else
                            state.comp_start[] + i - 1
                        end

                        if comp_idx <= state.total_components
                            component_data = state.components[comp_idx, :]
                            if state.invert_scale[]
                                component_data = -component_data
                            end
                            state.lines_obs[i][] = component_data
                        end
                    end
                else
                    # With shift - adjust ONLY channel scale without changing axis limits
                    if event.key == Keyboard.up && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] * 1.1
                    elseif event.key == Keyboard.down && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] / 1.1
                    end
                    # We don't modify any axis limits here, only the scaling factor
                    # This will affect how the channel data is plotted through the Observable
                end

            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                # Only handle page up/down if we're not using specific components
                if isnothing(state.specific_components)
                    # Handle component scrolling
                    current_start = state.comp_start[]
                    if event.key == Keyboard.page_up
                        new_start = max(1, current_start - state.n_visible_components)
                    else  # page_down
                        new_start = min(
                            state.total_components - state.n_visible_components + 1,
                            current_start + state.n_visible_components,
                        )
                    end

                    if new_start != current_start
                        state.comp_start[] = new_start
                        update_components!(state)
                    end
                end
            end
        end
    end
end

# Update component data when navigating
function update_components!(state)
    # Calculate global min/max if using global scale
    global_min, global_max = nothing, nothing
    if state.use_global_scale[]
        all_values = Float64[]
        # --- Use gridscale from state ---
        gridscale = state.topo_gridscale
        # --- End Use gridscale ---
        for i = 1:state.n_visible_components
            comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
                state.specific_components[i]
            else
                state.comp_start[] + i - 1
            end

            if comp_idx <= state.total_components
                tmp_layout = state.dat.layout[(in.(state.dat.layout.label, Ref(state.ica_result.data_label))), :]
                data = data_interpolation_topo(
                    state.ica_result.mixing[:, comp_idx],
                    permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
                    gridscale, # Use state value
                )
                append!(all_values, [v for v in data if !isnan(v)])
            end
        end

        if !isempty(all_values)
            global_min = minimum(all_values)
            global_max = maximum(all_values)
        end
    end

    for i = 1:state.n_visible_components
        comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
            state.specific_components[i]
        else
            state.comp_start[] + i - 1
        end

        if comp_idx <= state.total_components
            # Update component time series data with possible inversion
            component_data = state.components[comp_idx, :]
            if state.invert_scale[]
                component_data = -component_data
            end
            # Check if lines_obs exists for this index before updating
            if i <= length(state.lines_obs)
                state.lines_obs[i][] = component_data
            else
                println("Warning: Trying to update non-existent lines_obs at index $i")
            end

            # Clear and redraw topography using the viewer function
            if i <= length(state.topo_axs)
                topo_ax = state.topo_axs[i]
                # --- Pass parameters from state ---
                plot_topoplot_in_viewer!(
                    topo_ax.parent, # Pass the figure associated with the axis
                    topo_ax,
                    state.ica_result,
                    comp_idx,
                    state.dat.layout; # Semicolon for kwargs
                    use_global_scale = state.use_global_scale[],
                    global_min = global_min,
                    global_max = global_max,
                    # Pass stored kwargs from state
                    gridscale = state.topo_gridscale,
                    radius = state.topo_radius,
                    colormap = state.topo_colormap,
                    num_levels = state.topo_num_levels,
                    nan_color = state.topo_nan_color,
                    head_kwargs = state.topo_head_kwargs,
                    point_kwargs = state.topo_point_kwargs,
                    label_kwargs = state.topo_label_kwargs,
                )
                # --- End Pass parameters ---
                # Update title
                topo_ax.title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
            else
                println("Warning: Trying to update non-existent topo_axs at index $i")
            end
            # Update y-axis label
            if i <= length(state.axs)
                state.axs[i].ylabel = @sprintf("IC %d", comp_idx)
            end
        else
            # Handle invalid comp_idx: clear plots and labels
            if i <= length(state.lines_obs)
                state.lines_obs[i][] = zeros(Float64, size(state.components, 2)) # Clear line
            end
            if i <= length(state.topo_axs)
                empty!(state.topo_axs[i])
                state.topo_axs[i].title = @sprintf("Invalid IC %d", comp_idx)
            end
            if i <= length(state.axs)
                state.axs[i].ylabel = ""
            end
        end
    end
end

"""
    identify_eye_components(ica_result::InfoIca, dat::ContinuousData;
                          v_eog_channel::Symbol=:vEOG,
                          h_eog_channel::Symbol=:hEOG,
                          z_threshold::Float64=3.0,
                          exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value])

Identify ICA components potentially related to eye movements based on z-scored correlation.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data containing EOG channels.

# Keyword Arguments
- `v_eog_channel::Symbol`: Name of the vertical EOG channel (default: :vEOG).
- `h_eog_channel::Symbol`: Name of the horizontal EOG channel (default: :hEOG).
- `z_threshold::Float64`: Absolute Z-score threshold for identification (default: 3.0).
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.

# Returns
- `Dict{Symbol, Vector{Int}}`: Dictionary containing:
  - `:vertical_eye`: Vector of indices identified for vertical eye movements.
  - `:horizontal_eye`: Vector of indices identified for horizontal eye movements.
- `DataFrame`: DataFrame containing detailed correlation metrics per component.
"""
function identify_eye_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    v_eog_channel::Symbol = :vEOG,
    h_eog_channel::Symbol = :hEOG,
    z_threshold::Float64 = 3.0,
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value]
)
    # --- Input Validation ---
    if !(v_eog_channel in propertynames(dat.data))
        error("Vertical EOG channel $v_eog_channel not found in data")
    end
    if !(h_eog_channel in propertynames(dat.data))
        error("Horizontal EOG channel $h_eog_channel not found in data")
    end

    # --- Data Preparation ---
    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vertical_eye => Int[], :horizontal_eye => Int[]), DataFrame()
    end

    # Get EOG signals for valid samples only
    v_eog = dat.data[samples_to_use, v_eog_channel]
    h_eog = dat.data[samples_to_use, h_eog_channel]

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate components for valid samples
    components = ica_result.unmixing * dat_matrix
    n_components = size(components, 1)

    # Function to calculate correlations for all components
    function calculate_correlations(eog_signal)
        corrs = zeros(n_components)
        for comp_idx = 1:n_components 
             corrs[comp_idx] = abs(cor(components[comp_idx, :], eog_signal))
        end
        return corrs
    end

    identified_vEOG = Int[]
    identified_hEOG = Int[]
    v_eog_corr_z = Float64[]
    h_eog_corr_z = Float64[]
    v_eog_corrs = Float64[]
    h_eog_corrs = Float64[]

    # --- Vertical EOG --- 
    v_eog_corrs = calculate_correlations(v_eog)
    v_eog_corr_z = StatsBase.zscore(v_eog_corrs)
    identified_vEOG = findall(abs.(v_eog_corr_z) .> z_threshold)

    # --- Horizontal EOG --- 
    h_eog_corrs = calculate_correlations(h_eog)
    h_eog_corr_z = StatsBase.zscore(h_eog_corrs)
    identified_hEOG = findall(abs.(h_eog_corr_z) .> z_threshold)

    # --- Construct Results --- 
    sort!(identified_vEOG)
    sort!(identified_hEOG)

    result_dict = Dict{Symbol,Vector{Int}}( 
        :vertical_eye => identified_vEOG,
        :horizontal_eye => identified_hEOG,
    )
    
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :v_eog_channel_corr => v_eog_corrs,
        :v_eog_channel_zscore => v_eog_corr_z,
        :h_eog_channel_corr => h_eog_corrs,
        :h_eog_channel_zscore => h_eog_corr_z
    )

    return result_dict, metrics_df 
end

"""
    plot_eye_component_features(identified_comps::Dict, metrics_df::DataFrame; z_threshold::Float64=3.0)

Plot z-scores of EOG correlations from the metrics DataFrame and highlight identified components.

Uses the results from `identify_eye_components`.

# Arguments
- `identified_comps::Dict`: Dictionary returned by `identify_eye_components` (containing `:vertical_eye`, `:horizontal_eye`).
- `metrics_df::DataFrame`: DataFrame returned by `identify_eye_components`. Expected to have columns like `:vEOG_zscore`, `:hEOG_zscore` (based on the `v_eog_channel` and `h_eog_channel` arguments) and `:Component`.
- `v_eog_channel::Symbol`: Symbol of the vertical EOG channel used. Defaults to `:vEOG`. (Used to find column `Symbol("\$(v_eog_channel)_zscore")` in `metrics_df`).
- `h_eog_channel::Symbol`: Symbol of the horizontal EOG channel used. Defaults to `:hEOG`. (Used to find column `Symbol("\$(h_eog_channel)_zscore")` in `metrics_df`).
- `z_threshold::Float64`: Z-score threshold to draw lines on the plot (default: 3.0).

# Returns
- `fig::Figure`: The Makie Figure containing the z-score plots.
"""
function plot_eye_component_features(identified_comps::Dict, metrics_df::DataFrame; z_threshold::Float64=3.0)

    # Extract data from inputs
    v_eog_corr_z = metrics_df.v_eog_channel_zscore
    h_eog_corr_z = metrics_df.h_eog_channel_zscore
    final_vEOG = identified_comps[:vertical_eye] 
    final_hEOG = identified_comps[:horizontal_eye]
    
    # Check if data is empty
    if isempty(metrics_df) || isempty(v_eog_corr_z) || isempty(h_eog_corr_z)
        println("Warning: Could not plot eye component features, input DataFrame or z-scores are empty.")
        return Figure() # Return empty figure
    end
    
    n_components = nrow(metrics_df)

    fig = Figure(size=(800, 400)) 

    # Plot Vertical EOG Correlation Z-Scores
    ax_v = Axis(
        fig[1, 1],
        xlabel = "Component Number",
        ylabel = "Z-Score",
        title = "Vertical EOG Correlation Z-Scores"
    )
    # Use component indices from DataFrame for x-axis
    scatter!(ax_v, metrics_df.Component, v_eog_corr_z, color = :gray, markersize = 5)
    hlines!(ax_v, [z_threshold, -z_threshold], color = :gray, linestyle = :dash)
    # Highlight all identified components
    if !isempty(final_vEOG)
         # Get z-scores only for the identified components
        v_z_scores_highlight = metrics_df[in.(metrics_df.Component, Ref(final_vEOG)), :v_eog_channel_zscore]
        scatter!(ax_v, final_vEOG, v_z_scores_highlight, color = :blue, markersize = 8)
        for comp_idx in final_vEOG # Annotate each identified component
            text!(ax_v, comp_idx, metrics_df[comp_idx, :v_eog_channel_zscore], text=string(comp_idx), color=:blue, align=(:center,:bottom), fontsize=10)
        end
    end

    # Plot Horizontal EOG Correlation Z-Scores
    ax_h = Axis(
         fig[1, 2],
        xlabel = "Component Number",
        ylabel = "Z-Score",
        title = "Horizontal EOG Correlation Z-Scores"
    )
     # Use component indices from DataFrame for x-axis
    scatter!(ax_h, metrics_df.Component, h_eog_corr_z, color = :gray, markersize = 5)
    hlines!(ax_h, [z_threshold, -z_threshold], color = :gray, linestyle = :dash)
     # Highlight all identified components
    if !isempty(final_hEOG)
        # Get z-scores only for the identified components
        h_z_scores_highlight = metrics_df[in.(metrics_df.Component, Ref(final_hEOG)), :h_eog_channel_zscore]
        scatter!(ax_h, final_hEOG, h_z_scores_highlight, color = :blue, markersize = 8)
         for comp_idx in final_hEOG # Annotate each identified component
            text!(ax_h, comp_idx, metrics_df[comp_idx, :h_eog_channel_zscore], text=string(comp_idx), color=:blue, align=(:center,:bottom), fontsize=10)
        end
    end

    return fig
end

# --- EKG Artifact Identification ---

# Helper function (assuming similar logic needed elsewhere, or define inline)
function _get_samples_to_use(dat::ContinuousData, include_list, exclude_list)
    n_samples = size(dat.data, 1)
    keep = trues(n_samples) # Start with all true

    # Apply exclusions
    if !isnothing(exclude_list)
        for col in exclude_list
            if col in propertynames(dat.data) && eltype(dat.data[!, col]) == Bool
                keep .&= .!dat.data[!, col] # Set to false where exclude is true
            else
                @warn "Exclude column '$col' not found or not Bool type."
            end
        end
    end

    # Apply inclusions (if specified, overrides excludes implicitly)
    if !isnothing(include_list)
        include_mask = falses(n_samples)
        for col in include_list
             if col in propertynames(dat.data) && eltype(dat.data[!, col]) == Bool
                include_mask .|= dat.data[!, col] # Set to true where include is true
            else
                @warn "Include column '$col' not found or not Bool type."
            end
        end
        keep .&= include_mask # Keep only where include is true (within remaining)
    end

    return findall(keep)
end

"""
    _simple_findpeaks(data::AbstractVector; min_prominence_std::Real=2.0)

Basic peak finder. Finds indices where data point is greater than its
immediate neighbors and exceeds mean + min_prominence_std * std(data).
Returns indices of peaks.
"""
function _simple_findpeaks(data::AbstractVector; min_prominence_std::Real=2.0)
    if length(data) < 3
        return Int[]
    end
    peaks = Int[]
    mean_val = mean(data)
    std_val = std(data)
    # Handle zero std deviation case
    threshold = (std_val ≈ 0) ? mean_val : mean_val + min_prominence_std * std_val

    for i in 2:(length(data)-1)
        if data[i] > data[i-1] && data[i] > data[i+1] && data[i] > threshold
            push!(peaks, i)
        end
    end
    return peaks
end

"""
    identify_ekg_components(ica_result::InfoIca, dat::ContinuousData;
                              min_bpm::Real=40, max_bpm::Real=120,
                              min_prominence_std::Real=2.5,
                              min_peaks::Int=10,
                              max_ibi_std_s::Real=0.05,
                              include_samples::Union{Nothing,Vector{Symbol}} = nothing,
                              exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value])

Identify ICA components potentially related to EKG artifacts based on peak detection
and interval regularity, using only samples consistent with ICA calculation.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data (needed for sampling rate `fs` and sample selection columns).

# Keyword Arguments
- `min_bpm::Real`: Minimum plausible heart rate in beats per minute (default: 40).
- `max_bpm::Real`: Maximum plausible heart rate in beats per minute (default: 120).
- `min_prominence_std::Real`: Minimum peak prominence in standard deviations above mean (default: 2.5).
- `min_peaks::Int`: Minimum number of prominent peaks required within plausible heart rate range (default: 10).
- `max_ibi_std_s::Real`: Maximum standard deviation of the inter-beat intervals (in seconds) for component to be flagged (default: 0.05).
- `include_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to *include*. Defaults to `nothing` (include all unless excluded).
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to *exclude*. Defaults to `[:is_extreme_value]`.

# Returns
- `Vector{Int}`: Sorted vector of indices identified as potential EKG components.
- `DataFrame`: DataFrame containing metrics for each component (calculated on the included samples):
  - `:Component`: Component index (1 to n).
  - `:num_peaks`: Number of detected prominent peaks.
  - `:num_valid_ibis`: Number of inter-beat intervals within the plausible BPM range.
  - `:mean_ibi_s`: Mean inter-beat interval in seconds (if num_valid_ibis > 0).
  - `:std_ibi_s`: Standard deviation of inter-beat intervals in seconds (if num_valid_ibis > 1).
  - `:is_ekg_artifact`: Boolean flag indicating if component met the criteria.
"""
function identify_ecg_components( # Renamed from identify_ecg_components
    ica_result::InfoIca,
    dat::ContinuousData;
    min_bpm::Real=40,
    max_bpm::Real=120,
    min_prominence_std::Real=2.5,
    min_peaks::Int=10,
    max_ibi_std_s::Real=0.15,
    include_samples::Union{Nothing,Vector{Symbol}} = nothing, # New arg
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value] # New arg
)
    # --- Data Preparation ---

    # 1. Determine samples to use, consistent with run_ica
    samples_to_use = _get_samples_to_use(dat, include_samples, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying include/exclude criteria. Cannot identify EKG components."
        # Return empty results
        return Int[], DataFrame(Component=Int[], num_peaks=Int[], num_valid_ibis=Int[], mean_ibi_s=Float64[], std_ibi_s=Float64[], is_ekg_artifact=Bool[])
    end

    # 2. Extract relevant data *only for these samples*
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]

    # 3. Prepare matrix for unmixing (on the subset)
    dat_matrix_subset = permutedims(Matrix(data_subset_df)) 
    dat_matrix_subset .-= mean(dat_matrix_subset, dims=2)
    dat_matrix_subset ./= ica_result.scale

    # 4. Calculate components activations *only for these samples*
    components_subset = ica_result.unmixing * dat_matrix_subset 
    # Note: n_components is the total number of components from ICA result
    n_components_total = size(ica_result.unmixing, 1)
    fs = dat.sample_rate 

    # Convert BPM to plausible IBI range in seconds
    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    # Store results
    metrics = []
    identified_ekg = Int[]

    # Loop through components, analyzing the activations *from the subset*
    for comp_idx in 1:n_components_total # Iterate through all component indices
        ts = components_subset[comp_idx, :] # Time series from *filtered* samples

        # Find prominent peaks in the absolute signal magnitude
        peak_indices = _simple_findpeaks(abs.(ts); min_prominence_std=min_prominence_std)
        
        num_peaks = length(peak_indices)
        mean_ibi = NaN
        std_ibi = NaN
        num_valid_ibis = 0
        is_ekg = false

        if num_peaks >= 2
            # Calculate Inter-Beat Intervals (IBIs) in seconds
            ibis_s = diff(peak_indices) ./ fs

            # Filter IBIs based on plausible heart rate range
            valid_ibi_mask = (ibis_s .>= min_ibi_s) .& (ibis_s .<= max_ibi_s)
            valid_ibis = ibis_s[valid_ibi_mask]
            num_valid_ibis = length(valid_ibis)

            if num_valid_ibis > 1 
                mean_ibi = mean(valid_ibis)
                std_ibi = std(valid_ibis)

                # Check criteria
                if num_valid_ibis >= (min_peaks - 1) && 
                   std_ibi <= max_ibi_std_s
                    is_ekg = true
                    push!(identified_ekg, comp_idx)
                end
            elseif num_valid_ibis == 1 
                 mean_ibi = valid_ibis[1]
                 std_ibi = 0.0 
                 if num_valid_ibis >= (min_peaks - 1) && std_ibi <= max_ibi_std_s
                     is_ekg = true
                     push!(identified_ekg, comp_idx)
                 end
            end
        end

        # Store metrics for this component
        push!(metrics, (
            Component=comp_idx,
            num_peaks=num_peaks,
            num_valid_ibis=num_valid_ibis,
            mean_ibi_s=mean_ibi,
            std_ibi_s=std_ibi,
            is_ekg_artifact=is_ekg
        ))
    end

    # --- Construct Results ---
    metrics_df = DataFrame(metrics)
    # Ensure Component column matches 1:n_components_total if metrics is empty
    if isempty(metrics)
         metrics_df = DataFrame(Component=1:n_components_total, num_peaks=0, num_valid_ibis=0, mean_ibi_s=NaN, std_ibi_s=NaN, is_ekg_artifact=false)
    end
    sort!(identified_ekg)

    return identified_ekg, metrics_df
end



"""
    identify_spatial_kurtosis_components(ica_result::InfoIca, dat::ContinuousData;
                                      exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                      z_threshold::Float64 = 3.0)

Identify ICA components with high spatial kurtosis (localized, spot-like activity).

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `z_threshold::Float64`: Z-score threshold for identifying high spatial kurtosis components (default: 3.0).

# Returns
- `Vector{Int}`: Indices of components with high spatial kurtosis.
- `DataFrame`: DataFrame containing spatial kurtosis values and z-scores for all components.
"""
function identify_spatial_kurtosis_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    z_threshold::Float64 = 3.0
)
    # Calculate spatial kurtosis for each component's weights
    n_components = size(ica_result.mixing, 2)
    spatial_kurtosis = Float64[]
    
    for i in 1:n_components
        # Get component weights
        weights = ica_result.mixing[:, i]
        # Calculate kurtosis of the weights
        k = kurtosis(weights)
        push!(spatial_kurtosis, k)
    end

    # Calculate z-scores of spatial kurtosis values
    spatial_kurtosis_z = StatsBase.zscore(spatial_kurtosis)

    # Identify components with high spatial kurtosis (using z-scores)
    high_kurtosis_comps = findall(spatial_kurtosis_z .> z_threshold)  # Only positive deviations (localized activity)
    sort!(high_kurtosis_comps)

    # Create metrics DataFrame
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :SpatialKurtosis => spatial_kurtosis,
        :SpatialKurtosisZScore => spatial_kurtosis_z
    )

    return high_kurtosis_comps, metrics_df
end

"""
    plot_spatial_kurtosis_components(ica_result::InfoIca, dat::ContinuousData;
                                   exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                   z_threshold::Float64 = 3.0)

Plot spatial kurtosis z-scores for all ICA components and highlight those exceeding the threshold.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `z_threshold::Float64`: Z-score threshold for identifying high spatial kurtosis components (default: 3.0).

# Returns
- `fig::Figure`: The Makie Figure containing the spatial kurtosis plot.
"""
function plot_spatial_kurtosis_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    z_threshold::Float64 = 3.0
)
    # Get spatial kurtosis values
    high_kurtosis_comps, metrics_df = identify_spatial_kurtosis_components(
        ica_result, dat;
        exclude_samples = exclude_samples,
        z_threshold = z_threshold
    )

    # Create figure
    fig = Figure(size=(800, 400))
    
    # Plot spatial kurtosis z-scores
    ax = Axis(
        fig[1, 1],
        xlabel = "Component",
        ylabel = "Spatial Kurtosis Z-Score",
        title = "Component Spatial Kurtosis Z-Scores"
    )
    
    # Plot all components
    scatter!(ax, metrics_df.Component, metrics_df.SpatialKurtosisZScore, color=:gray)
    
    # Highlight high kurtosis components
    if !isempty(high_kurtosis_comps)
        high_kurtosis_values = metrics_df[in.(metrics_df.Component, Ref(high_kurtosis_comps)), :SpatialKurtosisZScore]
        scatter!(ax, high_kurtosis_comps, high_kurtosis_values, color=:red, markersize=8)
        
        # Add labels for high kurtosis components
        for (i, comp) in enumerate(high_kurtosis_comps)
            text!(ax, comp, high_kurtosis_values[i], text=string(comp), 
                  color=:red, align=(:center,:bottom), fontsize=10)
        end
    end
    
    # Add threshold line
    hlines!(ax, [z_threshold], color=:red, linestyle=:dash)
    
    # Add reference line for mean
    hlines!(ax, [0.0], color=:gray, linestyle=:dot, label="Mean")
    
    # Add legend with correct position specification
    axislegend(ax, position=(1.0, 1.0))
    
    return fig
end

"""
    plot_ecg_component_features(identified_comps::Vector{Int}, metrics_df::DataFrame;
                              min_bpm::Real=40, max_bpm::Real=120)

Plot metrics used for ECG component identification.

# Arguments
- `identified_comps::Vector{Int}`: Vector of identified ECG component indices.
- `metrics_df::DataFrame`: DataFrame containing ECG metrics (from identify_ecg_components).
- `min_bpm::Real`: Minimum plausible heart rate in BPM (default: 40).
- `max_bpm::Real`: Maximum plausible heart rate in BPM (default: 120).

# Returns
- `fig::Figure`: The Makie Figure containing the ECG metrics plots.
"""
function plot_ecg_component_features(
    identified_comps::Vector{Int},
    metrics_df::DataFrame;
    min_bpm::Real=40,
    max_bpm::Real=120
)
    # Create figure with two subplots
    fig = Figure(size=(1000, 400))
    
    # Plot 1: Number of peaks and valid IBIs
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Component",
        ylabel = "Count",
        title = "Peak and IBI Counts"
    )
    
    # Plot number of peaks
    scatter!(ax1, metrics_df.Component, metrics_df.num_peaks, 
             color=:gray, label="Total Peaks")
    
    # Plot number of valid IBIs
    scatter!(ax1, metrics_df.Component, metrics_df.num_valid_ibis,
             color=:blue, label="Valid IBIs")
    
    # Highlight identified components
    if !isempty(identified_comps)
        # Get metrics for identified components
        identified_metrics = metrics_df[in.(metrics_df.Component, Ref(identified_comps)), :]
        
        # Plot identified components with larger markers
        scatter!(ax1, identified_metrics.Component, identified_metrics.num_peaks,
                color=:red, markersize=8, label="Identified Components (Peaks)")
        scatter!(ax1, identified_metrics.Component, identified_metrics.num_valid_ibis,
                color=:red, markersize=8, label="Identified Components (IBIs)")
        
        # Add component numbers as labels
        for (i, comp) in enumerate(identified_comps)
            row = metrics_df[metrics_df.Component .== comp, :]
            text!(ax1, comp, row.num_peaks[1], text=string(comp),
                  color=:red, align=(:center,:bottom), fontsize=10)
        end
    end
    
    # Add legend
    axislegend(ax1, position=(1.0, 1.0))
    
    # Plot 2: IBI Statistics
    ax2 = Axis(
        fig[1, 2],
        xlabel = "Component",
        ylabel = "Seconds",
        title = "Inter-Beat Interval Statistics"
    )
    
    # Plot mean IBI
    scatter!(ax2, metrics_df.Component, metrics_df.mean_ibi_s,
             color=:blue, label="Mean IBI")
    
    # Plot IBI standard deviation
    scatter!(ax2, metrics_df.Component, metrics_df.std_ibi_s,
             color=:green, label="IBI Std Dev")
    
    # Add reference lines for plausible heart rate range
    min_ibi = 60.0 / max_bpm  # Convert BPM to seconds
    max_ibi = 60.0 / min_bpm
    hlines!(ax2, [min_ibi, max_ibi], color=:gray, linestyle=:dash,
            label="Plausible IBI Range")
    
    # Highlight identified components
    if !isempty(identified_comps)
        identified_metrics = metrics_df[in.(metrics_df.Component, Ref(identified_comps)), :]
        
        # Plot identified components with larger markers
        scatter!(ax2, identified_metrics.Component, identified_metrics.mean_ibi_s,
                color=:red, markersize=8, label="Identified Components (Mean)")
        scatter!(ax2, identified_metrics.Component, identified_metrics.std_ibi_s,
                color=:red, markersize=8, label="Identified Components (Std)")
        
        # Add component numbers as labels
        for (i, comp) in enumerate(identified_comps)
            row = metrics_df[metrics_df.Component .== comp, :]
            text!(ax2, comp, row.mean_ibi_s[1], text=string(comp),
                  color=:red, align=(:center,:bottom), fontsize=10)
        end
    end
    
    # Add legend
    axislegend(ax2, position=(1.0, 1.0))
    
    return fig
end

"""
    identify_line_noise_components(ica_result::InfoIca, dat::ContinuousData;
                                 exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Float64=3.0,
                                 min_harmonic_power::Real=0.5)

Identify ICA components with strong line noise characteristics.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz (default: 50.0 for European power).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum power ratio of harmonics relative to fundamental (default: 0.5).

# Returns
- `Vector{Int}`: Indices of components with strong line noise characteristics.
- `DataFrame`: DataFrame containing spectral metrics for all components.
"""
function identify_line_noise_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    z_threshold::Float64=3.0,
    min_harmonic_power::Real=0.5
)
    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], DataFrame()
    end

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate components for valid samples
    components = ica_result.unmixing * dat_matrix
    n_components = size(components, 1)
    fs = dat.sample_rate

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, but not too large)
    nfft = min(nextpow(2, size(components, 2)), 2^16)  # Cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, fs)
    psd = zeros(length(freqs), n_components)
    
    for i in 1:n_components
        # Zero-pad or truncate to nfft points
        signal = components[i, :]
        if length(signal) > nfft
            signal = signal[1:nfft]
        elseif length(signal) < nfft
            signal = [signal; zeros(nfft - length(signal))]
        end
        psd[:, i] = abs2.(FFTW.rfft(signal))
    end

    # Find indices for line frequency and harmonics
    line_idx = findmin(abs.(freqs .- line_freq))[2]
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)
    
    # Calculate metrics for each component
    metrics = []
    for i in 1:n_components
        # Get power at line frequency and surrounding band
        line_power = mean(psd[line_band, i])
        
        # Calculate power in surrounding bands (excluding line frequency)
        surrounding_bands = setdiff(1:length(freqs), line_band)
        surrounding_power = mean(psd[surrounding_bands, i])
        
        # Calculate power ratio
        power_ratio = line_power / (surrounding_power + eps())
        
        # Check for harmonics (2x and 3x line frequency)
        harmonic_powers = Float64[]
        for h in 2:3
            harmonic_freq = line_freq * h
            harmonic_idx = findmin(abs.(freqs .- harmonic_freq))[2]
            harmonic_band = findall(abs.(freqs .- harmonic_freq) .<= freq_bandwidth)
            harmonic_power = mean(psd[harmonic_band, i])
            push!(harmonic_powers, harmonic_power / line_power)
        end
        
        # Store metrics
        push!(metrics, (
            Component=i,
            LinePower=line_power,
            SurroundingPower=surrounding_power,
            PowerRatio=power_ratio,
            Harmonic2Ratio=harmonic_powers[1],
            Harmonic3Ratio=harmonic_powers[2]
        ))
    end

    # Create metrics DataFrame
    metrics_df = DataFrame(metrics)
    
    # Calculate z-scores of power ratios
    power_ratio_z = StatsBase.zscore(metrics_df.PowerRatio)
    metrics_df[!, :PowerRatioZScore] = power_ratio_z

    # Identify components with strong line noise characteristics
    line_noise_comps = findall(power_ratio_z .> z_threshold)
    
    # Additional check for harmonics
    if !isempty(line_noise_comps)
        harmonic_mask = (metrics_df.Harmonic2Ratio .> min_harmonic_power) .| 
                       (metrics_df.Harmonic3Ratio .> min_harmonic_power)
        line_noise_comps = intersect(line_noise_comps, findall(harmonic_mask))
    end
    
    sort!(line_noise_comps)

    return line_noise_comps, metrics_df
end

"""
    plot_line_noise_components(ica_result::InfoIca, dat::ContinuousData;
                             exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                             line_freq::Real=50.0,
                             freq_bandwidth::Real=1.0,
                             z_threshold::Float64=3.0,
                             min_harmonic_power::Real=0.5)

Plot spectral metrics used for line noise component identification.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum power ratio of harmonics relative to fundamental (default: 0.5).

# Returns
- `fig::Figure`: The Makie Figure containing the line noise metrics plots.
"""
function plot_line_noise_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    z_threshold::Float64=3.0,
    min_harmonic_power::Real=0.5
)
    # Get line noise components and metrics
    line_noise_comps, metrics_df = identify_line_noise_components(
        ica_result, dat;
        exclude_samples=exclude_samples,
        line_freq=line_freq,
        freq_bandwidth=freq_bandwidth,
        z_threshold=z_threshold,
        min_harmonic_power=min_harmonic_power
    )

    # Create figure with two subplots
    fig = Figure(size=(1000, 400))
    
    # Plot 1: Power Ratio Z-Scores
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Component",
        ylabel = "Power Ratio Z-Score",
        title = "Line Frequency Power Ratio Z-Scores"
    )
    
    # Plot all components with label
    scatter!(ax1, metrics_df.Component, metrics_df.PowerRatioZScore, 
             color=:gray, label="All Components")
    
    # Highlight identified components with label
    if !isempty(line_noise_comps)
        identified_metrics = metrics_df[in.(metrics_df.Component, Ref(line_noise_comps)), :]
        scatter!(ax1, identified_metrics.Component, identified_metrics.PowerRatioZScore,
                color=:red, markersize=8, label="Line Noise Components")
        
        # Add component numbers as labels
        for (i, comp) in enumerate(line_noise_comps)
            row = metrics_df[metrics_df.Component .== comp, :]
            text!(ax1, comp, row.PowerRatioZScore[1], text=string(comp),
                  color=:red, align=(:center,:bottom), fontsize=10)
        end
    end
    
    # Add threshold line with label
    hlines!(ax1, [z_threshold], color=:red, linestyle=:dash, label="Threshold")
    
    # Add legend
    axislegend(ax1, position=(1.0, 1.0))
    
    # Plot 2: Harmonic Ratios
    ax2 = Axis(
        fig[1, 2],
        xlabel = "Component",
        ylabel = "Power Ratio",
        title = "Harmonic Power Ratios"
    )
    
    # Plot harmonic ratios with labels
    scatter!(ax2, metrics_df.Component, metrics_df.Harmonic2Ratio,
             color=:blue, label="2nd Harmonic")
    scatter!(ax2, metrics_df.Component, metrics_df.Harmonic3Ratio,
             color=:green, label="3rd Harmonic")
    
    # Add reference line for minimum harmonic power with label
    hlines!(ax2, [min_harmonic_power], color=:gray, linestyle=:dash,
            label="Min Harmonic Power")
    
    # Highlight identified components with label
    if !isempty(line_noise_comps)
        identified_metrics = metrics_df[in.(metrics_df.Component, Ref(line_noise_comps)), :]
        scatter!(ax2, identified_metrics.Component, identified_metrics.Harmonic2Ratio,
                color=:red, markersize=8, label="Line Noise Components")
        scatter!(ax2, identified_metrics.Component, identified_metrics.Harmonic3Ratio,
                color=:red, markersize=8)
    end
    
    # Add legend
    axislegend(ax2, position=(1.0, 1.0))
    
    return fig
end

"""
    plot_component_spectrum(ica_result::InfoIca, dat::ContinuousData, comp_idx::Int;
                          exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                          line_freq::Real=50.0,
                          freq_bandwidth::Real=1.0,
                          window_size::Int=1024,
                          overlap::Real=0.5,
                          max_freq::Real=100.0)

Plot the power spectrum of a specific ICA component.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.
- `comp_idx::Int`: Index of the component to plot.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
"""
function plot_component_spectrum(
    ica_result::InfoIca,
    dat::ContinuousData,
    comp_idx::Int;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0
)
    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot plot component spectrum."
        return Figure()
    end

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate component activation
    components = ica_result.unmixing * dat_matrix
    fs = dat.sample_rate

    # Get the component time series
    signal = components[comp_idx, :]

    # Calculate power spectrum using Welch's method
    noverlap = Int(round(window_size * overlap))
    pgram = DSP.welch_pgram(signal, window_size, noverlap; fs=fs)
    freqs = DSP.freq(pgram)
    psd = DSP.power(pgram)

    # Create figure
    fig = Figure(size=(800, 400))
    
    # Plot power spectrum
    ax = Axis(
        fig[1, 1],
        xlabel = "Frequency (Hz)",
        ylabel = "Power Spectral Density (μV²/Hz)",
        title = "Power Spectrum of Component $comp_idx"
    )
    
    # Plot the spectrum
    lines!(ax, freqs, psd, color=:black, label="Power Spectrum")
    
    # Highlight line frequency and harmonics
    for h in 1:3
        freq = line_freq * h
        if freq <= max_freq
            # Add vertical line
            vlines!(ax, [freq], color=:red, linestyle=:dash, 
                   label=h==1 ? "Line Frequency" : "Harmonic")
            
            # Add shaded region around the frequency
            band_x = [freq-freq_bandwidth, freq+freq_bandwidth]
            band_y = [0, maximum(psd)]
            poly!(ax, [Point2f(band_x[1], band_y[1]), 
                      Point2f(band_x[2], band_y[1]),
                      Point2f(band_x[2], band_y[2]),
                      Point2f(band_x[1], band_y[2])],
                  color=(:red, 0.1))
            
            # Add frequency label
            text!(ax, freq, maximum(psd), text="$freq Hz",
                  color=:red, align=(:center, :bottom))
        end
    end
    
    # Set x-axis limits
    xlims!(ax, (0, max_freq))
    
    # Add legend
    axislegend(ax, position=(1.0, 1.0))
    
    return fig
end

"""
    plot_channel_spectrum(dat::ContinuousData, channel::Symbol;
                         exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                         line_freq::Real=50.0,
                         freq_bandwidth::Real=1.0,
                         window_size::Int=1024,
                         overlap::Real=0.5,
                         max_freq::Real=100.0)

Plot the power spectrum of a specific EEG channel.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `channel::Symbol`: The channel to plot (must be a column name in dat.data).

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
"""
function plot_channel_spectrum(
    dat::ContinuousData,
    channel::Symbol;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0
)
    # Check if channel exists
    if !(channel in propertynames(dat.data))
        error("Channel $channel not found in data")
    end

    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot plot channel spectrum."
        return Figure()
    end

    # Get the channel data
    signal = dat.data[samples_to_use, channel]
    fs = dat.sample_rate

    # Calculate power spectrum using Welch's method
    noverlap = Int(round(window_size * overlap))
    pgram = DSP.welch_pgram(signal, window_size, noverlap; fs=fs)
    freqs = DSP.freq(pgram)
    psd = DSP.power(pgram)

    # Create figure
    fig = Figure(size=(800, 400))
    
    # Plot power spectrum
    ax = Axis(
        fig[1, 1],
        xlabel = "Frequency (Hz)",
        ylabel = "Power Spectral Density (μV²/Hz)",
        title = "Power Spectrum of Channel $channel"
    )
    
    # Plot the spectrum
    lines!(ax, freqs, psd, color=:black, label="Power Spectrum")
    
    # Highlight line frequency and harmonics
    for h in 1:3
        freq = line_freq * h
        if freq <= max_freq
            # Add vertical line
            vlines!(ax, [freq], color=:red, linestyle=:dash, 
                   label=h==1 ? "Line Frequency" : "Harmonic")
            
            # Add shaded region around the frequency
            band_x = [freq-freq_bandwidth, freq+freq_bandwidth]
            band_y = [0, maximum(psd)]
            poly!(ax, [Point2f(band_x[1], band_y[1]), 
                      Point2f(band_x[2], band_y[1]),
                      Point2f(band_x[2], band_y[2]),
                      Point2f(band_x[1], band_y[2])],
                  color=(:red, 0.1))
            
            # Add frequency label
            text!(ax, freq, maximum(psd), text="$freq Hz",
                  color=:red, align=(:center, :bottom))
        end
    end
    
    # Set x-axis limits
    xlims!(ax, (0, max_freq))
    
    # Add legend
    axislegend(ax, position=(1.0, 1.0))
    
    return fig
end

"""
    plot_channel_spectrum(dat::ContinuousData, channel::Union{Symbol,Nothing}=nothing;
                         exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                         line_freq::Real=50.0,
                         freq_bandwidth::Real=1.0,
                         window_size::Int=1024,
                         overlap::Real=0.5,
                         max_freq::Real=100.0)

Plot the power spectrum of EEG channel(s).

# Arguments
- `dat::ContinuousData`: The continuous data.
- `channel::Union{Symbol,Nothing}`: The channel to plot (must be a column name in dat.data). If nothing, plots all EEG channels.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot(s).
"""
function plot_channel_spectrum(
    dat::ContinuousData,
    channel::Union{Symbol,Nothing}=nothing;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0
)
    println("Plotting channel spectrum")
    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot plot channel spectrum."
        return Figure()
    end

    # Get channels to plot
    channels_to_plot = if isnothing(channel)
        dat.layout.label
    else
        # Check if specified channel exists
        if !(channel in propertynames(dat.data))
            error("Channel $channel not found in data")
        end
        [channel]
    end

    # Calculate power spectra for all channels
    fs = dat.sample_rate
    noverlap = Int(round(window_size * overlap))
    spectra = Dict{Symbol, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    for ch in channels_to_plot
        signal = dat.data[samples_to_use, ch]
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs=fs)
        spectra[ch] = (DSP.freq(pgram), DSP.power(pgram))
    end

    # Create figure
    fig = Figure(size=(1000, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = "Frequency (Hz)",
        ylabel = "Power Spectral Density (μV²/Hz)",
        title = isnothing(channel) ? "Power Spectra of All Channels" : "Power Spectrum of Channel $channel"
    )
    
    # Plot spectra for all channels
    for (i, ch) in enumerate(channels_to_plot)
        freqs, psd = spectra[ch]
        if length(channels_to_plot) > 1
            lines!(ax, freqs, psd, label=string(ch))
        else
            lines!(ax, freqs, psd)
        end
    end
    
    # Highlight line frequency and harmonics
    for h in 1:3
        freq = line_freq * h
        if freq <= max_freq
            # Add vertical line without label
            vlines!(ax, [freq], color=:red, linestyle=:dash, label="")
            
            # Add shaded region around the frequency
            band_x = [freq-freq_bandwidth, freq+freq_bandwidth]
            band_y = [0, maximum([maximum(psd) for (_, psd) in values(spectra)])]
            poly!(ax, [Point2f(band_x[1], band_y[1]), 
                      Point2f(band_x[2], band_y[1]),
                      Point2f(band_x[2], band_y[2]),
                      Point2f(band_x[1], band_y[2])],
                  color=(:red, 0.1))
            
            # Add frequency label
            text!(ax, freq, band_y[2], text="$freq Hz",
                  color=:red, align=(:center, :bottom))
        end
    end
    
    # Set x-axis limits
    xlims!(ax, (0, max_freq))
    
    # Add legend only if we have multiple channels
    if length(channels_to_plot) > 1
        axislegend(ax, position=(1.0, 1.0))
    end
    
    return fig
end

# Helper function to plot spectrum in an axis
function _plot_spectrum!(ax::Axis, freqs::Vector{Float64}, psd::Vector{Float64}, 
                        line_freq::Real, freq_bandwidth::Real, max_freq::Real)
    # Plot the spectrum
    lines!(ax, freqs, psd, color=:black, label="Power Spectrum")
    
    # Highlight line frequency and harmonics
    for h in 1:3
        freq = line_freq * h
        if freq <= max_freq
            # Add vertical line without label
            vlines!(ax, [freq], color=:red, linestyle=:dash)
            
            # Add shaded region around the frequency
            band_x = [freq-freq_bandwidth, freq+freq_bandwidth]
            band_y = [0, maximum(psd)]
            poly!(ax, [Point2f(band_x[1], band_y[1]), 
                      Point2f(band_x[2], band_y[1]),
                      Point2f(band_x[2], band_y[2]),
                      Point2f(band_x[1], band_y[2])],
                  color=(:red, 0.1))
            
            # Add frequency label
            text!(ax, freq, maximum(psd), text="$freq Hz",
                  color=:red, align=(:center, :bottom))
        end
    end
    
    # Set x-axis limits
    xlims!(ax, (0, max_freq))
end

"""
    plot_selected_spectrum(selected_data::DataFrame, channel::Symbol; 
        line_freq::Real=50.0,
        freq_bandwidth::Real=1.0,
        window_size::Int=1024,
        overlap::Real=0.5,
        max_freq::Real=100.0
    )

Plot the power spectrum of a selected time region for a specific channel.

# Arguments
- `selected_data::DataFrame`: The selected time region data
- `channel::Symbol`: The channel to plot
- `line_freq::Real`: Line frequency (default: 50.0 Hz)
- `freq_bandwidth::Real`: Bandwidth around line frequency (default: 1.0 Hz)
- `window_size::Int`: Window size for Welch's method (default: 1024)
- `overlap::Real`: Overlap between windows (default: 0.5)
- `max_freq::Real`: Maximum frequency to plot (default: 100.0 Hz)

# Returns
- `fig::Figure`: The figure containing the plot
"""
function plot_selected_spectrum(
    selected_data::DataFrame,
    channel::Vector{Symbol};
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0
)
    # Get sampling rate from time column
    fs = 1 / mean(diff(selected_data.time))
    
    # Calculate minimum required samples
    min_samples = window_size  # Minimum samples needed for one window
    
    # Check if we have enough data
    if size(selected_data, 1) < min_samples
        @warn "Selected region is too short for the specified window size. Adjusting window size..."
        # Adjust window size to be half the data length, but ensure it's a power of 2
        window_size = 2^floor(Int, log2(size(selected_data, 1) / 2))
        if window_size < 32  # Set a minimum window size
            @warn "Selected region is too short for meaningful spectral analysis. Minimum window size of 32 samples required."
            return Figure(), Axis()
        end
    end
    
    # Get signal data for each channel
    spectra = Dict{Symbol, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    for ch in channel
        # Extract vector data from DataFrame
        signal = Vector{Float64}(selected_data[!, ch])
        
        # Calculate power spectrum using Welch's method
        noverlap = Int(round(window_size * overlap))
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs=fs)
        spectra[ch] = (DSP.freq(pgram), DSP.power(pgram))
    end
    
    # Create figure
    fig = Figure(size=(800, 400))
    ax = Axis(fig[1, 1],
        xlabel="Frequency (Hz)",
        ylabel="Power Spectral Density (μV²/Hz)",
        title="Power Spectrum of Selected Region"
    )
    
    # Plot spectra for all channels
    for (i, ch) in enumerate(channel)
        freqs, psd = spectra[ch]
        if length(channel) > 1
            lines!(ax, freqs, psd, label=string(ch))
        else
            lines!(ax, freqs, psd)
        end
    end
    
    # Add vertical lines for line frequency and harmonics
    for i in 1:3
        freq = i * line_freq
        if freq <= max_freq
            # Add vertical line without label
            vlines!(ax, [freq], color=:red, linestyle=:dash)
            
            # Add shaded region
            band_x = [freq - freq_bandwidth, freq + freq_bandwidth]
            band_y = [0, maximum([maximum(psd) for (_, psd) in values(spectra)])]
            vertices = [
                Point2f(band_x[1], band_y[1]),
                Point2f(band_x[2], band_y[1]),
                Point2f(band_x[2], band_y[2]),
                Point2f(band_x[1], band_y[2])
            ]
            poly!(ax, vertices, color=(:red, 0.2))
            
            # Add frequency label
            text!(ax, freq, band_y[2], text="$(Int(freq)) Hz", align=(:center, :bottom))
        end
    end
    
    # Set x-axis limits
    xlims!(ax, (0, max_freq))
    
    # Add legend if multiple channels
    if length(channel) > 1
        axislegend(ax, position=(1.0, 1.0))
    end
    display(GLMakie.Screen(), fig)
    return fig, ax
end


"""
    plot_selected_spectrum(selected_data::DataFrame, channel::Symbol;
                          line_freq::Real=50.0,
                          freq_bandwidth::Real=1.0,
                          window_size::Int=1024,
                          overlap::Real=0.5,
                          max_freq::Real=100.0)

Plot the power spectrum of a selected time region for a specific channel.

# Arguments
- `selected_data::DataFrame`: The selected time region data.
- `channel::Symbol`: The channel to plot (must be a column name in selected_data).

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
"""
function plot_selected_spectrum(
    selected_data::DataFrame,
    channel::Symbol;
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0
)
    # Check if channel exists
    if !(channel in propertynames(selected_data))
        error("Channel $channel not found in data")
    end

    # Get the channel data as a vector
    signal = Vector{Float64}(selected_data[!, channel])
    
    # Calculate sampling rate from time column
    fs = 1 / mean(diff(selected_data.time))

    # Calculate power spectrum using Welch's method
    noverlap = Int(round(window_size * overlap))
    pgram = DSP.welch_pgram(signal, window_size, noverlap; fs=fs)
    freqs = DSP.freq(pgram)
    psd = DSP.power(pgram)

    # Create figure
    fig = Figure(size=(800, 400))
    
    # Plot power spectrum
    ax = Axis(
        fig[1, 1],
        xlabel = "Frequency (Hz)",
        ylabel = "Power Spectral Density (μV²/Hz)",
        title = "Power Spectrum of Channel $channel"
    )
    
    # Plot the spectrum
    lines!(ax, freqs, psd, color=:black)
    
    # Highlight line frequency and harmonics
    for h in 1:3
        freq = line_freq * h
        if freq <= max_freq
            # Add vertical line without label
            vlines!(ax, [freq], color=:red, linestyle=:dash)
            
            # Add shaded region around the frequency
            band_x = [freq-freq_bandwidth, freq+freq_bandwidth]
            band_y = [0, maximum(psd)]
            poly!(ax, [Point2f(band_x[1], band_y[1]), 
                      Point2f(band_x[2], band_y[1]),
                      Point2f(band_x[2], band_y[2]),
                      Point2f(band_x[1], band_y[2])],
                  color=(:red, 0.1))
            
            # Add frequency label
            text!(ax, freq, maximum(psd), text="$freq Hz",
                  color=:red, align=(:center, :bottom))
        end
    end
    
    # Set x-axis limits
    xlims!(ax, (0, max_freq))
    display(GLMakie.Screen(), fig)
    return fig, ax
end


