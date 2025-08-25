"""
    plot_ica_topoplot(ica; ...)

Plot multiple ICA component topographies in a grid layout within a new Figure.

# Arguments
- `ica::InfoIca`: The ICA result object (contains layout information).

# Keyword Arguments
- `component_selection::Function`: Function that returns boolean vector for component filtering (default: all components).
- `dims=nothing`: Tuple specifying grid dimensions (rows, cols). If `nothing`, calculates best square-ish grid.
- `head_kwargs::Dict`: Controls the head outline for all plots.
- `point_kwargs::Dict`: Controls channel markers for all plots.
- `label_kwargs::Dict`: Controls channel labels for all plots.
- `topo_kwargs::Dict`: Controls the topography plots (e.g., `colormap`, `gridscale`, `num_levels`, `nan_color`).
- `colorbar_kwargs::Dict`: Controls colorbars (e.g., `plot_colorbar`, `width`, `colorbar_plot_numbers`). `colorbar_plot_numbers` specifies which plot indices (1-based) should have a visible colorbar.
- `use_global_scale::Bool=false`: If `true`, all topoplots share the same color scale based on the min/max across all plotted components. If `false`, each plot is scaled independently.

# Returns
- `fig::Figure`: The generated Makie Figure containing the grid of topoplots.

# Examples

## Basic Usage
```julia
# Plot first 10 components (default)
fig = plot_ica_topoplot(ica_result)

# Plot specific range of components
fig = plot_ica_topoplot(ica_result, component_selection = components(5:15))

# Plot specific components
fig = plot_ica_topoplot(ica_result, component_selection = components([1, 3, 5, 7]))

# Plot all components (if screen can handle it)
fig = plot_ica_topoplot(ica_result, component_selection = components())
```

## Advanced Selection
```julia
# Plot components with custom selection
fig = plot_ica_topoplot(ica_result, 
    component_selection = components(1:10)  # First 10 components
)

# Plot even-numbered components
fig = plot_ica_topoplot(ica_result, 
    component_selection = components(2:2:20)  # Even components 2, 4, 6, ..., 20
)
```
"""
function plot_ica_topoplot(
    ica;
    component_selection::Function = components(),
    dims = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    use_global_scale = false,
    method = :multiquadratic,  # :multiquadratic or :spherical_spline
    display_plot = true,
)


    # ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(ica.layout)
    _ensure_coordinates_3d!(ica.layout)


    # default kwargs and merged kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 100, :num_levels => 20, :nan_color => :transparent)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = topo_kwargs[:gridscale]
    num_levels = topo_kwargs[:num_levels]
    topo_kwargs = filter(p -> p.first ∉ (:gridscale, :num_levels, :levels), topo_kwargs)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30, :colorbar_plot_numbers => [])
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = colorbar_kwargs[:plot_colorbar]
    colorbar_plot_numbers = colorbar_kwargs[:colorbar_plot_numbers]
    colorbar_kwargs = filter(p -> p.first ∉ (:plot_colorbar, :colorbar_plot_numbers), colorbar_kwargs)

    # Get selected components using the helper function
    comps = get_selected_components(ica, component_selection)

    # Create figure with reduced margins
    fig = Figure(
        figure_padding = (0, 60, 0, 0), # Keep padding for right margin
        backgroundcolor = :white,
    )

    # Deal with plot dimensions
    if isnothing(dims)
        dims = best_rect(length(comps))
    end

    # Validate dimensions
    if length(dims) != 2 || any(dims .<= 0)
        throw(ArgumentError("Invalid dimensions: $dims. Expected [rows, cols] with positive values."))
    end

    # Ensure we have enough grid cells
    total_cells = dims[1] * dims[2]
    if total_cells < length(comps)
        throw(ArgumentError("Grid dimensions $dims provide $total_cells cells but need $(length(comps))."))
    end

    # Calculate all topo data first
    n_components = length(comps)
    all_data = Vector{Matrix{Float64}}(undef, n_components)
    for i in eachindex(comps)
        all_data[i] = _prepare_ica_topo_data(ica, comps[i], method, gridscale)
    end

    # Calculate levels for all plots (global or local)
    levels_result = _calculate_ica_topo_levels(all_data, use_global_scale, num_levels)

    # Create grid layouts and add plots in a single pass
    grids = Vector{GridLayout}(undef, n_components)
    for i in eachindex(comps)

        row, col = divrem(i - 1, dims[2]) .+ (1, 1)
        grids[i] = fig[row, col] = GridLayout()

        # grid layout and axis
        ax = Axis(grids[i][1, 1], title = @sprintf("IC %d (%.1f%%)", comps[i], ica.variance[comps[i]] * 100))

        # Get pre-calculated data and levels
        data = all_data[i]
        levels_to_use = use_global_scale ? levels_result : levels_result[i]

        co = _plot_ica_topo_on_axis!(
            ax,
            fig,
            data,
            ica,
            levels_to_use;
            gridscale = gridscale,
            colormap = topo_kwargs[:colormap],
            nan_color = topo_kwargs[:nan_color],
            head_kwargs = head_kwargs,
            point_kwargs = point_kwargs,
            label_kwargs = label_kwargs,
            method = method,
        )

        # Do we want to add a colourbar?
        if plot_colorbar

            width = get(colorbar_kwargs, :width, 10)
            height = get(colorbar_kwargs, :height, Relative(0.8))
            ticksize = get(colorbar_kwargs, :ticklabelsize, 10)

            # Create colorbar or placeholder in second column
            # TODO: I do not really follow colourbar logic in Makie
            # This seems a bit awkward!
            if i in colorbar_plot_numbers # Use accessed list
                Colorbar(
                    grids[i][1, 2],
                    co;
                    width = width,
                    height = height,
                    ticklabelsize = ticksize,
                    colorbar_kwargs...,
                )
            else
                Box(
                    grids[i][1, 2];
                    width = width,
                    height = height,
                    color = :transparent,
                    strokewidth = 0,
                    colorbar_kwargs...,
                )
            end

            # Set column sizes after creating the colorbar column
            colsize!(grids[i], 1, Relative(0.85))  # Plot column
            colsize!(grids[i], 2, Relative(0.15))  # Colorbar column

        end

    end

    if display_plot
        display_figure(fig)
    end
    return fig

end


# Create a state structure to hold the visualization state
mutable struct IcaComponentState
    # Data
    dat::ContinuousData
    ica_result::InfoIca
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

    #  Topoplot Kwargs 
    topo_gridscale::Int
    topo_colormap::Symbol
    topo_num_levels::Int
    topo_nan_color::Union{Symbol,Makie.Colorant}
    topo_method::Symbol  # :multiquadratic or :spherical_spline
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
        component_selection,
        n_visible_components,
        window_size;
        topo_kwargs = Dict(),
        head_kwargs = Dict(),
        point_kwargs = Dict(),
        label_kwargs = Dict(),
        method = :multiquadratic,  # :multiquadratic or :spherical_spline
    )
        # Prepare data matrix
        dat_matrix = prepare_ica_data_matrix(dat, ica_result)
        component_data = ica_result.unmixing * dat_matrix
        total_components = size(component_data, 1)

        # Use the layout from the ICA result (already subsetted to match the channels used in ICA)
        subsetted_layout = ica_result.layout

        # Create observables
        comp_start = Observable(1)
        use_global_scale = Observable(false)
        invert_scale = Observable(false)

        # Find index closest to time 0 to center the initial view
        time_zero_idx = findmin(abs.(dat.data.time))[2]
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

        # --- Process and store Topo Kwargs ---
        topo_defaults = Dict(:gridscale => 300, :colormap => :jet, :num_levels => 20, :nan_color => :transparent)
        processed_topo_kwargs = merge(topo_defaults, topo_kwargs)

        head_defaults = Dict(:color => :black, :linewidth => 2)
        processed_head_kwargs = merge(head_defaults, head_kwargs)

        point_defaults = Dict(:plot_points => false)
        processed_point_kwargs = merge(point_defaults, point_kwargs)

        label_defaults = Dict(:plot_labels => false)
        processed_label_kwargs = merge(label_defaults, label_kwargs)

        # Store processed values
        topo_gridscale = processed_topo_kwargs[:gridscale]
        topo_colormap = processed_topo_kwargs[:colormap]
        topo_num_levels = processed_topo_kwargs[:num_levels]
        topo_nan_color = processed_topo_kwargs[:nan_color]
        topo_method = method
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
            # Pass stored kwargs
            topo_gridscale,
            topo_colormap,
            topo_num_levels,
            topo_nan_color,
            topo_method,
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
- `component_selection=components()`: Component selection predicate (see `components()`). Defaults to all available components.
- `n_visible_components::Int=10`: Number of components visible vertically at once (controls the "window" size for scrolling).
- `window_size::Int=2000`: Initial time window size in samples.
- `topo_kwargs::Dict=Dict()`: Keyword arguments passed down for topography plots (see `_plot_topo_on_axis!`).
- `head_kwargs::Dict=Dict()`: Keyword arguments passed down for head outlines.
- `point_kwargs::Dict=Dict()`: Keyword arguments passed down for channel markers.
- `label_kwargs::Dict=Dict()`: Keyword arguments passed down for channel labels.

# Returns
- `fig::Figure`: The Makie Figure object containing the interactive plot.
"""
function plot_ica_component_activation(
    dat::ContinuousData,
    ica_result::InfoIca;
    component_selection = components(),
    n_visible_components::Int = 10,
    window_size::Int = 2000,
    topo_kwargs = Dict(),
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    method = :multiquadratic,  # :multiquadratic or :spherical_spline
)
    # Pass kwargs to constructor
    state = IcaComponentState(
        dat,
        ica_result,
        component_selection,
        n_visible_components,
        window_size;
        topo_kwargs = topo_kwargs,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
        method = method,
    )

    # Create figure with padding on the right for margin
    fig = Figure(figure_padding = (0, 20, 0, 0))

    # Setup GUI interface
    create_component_plots!(fig, state)
    add_navigation_controls!(fig, state)
    add_navigation_sliders!(fig, state)
    add_channel_menu!(fig, state)
    setup_keyboard_interactions!(fig, state)

    colsize!(fig.layout, 1, Relative(0.15))
    colsize!(fig.layout, 2, Relative(0.85))
    rowgap!(fig.layout, state.n_visible_components, 30)

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
    dat_matrix = permutedims(Matrix(dat.data[!, ica_result.layout.data.label]))
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
    layout::Layout,
    method::Symbol,
    levels;
    gridscale = 200,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    kwargs...,
)

    # Validate gridscale
    gridscale <= 0 && throw(ArgumentError("gridscale must be positive, got $gridscale"))

    # Calculate contour/coord range based on interpolation method
    contour_range_multiplier = method == :spherical_spline ? 4.0 : 2.0
    contour_range = DEFAULT_HEAD_RADIUS * contour_range_multiplier
    coord_range = range(-contour_range, contour_range, length = gridscale)

    co = contourf!(ax, coord_range, coord_range, data; levels = levels, kwargs...)
    plot_layout_2d!(
        fig,
        ax,
        layout;
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    hidedecorations!(ax, grid = false)

    return co

end

# Create a simpler override specifically for the component viewer
# Internal function for viewer integration
function _plot_ica_topo_in_viewer!(
    fig,
    topo_ax,
    ica,
    comp_idx;
    use_global_scale = false,
    gridscale = 100,
    colormap = :jet,
    num_levels = 20,
    nan_color = :transparent,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    method = :spherical_spline,
    pre_calculated_levels = nothing,
)
    # Prepare data using the new internal function
    data = _prepare_ica_topo_data(ica, comp_idx, method, gridscale)

    # Use pre-calculated levels if provided, otherwise calculate levels
    if isnothing(pre_calculated_levels)
        levels = _calculate_topo_levels(data; use_global_scale = use_global_scale, num_levels = num_levels)
    else
        levels = pre_calculated_levels
    end

    # Plot using the new internal function
    co = _plot_ica_topo_on_axis!(
        topo_ax,
        fig,
        data,
        ica,
        levels;
        gridscale = gridscale,
        colormap = colormap,
        nan_color = nan_color,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
        method = method,
    )

    return co
end

function create_component_plots!(fig, state)

    # Number of plots
    if length(state.components) == size(state.component_data, 1)
        num_plots = state.n_visible_components
    else
        num_plots = length(state.components)  # Show all components in the subset
    end

    for i = 1:num_plots

        topo_ax = Axis(fig[i, 1])
        push!(state.topo_axs, topo_ax)

        # Get the actual component number
        comp_idx = _get_component_index(state, i)

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

        linkaxes!(ax, ax_channel)

        # Channel overlay plot
        lines!(
            ax_channel,
            @lift(state.dat.data.time[$(state.xrange)]),
            @lift(if $(state.show_channel)
                xrange = $(state.xrange)
                if first(xrange) >= 1 && last(xrange) <= length($(state.channel_data))
                    $(state.channel_data)[xrange] .* $(state.channel_yscale)
                else
                    zeros(Float64, length(xrange))
                end
            else
                zeros(Float64, length($(state.xrange)))
            end),
            color = :grey,
        )

        # Set initial x-axis limits
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
            _plot_ica_topo_in_viewer!(
                fig,
                topo_ax,
                state.ica_result,
                comp_idx;
                use_global_scale = state.use_global_scale[],
                gridscale = state.topo_gridscale,
                colormap = state.topo_colormap,
                num_levels = state.topo_num_levels,
                nan_color = state.topo_nan_color,
                head_kwargs = state.topo_head_kwargs,
                point_kwargs = state.topo_point_kwargs,
                label_kwargs = state.topo_label_kwargs,
                method = state.topo_method,
            )
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
            size(state.component_data, 1) - state.n_visible_components + 1,
            state.comp_start[] + state.n_visible_components,
        )
        state.comp_start[] = new_start
        update_components!(state)
    end

    # Connect apply button
    on(apply_button.clicks) do _
        text_value = text_input.displayed_string[]
        if !isempty(text_value)
            comps = parse_string_to_ints(text_value, size(state.component_data, 1))
            if !isempty(comps)
                @info "Creating new plot with components: $comps"
                current_channel = state.channel_data[]
                show_channel = state.show_channel[]
                use_global = state.use_global_scale[]
                invert = state.invert_scale[]
                new_fig =
                    plot_ica_component_activation(state.dat, state.ica_result, component_selection = components(comps))
            end
        end
    end
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
                visible_times = true_times[true_times .>= state.dat.data.time[first(
                    current_range,
                )].&&true_times .<= state.dat.data.time[last(current_range)]]

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
        if event.action in (Keyboard.press, Keyboard.repeat)
            if event.key == Keyboard.left || event.key == Keyboard.right
                # Handle x-axis scrolling
                current_range = state.xrange[]
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - state.window_size)
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                else  # right
                    new_start = min(
                        size(state.component_data, 2) - state.window_size + 1,
                        first(current_range) + state.window_size,
                    )
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
                    # For subsets, update all components in the subset
                    num_plots = if length(state.components) == size(state.component_data, 1)
                        state.n_visible_components
                    else
                        length(state.components)  # Update all components in the subset
                    end

                    for i = 1:num_plots
                        # Get the correct component index
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
                # Only handle page up/down if we're not using fixed components
                if length(state.components) == size(state.component_data, 1)
                    # Handle component scrolling
                    current_start = state.comp_start[]
                    if event.key == Keyboard.page_up
                        new_start = max(1, current_start - state.n_visible_components)
                    else  # page_down
                        new_start = min(
                            size(state.component_data, 1) - state.n_visible_components + 1,
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
    # Number of plots
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
            all_data[data_count] =
                _prepare_ica_topo_data(state.ica_result, comp_idx, state.topo_method, state.topo_gridscale)
        end
    end

    # Resize to actual count
    resize!(all_data, data_count)

    # Calculate levels using the simplified function
    levels_result = _calculate_ica_topo_levels(all_data, state.use_global_scale[], state.topo_num_levels)

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
                @warn "Trying to update non-existent lines_obs at index $i"
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
                    state.ica_result,
                    comp_idx;
                    use_global_scale = state.use_global_scale[],
                    gridscale = state.topo_gridscale,
                    colormap = state.topo_colormap,
                    num_levels = state.topo_num_levels,
                    nan_color = state.topo_nan_color,
                    method = state.topo_method,
                    head_kwargs = state.topo_head_kwargs,
                    point_kwargs = state.topo_point_kwargs,
                    label_kwargs = state.topo_label_kwargs,
                    pre_calculated_levels = level_to_use,
                )
                # Update title
                topo_ax.title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
            else
                @warn "Trying to update non-existent topo_axs at index $i"
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
                state.topo_axs[i].title = @sprintf("Invalid IC %d", comp_idx)
            end
            if i <= length(state.axs)
                state.axs[i].ylabel = ""
            end
        end
    end
end




# Internal function to prepare ICA topoplot data (interpolation only)
function _prepare_ica_topo_data(ica::InfoIca, comp_idx::Int, method::Symbol, gridscale::Int)
    if method == :spherical_spline
        return _data_interpolation_topo_spherical_spline(ica.mixing[:, comp_idx], ica.layout, gridscale)
    elseif method == :multiquadratic
        return _data_interpolation_topo_multiquadratic(ica.mixing[:, comp_idx], ica.layout, gridscale)
    else
        throw(ArgumentError("Unknown interpolation method: $method"))
    end
end

# Internal function to plot ICA topoplot on an axis
function _plot_ica_topo_on_axis!(
    topo_ax,
    fig,
    data::Matrix{Float64},
    ica::InfoIca,
    levels;
    gridscale::Int = 100,
    colormap::Symbol = :jet,
    nan_color::Symbol = :transparent,
    head_kwargs::Dict = Dict(),
    point_kwargs::Dict = Dict(),
    label_kwargs::Dict = Dict(),
    method::Symbol = :spherical_spline,
)
    # Clear the axis
    empty!(topo_ax)

    # Plot using the existing low-level function
    co = _plot_topo_on_axis!(
        topo_ax,
        fig,
        data,
        ica.layout,
        method,
        levels;
        gridscale = gridscale,
        colormap = colormap,
        nan_color = nan_color,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    hidedecorations!(topo_ax, grid = false)

    return co
end

# Internal function to calculate ICA topoplot levels (global or local)
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

# Helper function to get component index based on state
function _get_component_index(state, i::Int)
    if length(state.components) == size(state.component_data, 1)
        return state.comp_start[] + i - 1
    else
        return state.components[i]
    end
end







################################################################################




"""
    plot_eog_component_features(identified_comps::Dict, metrics_df::DataFrame; z_threshold::Float64=3.0)

Plot z-scores of EOG correlations from the metrics DataFrame and highlight identified components.

Uses the results from `identify_eye_components`.

# Arguments
- `identified_comps::Dict`: Dictionary returned by `identify_eye_components` (containing `:vertical_eye`, `:horizontal_eye`).
- `metrics_df::DataFrame`: DataFrame returned by `identify_eye_components`. Expected to have columns like `:vEOG_zscore`, `:hEOG_zscore` (based on the `v_eog_channel` and `h_eog_channel` arguments) and `:Component`.
- `vEOG_channel::Symbol`: Symbol of the vertical EOG channel used. Defaults to `:vEOG`. (Used to find column `Symbol("\$(v_eog_channel)_zscore")` in `metrics_df`).
- `hEOG_channel::Symbol`: Symbol of the horizontal EOG channel used. Defaults to `:hEOG`. (Used to find column `Symbol("\$(h_eog_channel)_zscore")` in `metrics_df`).
- `z_threshold::Float64`: Z-score threshold to draw lines on the plot (default: 3.0).

# Returns
- `fig::Figure`: The Makie Figure containing the z-score plots.
"""
function plot_eog_component_features(
    identified_comps::Dict,
    metrics_df::DataFrame;
    z_threshold::Float64 = 3.0,
    display_plot::Bool = true,
)

    # Extract data from inputs
    vEOG_corr_z = metrics_df.vEOG_zscore
    hEOG_corr_z = metrics_df.hEOG_zscore
    final_vEOG = identified_comps[:vEOG]
    final_hEOG = identified_comps[:hEOG]

    # Check if data is empty
    if isempty(metrics_df) || isempty(vEOG_corr_z) || isempty(hEOG_corr_z)
        @warn "Could not plot eye component features, input DataFrame or z-scores are empty."
        return Figure() # Return empty figure
    end

    n_components = nrow(metrics_df)

    # Plot vEOG/hEOG Correlation Z-Scores
    fig = Figure()
    ax_v = Axis(fig[1, 1], xlabel = "Component Number", ylabel = "Z-Score", title = "vEOG Correlation Z-Scores")
    # Use component indices from DataFrame for x-axis
    scatter!(ax_v, metrics_df.Component, vEOG_corr_z, color = :gray, markersize = 5)
    hlines!(ax_v, [z_threshold, -z_threshold], color = :gray, linestyle = :dash)
    # Highlight all identified components
    if !isempty(final_vEOG)
        # Get z-scores only for the identified components
        vEOG_z_scores_highlight = metrics_df[in.(metrics_df.Component, Ref(final_vEOG)), :vEOG_zscore]
        scatter!(ax_v, final_vEOG, vEOG_z_scores_highlight, color = :blue, markersize = 8)
        for comp_idx in final_vEOG # Annotate each identified component
            text!(
                ax_v,
                comp_idx,
                metrics_df[comp_idx, :vEOG_zscore],
                text = string(comp_idx),
                color = :blue,
                align = (:center, :bottom),
                fontsize = 10,
            )
        end
    end

    ax_h = Axis(fig[1, 2], xlabel = "Component Number", ylabel = "Z-Score", title = "hEOG Correlation Z-Scores")
    # Use component indices from DataFrame for x-axis
    scatter!(ax_h, metrics_df.Component, hEOG_corr_z, color = :gray, markersize = 5)
    hlines!(ax_h, [z_threshold, -z_threshold], color = :gray, linestyle = :dash)
    # Highlight all identified components
    if !isempty(final_hEOG)
        # Get z-scores only for the identified components
        hEOG_z_scores_highlight = metrics_df[in.(metrics_df.Component, Ref(final_hEOG)), :hEOG_zscore]
        scatter!(ax_h, final_hEOG, hEOG_z_scores_highlight, color = :blue, markersize = 8)
        for comp_idx in final_hEOG # Annotate each identified component
            text!(
                ax_h,
                comp_idx,
                metrics_df[comp_idx, :hEOG_zscore],
                text = string(comp_idx),
                color = :blue,
                align = (:center, :bottom),
                fontsize = 10,
            )
        end
    end

    if display_plot
        display_figure(fig)
    end

    return fig, (ax_v, ax_h)

end

# ECG Component Identification 
"""
    _find_peaks(data::AbstractVector; min_prominence_std::Real=2.0, window_size::Int=1)

Enhanced peak finder. Finds indices where data point is greater than its
neighbors within a window (positive peaks) OR less than its neighbors within a window (negative peaks)
and exceeds the threshold. Returns indices of peaks.

# Arguments
- `data::AbstractVector`: Input data vector
- `min_prominence_std::Real`: Minimum prominence threshold in standard deviations (default: 2.0)
- `window_size::Int`: Number of samples to look on each side for comparison (default: 1)
"""
function _find_peaks(data::AbstractVector; min_prominence_std::Real = 2.0, window_size::Int = 1)
    if length(data) < 2 * window_size + 1
        return Int[]
    end

    peaks = Int[]
    mean_val = mean(data)
    std_val = std(data)
    # Handle zero std deviation case
    threshold_pos = (std_val ≈ 0) ? mean_val : mean_val + min_prominence_std * std_val
    threshold_neg = (std_val ≈ 0) ? mean_val : mean_val - min_prominence_std * std_val

    for i = (window_size+1):(length(data)-window_size)
        # Get the window around the current point
        left_window = data[(i-window_size):(i-1)]
        right_window = data[(i+1):(i+window_size)]

        # Positive peaks (greater than all neighbors in window and above threshold)
        if all(data[i] .> left_window) && all(data[i] .> right_window) && data[i] > threshold_pos
            push!(peaks, i)
            # Negative peaks (less than all neighbors in window and below threshold)
        elseif all(data[i] .< left_window) && all(data[i] .< right_window) && data[i] < threshold_neg
            push!(peaks, i)
        end
    end
    return peaks
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
    kurtosis_comps::Vector{Int},
    metrics_df::DataFrame,
    z_threshold::Float64 = 3.0,
    display_plot::Bool = true,
)

    # Create figure
    fig = Figure()

    # Plot spatial kurtosis z-scores
    ax = Axis(
        fig[1, 1],
        xlabel = "Component",
        ylabel = "Spatial Kurtosis Z-Score",
        title = "Component Spatial Kurtosis Z-Scores",
    )

    # Plot all components
    scatter!(ax, metrics_df.Component, metrics_df.SpatialKurtosisZScore, color = :gray)

    # Highlight high kurtosis components
    if !isempty(kurtosis_comps)
        kurtosis_values = metrics_df[in.(metrics_df.Component, Ref(kurtosis_comps)), :SpatialKurtosisZScore]
        scatter!(ax, kurtosis_comps, kurtosis_values, color = :red, markersize = 8)

        # Add labels for high kurtosis components
        for (i, comp) in enumerate(kurtosis_comps)
            text!(
                ax,
                comp,
                kurtosis_values[i],
                text = string(comp),
                color = :red,
                align = (:center, :bottom),
                fontsize = 10,
            )
        end
    end

    # Add reference lines
    hlines!(ax, [z_threshold], color = :red, linestyle = :dash)
    hlines!(ax, [0.0], color = :gray, linestyle = :dot)

    if display_plot
        display_figure(fig)
    end

    return fig, ax
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
function plot_ecg_component_features_(
    identified_comps::Vector{Int},
    metrics_df::DataFrame;
    min_bpm::Real = 40,
    max_bpm::Real = 120,
    max_ibi_std_s::Real = 0.2,
    min_peak_ratio::Real = 0.7,
    display_plot::Bool = true,
)
    # Create figure with two panels
    fig = Figure()

    # Calculate heart rates
    heart_rates = [isnan(ibi) || ibi <= 0 ? NaN : 60.0/ibi for ibi in metrics_df.mean_ibi_s]
    metrics_df[!, :heart_rate_bpm] = heart_rates

    # Left panel: Heart Rate vs Peak Ratio
    ax1 = Axis(fig[1, 1], xlabel = "Heart Rate (BPM)", ylabel = "Peak Ratio (valid/total)")
    # Right panel: Heart Rate vs IBI Regularity (std)
    ax2 = Axis(fig[1, 2], xlabel = "Heart Rate (BPM)", ylabel = "IBI Std Dev (seconds)")

    # Plot non-ECG components
    non_ecg_idx = setdiff(1:nrow(metrics_df), identified_comps)
    non_ecg_df = metrics_df[non_ecg_idx, :]

    # Filter out NaNs for plotting
    valid_non_ecg = findall(.!isnan.(non_ecg_df.heart_rate_bpm) .& .!isnan.(non_ecg_df.peak_ratio))
    if !isempty(valid_non_ecg)
        scatter!(
            ax1,
            non_ecg_df.heart_rate_bpm[valid_non_ecg],
            non_ecg_df.peak_ratio[valid_non_ecg],
            color = :gray,
            markersize = 10,
        )
    end

    # Filter valid points for second plot
    valid_non_ecg2 = findall(.!isnan.(non_ecg_df.heart_rate_bpm) .& .!isnan.(non_ecg_df.std_ibi_s))
    if !isempty(valid_non_ecg2)
        scatter!(
            ax2,
            non_ecg_df.heart_rate_bpm[valid_non_ecg2],
            non_ecg_df.std_ibi_s[valid_non_ecg2],
            color = :gray,
            markersize = 10,
        )
    end

    # Plot ECG components
    ecg_df = metrics_df[in.(metrics_df.Component, Ref(identified_comps)), :]

    # Filter out NaNs
    valid_ecg = findall(.!isnan.(ecg_df.heart_rate_bpm) .& .!isnan.(ecg_df.peak_ratio))
    if !isempty(valid_ecg)
        scatter!(ax1, ecg_df.heart_rate_bpm[valid_ecg], ecg_df.peak_ratio[valid_ecg], color = :black, markersize = 16)

        # Add component labels
        for i in valid_ecg
            text!(
                ax1,
                ecg_df.heart_rate_bpm[i],
                ecg_df.peak_ratio[i],
                text = string(ecg_df.Component[i]),
                align = (:center, :bottom),
                offset = (0, 3),
                fontsize = 12,
            )
        end
    end

    # Plot ECG components in second panel
    valid_ecg2 = findall(.!isnan.(ecg_df.heart_rate_bpm) .& .!isnan.(ecg_df.std_ibi_s))
    if !isempty(valid_ecg2)
        scatter!(ax2, ecg_df.heart_rate_bpm[valid_ecg2], ecg_df.std_ibi_s[valid_ecg2], color = :black, markersize = 16)

        # Add component labels
        for i in valid_ecg2
            text!(
                ax2,
                ecg_df.heart_rate_bpm[i],
                ecg_df.std_ibi_s[i],
                text = string(ecg_df.Component[i]),
                align = (:center, :bottom),
                offset = (0, 3),
                fontsize = 12,
            )
        end
    end

    # Add reference ranges for normal heart boundary and selection criterion
    vlines!(ax1, [min_bpm, max_bpm], color = (:black, 0.5), linestyle = :dash)
    vlines!(ax2, [min_bpm, max_bpm], color = (:black, 0.5), linestyle = :dash)
    hlines!(ax1, [min_peak_ratio], color = (:black, 0.5), linestyle = :dash)
    hlines!(ax2, [max_ibi_std_s], color = (:black, 0.5), linestyle = :dash)

    if display_plot
        display_figure(fig)
    end

    return fig, (ax1, ax2)

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
    line_noise_comps::Vector{Int},
    metrics_df::DataFrame,
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    z_threshold::Float64 = 3.0,
    min_harmonic_power::Real = 0.5,
    display_plot::Bool = true,
)

    # Create figure with two subplots
    fig = Figure(size = (1000, 400))

    # Plot 1: Power Ratio Z-Scores
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Component",
        ylabel = "Power Ratio Z-Score",
        title = "Line Frequency Power Ratio Z-Scores",
    )

    # Plot all components with label
    scatter!(ax1, metrics_df.Component, metrics_df.PowerRatioZScore, color = :gray, label = "All Components")

    # Highlight identified components with label
    if !isempty(line_noise_comps)
        identified_metrics = metrics_df[in.(metrics_df.Component, Ref(line_noise_comps)), :]
        scatter!(
            ax1,
            identified_metrics.Component,
            identified_metrics.PowerRatioZScore,
            color = :red,
            markersize = 8,
            label = "Line Noise Components",
        )

        # Add component numbers as labels
        for (i, comp) in enumerate(line_noise_comps)
            row = metrics_df[metrics_df.Component .== comp, :]
            text!(
                ax1,
                comp,
                row.PowerRatioZScore[1],
                text = string(comp),
                color = :red,
                align = (:center, :bottom),
                fontsize = 10,
            )
        end
    end

    # Add threshold line with label
    hlines!(ax1, [z_threshold], color = :red, linestyle = :dash, label = "Threshold")

    # Add legend
    axislegend(ax1, position = (1.0, 1.0))

    # Plot 2: Harmonic Ratios
    ax2 = Axis(fig[1, 2], xlabel = "Component", ylabel = "Power Ratio", title = "Harmonic Power Ratios")

    # Plot harmonic ratios with labels
    scatter!(ax2, metrics_df.Component, metrics_df.Harmonic2Ratio, color = :blue, label = "2nd Harmonic")
    scatter!(ax2, metrics_df.Component, metrics_df.Harmonic3Ratio, color = :green, label = "3rd Harmonic")

    # Add reference line for minimum harmonic power with label
    hlines!(ax2, [min_harmonic_power], color = :gray, linestyle = :dash, label = "Min Harmonic Power")

    # Highlight identified components with label
    if !isempty(line_noise_comps)
        identified_metrics = metrics_df[in.(metrics_df.Component, Ref(line_noise_comps)), :]
        scatter!(
            ax2,
            identified_metrics.Component,
            identified_metrics.Harmonic2Ratio,
            color = :red,
            markersize = 8,
            label = "Line Noise Components",
        )
        scatter!(ax2, identified_metrics.Component, identified_metrics.Harmonic3Ratio, color = :red, markersize = 8)
    end

    # Add legend
    axislegend(ax2, position = (1.0, 1.0))

    if display_plot
        display_figure(fig)
    end

    return fig, (ax1, ax2)
end


"""
    plot_ecg_component_features(identified_comps::Vector{Int64}, metrics_df::DataFrame)

Create a simplified visualization of ECG component detection metrics.

# Arguments
- `identified_comps::Vector{Int64}`: Vector of component indices identified as ECG artifacts
- `metrics_df::DataFrame`: DataFrame with component metrics

# Returns
- `fig::Figure`: The Makie Figure containing the plot
"""
function plot_ecg_component_features(identified_comps::Vector{Int64}, metrics_df::DataFrame)
    # Create figure with two panels
    fig = Figure(size = (1000, 600))

    # Calculate heart rates
    heart_rates = [isnan(ibi) || ibi <= 0 ? NaN : 60.0/ibi for ibi in metrics_df.mean_ibi_s]
    metrics_df[!, :heart_rate_bpm] = heart_rates

    # Left panel: Heart Rate vs Peak Ratio
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Heart Rate (BPM)",
        ylabel = "Peak Ratio (valid/total)",
        title = "ECG Detection Metrics",
    )

    # Right panel: Heart Rate vs IBI Regularity (std)
    ax2 =
        Axis(fig[1, 2], xlabel = "Heart Rate (BPM)", ylabel = "IBI Std Dev (seconds)", title = "Heart Rate Regularity")

    # Plot non-ECG components
    non_ecg_idx = setdiff(1:nrow(metrics_df), identified_comps)
    non_ecg_df = metrics_df[non_ecg_idx, :]

    # Filter out NaNs for plotting
    valid_non_ecg = findall(.!isnan.(non_ecg_df.heart_rate_bpm) .& .!isnan.(non_ecg_df.peak_ratio))
    if !isempty(valid_non_ecg)
        scatter!(
            ax1,
            non_ecg_df.heart_rate_bpm[valid_non_ecg],
            non_ecg_df.peak_ratio[valid_non_ecg],
            color = :gray,
            markersize = 8,
            label = "Non-ECG",
        )
    end

    # Filter valid points for second plot
    valid_non_ecg2 = findall(.!isnan.(non_ecg_df.heart_rate_bpm) .& .!isnan.(non_ecg_df.std_ibi_s))
    if !isempty(valid_non_ecg2)
        scatter!(
            ax2,
            non_ecg_df.heart_rate_bpm[valid_non_ecg2],
            non_ecg_df.std_ibi_s[valid_non_ecg2],
            color = :gray,
            markersize = 8,
            label = "Non-ECG",
        )
    end

    # Plot ECG components
    ecg_df = metrics_df[in.(metrics_df.Component, Ref(identified_comps)), :]

    # Filter out NaNs
    valid_ecg = findall(.!isnan.(ecg_df.heart_rate_bpm) .& .!isnan.(ecg_df.peak_ratio))
    if !isempty(valid_ecg)
        ecg_scatter1 = scatter!(
            ax1,
            ecg_df.heart_rate_bpm[valid_ecg],
            ecg_df.peak_ratio[valid_ecg],
            color = :red,
            markersize = 12,
            marker = :diamond,
            label = "ECG",
        )

        # Add component labels
        for i in valid_ecg
            text!(
                ax1,
                ecg_df.heart_rate_bpm[i],
                ecg_df.peak_ratio[i],
                text = string(ecg_df.Component[i]),
                align = (:center, :bottom),
                offset = (0, 3),
                fontsize = 12,
            )
        end
    end

    # Plot ECG components in second panel
    valid_ecg2 = findall(.!isnan.(ecg_df.heart_rate_bpm) .& .!isnan.(ecg_df.std_ibi_s))
    if !isempty(valid_ecg2)
        scatter!(
            ax2,
            ecg_df.heart_rate_bpm[valid_ecg2],
            ecg_df.std_ibi_s[valid_ecg2],
            color = :red,
            markersize = 12,
            marker = :diamond,
            label = "ECG",
        )

        # Add component labels
        for i in valid_ecg2
            text!(
                ax2,
                ecg_df.heart_rate_bpm[i],
                ecg_df.std_ibi_s[i],
                text = string(ecg_df.Component[i]),
                align = (:center, :bottom),
                offset = (0, 3),
                fontsize = 12,
            )
        end
    end

    # Add reference ranges using pre-defined values 
    # Normal heart rate range (typical values)
    vlines!(ax1, [60, 100], color = (:green, 0.5), linestyle = :dash, label = "Normal HR Range")
    vlines!(ax2, [60, 100], color = (:green, 0.5), linestyle = :dash)

    # Get threshold values from actual data when possible
    min_peak_ratio = 0.7  # Default if no components found
    max_std = 0.12        # Default if no components found

    if any(ecg_df.is_ecg_artifact)
        # Get actual values from data
        min_peak_ratio = minimum(ecg_df.peak_ratio[ecg_df.is_ecg_artifact])
        max_std = maximum(ecg_df.std_ibi_s[ecg_df.is_ecg_artifact])
    end

    # Add threshold lines
    hlines!(ax1, [min_peak_ratio], color = (:red, 0.5), linestyle = :dash, label = "Min Peak Ratio")
    hlines!(ax2, [max_std], color = (:red, 0.5), linestyle = :dash, label = "Max StdDev")

    # Add legends
    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rt)

    # Calculate BPM range from actual data
    min_hr = 40   # Default minimum heart rate
    max_hr = 120  # Default maximum heart rate
    if !isempty(ecg_df) && any(.!isnan.(ecg_df.heart_rate_bpm))
        # Use the actual range from identified components
        valid_hrs = filter(!isnan, ecg_df.heart_rate_bpm)
        if !isempty(valid_hrs)
            min_hr = floor(Int, minimum(valid_hrs))
            max_hr = ceil(Int, maximum(valid_hrs))
        end
    end

    # Display summary text
    Label(
        fig[2, 1:2],
        "Found $(length(identified_comps)) ECG components: $(join(identified_comps, ", "))\n" *
        "Criteria: Heart rate $min_hr-$max_hr BPM, StdDev ≤ $(round(max_std, digits=3))s, Peak Ratio ≥ $(round(min_peak_ratio, digits=2))",
        fontsize = 14,
        tellwidth = false,
    )

    return fig
end
