# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ERP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    :figure_title => ("ERP Plot", "Title for the plot window"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Whether to show the title"),

    # Line styling
    :linewidth => (2, "Line width for ERP traces"),
    :color => (:black, "Color for ERP traces (can be a single color or a vector of colors, one per dataset)"),
    :linestyle =>
        (:solid, "Line style for ERP traces (can be a single style or a vector of styles, one per dataset)"),
    :colormap => (:jet, "Colormap for multi-channel plots"),

    # Plot configuration
    :dims => (nothing, "Grid dimensions as (rows, cols). If nothing, automatically determined"),
    :yreversed => (false, "Whether to reverse the y-axis"),
    :average_channels => (false, "Whether to average across channels"),
    :interactive => (true, "Whether to enable interactive features"),

    # Legend parameters - get all Legend attributes with their actual defaults
    # This allows users to control any Legend parameter
    [
        Symbol("legend_$(attr)") => (get(LEGEND_DEFAULTS, attr, nothing), "Legend $(attr) parameter") for
        attr in propertynames(Legend)
    ]...,

    # Override specific legend parameters with custom defaults
    :legend => (true, "Whether to show the legend"),
    :legend_label => ("", "Title for the legend"),
    :legend_framevisible => (true, "Whether to show the frame of the legend"),
    :legend_position => (
        :lt,
        "Position of the legend for axislegend() (symbol like :lt, :rt, :lb, :rb, or tuple like (:left, :top), or (0.5, 0.5))",
    ),
    :legend_channel => ([], "If plotting multiple plots, within channel to put the legend on."),
    :legend_labels => ([], "If plotting multiple plots, within channel to put the legend on."),
    :legend_nbanks => (nothing, "Number of columns for the legend. If nothing, automatically determined."),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Origin lines
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),
)

"""
    plot_erp(filepath::String; 
             layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
             condition_selection::Function = conditions(),
             channel_selection::Function = channels(),
             sample_selection::Function = samples(),
             baseline_interval::BaselineInterval = nothing,
             kwargs...)

Load ERP data from a JLD2 file and create plots.

# Arguments
- `filepath::String`: Path to JLD2 file containing ErpData or Vector{ErpData}
- `layout`: Layout specification (see main plot_erp documentation)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering
- `sample_selection::Function`: Function that returns boolean vector for sample filtering
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Load and plot from file
plot_erp("grand_average_erps_good.jld2")

# With channel selection
plot_erp("grand_average_erps_good.jld2", channel_selection = channels([:PO7, :PO8]), layout = :grid)
```
"""
function plot_erp(
    filepath::String;
    layout::Union{Symbol,PlotLayout,Vector{Int}} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_interval::BaselineInterval = nothing,
    kwargs...,
)
    # Load data from file
    data = load_data(filepath)
    if isnothing(data)
        @minimal_error_throw "No data found in file: $filepath"
    end

    # Dispatch will handle ErpData vs Vector{ErpData} automatically
    return plot_erp(
        data;
        layout = layout,
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

"""
    plot_erp(dat::ErpData; 
             layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
             channel_selection::Function = channels(),
             sample_selection::Function = samples(),
             baseline_interval::BaselineInterval = nothing,
             kwargs...)

Create ERP plots with flexible layout options.

# Arguments
- `dat::ErpData`: ERP data structure
- `layout`: Layout specification:
  - `:single` (default): Single plot with all channels
  - `:grid`: Auto-calculated grid layout
  - `:topo`: Topographic layout based on channel positions
  - `PlotLayout`: Custom layout object
  - `Vector{Int}`: Custom grid dimensions [rows, cols]
  
- `channel_selection::Function`: Function that returns boolean vector for channel filtering
- `sample_selection::Function`: Function that returns boolean vector for sample filtering
- `baseline_interval::BaselineInterval`: Baseline correction interval. Can be `nothing` (no baseline), tuple like `(-0.2, 0.0)`, `IntervalTime`, or `IntervalIndex`. Default `nothing` means no baseline correction.
- `kwargs`: Additional keyword arguments

$(generate_kwargs_doc(PLOT_ERP_KWARGS))

# Examples
```julia
# Single plot with all channels
plot_erp(dat)

# Grid layout
plot_erp(dat, layout = :grid)

# Custom grid dimensions
plot_erp(dat, layout = [3, 4])

# Topographic layout
plot_erp(dat, layout = :topo)

# Custom layout object
layout = create_grid_layout(channels(dat), rows = 2, cols = 3)
plot_erp(dat, layout = layout)

# Disable interactivity
plot_erp(dat, interactive = false)
```

# Interactive Controls
When `interactive = true` (default):
- **Up Arrow**: Zoom in on Y-axis (compress Y limits)
- **Down Arrow**: Zoom out on Y-axis (expand Y limits)
- **Left Arrow**: Zoom in on X-axis (compress time range)
- **Right Arrow**: Zoom out on X-axis (expand time range)

# Mouse Selection
- **Shift + Left Click + Drag**: Select a time region (blue rectangle)
- **Right Click on Selection**: Open context menu with plot options
  - Topoplot (multiquadratic)
  - Topoplot (spherical_spline)
"""
function plot_erp(
    dat::ErpData;
    layout::Union{Symbol,PlotLayout,Vector{Int}} = :single,
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_interval::BaselineInterval = nothing,
    kwargs...,
)
    # For single ErpData, condition_selection doesn't apply (there's only one condition)
    return plot_erp(
        [dat];
        layout = layout,
        condition_selection = conditions(),  # Always select all (just the one condition)
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

"""
    plot_erp(datasets::Vector{ErpData}; 
             layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
             condition_selection::Function = conditions(),
             channel_selection::Function = channels(), 
             sample_selection::Function = samples(),
             baseline_interval::BaselineInterval = nothing,
             kwargs...)

Plot multiple ERP datasets on the same axis (e.g., conditions).
"""
function plot_erp(
    datasets::Vector{ErpData};
    layout::Union{Symbol,PlotLayout,Vector{Int}} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_interval::BaselineInterval = nothing,
    kwargs...,
)

    # Prepare kwargs and data
    plot_kwargs, _ = _prepare_plot_kwargs(kwargs)
    dat_subset, all_plot_channels, original_channels = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        baseline_interval = baseline_interval,
    )

    # set default plot title only for single layouts
    # For grid/topo layouts, we want individual channel names, not a global title
    if plot_kwargs[:show_title] && plot_kwargs[:title] == "" && layout == :single
        plot_kwargs[:title] =
            length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"
        if plot_kwargs[:average_channels]
            plot_kwargs[:title] = "Avg: $(print_vector(original_channels))"
        end
    end

    # Generate window title from datasets
    title_str = _generate_window_title(dat_subset)
    set_window_title(title_str)

    # Create figure and apply layout system
    fig = Figure(title = plot_kwargs[:figure_title])
    plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout)
    axes, channels = apply_layout!(fig, plot_layout; plot_kwargs...)

    # Store line references for control panel (if interactive)
    # Structure: line_refs[ax_idx][dataset_idx][channel_idx] = line
    line_refs = plot_kwargs[:interactive] ? [Dict{Int,Dict{Symbol,Any}}() for _ in axes] : nothing

    # Store legend references for linked interactions (if interactive)
    legend_refs = plot_kwargs[:interactive] ? Vector{Union{Legend,Nothing}}(undef, length(axes)) : nothing

    # Now do the actual plotting for each axis
    for (ax_idx, (ax, channel)) in enumerate(zip(axes, channels))
        channels_to_plot = plot_layout.type == :single ? all_plot_channels : [channel]
        @info "plot_erp ($layout): $(print_vector(channels_to_plot))"
        ax_line_refs = plot_kwargs[:interactive] ? line_refs[ax_idx] : nothing
        ax_result, leg = _plot_erp!(ax, dat_subset, channels_to_plot; line_refs = ax_line_refs, plot_kwargs...)
        if plot_kwargs[:interactive] && legend_refs !== nothing
            legend_refs[ax_idx] = leg
        end
    end

    # Apply our axis stuff
    _apply_axis_properties!.(axes; plot_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...) # slightly different for grid and topo layouts

    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)

    # Add keyboard interactivity if enabled
    if plot_kwargs[:interactive]
        _setup_shared_interactivity!(fig, axes, :erp)

        # Disable default interactions that conflict with our custom selection (all axes)
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end

        # Set up selection system for all axes (will work with linked axes)
        selection_state = SharedSelectionState(axes)

        # Set up selection system that works for all layouts
        _setup_unified_selection!(fig, axes, selection_state, datasets, plot_layout, _handle_erp_right_click!)

        # Set up channel selection events for topo and grid layouts
        if plot_layout.type in (:topo, :grid)
            _setup_channel_selection_events!(fig, selection_state, plot_layout, datasets, axes, plot_layout.type)
        end

        # Set up control panel (press 'c' to open)
        _setup_erp_control_panel!(fig, dat_subset, axes, baseline_interval, line_refs)

    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    # reset default title
    set_window_title("Makie")
    return fig, axes
end


"""
    plot_erp!(fig::Figure, ax::Axis, dat::ErpData; kwargs...)

Plot ERP data on an existing axis, mutating the figure and axis.

# Arguments
- `fig::Figure`: The figure to plot on
- `ax::Axis`: The axis to plot on  
- `dat::ErpData`: The ERP data to plot
- `kwargs...`: Additional plotting arguments (see PLOT_ERP_KWARGS)

# Returns
- `ax::Axis`: The axis that was plotted on
"""
plot_erp!(fig::Figure, ax::Axis, dat::ErpData; kwargs...) = plot_erp!(fig, ax, [dat]; kwargs...)

"""
    plot_erp!(fig::Figure, ax::Axis, datasets::Vector{ErpData}; kwargs...)

Plot multiple ERP datasets on an existing axis, mutating the figure and axis.

# Arguments
- `fig::Figure`: The figure to plot on
- `ax::Axis`: The axis to plot on
- `datasets::Vector{ErpData}`: The ERP datasets to plot
- `kwargs...`: Additional plotting arguments (see PLOT_ERP_KWARGS)

# Returns
- `ax::Axis`: The axis that was plotted on
"""
function plot_erp!(fig::Figure, ax::Axis, datasets::Vector{ErpData}; kwargs...)
    plot_kwargs, _ = _prepare_plot_kwargs(kwargs)
    baseline_interval = get(kwargs, :baseline_interval, nothing)
    dat_subset, all_plot_channels, _ = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = conditions(),
        channel_selection = get(plot_kwargs, :channel_selection, channels()),
        sample_selection = get(plot_kwargs, :sample_selection, samples()),
        baseline_interval = baseline_interval,
    )
    _plot_erp!(ax, dat_subset, all_plot_channels; plot_kwargs...)
    return ax
end


"""
    _plot_erp!(ax::Axis, datasets::Vector{ErpData}, channels::Vector{Symbol}; kwargs...)

Internal function to plot ERP data on an axis.
Handles both single and multiple datasets.
Note: datasets should already be subset based on channel_selection and sample_selection.

# Keyword Arguments
- `condition_mask::Vector{Bool}`: Mask to set visibility of each dataset. 
  `condition_mask[i]` controls visibility of `datasets[i]`. Defaults to all `true`.
- `line_refs::Union{Dict, Nothing}`: Dictionary to store line references for interactive updates.
  If provided, stores lines as line_refs[dataset_idx][channel] = line.
"""
function _plot_erp!(
    ax::Axis,
    datasets::Vector{ErpData},
    channels::Vector{Symbol};
    condition_mask::Vector{Bool} = Bool[],
    line_refs = nothing,
    kwargs...,
)

    kwargs = merge(PLOT_ERP_KWARGS, kwargs)

    # Default condition_mask to all true if not provided
    if isempty(condition_mask)
        condition_mask = fill(true, length(datasets))
    end

    # Compute colors and linestyles for each dataset
    all_colors = _compute_dataset_colors(
        kwargs[:color],
        length(datasets),
        length(channels),
        kwargs[:colormap],
        kwargs[:_color_explicitly_set],
    )
    all_linestyles = _compute_dataset_linestyles(kwargs[:linestyle], length(datasets))

    # Plot each dataset for ALL channels in this subplot
    for (dataset_idx, dat) in enumerate(datasets)

        for (channel_idx, channel) in enumerate(channels)

            # axis label
            label = isempty(kwargs[:legend_labels]) ? dat.condition_name : kwargs[:legend_labels][dataset_idx]
            if length(channels) > 1 # More than one channel in this subplot
                label *= " ($channel)"
            end

            color_idx = (dataset_idx - 1) * length(channels) + channel_idx

            # Set visibility based on condition_mask
            condition_visible = dataset_idx <= length(condition_mask) ? condition_mask[dataset_idx] : true

            # Always use Observable for y-data (allows updates for baseline changes and linked legend interactions)
            y_obs = Observable(dat.data[!, channel])
            line = lines!(
                ax,
                dat.data[!, :time],
                y_obs,
                linewidth = kwargs[:linewidth],
                color = all_colors[color_idx],
                linestyle = all_linestyles[dataset_idx],
                label = label,
                visible = condition_visible,
            )

            # Store line and y Observable if references are requested
            if line_refs !== nothing
                if !haskey(line_refs, dataset_idx)
                    line_refs[dataset_idx] = Dict{Symbol,Tuple{Any,Observable}}()
                end
                line_refs[dataset_idx][channel] = (line, y_obs)
            end
        end
    end

    _set_origin_lines!(ax; add_xy_origin = kwargs[:add_xy_origin])
    leg = _add_legend!(ax, channels, datasets, kwargs)

    return ax, leg  # Return both axis and legend
end


"""
    _prepare_plot_kwargs(kwargs)

Prepare plot kwargs by merging with defaults and tracking color setting.
Returns (plot_kwargs, color_explicitly_set).
"""
function _prepare_plot_kwargs(kwargs)
    color_explicitly_set = haskey(kwargs, :color)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)
    plot_kwargs[:_color_explicitly_set] = color_explicitly_set
    return plot_kwargs, color_explicitly_set
end

"""
    _prepare_erp_data(datasets, plot_kwargs; condition_selection, channel_selection, sample_selection)

Prepare ERP data for plotting: subset, average channels if needed, extract channel labels.
Returns (dat_subset, all_plot_channels, original_channels).
"""
function _prepare_erp_data(
    datasets::Vector{ErpData},
    plot_kwargs;
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    baseline_interval::BaselineInterval = nothing,
)
    # Data subsetting
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = true,
    )

    # Check if subsetting resulted in empty data
    if isempty(dat_subset)
        n_conditions = length(datasets)
        if n_conditions == 0
            @minimal_error_throw "No data available (empty dataset)"
        else
            @minimal_error_throw "No data matched the selection criteria. Available condition indices: 1:$n_conditions"
        end
    end

    # Apply baseline correction if requested
    if baseline_interval !== nothing
        @info "Applying baseline correction to $(length(dat_subset)) datasets"
        baseline!.(dat_subset, Ref(baseline_interval))
    end

    # Channel averaging if requested
    original_channels = nothing
    if plot_kwargs[:average_channels]
        @info "Averaging channels for $(length(dat_subset)) datasets"
        original_channels = channel_labels(dat_subset)
        dat_subset = channel_average(dat_subset, reduce = true)
    end

    # Extract channel labels
    selected_channels = channel_labels(dat_subset)
    extra_channels = extra_labels(dat_subset)
    all_plot_channels = vcat(selected_channels, extra_channels)

    return dat_subset, all_plot_channels, original_channels
end


function _handle_erp_right_click!(selection_state, mouse_x, data)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_erp_context_menu!(selection_state, data)
    end
end

function _show_erp_context_menu!(selection_state, data)

    menu_fig = Figure()
    plot_types = ["Topoplot (multiquadratic)", "Topoplot (spherical_spline)"]

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            original_data, x_min, x_max = _get_erp_selection_bounds(selection_state, data)

            # Create time-based sample selection for the topo plot
            time_sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max)

            if btn.label[] == "Topoplot (multiquadratic)"
                plot_topography(original_data, sample_selection = time_sample_selection, method = :multiquadratic)
            elseif btn.label[] == "Topoplot (spherical_spline)"
                plot_topography(original_data, sample_selection = time_sample_selection, method = :spherical_spline)
            end
        end
    end

    new_screen = getfield(Main, :GLMakie).Screen()
    display(new_screen, menu_fig)
end

"""
    _get_erp_selection_bounds(selection_state, data)

Extract time bounds from selection state and return original data with bounds.
Does not subset the data - preserves all electrodes for topo plots.
Returns (data, x_min, x_max).
"""
function _get_erp_selection_bounds(selection_state, data)
    x_min, x_max = minmax(selection_state.bounds[]...)
    return (data, x_min, x_max)
end

"""
    _compute_dataset_colors(color_val, n_datasets, n_channels, colormap, color_explicitly_set)

Compute colors for each dataset-channel combination.
Returns a vector of colors with length n_datasets * n_channels.
Colors cycle across all channel-dataset combinations.
"""
function _compute_dataset_colors(color_val, n_datasets::Int, n_channels::Int, colormap, color_explicitly_set::Bool)
    n_total = n_datasets * n_channels
    if color_val isa Vector # User specified colors - cycle through them for all channel-dataset combinations
        return [color_val[(i-1)%length(color_val)+1] for i = 1:n_total]
    elseif n_total > 1 && !color_explicitly_set # Default: use colormap for all channel-dataset combinations
        return Makie.cgrad(colormap, n_total, categorical = true)
    else # Single color - use for all channel-dataset combinations
        return [color_val for _ = 1:n_total]
    end
end

"""
    _compute_dataset_linestyles(linestyle_val, n_datasets)

Compute linestyles for each dataset based on user input.
Returns a vector of linestyles, one per dataset.
"""
function _compute_dataset_linestyles(linestyle_val, n_datasets::Int)
    if linestyle_val isa Vector # User specified linestyles per dataset - wrap if needed
        return [linestyle_val[(i-1)%length(linestyle_val)+1] for i = 1:n_datasets]
    else # Single linestyle - use for all datasets
        return [linestyle_val for _ = 1:n_datasets]
    end
end


"""
    _add_legend!(ax::Axis, channels::Vector{Symbol}, datasets::Vector{ErpData}, kwargs::Dict)

Add legend to axis if conditions are met.
Handles legend_channel filtering, position, and other legend attributes.
"""
function _add_legend!(ax::Axis, channels::Vector{Symbol}, datasets::Vector{ErpData}, kwargs::Dict)

    # Check if legend should be shown
    if !kwargs[:legend] || (length(channels) == 1 && length(datasets) == 1)
        return nothing
    end
    if !isempty(kwargs[:legend_channel]) && isempty(intersect(kwargs[:legend_channel], channels))
        return nothing
    end

    # Extract legend parameters
    legend_label = kwargs[:legend_label]
    legend_position = kwargs[:legend_position]
    if kwargs[:legend_nbanks] === nothing
        kwargs[:legend_nbanks] = length(channels) > 10 ? cld(length(channels), 10) : 1
    end
    legend_kwargs = _extract_legend_kwargs(kwargs)

    # Add legend with position and optional label
    if legend_label != ""
        leg = axislegend(ax, legend_label; position = legend_position, legend_kwargs...)
    else
        leg = axislegend(ax; position = legend_position, legend_kwargs...)
    end

    return leg
end

"""
    _setup_linked_legend_interactions!(line_refs::Vector{<:Dict})

Set up linked legend interactions so clicking a legend entry in one plot
toggles visibility of the corresponding condition in all plots.
"""
function _setup_linked_legend_interactions!(line_refs::Vector{<:Dict})
    # Create a mapping: dataset_idx -> all lines across all axes for that dataset
    dataset_lines = Dict{Int,Vector{Any}}()

    # Collect all lines for each dataset
    for ax_line_refs in line_refs
        for (dataset_idx, channel_lines) in ax_line_refs
            lines = get!(dataset_lines, dataset_idx, Any[])
            append!(lines, [line_data[1] for line_data in values(channel_lines)])
        end
    end

    # When any line's visibility changes, update all other lines for that dataset
    for lines in values(dataset_lines)
        length(lines) > 1 || continue  # Only need syncing if there are multiple lines

        # Create a flag to prevent infinite loops
        syncing = Ref(false)

        for line in lines
            other_lines = [l for l in lines if l !== line]
            on(line.visible) do visible_val
                syncing[] && return  # Skip if already syncing
                syncing[] = true
                for other_line in other_lines
                    other_line.visible = visible_val
                end
                syncing[] = false
            end
        end
    end

end

"""
    _setup_erp_control_panel!(fig::Figure, dat_subset::Vector{ErpData}, axes::Vector{Axis}, 
                               baseline_interval::BaselineInterval,
                               line_refs::Union{Vector{<:Dict},Nothing} = nothing)

Set up a control panel that opens when 'c' key is pressed.
Allows adjusting baseline and toggling conditions.
"""
function _setup_erp_control_panel!(
    fig::Figure,
    dat_subset::Vector{ErpData},
    axes::Vector{Axis},
    baseline_interval::BaselineInterval,
    line_refs::Union{Vector{<:Dict},Nothing} = nothing,
)

    control_fig = Ref{Union{Figure,Nothing}}(nothing)

    # Set up linked legend interactions
    if line_refs !== nothing
        _setup_linked_legend_interactions!(line_refs)
    end

    # State: baseline values and condition selections
    start_val, stop_val = _extract_baseline_values(baseline_interval)
    baseline_start_obs = Observable(start_val === nothing ? "" : string(start_val))
    baseline_stop_obs = Observable(stop_val === nothing ? "" : string(stop_val))
    condition_checked = [Observable(true) for _ in dat_subset]

    # Track previous baseline to avoid unnecessary updates
    previous_baseline = Ref{Union{Tuple{Float64,Float64},Nothing}}(nothing)

    function update_plot!()
        # Parse baseline values from observables
        start_str, stop_str = baseline_start_obs[], baseline_stop_obs[]
        start_val, stop_val = _parse_baseline_values(start_str, stop_str)

        # Convert to tuple if valid (baseline! accepts tuples and converts internally)
        baseline_interval_new = (start_val !== nothing && stop_val !== nothing) ? (start_val, stop_val) : nothing

        # Check if baseline actually changed
        baseline_changed = baseline_interval_new !== previous_baseline[]

        # Apply baseline if it changed
        if baseline_changed && baseline_interval_new !== nothing
            baseline!.(dat_subset, Ref(baseline_interval_new))
            previous_baseline[] = baseline_interval_new
        end

        # Build condition mask
        condition_mask = [checked[] for checked in condition_checked]

        # Update existing lines (y-data if baseline changed, visibility always)
        for (ax_idx, ax_line_refs) in enumerate(line_refs)
            ax_idx > length(axes) && continue
            for (dataset_idx, channel_lines) in ax_line_refs
                dataset_idx > length(dat_subset) && continue
                dat = dat_subset[dataset_idx]

                for (channel, line_data) in channel_lines
                    line, y_obs = line_data

                    # Update y-data if baseline changed
                    if baseline_changed
                        y_obs[] = dat.data[!, channel]
                    end

                    # Update visibility based on condition_mask
                    line.visible = dataset_idx <= length(condition_mask) ? condition_mask[dataset_idx] : true
                end
            end
        end
    end

    # Keyboard handler for 'c' key
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.c

            # Create new control panel
            control_fig[] = Figure(title = "ERP Control Panel", size = (300, 400))
            layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)

            # Baseline section
            Label(layout[1, 1], "Baseline Interval", fontsize = 14, font = :bold)
            baseline_layout = GridLayout(layout[2, 1], tellwidth = false, colgap = 10)

            _create_baseline_textbox(baseline_layout, 1, "Start (ms):", baseline_start_obs, " ", 100)
            _create_baseline_textbox(baseline_layout, 2, "End (ms):", baseline_stop_obs, " ", 100)

            # Apply button
            apply_btn = Button(layout[3, 1], label = "Apply Baseline", width = 200)
            on(apply_btn.clicks) do _
                update_plot!()
            end

            # Conditions section
            Label(layout[4, 1], "Conditions", fontsize = 14, font = :bold)
            conditions_layout = GridLayout(layout[5, 1], tellwidth = false, rowgap = 5)

            for (i, dat) in enumerate(dat_subset)
                cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][])
                Label(conditions_layout[i, 2], dat.condition_name)
                connect!(condition_checked[i], cb.checked)
            end

            # Auto-update on condition changes (re-plot)
            for checked in condition_checked
                on(checked) do _
                    update_plot!()
                end
            end

            display(control_fig[])

        end
    end
end
