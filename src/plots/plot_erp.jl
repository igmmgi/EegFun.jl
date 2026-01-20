# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ERP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("ERP Plot", "Title for the plot window"),
    :interactive => (true, "Enable interactive features (true/false)"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    :yreversed => (false, "Whether to reverse the y-axis"),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Show title (true/false)"),

    # Line styling
    :linewidth => (2, "Line width for ERPs"),
    :color => (:black, "Color for ERPs (single color or a vector of colors, one per dataset)"),
    :linestyle => (:solid, "Line style for ERPs (single style or a vector of styles, one per dataset)"),
    :colormap => (:jet, "Colormap for multi-channel plots"),

    # Plot configuration
    :yreversed => (false, "Reverse the y-axis (true/false)"),
    :average_channels => (false, "Average across channels (true/false)"),

    # Legend parameters - get all Legend attributes with their actual defaults
    # This allows users to control any Legend parameter
    [
        Symbol("legend_$(attr)") => (get(LEGEND_DEFAULTS, attr, nothing), "Legend $(attr) parameter") for
        attr in propertynames(Legend)
    ]...,

    # Override specific legend parameters with custom defaults
    :legend => (true, "Show the legend (true/false)"),
    :legend_label => ("", "Title for the legend"),
    :legend_framevisible => (true, "Legend frame visible (true/false)"),
    :legend_position => (:lt, "Legend position (:lt, :rt, :lb, :rb, or tuple like (:left, :top), or (0.5, 0.5))"),
    :legend_channel => ([], "Which channel to put the legend on."),
    :legend_labels => ([], "Legend labels."),
    :legend_nbanks => (nothing, "Number of legend columns."),

    # Individual axes grids
    :xgrid => (false, "Show x-axis grid (true/false)"),
    :ygrid => (false, "Show y-axis grid (true/false)"),
    :xminorgrid => (false, "Show x-axis minor grid (true/false)"),
    :yminorgrid => (false, "Show y-axis minor grid (true/false)"),

    # Origin lines
    :add_xy_origin => (true, "Add origin lines at x=0 and y=0 (true/false)"),

    # Layout parameters (for topo and other layouts)
    :layout_topo_plot_width => (0.05, "Width of individual plots (fraction of figure width)"),
    :layout_topo_plot_height => (0.05, "Height of individual plots (fraction of figure height)"),
    :layout_topo_scale_offset => (0.1, "Offset factor for scale plot position"),
    :layout_topo_scale_pos => ((0.8, -0.8), "Fallback position for scale plot in topo layout as (x, y) tuple"),

    # Grid layout parameters
    :layout_grid_rowgap => (10, "Gap between rows (in pixels)"),
    :layout_grid_colgap => (10, "Gap between columns (in pixels)"),
    :layout_grid_dims =>
        (nothing, "Grid dimensions as (rows, cols) tuple for grid layouts. If nothing, automatically determined"),
    :layout_grid_skip_positions =>
        (nothing, "Positions to skip in grid layout as vector of (row, col) tuples, e.g., [(2,1), (2,3)]"),

    # General layout parameters
    :figure_padding =>
        ((10, 10, 10, 10), "Padding around entire figure as (left, right, top, bottom) tuple (in pixels)"),
)

"""
    plot_erp(filepath::String; 
             layout::Union{Symbol, PlotLayout} = :single,
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
    layout::Union{Symbol,PlotLayout} = :single,
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
             layout::Union{Symbol, PlotLayout} = :single,
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
plot_erp(dat, layout = :grid, layout_grid_dims = (3, 4))

# Topographic layout
plot_erp(dat, layout = :topo)

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
  - Topoplot 
"""
function plot_erp(
    dat::ErpData;
    layout::Union{Symbol,PlotLayout} = :single,
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
             layout::Union{Symbol, PlotLayout} = :single,
             condition_selection::Function = conditions(),
             channel_selection::Function = channels(), 
             sample_selection::Function = samples(),
             baseline_interval::BaselineInterval = nothing,
             kwargs...)

Plot multiple ERP datasets on the same axis (e.g., conditions).
"""
function plot_erp(
    datasets::Vector{ErpData};
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    baseline_interval::BaselineInterval = nothing,
    return_line_refs::Bool = false,  # Internal parameter, not in PLOT_ERP_KWARGS
    kwargs...,
)

    # Prepare kwargs and data
    plot_kwargs, user_provided_color = _prepare_plot_kwargs(kwargs)
    dat_subset, all_channels, channel_selection_func, original_channels = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        baseline_interval = baseline_interval,
    )

    # Apply channel_selection to determine which channels to plot
    # dat_subset has all channels, but we only plot the selected ones
    selected_channels =
        get_selected_channels(first(dat_subset), channel_selection_func; include_meta = false, include_extra = true)
    # Preserve order from selected_channels (user's channel_selection order)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Check if any channels remain after filtering
    if isempty(all_plot_channels)
        @minimal_error_throw "No valid channels found. Selected channels: $selected_channels, Available channels: $all_channels"
    end

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

    # Extract layout_* parameters, remove prefix, and pass to create_layout
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    # Create figure and apply layout system
    fig = Figure(title = plot_kwargs[:figure_title], figure_padding = plot_kwargs[:figure_padding])

    plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout; layout_kwargs...)

    # For :topo layout, set default legend_channel to last channel if not explicitly set
    if plot_layout.type == :topo && isempty(plot_kwargs[:legend_channel]) && !isempty(all_plot_channels)
        plot_kwargs[:legend_channel] = [all_plot_channels[end]]
        plot_kwargs[:legend_position] = (5, -2) # and put it a bit outside the plot
    end
    axes, channels = _apply_layout!(fig, plot_layout; plot_kwargs...)

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
        leg = _plot_erp!(
            ax,
            dat_subset,
            channels_to_plot;
            line_refs = ax_line_refs,
            user_provided_color = user_provided_color,
            plot_kwargs...,
        )
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

        # Set up control panel (press 'c' to open) - must be before selection to capture condition_checked
        condition_checked_ref = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing)
        _setup_erp_control_panel!(fig, dat_subset, axes, baseline_interval, line_refs, condition_checked_ref)

        # Create right-click handler that has access to condition visibility
        right_click_handler =
            (selection_state, mouse_x, data) ->
                _handle_erp_right_click!(selection_state, mouse_x, data, condition_checked_ref)

        # Set up selection system that works for all layouts
        _setup_unified_selection!(fig, axes, selection_state, dat_subset, plot_layout, right_click_handler)

        # Set up channel selection events for topo and grid layouts
        if plot_layout.type in (:topo, :grid)
            _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes, plot_layout.type)
        end

    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    # reset default title
    set_window_title("Makie")
    # Optionally return line references if requested
    if return_line_refs
        return fig, axes, line_refs
    else
        return fig, axes
    end
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
    plot_kwargs, user_provided_color = _prepare_plot_kwargs(kwargs)
    baseline_interval = get(kwargs, :baseline_interval, nothing)
    channel_selection_func = get(plot_kwargs, :channel_selection, channels())
    dat_subset, all_channels, _, _ = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = conditions(),
        channel_selection = channel_selection_func,
        sample_selection = get(plot_kwargs, :sample_selection, samples()),
        baseline_interval = baseline_interval,
    )
    # Apply channel_selection to determine which channels to plot
    selected_channels =
        get_selected_channels(first(dat_subset), channel_selection_func; include_meta = false, include_extra = true)
    # Preserve order from selected_channels (user's channel_selection order)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]
    _plot_erp!(ax, dat_subset, all_plot_channels; user_provided_color = user_provided_color, plot_kwargs...)
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
    user_provided_color::Bool = false,
    kwargs...,
)

    # Merge kwargs with defaults (kwargs are already partially merged from calling function)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)

    # Default condition_mask to all true if not provided
    if isempty(condition_mask)
        condition_mask = fill(true, length(datasets))
    end

    # Compute colors and linestyles for each dataset
    all_colors = _compute_dataset_colors(
        plot_kwargs[:color],
        length(datasets),
        length(channels),
        plot_kwargs[:colormap],
        user_provided_color,
    )
    all_linestyles = _compute_dataset_linestyles(plot_kwargs[:linestyle], length(datasets))

    # Plot each dataset for ALL channels in this subplot
    for (dataset_idx, dat) in enumerate(datasets)

        for (channel_idx, channel) in enumerate(channels)

            # axis label
            label = isempty(plot_kwargs[:legend_labels]) ? dat.condition_name : plot_kwargs[:legend_labels][dataset_idx]
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
                linewidth = plot_kwargs[:linewidth],
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

    _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
    leg = _add_legend!(ax, channels, datasets, plot_kwargs)

    return leg
end


"""
    _prepare_plot_kwargs(kwargs)

Prepare plot kwargs by merging with defaults and tracking whether user provided color.
Returns (plot_kwargs, user_provided_color).
"""
function _prepare_plot_kwargs(kwargs)
    user_provided_color = haskey(kwargs, :color)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)
    return plot_kwargs, user_provided_color
end

"""
    _prepare_erp_data(datasets, plot_kwargs; condition_selection, channel_selection, sample_selection)

Prepare ERP data for plotting: subset by condition and sample only (NOT channels).
Channel selection is applied later at plot time to keep all channels available for topoplots.
Returns (dat_subset_unfiltered_channels, all_channels, channel_selection_func, original_channels).
"""
function _prepare_erp_data(
    datasets::Vector{ErpData},
    plot_kwargs;
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    baseline_interval::BaselineInterval = nothing,
)
    # Data subsetting - ONLY by condition and sample, NOT by channel
    # This keeps all channels available for right-click topoplots
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channels(),  # Select ALL channels (no filtering)
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

    # Extract ALL channel labels (not filtered by channel_selection)
    # We need this before averaging to know which channels to average
    all_channels = channel_labels(dat_subset)
    extra_channels = extra_labels(dat_subset)
    all_channels = vcat(all_channels, extra_channels)

    # Apply channel_selection to determine which channels to plot/average
    selected_channels =
        get_selected_channels(first(dat_subset), channel_selection; include_meta = false, include_extra = true)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Channel averaging if requested - average only the selected channels
    original_channels = nothing
    if plot_kwargs[:average_channels]
        @info "Averaging channels for $(length(dat_subset)) datasets"
        original_channels = all_plot_channels
        # Average only the selected channels, similar to plot_epochs
        for dat in dat_subset
            channel_average!(
                dat;
                channel_selections = [channels(all_plot_channels)],
                output_labels = [:avg],
                reduce = false,
            )
        end
        # After averaging, update all_channels to include the averaged channel
        all_channels = [:avg]
        # Update channel_selection to select :avg instead of original channels
        channel_selection = channels(:avg)
    end

    return dat_subset, all_channels, channel_selection, original_channels
end


function _handle_erp_right_click!(selection_state, mouse_x, data, condition_checked_ref)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_erp_context_menu!(selection_state, data, condition_checked_ref)
    end
end

function _show_erp_context_menu!(selection_state, data, condition_checked_ref)

    menu_fig = Figure(size = (300, 300))

    # Filter by visible conditions to determine if we have multiple visible conditions
    data_to_plot = _filter_visible_conditions(data, condition_checked_ref)
    has_multiple_conditions = data_to_plot isa Vector{ErpData} && length(data_to_plot) > 1

    plot_types = ["Topoplot"]

    # Only add average options if multiple visible conditions
    if has_multiple_conditions
        push!(plot_types, "Topoplot (average)")
    end

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            original_data, x_min, x_max = _get_erp_selection_bounds(selection_state, data)

            # Create time-based sample selection for the topo plot
            time_sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max)

            # Filter by visible conditions if condition_checked is available (already done above, but do again for consistency)
            data_to_plot = _filter_visible_conditions(original_data, condition_checked_ref)

            if btn.label[] == "Topoplot"
                plot_topography(data_to_plot, sample_selection = time_sample_selection)
            elseif btn.label[] == "Topoplot (average)"
                avg_data = _average_conditions(data_to_plot)
                plot_topography(avg_data, sample_selection = time_sample_selection)
            end
        end
    end

    new_screen = GLMakie.Screen(size = (300, 300))
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
    _filter_visible_conditions(data, condition_checked_ref)

Filter data by visible conditions from the control panel.
Returns filtered Vector{ErpData} or single ErpData if only one condition.
"""
function _filter_visible_conditions(data, condition_checked_ref)
    # If no condition_checked available or data is not Vector, return as-is
    if condition_checked_ref[] === nothing || !(data isa Vector{ErpData})
        return data
    end

    condition_checked = condition_checked_ref[]
    if length(condition_checked) != length(data)
        return data  # Mismatch, return as-is
    end

    # Filter by visible conditions
    visible_data = [data[i] for i in eachindex(data) if condition_checked[i][]]

    # Return single ErpData if only one visible, otherwise Vector
    return length(visible_data) == 1 ? visible_data[1] : visible_data
end

"""
    _average_conditions(erps::Vector{ErpData}) -> ErpData

Average multiple ErpData conditions together into a single ErpData.
Reuses _create_grand_average which averages ErpData objects together.
"""
function _average_conditions(erps::Vector{ErpData})
    if isempty(erps)
        @minimal_error_throw("Cannot average empty ERP list")
    end
    if length(erps) == 1
        return erps[1]  # Nothing to average
    end

    # Reuse _create_grand_average - it averages ErpData together (same logic for conditions or participants)
    # Use first condition number as cond_num (doesn't matter for averaging)
    avg_erp = _create_grand_average(erps, first(erps).condition)

    # Update condition name to reflect averaging across conditions
    avg_cond_name = "avg_" * join([erp.condition_name for erp in erps], "_")
    return ErpData(
        avg_erp.file,
        avg_erp.condition,
        avg_cond_name,
        avg_erp.data,
        avg_erp.layout,
        avg_erp.sample_rate,
        avg_erp.analysis_info,
        avg_erp.n_epochs,
    )
end

"""
    _average_conditions(erp::ErpData) -> ErpData

Single ErpData case - just return it.
"""
_average_conditions(erp::ErpData) = erp

"""
    _compute_dataset_colors(color_val, n_datasets, n_channels, colormap, color_explicitly_set)

Compute colors for each dataset-channel combination.
Returns a vector of colors with length n_datasets * n_channels.
Colors cycle across all channel-dataset combinations.
"""
function _compute_dataset_colors(color_val, n_datasets::Int, n_channels::Int, colormap, user_provided_color::Bool)
    n_total = n_datasets * n_channels

    # If user provided a vector of colors, use those (cycle if needed)
    if color_val isa Vector
        return [color_val[(i-1)%length(color_val)+1] for i = 1:n_total]
    end

    # If user provided a single color, use it for all items
    if user_provided_color
        return [color_val for _ = 1:n_total]
    end

    # User didn't provide color: use colormap for multiple items, default color for single item
    if n_total > 1
        # Makie.cgrad returns a gradient, convert to vector of colors
        gradient = Makie.cgrad(colormap, n_total, categorical = true)
        return [gradient[i] for i = 1:n_total]
    else
        return [color_val for _ = 1:n_total]  # Use default color for single item
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
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}} = Ref{Union{Vector{Observable{Bool}},Nothing}}(
        nothing,
    ),
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
    condition_checked_ref[] = condition_checked  # Store for access by right-click handler

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

            for (idx, dat) in enumerate(dat_subset)
                cb = Checkbox(conditions_layout[idx, 1], checked = condition_checked[idx][])
                Label(conditions_layout[idx, 2], dat.condition_name)
                connect!(condition_checked[idx], cb.checked)
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
