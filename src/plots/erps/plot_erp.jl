# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ERP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("ERP Plot", "Title for the plot window"),
    :interactive => (true, "Enable interactive features (true/false)"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :selection_color => (:blue, "Color for interactive selection rectangles"),
    :selection_alpha => (0.3, "Alpha (transparency) for interactive selection rectangles"),
    :theme_fontsize => (24, "Font size for theme"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    :yreversed => (false, "Whether to reverse the y-axis"),
    :xticks => (nothing, "X-axis tick positions (e.g., -0.1:0.1:0.8 or [0, 0.2, 0.4]). If nothing, automatically determined"),
    :yticks => (nothing, "Y-axis tick positions (e.g., -4:2:4 or [-2, 0, 2]). If nothing, automatically determined"),
    :time_unit => (
        :s,
        "Time unit for x-axis display (:s or :ms). Only affects axis labels and tick formatting — all intervals remain in seconds.",
    ),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Show title (true/false)"),

    # Line styling
    :linewidth => (2, "Line width for ERPs"),
    :color => (:black, "Color for ERPs (single color or a vector of colors, one per dataset)"),
    :linestyle => (:solid, "Line style for ERPs (single style or a vector of styles, one per dataset)"),
    :colormap => (:jet, "Colormap for multi-channel plots"),

    # Plot configuration
    :average_channels => (false, "Average across channels (true/false)"),
    :error_bars => (:none, "Error bars type to plot: :none, :sem, :within_sem, :ci95"),

    # Legend parameters - get all Legend attributes with their actual defaults
    # This allows users to control any Legend parameter
    [Symbol("legend_$(attr)") => (get(LEGEND_DEFAULTS, attr, nothing), "Legend $(attr) parameter") for attr in propertynames(Legend)]...,

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

    # Origin lines and scale indicator
    :add_xy_origin => (true, "Add origin lines at x=0 and y=0 (true/false)"),
    :axis_type => (:standard, "Type of axis to draw (:standard or :origin)"),
    :scale_x_value => (nothing, "X-axis scale value/step size (e.g. 0.1 for 100 ms)"),
    :scale_y_value => (nothing, "Y-axis scale value/step size (e.g. 5.0 for 5 μV)"),

    # Layout parameters (for topo and other layouts)
    :layout_topo_plot_width => (0.05, "Width of individual plots (fraction of figure width)"),
    :layout_topo_plot_height => (0.05, "Height of individual plots (fraction of figure height)"),
    :layout_topo_scale_offset => (0.1, "Offset factor for scale plot position"),
    :layout_topo_scale_pos => ((0.8, -0.8), "Fallback position for scale plot in topo layout as (x, y) tuple"),

    # Grid layout parameters
    :layout_grid_rowgap => (10, "Gap between rows (in pixels)"),
    :layout_grid_colgap => (10, "Gap between columns (in pixels)"),
    :layout_grid_dims => (nothing, "Grid dimensions as (rows, cols) tuple for grid layouts. If nothing, automatically determined"),
    :layout_grid_skip_positions => (nothing, "Positions to skip in grid layout as vector of (row, col) tuples, e.g., [(2,1), (2,3)]"),

    # General layout parameters
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),

    # Highlight regions
    :highlight_regions => (
        nothing,
        "Highlight regions as a NamedTuple or Vector of NamedTuples. Each region: (x1, x2, y1=-Inf, y2=Inf, color=:gray, alpha=0.3)",
    ),
)

"""
    plot_erp(filepath::String; 
             input_dir::String = pwd(),
             participant_selection::Function = participants(),
             layout::Union{Symbol, PlotLayout} = :single,
             condition_selection::Function = conditions(),
             channel_selection::Function = channels(),
             channel_plot_order::Union{Nothing, Vector{Symbol}} = nothing,
             sample_selection::Function = samples(),
             interval_selection::Interval = times(),
             baseline_interval::Interval = nothing,
             kwargs...)

Load ERP data and create plots. Accepts either a direct `.jld2` filepath or a filename 
pattern to discover and plot all matching files.

# Arguments
- `filepath::String`: Either a `.jld2` file path, or a pattern string (e.g. `"erps_final"`) 
  to match against files in `input_dir`
- `input_dir::String`: Directory to search for pattern-matched files (default: `pwd()`)
- `participant_selection::Function`: Participant filter for pattern mode (default: `participants()`)
- `layout`: Layout specification (see main plot_erp documentation)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering
- `sample_selection::Function`: Function that returns boolean vector for sample filtering
- `kwargs`: Additional keyword arguments

$(_generate_kwargs_doc(PLOT_ERP_KWARGS))

# Examples
```julia
# Load and plot from file
plot_erp("grand_average_erps_final.jld2")

# Plot all files matching pattern in current directory
plot_erp("erps_final")

# Plot specific participant
plot_erp("erps_final", participant_selection = participants(1))

# With channel selection
plot_erp("erps_final", channel_selection = channels([:PO7, :PO8]), layout = :grid)
```
"""
function plot_erp(
    filepath::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    # Detect: filepath (ends with .jld2) vs pattern
    if endswith(filepath, ".jld2")
        # Direct file path — existing behavior
        data = read_data(filepath)
        if isnothing(data)
            @minimal_error "No data found in file: $filepath"
        end

        # Dispatch will handle ErpData vs Vector{ErpData} automatically
        return plot_erp(
            data;
            layout = layout,
            condition_selection = condition_selection,
            channel_selection = channel_selection,
            channel_plot_order = channel_plot_order,
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            baseline_interval = baseline_interval,
            kwargs...,
        )
    else
        # Pattern-based discovery — one plot per file
        files = _find_batch_files(filepath, input_dir, participant_selection)
        if isempty(files)
            all_matching = _find_batch_files(filepath, input_dir)  # without participant filter
            if isempty(all_matching)
                @minimal_error "No ERP files matching pattern '$filepath' in $input_dir"
            else
                avail_ids = sort(unique(_extract_participant_id.(all_matching)))
                @minimal_error "Pattern '$filepath' matched $(length(all_matching)) file(s), but none passed the participant selection. Available participant IDs: $avail_ids"
            end
        end

        results = NamedTuple[]
        for file in sort(files, by = _natural_sort_key)
            file_path = joinpath(input_dir, file)
            @info "Plotting: $file"
            data = read_data(file_path)
            isnothing(data) && continue

            result = plot_erp(
                data;
                layout = layout,
                condition_selection = condition_selection,
                channel_selection = channel_selection,
                channel_plot_order = channel_plot_order,
                sample_selection = sample_selection,
                interval_selection = interval_selection,
                baseline_interval = baseline_interval,
                kwargs...,
            )
            push!(results, result)
        end
        return results
    end
end

function plot_erp(
    dat::ErpData;
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    # For single ErpData, condition_selection doesn't apply (there's only one condition)
    return plot_erp(
        [dat];
        layout = layout,
        condition_selection = conditions(),  # Always select all (just the one condition)
        channel_selection = channel_selection,
        channel_plot_order = channel_plot_order,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

function plot_erp(
    datasets::Vector{ErpData};
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    # If there are datasets from multiple participants, automatically create separate windows
    file_ids = unique(d.file for d in datasets)
    if length(file_ids) > 1
        results = NamedTuple[]
        for file_id in file_ids
            participant_data = filter(d -> d.file == file_id, datasets)
            result = plot_erp(
                participant_data;
                layout = layout,
                condition_selection = condition_selection,
                channel_selection = channel_selection,
                channel_plot_order = channel_plot_order,
                sample_selection = sample_selection,
                interval_selection = interval_selection,
                baseline_interval = baseline_interval,
                kwargs...,
            )
            push!(results, result)
        end
        return results
    end

    # Prepare kwargs and data
    plot_kwargs, user_provided_color = _prepare_plot_kwargs(kwargs)
    dat_subset, all_channels, channel_selection_func, original_channels = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
    )

    variances_dict = Dict{Int,Tuple{DataFrames.DataFrame,Int}}()

    return _plot_erp_core(
        dat_subset,
        all_channels,
        channel_selection_func,
        original_channels,
        variances_dict,
        plot_kwargs,
        user_provided_color,
        layout,
        baseline_interval,
        channel_plot_order,
    )
end

function plot_erp_errorbar(
    filepath::String;
    input_dir::String = "",
    participant_selection::Function = participants(),
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    # Pattern-based discovery — one plot per file? No, plot_erp_errorbar groups them!
    files = _find_batch_files(filepath, input_dir, participant_selection)
    if isempty(files)
        all_matching = _find_batch_files(filepath, input_dir)  # without participant filter
        if isempty(all_matching)
            @minimal_error "No ERP files matching pattern '$filepath' in $input_dir"
        else
            avail_ids = sort(unique(_extract_participant_id.(all_matching)))
            @minimal_error "Pattern '$filepath' matched $(length(all_matching)) file(s), but none passed the participant selection. Available participant IDs: $avail_ids"
        end
    end

    datasets = ErpData[]
    for file in sort(files, by = _natural_sort_key)
        file_path = joinpath(input_dir, file)
        @info "Loading: $file"
        data = read_data(file_path)
        isnothing(data) && continue

        if data isa Vector{ErpData}
            append!(datasets, data)
        else
            push!(datasets, data)
        end
    end

    if isempty(datasets)
        @minimal_error "No data loaded from pattern: $filepath"
    end

    return plot_erp_errorbar(
        datasets;
        layout = layout,
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        channel_plot_order = channel_plot_order,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

function plot_erp_errorbar(
    datasets::Vector{ErpData};
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
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
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
    )

    @info "Averaging across participants dynamically..."
    dat_subset, variances_dict = _compute_dynamic_grand_averages(dat_subset, plot_kwargs[:error_bars])



    return _plot_erp_core(
        dat_subset,
        all_channels,
        channel_selection_func,
        original_channels,
        variances_dict,
        plot_kwargs,
        user_provided_color,
        layout,
        baseline_interval,
        channel_plot_order,
    )
end


function _plot_erp_core(
    dat_subset,
    all_channels,
    channel_selection_func,
    original_channels,
    variances_dict,
    plot_kwargs,
    user_provided_color,
    layout,
    baseline_interval,
    channel_plot_order,
)
    # Apply channel_selection to determine which channels to plot
    # dat_subset has all channels, but we only plot the selected ones
    selected_channels = get_selected_channels(first(dat_subset), channel_selection_func; include_meta = false, include_extra = true)
    # Preserve order from selected_channels (user's channel_selection order)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Apply user-specified plotting order if provided
    if !isnothing(channel_plot_order)
        ordered = Symbol[ch for ch in channel_plot_order if ch in all_plot_channels]
        if !isempty(ordered)
            all_plot_channels = ordered
        end
    end

    # Check if any channels remain after filtering
    if isempty(all_plot_channels)
        @minimal_error "No valid channels found. Selected channels: $selected_channels, Available channels: $all_channels"
    end

    # set default plot title only for single layouts
    # For grid/topo layouts, we want individual channel names, not a global title
    if plot_kwargs[:show_title] && plot_kwargs[:title] == "" && layout == :single
        plot_kwargs[:title] = length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(_print_vector(all_plot_channels))"
        if plot_kwargs[:average_channels]
            plot_kwargs[:title] = "Avg: $(_print_vector(original_channels))"
        end
    end

    # Generate window title from datasets
    title_str = _generate_window_title(dat_subset)
    _set_window_title(title_str)

    # Extract layout_* parameters, remove prefix, and pass to create_layout
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    # Create figure and apply layout system
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])
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
        @info "plot_erp ($layout): $(_print_vector(channels_to_plot))"
        ax_line_refs = plot_kwargs[:interactive] ? line_refs[ax_idx] : nothing
        leg = _plot_erp!(
            ax,
            dat_subset,
            channels_to_plot;
            line_refs = ax_line_refs,
            user_provided_color = user_provided_color,
            variances_dict = variances_dict,
            plot_kwargs...,
        )
        if plot_kwargs[:interactive] && !isnothing(legend_refs)
            legend_refs[ax_idx] = leg
        end
    end

    # Draw highlight regions on all axes (before axis properties so they render behind)
    if !isnothing(plot_kwargs[:highlight_regions])
        for ax in axes
            _draw_highlight_regions!(ax, plot_kwargs[:highlight_regions])
        end
    end

    # Apply our axis stuff
    _apply_axis_properties!.(axes; plot_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...) # slightly different for grid and topo layouts

    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)

    # Add keyboard interactivity if enabled
    if plot_kwargs[:interactive]
        _setup_shared_interactivity!(fig, axes, :erp; zoom_step = plot_kwargs[:zoom_step])

        # Disable default interactions that conflict with our custom selection (all axes)
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end

        # Set up selection system for all axes (will work with linked axes)
        selection_state =
            SharedSelectionState(axes; selection_color = plot_kwargs[:selection_color], selection_alpha = plot_kwargs[:selection_alpha])

        # Set up control panel (press 'c' to open) - must be before selection to capture condition_checked
        condition_checked_ref = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing)
        _setup_erp_control_panel!(fig, dat_subset, axes, baseline_interval, line_refs, condition_checked_ref)

        # Create right-click handler that has access to condition visibility
        right_click_handler =
            (selection_state, mouse_x, data) -> _handle_erp_right_click!(selection_state, mouse_x, data, condition_checked_ref)

        # Set up selection system that works for all layouts
        _setup_unified_selection!(fig, axes, selection_state, dat_subset, right_click_handler)

        # Set up channel selection events for topo and grid layouts
        if plot_layout.type in (:topo, :grid)
            channel_rc_handler = (selected_chs, data) -> _show_channel_average_menu!(selected_chs, data, condition_checked_ref)
            _setup_channel_selection_events!(
                fig,
                selection_state,
                plot_layout,
                dat_subset,
                axes;
                channel_right_click_handler = channel_rc_handler,
            )
        end

    end

    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end

    # reset default title
    _set_window_title("Makie")
    # Return named tuple - line_refs available for advanced use (e.g., plot_erp_measurements)
    # Most users can ignore it via: (; fig, axes) = plot_erp(...)
    return (fig = fig, axes = axes, line_refs = line_refs)
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
function plot_erp!(fig::Figure, ax::Axis, datasets::Vector{ErpData}; kwargs...)
    # Extract special parameters before validation (like plot_erp does)
    baseline_interval = get(kwargs, :baseline_interval, nothing)
    channel_selection_func = get(kwargs, :channel_selection, channels())
    channel_plot_order = get(kwargs, :channel_plot_order, nothing)

    # Remove them from kwargs before validation
    filtered_kwargs = pairs(NamedTuple(filter(p -> p[1] ∉ [:baseline_interval, :channel_selection, :channel_plot_order], pairs(kwargs))))

    plot_kwargs, user_provided_color = _prepare_plot_kwargs(filtered_kwargs)

    dat_subset, all_channels, _, _ = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = conditions(),
        channel_selection = channel_selection_func,
        sample_selection = get(plot_kwargs, :sample_selection, samples()),
        baseline_interval = baseline_interval,
    )
    variances_dict = Dict{Int,Tuple{DataFrames.DataFrame,Int}}()

    # Apply channel_selection to determine which channels to plot
    selected_channels = get_selected_channels(first(dat_subset), channel_selection_func; include_meta = false, include_extra = true)
    # Preserve order from selected_channels (user's channel_selection order)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Apply user-specified plotting order if provided
    if !isnothing(channel_plot_order)
        ordered = Symbol[ch for ch in channel_plot_order if ch in all_plot_channels]
        if !isempty(ordered)
            all_plot_channels = ordered
        end
    end

    _plot_erp!(
        ax,
        dat_subset,
        all_plot_channels;
        user_provided_color = user_provided_color,
        variances_dict = variances_dict,
        plot_kwargs...,
    )
    return ax
end

"""
    plot_erp_errorbar!(fig::Figure, ax::Axis, datasets::Vector{ErpData}; kwargs...)

Plot grand-averaged ERP data with error ribbons on an existing axis, mutating the figure and axis.
"""
plot_erp_errorbar!(fig::Figure, ax::Axis, dat::ErpData; kwargs...) = plot_erp_errorbar!(fig, ax, [dat]; kwargs...)
function plot_erp_errorbar!(fig::Figure, ax::Axis, datasets::Vector{ErpData}; kwargs...)
    # Extract special parameters before validation (like plot_erp does)
    baseline_interval = get(kwargs, :baseline_interval, nothing)
    channel_selection_func = get(kwargs, :channel_selection, channels())
    channel_plot_order = get(kwargs, :channel_plot_order, nothing)

    # Remove them from kwargs before validation
    filtered_kwargs = pairs(NamedTuple(filter(p -> p[1] ∉ [:baseline_interval, :channel_selection, :channel_plot_order], pairs(kwargs))))

    plot_kwargs, user_provided_color = _prepare_plot_kwargs(filtered_kwargs)

    dat_subset, all_channels, _, _ = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = conditions(),
        channel_selection = channel_selection_func,
        sample_selection = get(plot_kwargs, :sample_selection, samples()),
        baseline_interval = baseline_interval,
    )

    @info "Averaging across participants dynamically..."
    dat_subset, variances_dict = _compute_dynamic_grand_averages(dat_subset, plot_kwargs[:error_bars])

    # Apply channel_selection to determine which channels to plot
    selected_channels = get_selected_channels(first(dat_subset), channel_selection_func; include_meta = false, include_extra = true)
    # Preserve order from selected_channels (user's channel_selection order)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Apply user-specified plotting order if provided
    if !isnothing(channel_plot_order)
        ordered = Symbol[ch for ch in channel_plot_order if ch in all_plot_channels]
        if !isempty(ordered)
            all_plot_channels = ordered
        end
    end

    _plot_erp!(
        ax,
        dat_subset,
        all_plot_channels;
        user_provided_color = user_provided_color,
        variances_dict = variances_dict,
        plot_kwargs...,
    )
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
    variances_dict::Dict{Int,Tuple{DataFrames.DataFrame,Int}} = Dict{Int,Tuple{DataFrames.DataFrame,Int}}(),
    kwargs...,
)

    # Merge kwargs with defaults (kwargs are already partially merged from calling function)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)

    # Default condition_mask to all true if not provided
    if isempty(condition_mask)
        condition_mask = fill(true, length(datasets))
    end

    # Compute colors and linestyles for each dataset
    all_colors =
        _compute_dataset_colors(plot_kwargs[:color], length(datasets), length(channels), plot_kwargs[:colormap], user_provided_color)
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

            error_type = plot_kwargs[:error_bars]
            if error_type != :none && haskey(variances_dict, dataset_idx)
                var_df, n_stat = variances_dict[dataset_idx]
                if hasproperty(var_df, channel)
                    var_col = var_df[!, channel]
                    sem = sqrt.(var_col) ./ sqrt(n_stat)
                    err_val = error_type == :ci95 ? sem .* 1.96 : sem

                    # Create bounds tied to the Observable mean so they shift correctly if baseline changes interactively
                    lower_bound = lift(y -> y .- err_val, y_obs)
                    upper_bound = lift(y -> y .+ err_val, y_obs)

                    band_color = to_color(all_colors[color_idx])
                    band_color_alpha = RGBAf(band_color.r, band_color.g, band_color.b, 0.2)

                    band!(ax, dat.data[!, :time], lower_bound, upper_bound, color = band_color_alpha, visible = condition_visible)
                end
            end

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
            if !isnothing(line_refs)
                if !haskey(line_refs, dataset_idx)
                    line_refs[dataset_idx] = Dict{Symbol,Tuple{Any,Observable}}()
                end
                line_refs[dataset_idx][channel] = (line, y_obs)
            end
        end
    end

    _set_axis_properties!(
        ax;
        xlim = plot_kwargs[:xlim],
        ylim = plot_kwargs[:ylim],
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
        yreversed = plot_kwargs[:yreversed],
        scale_x_value = plot_kwargs[:scale_x_value],
        scale_y_value = plot_kwargs[:scale_y_value],
    )
    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )
    _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
    _add_origin_scale_indicator!(
        ax;
        axis_type = plot_kwargs[:axis_type],
        scale_x_value = plot_kwargs[:scale_x_value],
        scale_y_value = plot_kwargs[:scale_y_value],
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
    )
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
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
)
    # Data subsetting - ONLY by condition and sample, NOT by channel
    # This keeps all channels available for right-click topoplots
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channels(),  # Select ALL channels (no filtering)
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = true,
    )

    # Check if subsetting resulted in empty data
    if isempty(dat_subset)
        n_conditions = length(datasets)
        if n_conditions == 0
            @minimal_error "No data available (empty dataset)"
        else
            @minimal_error "No data matched the selection criteria. Available condition indices: 1:$n_conditions"
        end
    end

    # Apply baseline correction if requested
    if !isnothing(baseline_interval)
        @info "Applying baseline correction to $(length(dat_subset)) datasets"
        baseline!.(dat_subset, Ref(baseline_interval))
    end

    # Extract ALL channel labels (not filtered by channel_selection)
    # We need this before averaging to know which channels to average
    all_channels = channel_labels(dat_subset)
    extra_channels = extra_labels(dat_subset)
    all_channels = vcat(all_channels, extra_channels)

    # Apply channel_selection to determine which channels to plot/average
    selected_channels = get_selected_channels(first(dat_subset), channel_selection; include_meta = false, include_extra = true)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Channel averaging if requested - average only the selected channels
    original_channels = nothing
    if plot_kwargs[:average_channels]
        @info "Averaging channels for $(length(dat_subset)) datasets"
        original_channels = all_plot_channels
        # Average only the selected channels, similar to plot_epochs
        for dat in dat_subset
            channel_average!(dat; channel_selections = [channels(all_plot_channels)], output_labels = [:avg], reduce = false)
        end
        # After averaging, update all_channels to include the averaged channel
        all_channels = [:avg]
        # Update channel_selection to select :avg instead of original channels
        channel_selection = channels(:avg)
    end

    return dat_subset, all_channels, channel_selection, original_channels
end

function _compute_dynamic_grand_averages(datasets::Vector{ErpData}, error_bars::Symbol)
    # Check if they are already grand averaged
    if all(dat -> dat.file == "grand_avg", datasets)
        @minimal_warning "Datasets are already grand averages. Skipping dynamic participant averaging."
        return datasets, Dict{Int,Tuple{DataFrames.DataFrame,Int}}()
    end

    # Group by condition
    erps_by_condition = OrderedDict{Int,Vector{ErpData}}()
    for dat in datasets
        cond = dat.condition
        if !haskey(erps_by_condition, cond)
            erps_by_condition[cond] = ErpData[]
        end
        push!(erps_by_condition[cond], dat)
    end

    grand_averages = ErpData[]
    variances_dict = Dict{Int,Tuple{DataFrames.DataFrame,Int}}()

    # Pre-compute participant global means if Cousineau-Morey is requested
    participant_global_means = Dict{String,DataFrames.DataFrame}()
    grand_global_mean = DataFrames.DataFrame()
    if error_bars == :within_sem
        participant_erps = Dict{String,Vector{ErpData}}()
        for erp in datasets
            pid = erp.file
            if !haskey(participant_erps, pid)
                participant_erps[pid] = ErpData[]
            end
            push!(participant_erps[pid], erp)
        end

        first_erp = datasets[1]
        metadata_cols = meta_labels(first_erp)
        all_channel_sets = [setdiff(propertynames(erp.data), metadata_cols) for erp in datasets]
        eeg_channels = collect(intersect(all_channel_sets...))
        n_points = nrow(first_erp.data)

        for (pid, erplist) in participant_erps
            df_mean = DataFrames.DataFrame()
            n_erps_p = length(erplist)
            for ch in eeg_channels
                avg_buf = zeros(Float64, n_points)
                for erp in erplist
                    col_data = erp.data[!, ch]::Vector{Float64}
                    @inbounds @simd for i = 1:n_points
                        avg_buf[i] += col_data[i]
                    end
                end
                avg_buf ./= n_erps_p
                df_mean[!, ch] = avg_buf
            end
            participant_global_means[pid] = df_mean
        end

        n_participants = length(participant_global_means)
        if n_participants > 0
            for ch in eeg_channels
                avg_buf = zeros(Float64, n_points)
                for df_mean in values(participant_global_means)
                    col_data = df_mean[!, ch]::Vector{Float64}
                    @inbounds @simd for i = 1:n_points
                        avg_buf[i] += col_data[i]
                    end
                end
                avg_buf ./= n_participants
                grand_global_mean[!, ch] = avg_buf
            end
        end
    end

    K_conditions = length(erps_by_condition)
    cm_correction = K_conditions > 1 ? (K_conditions / (K_conditions - 1)) : 1.0

    dataset_idx = 1
    for cond_num in sort(collect(keys(erps_by_condition)))
        erps = erps_by_condition[cond_num]
        first_erp = erps[1]
        metadata_cols = meta_labels(first_erp)
        all_channel_sets = [setdiff(propertynames(erp.data), metadata_cols) for erp in erps]
        eeg_channels = collect(intersect(all_channel_sets...))

        grand_avg_data = copy(first_erp.data)
        cols_to_remove = [:condition, :condition_name, :n_epochs]
        for col in cols_to_remove
            if hasproperty(grand_avg_data, col)
                select!(grand_avg_data, Not(col))
            end
        end

        var_df = DataFrames.DataFrame()
        for col in metadata_cols
            if hasproperty(grand_avg_data, col)
                var_df[!, col] = grand_avg_data[!, col]
            end
        end

        n_points = nrow(grand_avg_data)
        n_erps = length(erps)

        avg_buf = Vector{Float64}(undef, n_points)
        for ch in eeg_channels
            fill!(avg_buf, 0.0)
            for erp in erps
                col_data = erp.data[!, ch]::Vector{Float64}
                @inbounds @simd for i = 1:n_points
                    avg_buf[i] += col_data[i]
                end
            end
            @inbounds @simd for i = 1:n_points
                avg_buf[i] /= n_erps
            end
            grand_avg_data[!, ch] = copy(avg_buf)

            if n_erps > 1 && error_bars != :none
                var_buf = zeros(Float64, n_points)
                if error_bars == :within_sem
                    for erp in erps
                        pid = erp.file
                        col_data = erp.data[!, ch]::Vector{Float64}
                        p_mean = participant_global_means[pid][!, ch]::Vector{Float64}
                        g_mean = grand_global_mean[!, ch]::Vector{Float64}

                        @inbounds @simd for i = 1:n_points
                            y_corr = col_data[i] - p_mean[i] + g_mean[i]
                            diff = y_corr - avg_buf[i]
                            var_buf[i] += diff * diff
                        end
                    end
                    @inbounds @simd for i = 1:n_points
                        var_buf[i] = (var_buf[i] / (n_erps - 1)) * cm_correction
                    end
                else
                    for erp in erps
                        col_data = erp.data[!, ch]::Vector{Float64}
                        @inbounds @simd for i = 1:n_points
                            diff = col_data[i] - avg_buf[i]
                            var_buf[i] += diff * diff
                        end
                    end
                    @inbounds @simd for i = 1:n_points
                        var_buf[i] /= (n_erps - 1)
                    end
                end
                var_df[!, ch] = var_buf
            end
        end

        total_epochs = sum(erp.n_epochs for erp in erps)
        grand_avg = ErpData(
            "grand_avg",
            cond_num,
            "grand_avg_$(first_erp.condition_name)",
            grand_avg_data,
            copy(first_erp.layout),
            first_erp.sample_rate,
            copy(first_erp.analysis_info),
            total_epochs,
        )

        push!(grand_averages, grand_avg)
        if error_bars != :none && n_erps > 1
            variances_dict[dataset_idx] = (var_df, n_erps)
        end
        dataset_idx += 1
    end

    return grand_averages, variances_dict
end


function _handle_erp_right_click!(selection_state, mouse_x, data, condition_checked_ref)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_erp_context_menu!(selection_state, data, condition_checked_ref)
    end
end

function _show_erp_context_menu!(selection_state, data, condition_checked_ref)

    _set_window_title("ERP Context Menu")
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
    if isnothing(condition_checked_ref[]) || !(data isa Vector{ErpData})
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
        @minimal_error("Cannot average empty ERP list")
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
    _show_channel_average_menu!(selected_channels, data, condition_checked_ref)

Show a popup context menu for channel selection right-click.
Offers "Plot Average of Selected" which opens a new ERP plot averaging the selected channels.
"""
function _show_channel_average_menu!(selected_channels, data, condition_checked_ref)
    _set_window_title("Channel Selection Menu")
    menu_fig = Figure(size = (300, 150))

    data_to_plot = _filter_visible_conditions(data, condition_checked_ref)

    btn = Button(menu_fig[1, 1], label = "Plot Average of Selected ($(length(selected_channels)) ch)")
    on(btn.clicks) do _
        plot_erp(data_to_plot; channel_selection = channels(selected_channels), average_channels = true)
    end

    new_screen = GLMakie.Screen(size = (300, 150))
    display(new_screen, menu_fig)
end


"""
    _draw_highlight_regions!(ax::Axis, regions)

Draw shaded rectangular highlight regions on an axis.
Accepts a single NamedTuple or a Vector of NamedTuples.
Each region supports: x1, x2 (required), y1 (default -Inf), y2 (default Inf),
color (default :gray), alpha (default 0.3).
"""
function _draw_highlight_regions!(ax::Axis, regions)
    # Normalize single region to vector
    region_list = regions isa AbstractVector ? regions : [regions]

    for region in region_list
        x1 = region.x1
        x2 = region.x2
        color = get(region, :color, :gray)
        alpha = get(region, :alpha, 0.3)
        has_y1 = haskey(region, :y1)
        has_y2 = haskey(region, :y2)

        if has_y1 || has_y2
            # Explicit y bounds — use poly! for a bounded rectangle
            y1 = get(region, :y1, -1e10)
            y2 = get(region, :y2, 1e10)
            rect = [Point2f(x1, y1), Point2f(x2, y1), Point2f(x2, y2), Point2f(x1, y2)]
            poly!(ax, rect, color = (color, alpha), strokewidth = 0)
        else
            # No y bounds — use vspan! which spans full axis height without affecting ylim
            vspan!(ax, x1, x2, color = (color, alpha))
        end
    end
end

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
    if isnothing(kwargs[:legend_nbanks])
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
                               baseline_interval::Interval,
                               line_refs::Union{Vector{<:Dict},Nothing} = nothing)

Set up a control panel that opens when 'c' key is pressed.
Allows adjusting baseline and toggling conditions.
"""
function _setup_erp_control_panel!(
    fig::Figure,
    dat_subset::Vector{ErpData},
    axes::Vector{Axis},
    baseline_interval::Interval,
    line_refs::Union{Vector{<:Dict},Nothing} = nothing,
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}} = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing),
)

    control_fig = Ref{Union{Figure,Nothing}}(nothing)

    # Set up linked legend interactions
    if !isnothing(line_refs)
        _setup_linked_legend_interactions!(line_refs)
    end

    # State: baseline values and condition selections
    start_val, stop_val = _extract_baseline_values(baseline_interval)
    baseline_start_obs = Observable(isnothing(start_val) ? "" : string(start_val))
    baseline_stop_obs = Observable(isnothing(stop_val) ? "" : string(stop_val))
    condition_checked = [Observable(true) for _ in dat_subset]
    condition_checked_ref[] = condition_checked  # Store for access by right-click handler

    # Track previous baseline to avoid unnecessary updates
    previous_baseline = Ref{Union{Tuple{Float64,Float64},Nothing}}(nothing)

    function update_plot!()
        # Parse baseline values from observables
        start_str, stop_str = baseline_start_obs[], baseline_stop_obs[]
        start_val, stop_val = _parse_baseline_values(start_str, stop_str)

        # Convert to tuple if valid (baseline! accepts tuples and converts internally)
        baseline_interval_new = (!isnothing(start_val) && !isnothing(stop_val)) ? (start_val, stop_val) : nothing

        # Check if baseline actually changed
        baseline_changed = baseline_interval_new !== previous_baseline[]

        # Apply baseline if it changed
        if baseline_changed && !isnothing(baseline_interval_new)
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
