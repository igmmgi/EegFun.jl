# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_EPOCHS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Data processing
    :average_channels => (false, "Whether to average across channels"),
    :plot_avg_trials => (true, "Whether to draw ERP average overlay"),

    # Axis limits, labels, and direction
    :title => (nothing, "Plot title. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :yreversed => (false, "Whether to reverse the y-axis"),

    # Line styling
    :linewidth => (1, "Line width for epoch traces"),
    :avg_linewidth_multiplier => (2.0, "Multiplier for average line width."),
    :trial_alpha => (0.25, "Alpha (transparency) for individual trial traces"),
    :color => (:black, "Color for epoch traces (can be a single color or a vector of colors, one per condition)"),
    :colormap => (:jet, "Colormap for multi-condition plots"),

    # Layout configuration
    :layout => (:single, "Layout type: :single, :grid, or :topo"),

    # Layout parameters (for topo and other layouts)
    :layout_topo_plot_width => (0.10, "Width of individual plots (fraction of figure width)"),
    :layout_topo_plot_height => (0.10, "Height of individual plots (fraction of figure height)"),
    :layout_topo_margin => (0.12, "Margin between plots"),
    :layout_topo_scale_offset => (0.1, "Offset factor for scale plot position"),
    :layout_topo_scale_pos => ((0.8, -0.8), "Fallback position for scale plot in topo layout as (x, y) tuple"),

    # Grid layout parameters
    :layout_grid_rowgap => (10, "Gap between rows (in pixels)"),
    :layout_grid_colgap => (10, "Gap between columns (in pixels)"),
    :layout_grid_dims =>
        (nothing, "Grid dimensions as (rows, cols) tuple for grid layouts. If nothing, automatically determined"),
    :layout_grid_skip_positions =>
        (nothing, "Positions to skip in grid layout as vector of (row, col) tuples, e.g., [(2,1), (2,3)]"),

    # Display options
    :theme_fontsize => (24, "Font size for theme"),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Origin lines
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),

    # Interactive features
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
    :legend_labels => ([], "If plotting multiple plots, custom labels for conditions."),
    :legend_nbanks => (nothing, "Number of columns for the legend. If nothing, automatically determined."),

    # General layout parameters
    :figure_padding =>
        ((10, 10, 10, 10), "Padding around entire figure as (left, right, top, bottom) tuple (in pixels)"),
)

"""
    plot_epochs(filepath::String; 
               channel_selection::Function = channels(),
               sample_selection::Function = samples(), 
               epoch_selection::Function = epochs(),
               include_extra::Bool = false,
               layout = :single,
               kwargs...)

Load epoch data from a JLD2 file and create plots.

# Arguments
- `filepath::String`: Path to JLD2 file containing EpochData
- `channel_selection::Function`: Function that returns boolean vector for channel filtering
- `sample_selection::Function`: Function that returns boolean vector for sample filtering
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering
- `include_extra::Bool`: Whether to include extra channels
- `layout`: Layout specification (see main plot_epochs documentation)
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Load and plot from file
plot_epochs("Flank_C_3_epochs_original.jld2")

# With channel selection
plot_epochs("Flank_C_3_epochs_original.jld2", channel_selection = channels([:PO7, :PO8]), layout = :grid)
```
"""
function plot_epochs(
    filepath::String;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)

    data = load_data(filepath)
    isnothing(data) && @minimal_error_throw "No data found in file: $filepath"

    # Dispatch to main plot_epochs function (handles both EpochData and Vector{EpochData})
    return plot_epochs(
        data;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
        layout = layout,
        kwargs...,
    )
end

"""
    plot_epochs(datasets::Vector{EpochData}; 
                condition_selection::Function = conditions(),
                channel_selection::Function = channels(),
                sample_selection::Function = samples(), 
                epoch_selection::Function = epochs(),
                include_extra::Bool = false,
                layout = :single,
                kwargs...)

Plot multiple epoch datasets (conditions) with flexible layout options.
Conditions are overlaid on the same plot with different colors.

# Arguments
- `datasets::Vector{EpochData}`: Vector of epoch data structures (one per condition)
- `condition_selection::Function`: Function that returns boolean vector for condition filtering (default: `conditions()`)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: `channels()`)
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: `samples()`)
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering (default: `epochs()`)
- `include_extra::Bool`: Whether to include extra channels (default: `false`)
- `layout`: Layout specification (see single EpochData method documentation)
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot all conditions
fig, ax = plot_epochs([dat1, dat2])

# Plot specific conditions
fig, ax = plot_epochs([dat1, dat2], condition_selection = conditions([1, 2]))
```
"""
function plot_epochs(
    datasets::Vector{EpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)::Tuple{Figure,Union{Axis,Vector{Axis}}}

    user_provided_color = haskey(kwargs, :color)
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)

    # Use subset to filter by condition, sample, and epoch - NOT by channel
    # Channel selection is applied later at plot time to keep all channels available for topoplots
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channels(),  # Select ALL channels (no filtering)
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )

    # Generate window title from datasets
    title_str = _generate_window_title(dat_subset)
    set_window_title(title_str)

    # Check if subsetting resulted in empty data
    isempty(dat_subset) && @minimal_error_throw "No data matched the selection criteria!"

    # Get ALL channels from first dataset (not filtered by channel_selection)
    all_channels = channel_labels(dat_subset[1])
    extra_channels = extra_labels(dat_subset[1])
    all_channels = vcat(all_channels, extra_channels)

    # Apply channel_selection to determine which channels to plot
    selected_channels =
        get_selected_channels(first(dat_subset), channel_selection; include_meta = false, include_extra = include_extra)
    # Preserve order from selected_channels (user's channel_selection order)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    # Validate we have channels to plot
    isempty(all_plot_channels) && @minimal_error_throw "No channels selected for plotting"

    # If average_channels is requested, apply channel_average! directly to dat_subset (mutate in place)
    if plot_kwargs[:average_channels]
        for dat in dat_subset
            channel_average!(
                dat;
                channel_selections = [channels(all_plot_channels)],
                output_labels = [:avg],
                reduce = false,
            )
        end
    end

    # Create dat_subset_avg: average over trials (ERP data) for the average overlay
    dat_subset_avg = [average_epochs(dat) for dat in dat_subset]

    # Compute colors - need to account for number of channels when in :single layout
    n_conditions = length(dat_subset)
    n_channels_for_colors = (layout == :single && length(all_plot_channels) > 1) ? length(all_plot_channels) : 1
    condition_colors_list = _compute_dataset_colors(
        plot_kwargs[:color],
        n_conditions,
        n_channels_for_colors,
        plot_kwargs[:colormap],
        user_provided_color,
    )

    @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(n_conditions) conditions"

    # Extract layout_* parameters, remove prefix, and pass to create_layout
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    # Create plot_layout object for unified selection system
    plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout; layout_kwargs...)

    # Create figure with padding (guaranteed to be in plot_kwargs)
    fig = Figure(figure_padding = plot_kwargs[:figure_padding])
    axes = Axis[]

    # Apply theme font size early
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])

    # Initialize line references for control panel if interactive
    line_refs = nothing
    if plot_kwargs[:interactive]
        n_axes = (layout == :grid || layout == :topo) ? length(all_plot_channels) : 1
        line_refs = [Dict{Int,Dict{Symbol,Any}}() for _ = 1:n_axes]
    end

    if plot_kwargs[:average_channels]
        ax = Axis(fig[1, 1])
        push!(axes, ax)
        ax_line_refs = plot_kwargs[:interactive] ? line_refs[1] : nothing

        # For each condition, plot using dat_subset 
        for (cond_idx, dat) in enumerate(dat_subset)
            cond_plot_kwargs = merge(plot_kwargs, Dict(:color => condition_colors_list[cond_idx]))
            label = _get_condition_label(dat, cond_idx, plot_kwargs[:legend_labels])

            trial_line, trial_y_obs =
                _plot_epochs!(ax, dat, [:avg], cond_plot_kwargs; label = label, line_refs = ax_line_refs)

            if plot_kwargs[:interactive]
                _store_line_ref!(ax_line_refs, cond_idx, :avg, trial_line, trial_y_obs)
            end

            if plot_kwargs[:plot_avg_trials]
                erp_dat = dat_subset_avg[cond_idx]
                avg_label = label * " (avg)"
                avg_line, y_obs = _plot_erp_average!(
                    ax,
                    erp_dat,
                    [:avg],
                    cond_plot_kwargs;
                    label = avg_label,
                    line_refs = ax_line_refs,
                )
                if plot_kwargs[:interactive]
                    _store_avg_line_ref!(ax_line_refs, cond_idx, :avg, avg_line, y_obs)
                end
            end
        end

        # Add legend
        _add_epochs_legend!(ax, [:avg], dat_subset, plot_kwargs)

        # Set axis properties
        if length(all_plot_channels) == 1
            ax.title = string(all_plot_channels[1])
        else
            ax.title = "Avg: $(print_vector(all_plot_channels, max_length = 8, n_ends = 3))"
        end

        _set_axis_properties!(
            ax;
            xlim = plot_kwargs[:xlim],
            ylim = plot_kwargs[:ylim],
            xlabel = plot_kwargs[:xlabel],
            ylabel = plot_kwargs[:ylabel],
            yreversed = plot_kwargs[:yreversed],
        )
        _set_axis_grid!(
            ax;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )
        _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])

    else # Multi-channel layout (single, grid, topo)

        # Determine layout type
        if layout == :single
            ax = Axis(fig[1, 1])
            push!(axes, ax)
            ax_line_refs = plot_kwargs[:interactive] ? line_refs[1] : nothing

            # Plot all conditions for each channel
            for (ch_idx, ch) in enumerate(all_plot_channels)
                for (cond_idx, dat) in enumerate(dat_subset)
                    channel_suffix = length(all_plot_channels) > 1 ? " ($ch)" : ""
                    label = _get_condition_label(dat, cond_idx, plot_kwargs[:legend_labels], channel_suffix)

                    color_idx =
                        length(all_plot_channels) > 1 ? ((cond_idx - 1) * length(all_plot_channels) + ch_idx) : cond_idx
                    cond_plot_kwargs = merge(plot_kwargs, Dict(:color => condition_colors_list[color_idx]))

                    trial_line, trial_y_obs =
                        _plot_epochs!(ax, dat, [ch], cond_plot_kwargs; label = label, line_refs = ax_line_refs)

                    if plot_kwargs[:interactive]
                        _store_line_ref!(ax_line_refs, cond_idx, ch, trial_line, trial_y_obs)
                    end

                    if plot_kwargs[:plot_avg_trials]
                        erp_dat = average_epochs(dat)
                        avg_label = label * " (avg)"
                        avg_line, y_obs = _plot_erp_average!(
                            ax,
                            erp_dat,
                            [ch],
                            cond_plot_kwargs;
                            label = avg_label,
                            line_refs = ax_line_refs,
                        )
                        if plot_kwargs[:interactive]
                            _store_avg_line_ref!(ax_line_refs, cond_idx, ch, avg_line, y_obs)
                        end
                    end
                end
            end

            # Add legend
            _add_epochs_legend!(ax, all_plot_channels, dat_subset, plot_kwargs)

            # Set axis properties
            _set_axis_properties!(
                ax;
                xlim = plot_kwargs[:xlim],
                ylim = plot_kwargs[:ylim],
                xlabel = plot_kwargs[:xlabel],
                ylabel = plot_kwargs[:ylabel],
                yreversed = plot_kwargs[:yreversed],
            )
            _set_axis_grid!(
                ax;
                xgrid = plot_kwargs[:xgrid],
                ygrid = plot_kwargs[:ygrid],
                xminorgrid = plot_kwargs[:xminorgrid],
                yminorgrid = plot_kwargs[:yminorgrid],
            )
            _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])

        elseif layout == :grid
            # Use layout_grid_dims if provided, otherwise calculate best rectangle
            grid_dims =
                plot_kwargs[:layout_grid_dims] !== nothing ? plot_kwargs[:layout_grid_dims] :
                best_rect(length(all_plot_channels))
            _plot_epochs_grid_multi!(
                fig,
                axes,
                dat_subset,
                all_plot_channels,
                grid_dims,
                condition_colors_list,
                plot_kwargs,
                line_refs,
            )

        elseif layout == :topo
            _plot_epochs_topo_multi!(
                fig,
                axes,
                dat_subset,
                all_plot_channels,
                condition_colors_list,
                plot_kwargs,
                line_refs,
            )
        end
    end

    # Setup interactivity if requested
    if plot_kwargs[:interactive]
        _setup_shared_interactivity!(fig, axes, :epochs)

        # Disable default interactions that conflict with our custom selection (all axes)
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end

        # Set up selection system for all axes (will work with linked axes)
        selection_state = SharedSelectionState(axes)

        # Set up control panel (press 'c' to open) - must be before selection to capture condition_checked
        condition_checked_ref = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing)
        _setup_epochs_control_panel!(
            fig,
            dat_subset,
            dat_subset_avg,
            axes,
            get(plot_kwargs, :baseline_interval, nothing),
            line_refs,
            condition_checked_ref,
        )

        # Create right-click handler that has access to condition visibility
        right_click_handler =
            (selection_state, mouse_x, data) ->
                _handle_epochs_right_click!(selection_state, mouse_x, data, condition_checked_ref)

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

    set_window_title("Makie")
    return fig, axes
end

"""
    plot_epochs(dat::EpochData; 
                channel_selection::Function = channels(),
                sample_selection::Function = samples(), 
                epoch_selection::Function = epochs(),
                include_extra::Bool = false,
                layout = :single,
                kwargs...)

Plot epoch data with flexible layout options.

# Arguments
- `dat::EpochData`: Epoch data structure containing multiple trials
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: `channels()`)
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: `samples()`)
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering (default: `epochs()`)
- `include_extra::Bool`: Whether to include extra channels (default: `false`)
- `layout`: Layout specification:
  - `:single` (default): Single plot with all channels
  - `:grid`: Auto-calculated grid layout
  - `:topo`: Topographic layout based on channel positions
  - `Vector{Int}`: Custom grid dimensions [rows, cols]

$(generate_kwargs_doc(PLOT_EPOCHS_KWARGS))

# Returns
- `Figure`: The Makie Figure object
- `Union{Axis, Vector{Axis}}`: Single axis for single layout, or vector of axes for grid/topo layouts

# Examples
```julia
# Single plot with all channels
fig, ax = plot_epochs(dat)

# Grid layout
fig, axes = plot_epochs(dat, layout = :grid)

# Custom grid dimensions
fig, axes = plot_epochs(dat, layout = [2, 3])

# Don't display plot
fig, ax = plot_epochs(dat; display_plot = false)

# Custom styling
fig, ax = plot_epochs(dat; 
    color = [:blue, :red],
    linewidth = [1, 3],
    title = "Custom Epochs"
)
```
"""
function plot_epochs(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)::Tuple{Figure,Union{Axis,Vector{Axis}}}
    return plot_epochs(
        [dat];
        condition_selection = conditions(),  # Always select all (just the one condition)
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
        layout = layout,
        kwargs...,
    )
end

function _plot_epochs!(
    ax,
    dat,
    channels,
    plot_kwargs;
    label::Union{String,Nothing} = nothing,
    line_refs = nothing,
)::Tuple{Lines,Observable}

    @info "plot_epochs: $(print_vector(channels))"
    @assert length(channels) == 1 "_plot_epochs! expects a single channel"

    # Cache time vector and styles
    time_vec = dat.data[1][!, :time]
    # Color is passed as single value or from color cycle, not array
    trial_color = plot_kwargs[:color] isa Vector ? plot_kwargs[:color][1] : plot_kwargs[:color]
    trial_linewidth = plot_kwargs[:linewidth] isa Vector ? plot_kwargs[:linewidth][1] : plot_kwargs[:linewidth]

    # Concatenate trials initially
    y_cat = _concatenate_trials(dat.data, channels[1], time_vec)
    time_cat = similar(y_cat)
    pos = 1
    @inbounds for t = 1:length(dat.data)
        for i = 1:length(time_vec)
            time_cat[pos] = time_vec[i]
            pos += 1
        end
        if t != length(dat.data)
            time_cat[pos] = NaN
            pos += 1
        end
    end

    # Use Observable for y-data to allow updates (e.g., for baseline changes)
    y_obs = Observable(y_cat)

    line = lines!(
        ax,
        time_cat,
        y_obs,
        color = trial_color,
        linewidth = trial_linewidth,
        alpha = plot_kwargs[:trial_alpha],
        label = label,  # Should be nothing for trial lines, label for average lines
    )

    return line, y_obs
end


function _plot_erp_average!(
    ax,
    erp_dat::ErpData,
    channels::Vector{Symbol},
    plot_kwargs;
    label::Union{String,Nothing} = nothing,
    line_refs = nothing,
)::Tuple{Lines,Observable}
    @assert length(channels) == 1 "_plot_erp_average! expects a single channel"
    ch = channels[1]

    time_vec = erp_dat.data[!, :time]
    # Color is passed as single value or from color cycle, not array
    avg_color = plot_kwargs[:color] isa Vector ? plot_kwargs[:color][1] : plot_kwargs[:color]
    base_linewidth = plot_kwargs[:linewidth] isa Vector ? plot_kwargs[:linewidth][1] : plot_kwargs[:linewidth]
    avg_linewidth = base_linewidth * plot_kwargs[:avg_linewidth_multiplier]

    # Use Observable for y-data to allow updates for baseline changes
    y_obs = Observable(erp_dat.data[!, ch])
    line = lines!(ax, time_vec, y_obs, color = avg_color, linewidth = avg_linewidth, alpha = 1.0, label = label)
    return line, y_obs
end

"""
    _plot_epochs_grid_multi!(fig, axes, datasets, all_plot_channels, grid_dims, condition_colors, plot_kwargs)

Create a grid layout for plotting epochs from multiple conditions.
"""
function _plot_epochs_grid_multi!(
    fig::Figure,
    axes::Vector{Axis},
    datasets::Vector{EpochData},
    all_plot_channels::Vector{Symbol},
    grid_dims::Tuple{Int,Int},
    condition_colors::Vector,
    plot_kwargs::Dict,
    line_refs = nothing,
)
    rows, cols = grid_dims

    # Calculate y-range if not provided (use first dataset as reference)
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(datasets[1]; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    for (idx, channel) in enumerate(all_plot_channels)
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)

        # Get line_refs for this axis
        ax_line_refs = line_refs !== nothing && idx <= length(line_refs) ? line_refs[idx] : nothing

        # Plot all conditions for this channel
        for (cond_idx, dat) in enumerate(datasets)
            label = _get_condition_label(dat, cond_idx, plot_kwargs[:legend_labels])
            cond_plot_kwargs = merge(plot_kwargs, Dict(:color => condition_colors[cond_idx]))
            trial_line, trial_y_obs =
                _plot_epochs!(ax, dat, [channel], cond_plot_kwargs; label = label, line_refs = ax_line_refs)

            if ax_line_refs !== nothing
                _store_line_ref!(ax_line_refs, cond_idx, channel, trial_line, trial_y_obs)
            end

            if plot_kwargs[:plot_avg_trials]
                erp_dat = average_epochs(dat)
                avg_label = label * " (avg)"
                avg_line, y_obs = _plot_erp_average!(
                    ax,
                    erp_dat,
                    [channel],
                    cond_plot_kwargs;
                    label = avg_label,
                    line_refs = ax_line_refs,
                )
                if ax_line_refs !== nothing
                    _store_avg_line_ref!(ax_line_refs, cond_idx, channel, avg_line, y_obs)
                end
            end
        end

        # Add legend for this channel
        _add_epochs_legend!(ax, [channel], datasets, plot_kwargs)

        # Set axis properties with ylim
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim))

        # Only add x and y labels to outer left column and bottom row
        if col != 1
            axis_kwargs = merge(axis_kwargs, Dict(:ylabel => ""))
            ax.yticklabelsvisible = false
        end
        if row != rows
            axis_kwargs = merge(axis_kwargs, Dict(:xlabel => ""))
            ax.xticklabelsvisible = false
        end

        ax.title = string(channel)
        _set_axis_properties!(
            ax;
            xlim = axis_kwargs[:xlim],
            ylim = axis_kwargs[:ylim],
            xlabel = axis_kwargs[:xlabel],
            ylabel = axis_kwargs[:ylabel],
            yreversed = axis_kwargs[:yreversed],
        )
        _set_axis_grid!(
            ax;
            xgrid = axis_kwargs[:xgrid],
            ygrid = axis_kwargs[:ygrid],
            xminorgrid = axis_kwargs[:xminorgrid],
            yminorgrid = axis_kwargs[:yminorgrid],
        )
        _set_origin_lines!(ax; add_xy_origin = axis_kwargs[:add_xy_origin])
    end

    # Link axes for synchronized zooming
    if length(axes) > 1
        Makie.linkaxes!(axes...)
    end
end

"""
    _plot_epochs_topo_multi!(fig, axes, datasets, all_plot_channels, condition_colors, plot_kwargs)

Create a topographic layout for plotting epochs from multiple conditions.
"""
function _plot_epochs_topo_multi!(
    fig::Figure,
    axes::Vector{Axis},
    datasets::Vector{EpochData},
    all_plot_channels::Vector{Symbol},
    condition_colors::Vector,
    plot_kwargs::Dict,
    line_refs = nothing,
)
    # Use first dataset for layout (all should have same layout after subsetting)
    dat = datasets[1]

    # Ensure 2D coordinates exist
    if !all(in.([:x2, :y2], Ref(propertynames(dat.layout.data))))
        polar_to_cartesian_xy!(dat.layout)
    end

    # Determine global y-lims for consistency across small axes
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(dat; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    # Normalize positions to [0,1]
    x2 = dat.layout.data.x2
    y2 = dat.layout.data.y2
    minx, maxx = extrema(x2)
    miny, maxy = extrema(y2)
    xrange = maxx - minx
    yrange = maxy - miny
    xrange = xrange == 0 ? 1.0 : xrange
    yrange = yrange == 0 ? 1.0 : yrange

    plot_w = plot_kwargs[:layout_topo_plot_width]
    plot_h = plot_kwargs[:layout_topo_plot_height]
    margin = plot_kwargs[:layout_topo_margin]

    # Map channel -> position
    pos_map = Dict{Symbol,Tuple{Float64,Float64}}()
    for (lab, x, y) in zip(dat.layout.data.label, x2, y2)
        nx = (x - minx) / xrange
        ny = (y - miny) / yrange
        pos_map[Symbol(lab)] = (nx, ny)
    end

    # Create axes at positions
    for (ch_idx, ch) in enumerate(all_plot_channels)
        pos = get(pos_map, ch, (0.5, 0.5))
        ax = Axis(
            fig[1, 1],
            width = Relative(plot_w),
            height = Relative(plot_h),
            halign = clamp(pos[1], margin, 1 - margin),
            valign = clamp(pos[2], margin, 1 - margin),
        )
        push!(axes, ax)

        # Get line_refs for this axis
        ax_line_refs = line_refs !== nothing && ch_idx <= length(line_refs) ? line_refs[ch_idx] : nothing

        # Plot all conditions for this channel
        for (cond_idx, dat_cond) in enumerate(datasets)
            label = _get_condition_label(dat_cond, cond_idx, plot_kwargs[:legend_labels])
            cond_plot_kwargs = merge(plot_kwargs, Dict(:color => condition_colors[cond_idx]))
            trial_line, trial_y_obs =
                _plot_epochs!(ax, dat_cond, [ch], cond_plot_kwargs; label = label, line_refs = ax_line_refs)

            if ax_line_refs !== nothing
                _store_line_ref!(ax_line_refs, cond_idx, ch, trial_line, trial_y_obs)
            end

            if plot_kwargs[:plot_avg_trials]
                erp_dat = average_epochs(dat_cond)
                avg_label = label * " (avg)"
                avg_line, y_obs =
                    _plot_erp_average!(ax, erp_dat, [ch], cond_plot_kwargs; label = avg_label, line_refs = ax_line_refs)
                if ax_line_refs !== nothing
                    _store_avg_line_ref!(ax_line_refs, cond_idx, ch, avg_line, y_obs)
                end
            end
        end

        # Add legend for this channel
        _add_epochs_legend!(ax, [ch], datasets, plot_kwargs)

        # Suppress axis labels on all but the final axis
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim, :xlabel => "", :ylabel => ""))
        ax.title = string(ch)
        _set_axis_properties!(
            ax;
            xlim = axis_kwargs[:xlim],
            ylim = axis_kwargs[:ylim],
            xlabel = axis_kwargs[:xlabel],
            ylabel = axis_kwargs[:ylabel],
            yreversed = axis_kwargs[:yreversed],
        )
        _set_axis_grid!(
            ax;
            xgrid = axis_kwargs[:xgrid],
            ygrid = axis_kwargs[:ygrid],
            xminorgrid = axis_kwargs[:xminorgrid],
            yminorgrid = axis_kwargs[:yminorgrid],
        )
        _set_origin_lines!(ax; add_xy_origin = axis_kwargs[:add_xy_origin])
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticksvisible = false
        hidespines!(ax)
    end

    # Optional extra scale axis in bottom-right
    # Note: This is a custom feature for plot_epochs, not present in plot_erp
    # Using layout_topo_scale_pos for consistency with the layout system
    if haskey(plot_kwargs, :layout_show_scale) && plot_kwargs[:layout_show_scale]
        scale_pos = plot_kwargs[:layout_topo_scale_pos]
        scale_ax = Axis(
            fig[1, 1],
            width = Relative(plot_kwargs[:layout_topo_plot_width]),
            height = Relative(plot_kwargs[:layout_topo_plot_height]),
            halign = scale_pos[1],
            valign = scale_pos[2],
        )
        push!(axes, scale_ax)
        # No data in this axis; just show labels and limits
        tmin, tmax = (dat.data[1].time[1], dat.data[1].time[end])
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim, :xlim => (tmin, tmax)))
        scale_ax.title = ""
        _set_axis_properties!(
            scale_ax;
            xlim = axis_kwargs[:xlim],
            ylim = axis_kwargs[:ylim],
            xlabel = axis_kwargs[:xlabel],
            ylabel = axis_kwargs[:ylabel],
            yreversed = axis_kwargs[:yreversed],
        )
        _set_axis_grid!(
            scale_ax;
            xgrid = axis_kwargs[:xgrid],
            ygrid = axis_kwargs[:ygrid],
            xminorgrid = axis_kwargs[:xminorgrid],
            yminorgrid = axis_kwargs[:yminorgrid],
        )
        _set_origin_lines!(scale_ax; add_xy_origin = axis_kwargs[:add_xy_origin])
    end
end

"""
    _add_epochs_legend!(ax::Axis, channels::Vector{Symbol}, datasets::Vector{EpochData}, kwargs::Dict)

Add legend to axis for epochs plot if conditions are met.
Similar to _add_legend! for ERPs but adapted for EpochData.
"""
function _add_epochs_legend!(ax::Axis, channels::Vector{Symbol}, datasets::Vector{EpochData}, kwargs::Dict)
    # Check if legend should be shown
    # Do not show if requested false, or single channel + single dataset
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

# Helper functions for plot_epochs

"""
    _get_condition_label(dat, cond_idx, legend_labels, channel_suffix = "")

Get the label for a condition, with optional channel suffix.
"""
function _get_condition_label(dat, cond_idx, legend_labels, channel_suffix = "")
    label =
        isempty(legend_labels) ? dat.condition_name :
        (cond_idx <= length(legend_labels) ? legend_labels[cond_idx] : dat.condition_name)
    return channel_suffix == "" ? label : label * channel_suffix
end

"""
    _store_line_ref!(ax_line_refs, cond_idx, ch, trial_line, trial_y_obs)

Store trial line reference in the line_refs structure.
"""
function _store_line_ref!(ax_line_refs, cond_idx, ch, trial_line, trial_y_obs)
    if !haskey(ax_line_refs, cond_idx)
        ax_line_refs[cond_idx] = Dict{Symbol,Any}()
    end
    ax_line_refs[cond_idx][ch] = Dict{Symbol,Any}(:trials => (trial_line, trial_y_obs), :average => nothing)
end

"""
    _store_avg_line_ref!(ax_line_refs, cond_idx, ch, avg_line, y_obs)

Store average line reference in the line_refs structure.
"""
function _store_avg_line_ref!(ax_line_refs, cond_idx, ch, avg_line, y_obs)
    ax_line_refs[cond_idx][ch][:average] = (avg_line, y_obs)
end

# Helper to concatenate trials with NaN separators
function _concatenate_trials(trials, ch, time_vec)
    m = length(trials)
    n = length(time_vec)
    total_len = m * n + (m - 1)
    y_cat = Vector{Float64}(undef, total_len)
    pos = 1
    @inbounds for t = 1:m
        df = trials[t]
        y = df[!, ch]
        for i = 1:n
            y_cat[pos] = y[i]
            pos += 1
        end
        if t != m
            y_cat[pos] = NaN
            pos += 1
        end
    end
    return y_cat
end

"""
    _setup_linked_legend_interactions_epochs!(line_refs::Vector{<:Dict})

Set up linked legend interactions so clicking a legend entry in one plot
toggles visibility of the corresponding condition in all plots.
Structure: line_refs[ax_idx][condition_idx][ch][:trials] and [:average]

Note: Both trial and average lines appear in the legend (trial with condition name, 
average with condition name + " (avg)"). We sync both types across plots.
"""
function _setup_linked_legend_interactions_epochs!(line_refs::Vector{<:Dict})
    # Create separate mappings: condition_idx -> trial lines and condition_idx -> average lines
    condition_trial_lines = Dict{Int,Vector{Any}}()
    condition_avg_lines = Dict{Int,Vector{Any}}()

    # Collect all lines for each condition, separating trials and averages
    # Structure: line_refs[ax_idx][condition_idx][ch][:trials] and [:average]
    for ax_line_refs in line_refs
        for (cond_idx, cond_line_data) in ax_line_refs
            # Iterate over channels
            for (ch, ch_line_data) in cond_line_data
                if !isa(ch_line_data, Dict)
                    continue
                end
                # Collect trial lines
                if haskey(ch_line_data, :trials) && ch_line_data[:trials] !== nothing
                    trial_data = ch_line_data[:trials]
                    trial_line = nothing
                    if trial_data isa Tuple && length(trial_data) == 2
                        trial_line = trial_data[1]  # Get the line, not the tuple
                    elseif trial_data isa Lines
                        trial_line = trial_data
                    end
                    if trial_line !== nothing
                        trial_lines = get!(condition_trial_lines, cond_idx, Any[])
                        push!(trial_lines, trial_line)
                    end
                end

                # Collect average lines
                if haskey(ch_line_data, :average) && ch_line_data[:average] !== nothing
                    avg_data = ch_line_data[:average]
                    avg_line = nothing
                    if avg_data isa Tuple && length(avg_data) == 2
                        avg_line = avg_data[1]  # Get the line, not the tuple
                    elseif avg_data isa Lines
                        avg_line = avg_data
                    end
                    if avg_line !== nothing
                        avg_lines = get!(condition_avg_lines, cond_idx, Any[])
                        push!(avg_lines, avg_line)
                    end
                end
            end
        end
    end

    # Sync trial lines independently across plots
    for lines in values(condition_trial_lines)
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

    # Sync average lines independently across plots
    for lines in values(condition_avg_lines)
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
    _setup_epochs_control_panel!(fig::Figure, dat_subset::Vector{EpochData}, dat_subset_avg::Vector{ErpData}, axes::Vector{Axis}, 
                                 baseline_interval::BaselineInterval,
                                 line_refs::Vector{<:Dict})

Set up a control panel that opens when 'c' key is pressed.
Allows adjusting baseline and toggling conditions.
"""
function _setup_epochs_control_panel!(
    fig::Figure,
    dat_subset::Vector{EpochData},
    dat_subset_avg::Vector{ErpData},
    axes::Vector{Axis},
    baseline_interval::BaselineInterval,
    line_refs::Vector{<:Dict},
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}} = Ref{Union{Vector{Observable{Bool}},Nothing}}(
        nothing,
    ),
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)

    # Set up linked legend interactions
    _setup_linked_legend_interactions_epochs!(line_refs)

    # State: baseline values and condition selections
    start_val, stop_val = _extract_baseline_values(baseline_interval)
    baseline_start_obs = Observable(start_val === nothing ? "" : string(start_val))
    baseline_stop_obs = Observable(stop_val === nothing ? "" : string(stop_val))
    condition_checked = [Observable(true) for _ in dat_subset]
    condition_checked_ref[] = condition_checked  # Store for access by right-click handler

    # Track previous baseline to avoid unnecessary updates
    previous_baseline = Ref{Union{Tuple{Float64,Float64},Nothing}}(nothing)

    # Update plot (toggle visibility instead of re-plotting)
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
            baseline!.(dat_subset_avg, Ref(baseline_interval_new))
            previous_baseline[] = baseline_interval_new
        end

        # Build condition mask
        condition_mask = [checked[] for checked in condition_checked]

        # Update line visibility and y-data
        # Structure: line_refs[ax_idx][cond_idx][ch][:trials] and [:average]
        for (ax_idx, ax_line_refs) in enumerate(line_refs)
            ax_idx > length(axes) && continue
            for (cond_idx, cond_line_data) in ax_line_refs
                cond_idx > length(dat_subset) && continue
                visible = cond_idx <= length(condition_mask) ? condition_mask[cond_idx] : true
                dat = dat_subset[cond_idx]

                # Iterate over channels
                for (ch, line_data) in cond_line_data
                    if !isa(line_data, Dict)
                        continue  # Skip if not the expected structure
                    end

                    # Update trial line visibility and y-data
                    if haskey(line_data, :trials) && line_data[:trials] !== nothing
                        trial_line, trial_y_obs = line_data[:trials]
                        trial_line.visible = visible
                        if baseline_changed
                            time_vec = dat.data[1][!, :time]
                            trial_y_obs[] = _concatenate_trials(dat.data, ch, time_vec)
                        end
                    end

                    # Update average line visibility and y-data (if present)
                    if haskey(line_data, :average) && line_data[:average] !== nothing
                        avg_line, y_obs = line_data[:average]
                        avg_line.visible = visible
                        if baseline_changed && cond_idx <= length(dat_subset_avg)
                            y_obs[] = dat_subset_avg[cond_idx].data[!, ch]
                        end
                    end
                end
            end
        end

    end

    # Keyboard handler for 'c' key
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.c

            # Create new control panel
            control_fig[] = Figure(title = "Epochs Control Panel", size = (300, 400))
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

            # Auto-update on condition changes
            for checked in condition_checked
                on(checked) do _
                    update_plot!()
                end
            end

            display(control_fig[])

        end
    end
end

# =============================================================================
# RIGHT-CLICK HANDLERS FOR EPOCHS
# =============================================================================

function _handle_epochs_right_click!(selection_state, mouse_x, data, condition_checked_ref)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_epochs_context_menu!(selection_state, data, condition_checked_ref)
    end
end

function _show_epochs_context_menu!(selection_state, data, condition_checked_ref)

    menu_fig = Figure()

    # Filter by visible conditions to determine if we have multiple visible conditions
    data_to_plot = _filter_visible_conditions_epochs(data, condition_checked_ref)
    has_multiple_conditions = data_to_plot isa Vector{EpochData} && length(data_to_plot) > 1

    plot_types = ["Topoplot"]

    # Only add average options if multiple visible conditions
    if has_multiple_conditions
        push!(plot_types, "Topoplot (average)")
    end

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            original_data, x_min, x_max = _get_epochs_selection_bounds(selection_state, data)

            # Create time-based sample selection for the plots
            time_sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max)

            # Filter by visible conditions if condition_checked is available (already done above, but do again for consistency)
            data_to_plot = _filter_visible_conditions_epochs(original_data, condition_checked_ref)

            # Convert EpochData to ErpData for topoplot (average epochs first)
            if btn.label[] in [
                "Topoplot",
                "Topoplot (average)",
            ]
                # Convert EpochData to ErpData by averaging epochs
                erp_data =
                    data_to_plot isa Vector{EpochData} ? [average_epochs(dat) for dat in data_to_plot] :
                    [average_epochs(data_to_plot)]

                # Handle single vs vector
                erp_data = length(erp_data) == 1 ? erp_data[1] : erp_data

                if btn.label[] == "Topoplot"
                    plot_topography(erp_data, sample_selection = time_sample_selection)
                elseif btn.label[] == "Topoplot (average)"
                    avg_data = erp_data isa Vector{ErpData} ? _average_conditions(erp_data) : erp_data
                    plot_topography(avg_data, sample_selection = time_sample_selection)
                end
            end
        end
    end

    new_screen = GLMakie.Screen()
    display(new_screen, menu_fig)
end

"""
    _get_epochs_selection_bounds(selection_state, data)

Extract time bounds from selection state and return original data with bounds.
Does not subset the data - preserves all electrodes for topo plots.
Returns (data, x_min, x_max).
"""
function _get_epochs_selection_bounds(selection_state, data)
    x_min, x_max = minmax(selection_state.bounds[]...)
    return (data, x_min, x_max)
end

"""
    _filter_visible_conditions_epochs(data, condition_checked_ref)

Filter EpochData by visible conditions from the control panel.
Returns filtered Vector{EpochData} or single EpochData if only one condition.
"""
function _filter_visible_conditions_epochs(data, condition_checked_ref)
    # If no condition_checked available or data is not Vector, return as-is
    if condition_checked_ref[] === nothing || !(data isa Vector{EpochData})
        return data
    end

    condition_checked = condition_checked_ref[]
    if length(condition_checked) != length(data)
        return data  # Mismatch, return as-is
    end

    # Filter by visible conditions
    visible_data = [data[i] for i in eachindex(data) if condition_checked[i][]]

    # Return single EpochData if only one visible, otherwise Vector
    return length(visible_data) == 1 ? visible_data[1] : visible_data
end
