# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_EPOCHS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Data processing
    :average_channels => (false, "Whether to average across channels"),
    :plot_avg_trials => (true, "Whether to draw ERP average overlay"),
    :baseline_interval => (nothing, "Baseline interval as (start, stop) tuple or range in seconds"),

    # Axis limits, labels, and direction
    :plot_title => (nothing, "Plot title. If nothing, automatically determined"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :yreversed => (false, "Whether to reverse the y-axis"),

    # Line styling
    :linewidth => (nothing, "Line width for epoch traces"),
    :avg_linewidth_multiplier => (2.0, "Multiplier for average line width."),
    :trial_alpha => (0.25, "Alpha (transparency) for individual trial traces"),
    :color => (nothing, "Color for epoch traces (can be a single color or a vector of colors, one per condition)"),
    :colormap => (nothing, "Colormap for multi-condition plots"),

    # Layout configuration
    :layout => (:single, "Layout type: :single, :grid, or :topo"),

    # Layout parameters
    :layout_kwargs => (
        NamedTuple(),
        "Additional parameters to configure grid or topographical layouts (e.g., grid_dims, topo_plot_width, grid_rowgap)",
    ),

    # Display options
    :theme_fontsize => (nothing, "Font size for theme"),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Origin lines
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),
    :axis_type => (:standard, "Type of axis to draw (:standard or :origin)"),
    :scale_x_value => (nothing, "X-axis scale value/step size (e.g. 0.1 for 100 ms)"),
    :scale_y_value => (nothing, "Y-axis scale value/step size (e.g. 50.0 for 50 μV)"),

    # Interactive features
    :interactive => (true, "Whether to enable interactive features"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :selection_color => (:blue, "Color for interactive selection rectangles"),
    :selection_alpha => (0.3, "Alpha (transparency) for interactive selection rectangles"),

    # Legend parameters
    :legend_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Legend block (see Makie's [Legend documentation](https://docs.makie.org/stable/reference/blocks/legend/) for available attributes).",
    ),

    # Override specific legend parameters with custom defaults
    :legend => (true, "Whether to show the legend"),
    :legend_label => ("", "Title for the legend"),
    :legend_framevisible => (true, "Whether to show the frame of the legend"),
    :legend_position =>
        (:lt, "Position of the legend for axislegend() (symbol like :lt, :rt, :lb, :rb, or tuple like (:left, :top), or (0.5, 0.5))"),
    :legend_channel => ([], "If plotting multiple plots, within channel to put the legend on."),
    :legend_labels => ([], "If plotting multiple plots, custom labels for conditions."),
    :legend_nbanks => (nothing, "Number of columns for the legend. If nothing, automatically determined."),

    # General layout parameters
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),
)

"""
    plot_epochs(filename::String; 
               input_dir::String = pwd(),
               participant_selection::Function = participants(),
               channel_selection::Function = channels(),
               sample_selection::Function = samples(),
    interval_selection::Interval = times(), 
               epoch_selection::Function = epochs(),
               include_extra::Bool = false,
               layout = :single,
               kwargs...)

Load epoch data and create plots. Accepts either a direct `.jld2` filepath or a pattern
string to discover and plot all matching files (one plot per file).

$(_generate_kwargs_doc(PLOT_EPOCHS_KWARGS))

# Examples
```julia
plot_epochs("epochs_raw.jld2")
plot_epochs("epochs_raw")
plot_epochs("epochs_raw", participant_selection = participants(1))
```
"""
function plot_epochs(
    filename::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)
    if endswith(filename, ".jld2")
        data = read_data(filename)
        isnothing(data) && @minimal_error "No data found in file: $filename"
        return plot_epochs(
            data;
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            epoch_selection = epoch_selection,
            include_extra = include_extra,
            layout = layout,
            kwargs...,
        )
    else
        files = _find_batch_files(filename, input_dir, participant_selection)
        if isempty(files)
            all_matching = _find_batch_files(filename, input_dir)  # without participant filter
            if isempty(all_matching)
                @minimal_error "No files matching pattern '$filename' in $input_dir"
            else
                avail_ids = sort(unique(_extract_participant_id.(all_matching)))
                @minimal_error "Pattern '$filename' matched $(length(all_matching)) file(s), but none passed the participant selection. Available participant IDs: $avail_ids"
            end
        end

        results = NamedTuple[]
        for file in sort(files, by = _natural_sort_key)
            file_path = joinpath(input_dir, file)
            @info "Plotting: $file"
            data = read_data(file_path)
            isnothing(data) && continue
            result = plot_epochs(
                data;
                channel_selection = channel_selection,
                sample_selection = sample_selection,
                interval_selection = interval_selection,
                epoch_selection = epoch_selection,
                include_extra = include_extra,
                layout = layout,
                kwargs...,
            )
            push!(results, result)
        end
        return results
    end
end
function plot_epochs(
    datasets::AbstractVector{<:EpochData};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)::NamedTuple{(:fig, :axes),Tuple{Figure,Vector{Axis}}}

    user_provided_color = haskey(kwargs, :color)
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)

    # Use subset to filter by condition, sample, and epoch - NOT by channel
    # Channel selection is applied later at plot time to keep all channels available for topoplots
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channels(),  # Select ALL channels (no filtering)
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )

    # Check if subsetting resulted in empty data before generating window title
    isempty(dat_subset) && @minimal_error "No data matched the selection criteria!"

    # Apply initial baseline if requested
    if !isnothing(plot_kwargs[:baseline_interval])
        baseline!(dat_subset, plot_kwargs[:baseline_interval])
    end

    # Generate window title from datasets
    title_str = _generate_window_title(dat_subset)
    _set_window_title(title_str)

    try
        # Get ALL channels from first dataset (not filtered by channel_selection)
        all_channels = channel_labels(dat_subset[1])
        extra_channels = extra_labels(dat_subset[1])
        all_channels = vcat(all_channels, extra_channels)

        # Apply channel_selection to determine which channels to plot
        selected_channels = get_selected_channels(first(dat_subset), channel_selection; include_meta = false, include_extra = include_extra)
        # Preserve order from selected_channels (user's channel_selection order)
        all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

        # Validate we have channels to plot
        isempty(all_plot_channels) && @minimal_error "No channels selected for plotting"

        # If average_channels is requested, apply channel_average! directly to dat_subset (mutate in place)
        if plot_kwargs[:average_channels]
            for dat in dat_subset
                channel_average!(dat; channel_selections = [channels(all_plot_channels)], output_labels = [:avg], reduce = false)
            end
        end

        # Create dat_subset_avg: average over trials (ERP data) for the average overlay
        dat_subset_avg = [average_epochs(dat) for dat in dat_subset]

        # Precompute time_cat vectors once per condition for all channels
        time_cat_list = [_build_time_cat(n_epochs(dat), dat.data[1][!, :time]) for dat in dat_subset]

        # Compute colors - need to account for number of channels when in :single layout (only if not averaging channels)
        n_conditions = length(dat_subset)
        n_channels_for_colors =
            (layout == :single && length(all_plot_channels) > 1 && !plot_kwargs[:average_channels]) ? length(all_plot_channels) : 1
        color_cycle =
            haskey(Makie.current_default_theme(), :palette) ? Makie.current_default_theme()[:palette][:color][] : Makie.wong_colors()
        resolved_colormap =
            isnothing(plot_kwargs[:colormap]) ?
            (haskey(Makie.current_default_theme(), :colormap) ? Makie.current_default_theme()[:colormap][] : nothing) :
            plot_kwargs[:colormap]
        condition_colors_list = _compute_dataset_colors(
            plot_kwargs[:color],
            n_conditions,
            n_channels_for_colors,
            resolved_colormap,
            user_provided_color,
            color_cycle,
        )

        @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(n_conditions) conditions"

        # Extract layout_* parameters, remove prefix, and pass to create_layout
        layout_kwargs = pop!(plot_kwargs, :layout_kwargs, (;))

        # Create plot_layout object for unified selection system
        plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout; layout_kwargs...)

        # Create figure with padding (guaranteed to be in plot_kwargs)
        fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
        fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
        axes = Axis[]

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

            _plot_epochs_conditions_on_axis!(
                ax,
                dat_subset,
                :avg,
                condition_colors_list,
                plot_kwargs,
                ax_line_refs;
                dat_subset_avg = dat_subset_avg,
                time_cat_list = time_cat_list,
            )

            # Add legend
            _add_epochs_legend!(ax, [:avg], dat_subset, plot_kwargs)

            # Set axis properties
            if length(all_plot_channels) == 1
                ax.title = string(all_plot_channels[1])
            else
                ax.title = "Avg: $(_print_vector(all_plot_channels, max_length = 8, n_ends = 3))"
            end
            _finalize_epochs_axis!(ax, plot_kwargs)

        else # Multi-channel layout (single, grid, topo)

            # Determine layout type
            if layout == :single
                ax = Axis(fig[1, 1])
                push!(axes, ax)
                ax_line_refs = plot_kwargs[:interactive] ? line_refs[1] : nothing

                # Plot all conditions for each channel
                for (ch_idx, ch) in enumerate(all_plot_channels)
                    ch_colors = [
                        condition_colors_list[length(all_plot_channels) > 1 ? ((ci - 1) * length(all_plot_channels) + ch_idx) : ci] for
                        ci = 1:length(dat_subset)
                    ]
                    channel_suffix = length(all_plot_channels) > 1 ? " ($ch)" : ""
                    _plot_epochs_conditions_on_axis!(
                        ax,
                        dat_subset,
                        ch,
                        ch_colors,
                        plot_kwargs,
                        ax_line_refs;
                        dat_subset_avg = dat_subset_avg,
                        time_cat_list = time_cat_list,
                        channel_suffix = channel_suffix,
                    )
                end

                # Add legend
                _add_epochs_legend!(ax, all_plot_channels, dat_subset, plot_kwargs)

                # Set axis properties
                _finalize_epochs_axis!(ax, plot_kwargs)

            elseif layout == :grid
                # Use layout_grid_dims if provided, otherwise calculate best rectangle
                grid_dims = plot_layout.type == :grid ? Tuple(plot_layout.dims) : _best_rect(length(all_plot_channels))
                _plot_epochs_grid_multi!(
                    fig,
                    axes,
                    dat_subset,
                    all_plot_channels,
                    grid_dims,
                    condition_colors_list,
                    plot_kwargs,
                    line_refs;
                    dat_subset_avg = dat_subset_avg,
                    time_cat_list = time_cat_list,
                )

            elseif layout == :topo
                _plot_epochs_topo_multi!(
                    fig,
                    axes,
                    dat_subset,
                    all_plot_channels,
                    condition_colors_list,
                    plot_kwargs,
                    line_refs;
                    dat_subset_avg = dat_subset_avg,
                    time_cat_list = time_cat_list,
                    plot_layout = plot_layout,
                )
            end
        end

        # Setup interactivity if requested
        if plot_kwargs[:interactive]
            _setup_shared_interactivity!(fig, axes, :epochs; zoom_step = plot_kwargs[:zoom_step])

            # Disable default interactions that conflict with our custom selection (all axes)
            for ax in axes
                deregister_interaction!(ax, :rectanglezoom)
            end

            # Set up selection system for all axes (will work with linked axes)
            selection_state =
                SharedSelectionState(axes; selection_color = plot_kwargs[:selection_color], selection_alpha = plot_kwargs[:selection_alpha])

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
                (selection_state, mouse_x, data) -> _handle_epochs_right_click!(selection_state, mouse_x, data, condition_checked_ref)

            # Set up selection system that works for all layouts
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, right_click_handler)

            # Set up channel selection events for topo and grid layouts
            if plot_layout.type in (:topo, :grid)
                _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes)
            end
        end

        if plot_kwargs[:display_plot]
            _display_figure(fig)
        end

        return (fig = fig, axes = axes)
    finally
        _set_window_title("Makie")
    end
end
function plot_epochs(
    dat::EpochData;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)::NamedTuple{(:fig, :axes),Tuple{Figure,Union{Axis,Vector{Axis}}}}
    result = plot_epochs(
        [dat];
        condition_selection = conditions(),  # Always select all (just the one condition)
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
        layout = layout,
        kwargs...,
    )
    # Return NamedTuple to match signature
    return (fig = result.fig, axes = result.axes)
end

"""Finalize standard axis properties, grids, origin lines, and scale indicator."""
function _finalize_epochs_axis!(ax::Axis, kw::Dict)
    _set_axis_properties!(
        ax;
        xlim = kw[:xlim],
        ylim = kw[:ylim],
        xlabel = kw[:xlabel],
        ylabel = kw[:ylabel],
        yreversed = kw[:yreversed],
        scale_x_value = kw[:scale_x_value],
        scale_y_value = kw[:scale_y_value],
    )
    _set_axis_grid!(ax; xgrid = kw[:xgrid], ygrid = kw[:ygrid], xminorgrid = kw[:xminorgrid], yminorgrid = kw[:yminorgrid])
    _set_origin_lines!(ax; add_xy_origin = kw[:add_xy_origin])
    _add_origin_scale_indicator!(
        ax;
        axis_type = kw[:axis_type],
        scale_x_value = kw[:scale_x_value],
        scale_y_value = kw[:scale_y_value],
        xlabel = kw[:xlabel],
        ylabel = kw[:ylabel],
    )
end

"""Plot all conditions for a given channel on an axis, with optional average overlay and line tracking."""
function _plot_epochs_conditions_on_axis!(
    ax::Axis,
    datasets::Vector{EpochData},
    channel::Symbol,
    colors::Vector,
    plot_kwargs::Dict,
    ax_line_refs;
    dat_subset_avg::Union{Vector{ErpData},Nothing} = nothing,
    time_cat_list::Union{Vector{<:AbstractVector},Nothing} = nothing,
    channel_suffix::String = "",
)
    for (cond_idx, dat) in enumerate(datasets)
        label = _get_condition_label(dat, cond_idx, plot_kwargs[:legend_labels], channel_suffix)
        cond_plot_kwargs = merge(plot_kwargs, Dict(:color => colors[cond_idx]))
        tc = isnothing(time_cat_list) ? nothing : time_cat_list[cond_idx]

        trial_line, trial_y_obs =
            _plot_epochs!(ax, dat, [channel], cond_plot_kwargs; time_cat = tc, label = label, line_refs = ax_line_refs)

        if !isnothing(ax_line_refs)
            _store_line_ref!(ax_line_refs, cond_idx, channel, trial_line, trial_y_obs)
        end

        if plot_kwargs[:plot_avg_trials]
            erp_dat = !isnothing(dat_subset_avg) ? dat_subset_avg[cond_idx] : average_epochs(dat)
            avg_label = label * " (avg)"
            avg_line, y_obs = _plot_erp_average!(ax, erp_dat, [channel], cond_plot_kwargs; label = avg_label, line_refs = ax_line_refs)
            if !isnothing(ax_line_refs)
                _store_avg_line_ref!(ax_line_refs, cond_idx, channel, avg_line, y_obs)
            end
        end
    end
end

"""Plot single-channel epoch traces (NaN-separated trials) onto an axis."""
function _plot_epochs!(
    ax,
    dat,
    channels,
    plot_kwargs;
    time_cat::Union{AbstractVector,Nothing} = nothing,
    label::Union{String,Nothing} = nothing,
    line_refs = nothing,
)::Tuple{Lines,Observable}

    @debug "plot_epochs: $(_print_vector(channels))"
    length(channels) == 1 || @minimal_error "_plot_epochs! expects a single channel"

    # Cache time vector and styles
    time_vec = dat.data[1][!, :time]
    # Color is passed as single value or from color cycle, not array
    trial_color = plot_kwargs[:color] isa Vector ? plot_kwargs[:color][1] : plot_kwargs[:color]
    actual_lw = _resolve_theme_linewidth(ax, plot_kwargs[:linewidth], 1)
    trial_linewidth = actual_lw isa Vector ? actual_lw[1] : actual_lw

    # Concatenate trials efficiently
    tc = isnothing(time_cat) ? _build_time_cat(n_epochs(dat), time_vec) : time_cat
    y_cat = _concatenate_trials(dat.data, channels[1], time_vec)

    # Use Observable for y-data to allow updates (e.g., for baseline changes)
    y_obs = Observable(y_cat)

    line = lines!(
        ax,
        tc,
        y_obs,
        color = trial_color,
        linewidth = trial_linewidth,
        alpha = plot_kwargs[:trial_alpha],
        label = label,  # Should be nothing for trial lines, label for average lines
    )

    return line, y_obs
end


"""Overlay the trial-averaged ERP waveform on an existing axis."""
function _plot_erp_average!(
    ax,
    erp_dat::ErpData,
    channels::Vector{Symbol},
    plot_kwargs;
    label::Union{String,Nothing} = nothing,
    line_refs = nothing,
)::Tuple{Lines,Observable}
    length(channels) == 1 || @minimal_error "_plot_erp_average! expects a single channel"
    ch = channels[1]

    time_vec = erp_dat.data[!, :time]
    # Color is passed as single value or from color cycle, not array
    avg_color = plot_kwargs[:color] isa Vector ? plot_kwargs[:color][1] : plot_kwargs[:color]
    actual_lw = _resolve_theme_linewidth(ax, plot_kwargs[:linewidth], 1)
    base_linewidth = actual_lw isa Vector ? actual_lw[1] : actual_lw
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
    grid_dims::Union{Tuple{Int,Int},Vector{Int}},
    condition_colors::Vector,
    plot_kwargs::Dict,
    line_refs = nothing;
    dat_subset_avg::Union{Vector{ErpData},Nothing} = nothing,
    time_cat_list::Union{Vector{Vector{Float64}},Nothing} = nothing,
)
    rows, cols = grid_dims

    # Calculate y-range if not provided across all datasets
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        ylim = _ylimits(datasets; channel_selection = channels(all_plot_channels))
    end

    for (idx, channel) in enumerate(all_plot_channels)
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)

        # Get line_refs for this axis
        ax_line_refs = !isnothing(line_refs) && idx <= length(line_refs) ? line_refs[idx] : nothing

        # Plot all conditions for this channel
        _plot_epochs_conditions_on_axis!(
            ax,
            datasets,
            channel,
            condition_colors,
            plot_kwargs,
            ax_line_refs;
            dat_subset_avg = dat_subset_avg,
            time_cat_list = time_cat_list,
        )

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
        _finalize_epochs_axis!(ax, axis_kwargs)
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
    line_refs = nothing;
    dat_subset_avg::Union{Vector{ErpData},Nothing} = nothing,
    time_cat_list::Union{Vector{Vector{Float64}},Nothing} = nothing,
    plot_layout = nothing,
)
    # Use first dataset for layout (all should have same layout after subsetting)
    dat = datasets[1]

    # Ensure 2D coordinates exist
    if !all(in.([:x2, :y2], Ref(propertynames(dat.layout.data))))
        polar_to_cartesian_xy!(dat.layout)
    end

    # Determine global y-lims across all datasets for consistency across small axes
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        ylim = _ylimits(datasets; channel_selection = channels(all_plot_channels))
    end

    plot_w = get(plot_kwargs, :layout_topo_plot_width, !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_plot_width, 0.1) : 0.1)
    plot_h = get(plot_kwargs, :layout_topo_plot_height, !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_plot_height, 0.1) : 0.1)

    # Normalize positions to [0,1] based on ALL channels in the layout
    # This ensures positions match the layout file regardless of which channels are selected
    # This matches how _apply_layout! handles topo layouts
    x2 = dat.layout.data.x2
    y2 = dat.layout.data.y2
    minx, maxx = extrema(x2)
    miny, maxy = extrema(y2)
    xrange = maxx - minx == 0 ? 1.0 : maxx - minx
    yrange = maxy - miny == 0 ? 1.0 : maxy - miny

    # Create axes at positions
    # Normalize positions to [0,1] - no margin needed since halign/valign centers the plot at the position
    for (ch_idx, ch) in enumerate(all_plot_channels)
        idx = findfirst(==(ch), dat.layout.data.label)
        if !isnothing(idx)
            x = dat.layout.data.x2[idx]
            y = dat.layout.data.y2[idx]
            halign = (x - minx) / xrange
            valign = (y - miny) / yrange
        else
            # Channel not in layout - use center position
            halign = 0.5
            valign = 0.5
        end

        ax = Axis(fig[1, 1], width = Relative(plot_w), height = Relative(plot_h), halign = halign, valign = valign)
        push!(axes, ax)

        # Get line_refs for this axis
        ax_line_refs = !isnothing(line_refs) && ch_idx <= length(line_refs) ? line_refs[ch_idx] : nothing

        # Plot all conditions for this channel
        _plot_epochs_conditions_on_axis!(
            ax,
            datasets,
            ch,
            condition_colors,
            plot_kwargs,
            ax_line_refs;
            dat_subset_avg = dat_subset_avg,
            time_cat_list = time_cat_list,
        )

        # Add legend for this channel
        _add_epochs_legend!(ax, [ch], datasets, plot_kwargs)

        # Suppress axis labels on all but the final axis
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim, :xlabel => "", :ylabel => ""))
        ax.title = string(ch)
        _finalize_epochs_axis!(ax, axis_kwargs)
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticksvisible = false
        hidespines!(ax)
    end

    show_scale =
        get(plot_kwargs, :layout_topo_show_scale, !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_show_scale, false) : false)
    if show_scale
        @info "Showing scale axis"
        scale_pos = get(
            plot_kwargs,
            :layout_topo_scale_pos,
            !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_scale_pos, (0.8, -0.8)) : (0.8, -0.8),
        )
        scale_halign = 0.0 <= scale_pos[1] <= 1.0 ? scale_pos[1] : clamp((scale_pos[1] + 1) / 2, 0.0, 1.0)
        scale_valign = 0.0 <= scale_pos[2] <= 1.0 ? scale_pos[2] : clamp((scale_pos[2] + 1) / 2, 0.0, 1.0)
        scale_ax = Axis(fig[1, 1], width = Relative(plot_w), height = Relative(plot_h), halign = scale_halign, valign = scale_valign)
        push!(axes, scale_ax)
        # No data in this axis; just show labels and limits
        tmin, tmax = (dat.data[1].time[1], dat.data[1].time[end])
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim, :xlim => (tmin, tmax)))
        scale_ax.title = ""
        _finalize_epochs_axis!(scale_ax, axis_kwargs)
        # Scale axis should show labels and ticks (unlike channel axes which are hidden)
        scale_ax.xticklabelsvisible = true
        scale_ax.yticklabelsvisible = true
        scale_ax.xticksvisible = true
        scale_ax.yticksvisible = true
    end

    # Link axes for synchronized zooming (arrow keys control all axes together)
    if length(axes) > 1
        Makie.linkaxes!(axes...)
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
    nbanks = something(get(kwargs, :legend_nbanks, nothing), length(channels) > 10 ? cld(length(channels), 10) : 1)
    legend_kwargs = get(kwargs, :legend_kwargs, (;))

    # Add legend with position and optional label
    if isempty(legend_label)
        leg = axislegend(ax; position = legend_position, nbanks = nbanks, legend_kwargs...)
    else
        leg = axislegend(ax, legend_label; position = legend_position, nbanks = nbanks, legend_kwargs...)
    end

    return leg
end

# Helper functions for plot_epochs

"""
    _get_condition_label(dat, cond_idx, legend_labels, channel_suffix = "")

Get the label for a condition, with optional channel suffix.
"""
function _get_condition_label(dat, cond_idx, legend_labels, channel_suffix = "")
    label = isempty(legend_labels) ? dat.condition_name : (cond_idx <= length(legend_labels) ? legend_labels[cond_idx] : dat.condition_name)
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

"""Precompute NaN-separated time vector for epoched trial overlays."""
function _build_time_cat(n_trials::Int, time_vec::AbstractVector)
    n = length(time_vec)
    total_len = n_trials * n + (n_trials - 1)
    T = eltype(time_vec)
    time_cat = Vector{T}(undef, total_len)
    pos = 1
    @inbounds for t = 1:n_trials
        copyto!(time_cat, pos, time_vec, 1, n)
        pos += n
        if t != n_trials
            time_cat[pos] = NaN
            pos += 1
        end
    end
    return time_cat
end

"""Concatenate trial data for one channel into a single NaN-separated vector."""
function _concatenate_trials(trials, ch, time_vec)
    m = length(trials)
    n = length(time_vec)
    total_len = m * n + (m - 1)
    col_idx = ch isa Int ? ch : DataFrames.columnindex(trials[1], ch)
    T = eltype(trials[1][!, col_idx])
    y_cat = Vector{T}(undef, total_len)
    pos = 1
    @inbounds for t = 1:m
        y = trials[t][!, col_idx]
        copyto!(y_cat, pos, y, 1, n)
        pos += n
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
                if haskey(ch_line_data, :trials) && !isnothing(ch_line_data[:trials])
                    trial_data = ch_line_data[:trials]
                    trial_line = nothing
                    if trial_data isa Tuple && length(trial_data) == 2
                        trial_line = trial_data[1]  # Get the line, not the tuple
                    elseif trial_data isa Lines
                        trial_line = trial_data
                    end
                    if !isnothing(trial_line)
                        trial_lines = get!(condition_trial_lines, cond_idx, Any[])
                        push!(trial_lines, trial_line)
                    end
                end

                # Collect average lines
                if haskey(ch_line_data, :average) && !isnothing(ch_line_data[:average])
                    avg_data = ch_line_data[:average]
                    avg_line = nothing
                    if avg_data isa Tuple && length(avg_data) == 2
                        avg_line = avg_data[1]  # Get the line, not the tuple
                    elseif avg_data isa Lines
                        avg_line = avg_data
                    end
                    if !isnothing(avg_line)
                        avg_lines = get!(condition_avg_lines, cond_idx, Any[])
                        push!(avg_lines, avg_line)
                    end
                end
            end
        end
    end

    # Sync trial and average lines independently across plots
    for lines in values(condition_trial_lines)
        _sync_lines!(lines)
    end
    for lines in values(condition_avg_lines)
        _sync_lines!(lines)
    end
end

"""Helper to sync visibility across lines in a group."""
function _sync_lines!(lines::Vector)
    length(lines) > 1 || return
    syncing = Ref(false)
    for line in lines
        other_lines = [l for l in lines if l !== line]
        on(line.visible) do visible_val
            syncing[] && return
            syncing[] = true
            for other_line in other_lines
                if other_line.visible[] != visible_val
                    other_line.visible = visible_val
                end
            end
            syncing[] = false
        end
    end
end

"""
    _setup_epochs_control_panel!(fig::Figure, dat_subset::Vector{EpochData}, dat_subset_avg::Vector{ErpData}, axes::Vector{Axis}, 
                                 baseline_interval::Interval,
                                 line_refs::Vector{<:Dict})

Set up a control panel that opens when 'c' key is pressed.
Allows adjusting baseline and toggling conditions.
"""
function _setup_epochs_control_panel!(
    fig::Figure,
    dat_subset::Vector{EpochData},
    dat_subset_avg::Vector{ErpData},
    axes::Vector{Axis},
    baseline_interval::Interval,
    line_refs::Vector{<:Dict},
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}} = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing),
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)
    control_screen = Ref{Any}(nothing)

    # Set up linked legend interactions
    _setup_linked_legend_interactions_epochs!(line_refs)

    # State: baseline values and condition selections
    start_val, stop_val = _extract_baseline_values(baseline_interval)
    baseline_start_obs = Observable(isnothing(start_val) ? "" : string(start_val))
    baseline_stop_obs = Observable(isnothing(stop_val) ? "" : string(stop_val))
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
        baseline_interval_new = (!isnothing(start_val) && !isnothing(stop_val)) ? (start_val, stop_val) : nothing

        # Check if baseline actually changed
        baseline_changed = baseline_interval_new != previous_baseline[]

        # Apply baseline if it changed
        if baseline_changed && !isnothing(baseline_interval_new)
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
                    if haskey(line_data, :trials) && !isnothing(line_data[:trials])
                        trial_line, trial_y_obs = line_data[:trials]
                        if trial_line.visible[] != visible
                            trial_line.visible = visible
                        end
                        if baseline_changed
                            time_vec = dat.data[1][!, :time]
                            trial_y_obs[] = _concatenate_trials(dat.data, ch, time_vec)
                        end
                    end

                    # Update average line visibility and y-data (if present)
                    if haskey(line_data, :average) && !isnothing(line_data[:average])
                        avg_line, y_obs = line_data[:average]
                        if avg_line.visible[] != visible
                            avg_line.visible = visible
                        end
                        if baseline_changed && cond_idx <= length(dat_subset_avg)
                            y_obs[] = dat_subset_avg[cond_idx].data[!, ch]
                        end
                    end
                end
            end
        end
    end

    # Auto-update on condition changes (registered once outside keyboard handler)
    for checked in condition_checked
        on(checked) do _
            update_plot!()
        end
    end

    # Keyboard handler for 'c' key
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.c
            # If control panel is already open, do not spawn another
            if !isnothing(control_screen[]) && hasmethod(isopen, Tuple{typeof(control_screen[])}) && isopen(control_screen[])
                return
            end

            # Create new control panel
            control_fig[] = _eegfun_figure(title = "Epochs Control Panel", size = (300, 400))
            panel_layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)

            # Baseline section
            Label(panel_layout[1, 1], "Baseline Interval", fontsize = 14, font = :bold)
            baseline_layout = GridLayout(panel_layout[2, 1], tellwidth = false, colgap = 10)

            _create_baseline_textbox(baseline_layout, 1, "Start (s):", baseline_start_obs, " ", 100)
            _create_baseline_textbox(baseline_layout, 2, "End (s):", baseline_stop_obs, " ", 100)

            # Apply button
            apply_btn = Button(panel_layout[3, 1], label = "Apply Baseline", width = 200)
            on(apply_btn.clicks) do _
                update_plot!()
            end

            # Conditions section
            Label(panel_layout[4, 1], "Conditions", fontsize = 14, font = :bold)
            conditions_layout = GridLayout(panel_layout[5, 1], tellwidth = false, rowgap = 5)

            for (i, dat) in enumerate(dat_subset)
                cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][])
                Label(conditions_layout[i, 2], dat.condition_name)
                connect!(condition_checked[i], cb.checked)
            end

            control_screen[] = _display_in_screen(control_fig[])
        end
    end
end

# === RIGHT-CLICK HANDLERS FOR EPOCHS ===

"""Dispatch right-click to the epochs context menu if inside the selection."""
function _handle_epochs_right_click!(selection_state, mouse_x, data, condition_checked_ref)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_epochs_context_menu!(selection_state, data, condition_checked_ref)
    end
end

"""Show a context menu with topoplot options for the current epoch selection."""
function _show_epochs_context_menu!(selection_state, data, condition_checked_ref)

    menu_fig = _eegfun_figure()

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
            if btn.label[] in ["Topoplot", "Topoplot (average)"]
                # Convert EpochData to ErpData by averaging epochs
                erp_data =
                    data_to_plot isa Vector{EpochData} ? [average_epochs(dat) for dat in data_to_plot] : [average_epochs(data_to_plot)]

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

    _display_in_screen(menu_fig)
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
    if isnothing(condition_checked_ref[]) || !(data isa Vector{EpochData})
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
