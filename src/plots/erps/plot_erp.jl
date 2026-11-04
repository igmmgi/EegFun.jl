const PLOT_ERP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :window_title => ("", "Title for the OS plot window. If empty, it's generated automatically."),
    :figure_title => ("", "Title drawn at the top of the entire figure canvas"),
    :figure_title_fontsize => (24, "Font size for figure title"),
    :interactive => (true, "Enable interactive features (true/false)"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :selection_color => (:blue, "Color for interactive selection rectangles"),
    :selection_alpha => (0.3, "Alpha (transparency) for interactive selection rectangles"),
    :theme_fontsize => (nothing, "Font size for theme"),

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
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Show plot title (true/false)"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Line styling
    :linewidth => (nothing, "Line width for ERPs"),
    :color => (nothing, "Color for ERPs (single color or a vector of colors, one per dataset)"),
    :linestyle => (nothing, "Line style for ERPs (single style or a vector of styles, one per dataset)"),
    :colormap => (nothing, "Colormap for multi-channel plots"),

    # Plot configuration
    :average_channels => (false, "Average across channels (true/false)"),
    :error_bars => (:none, "Error bars type to plot: :none, :sem, :within_sem, :ci95"),

    # Legend parameters
    :legend_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Legend block (see Makie's [Legend documentation](https://docs.makie.org/stable/reference/blocks/legend/) for available attributes).",
    ),
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

    # Layout parameters
    :layout_kwargs => (
        NamedTuple(),
        "Additional parameters to configure grid or topographical layouts (e.g., grid_dims, topo_plot_width, grid_rowgap)",
    ),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),

    # Highlight regions
    :highlight_regions => (
        nothing,
        "Highlight regions as a NamedTuple or Vector of NamedTuples. Each region: (x1, x2, y1=-Inf, y2=Inf, color=:gray, alpha=0.3)",
    ),
)

function _order_plot_channels(selected_channels, all_channels, channel_plot_order)
    plot_chs = [ch for ch in selected_channels if ch in all_channels]
    if !isnothing(channel_plot_order)
        ordered = Symbol[Symbol(ch) for ch in channel_plot_order if Symbol(ch) in plot_chs]
        !isempty(ordered) && (plot_chs = ordered)
    end
    return plot_chs
end

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
plot_erp("grand_average_erps_final.jld2")
plot_erp("erps_final")
plot_erp("erps_final", participant_selection = participants(1))
plot_erp("erps_final", channel_selection = channels([:PO7, :PO8]), layout = :grid)
```
"""
function plot_erp(
    filepath::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection = conditions(),
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
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
    end

    files = _require_batch_files(filepath, input_dir, participant_selection; kind = "ERP")
    results = NamedTuple[]
    for file in files
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(
            results,
            plot_erp(
                data;
                layout = layout,
                condition_selection = condition_selection,
                channel_selection = channel_selection,
                channel_plot_order = channel_plot_order,
                sample_selection = sample_selection,
                interval_selection = interval_selection,
                baseline_interval = baseline_interval,
                kwargs...,
            ),
        )
    end
    return results
end

function plot_erp(
    dat::ErpData;
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    return plot_erp(
        [dat];
        layout = layout,
        condition_selection = conditions(),
        channel_selection = channel_selection,
        channel_plot_order = channel_plot_order,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

function plot_erp(
    datasets::AbstractVector{<:ErpData};
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection = conditions(),
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    file_ids = unique(d.file for d in datasets)
    if length(file_ids) > 1
        return NamedTuple[
            plot_erp(
                filter(d -> d.file == file_id, datasets);
                layout = layout,
                condition_selection = condition_selection,
                channel_selection = channel_selection,
                channel_plot_order = channel_plot_order,
                sample_selection = sample_selection,
                interval_selection = interval_selection,
                baseline_interval = baseline_interval,
                kwargs...,
            ) for file_id in file_ids
        ]
    end

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

    return _plot_erp_core(
        dat_subset,
        all_channels,
        channel_selection_func,
        original_channels,
        Dict{Int,Tuple{DataFrames.DataFrame,Int}}(),
        plot_kwargs,
        user_provided_color,
        layout,
        baseline_interval,
        channel_plot_order,
    )
end

"""
    plot_erp_errorbar(filepath::String; kwargs...)
    plot_erp_errorbar(dataset::ErpData; kwargs...)
    plot_erp_errorbar(datasets::AbstractVector{<:ErpData}; kwargs...)

Plot grand-averaged ERP waveforms with across-participant error ribbons (`:sem`, `:within_sem`, or `:ci95`).
"""
function plot_erp_errorbar(
    filepath::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection = conditions(),
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    files = _require_batch_files(filepath, input_dir, participant_selection; kind = "ERP")
    datasets = ErpData[]
    for file in files
        @info "Loading: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        data isa Vector{ErpData} ? append!(datasets, data) : push!(datasets, data)
    end
    isempty(datasets) && @minimal_error "No data loaded from pattern: $filepath"

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
    dat::ErpData;
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
    return plot_erp_errorbar(
        [dat];
        layout = layout,
        condition_selection = conditions(),
        channel_selection = channel_selection,
        channel_plot_order = channel_plot_order,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

function plot_erp_errorbar(
    datasets::AbstractVector{<:ErpData};
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection = conditions(),
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    kwargs...,
)
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
    selected_channels = get_selected_channels(first(dat_subset), channel_selection_func; include_meta = false, include_extra = true)
    all_plot_channels = _order_plot_channels(selected_channels, all_channels, channel_plot_order)
    isempty(all_plot_channels) &&
        @minimal_error "No valid channels found. Selected channels: $selected_channels, Available channels: $all_channels"

    if isnothing(plot_kwargs[:plot_title]) && layout == :single
        plot_kwargs[:plot_title] = if plot_kwargs[:average_channels]
            "Avg: $(_print_vector(original_channels))"
        elseif length(all_plot_channels) == 1
            string(all_plot_channels[1])
        else
            "$(_print_vector(all_plot_channels))"
        end
    end

    title_str = !isempty(plot_kwargs[:window_title]) ? plot_kwargs[:window_title] : _generate_window_title(dat_subset)
    _set_window_title(title_str)

    try
        layout_kwargs = get(plot_kwargs, :layout_kwargs, (;))
        fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
        fig = _eegfun_figure(; title = title_str, figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)

        plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout; layout_kwargs...)

        if plot_layout.type == :topo && isempty(plot_kwargs[:legend_channel]) && !isempty(all_plot_channels)
            plot_kwargs[:legend_channel] = [all_plot_channels[end]]
            plot_kwargs[:legend_position] = (5, -2)
        end
        axes, channels = _apply_layout!(fig, plot_layout; plot_kwargs...)

        is_interactive = plot_kwargs[:interactive]
        line_refs = is_interactive ? [Dict{Int,Dict{Symbol,Any}}() for _ in axes] : nothing
        legend_refs = is_interactive ? Vector{Union{Legend,Nothing}}(undef, length(axes)) : nothing

        @info "plot_erp ($layout): $(_print_vector(all_plot_channels))"

        for (ax_idx, (ax, channel)) in enumerate(zip(axes, channels))
            channels_to_plot = plot_layout.type == :single ? all_plot_channels : [channel]
            ax_line_refs = is_interactive ? line_refs[ax_idx] : nothing
            leg = _plot_erp!(
                ax,
                dat_subset,
                channels_to_plot,
                plot_kwargs;
                line_refs = ax_line_refs,
                user_provided_color = user_provided_color,
                variances_dict = variances_dict,
            )
            if is_interactive && !isnothing(legend_refs)
                legend_refs[ax_idx] = leg
            end
        end

        if !isnothing(plot_kwargs[:highlight_regions])
            for ax in axes
                _draw_highlight_regions!(ax, plot_kwargs[:highlight_regions])
            end
        end

        _apply_axis_properties!.(axes; plot_kwargs...)
        _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...)
        length(axes) > 1 && linkaxes!(axes...)

        if is_interactive
            _setup_shared_interactivity!(fig, axes, :erp; zoom_step = plot_kwargs[:zoom_step])
            for ax in axes
                deregister_interaction!(ax, :rectanglezoom)
            end

            selection_state =
                SharedSelectionState(axes; selection_color = plot_kwargs[:selection_color], selection_alpha = plot_kwargs[:selection_alpha])

            condition_checked_ref = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing)
            _setup_erp_control_panel!(fig, dat_subset, axes, baseline_interval, line_refs, condition_checked_ref, legend_refs)

            right_click_handler = (sel, mouse_x, data) -> _handle_erp_right_click!(sel, mouse_x, data, condition_checked_ref)
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, right_click_handler)

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

        if !isempty(plot_kwargs[:figure_title])
            Label(fig[0, :], plot_kwargs[:figure_title], fontsize = plot_kwargs[:figure_title_fontsize], font = :bold, tellwidth = false)
        end

        plot_kwargs[:display_plot] && _display_figure(fig)
        return (fig = fig, axes = axes, line_refs = line_refs)
    finally
        _set_window_title("Makie")
    end
end

function _plot_erp_on_existing_axis!(ax::Axis, datasets::AbstractVector{<:ErpData}; dynamic_grand_avg::Bool = false, kwargs...)
    baseline_interval = get(kwargs, :baseline_interval, nothing)
    channel_selection_func = get(kwargs, :channel_selection, channels())
    channel_plot_order = get(kwargs, :channel_plot_order, nothing)

    filtered_kwargs = pairs(NamedTuple(filter(p -> p[1] ∉ (:baseline_interval, :channel_selection, :channel_plot_order), pairs(kwargs))))
    plot_kwargs, user_provided_color = _prepare_plot_kwargs(filtered_kwargs)

    dat_subset, all_channels, updated_channel_selection, _ = _prepare_erp_data(
        datasets,
        plot_kwargs;
        condition_selection = conditions(),
        channel_selection = channel_selection_func,
        sample_selection = get(plot_kwargs, :sample_selection, samples()),
        baseline_interval = baseline_interval,
    )

    if dynamic_grand_avg
        @info "Averaging across participants dynamically..."
        dat_subset, variances_dict = _compute_dynamic_grand_averages(dat_subset, plot_kwargs[:error_bars])
        sel_func = channel_selection_func
    else
        variances_dict = Dict{Int,Tuple{DataFrames.DataFrame,Int}}()
        sel_func = updated_channel_selection
    end

    selected_channels = get_selected_channels(first(dat_subset), sel_func; include_meta = false, include_extra = true)
    all_plot_channels = _order_plot_channels(selected_channels, all_channels, channel_plot_order)

    _plot_erp!(ax, dat_subset, all_plot_channels, plot_kwargs; user_provided_color = user_provided_color, variances_dict = variances_dict)
    return ax
end

"""
    plot_erp!(fig::Figure, ax::Axis, dat::ErpData; kwargs...)
    plot_erp!(fig::Figure, ax::Axis, datasets::AbstractVector{<:ErpData}; kwargs...)

Plot ERP data on an existing axis, mutating the figure and axis. Returns `ax`.
"""
plot_erp!(fig::Figure, ax::Axis, dat::ErpData; kwargs...) = plot_erp!(fig, ax, [dat]; kwargs...)
plot_erp!(fig::Figure, ax::Axis, datasets::AbstractVector{<:ErpData}; kwargs...) =
    _plot_erp_on_existing_axis!(ax, datasets; dynamic_grand_avg = false, kwargs...)

"""
    plot_erp_errorbar!(fig::Figure, ax::Axis, dat::ErpData; kwargs...)
    plot_erp_errorbar!(fig::Figure, ax::Axis, datasets::AbstractVector{<:ErpData}; kwargs...)

Plot grand-averaged ERP data with error ribbons on an existing axis, mutating the figure and axis. Returns `ax`.
"""
plot_erp_errorbar!(fig::Figure, ax::Axis, dat::ErpData; kwargs...) = plot_erp_errorbar!(fig, ax, [dat]; kwargs...)
plot_erp_errorbar!(fig::Figure, ax::Axis, datasets::AbstractVector{<:ErpData}; kwargs...) =
    _plot_erp_on_existing_axis!(ax, datasets; dynamic_grand_avg = true, kwargs...)

function _plot_erp!(
    ax::Axis,
    datasets::AbstractVector{<:ErpData},
    channels::AbstractVector{Symbol},
    plot_kwargs::Dict{Symbol,Any};
    condition_mask::AbstractVector{Bool} = Bool[],
    line_refs = nothing,
    user_provided_color::Bool = false,
    variances_dict::Dict{Int,Tuple{DataFrames.DataFrame,Int}} = Dict{Int,Tuple{DataFrames.DataFrame,Int}}(),
)
    def_theme = Makie.current_default_theme()
    color_cycle = haskey(def_theme, :palette) ? def_theme[:palette][:color][] : Makie.wong_colors()
    resolved_colormap = _resolve_theme_colormap(ax, plot_kwargs[:colormap], nothing)
    n_datasets = length(datasets)
    n_channels = length(channels)
    all_colors = _compute_dataset_colors(plot_kwargs[:color], n_datasets, n_channels, resolved_colormap, user_provided_color, color_cycle)
    all_linestyles = _compute_dataset_linestyles(plot_kwargs[:linestyle], n_datasets)
    lw = _resolve_theme_linewidth(ax, plot_kwargs[:linewidth], 2)
    legend_labels = plot_kwargs[:legend_labels]
    error_type = plot_kwargs[:error_bars]

    for (dataset_idx, dat) in enumerate(datasets)
        time_col = dat.data[!, :time]
        condition_visible = isempty(condition_mask) || dataset_idx > length(condition_mask) ? true : condition_mask[dataset_idx]
        base_label = isempty(legend_labels) ? dat.condition_name : legend_labels[dataset_idx]
        has_var = error_type != :none && haskey(variances_dict, dataset_idx)

        for (channel_idx, channel) in enumerate(channels)
            label = n_channels > 1 ? "$base_label ($channel)" : base_label
            color_idx = (dataset_idx - 1) * n_channels + channel_idx
            y_obs = Observable(dat.data[!, channel])

            if has_var
                var_df, n_stat = variances_dict[dataset_idx]
                if hasproperty(var_df, channel)
                    var_col = var_df[!, channel]
                    inv_sqrt_n = (error_type == :ci95 ? 1.96 : 1.0) / sqrt(n_stat)
                    err_val = [sqrt(v) * inv_sqrt_n for v in var_col]

                    lower_bound = lift(y -> y .- err_val, y_obs)
                    upper_bound = lift(y -> y .+ err_val, y_obs)
                    band_color = to_color(all_colors[color_idx])
                    band!(
                        ax,
                        time_col,
                        lower_bound,
                        upper_bound,
                        color = RGBAf(band_color.r, band_color.g, band_color.b, 0.2),
                        visible = condition_visible,
                    )
                end
            end

            line = lines!(
                ax,
                time_col,
                y_obs,
                linewidth = lw,
                color = all_colors[color_idx],
                linestyle = all_linestyles[dataset_idx],
                label = label,
                visible = condition_visible,
            )

            if !isnothing(line_refs)
                ds_refs = get!(line_refs, dataset_idx) do
                    Dict{Symbol,Tuple{Any,Observable}}()
                end
                ds_refs[channel] = (line, y_obs)
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
    return _add_legend!(ax, channels, datasets, plot_kwargs)
end

function _prepare_plot_kwargs(kwargs)
    kwargs_norm = _normalize_plot_kwargs(kwargs)
    user_provided_color = haskey(kwargs_norm, :color)
    return _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs_norm), user_provided_color
end

function _prepare_erp_data(
    datasets::AbstractVector{<:ErpData},
    plot_kwargs;
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
)
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channels(),
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = true,
    )

    if isempty(dat_subset)
        n_conditions = length(datasets)
        if n_conditions == 0
            @minimal_error "No data available (empty dataset)"
        else
            @minimal_error "No data matched the selection criteria. Available condition indices: 1:$n_conditions"
        end
    end

    if !isnothing(baseline_interval)
        @info "Applying baseline correction to $(length(dat_subset)) datasets"
        baseline!.(dat_subset, Ref(baseline_interval))
    end

    all_channels = vcat(channel_labels(dat_subset), extra_labels(dat_subset))
    selected_channels = get_selected_channels(first(dat_subset), channel_selection; include_meta = false, include_extra = true)
    all_plot_channels = [ch for ch in selected_channels if ch in all_channels]

    original_channels = nothing
    if plot_kwargs[:average_channels]
        @info "Averaging channels for $(length(dat_subset)) datasets"
        original_channels = all_plot_channels
        ch_sel = [channels(all_plot_channels)]
        for dat in dat_subset
            channel_average!(dat; channel_selections = ch_sel, output_labels = [:avg], reduce = false)
        end
        all_channels = [:avg]
        channel_selection = channels(:avg)
    end

    return dat_subset, all_channels, channel_selection, original_channels
end

@inline function _accum_col!(buf::Vector{Float64}, col::AbstractVector)
    @inbounds @simd for i in eachindex(buf)
        buf[i] += col[i]
    end
end

@inline function _accum_sqdiff!(var_buf::Vector{Float64}, col::AbstractVector, avg_buf::Vector{Float64})
    @inbounds @simd for i in eachindex(var_buf)
        diff = col[i] - avg_buf[i]
        var_buf[i] += diff * diff
    end
end

@inline function _accum_within_sqdiff!(
    var_buf::Vector{Float64},
    col::AbstractVector,
    p_mean::Vector{Float64},
    g_mean::Vector{Float64},
    avg_buf::Vector{Float64},
)
    @inbounds @simd for i in eachindex(var_buf)
        diff = (col[i] - p_mean[i] + g_mean[i]) - avg_buf[i]
        var_buf[i] += diff * diff
    end
end

function _compute_dynamic_grand_averages(datasets::Vector{ErpData}, error_bars::Symbol)
    if all(dat -> dat.file == "grand_avg", datasets)
        @minimal_warning "Datasets are already grand averages. Skipping dynamic participant averaging."
        return datasets, Dict{Int,Tuple{DataFrames.DataFrame,Int}}()
    end

    erps_by_condition = OrderedDict{Int,Vector{ErpData}}()
    for dat in datasets
        push!(get!(erps_by_condition, dat.condition, ErpData[]), dat)
    end

    grand_averages = ErpData[]
    variances_dict = Dict{Int,Tuple{DataFrames.DataFrame,Int}}()

    participant_global_means = Dict{String,DataFrames.DataFrame}()
    grand_global_mean = DataFrames.DataFrame()
    if error_bars == :within_sem
        participant_erps = Dict{String,Vector{ErpData}}()
        for erp in datasets
            push!(get!(participant_erps, erp.file, ErpData[]), erp)
        end

        first_erp = datasets[1]
        metadata_cols = meta_labels(first_erp)
        eeg_channels = collect(intersect((setdiff(propertynames(erp.data), metadata_cols) for erp in datasets)...))
        n_points = n_samples(first_erp)

        for (pid, erplist) in participant_erps
            df_mean = DataFrames.DataFrame()
            inv_n_p = 1.0 / length(erplist)
            for ch in eeg_channels
                avg_buf = zeros(Float64, n_points)
                for erp in erplist
                    _accum_col!(avg_buf, erp.data[!, ch])
                end
                avg_buf .*= inv_n_p
                df_mean[!, ch] = avg_buf
            end
            participant_global_means[pid] = df_mean
        end

        n_participants = length(participant_global_means)
        if n_participants > 0
            inv_np = 1.0 / n_participants
            for ch in eeg_channels
                avg_buf = zeros(Float64, n_points)
                for df_mean in values(participant_global_means)
                    _accum_col!(avg_buf, df_mean[!, ch]::Vector{Float64})
                end
                avg_buf .*= inv_np
                grand_global_mean[!, ch] = avg_buf
            end
        end
    end

    K_conditions = length(erps_by_condition)
    cm_correction = K_conditions > 1 ? (K_conditions / (K_conditions - 1)) : 1.0

    for (dataset_idx, cond_num) in enumerate(sort!(collect(keys(erps_by_condition))))
        erps = erps_by_condition[cond_num]
        first_erp = erps[1]
        metadata_cols = meta_labels(first_erp)
        eeg_channels = collect(intersect((setdiff(propertynames(erp.data), metadata_cols) for erp in erps)...))

        grand_avg_data = copy(first_erp.data)
        var_df = DataFrames.DataFrame()
        for col in metadata_cols
            hasproperty(grand_avg_data, col) && (var_df[!, col] = grand_avg_data[!, col])
        end

        n_points = nrow(grand_avg_data)
        n_erps = length(erps)
        inv_n = 1.0 / n_erps
        avg_buf = Vector{Float64}(undef, n_points)

        for ch in eeg_channels
            fill!(avg_buf, 0.0)
            for erp in erps
                _accum_col!(avg_buf, erp.data[!, ch])
            end
            avg_buf .*= inv_n
            grand_avg_data[!, ch] = copy(avg_buf)

            if n_erps > 1 && error_bars != :none
                var_buf = zeros(Float64, n_points)
                if error_bars == :within_sem
                    g_mean = grand_global_mean[!, ch]::Vector{Float64}
                    for erp in erps
                        _accum_within_sqdiff!(
                            var_buf,
                            erp.data[!, ch],
                            participant_global_means[erp.file][!, ch]::Vector{Float64},
                            g_mean,
                            avg_buf,
                        )
                    end
                    var_buf .*= (cm_correction / (n_erps - 1))
                else
                    for erp in erps
                        _accum_sqdiff!(var_buf, erp.data[!, ch], avg_buf)
                    end
                    var_buf ./= (n_erps - 1)
                end
                var_df[!, ch] = var_buf
            end
        end

        push!(
            grand_averages,
            ErpData(
                "grand_avg",
                cond_num,
                "grand_avg_$(first_erp.condition_name)",
                grand_avg_data,
                copy(first_erp.layout),
                first_erp.sample_rate,
                copy(first_erp.analysis_info),
                sum(erp.n_epochs for erp in erps),
            ),
        )
        error_bars != :none && n_erps > 1 && (variances_dict[dataset_idx] = (var_df, n_erps))
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
    menu_fig = _eegfun_figure(size = (300, 300))

    data_to_plot = _filter_visible_conditions(data, condition_checked_ref)
    plot_types = data_to_plot isa Vector{ErpData} && length(data_to_plot) > 1 ? ["Topoplot", "Topoplot (average)"] : ["Topoplot"]

    for (idx, plot_type) in enumerate(plot_types)
        btn = Button(menu_fig[idx, 1], label = plot_type)
        on(btn.clicks) do _
            original_data, x_min, x_max = _get_erp_selection_bounds(selection_state, data)
            time_sel = x -> (x.time .>= x_min) .& (x.time .<= x_max)
            vis_data = _filter_visible_conditions(original_data, condition_checked_ref)

            if btn.label[] == "Topoplot"
                plot_topography(vis_data, sample_selection = time_sel)
            elseif btn.label[] == "Topoplot (average)"
                plot_topography(_average_conditions(vis_data), sample_selection = time_sel)
            end
        end
    end

    _display_in_screen(menu_fig; size = (300, 300))
end

function _get_erp_selection_bounds(selection_state, data)
    x_min, x_max = minmax(selection_state.bounds[]...)
    return (data, x_min, x_max)
end

function _average_conditions(erps::Vector{ErpData})
    isempty(erps) && @minimal_error("Cannot average empty ERP list")
    length(erps) == 1 && return erps[1]

    avg_erp = _create_grand_average(erps, first(erps).condition)
    avg_cond_name = "avg_" * join((erp.condition_name for erp in erps), "_")
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

_average_conditions(erp::ErpData) = erp

function _show_channel_average_menu!(selected_channels, data, condition_checked_ref)
    _set_window_title("Channel Selection Menu")
    menu_fig = _eegfun_figure(size = (300, 150))
    data_to_plot = _filter_visible_conditions(data, condition_checked_ref)

    btn = Button(menu_fig[1, 1], label = "Plot Average of Selected ($(length(selected_channels)) ch)")
    on(btn.clicks) do _
        plot_erp(data_to_plot; channel_selection = channels(selected_channels), average_channels = true)
    end
    _display_in_screen(menu_fig; size = (300, 150))
end

function _draw_highlight_regions!(ax::Axis, regions)
    region_list = regions isa AbstractVector ? regions : (regions,)
    for region in region_list
        x1, x2 = region.x1, region.x2
        color = get(region, :color, :gray)
        alpha = get(region, :alpha, 0.3)
        if haskey(region, :y1) || haskey(region, :y2)
            y1 = get(region, :y1, -1e10)
            y2 = get(region, :y2, 1e10)
            poly!(ax, Point2f[(x1, y1), (x2, y1), (x2, y2), (x1, y2)], color = (color, alpha), strokewidth = 0)
        else
            vspan!(ax, x1, x2, color = (color, alpha))
        end
    end
end

function _compute_dataset_colors(color_val, n_datasets::Int, n_channels::Int, colormap, user_provided_color::Bool, color_cycle)
    n_total = n_datasets * n_channels
    if color_val isa Vector
        return [color_val[mod1(i, length(color_val))] for i = 1:n_total]
    elseif user_provided_color && !isnothing(color_val)
        return fill(color_val, n_total)
    elseif n_total > 1 && !isnothing(colormap) && colormap != Makie.automatic
        gradient = Makie.cgrad(colormap, n_total, categorical = true)
        return [gradient[i] for i = 1:n_total]
    else
        return [color_cycle[mod1(i, length(color_cycle))] for i = 1:n_total]
    end
end

function _compute_dataset_linestyles(linestyle_val, n_datasets::Int)
    if linestyle_val isa Vector
        return [linestyle_val[mod1(i, length(linestyle_val))] for i = 1:n_datasets]
    elseif isnothing(linestyle_val) || linestyle_val == Makie.automatic
        return fill(nothing, n_datasets)
    else
        return fill(linestyle_val, n_datasets)
    end
end

function _setup_linked_legend_interactions!(line_refs::Vector{<:Dict})
    length(line_refs) > 1 || return

    dataset_lines = Dict{Any,Vector{Any}}()
    for ax_line_refs in line_refs
        for (dataset_idx, channel_lines) in ax_line_refs
            if length(channel_lines) == 1
                push!(get!(dataset_lines, dataset_idx, Any[]), first(values(channel_lines))[1])
            elseif length(channel_lines) > 1
                for (ch, line_data) in channel_lines
                    push!(get!(dataset_lines, (dataset_idx, ch), Any[]), line_data[1])
                end
            end
        end
    end

    foreach(_sync_lines!, values(dataset_lines))
end

function _setup_erp_control_panel!(
    fig::Figure,
    dat_subset::Vector{ErpData},
    axes::Vector{Axis},
    baseline_interval::Interval,
    line_refs::Union{Vector{<:Dict},Nothing} = nothing,
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}} = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing),
    legend_refs = nothing,
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)
    !isnothing(line_refs) && _setup_linked_legend_interactions!(line_refs)

    start_val, stop_val = _extract_baseline_values(baseline_interval)
    baseline_start_obs = Observable(isnothing(start_val) ? "" : string(start_val))
    baseline_stop_obs = Observable(isnothing(stop_val) ? "" : string(stop_val))
    condition_checked = [Observable(true) for _ in dat_subset]
    condition_checked_ref[] = condition_checked
    previous_baseline = Ref{Union{Tuple{Float64,Float64},Nothing}}(nothing)

    function update_plot!()
        s_val, e_val = _parse_baseline_values(baseline_start_obs[], baseline_stop_obs[])
        baseline_interval_new = (!isnothing(s_val) && !isnothing(e_val)) ? (s_val, e_val) : nothing
        baseline_changed = baseline_interval_new != previous_baseline[]

        if baseline_changed && !isnothing(baseline_interval_new)
            baseline!.(dat_subset, Ref(baseline_interval_new))
            previous_baseline[] = baseline_interval_new
        end

        condition_mask = [checked[] for checked in condition_checked]
        isnothing(line_refs) && return

        for (ax_idx, ax_line_refs) in enumerate(line_refs)
            ax_idx > length(axes) && continue
            for (dataset_idx, channel_lines) in ax_line_refs
                dataset_idx > length(dat_subset) && continue
                dat = dat_subset[dataset_idx]
                new_vis = dataset_idx <= length(condition_mask) ? condition_mask[dataset_idx] : true

                for (channel, (line, y_obs)) in channel_lines
                    baseline_changed && (y_obs[] = dat.data[!, channel])
                    line.visible[] != new_vis && (line.visible = new_vis)
                end
            end
        end
    end

    for checked in condition_checked
        on(_ -> update_plot!(), checked)
    end

    control_open = Ref(false)

    on(events(fig).keyboardbutton) do event
        event.action == Keyboard.press || return
        if event.key == Keyboard.c
            if control_open[] && !isnothing(control_fig[])
                _display_in_screen(control_fig[])
                return
            end

            control_fig[] = _eegfun_figure(title = "ERP Control Panel", size = (320, 420))
            control_open[] = true
            layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)

            Label(layout[1, 1], "Baseline Interval", fontsize = 14, font = :bold)
            baseline_layout = GridLayout(layout[2, 1], tellwidth = false, colgap = 10)
            _create_baseline_textbox(baseline_layout, 1, "Start (s):", baseline_start_obs, " ", 100)
            _create_baseline_textbox(baseline_layout, 2, "End (s):", baseline_stop_obs, " ", 100)

            apply_btn = Button(layout[3, 1], label = "Apply Baseline", width = 200)
            on(_ -> update_plot!(), apply_btn.clicks)

            Label(layout[4, 1], "Conditions", fontsize = 14, font = :bold)
            conditions_layout = GridLayout(layout[5, 1], tellwidth = false, rowgap = 5)
            for (idx, dat) in enumerate(dat_subset)
                cb = Checkbox(conditions_layout[idx, 1], checked = condition_checked[idx][])
                Label(conditions_layout[idx, 2], dat.condition_name)
                on(val -> (condition_checked[idx][] = val), cb.checked)
            end

            leg_btn = Button(layout[6, 1], label = "Toggle Legend (L)", width = 200)
            on(_ -> _toggle_legends!(legend_refs), leg_btn.clicks)

            _display_in_screen(control_fig[])
        elseif event.key == Keyboard.l
            _toggle_legends!(legend_refs)
        end
    end
end
