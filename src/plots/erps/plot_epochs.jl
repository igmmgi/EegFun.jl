const PLOT_EPOCHS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Data processing
    :average_channels => (false, "Whether to average across channels"),
    :plot_avg_trials => (true, "Whether to draw ERP average overlay"),
    :baseline_interval => (nothing, "Baseline interval as (start, stop) tuple or range in seconds"),

    # Axis limits, labels, and direction
    :plot_title => (nothing, "Plot title. If nothing, automatically determined"),
    :show_title => (true, "Whether to show the plot title"),
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

function _require_batch_files(pattern::String, input_dir::String, participant_selection; kind::String = "")
    files = _find_batch_files(pattern, input_dir, participant_selection)
    if isempty(files)
        all_matching = _find_batch_files(pattern, input_dir)
        label = isempty(kind) ? "files" : "$kind files"
        if isempty(all_matching)
            @minimal_error "No $label matching pattern '$pattern' in $input_dir"
        else
            avail_ids = sort(unique(_extract_participant_id.(all_matching)))
            @minimal_error "Pattern '$pattern' matched $(length(all_matching)) file(s), but none passed the participant selection. Available participant IDs: $avail_ids"
        end
    end
    return sort(files, by = _natural_sort_key)
end

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
    end

    files = _require_batch_files(filename, input_dir, participant_selection)
    results = NamedTuple[]
    for file in files
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(
            results,
            plot_epochs(
                data;
                channel_selection = channel_selection,
                sample_selection = sample_selection,
                interval_selection = interval_selection,
                epoch_selection = epoch_selection,
                include_extra = include_extra,
                layout = layout,
                kwargs...,
            ),
        )
    end
    return results
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
)
    user_provided_color = haskey(kwargs, :color)
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)

    # Subset by condition, sample, and epoch (keep all channels available for interactive topoplots)
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channels(),
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )
    isempty(dat_subset) && @minimal_error "No data matched the selection criteria!"

    !isnothing(plot_kwargs[:baseline_interval]) && baseline!(dat_subset, plot_kwargs[:baseline_interval])

    _set_window_title(_generate_window_title(dat_subset))

    try
        all_channels = vcat(channel_labels(dat_subset[1]), extra_labels(dat_subset[1]))
        selected_channels = get_selected_channels(first(dat_subset), channel_selection; include_meta = false, include_extra = include_extra)
        all_plot_channels = [ch for ch in selected_channels if ch in all_channels]
        isempty(all_plot_channels) && @minimal_error "No channels selected for plotting"

        if plot_kwargs[:average_channels]
            ch_sel = [channels(all_plot_channels)]
            for dat in dat_subset
                channel_average!(dat; channel_selections = ch_sel, output_labels = [:avg], reduce = false)
            end
        end

        dat_subset_avg = [average_epochs(dat) for dat in dat_subset]
        time_cat_list = [_build_time_cat(n_epochs(dat), dat.data[1][!, :time]) for dat in dat_subset]

        n_conditions = length(dat_subset)
        n_channels_for_colors =
            (layout == :single && length(all_plot_channels) > 1 && !plot_kwargs[:average_channels]) ? length(all_plot_channels) : 1
        def_theme = Makie.current_default_theme()
        color_cycle = haskey(def_theme, :palette) ? def_theme[:palette][:color][] : Makie.wong_colors()
        resolved_colormap =
            isnothing(plot_kwargs[:colormap]) ? (haskey(def_theme, :colormap) ? def_theme[:colormap][] : nothing) : plot_kwargs[:colormap]
        condition_colors_list = _compute_dataset_colors(
            plot_kwargs[:color],
            n_conditions,
            n_channels_for_colors,
            resolved_colormap,
            user_provided_color,
            color_cycle,
        )

        @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(n_conditions) conditions"

        layout_kwargs = pop!(plot_kwargs, :layout_kwargs, (;))
        plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout; layout_kwargs...)

        fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
        fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
        axes = Axis[]
        legend_refs = Union{Legend,Nothing}[]

        line_refs = if plot_kwargs[:interactive]
            n_axes = (layout == :grid || layout == :topo) ? length(all_plot_channels) : 1
            [Dict{Int,Dict{Symbol,Any}}() for _ = 1:n_axes]
        else
            nothing
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

            leg = _add_legend!(ax, [:avg], dat_subset, plot_kwargs)
            !isnothing(leg) && push!(legend_refs, leg)

            ax.title =
                length(all_plot_channels) == 1 ? string(all_plot_channels[1]) :
                "Avg: $(_print_vector(all_plot_channels, max_length = 8, n_ends = 3))"
            _finalize_epochs_axis!(ax, plot_kwargs)
        elseif layout == :single
            ax = Axis(fig[1, 1])
            push!(axes, ax)
            ax_line_refs = plot_kwargs[:interactive] ? line_refs[1] : nothing
            n_ch = length(all_plot_channels)

            for (ch_idx, ch) in enumerate(all_plot_channels)
                ch_colors = [condition_colors_list[n_ch > 1 ? ((ci - 1) * n_ch + ch_idx) : ci] for ci = 1:n_conditions]
                channel_suffix = n_ch > 1 ? " ($ch)" : ""
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

            leg = _add_legend!(ax, all_plot_channels, dat_subset, plot_kwargs)
            !isnothing(leg) && push!(legend_refs, leg)
            _finalize_epochs_axis!(ax, plot_kwargs)
        elseif layout == :grid
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

        if plot_kwargs[:interactive]
            _setup_shared_interactivity!(fig, axes, :epochs; zoom_step = plot_kwargs[:zoom_step])
            for ax in axes
                deregister_interaction!(ax, :rectanglezoom)
            end

            selection_state =
                SharedSelectionState(axes; selection_color = plot_kwargs[:selection_color], selection_alpha = plot_kwargs[:selection_alpha])

            condition_checked_ref = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing)
            _setup_epochs_control_panel!(
                fig,
                dat_subset,
                dat_subset_avg,
                axes,
                get(plot_kwargs, :baseline_interval, nothing),
                line_refs,
                condition_checked_ref,
                legend_refs,
            )

            right_click_handler = (sel, mouse_x, data) -> _handle_epochs_right_click!(sel, mouse_x, data, condition_checked_ref)
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, right_click_handler)

            if plot_layout.type in (:topo, :grid)
                _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes)
            end
        end

        plot_kwargs[:display_plot] && _display_figure(fig)
        return (fig = fig, axes = axes)
    finally
        _set_window_title("Makie")
    end
end

function plot_epochs(
    dat::EpochData;
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)
    return plot_epochs(
        [dat];
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
        layout = layout,
        kwargs...,
    )
end

function _finalize_epochs_axis!(ax::Axis, kw::Dict; xlim = kw[:xlim], ylim = kw[:ylim], xlabel = kw[:xlabel], ylabel = kw[:ylabel])
    if !get(kw, :show_title, true)
        ax.title = ""
    elseif !isnothing(kw[:plot_title]) && isempty(ax.title[])
        ax.title = string(kw[:plot_title])
    end
    _set_axis_properties!(
        ax;
        xlim = xlim,
        ylim = ylim,
        xlabel = xlabel,
        ylabel = ylabel,
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
        xlabel = xlabel,
        ylabel = ylabel,
    )
end

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
    actual_lw = _resolve_theme_linewidth(ax, plot_kwargs[:linewidth], 1)
    base_lw = actual_lw isa Vector ? actual_lw[1] : actual_lw
    avg_lw = base_lw * plot_kwargs[:avg_linewidth_multiplier]
    trial_alpha = plot_kwargs[:trial_alpha]
    plot_avg = plot_kwargs[:plot_avg_trials]
    legend_labels = plot_kwargs[:legend_labels]

    for (cond_idx, dat) in enumerate(datasets)
        label = _get_condition_label(dat, cond_idx, legend_labels, channel_suffix)
        cond_color = colors[cond_idx]
        line_color = cond_color isa Vector ? cond_color[1] : cond_color
        tc = isnothing(time_cat_list) ? nothing : time_cat_list[cond_idx]

        trial_line, trial_y_obs =
            _plot_epochs!(ax, dat, channel; color = line_color, linewidth = base_lw, alpha = trial_alpha, time_cat = tc, label = label)
        !isnothing(ax_line_refs) && _store_line_ref!(ax_line_refs, cond_idx, channel, trial_line, trial_y_obs)

        if plot_avg
            erp_dat = !isnothing(dat_subset_avg) ? dat_subset_avg[cond_idx] : average_epochs(dat)
            avg_line, y_obs = _plot_erp_average!(ax, erp_dat, channel; color = line_color, linewidth = avg_lw, label = label * " (avg)")
            !isnothing(ax_line_refs) && _store_avg_line_ref!(ax_line_refs, cond_idx, channel, avg_line, y_obs)
        end
    end
end

function _plot_epochs!(
    ax,
    dat,
    ch::Symbol;
    color = nothing,
    linewidth = 1,
    alpha = 0.25,
    time_cat::Union{AbstractVector,Nothing} = nothing,
    label::Union{String,Nothing} = nothing,
)
    time_vec = dat.data[1][!, :time]
    tc = isnothing(time_cat) ? _build_time_cat(n_epochs(dat), time_vec) : time_cat
    y_obs = Observable(_concatenate_trials(dat.data, ch, time_vec))
    line = lines!(ax, tc, y_obs, color = color, linewidth = linewidth, alpha = alpha, label = label)
    return line, y_obs
end

function _plot_erp_average!(ax, erp_dat::ErpData, ch::Symbol; color = nothing, linewidth = 2, label::Union{String,Nothing} = nothing)
    y_obs = Observable(erp_dat.data[!, ch])
    line = lines!(ax, erp_dat.data[!, :time], y_obs, color = color, linewidth = linewidth, alpha = 1.0, label = label)
    return line, y_obs
end

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
    time_cat_list::Union{Vector{<:AbstractVector},Nothing} = nothing,
)
    rows, cols = grid_dims
    ylim = something_or_calc_ylim(plot_kwargs[:ylim], datasets, all_plot_channels)

    for (idx, channel) in enumerate(all_plot_channels)
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)

        ax_line_refs = !isnothing(line_refs) && idx <= length(line_refs) ? line_refs[idx] : nothing
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
        _add_legend!(ax, [channel], datasets, plot_kwargs)

        col != 1 && (ax.yticklabelsvisible = false)
        row != rows && (ax.xticklabelsvisible = false)
        ax.title = string(channel)
        _finalize_epochs_axis!(
            ax,
            plot_kwargs;
            ylim = ylim,
            xlabel = row == rows ? plot_kwargs[:xlabel] : "",
            ylabel = col == 1 ? plot_kwargs[:ylabel] : "",
        )
    end

    length(axes) > 1 && Makie.linkaxes!(axes...)
end

something_or_calc_ylim(ylim, datasets, chs) = isnothing(ylim) ? _ylimits(datasets; channel_selection = channels(chs)) : ylim

function _plot_epochs_topo_multi!(
    fig::Figure,
    axes::Vector{Axis},
    datasets::Vector{EpochData},
    all_plot_channels::Vector{Symbol},
    condition_colors::Vector,
    plot_kwargs::Dict,
    line_refs = nothing;
    dat_subset_avg::Union{Vector{ErpData},Nothing} = nothing,
    time_cat_list::Union{Vector{<:AbstractVector},Nothing} = nothing,
    plot_layout = nothing,
)
    dat = datasets[1]
    if !all(in.([:x2, :y2], Ref(propertynames(dat.layout.data))))
        polar_to_cartesian_xy!(dat.layout)
    end

    ylim = something_or_calc_ylim(plot_kwargs[:ylim], datasets, all_plot_channels)
    plot_w = get(plot_kwargs, :layout_topo_plot_width, !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_plot_width, 0.1) : 0.1)
    plot_h = get(plot_kwargs, :layout_topo_plot_height, !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_plot_height, 0.1) : 0.1)

    x2 = dat.layout.data[!, :x2]
    y2 = dat.layout.data[!, :y2]
    labels = dat.layout.data[!, :label]
    minx, maxx = extrema(x2)
    miny, maxy = extrema(y2)
    xrange = maxx == minx ? 1.0 : maxx - minx
    yrange = maxy == miny ? 1.0 : maxy - miny

    for (ch_idx, ch) in enumerate(all_plot_channels)
        idx = findfirst(==(ch), labels)
        halign, valign = isnothing(idx) ? (0.5, 0.5) : ((x2[idx] - minx) / xrange, (y2[idx] - miny) / yrange)

        ax = Axis(fig[1, 1], width = Relative(plot_w), height = Relative(plot_h), halign = halign, valign = valign)
        push!(axes, ax)

        ax_line_refs = !isnothing(line_refs) && ch_idx <= length(line_refs) ? line_refs[ch_idx] : nothing
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
        _add_legend!(ax, [ch], datasets, plot_kwargs)

        ax.title = string(ch)
        _finalize_epochs_axis!(ax, plot_kwargs; ylim = ylim, xlabel = "", ylabel = "")
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticksvisible = false
        hidespines!(ax)
    end

    show_scale =
        get(plot_kwargs, :layout_topo_show_scale, !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_show_scale, false) : false)
    if show_scale
        scale_pos = get(
            plot_kwargs,
            :layout_topo_scale_pos,
            !isnothing(plot_layout) ? get(plot_layout.metadata, :topo_scale_pos, (0.8, -0.8)) : (0.8, -0.8),
        )
        scale_halign = 0.0 <= scale_pos[1] <= 1.0 ? scale_pos[1] : clamp((scale_pos[1] + 1) / 2, 0.0, 1.0)
        scale_valign = 0.0 <= scale_pos[2] <= 1.0 ? scale_pos[2] : clamp((scale_pos[2] + 1) / 2, 0.0, 1.0)
        scale_ax = Axis(fig[1, 1], width = Relative(plot_w), height = Relative(plot_h), halign = scale_halign, valign = scale_valign)
        push!(axes, scale_ax)
        t_col = dat.data[1][!, :time]
        scale_ax.title = ""
        _finalize_epochs_axis!(scale_ax, plot_kwargs; ylim = ylim, xlim = (t_col[1], t_col[end]))
        scale_ax.xticklabelsvisible = true
        scale_ax.yticklabelsvisible = true
        scale_ax.xticksvisible = true
        scale_ax.yticksvisible = true
    end

    length(axes) > 1 && Makie.linkaxes!(axes...)
end

function _add_legend!(ax::Axis, channels::AbstractVector{Symbol}, datasets::AbstractVector, kwargs::Dict)
    kwargs[:legend] || return nothing
    if !isempty(kwargs[:legend_channel]) && isempty(intersect(kwargs[:legend_channel], channels))
        return nothing
    end

    legend_label = kwargs[:legend_label]
    legend_position = kwargs[:legend_position]
    nbanks = something(get(kwargs, :legend_nbanks, nothing), length(channels) > 10 ? cld(length(channels), 10) : 1)
    legend_kwargs = get(kwargs, :legend_kwargs, (;))
    framevisible = get(kwargs, :legend_framevisible, true)

    leg = if isempty(legend_label)
        axislegend(ax; position = legend_position, nbanks = nbanks, framevisible = framevisible, legend_kwargs...)
    else
        axislegend(ax, legend_label; position = legend_position, nbanks = nbanks, framevisible = framevisible, legend_kwargs...)
    end

    if length(channels) == 1 && length(datasets) == 1
        leg.blockscene.visible[] = false
    end
    return leg
end

function _get_condition_label(dat, cond_idx, legend_labels, channel_suffix = "")
    label = isempty(legend_labels) ? dat.condition_name : (cond_idx <= length(legend_labels) ? legend_labels[cond_idx] : dat.condition_name)
    return isempty(channel_suffix) ? label : label * channel_suffix
end

function _store_line_ref!(ax_line_refs, cond_idx, ch, trial_line, trial_y_obs)
    cond_dict = get!(ax_line_refs, cond_idx) do
        Dict{Symbol,Any}()
    end
    cond_dict[ch] = Dict{Symbol,Any}(:trials => (trial_line, trial_y_obs), :average => nothing)
end

_store_avg_line_ref!(ax_line_refs, cond_idx, ch, avg_line, y_obs) = (ax_line_refs[cond_idx][ch][:average] = (avg_line, y_obs))

function _build_time_cat(n_trials::Int, time_vec::AbstractVector{T}) where {T}
    n = length(time_vec)
    time_cat = Vector{T}(undef, n_trials * n + (n_trials - 1))
    pos = 1
    @inbounds for t = 1:n_trials
        copyto!(time_cat, pos, time_vec, 1, n)
        pos += n
        if t != n_trials
            time_cat[pos] = T(NaN)
            pos += 1
        end
    end
    return time_cat
end

function _concatenate_trials_col(::AbstractVector{T}, trials::Vector{DataFrame}, col_idx::Int, n::Int) where {T}
    m = length(trials)
    y_cat = Vector{T}(undef, m * n + (m - 1))
    pos = 1
    @inbounds for t = 1:m
        col = trials[t][!, col_idx]
        if col isa Vector{T}
            copyto!(y_cat, pos, col, 1, n)
        else
            copyto!(y_cat, pos, col::AbstractVector{T}, 1, n)
        end
        pos += n
        if t != m
            y_cat[pos] = T(NaN)
            pos += 1
        end
    end
    return y_cat
end

function _concatenate_trials(trials, ch, time_vec)
    col_idx = ch isa Int ? ch : DataFrames.columnindex(trials[1], ch)
    return _concatenate_trials_col(trials[1][!, col_idx], trials, col_idx, length(time_vec))
end

function _setup_linked_legend_interactions_epochs!(line_refs::Vector{<:Dict})
    length(line_refs) > 1 || return

    condition_trial_lines = Dict{Any,Vector{Any}}()
    condition_avg_lines = Dict{Any,Vector{Any}}()

    for ax_line_refs in line_refs
        for (cond_idx, cond_line_data) in ax_line_refs
            for (ch, ch_line_data) in cond_line_data
                ch_line_data isa Dict || continue
                group_key = length(cond_line_data) == 1 ? cond_idx : (cond_idx, ch)

                for (slot, target_dict) in ((:trials, condition_trial_lines), (:average, condition_avg_lines))
                    entry = get(ch_line_data, slot, nothing)
                    isnothing(entry) && continue
                    line_obj = entry isa Tuple ? entry[1] : (entry isa Lines ? entry : nothing)
                    !isnothing(line_obj) && push!(get!(target_dict, group_key, Any[]), line_obj)
                end
            end
        end
    end

    foreach(_sync_lines!, values(condition_trial_lines))
    foreach(_sync_lines!, values(condition_avg_lines))
end

function _sync_lines!(lines::Vector)
    length(lines) > 1 || return
    syncing = Ref(false)
    for line in lines
        other_lines = [l for l in lines if l !== line]
        on(line.visible) do visible_val
            syncing[] && return
            syncing[] = true
            for other_line in other_lines
                other_line.visible[] != visible_val && (other_line.visible = visible_val)
            end
            syncing[] = false
        end
    end
end

function _toggle_legends!(legend_refs)
    isnothing(legend_refs) && return
    for leg in legend_refs
        if !isnothing(leg) && isdefined(leg, :blockscene)
            leg.blockscene.visible[] = !leg.blockscene.visible[]
        end
    end
end

function _setup_epochs_control_panel!(
    fig::Figure,
    dat_subset::Vector{EpochData},
    dat_subset_avg::Vector{ErpData},
    axes::Vector{Axis},
    baseline_interval::Interval,
    line_refs::Vector{<:Dict},
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}} = Ref{Union{Vector{Observable{Bool}},Nothing}}(nothing),
    legend_refs::Vector = Union{Legend,Nothing}[],
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)
    control_screen = Ref{Any}(nothing)

    _setup_linked_legend_interactions_epochs!(line_refs)

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
            baseline!.(dat_subset_avg, Ref(baseline_interval_new))
            previous_baseline[] = baseline_interval_new
        end

        condition_mask = [checked[] for checked in condition_checked]

        for (ax_idx, ax_line_refs) in enumerate(line_refs)
            ax_idx > length(axes) && continue
            for (cond_idx, cond_line_data) in ax_line_refs
                cond_idx > length(dat_subset) && continue
                visible = cond_idx <= length(condition_mask) ? condition_mask[cond_idx] : true
                dat = dat_subset[cond_idx]

                for (ch, line_data) in cond_line_data
                    line_data isa Dict || continue

                    trials_entry = get(line_data, :trials, nothing)
                    if !isnothing(trials_entry)
                        trial_line, trial_y_obs = trials_entry
                        trial_line.visible[] != visible && (trial_line.visible = visible)
                        if baseline_changed
                            trial_y_obs[] = _concatenate_trials(dat.data, ch, dat.data[1][!, :time])
                        end
                    end

                    avg_entry = get(line_data, :average, nothing)
                    if !isnothing(avg_entry)
                        avg_line, y_obs = avg_entry
                        avg_line.visible[] != visible && (avg_line.visible = visible)
                        if baseline_changed && cond_idx <= length(dat_subset_avg)
                            y_obs[] = dat_subset_avg[cond_idx].data[!, ch]
                        end
                    end
                end
            end
        end
    end

    for checked in condition_checked
        on(_ -> update_plot!(), checked)
    end

    on(events(fig).keyboardbutton) do event
        event.action == Keyboard.press || return
        if event.key == Keyboard.c
            if !isnothing(control_screen[]) && hasmethod(isopen, Tuple{typeof(control_screen[])}) && isopen(control_screen[])
                return
            end

            control_fig[] = _eegfun_figure(title = "Epochs Control Panel", size = (300, 400))
            panel_layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)

            Label(panel_layout[1, 1], "Baseline Interval", fontsize = 14, font = :bold)
            baseline_layout = GridLayout(panel_layout[2, 1], tellwidth = false, colgap = 10)
            _create_baseline_textbox(baseline_layout, 1, "Start (s):", baseline_start_obs, " ", 100)
            _create_baseline_textbox(baseline_layout, 2, "End (s):", baseline_stop_obs, " ", 100)

            apply_btn = Button(panel_layout[3, 1], label = "Apply Baseline", width = 200)
            on(_ -> update_plot!(), apply_btn.clicks)

            Label(panel_layout[4, 1], "Conditions", fontsize = 14, font = :bold)
            conditions_layout = GridLayout(panel_layout[5, 1], tellwidth = false, rowgap = 5)
            for (i, dat) in enumerate(dat_subset)
                cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][])
                Label(conditions_layout[i, 2], dat.condition_name)
                connect!(condition_checked[i], cb.checked)
            end

            leg_btn = Button(panel_layout[6, 1], label = "Toggle Legend (L)", width = 200)
            on(_ -> _toggle_legends!(legend_refs), leg_btn.clicks)

            control_screen[] = _display_in_screen(control_fig[])
        elseif event.key == Keyboard.l
            _toggle_legends!(legend_refs)
        end
    end
end

function _handle_epochs_right_click!(selection_state, mouse_x, data, condition_checked_ref)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_epochs_context_menu!(selection_state, data, condition_checked_ref)
    end
end

function _show_epochs_context_menu!(selection_state, data, condition_checked_ref)
    menu_fig = _eegfun_figure()
    data_to_plot = _filter_visible_conditions(data, condition_checked_ref)
    plot_types = data_to_plot isa Vector{EpochData} && length(data_to_plot) > 1 ? ["Topoplot", "Topoplot (average)"] : ["Topoplot"]

    for (idx, plot_type) in enumerate(plot_types)
        btn = Button(menu_fig[idx, 1], label = plot_type)
        on(btn.clicks) do _
            original_data, x_min, x_max = _get_epochs_selection_bounds(selection_state, data)
            time_sel = x -> (x.time .>= x_min) .& (x.time .<= x_max)
            vis_data = _filter_visible_conditions(original_data, condition_checked_ref)
            erps = vis_data isa Vector{EpochData} ? [average_epochs(dat) for dat in vis_data] : [average_epochs(vis_data)]
            erp_data = length(erps) == 1 ? erps[1] : erps

            if btn.label[] == "Topoplot"
                plot_topography(erp_data, sample_selection = time_sel)
            elseif btn.label[] == "Topoplot (average)"
                plot_topography(erp_data isa Vector{ErpData} ? _average_conditions(erp_data) : erp_data, sample_selection = time_sel)
            end
        end
    end

    _display_in_screen(menu_fig)
end

function _get_epochs_selection_bounds(selection_state, data)
    x_min, x_max = minmax(selection_state.bounds[]...)
    return (data, x_min, x_max)
end

function _filter_visible_conditions(data, condition_checked_ref)
    if isnothing(condition_checked_ref[]) || !(data isa AbstractVector)
        return data
    end
    condition_checked = condition_checked_ref[]
    length(condition_checked) == length(data) || return data
    visible_data = [data[i] for i in eachindex(data) if condition_checked[i][]]
    return length(visible_data) == 1 ? visible_data[1] : visible_data
end
