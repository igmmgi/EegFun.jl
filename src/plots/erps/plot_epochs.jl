# Multi-trial and butterfly ERP/epoch waveform plotting functions.

######################################
######### Types & Constants ##########
######################################

const PLOT_EPOCHS_KWARGS = merge(
    _BASE_WAVEFORM_KWARGS,
    Dict{Symbol,Tuple{Any,String}}(
        :plot_avg_trials => (true, "Whether to draw ERP average overlay"),
        :baseline_interval => (nothing, "Baseline interval as (start, stop) tuple or range in seconds"),
        :avg_linewidth_multiplier => (2.0, "Multiplier for average line width."),
        :trial_alpha => (0.25, "Alpha (transparency) for individual trial traces"),
        :layout => (:single, "Layout type: :single, :grid, or :topo"),
    ),
)

######################################
########### User Functions ###########
######################################

"""
    plot_epochs(datasets::AbstractVector{<:EpochData}; condition_selection = conditions(), channel_selection = channels(), sample_selection = samples(), interval_selection = times(), epoch_selection = epochs(), include_extra = false, layout = :single, kwargs...)
    plot_epochs(dat::EpochData; kwargs...)
    plot_epochs(filename::String; input_dir = pwd(), participant_selection = participants(), kwargs...)

Plot single-trial epoch waveforms across `:single`, `:grid`, `:topo`, or custom `PlotLayout` layouts
from in-memory `EpochData` or a `.jld2` filepath / filename pattern.

$(_generate_kwargs_doc(PLOT_EPOCHS_KWARGS))

# Examples
```julia
plot_epochs(epochs_data; layout = :grid)
plot_epochs("epochs_raw.jld2")
plot_epochs("epochs_raw", participant_selection = participants(1))
```
"""
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
            n_axes = plot_layout.type in (:grid, :topo, :custom) ? length(plot_layout.channels) : 1
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
        elseif plot_layout.type == :single
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
        elseif plot_layout.type in (:grid, :topo, :custom)
            created_axes, layout_channels = _apply_layout!(fig, plot_layout; plot_kwargs...)
            append!(axes, created_axes)
            ylim = _something_or_calc_ylim(plot_kwargs[:ylim], dat_subset, layout_channels)
            for (idx, (ax, ch)) in enumerate(zip(axes, layout_channels))
                ax_line_refs = !isnothing(line_refs) && idx <= length(line_refs) ? line_refs[idx] : nothing
                _plot_epochs_conditions_on_axis!(
                    ax,
                    dat_subset,
                    ch,
                    condition_colors_list,
                    plot_kwargs,
                    ax_line_refs;
                    dat_subset_avg = dat_subset_avg,
                    time_cat_list = time_cat_list,
                )
                leg = _add_legend!(ax, [ch], dat_subset, plot_kwargs)
                !isnothing(leg) && push!(legend_refs, leg)
                _finalize_epochs_axis!(ax, plot_kwargs; ylim = ylim)
            end
            _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...)
            if plot_layout.type == :topo && get(plot_layout.metadata, :topo_show_scale, false)
                plot_w = get(plot_layout.metadata, :topo_plot_width, 0.1)
                plot_h = get(plot_layout.metadata, :topo_plot_height, 0.1)
                scale_pos = get(plot_layout.metadata, :topo_scale_pos, (0.8, -0.8))
                scale_hx = 0.0 <= scale_pos[1] <= 1.0 ? scale_pos[1] : clamp((scale_pos[1] + 1) / 2, 0.0, 1.0)
                scale_vy = 0.0 <= scale_pos[2] <= 1.0 ? scale_pos[2] : clamp((scale_pos[2] + 1) / 2, 0.0, 1.0)
                scale_ax = Axis(fig[1, 1], width = Relative(plot_w), height = Relative(plot_h), halign = scale_hx, valign = scale_vy)
                push!(axes, scale_ax)
                t_col = dat_subset[1].data[1][!, :time]
                _finalize_epochs_axis!(scale_ax, plot_kwargs; ylim = ylim, xlim = (t_col[1], t_col[end]))
            end
            length(axes) > 1 && Makie.linkaxes!(axes...)
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

            right_click_handler = (sel, mouse_x, data) -> _handle_waveform_right_click!(sel, mouse_x, data, condition_checked_ref)
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, right_click_handler)

            if plot_layout.type in (:topo, :grid, :custom)
                _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes)
            end
        end

        plot_kwargs[:display_plot] && _display_figure(fig)
        return (fig = fig, axes = axes)
    finally
        _set_window_title("Makie")
    end
end

plot_epochs(dat::EpochData; kwargs...) = plot_epochs([dat]; kwargs...)

function plot_epochs(filename::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    return _dispatch_plot_filepath(
        data -> plot_epochs(data; kwargs...),
        filename;
        input_dir = input_dir,
        participant_selection = participant_selection,
    )
end

######################################
######### Internal Functions #########
######################################

_something_or_calc_ylim(ylim, datasets, chs) = isnothing(ylim) ? _ylimits(datasets; channel_selection = channels(chs)) : ylim

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
    _setup_linked_legend_interactions_epochs!(line_refs)
    _setup_waveform_control_panel!(
        fig,
        dat_subset,
        baseline_interval,
        condition_checked_ref,
        legend_refs;
        panel_title = "Epochs Control Panel",
    ) do baseline_changed, baseline_interval_new, condition_mask
        if baseline_changed
            baseline!.(dat_subset, Ref(baseline_interval_new))
            baseline!.(dat_subset_avg, Ref(baseline_interval_new))
        end
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
                        baseline_changed && (trial_y_obs[] = _concatenate_trials(dat.data, ch, dat.data[1][!, :time]))
                    end
                    avg_entry = get(line_data, :average, nothing)
                    if !isnothing(avg_entry)
                        avg_line, y_obs = avg_entry
                        avg_line.visible[] != visible && (avg_line.visible = visible)
                        baseline_changed && cond_idx <= length(dat_subset_avg) && (y_obs[] = dat_subset_avg[cond_idx].data[!, ch])
                    end
                end
            end
        end
    end
end
