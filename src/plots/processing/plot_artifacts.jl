# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_ARTIFACT_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Display plot."),
    :plot_title => (nothing, "Custom plot title override."),
    :show_title => (true, "Whether to show the plot title."),
    :legend => (true, "Whether to display channel legend."),
    :colormap => (nothing, "Colormap (e.g., :Set1_9, :tab20)."),
    :linewidth_normal => (1, "Line width for normal channels."),
    :linewidth_rejected => (2, "Line width for rejected channels."),
    :alpha_normal => (0.2, "Transparency for normal channels."),
    :alpha_rejected => (1.0, "Transparency for rejected channels."),
    :selection_threshold => (50, "Distance threshold for channel selection."),
    :xlim => (nothing, "X-axis limits as (min, max) tuple or nothing for auto-scaling"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple or nothing for auto-scaling"),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),
)

const _EMPTY_CHANNEL_SET = Set{Symbol}()

function _create_artifact_controls(layout_pos, n_epochs::Int, n_artifact_epochs::Int)
    controls_layout = GridLayout(layout_pos, tellheight = false, tellwidth = false)
    back_button = Button(controls_layout[1, 1], label = "◀ Previous")
    epoch_label = Label(controls_layout[1, 2], "Epoch 1 / $n_epochs", justification = :center)
    forward_button = Button(controls_layout[1, 3], label = "Next ▶")
    back_artifact_button = Button(controls_layout[1, 4], label = "◀ Previous Artifact")
    art_str = n_artifact_epochs > 0 ? "Artifact 1 / $n_artifact_epochs" : "No artifacts found"
    artifact_label = Label(controls_layout[1, 5], art_str, justification = :center)
    forward_artifact_button = Button(controls_layout[1, 6], label = "Next Artifact ▶")
    return (; back_button, epoch_label, forward_button, back_artifact_button, artifact_label, forward_artifact_button)
end

function _index_rejected_artifacts(artifacts::EpochRejectionInfo, colormap_name)
    rejected_by_epoch = Dict{Int,Set{Symbol}}()
    rejected_color_map = Dict{Symbol,Any}()
    colors = Makie.to_colormap(colormap_name)
    n_colors = length(colors)
    for (i, r) in enumerate(artifacts.rejected)
        push!(get!(Set{Symbol}, rejected_by_epoch, r.epoch), r.channel)
        rejected_color_map[r.channel] = colors[((i-1)%n_colors)+1]
    end
    epochs_with_artifacts = sort!(collect(keys(rejected_by_epoch)))
    return rejected_by_epoch, epochs_with_artifacts, rejected_color_map
end

function _plot_artifact_channels!(
    ax::Axis,
    epoch::DataFrame,
    selected_channels::Vector{Symbol},
    rejected_set::Set{Symbol},
    selected_set::Set{Symbol},
    rejected_color_map::Dict{Symbol,Any},
    default_color,
    plot_kwargs::Dict{Symbol,Any};
    is_repaired::Bool = false,
)
    t = epoch.time
    alpha_rej = plot_kwargs[:alpha_rejected]
    alpha_norm = plot_kwargs[:alpha_normal]
    lw_rej = plot_kwargs[:linewidth_rejected]
    lw_norm = plot_kwargs[:linewidth_normal]

    for ch in selected_channels
        is_rej = ch in rejected_set
        is_sel = ch in selected_set
        col = is_rej ? rejected_color_map[ch] : default_color
        alpha = is_rej ? alpha_rej : alpha_norm
        lw = is_rej ? lw_rej : lw_norm
        is_sel && (lw *= 2)

        lbl_str = is_rej ? (is_repaired ? "$ch (repaired)" : "$ch (rejected)") : string(ch)
        lbl = is_sel ? rich(lbl_str, font = :bold) : lbl_str
        lines!(ax, t, epoch[!, ch]; color = col, alpha = alpha, linewidth = lw, label = lbl)
    end
end

function _reset_artifact_legend!(legend_ref::Base.RefValue{Union{Nothing,Legend}})
    if !isnothing(legend_ref[])
        delete!(legend_ref[])
        legend_ref[] = nothing
    end
end

function _set_artifact_legend!(ax::Axis, n_channels::Int, legend_ref::Base.RefValue{Union{Nothing,Legend}})
    n_cols = n_channels > 10 ? cld(n_channels, 20) : 1
    legend_ref[] = axislegend(ax, position = (:right, :top), nbanks = n_cols)
end

function _run_artifact_viewer(
    epochs_raw::EpochData,
    artifacts::EpochRejectionInfo,
    channel_selection,
    n_axes::Int,
    render_epoch!::Function;
    kwargs...,
)
    _set_window_title(_generate_window_title(epochs_raw))
    plot_kwargs = _merge_plot_kwargs(PLOT_ARTIFACT_KWARGS, kwargs)
    selected_channels = get_selected_channels(epochs_raw, channel_selection, include_meta = false, include_extra = false)

    fig = _eegfun_figure()
    axes = [Axis(fig[i, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)") for i = 1:n_axes]
    for ax in axes
        _set_axis_grid!(
            ax;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )
        _set_axis_properties!(ax; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])
        _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
    end
    n_axes == 2 && linkxaxes!(axes[1], axes[2])

    actual_cmap = _resolve_theme_colormap(axes[1], plot_kwargs[:colormap])
    rejected_by_epoch, epochs_with_artifacts, rejected_color_map = _index_rejected_artifacts(artifacts, actual_cmap)
    default_color = _resolve_theme_linecolor(axes[1], nothing)

    n_eps = n_epochs(epochs_raw)
    controls = _create_artifact_controls(fig[n_axes+1, 1], n_eps, length(epochs_with_artifacts))
    plot_rel = n_axes == 1 ? 0.9 : 0.45
    for i = 1:n_axes
        rowsize!(fig.layout, i, Relative(plot_rel))
    end
    rowsize!(fig.layout, n_axes + 1, Relative(0.1))

    epoch_idx = Observable(1)
    artifact_idx = Observable(1)
    legends = [Ref{Union{Nothing,Legend}}(nothing) for _ = 1:n_axes]
    selected_set = Set{Symbol}()

    refresh! =
        (idx::Int) -> begin
            for j = 1:n_axes
                empty!(axes[j])
                _reset_artifact_legend!(legends[j])
            end
            render_epoch!(
                idx,
                axes,
                legends,
                selected_channels,
                selected_set,
                rejected_by_epoch,
                rejected_color_map,
                default_color,
                plot_kwargs,
            )
        end

    sel_thresh = plot_kwargs[:selection_threshold]
    for ax in axes
        on(events(ax).mousebutton, priority = 0) do event
            if event.button == Mouse.left && event.action == Mouse.press && _is_shift_held(fig)
                mouse_pos = _get_mouse_position(ax)
                isnothing(mouse_pos) && return Consume(false)
                closest_ch, min_dist = _find_closest_channel(mouse_pos[1], mouse_pos[2], epochs_raw.data[epoch_idx[]], selected_channels)
                if !isnothing(closest_ch) && min_dist < sel_thresh
                    _toggle_channel_selection(closest_ch, selected_set)
                    refresh!(epoch_idx[])
                end
            end
            return Consume(false)
        end
    end

    on(controls.back_button.clicks) do _
        epoch_idx[] = max(1, epoch_idx[] - 1)
    end
    on(controls.forward_button.clicks) do _
        epoch_idx[] = min(n_eps, epoch_idx[] + 1)
    end
    on(controls.back_artifact_button.clicks) do _
        artifact_idx[] = max(1, artifact_idx[] - 1)
    end
    on(controls.forward_artifact_button.clicks) do _
        artifact_idx[] = min(length(epochs_with_artifacts), artifact_idx[] + 1)
    end

    on(epoch_idx) do idx
        refresh!(idx)
        controls.epoch_label.text = "Epoch $idx / $n_eps"
    end
    on(artifact_idx) do idx
        if !isempty(epochs_with_artifacts)
            epoch_idx[] = epochs_with_artifacts[idx]
            controls.artifact_label.text = "Artifact $idx / $(length(epochs_with_artifacts))"
        end
    end

    refresh!(epoch_idx[])
    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

"""
    plot_artifact_detection(epochs::EpochData, artifacts::EpochRejectionInfo; channel_selection=channels(), kwargs...)
    plot_artifact_detection(epochs::AbstractVector{<:EpochData}, artifacts::AbstractVector{<:EpochRejectionInfo}; kwargs...)

Interactive plot of epoch artifact detection results with epoch and artifact navigation controls.

$(_generate_kwargs_doc(PLOT_ARTIFACT_KWARGS))
"""
function plot_artifact_detection(epochs::EpochData, artifacts::EpochRejectionInfo; channel_selection = channels(), kwargs...)
    return _run_artifact_viewer(
        epochs,
        artifacts,
        channel_selection,
        1,
        (idx, axes, legends, sel_chs, sel_set, rej_by_ep, rej_cmap, def_col, pkw) -> begin
            ax = axes[1]
            _set_origin_lines!(ax; add_xy_origin = pkw[:add_xy_origin])
            rej_set = get(rej_by_ep, idx, _EMPTY_CHANNEL_SET)
            ax.title = !pkw[:show_title] ? "" : something(pkw[:plot_title], "Artifact Detection - Epoch $idx")
            _plot_artifact_channels!(ax, epochs.data[idx], sel_chs, rej_set, sel_set, rej_cmap, def_col, pkw)
            pkw[:legend] && _set_artifact_legend!(ax, length(sel_chs), legends[1])
        end;
        kwargs...,
    )
end

function plot_artifact_detection(
    epochs::AbstractVector{<:EpochData},
    artifacts::AbstractVector{<:EpochRejectionInfo};
    channel_selection = channels(),
    kwargs...,
)
    length(epochs) == length(artifacts) ||
        @minimal_error("epochs and artifacts must have the same length (got $(length(epochs)) and $(length(artifacts)))")
    return [plot_artifact_detection(epochs[i], artifacts[i]; channel_selection = channel_selection, kwargs...) for i in eachindex(epochs)]
end

"""
    plot_artifact_repair(epochs_raw::EpochData, epochs_repaired::EpochData, artifacts::EpochRejectionInfo; channel_selection=channels(), kwargs...)

Interactive before/after comparison of original and repaired epochs with navigation controls.

$(_generate_kwargs_doc(PLOT_ARTIFACT_KWARGS))
"""
function plot_artifact_repair(
    epochs_raw::EpochData,
    epochs_repaired::EpochData,
    artifacts::EpochRejectionInfo;
    channel_selection = channels(),
    kwargs...,
)
    return _run_artifact_viewer(
        epochs_raw,
        artifacts,
        channel_selection,
        2,
        (idx, axes, legends, sel_chs, sel_set, rej_by_ep, rej_cmap, def_col, pkw) -> begin
            ax1, ax2 = axes[1], axes[2]
            _set_origin_lines!(ax1; add_xy_origin = pkw[:add_xy_origin])
            _set_origin_lines!(ax2; add_xy_origin = pkw[:add_xy_origin])
            rej_set = get(rej_by_ep, idx, _EMPTY_CHANNEL_SET)
            ax1.title = pkw[:show_title] ? "Original Data - Epoch $idx" : ""
            ax2.title = pkw[:show_title] ? "Repaired Data - Epoch $idx" : ""
            _plot_artifact_channels!(ax1, epochs_raw.data[idx], sel_chs, rej_set, sel_set, rej_cmap, def_col, pkw)
            _plot_artifact_channels!(ax2, epochs_repaired.data[idx], sel_chs, rej_set, sel_set, rej_cmap, def_col, pkw; is_repaired = true)
            if pkw[:legend]
                _set_artifact_legend!(ax1, length(sel_chs), legends[1])
                _set_artifact_legend!(ax2, length(sel_chs), legends[2])
            end
        end;
        kwargs...,
    )
end

"""
    plot_artifact_rejection(epochs_raw::EpochData, epochs_rejected::EpochData, artifacts::EpochRejectionInfo; channel_selection=channels(), kwargs...)

Interactive before/after comparison between original and post-rejection epochs aligned by epoch number.
Rejected epochs display a blank lower panel with red axis spines.

$(_generate_kwargs_doc(PLOT_ARTIFACT_KWARGS))
"""
function plot_artifact_rejection(
    epochs_raw::EpochData,
    epochs_rejected::EpochData,
    artifacts::EpochRejectionInfo;
    channel_selection = channels(),
    kwargs...,
)
    epoch_number_to_idx = Dict{Int,Int}(ep_df.epoch[1] => i for (i, ep_df) in enumerate(epochs_rejected.data))
    return _run_artifact_viewer(
        epochs_raw,
        artifacts,
        channel_selection,
        2,
        (idx, axes, legends, sel_chs, sel_set, rej_by_ep, rej_cmap, def_col, pkw) -> begin
            ax1, ax2 = axes[1], axes[2]
            _set_origin_lines!(ax1; add_xy_origin = pkw[:add_xy_origin])
            ep_orig = epochs_raw.data[idx]
            ep_num = ep_orig.epoch[1]
            rej_set = get(rej_by_ep, ep_num, _EMPTY_CHANNEL_SET)
            rej_idx = get(epoch_number_to_idx, ep_num, nothing)

            ax1.title = pkw[:show_title] ? "Original Data - Epoch $ep_num" : ""
            ax2.title = pkw[:show_title] ? "Rejected Data - Epoch $ep_num" : ""

            _plot_artifact_channels!(ax1, ep_orig, sel_chs, rej_set, sel_set, rej_cmap, def_col, pkw)
            if !isnothing(rej_idx)
                _set_axis_spine_color!(ax1, def_col)
                _set_axis_spine_color!(ax2, def_col)
                _set_origin_lines!(ax2; add_xy_origin = pkw[:add_xy_origin])
                _plot_artifact_channels!(ax2, epochs_rejected.data[rej_idx], sel_chs, rej_set, sel_set, rej_cmap, def_col, pkw)
                if pkw[:legend]
                    _set_artifact_legend!(ax1, length(sel_chs), legends[1])
                    _set_artifact_legend!(ax2, length(sel_chs), legends[2])
                end
            else
                _set_axis_spine_color!(ax1, :red)
                _set_axis_spine_color!(ax2, :red)
                pkw[:legend] && _set_artifact_legend!(ax1, length(sel_chs), legends[1])
            end
        end;
        kwargs...,
    )
end
