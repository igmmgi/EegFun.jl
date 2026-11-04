"""
    plot_erp_stats(result::StatsResult;
                   layout::Union{Symbol,PlotLayout} = :single,
                   condition_selection = conditions(),
                   channel_selection = channels(),
                   channel_plot_order = nothing,
                   plot_erp::Bool = true,
                   plot_difference::Bool = false,
                   plot_tvalues::Bool = false,
                   plot_significance::Bool = false,
                   plot_critical_t::Bool = false,
                   plot_se::Bool = false,
                   difference_offset::Real = 0.0,
                   significance_position::Union{Symbol,Real} = :auto,
                   significance_color = (:gray, 0.6),
                   kwargs...)

Plot ERP waveforms and statistical results for `AnalyticResult` or `PermutationResult`
across `:single`, `:grid`, `:topo`, or custom `PlotLayout` layouts.

# Arguments
- `result::StatsResult`: Result from `analytic_test` or `permutation_test`.
- `layout`: `:single` (default), `:grid`, `:topo`, or a `PlotLayout`.
- `condition_selection`: Predicate to select conditions when `plot_erp = true`.
- `channel_selection`: Predicate, symbol(s), or name(s) selecting channels (default: `channels()`).
- `channel_plot_order`: Optional vector overriding channel plot order.
- `plot_erp`: Plot condition average waveforms (default: `true`).
- `plot_difference`: Plot difference wave `A - B` (default: `false`).
- `plot_tvalues`: Plot t-statistics (default: `false`; auto-enabled when `plot_critical_t = true`).
- `plot_significance`: Highlight significant time intervals (default: `false`).
- `plot_critical_t`: Plot critical t thresholds (default: `false`).
- `plot_se`: Plot ±1 SEM bands around condition and difference waveforms (default: `false`).
- `difference_offset`: Vertical offset for the difference wave (default: `0.0`).
- `significance_position`: `:auto` (default), `:zero`, `:bottom`, or custom `Real` y-coordinate.
- `significance_color`: Color for significance bars (default: `(:gray, 0.6)`).

$(_generate_kwargs_doc(PLOT_ERP_KWARGS))

# Returns
- Named tuple `(fig, axes)`
"""
function plot_erp_stats(
    result::StatsResult;
    layout::Union{Symbol,PlotLayout} = :single,
    condition_selection = conditions(),
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    plot_erp::Bool = true,
    plot_difference::Bool = false,
    plot_tvalues::Bool = false,
    plot_significance::Bool = false,
    plot_critical_t::Bool = false,
    plot_se::Bool = false,
    difference_offset::Real = 0.0,
    significance_position::Union{Symbol,Real} = :auto,
    significance_color = (:gray, 0.6),
    kwargs...,
)
    plot_tvalues = plot_tvalues || plot_critical_t
    (plot_erp || plot_difference || plot_tvalues) ||
        @minimal_error("At least one of plot_erp, plot_difference, or plot_tvalues must be true")

    kwargs_norm = _normalize_plot_kwargs(kwargs)
    haskey(kwargs_norm, :show_significance) && (plot_significance = Bool(pop!(kwargs_norm, :show_significance)))
    haskey(kwargs_norm, :significance_position) && (significance_position = pop!(kwargs_norm, :significance_position))
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs_norm)
    plot_kwargs[:window_title] = get(kwargs_norm, :window_title, "ERP Stats")

    default_stats_colors = Any[:blue, :red, :black, :purple]
    user_color = plot_kwargs[:color]
    stats_colors = if user_color isa Vector
        Any[i <= length(user_color) ? user_color[i] : default_stats_colors[i] for i = 1:4]
    elseif haskey(kwargs_norm, :color) && !isnothing(user_color)
        Any[user_color, user_color, :black, :purple]
    else
        default_stats_colors
    end

    all_electrodes = result.electrodes
    selected_channels = all_electrodes[_to_channel_func(channel_selection)(all_electrodes)]
    isempty(selected_channels) && @minimal_error("No channels matched the selection. Available channels: $(all_electrodes)")

    if !isnothing(channel_plot_order)
        ordered = Symbol[Symbol(ch) for ch in channel_plot_order if Symbol(ch) in selected_channels]
        isempty(ordered) && @minimal_error(
            "No channels in channel_plot_order matched the selection. Selected: $(selected_channels), Requested order: $(channel_plot_order)",
        )
        selected_channels = ordered
    end

    title_suffix = if result isa AnalyticResult
        "Analytic Test ($(result.test_info.correction_method))"
    else
        "Permutation Test ($(result.test_info.cluster_info.threshold_method))"
    end

    plot_kwargs[:ylabel] = if plot_tvalues && (plot_erp || plot_difference)
        "Amplitude (μV) / t-statistic"
    elseif plot_tvalues
        "t-statistic"
    else
        "Amplitude (μV)"
    end

    if layout == :single && isnothing(plot_kwargs[:plot_title])
        ch_label = length(selected_channels) == 1 ? string(selected_channels[1]) : _print_vector(selected_channels)
        plot_kwargs[:plot_title] = "$ch_label - $title_suffix"
    end

    layout_kwargs = pop!(plot_kwargs, :layout_kwargs, (;))
    fig = _eegfun_figure(size = (800, 600), title = plot_kwargs[:window_title], figure_padding = plot_kwargs[:figure_padding])

    plot_layout = create_layout(layout, selected_channels, result.data[1].layout; layout_kwargs...)
    axes, layout_channels = _apply_layout!(fig, plot_layout; plot_kwargs...)

    cond_indices = plot_erp ? get_selected_conditions(result.data, condition_selection) : Int[]
    elec_index_map = Dict{Symbol,Int}(ch => i for (i, ch) in enumerate(all_electrodes))

    for (ax, channel) in zip(axes, layout_channels)
        channels_to_plot = plot_layout.type == :single ? selected_channels : Symbol[channel]
        first_ch = first(channels_to_plot)
        for ch in channels_to_plot
            channel_idx = get(elec_index_map, ch, 0)
            channel_idx == 0 && continue
            _plot_erp_stats_channel!(
                ax,
                result,
                channel_idx,
                ch,
                cond_indices;
                plot_erp = plot_erp,
                plot_difference = plot_difference,
                plot_tvalues = plot_tvalues,
                plot_significance = plot_significance,
                plot_critical_t = plot_critical_t,
                plot_se = plot_se,
                colors = stats_colors,
                difference_offset = difference_offset,
                significance_position = significance_position,
                significance_color = significance_color,
                linewidth = plot_kwargs[:linewidth],
                legend_labels = plot_kwargs[:legend_labels],
                add_labels = (ch === first_ch),
            )
        end
    end

    for ax in axes
        _apply_axis_properties!(ax; plot_kwargs...)
    end
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...)

    if plot_layout.type in (:grid, :topo) && plot_kwargs[:show_title]
        for (ax, channel) in zip(axes, layout_channels)
            ax.title = string(channel)
        end
    end

    length(axes) > 1 && linkaxes!(axes...)

    if plot_kwargs[:add_xy_origin]
        for ax in axes
            _set_origin_lines!(ax)
        end
    end

    if !isempty(axes) && plot_kwargs[:legend]
        axislegend(axes[1], position = plot_kwargs[:legend_position], framevisible = plot_kwargs[:legend_framevisible])
    end

    if !isempty(plot_kwargs[:figure_title])
        Label(fig[0, :], plot_kwargs[:figure_title], fontsize = plot_kwargs[:figure_title_fontsize], font = :bold, tellwidth = false)
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = axes)
end

function _plot_erp_stats_channel!(
    ax::Axis,
    result::StatsResult,
    channel_idx::Int,
    channel_sym::Symbol,
    cond_indices::Vector{Int};
    plot_erp::Bool = true,
    plot_difference::Bool = false,
    plot_tvalues::Bool = false,
    plot_significance::Bool = false,
    plot_critical_t::Bool = false,
    plot_se::Bool = false,
    colors::Vector = Any[:blue, :red, :black, :purple],
    difference_offset::Real = 0.0,
    significance_position::Union{Symbol,Real} = :auto,
    significance_color = (:gray, 0.6),
    linewidth = 2,
    legend_labels::Vector = [],
    add_labels::Bool = true,
)
    time_points = result.time_points
    t_values = @view result.stat_matrix.t[channel_idx, :]
    lw = _resolve_theme_linewidth(ax, linewidth, 2)

    cond_A_erp = result.data[1]
    cond_B_erp = result.data[2]
    channel_sym ∉ cond_A_erp && @minimal_error("Channel $channel_sym not found in result data")
    cond_A_avg = cond_A_erp.data[!, channel_sym]
    cond_B_avg = cond_B_erp.data[!, channel_sym]
    erp_time_points = cond_A_erp.data[!, :time]

    critical_t_pos = nothing
    critical_t_neg = nothing
    if plot_critical_t && plot_tvalues
        crit = result.critical_t
        n_tp = length(time_points)
        if crit isa Real
            critical_t_pos = fill(Float64(crit), n_tp)
            critical_t_neg = fill(-Float64(crit), n_tp)
        elseif crit isa Tuple{<:Real,<:Real}
            critical_t_pos = fill(Float64(crit[1]), n_tp)
            critical_t_neg = fill(Float64(crit[2]), n_tp)
        elseif crit isa Tuple{<:AbstractMatrix,<:AbstractMatrix}
            critical_t_pos = Float64.(@view crit[1][channel_idx, :])
            critical_t_neg = Float64.(@view crit[2][channel_idx, :])
        elseif crit isa AbstractMatrix
            critical_t_pos = Float64.(@view crit[channel_idx, :])
            critical_t_neg = .-critical_t_pos
        end
    end

    if plot_erp
        for (i, orig_idx) in enumerate(cond_indices)
            erp = result.data[orig_idx]
            erp_name = (isempty(legend_labels) || i > length(legend_labels)) ? erp.condition_name : legend_labels[i]
            erp_avg = erp.data[!, channel_sym]
            lines!(ax, erp_time_points, erp_avg, color = colors[i], linewidth = lw, label = add_labels ? erp_name : nothing)
            if plot_se && (orig_idx == 1 || orig_idx == 2)
                se_mat = orig_idx == 1 ? result.se_cond1 : result.se_cond2
                se_data = @view se_mat[channel_idx, :]
                band!(ax, erp_time_points, erp_avg .- se_data, erp_avg .+ se_data, color = (colors[i], 0.15))
            end
        end
    end

    diff_wave_plot = nothing
    if plot_difference || plot_significance
        diff_wave = cond_A_avg .- cond_B_avg
        diff_wave_plot = difference_offset != 0.0 && plot_difference ? (diff_wave .+ difference_offset) : diff_wave
    end

    if plot_difference
        diff_label = if difference_offset != 0.0
            hlines!(ax, [difference_offset], color = (:gray, 0.7), linewidth = 1)
            text!(
                ax,
                erp_time_points[end] * 0.98,
                difference_offset,
                text = "0 μV",
                align = (:right, :center),
                color = :gray,
                fontsize = 18,
            )
            "Difference (offset=$(difference_offset))"
        else
            length(legend_labels) >= 3 ? legend_labels[3] : "Difference"
        end
        lines!(ax, erp_time_points, diff_wave_plot, color = colors[3], linewidth = lw, label = add_labels ? diff_label : nothing)
        if plot_se
            se_diff_channel = @view result.se_diff[channel_idx, :]
            band!(
                ax,
                erp_time_points,
                diff_wave_plot .- se_diff_channel,
                diff_wave_plot .+ se_diff_channel,
                color = (colors[3], 0.15),
                label = "±1 SEM",
            )
        end
    end

    plot_tvalues && lines!(ax, time_points, t_values, color = colors[4], linewidth = lw, label = "t-statistic")

    if plot_significance
        sig_pos = @view result.masks.positive[channel_idx, :]
        sig_neg = @view result.masks.negative[channel_idx, :]
        if any(sig_pos) || any(sig_neg)
            sig_regions = _find_continuous_regions(sig_pos .| sig_neg, time_points)
            bar_y, bar_height = _compute_significance_position(
                significance_position,
                plot_erp,
                plot_difference,
                plot_tvalues,
                cond_A_avg,
                cond_B_avg,
                diff_wave_plot,
                t_values,
            )
            for (t_start, t_end) in sig_regions
                poly!(
                    ax,
                    Point2f[
                        Point2f(t_start, bar_y),
                        Point2f(t_end, bar_y),
                        Point2f(t_end, bar_y + bar_height),
                        Point2f(t_start, bar_y + bar_height),
                    ],
                    color = significance_color,
                    strokewidth = 0,
                )
            end
        end
    end

    if plot_critical_t && plot_tvalues && !isnothing(critical_t_pos)
        any(!isnan, critical_t_pos) &&
            lines!(ax, time_points, critical_t_pos, color = :grey, linewidth = lw, linestyle = :dashdot, label = "Critical t+")
        !isnothing(critical_t_neg) &&
            any(!isnan, critical_t_neg) &&
            lines!(ax, time_points, critical_t_neg, color = :grey, linewidth = lw, linestyle = :dashdot, label = "Critical t-")
    end

    return nothing
end

function _nonnan_extrema(v::AbstractVector)
    lo = Inf
    hi = -Inf
    @inbounds for x in v
        val = Float64(x)
        if !isnan(val)
            val < lo && (lo = val)
            val > hi && (hi = val)
        end
    end
    return lo <= hi ? (lo, hi) : (0.0, 1.0)
end

function _amp_extrema(plot_erp::Bool, plot_difference::Bool, cond_A_avg, cond_B_avg, diff_wave_plot)
    if plot_erp
        a_min, a_max = extrema(cond_A_avg)
        b_min, b_max = extrema(cond_B_avg)
        return min(a_min, b_min), max(a_max, b_max)
    elseif plot_difference && !isnothing(diff_wave_plot)
        return extrema(diff_wave_plot)
    end
    return 0.0, 0.0
end

function _compute_significance_position(position, plot_erp, plot_difference, plot_tvalues, cond_A_avg, cond_B_avg, diff_wave_plot, t_values)
    data_range, data_min = if plot_tvalues && !plot_difference && !plot_erp
        t_min, t_max = _nonnan_extrema(t_values)
        t_max - t_min, t_min
    elseif plot_erp || plot_difference
        amp_min, amp_max = _amp_extrema(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
        amp_max - amp_min, amp_min
    else
        1.0, 0.0
    end

    bar_height = data_range * 0.05
    bar_height ≈ 0.0 && (bar_height = 0.05)

    bar_y = if position isa Real
        Float64(position)
    elseif position == :zero
        -bar_height / 2
    elseif position == :bottom
        data_min - data_range * 0.05
    elseif plot_erp || plot_difference
        amp_min, amp_max = _amp_extrema(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
        amp_min <= 0.0 <= amp_max ? (-bar_height / 2) : (amp_min - data_range * 0.05)
    else
        data_min - data_range * 0.05
    end
    return bar_y, bar_height
end

function _max_abs_nonnan(A::AbstractArray)
    m = 0.0
    @inbounds for x in A
        v = abs(Float64(x))
        if !isnan(v) && v > m
            m = v
        end
    end
    return m ≈ 0.0 ? 1.0 : m
end

function _transpose_f32(t_values::AbstractMatrix{<:Real})
    n_chans, n_times = size(t_values)
    out = Matrix{Float32}(undef, n_times, n_chans)
    @inbounds for c = 1:n_chans, t = 1:n_times
        out[t, c] = Float32(t_values[c, t])
    end
    return out
end

function _transpose_masked_f32(t_values::AbstractMatrix{<:Real}, mask::AbstractMatrix{Bool})
    n_chans, n_times = size(t_values)
    out = Matrix{Float32}(undef, n_times, n_chans)
    @inbounds for c = 1:n_chans, t = 1:n_times
        out[t, c] = mask[c, t] ? Float32(t_values[c, t]) : NaN32
    end
    return out
end

function _plot_stat_heatmap_core!(
    fig::Union{Figure,Nothing},
    ax::Axis,
    time_points::AbstractVector{<:Real},
    electrodes::AbstractVector,
    t_values::AbstractMatrix{<:Real},
    mask::AbstractMatrix{Bool};
    kwargs...,
)
    kwargs_norm = _normalize_plot_kwargs(kwargs)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs_norm; validate = false)

    time_unit = get(kwargs_norm, :time_unit, :s)
    t_scale = time_unit == :ms ? 1000.0 : 1.0
    x_times = time_points .* t_scale

    n_chans = length(electrodes)
    plot_kwargs[:ylabel] = get(kwargs_norm, :ylabel, "Electrode (1 to $(n_chans))")
    plot_kwargs[:xlabel] = get(kwargs_norm, :xlabel, time_unit == :ms ? "Time (ms)" : "Time (s)")
    pop!(plot_kwargs, :time_unit, nothing)

    if isnothing(get(kwargs_norm, :yticks, nothing))
        if n_chans > 40
            step = max(1, div(n_chans, 14))
            y_indices = [1; step:step:n_chans]
            plot_kwargs[:yticks] = (y_indices, ["$(i): $(electrodes[i])" for i in y_indices])
            ax.yminorticks = IntervalsBetween(step)
            ax.yminorticksvisible = true
        else
            plot_kwargs[:yticks] = (1:n_chans, String.(electrodes))
        end
    end

    baseline_window = get(kwargs_norm, :baseline_window, nothing)
    if !isnothing(baseline_window)
        vspan!(ax, baseline_window[1] * t_scale, baseline_window[2] * t_scale, color = (:gray85, 0.25))
        vlines!(ax, [0.0], color = (:black, 0.75), linestyle = :solid, linewidth = 1.2)
    end

    max_t = _max_abs_nonnan(t_values)
    cmap = _resolve_theme_colormap(ax, get(kwargs_norm, :colormap, :coolwarm))
    clims = something(get(kwargs_norm, :colorrange, nothing), (-max_t, max_t))
    alpha_subdued = get(kwargs_norm, :alpha_subdued, 0.20)
    mask_style = get(kwargs_norm, :mask_style, :alpha)

    hm = if mask_style == :alpha
        heatmap!(ax, x_times, 1:n_chans, _transpose_f32(t_values), colormap = cmap, colorrange = clims, alpha = alpha_subdued)
        heatmap!(
            ax,
            x_times,
            1:n_chans,
            _transpose_masked_f32(t_values, mask),
            colormap = cmap,
            colorrange = clims,
            nan_color = :transparent,
        )
    elseif mask_style == :hide
        heatmap!(ax, x_times, 1:n_chans, _transpose_masked_f32(t_values, mask), colormap = cmap, colorrange = clims, nan_color = :transparent)
    elseif mask_style == :contour
        h = heatmap!(ax, x_times, 1:n_chans, _transpose_f32(t_values), colormap = cmap, colorrange = clims)
        any(mask) && contour!(ax, x_times, 1:n_chans, _transpose_f32(mask), levels = [0.5], color = :black, linewidth = 1.8)
        h
    else
        heatmap!(ax, x_times, 1:n_chans, _transpose_f32(t_values), colormap = cmap, colorrange = clims)
    end

    xlims!(ax, extrema(x_times))
    ylims!(ax, 0.5, n_chans + 0.5)

    if get(kwargs_norm, :colorbar_plot, true) && !isnothing(fig)
        Colorbar(fig[1, 2], hm, label = get(kwargs_norm, :colorbar_label, "t-value"))
    end

    _apply_axis_properties!(ax; plot_kwargs...)
    return (hm = hm, ax = ax)
end

"""
    plot_stat_heatmap!(fig, ax, result::StatsResult; kwargs...)

Mutating variant of `plot_stat_heatmap` to draw a 2D statistical heatmap into an existing `Axis`.
"""
function plot_stat_heatmap!(fig::Union{Figure,Nothing}, ax::Axis, result::StatsResult; kwargs...)
    t_values = result.stat_matrix.t
    has_masks = hasproperty(result, :masks) && !isnothing(result.masks)
    mask = has_masks ? (result.masks.positive .| result.masks.negative) : trues(size(t_values))
    return _plot_stat_heatmap_core!(fig, ax, result.time_points, result.electrodes, t_values, mask; kwargs...)
end

"""
    plot_stat_heatmap(result::StatsResult; kwargs...)

Plot a 2D heatmap of t-statistics (`Channels × Time`) from an analytic or permutation test.
"""
function plot_stat_heatmap(result::StatsResult; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs; validate = false)
    window_title = get(kwargs, :window_title, "T-Statistic Heatmap")
    xlabel = get(kwargs, :xlabel, get(kwargs, :time_unit, :s) == :ms ? "Time (ms)" : "Time (s)")
    ylabel = get(kwargs, :ylabel, "Channels")

    fig = _eegfun_figure(size = (850, 600), title = window_title, figure_padding = plot_kwargs[:figure_padding])
    ax = Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel)

    plot_stat_heatmap!(fig, ax, result; kwargs...)

    if !isempty(plot_kwargs[:figure_title])
        Label(fig[0, :], plot_kwargs[:figure_title], fontsize = plot_kwargs[:figure_title_fontsize], font = :bold, tellwidth = false)
    end
    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax])
end
