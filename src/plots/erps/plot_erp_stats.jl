"""
Plotting functions for statistical test results (analytic and permutation tests).
"""


"""
    plot_erp_stats(result::StatsResult;
                    layout::Union{Symbol, PlotLayout} = :single,
                    channel_selection::Function = channels(),
                    channel_plot_order::Union{Nothing, Vector{Symbol}} = nothing,
                    plot_erp::Bool = true,
                    plot_difference::Bool = false,
                    plot_tvalues::Bool = false,
                    plot_significance::Bool = false,
                    plot_critical_t::Bool = false,
                    plot_se::Bool = false,
                    colors::Vector = [:blue, :red, :black, :purple],
                    difference_offset::Real = 0.0,
                    significance_position::Union{Symbol, Real} = :auto,
                    significance_color = (:gray, 0.6),
                    kwargs...)

Plot ERP waveforms and statistical results for analytic or permutation tests with flexible layout support.

Works with both `AnalyticResult` (from `analytic_test`) and `PermutationResult` (from `permutation_test`).

# Arguments
- `result::StatsResult`: Results from `analytic_test` or `permutation_test`
- `layout`: Layout specification:
  - `:single` (default): Single plot with selected channels overlaid
  - `:grid`: Grid layout with one subplot per selected channel
  - `:topo`: Topographic layout based on channel positions
  - `PlotLayout`: Custom layout object
- `channel_selection::Function`: Predicate to select channels (default: `channels()` - all channels)
- `channel_plot_order::Union{Nothing, Vector{Symbol}}`: Override the plotting order of selected channels (default: `nothing` — data order).
  When provided, channels are plotted in the specified order. Only channels present in both `channel_plot_order` and the selection are plotted
- `plot_erp::Bool`: Whether to plot ERP waveforms (condition averages) (default: true)
- `plot_difference::Bool`: Whether to plot difference wave (A-B) (default: false)
- `plot_tvalues::Bool`: Whether to plot t-statistics (default: false)
- `plot_significance::Bool`: Whether to highlight significant time points (default: false)
- `plot_critical_t::Bool`: Whether to plot critical t-values (default: false). Only relevant when `plot_tvalues=true`
- `difference_offset::Real`: Vertical offset for difference wave (default: 0.0). Set to non-zero to shift for visibility
- `plot_se::Bool`: Whether to plot ±1 SEM bands around waveforms (default: false). Bands are drawn around condition ERPs and difference wave
- `significance_position::Union{Symbol, Real}`: Position for significance bars (default: `:auto`). Options:
  - `:auto` - Automatically place at y=0 if visible, otherwise at bottom (default)
  - `:zero` - Always place at y=0
  - `:bottom` - Always place at bottom of plot/spine
  - `Real` - Custom y-position (e.g., `-5.0` to place at -5 μV)
- `significance_color`: Color for significance bars (default: `(:gray, 0.6)`). Can be any Makie color specification

$(_generate_kwargs_doc(PLOT_ERP_KWARGS))

# Returns
- Named tuple `(fig, axes)` 

# Examples
```julia
# With analytic test results - single channel
result = analytic_test(prepared, correction_method=:no)
plot_erp_stats(result, channel_selection=channels(:PO7), 
               plot_erp=true, plot_difference=true, plot_significance=true)

# Grid layout with multiple channels
plot_erp_stats(result, channel_selection=channels([:PO7, :PO8, :Oz, :Pz]),
               layout=:grid, plot_significance=true)

# With permutation test results
result_perm = permutation_test(prepared, n_permutations=1000)
plot_erp_stats(result_perm, channel_selection=channels(:PO7),
               plot_erp=true, plot_significance=true, plot_critical_t=true)
```
"""
function plot_erp_stats(
    result::StatsResult;
    layout::Union{Symbol,PlotLayout}=:single,
    condition_selection::Function=conditions(),
    channel_selection::Function=channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}}=nothing,
    plot_erp::Bool=true,
    plot_difference::Bool=false,
    plot_tvalues::Bool=false,
    plot_significance::Bool=false,
    plot_critical_t::Bool=false,
    plot_se::Bool=false,
    difference_offset::Real=0.0,
    significance_position::Union{Symbol,Real}=:auto,
    significance_color=(:gray, 0.6),
    kwargs...,
)
    # If plot_critical_t is requested, automatically enable t-value plotting
    plot_tvalues = plot_tvalues || plot_critical_t

    # Validate that at least something is being plotted
    if !plot_erp && !plot_difference && !plot_tvalues
        error("At least one of plot_erp, plot_difference, or plot_tvalues must be true")
    end

    # Prepare kwargs (reuse PLOT_ERP_KWARGS, override stats-specific defaults)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)
    plot_kwargs[:window_title] = get(kwargs, :window_title, "ERP Stats")

    # Compute colors for stats plot elements: [cond1, cond2, difference, t-values]
    default_stats_colors = [:blue, :red, :black, :purple]
    user_color = plot_kwargs[:color]
    stats_colors = if user_color isa Vector
        # Pad with defaults if user provides fewer than 4
        [i <= length(user_color) ? user_color[i] : default_stats_colors[i] for i = 1:4]
    elseif haskey(kwargs, :color)
        # User provided a single color — use it for both conditions
        [user_color, user_color, :black, :purple]
    else
        default_stats_colors
    end

    # Determine selected channels from result electrodes
    all_electrodes = result.electrodes
    selected_mask = channel_selection(all_electrodes)
    selected_channels = all_electrodes[selected_mask]

    if isempty(selected_channels)
        error("No channels matched the selection. Available channels: $(all_electrodes)")
    end

    # Apply user-specified plotting order if provided
    if !isnothing(channel_plot_order)
        ordered = Symbol[ch for ch in channel_plot_order if ch in selected_channels]
        if isempty(ordered)
            error(
                "No channels in channel_plot_order matched the selection. Selected: $(selected_channels), Requested order: $(channel_plot_order)",
            )
        end
        selected_channels = ordered
    end

    # Set default title
    title_suffix = if isa(result, AnalyticResult)
        "Analytic Test ($(result.test_info.correction_method))"
    else
        "Permutation Test ($(result.test_info.cluster_info.threshold_method))"
    end

    # Determine y-axis label
    if plot_tvalues && (plot_erp || plot_difference)
        plot_kwargs[:ylabel] = "Amplitude (μV) / t-statistic"
    elseif plot_tvalues
        plot_kwargs[:ylabel] = "t-statistic"
    else
        plot_kwargs[:ylabel] = "Amplitude (μV)"
    end

    # Set title for single layout
    if layout == :single && isnothing(plot_kwargs[:plot_title])
        if length(selected_channels) == 1
            plot_kwargs[:plot_title] = "$(selected_channels[1]) - $title_suffix"
        else
            plot_kwargs[:plot_title] = "$(_print_vector(selected_channels)) - $title_suffix"
        end
    end

    # Extract layout_* parameters for layout system
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    # Get the EEG layout from the result data (for topo positioning)
    eeg_layout = result.data[1].layout

    # Create figure and apply layout system
    fig = Figure(size=(800, 600), title=plot_kwargs[:window_title], figure_padding=plot_kwargs[:figure_padding])



    plot_layout = create_layout(layout, selected_channels, eeg_layout; layout_kwargs...)
    axes, layout_channels = _apply_layout!(fig, plot_layout; plot_kwargs...)

    # Render stats content on each axis
    for (ax_idx, (ax, channel)) in enumerate(zip(axes, layout_channels))
        channels_to_plot = plot_layout.type == :single ? selected_channels : [channel]

        for ch in channels_to_plot
            channel_idx = findfirst(==(ch), all_electrodes)
            isnothing(channel_idx) && continue

            _plot_erp_stats_channel!(
                ax,
                result,
                channel_idx,
                ch;
                condition_selection=condition_selection,
                plot_erp=plot_erp,
                plot_difference=plot_difference,
                plot_tvalues=plot_tvalues,
                plot_significance=plot_significance,
                plot_critical_t=plot_critical_t,
                plot_se=plot_se,
                colors=stats_colors,
                difference_offset=difference_offset,
                significance_color=significance_color,
                linewidth=plot_kwargs[:linewidth],
                legend_labels=plot_kwargs[:legend_labels],
                add_labels=(ch == first(channels_to_plot))
            )
        end
    end

    # Apply axis properties and layout properties
    _apply_axis_properties!.(axes; plot_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...)

    # Override grid/topo titles to include stats info
    if plot_layout.type in (:grid, :topo)
        for (ax_idx, (ax, channel)) in enumerate(zip(axes, layout_channels))
            ax.title = string(channel)
        end
    end

    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)

    # Add origin lines (x=0, y=0)
    if plot_kwargs[:add_xy_origin]
        for ax in axes
            _set_origin_lines!(ax)
        end
    end

    # Add a single legend (on the first axis)
    if !isempty(axes)
        axislegend(axes[1], position=plot_kwargs[:legend_position], framevisible=plot_kwargs[:legend_framevisible])
    end

    # Draw supertitle if user explicitly provided a figure_title
    if !isempty(plot_kwargs[:figure_title])
        Label(fig[0, :], plot_kwargs[:figure_title], fontsize=plot_kwargs[:figure_title_fontsize], font=:bold, tellwidth=false)
    end

    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end

    return (fig=fig, axes=axes)
end


"""
    _plot_erp_stats_channel!(ax, result, channel_idx, channel_sym; kwargs...)

Internal helper to render ERPs, difference wave, t-values, significance bars, 
and critical t-lines for a single channel on a single axis.
"""
function _plot_erp_stats_channel!(
    ax::Axis,
    result::StatsResult,
    channel_idx::Int,
    channel_sym::Symbol;
    condition_selection::Function=conditions(),
    plot_erp::Bool=true,
    plot_difference::Bool=false,
    plot_tvalues::Bool=false,
    plot_significance::Bool=false,
    plot_critical_t::Bool=false,
    plot_se::Bool=false,
    colors::Vector=[:blue, :red, :black, :purple],
    difference_offset::Real=0.0,
    significance_position::Union{Symbol,Real}=:auto,
    significance_color=(:gray, 0.6),
    linewidth::Int=2,
    legend_labels::Vector=[],
    add_labels::Bool=true,
)
    time_points = result.time_points
    t_values = result.stat_matrix.t[channel_idx, :]

    plotted_erps = subset(collect(result.data); condition_selection=condition_selection)

    # Get condition averages for this channel (index by Symbol, not integer position)
    cond_A_erp = result.data[1]
    cond_B_erp = result.data[2]
    !(channel_sym in channel_labels(cond_A_erp)) && error("Channel $channel_sym not found in result data")
    cond_A_avg = cond_A_erp.data[!, channel_sym]
    cond_B_avg = cond_B_erp.data[!, channel_sym]
    erp_time_points = cond_A_erp.data[!, :time]

    # Pre-compute critical t-values if needed
    critical_t_pos = nothing
    critical_t_neg = nothing
    if (plot_significance || plot_critical_t) && plot_tvalues
        if isa(result, AnalyticResult)
            critical_t_pos = fill(result.critical_t, length(time_points))
            critical_t_neg = fill(-result.critical_t, length(time_points))
        else
            critical_t_channel = result.critical_t[channel_idx, :]
            critical_t_pos = critical_t_channel
            critical_t_neg = -critical_t_channel
        end
    end

    # Plot condition averages (ERP waveforms) with optional SE bands
    if plot_erp
        for (i, erp) in enumerate(plotted_erps)
            erp_name = isempty(legend_labels) ? erp.condition_name : (i <= length(legend_labels) ? legend_labels[i] : erp.condition_name)
            erp_avg = erp.data[!, channel_sym]
            lines!(ax, erp_time_points, erp_avg, color=colors[i], linewidth=linewidth, label=add_labels ? erp_name : nothing)

            # SE bands around individual condition ERPs (full display interval)
            if plot_se
                orig_idx = findfirst(x -> x.condition == erp.condition, result.data)
                if !isnothing(orig_idx)
                    se_data = orig_idx == 1 ? result.se_cond1[channel_idx, :] : result.se_cond2[channel_idx, :]
                    band!(ax, erp_time_points, erp_avg .- se_data, erp_avg .+ se_data, color=(colors[i], 0.15))
                end
            end
        end
    end

    # Compute difference wave only when needed
    diff_wave = nothing
    diff_wave_plot = nothing
    if plot_difference || plot_significance
        diff_wave = cond_A_avg .- cond_B_avg
        diff_wave_plot = diff_wave
    end

    # Plot difference wave (if requested)
    if plot_difference
        diff_label = length(legend_labels) >= 3 ? legend_labels[3] : "Difference"
        if difference_offset != 0.0
            diff_wave_plot = diff_wave .+ difference_offset
            diff_label = "Difference (offset=$(difference_offset))"

            hlines!(ax, [difference_offset], color=(:gray, 0.7), linewidth=1)
            text!(
                ax,
                erp_time_points[end] * 0.98,
                difference_offset,
                text="0 μV",
                align=(:right, :center),
                color=:gray,
                fontsize=18,
            )
        else
            diff_label = "Difference"
        end

        lines!(ax, erp_time_points, diff_wave_plot, color=colors[3], linewidth=linewidth, label=add_labels ? diff_label : nothing)

        # SE band around difference wave (full display interval)
        if plot_se
            se_diff_channel = result.se_diff[channel_idx, :]
            band!(
                ax,
                erp_time_points,
                diff_wave_plot .- se_diff_channel,
                diff_wave_plot .+ se_diff_channel,
                color=(colors[3], 0.15),
                label="±1 SEM",
            )
        end
    end

    # Plot t-values (if requested)
    if plot_tvalues
        lines!(ax, time_points, t_values, color=colors[4], linewidth=linewidth, label="t-statistic")
    end

    # Show significance regions as bars
    if plot_significance
        sig_pos = result.masks.positive[channel_idx, :]
        sig_neg = result.masks.negative[channel_idx, :]
        sig_any = sig_pos .| sig_neg

        if any(sig_any)
            sig_regions = _find_continuous_regions(sig_any, time_points)
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
                rect_vertices = [
                    Point2f(t_start, bar_y),
                    Point2f(t_end, bar_y),
                    Point2f(t_end, bar_y + bar_height),
                    Point2f(t_start, bar_y + bar_height),
                ]
                poly!(ax, rect_vertices, color=significance_color, strokewidth=0)
            end
        end
    end

    # Show critical t-values if requested
    if plot_critical_t && plot_tvalues && !isnothing(critical_t_pos)
        lines!(ax, time_points, critical_t_pos, color=:grey, linewidth=linewidth, linestyle=:dashdot, label="Critical t+")
        lines!(ax, time_points, critical_t_neg, color=:grey, linewidth=linewidth, linestyle=:dashdot, label="Critical t-")
    end

    return nothing
end


"""
    _compute_significance_position(position, plot_erp, plot_difference, plot_tvalues,
                               cond_A_avg, cond_B_avg, diff_wave_plot, t_values)

Compute bar_y and bar_height for significance bar rendering. Uses 5% of data range for bar height.
"""
function _compute_significance_position(position, plot_erp, plot_difference, plot_tvalues, cond_A_avg, cond_B_avg, diff_wave_plot, t_values)
    bar_height_frac = 0.05  # 5% of range

    # Determine the relevant data range for bar sizing
    data_range, data_min =
        _get_significance_data_range(plot_erp, plot_difference, plot_tvalues, cond_A_avg, cond_B_avg, diff_wave_plot, t_values)
    bar_height = data_range * bar_height_frac

    # Fallback for zero range
    if bar_height ≈ 0.0
        bar_height = 0.05
    end

    # Determine bar_y based on position
    if position isa Real
        bar_y = Float64(position)
    elseif position == :zero
        bar_y = -bar_height / 2  # Center around zero
    elseif position == :bottom
        bar_y = data_min - data_range * 0.05
    else  # :auto (default)
        if plot_erp || plot_difference
            amp_min = _get_amp_min(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
            amp_max = _get_amp_max(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
            if amp_min <= 0.0 <= amp_max
                bar_y = -bar_height / 2  # Center around zero
            else
                bar_y = amp_min - data_range * 0.05
            end
        else
            bar_y = data_min - data_range * 0.05
        end
    end

    return bar_y, bar_height
end


"""
    _get_significance_data_range(...)

Returns `(range, min)` for the plotted data, used to size and position significance bars.
"""
function _get_significance_data_range(plot_erp, plot_difference, plot_tvalues, cond_A_avg, cond_B_avg, diff_wave_plot, t_values)
    if plot_tvalues && !plot_difference && !plot_erp
        valid = t_values[.!isnan.(t_values)]
        t_min, t_max = extrema(valid)
        return t_max - t_min, t_min
    elseif plot_erp || plot_difference
        amp_min = _get_amp_min(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
        amp_max = _get_amp_max(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
        return amp_max - amp_min, amp_min
    else
        return 1.0, 0.0
    end
end

function _get_amp_min(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
    if plot_erp
        min(minimum(cond_A_avg), minimum(cond_B_avg))
    elseif plot_difference && !isnothing(diff_wave_plot)
        minimum(diff_wave_plot)
    else
        0.0
    end
end

function _get_amp_max(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
    if plot_erp
        max(maximum(cond_A_avg), maximum(cond_B_avg))
    elseif plot_difference && !isnothing(diff_wave_plot)
        maximum(diff_wave_plot)
    else
        0.0
    end
end

"""
    plot_stat_heatmap(result::StatsResult; kwargs...)

Plots a 2D heatmap of the t-statistics (Channels x Time).
"""
function plot_stat_heatmap(result::StatsResult; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)
    plot_kwargs[:window_title] = get(kwargs, :window_title, "T-Statistic Heatmap")
    plot_kwargs[:xlabel] = get(kwargs, :xlabel, "Time (s)")
    plot_kwargs[:ylabel] = get(kwargs, :ylabel, "Channels")

    time_points = result.time_points
    electrodes = result.electrodes
    t_values = result.stat_matrix.t

    fig = Figure(size=(800, 600), title=plot_kwargs[:window_title], figure_padding=plot_kwargs[:figure_padding])

    if length(electrodes) > 40
        yticks = (1:5:length(electrodes), String.(electrodes)[1:5:end])
    else
        yticks = (1:length(electrodes), String.(electrodes))
    end

    ax = Axis(fig[1, 1], yticks=yticks)

    max_t = maximum(abs.(t_values))

    hm = heatmap!(
        ax,
        time_points,
        1:length(electrodes),
        transpose(t_values),
        colormap=get(kwargs, :colormap, :RdBu_r),
        colorrange=get(kwargs, :colorrange, (-max_t, max_t)),
    )

    Colorbar(fig[1, 2], hm, label="t-value")

    _apply_axis_properties!(ax; plot_kwargs...)

    # Draw supertitle if user explicitly provided a figure_title
    if !isempty(plot_kwargs[:figure_title])
        Label(fig[0, :], plot_kwargs[:figure_title], fontsize=plot_kwargs[:figure_title_fontsize], font=:bold, tellwidth=false)
    end

    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end

    return (fig=fig, axes=[ax])
end
