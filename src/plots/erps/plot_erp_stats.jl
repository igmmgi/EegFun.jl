"""
Plotting functions for statistical test results (analytic and permutation tests).
"""


"""
    plot_erp_stats(result::StatsResult;
                    layout::Union{Symbol, PlotLayout} = :single,
                    channel_selection::Function = channels(),
                    plot_erp::Bool = true,
                    plot_difference::Bool = false,
                    plot_tvalues::Bool = false,
                    plot_significance::Bool = false,
                    plot_critical_t::Bool = false,
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
- `plot_erp::Bool`: Whether to plot ERP waveforms (condition averages) (default: true)
- `plot_difference::Bool`: Whether to plot difference wave (A-B) (default: false)
- `plot_tvalues::Bool`: Whether to plot t-statistics (default: false)
- `plot_significance::Bool`: Whether to highlight significant time points (default: false)
- `plot_critical_t::Bool`: Whether to plot critical t-values (default: false). Only relevant when `plot_tvalues=true`
- `difference_offset::Real`: Vertical offset for difference wave (default: 0.0). Set to non-zero to shift for visibility
- `significance_position::Union{Symbol, Real}`: Position for significance bars (default: `:auto`). Options:
  - `:auto` - Automatically place at y=0 if visible, otherwise at bottom (default)
  - `:zero` - Always place at y=0
  - `:bottom` - Always place at bottom of plot/spine
  - `Real` - Custom y-position (e.g., `-5.0` to place at -5 μV)
- `significance_color`: Color for significance bars (default: `(:gray, 0.6)`). Can be any Makie color specification

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
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    plot_erp::Bool = true,
    plot_difference::Bool = false,
    plot_tvalues::Bool = false,
    plot_significance::Bool = false,
    plot_critical_t::Bool = false,
    difference_offset::Real = 0.0,
    significance_position::Union{Symbol,Real} = :auto,
    significance_color = (:gray, 0.6),
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
    plot_kwargs[:figure_title] = get(kwargs, :figure_title, "ERP Stats")

    # Determine selected channels from result electrodes
    all_electrodes = result.electrodes
    selected_mask = channel_selection(all_electrodes)
    selected_channels = all_electrodes[selected_mask]

    if isempty(selected_channels)
        error("No channels matched the selection. Available channels: $(all_electrodes)")
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
    if layout == :single && plot_kwargs[:show_title] && plot_kwargs[:title] == ""
        if length(selected_channels) == 1
            plot_kwargs[:title] = "$(selected_channels[1]) - $title_suffix"
        else
            plot_kwargs[:title] = "$(_print_vector(selected_channels)) - $title_suffix"
        end
    end

    # Extract layout_* parameters for layout system
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    # Get the EEG layout from the result data (for topo positioning)
    eeg_layout = result.data[1].layout

    # Create figure and apply layout system
    fig = Figure(size = (800, 600), title = plot_kwargs[:figure_title], figure_padding = plot_kwargs[:figure_padding])
    plot_layout = create_layout(layout, selected_channels, eeg_layout; layout_kwargs...)
    axes, layout_channels = _apply_layout!(fig, plot_layout; plot_kwargs...)

    # Render stats content on each axis
    for (ax_idx, (ax, channel)) in enumerate(zip(axes, layout_channels))
        channels_to_plot = plot_layout.type == :single ? selected_channels : [channel]

        for ch in channels_to_plot
            channel_idx = findfirst(==(ch), all_electrodes)
            channel_idx === nothing && continue

            _plot_erp_stats_channel!(
                ax,
                result,
                channel_idx,
                ch;
                plot_erp = plot_erp,
                plot_difference = plot_difference,
                plot_tvalues = plot_tvalues,
                plot_significance = plot_significance,
                plot_critical_t = plot_critical_t,
                difference_offset = difference_offset,
                significance_position = significance_position,
                significance_color = significance_color,
                linewidth = plot_kwargs[:linewidth],
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
        axislegend(axes[1], position = plot_kwargs[:legend_position])
    end

    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end

    return (fig = fig, axes = axes)
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
    plot_erp::Bool = true,
    plot_difference::Bool = false,
    plot_tvalues::Bool = false,
    plot_significance::Bool = false,
    plot_critical_t::Bool = false,
    difference_offset::Real = 0.0,
    significance_position::Union{Symbol,Real} = :auto,
    significance_color = (:gray, 0.6),
    linewidth::Int = 2,
)
    time_points = result.time_points
    t_values = result.stat_matrix.t[channel_idx, :]

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

    # Plot condition averages (ERP waveforms)
    if plot_erp
        cond_A_name = result.data[1].condition_name
        cond_B_name = result.data[2].condition_name
        lines!(ax, erp_time_points, cond_A_avg, color = :blue, linewidth = linewidth, label = cond_A_name)
        lines!(ax, erp_time_points, cond_B_avg, color = :red, linewidth = linewidth, label = cond_B_name)
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
        if difference_offset != 0.0
            diff_wave_plot = diff_wave .+ difference_offset
            diff_label = "Difference (offset=$(difference_offset))"

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
        else
            diff_label = "Difference"
        end

        lines!(ax, erp_time_points, diff_wave_plot, color = :black, linewidth = linewidth, label = diff_label)
    end

    # Plot t-values (if requested)
    if plot_tvalues
        lines!(ax, time_points, t_values, color = :purple, linewidth = linewidth, label = "t-statistic")
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
                poly!(ax, rect_vertices, color = significance_color, strokewidth = 0)
            end
        end
    end

    # Show critical t-values if requested
    if plot_critical_t && plot_tvalues && critical_t_pos !== nothing
        lines!(ax, time_points, critical_t_pos, color = :grey, linewidth = linewidth, linestyle = :dashdot, label = "Critical t+")
        lines!(ax, time_points, critical_t_neg, color = :grey, linewidth = linewidth, linestyle = :dashdot, label = "Critical t-")
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
    elseif plot_difference && diff_wave_plot !== nothing
        minimum(diff_wave_plot)
    else
        0.0
    end
end

function _get_amp_max(plot_erp, plot_difference, cond_A_avg, cond_B_avg, diff_wave_plot)
    if plot_erp
        max(maximum(cond_A_avg), maximum(cond_B_avg))
    elseif plot_difference && diff_wave_plot !== nothing
        maximum(diff_wave_plot)
    else
        0.0
    end
end
