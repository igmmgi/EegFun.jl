"""
Simple plotting functions for statistical test results.

This module provides visualization functions for permutation test and analytic t-test results.
"""

"""
    find_continuous_regions(mask::BitVector, time_points::Vector{Float64})

Find continuous regions where mask is true, returning start and end times.
"""
function find_continuous_regions(mask::BitVector, time_points::Vector{Float64})
    regions = Vector{Tuple{Float64,Float64}}()

    if isempty(mask) || !any(mask)
        return regions
    end

    in_region = false
    start_idx = 0

    for (i, is_sig) in enumerate(mask)
        if is_sig && !in_region
            # Start of a new region
            in_region = true
            start_idx = i
        elseif !is_sig && in_region
            # End of a region
            in_region = false
            push!(regions, (time_points[start_idx], time_points[i-1]))
        end
    end

    # Handle case where region extends to end
    if in_region
        push!(regions, (time_points[start_idx], time_points[end]))
    end

    return regions
end

"""
    plot_analytic_ttest(result::StatisticalTestResult;
                        channel::Symbol,
                        plot_erp::Bool = true,
                        plot_difference::Bool = false,
                        plot_tvalues::Bool = false,
                        show_significance::Bool = false,
                        show_critical_t::Bool = false,
                        shift_difference::Bool = false)

Plot ERP waveforms and statistical results for analytic t-test or cluster permutation test with flexible component control.

Works with both `AnalyticTTestResult` (from `analytic_ttest`) and `ClusterPermutationResult` (from `cluster_permutation_test`).

# Arguments
- `result::StatisticalTestResult`: Results from `analytic_ttest` or `cluster_permutation_test`
- `channel::Symbol`: Channel/electrode to plot
- `plot_erp::Bool`: Whether to plot ERP waveforms (condition averages) (default: true)
- `plot_difference::Bool`: Whether to plot difference wave (A-B) (default: false)
- `plot_tvalues::Bool`: Whether to plot t-statistics (default: false)
- `show_significance::Bool`: Whether to highlight significant time points (default: false)
- `show_critical_t::Bool`: Whether to plot critical t-values (default: false). Only relevant when `plot_tvalues=true`
- `shift_difference::Bool`: Whether to shift difference wave up for visibility (default: false). Only relevant when `plot_difference=true`
- `sig_bar_position::Union{Symbol, Float64}`: Position for significance bars (default: `:auto`). Options:
  - `:auto` - Automatically place at y=0 if visible, otherwise at bottom (default)
  - `:zero` - Always place at y=0
  - `:bottom` - Always place at bottom of plot/spine
  - `Float64` - Custom y-position (e.g., `-5.0` to place at -5 μV)
- `sig_bar_color`: Color for significance bars (default: `(:gray, 0.6)`). Can be any Makie color specification (e.g., `:red`, `(:blue, 0.5)`, `RGB(1, 0, 0)`)

# Returns
- `fig::Figure`: Makie figure

# Examples
```julia
# With analytic t-test results
result_analytic = analytic_ttest(prepared, correction_method=:no)
fig = plot_analytic_ttest(result_analytic, channel=:PO7, 
                         plot_erp=true, plot_difference=true, show_significance=true)

# With cluster permutation test results
result_cluster = cluster_permutation_test(prepared, n_permutations=1000)
fig = plot_analytic_ttest(result_cluster, channel=:PO7,
                         plot_erp=true, plot_difference=true, show_significance=true, show_critical_t=true)
```
"""
function plot_analytic_ttest(
    result::StatisticalTestResult;
    channel::Symbol,
    plot_erp::Bool = true,
    plot_difference::Bool = false,
    plot_tvalues::Bool = false,
    show_significance::Bool = false,
    show_critical_t::Bool = false,
    shift_difference::Bool = false,
    sig_bar_position::Union{Symbol,Float64} = :auto,
    sig_bar_color = (:gray, 0.6),
)
    # If show_critical_t is requested, automatically enable t-value plotting
    # (since critical t-values only make sense when plotting t-values)
    plot_tvalues = plot_tvalues || show_critical_t

    # Validate that at least something is being plotted
    if !plot_erp && !plot_difference && !plot_tvalues
        error("At least one of plot_erp, plot_difference, or plot_tvalues must be true")
    end

    # Find channel index
    channel_idx = findfirst(==(channel), result.electrodes)
    if channel_idx === nothing
        error("Channel $channel not found in results. Available channels: $(result.electrodes)")
    end

    # Get data for this channel
    time_points = result.time_points
    t_values = result.stat_matrix.t[channel_idx, :]

    # Get p-values if available (only for AnalyticTTestResult)
    p_values = if isa(result, AnalyticTTestResult)
        result.stat_matrix.p[channel_idx, :]
    else
        # For ClusterPermutationResult, we don't have p_matrix, but we don't need it for plotting
        similar(t_values, Float64)  # Dummy array, won't be used
    end

    # Get condition averages for this channel (both result types now have data)
    cond_A_erp = result.data[1]
    cond_B_erp = result.data[2]
    channel_col = findfirst(==(channel), channel_labels(cond_A_erp))
    channel_col === nothing && error("Channel $channel not found in result data")
    cond_A_avg = cond_A_erp.data[!, channel_col]
    cond_B_avg = cond_B_erp.data[!, channel_col]
    erp_time_points = cond_A_erp.data[!, :time]

    # Calculate difference: A - B
    # When A = B, difference = 0
    # When A > B, difference > 0 (positive)
    # When A < B, difference < 0 (negative)
    diff_wave = cond_A_avg .- cond_B_avg

    # Pre-compute critical t-values if needed for significance bars or critical t lines
    critical_t_pos = nothing
    critical_t_neg = nothing
    if (show_significance || show_critical_t) && plot_tvalues
        if isa(result, AnalyticTTestResult)
            # Use pre-computed critical t-value (uniform across all points)
            critical_t_pos = fill(result.critical_t, length(time_points))
            critical_t_neg = fill(-result.critical_t, length(time_points))
        else
            # For ClusterPermutationResult, use the pre-computed critical_t (varies by electrode/time)
            critical_t_channel = result.critical_t[channel_idx, :]
            critical_t_pos = critical_t_channel
            critical_t_neg = -critical_t_channel
        end
    end

    # Create figure
    fig = Figure(size = (800, 600))

    # Always use single axis (no dual axes)
    # Get title suffix based on result type
    title_suffix = if isa(result, AnalyticTTestResult)
        "Analytic t-test ($(result.test_info.correction_method))"
    else
        "Cluster permutation test ($(result.test_info.cluster_info.threshold_method))"
    end

    # Determine y-axis label
    if plot_tvalues && (plot_erp || plot_difference)
        ylabel_str = "Amplitude (μV) / t-statistic"
    elseif plot_tvalues
        ylabel_str = "t-statistic"
    else
        ylabel_str = "Amplitude (μV)"
    end

    ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = ylabel_str, title = "$channel - $title_suffix")

    # Plot condition averages (ERP waveforms)
    if plot_erp
        cond_A_name = result.data[1].condition_name
        cond_B_name = result.data[2].condition_name
        lines!(ax, erp_time_points, cond_A_avg, color = :blue, linewidth = 2, label = cond_A_name)
        lines!(ax, erp_time_points, cond_B_avg, color = :red, linewidth = 2, label = cond_B_name)
    end

    # Initialize diff_offset for use in significance markers
    diff_offset = 0.0
    diff_wave_plot = diff_wave

    # Plot difference wave (if requested)
    if plot_difference
        if shift_difference && plot_erp
            # Shift difference wave up for visibility (similar to MATLAB)
            # Calculate offset to place difference wave above the condition averages
            max_cond = maximum([maximum(cond_A_avg), maximum(cond_B_avg)])
            min_cond = minimum([minimum(cond_A_avg), minimum(cond_B_avg)])
            range_cond = max_cond - min_cond
            diff_offset = max_cond + range_cond * 0.3  # Offset by 30% of condition range

            # The difference wave: when conditions overlap (A = B), diff = 0, so shifted wave = offset
            diff_wave_plot = diff_wave .+ diff_offset
            diff_label = "Difference (A-B, shifted)"

            # Horizontal line at offset level = "zero difference" reference
            hlines!(ax, [diff_offset], color = (:gray, 0.7), linewidth = 1, linestyle = :dot)

            # Add text annotation to clarify: this line represents zero difference
            text!(
                ax,
                erp_time_points[end] * 0.98,
                diff_offset,
                text = "0 μV (A=B)",
                align = (:right, :center),
                color = :gray,
                fontsize = 9,
            )
        else
            # Plot difference wave at actual zero (not shifted)
            diff_wave_plot = diff_wave
            diff_offset = 0.0
            diff_label = "Difference (A-B)"
        end

        lines!(
            ax,
            erp_time_points,
            diff_wave_plot,
            color = :black,
            linewidth = 2,
            linestyle = :dash,
            label = diff_label,
        )

        # Add zero lines
        vlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)
        hlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)
    end

    # Plot t-values (if requested)
    if plot_tvalues
        # Plot t-value line on same axis
        lines!(ax, time_points, t_values, color = :purple, linewidth = 2, linestyle = :solid, label = "t-statistic")

        # Add zero line for t-values (only if not already added for difference wave)
        if !plot_difference
            hlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)
            vlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)
        end
    end

    # Show significance regions as grey bars at y=0 (or bottom spine)
    if show_significance
        # Use the result's significance masks, which already have the correction applied
        # For AnalyticTTestResult, these masks reflect the correction method (no correction or Bonferroni)
        # For ClusterPermutationResult, these masks reflect cluster-level significance
        sig_pos = result.masks.positive[channel_idx, :]
        sig_neg = result.masks.negative[channel_idx, :]
        sig_any = sig_pos .| sig_neg  # Any significance (positive or negative)

        if any(sig_any)
            # Find continuous significant regions
            sig_regions = find_continuous_regions(sig_any, time_points)

            # Always use the single axis for significance bars
            sig_ax = ax

            # Calculate bar position and height based on data range and user preference
            if sig_bar_position isa Float64
                # User specified custom position
                bar_y = sig_bar_position
                # Calculate height based on data range
                if plot_tvalues && !plot_difference && !plot_erp
                    t_range = maximum(t_values[.!isnan.(t_values)]) - minimum(t_values[.!isnan.(t_values)])
                    bar_height = t_range * 0.02  # 2% of range
                elseif plot_difference || plot_erp
                    if plot_erp
                        amp_range =
                            maximum([maximum(cond_A_avg), maximum(cond_B_avg)]) -
                            minimum([minimum(cond_A_avg), minimum(cond_B_avg)])
                    else
                        amp_range = maximum(diff_wave_plot) - minimum(diff_wave_plot)
                    end
                    bar_height = amp_range * 0.02  # 2% of range
                else
                    bar_height = 0.05
                end
            elseif sig_bar_position == :zero
                # Always place at y=0
                bar_y = 0.0
                if plot_tvalues && !plot_difference && !plot_erp
                    t_range = maximum(t_values[.!isnan.(t_values)]) - minimum(t_values[.!isnan.(t_values)])
                    bar_height = t_range * 0.02
                elseif plot_difference || plot_erp
                    if plot_erp
                        amp_range =
                            maximum([maximum(cond_A_avg), maximum(cond_B_avg)]) -
                            minimum([minimum(cond_A_avg), minimum(cond_B_avg)])
                    else
                        amp_range = maximum(diff_wave_plot) - minimum(diff_wave_plot)
                    end
                    bar_height = amp_range * 0.02
                else
                    bar_height = 0.05
                end
            elseif sig_bar_position == :bottom
                # Always place at bottom of plot
                if plot_tvalues && !plot_difference && !plot_erp
                    t_min = minimum(t_values[.!isnan.(t_values)])
                    t_max = maximum(t_values[.!isnan.(t_values)])
                    t_range = t_max - t_min
                    bar_y = t_min - t_range * 0.05  # 5% below minimum
                    bar_height = t_range * 0.02
                elseif plot_difference || plot_erp
                    if plot_erp
                        amp_min = minimum([minimum(cond_A_avg), minimum(cond_B_avg)])
                        amp_max = maximum([maximum(cond_A_avg), maximum(cond_B_avg)])
                    else
                        amp_min = minimum(diff_wave_plot)
                        amp_max = maximum(diff_wave_plot)
                    end
                    amp_range = amp_max - amp_min
                    bar_y = amp_min - amp_range * 0.05  # 5% below minimum
                    bar_height = amp_range * 0.02
                else
                    bar_y = -0.1
                    bar_height = 0.05
                end
            else  # :auto (default)
                # Automatic positioning: y=0 if visible, otherwise at bottom
                if plot_tvalues && !plot_difference && !plot_erp
                    # For t-value only plot, place bars at bottom of t-value range
                    t_min = minimum(t_values[.!isnan.(t_values)])
                    t_max = maximum(t_values[.!isnan.(t_values)])
                    t_range = t_max - t_min
                    bar_y = t_min - t_range * 0.05  # 5% below minimum
                    bar_height = t_range * 0.02  # 2% of range
                elseif plot_difference || plot_erp
                    # For amplitude plot, place bars at y=0 or slightly below
                    if plot_erp
                        amp_min = minimum([minimum(cond_A_avg), minimum(cond_B_avg)])
                        amp_max = maximum([maximum(cond_A_avg), maximum(cond_B_avg)])
                    else
                        amp_min = minimum(diff_wave_plot)
                        amp_max = maximum(diff_wave_plot)
                    end
                    amp_range = amp_max - amp_min
                    # Place bars at y=0 if it's visible, otherwise at bottom
                    if amp_min <= 0.0 <= amp_max
                        bar_y = 0.0
                        bar_height = amp_range * 0.02  # 2% of range
                    else
                        bar_y = amp_min - amp_range * 0.05  # 5% below minimum
                        bar_height = amp_range * 0.02  # 2% of range
                    end
                else
                    # Fallback: place at y=0
                    bar_y = 0.0
                    bar_height = 0.05
                end
            end

            # Plot grey bars for each significant region
            for (t_start, t_end) in sig_regions
                # Create rectangle vertices: bottom-left, bottom-right, top-right, top-left
                rect_vertices = [
                    Point2f(t_start, bar_y),
                    Point2f(t_end, bar_y),
                    Point2f(t_end, bar_y + bar_height),
                    Point2f(t_start, bar_y + bar_height),
                ]
                # Use poly! to draw filled rectangle
                poly!(sig_ax, rect_vertices, color = sig_bar_color, strokewidth = 0)
            end
        end
    end

    # Show critical t-values if requested (only relevant for t-value plots)
    if show_critical_t && plot_tvalues && critical_t_pos !== nothing
        # Plot critical t boundaries on same axis (symmetric around zero)
        # These are the actual critical t-values, not scaled to amplitude
        lines!(
            ax,
            time_points,
            critical_t_pos,
            color = :grey,
            linewidth = 2,
            linestyle = :dashdot,
            label = "Critical t+",
        )
        lines!(
            ax,
            time_points,
            critical_t_neg,
            color = :grey,
            linewidth = 2,
            linestyle = :dashdot,
            label = "Critical t-",
        )
    end

    # Add legend
    axislegend(ax, position = :rt)

    return fig
end
