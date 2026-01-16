"""
Plotting functions for MVPA/decoding results.

This module provides functions for visualizing decoding accuracy over time,
confusion matrices, and other decoding analysis results.
"""
# ==============================================================================
#   DEFAULT KEYWORD ARGUMENTS
# ==============================================================================

const PLOT_DECODING_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("Decoding Results", "Title for the plot window"),
    :interactive => (true, "Enable interactive features (true/false)"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("Classification Accuracy", "Label for y-axis"),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Show title (true/false)"),

    # Line styling
    :linewidth => (2, "Line width for decoding curve"),
    :color => (:blue, "Color for decoding curve"),
    :linestyle => (:solid, "Line style"),

    # Chance level
    :chance_color => (:gray, "Color for chance level line"),
    :chance_linestyle => (:dash, "Line style for chance level"),
    :chance_linewidth => (1, "Line width for chance level"),

    # Error bars
    :show_error => (true, "Show standard error bars (true/false)"),
    :error_color => (:blue, "Color for error bars"),
    :error_alpha => (0.3, "Transparency for error shading"),

    # Grid
    :xgrid => (true, "Show x-axis grid (true/false)"),
    :ygrid => (true, "Show y-axis grid (true/false)"),

    # Origin lines
    :add_xy_origin => (true, "Add origin lines at x=0 and y=chance (true/false)"),

    # Figure padding
    :figure_padding => ((10, 10, 10, 10), "Padding around entire figure as (left, right, top, bottom) tuple (in pixels)"),
)

# ==============================================================================
#   HELPER FUNCTIONS
# ==============================================================================

"""
    _add_decoding_origin_lines!(ax::Axis, chance_level::Float64, plot_kwargs::Dict)

Add origin lines for decoding plots: x=0 and y=chance_level.
Always plots chance level line with label.
"""
function _add_decoding_origin_lines!(ax::Axis, chance_level::Float64, plot_kwargs::Dict)
    
    # Add x=0 line if add_xy_origin is true
        vlines!(ax, 0, color = :black, linewidth = 1, linestyle = :dash)
    # Always add y=chance_level line with label
    hlines!(
        ax,
        chance_level,
        color = plot_kwargs[:chance_color],
        linewidth = plot_kwargs[:chance_linewidth],
        linestyle = plot_kwargs[:chance_linestyle],
        label = "Chance ($(round(chance_level, digits = 3)))",
    )
end

"""
    _plot_error_band!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64}, 
                      stderror::Union{Vector{Float64}, Nothing}, plot_kwargs::Dict)

Plot error band (standard error shading) on axis.
"""
function _plot_error_band!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64},Nothing},
    plot_kwargs::Dict,
)
    band!(
        ax,
        times,
        accuracy .- stderror,
        accuracy .+ stderror,
        color = (plot_kwargs[:error_color], plot_kwargs[:error_alpha]),
        label = "±1 SE",
    )
end

"""
    _plot_accuracy_curve!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64}, 
                          plot_kwargs::Dict; label::String = "Accuracy")

Plot main accuracy curve on axis.
"""
function _plot_accuracy_curve!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    plot_kwargs::Dict;
    label::String = "Accuracy",
)
    lines!(
        ax,
        times,
        accuracy,
        color = plot_kwargs[:color],
        linewidth = plot_kwargs[:linewidth],
        linestyle = plot_kwargs[:linestyle],
        label = label,
    )
end

"""
    _setup_axis_limits!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64},
                        stderror::Union{Vector{Float64}, Nothing}, plot_kwargs::Dict)

Setup axis limits for decoding plot.
"""
function _setup_axis_limits!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64}, Nothing},
    plot_kwargs::Dict,
)
    # X-axis limits
    xlim_val = plot_kwargs[:xlim]
    if !isnothing(xlim_val)
        xlims!(ax, xlim_val)
    else
        xlims!(ax, (times[1], times[end]))
    end

    # Y-axis limits
    ylim_val = plot_kwargs[:ylim]
    if !isnothing(ylim_val)
        ylims!(ax, ylim_val)
    else
        # Auto-determine y limits
        y_min = minimum(accuracy)
        y_max = maximum(accuracy)
        if plot_kwargs[:show_error] && !isnothing(stderror)
            y_min = min(y_min, minimum(accuracy .- stderror))
            y_max = max(y_max, maximum(accuracy .+ stderror))
        end
        y_range = y_max - y_min
        ylims!(ax, (y_min - 0.05 * y_range, y_max + 0.05 * y_range))
    end
end

"""
    _plot_decoding_to_axis!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64},
                            stderror::Union{Vector{Float64}, Nothing}, chance_level::Float64,
                            plot_kwargs::Dict; show_legend::Bool = true, curve_label::String = "Accuracy")

Base function that plots decoding data to an existing axis.

This function handles the core plotting logic:
- Setting up axis limits
- Adding origin lines (x=0 and y=chance_level)
- Plotting error band
- Plotting accuracy curve
- Optionally showing legend

# Arguments
- `ax::Axis`: The axis to plot to
- `times::Vector{Float64}`: Time points
- `accuracy::Vector{Float64}`: Accuracy values
- `stderror::Union{Vector{Float64}, Nothing}`: Standard error values (optional)
- `chance_level::Float64`: Chance level for reference line
- `plot_kwargs::Dict`: Plotting keyword arguments
- `show_legend::Bool`: Whether to show legend (default: true)
- `curve_label::String`: Label for accuracy curve (default: "Accuracy")
"""
function _plot_decoding_to_axis!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64}, Nothing},
    chance_level::Float64,
    plot_kwargs::Dict;
    show_legend::Bool = true,
    curve_label::String = "Accuracy",
)
    # Setup axis limits
    _setup_axis_limits!(ax, times, accuracy, stderror, plot_kwargs)

    # Add origin lines (x=0 and y=chance_level)
    if plot_kwargs[:add_xy_origin]
        _add_decoding_origin_lines!(ax, chance_level, plot_kwargs)
    end

    # Plot error band
    if plot_kwargs[:show_error] && !isnothing(stderror)
        _plot_error_band!(ax, times, accuracy, stderror, plot_kwargs)
    end

    # Plot main accuracy curve
    _plot_accuracy_curve!(ax, times, accuracy, plot_kwargs; label = curve_label)

    # Show legend if requested
    if show_legend
        axislegend(ax, position = :rt)
    end
end

# ==============================================================================
#   MAIN PLOTTING FUNCTIONS
# ==============================================================================

"""
    plot_decoding(decoded::DecodedData; kwargs...)

Plot decoding accuracy over time.

Creates a line plot showing classification accuracy at each time point,
with optional error bars and chance level reference line.

# Arguments
- `decoded::DecodedData`: DecodedData object containing decoding results
- `kwargs`: Additional keyword arguments (see PLOT_DECODING_KWARGS)

# Examples
```julia
# Basic plot
plot_decoding(decoded)

# Custom styling
plot_decoding(decoded, color=:red, linewidth=3, show_error=false)

# With custom title
plot_decoding(decoded, title="Face vs. Object Decoding")
```
"""
function plot_decoding(decoded::DecodedData; kwargs...)
    # Merge defaults with user kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)

    # Extract parameters
    display_plot = plot_kwargs[:display_plot]
    figure_title = plot_kwargs[:figure_title]
    times = decoded.times
    accuracy = decoded.average_score
    stderror = decoded.stderror
    chance_level = decoded.parameters.chance_level
    title_text = plot_kwargs[:title]
    show_title = plot_kwargs[:show_title]

    # Create figure
    fig = Figure(title = figure_title, size = (800, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
        xgridvisible = plot_kwargs[:xgrid],
        ygridvisible = plot_kwargs[:ygrid],
    )

    # Set title if requested
    if show_title && !isempty(title_text)
        ax.title = title_text
    elseif show_title && isempty(title_text)
        # Default title from decoded data
        ax.title = "Decoding: $(join(decoded.condition_names, " vs "))"
    end

    # Plot decoding data to axis
    _plot_decoding_to_axis!(ax, times, accuracy, stderror, chance_level, plot_kwargs)

    # Display if requested
    if display_plot
        display(fig)
    end

    return fig, ax
end

"""
    plot_decoding(decoded_list::Vector{DecodedData}; kwargs...)

Plot multiple decoding results in separate subplots (one per subject).

Creates a grid of subplots using `best_rect` to determine optimal layout,
with each subject's decoding results in its own subplot.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects, one per subject
- `kwargs`: Additional keyword arguments (see PLOT_DECODING_KWARGS)

# Examples
```julia
# Plot all subjects in separate subplots
plot_decoding(all_decoded, title = "Individual Subjects")
```
"""
function plot_decoding(decoded_list::Vector{DecodedData}; kwargs...)
    if isempty(decoded_list)
        @minimal_error_throw("Cannot plot empty decoded data list")
    end

    # Merge defaults with user kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)

    # Extract parameters
    display_plot = plot_kwargs[:display_plot]
    figure_title = plot_kwargs[:figure_title]
    title_text = plot_kwargs[:title]
    show_title = plot_kwargs[:show_title]
    
    # Determine optimal subplot layout using best_rect
    n_subjects = length(decoded_list)
    rows, cols = best_rect(n_subjects)
    
    # Create figure with appropriate size for subplots
    fig = Figure(title = figure_title, size = (400 * cols, 300 * rows))
    
    # Determine common time range and y limits across all subjects
    all_times = [d.times for d in decoded_list]
    time_min = minimum([t[1] for t in all_times])
    time_max = maximum([t[end] for t in all_times])
    
    all_accuracies = [d.average_score for d in decoded_list]
    y_min = minimum([minimum(acc) for acc in all_accuracies])
    y_max = maximum([maximum(acc) for acc in all_accuracies])
    
    # Account for error bars
    for d in decoded_list
        if !isnothing(d.stderror)
            y_min = min(y_min, minimum(d.average_score .- d.stderror))
            y_max = max(y_max, maximum(d.average_score .+ d.stderror))
        end
    end
    
    y_range = y_max - y_min
    y_lims = (y_min - 0.05 * y_range, y_max + 0.05 * y_range)
    
    # Get chance level (use first decoded's chance level)
    chance_level = decoded_list[1].parameters.chance_level
    
    # Create subplot for each subject
    for (idx, decoded) in enumerate(decoded_list)
        # Calculate row and column for this subplot
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        
        # Create axis for this subject
        ax = Axis(
            fig[row, col],
            xlabel = (row == rows) ? plot_kwargs[:xlabel] : "",  # Only show xlabel on bottom row
            ylabel = (col == 1) ? plot_kwargs[:ylabel] : "",    # Only show ylabel on left column
            xgridvisible = plot_kwargs[:xgrid],
            ygridvisible = plot_kwargs[:ygrid],
        )

        # Set title for this subplot
        if show_title
            if !isempty(title_text)
                ax.title = title_text
            else
                # Use subject identifier from file name or index
                subject_id = isnothing(decoded.file) ? "Subject $idx" : basename(decoded.file)
                ax.title = subject_id
            end
        end

        # Set axis limits
        xlims!(ax, (time_min, time_max))
        ylims!(ax, y_lims)

        # Get data for this subject
        times = decoded.times
        accuracy = decoded.average_score
        stderror = decoded.stderror

        # Plot decoding data to axis (no legend for individual subplots)
        _plot_decoding_to_axis!(ax, times, accuracy, stderror, chance_level, plot_kwargs; show_legend = false, curve_label = "")
    end
    
    # Display if requested
    if display_plot
        display(fig)
    end
    
    return fig
end

"""
    plot_decoding(decoded::DecodedData, stats::DecodingStatisticsResult; kwargs...)

Plot decoding accuracy over time with significance markers from statistical test results.

# Arguments
- `decoded::DecodedData`: DecodedData object containing decoding results
- `stats::DecodingStatisticsResult`: DecodingStatisticsResult from `test_against_chance` or `test_against_chance_cluster`
- `kwargs`: Additional keyword arguments (see PLOT_DECODING_KWARGS)
  - `show_significance::Bool`: Show significance markers (default: true)
  - `sig_color::ColorType`: Color for significance markers (default: :yellow)
  - `sig_alpha::Float64`: Transparency for significance markers (default: 0.3)
  - `sig_bar_position::Union{Float64, Symbol}`: Y-position for significance bars - `:bottom` (default), `:top`, or a Float64 value

# Examples
```julia
# Test and plot with significance
stats = test_against_chance(decoded_list, alpha=0.05, correction_method=:bonferroni)
plot_decoding(grand_avg, stats, show_significance=true)
```
"""
function plot_decoding(decoded::DecodedData, stats::DecodingStatisticsResult; kwargs...)
    # Merge defaults with user kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)
    
    # Extract parameters
    show_significance = get(kwargs, :show_significance, true)
    sig_color = get(kwargs, :sig_color, :yellow)
    sig_alpha = get(kwargs, :sig_alpha, 0.3)
    sig_bar_position = get(kwargs, :sig_bar_position, :bottom)
    
    # Validate that times match
    if decoded.times != stats.times
        @minimal_error_throw("DecodedData and DecodingStatisticsResult must have matching time vectors")
    end
    
    # Create base plot (without displaying yet)
    display_plot_orig = plot_kwargs[:display_plot]
    fig, ax = plot_decoding(decoded; kwargs..., display_plot=false, show_significance=false)
    
    # Compute y-limits from data for significance bar positioning
    accuracy = decoded.average_score
    stderror = decoded.stderror
    y_min = minimum(accuracy)
    y_max = maximum(accuracy)
    if !isnothing(stderror)
        y_min = min(y_min, minimum(accuracy .- stderror))
        y_max = max(y_max, maximum(accuracy .+ stderror))
    end
    y_range = y_max - y_min
    y_min -= 0.05 * y_range
    y_max += 0.05 * y_range
    y_range = y_max - y_min
    
    # Add significance markers if requested
    if show_significance && any(stats.significant_mask)
        # Find continuous significant regions
        sig_regions = _find_continuous_regions(stats.significant_mask, stats.times)
        
        # Determine bar position
        if sig_bar_position == :bottom
            bar_y = y_min + 0.02 * y_range
        elseif sig_bar_position == :top
            bar_y = y_max - 0.02 * y_range
        elseif isa(sig_bar_position, Number)
            bar_y = Float64(sig_bar_position)
        else
            bar_y = y_min + 0.02 * y_range
        end
        
        bar_height = y_range * 0.02  # 2% of y-range
        
        # Plot significance bars for each continuous region
        for region in sig_regions
            t_start, t_end = region
            # Use poly! to create a horizontal bar (rectangle)
            poly!(
                ax,
                Rect(t_start, bar_y, t_end - t_start, bar_height),
                color = (sig_color, sig_alpha),
                strokewidth = 0,
            )
        end
        
        # If clusters are available, add cluster labels
        if !isnothing(stats.clusters) && !isempty(stats.clusters)
            for cluster in stats.clusters
                if cluster.is_significant
                    # Add text label at cluster center
                    t_center = (cluster.time_range[1] + cluster.time_range[2]) / 2
                    text!(
                        ax,
                        t_center,
                        bar_y + bar_height * 1.5,
                        text = "p=$(round(cluster.p_value, digits=3))",
                        align = (:center, :bottom),
                        fontsize = 10,
                        color = sig_color,
                    )
                end
            end
        end
    end
    
    # Display if originally requested
    if display_plot_orig
        display(fig)
    end
    
    return fig
end

"""
    _find_continuous_regions(mask::BitVector, times::Vector{Float64})

Find continuous regions of true values in a boolean mask.

# Returns
- `regions::Vector{Tuple{Float64, Float64}}`: List of (start, end) time ranges
"""
function _find_continuous_regions(mask::BitVector, times::Vector{Float64})
    regions = Tuple{Float64, Float64}[]
    
    if !any(mask)
        return regions
    end
    
    in_region = false
    region_start = 0
    
    for (idx, is_sig) in enumerate(mask)
        if is_sig && !in_region # Start new region
            in_region = true
            region_start = idx
        elseif !is_sig && in_region # End current region
            in_region = false
            region_end = idx - 1
            push!(regions, (times[region_start], times[region_end]))
        end
    end
    
    # Handle region that extends to end
    if in_region
        push!(regions, (times[region_start], times[end]))
    end
    
    return regions
end

"""
    plot_confusion_matrix(decoded::DecodedData; time_point::Union{Float64, Int, Nothing} = nothing, kwargs...)

Plot confusion matrix for decoding results.

# Arguments
- `decoded::DecodedData`: DecodedData object containing confusion matrices
- `time_point::Union{Float64, Int, Nothing}`: Time point to plot (in seconds or index). If nothing, plots average across all time points
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot confusion matrix at a specific time point
plot_confusion_matrix(decoded, time_point=0.3)

# Plot average confusion matrix
plot_confusion_matrix(decoded)
```
"""
function plot_confusion_matrix(decoded::DecodedData; time_point::Union{Float64, Int, Nothing} = nothing, kwargs...)
    if isnothing(decoded.confusion_matrix)
        @minimal_error_throw("No confusion matrix data available in DecodedData")
    end

    # Determine which time point to plot
    if isnothing(time_point)
        # Average across all time points
        confusion = mean(decoded.confusion_matrix, dims = 1)[1, :, :]
        title_text = "Average Confusion Matrix"
    elseif time_point isa Float64
        # Find closest time point
        time_idx = argmin(abs.(decoded.times .- time_point))
        confusion = decoded.confusion_matrix[time_idx, :, :]
        title_text = "Confusion Matrix at $(round(time_point, digits=3)) s"
    else
        # Use as index
        confusion = decoded.confusion_matrix[time_point, :, :]
        title_text = "Confusion Matrix at $(round(decoded.times[time_point], digits=3)) s"
    end

    # Create figure
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1], title = title_text, xlabel = "Predicted", ylabel = "True")

    # Create heatmap
    heatmap!(
        ax,
        confusion,
        colormap = :viridis,
        colorrange = (0, 1),
    )

    # Add text labels
    n_classes = size(confusion, 1)
    for i in 1:n_classes
        for j in 1:n_classes
            text!(
                ax,
                j,
                i,
                text = string(round(confusion[i, j], digits = 2)),
                color = confusion[i, j] > 0.5 ? :white : :black,
                align = (:center, :center),
            )
        end
    end

    # Set ticks
    ax.xticks = (1:n_classes, decoded.condition_names)
    ax.yticks = (1:n_classes, decoded.condition_names)

    display(fig)
    return fig
end

