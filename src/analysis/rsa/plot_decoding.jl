"""
Plotting functions for MVPA/decoding results.

This module provides functions for visualizing decoding accuracy over time,
confusion matrices, and other decoding analysis results.
"""

using Makie
using CairoMakie

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
    :show_chance => (true, "Show chance level line (true/false)"),
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
    plot_kwargs = merge(PLOT_DECODING_KWARGS, Dict(kwargs))

    # Extract parameters (unwrap tuples if needed)
    _get_val(k) = isa(plot_kwargs[k], Tuple) ? plot_kwargs[k][1] : plot_kwargs[k]
    display_plot = _get_val(:display_plot)
    figure_title = _get_val(:figure_title)
    times = decoded.times
    accuracy = decoded.average_score
    stderror = decoded.stderror
    chance_level = decoded.chance_level
    show_chance = _get_val(:show_chance)
    show_error = _get_val(:show_error) && !isnothing(stderror)
    title_text = _get_val(:title)
    show_title = _get_val(:show_title)

    # Create figure
    fig = Figure(title = figure_title, resolution = (800, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = _get_val(:xlabel),
        ylabel = _get_val(:ylabel),
        xgridvisible = _get_val(:xgrid),
        ygridvisible = _get_val(:ygrid),
    )

    # Set title if requested
    if show_title && !isempty(title_text)
        ax.title = title_text
    elseif show_title && isempty(title_text)
        # Default title from decoded data
        method_str = string(decoded.method) |> uppercase
        ax.title = "$(method_str) Decoding: $(join(decoded.condition_names, " vs "))"
    end

    # Set axis limits
    xlim_val = _get_val(:xlim)
    if !isnothing(xlim_val)
        xlims!(ax, xlim_val)
    else
        xlims!(ax, (times[1], times[end]))
    end

    ylim_val = _get_val(:ylim)
    if !isnothing(ylim_val)
        ylims!(ax, ylim_val)
    else
        # Auto-determine y limits
        y_min = minimum(accuracy)
        y_max = maximum(accuracy)
        if show_error && !isnothing(stderror)
            y_min = min(y_min, minimum(accuracy .- stderror))
            y_max = max(y_max, maximum(accuracy .+ stderror))
        end
        y_range = y_max - y_min
        ylims!(ax, (y_min - 0.05 * y_range, y_max + 0.05 * y_range))
    end

    # Add origin lines if requested
    if _get_val(:add_xy_origin)
        vlines!(ax, 0, color = :black, linewidth = 1, linestyle = :dash)
        hlines!(ax, chance_level, color = _get_val(:chance_color), linewidth = _get_val(:chance_linewidth), linestyle = _get_val(:chance_linestyle))
    end

    # Plot chance level line
    if show_chance
        hlines!(
            ax,
            chance_level,
            color = _get_val(:chance_color),
            linewidth = _get_val(:chance_linewidth),
            linestyle = _get_val(:chance_linestyle),
            label = "Chance ($(round(chance_level, digits = 3)))",
        )
    end

    # Plot error shading if available
    if show_error
        band!(
            ax,
            times,
            accuracy .- stderror,
            accuracy .+ stderror,
            color = (_get_val(:error_color), _get_val(:error_alpha)),
            label = "±1 SE",
        )
    end

    # Plot main accuracy curve
    lines!(
        ax,
        times,
        accuracy,
        color = _get_val(:color),
        linewidth = _get_val(:linewidth),
        linestyle = _get_val(:linestyle),
        label = "Accuracy",
    )

    # Add legend if multiple elements
    if show_chance || show_error
        axislegend(ax, position = :rt)
    end

    # Display if requested
    if display_plot
        display(fig)
    end

    return fig
end

"""
    plot_decoding(decoded_list::Vector{DecodedData}; kwargs...)

Plot multiple decoding results on the same axes.

Useful for comparing decoding across different conditions, participants,
or analysis parameters.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects
- `kwargs`: Additional keyword arguments
- `colors::Union{Vector, Symbol}`: Colors for each curve (default: automatic)
- `labels::Vector{String}`: Labels for each curve (default: from condition_names)

# Examples
```julia
# Compare two decoding analyses
plot_decoding([decoded1, decoded2], colors=[:red, :blue], labels=["Condition A", "Condition B"])
```
"""
function plot_decoding(decoded_list::Vector{DecodedData}; kwargs...)
    if isempty(decoded_list)
        @minimal_error_throw("Cannot plot empty decoded data list")
    end

    # Merge defaults with user kwargs
    plot_kwargs = merge(PLOT_DECODING_KWARGS, Dict(kwargs))

    # Extract parameters (unwrap tuples if needed)
    _get_val(k) = isa(plot_kwargs[k], Tuple) ? plot_kwargs[k][1] : plot_kwargs[k]
    display_plot = _get_val(:display_plot)
    figure_title = _get_val(:figure_title)
    colors = get(kwargs, :colors, nothing)
    labels = get(kwargs, :labels, nothing)

    # Auto-generate colors if not provided
    if isnothing(colors)
        colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]
        # Cycle if needed
        colors = [colors[mod1(i, length(colors))] for i in 1:length(decoded_list)]
    elseif colors isa Symbol
        colors = fill(colors, length(decoded_list))
    end

    # Auto-generate labels if not provided
    if isnothing(labels)
        labels = [join(d.condition_names, " vs ") for d in decoded_list]
    end

    # Create figure
    fig = Figure(title = figure_title, resolution = (800, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = _get_val(:xlabel),
        ylabel = _get_val(:ylabel),
        xgridvisible = _get_val(:xgrid),
        ygridvisible = _get_val(:ygrid),
    )

    # Determine common time range and y limits
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
    ylims!(ax, (y_min - 0.05 * y_range, y_max + 0.05 * y_range))
    xlims!(ax, (time_min, time_max))

    # Get chance level (use first decoded's chance level, or average if different)
    chance_level = decoded_list[1].chance_level

    # Add origin lines
    if _get_val(:add_xy_origin)
        vlines!(ax, 0, color = :black, linewidth = 1, linestyle = :dash)
        hlines!(ax, chance_level, color = _get_val(:chance_color), linewidth = _get_val(:chance_linewidth), linestyle = _get_val(:chance_linestyle))
    end

    # Plot chance level
    if _get_val(:show_chance)
        hlines!(
            ax,
            chance_level,
            color = _get_val(:chance_color),
            linewidth = _get_val(:chance_linewidth),
            linestyle = _get_val(:chance_linestyle),
            label = "Chance ($(round(chance_level, digits = 3)))",
        )
    end

    # Plot each decoding result
    for (idx, decoded) in enumerate(decoded_list)
        times = decoded.times
        accuracy = decoded.average_score
        stderror = decoded.stderror

        # Error shading
        if _get_val(:show_error) && !isnothing(stderror)
            band!(
                ax,
                times,
                accuracy .- stderror,
                accuracy .+ stderror,
                color = (colors[idx], _get_val(:error_alpha)),
            )
        end

        # Main curve
        lines!(
            ax,
            times,
            accuracy,
            color = colors[idx],
            linewidth = _get_val(:linewidth),
            linestyle = _get_val(:linestyle),
            label = labels[idx],
        )
    end

    # Add legend
    axislegend(ax, position = :rt)

    # Display if requested
    if display_plot
        display(fig)
    end

    return fig
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
    fig = Figure(resolution = (600, 600))
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

