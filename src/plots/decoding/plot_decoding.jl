"""
Plotting function for MVPA/decoding results.
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

    # Significance markers
    :sig_color => (:black, "Color for significance markers"),
    :sig_alpha => (0.5, "Transparency for significance markers"),

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
"""
function _add_decoding_origin_lines!(ax::Axis, chance_level::Float64, plot_kwargs::Dict)
    # Add x=0 line and y=chance_level line
    vlines!(ax, 0, color = :black, linewidth = 1, linestyle = :solid)
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
    band!(ax, times, accuracy .- stderror, accuracy .+ stderror, color = (plot_kwargs[:error_color], plot_kwargs[:error_alpha]))
end


"""
    _plot_accuracy_curve!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64}, 
                          plot_kwargs::Dict; label::String = "Accuracy")

Plot main accuracy curve on axis.
"""
function _plot_accuracy_curve!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64}, plot_kwargs::Dict; label::String = "Accuracy")
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
    stderror::Union{Vector{Float64},Nothing},
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
                            plot_kwargs::Dict)

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
"""
function _plot_decoding_to_axis!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64},Nothing},
    chance_level::Float64,
    plot_kwargs::Dict;
)

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
    _plot_accuracy_curve!(ax, times, accuracy, plot_kwargs)

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

    plot_decoding(filepath::String; input_dir=pwd(), participant_selection=participants(), kwargs...)

Load decoded data and plot. Accepts either a `.jld2` filepath or a pattern
to discover and plot all matching files (one plot per file).

# Examples
```julia
plot_decoding("decoded.jld2")
plot_decoding("decoded")
```
"""
function plot_decoding(filepath::String; input_dir::String = pwd(), participant_selection::Function = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_decoding(data; kwargs...)
    else
        files = _find_batch_files(filepath, input_dir, participant_selection)
        isempty(files) && @minimal_error "No files matching pattern '$filepath' in $input_dir"

        results = NamedTuple[]
        for file in sort(files, by = _natural_sort_key)
            file_path = joinpath(input_dir, file)
            @info "Plotting: $file"
            data = read_data(file_path)
            isnothing(data) && continue
            result = plot_decoding(data; kwargs...)
            push!(results, result)
        end
        return results
    end
end

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

    _set_window_title("Decoding")
    if display_plot
        display(fig)
    end
    _set_window_title("Makie")

    return fig, ax
end
function plot_decoding(decoded_list::Vector{DecodedData}; kwargs...)
    if isempty(decoded_list)
        @minimal_error("Cannot plot empty decoded data list")
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
    rows, cols = _best_rect(n_subjects)

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
        _plot_decoding_to_axis!(ax, times, accuracy, stderror, chance_level, plot_kwargs)
    end

    _set_window_title("Decoding")
    if display_plot
        display(fig)
    end
    _set_window_title("Makie")

    return fig
end
function plot_decoding(decoded::DecodedData, stats::DecodingStatisticsResult; kwargs...)
    # Merge defaults with user kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)

    # Extract parameters
    show_significance = get(kwargs, :show_significance, true)
    sig_color = plot_kwargs[:sig_color]
    sig_alpha = plot_kwargs[:sig_alpha]
    sig_bar_position = get(kwargs, :sig_bar_position, :bottom)

    # Validate that times match
    if decoded.times != stats.times
        @minimal_error("DecodedData and DecodingStatisticsResult must have matching time vectors")
    end

    # Create base plot (without displaying yet)
    display_plot_orig = plot_kwargs[:display_plot]
    fig, ax = plot_decoding(decoded; kwargs..., display_plot = false)

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
            poly!(ax, Rect(t_start, bar_y, t_end - t_start, bar_height), color = (sig_color, sig_alpha), strokewidth = 0)
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

    _set_window_title("Decoding")
    if display_plot_orig
        display(fig)
    end
    _set_window_title("Makie")

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
function plot_confusion_matrix(decoded::DecodedData; time_point::Union{Float64,Int,Nothing} = nothing, display_plot::Bool = true, kwargs...)
    if isnothing(decoded.confusion_matrix)
        @minimal_error("No confusion matrix data available in DecodedData")
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
    heatmap!(ax, confusion, colormap = :jet, colorrange = (0, 1))

    # Add text labels
    n_classes = size(confusion, 1)
    for i = 1:n_classes
        for j = 1:n_classes
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

    _set_window_title("Confusion Matrix")
    if display_plot
        display(fig)
    end
    _set_window_title("Makie")
    return fig
end

