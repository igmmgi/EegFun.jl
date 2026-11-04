# Plotting functions for MVPA/decoding results.

const PLOT_DECODING_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("", "Title drawn at the top of the entire figure canvas"),
    :figure_title_fontsize => (24, "Font size for figure title"),
    :interactive => (true, "Enable interactive features (true/false)"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("Classification Accuracy", "Label for y-axis"),

    # Title
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Show plot title (true/false)"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Line styling
    :linewidth => (2, "Line width for decoding curve"),
    :color => (nothing, "Color for decoding curve"),
    :linestyle => (nothing, "Line style"),

    # Chance level
    :chance_color => (:gray, "Color for chance level line"),
    :chance_linestyle => (:dash, "Line style for chance level"),
    :chance_linewidth => (1, "Line width for chance level"),

    # Error bars
    :show_error => (true, "Show standard error bars (true/false)"),
    :error_color => (:blue, "Color for error bars"),
    :error_alpha => (0.3, "Transparency for error shading"),

    # Significance markers
    :show_significance => (true, "Show significance markers when statistics are provided (true/false)"),
    :significance_position => (:bottom, "Vertical position of significance bars (:bottom, :top, or numeric y value)"),
    :sig_color => (:black, "Color for significance markers"),
    :sig_alpha => (0.5, "Transparency for significance markers"),

    # Legend & Colorbar / Heatmap
    :legend => (false, "Show legend (true/false)"),
    :colormap => (nothing, "Colormap for confusion matrix"),
    :colorrange => (nothing, "Color scale limits for confusion matrix (default: (0, 1))"),
    :colorbar_plot => (true, "Show colorbar on confusion matrix (true/false)"),

    # Grid
    :xgrid => (true, "Show x-axis grid (true/false)"),
    :ygrid => (true, "Show y-axis grid (true/false)"),

    # Origin lines
    :add_xy_origin => (true, "Add origin lines at x=0 and y=chance (true/false)"),

    # Figure padding
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),
)

function _decoding_ylim_bounds(accuracy::Vector{Float64}, stderror::Union{Vector{Float64},Nothing}, show_error::Bool)
    y_min = Inf
    y_max = -Inf
    if show_error && !isnothing(stderror)
        @inbounds for i in eachindex(accuracy, stderror)
            a = accuracy[i]
            s = stderror[i]
            y_min = min(y_min, a - s)
            y_max = max(y_max, a + s)
        end
    else
        @inbounds for i in eachindex(accuracy)
            a = accuracy[i]
            y_min = min(y_min, a)
            y_max = max(y_max, a)
        end
    end
    return y_min, y_max
end

function _decoding_auto_ylim(accuracy::Vector{Float64}, stderror::Union{Vector{Float64},Nothing}, show_error::Bool)
    y_min, y_max = _decoding_ylim_bounds(accuracy, stderror, show_error)
    y_range = y_max - y_min
    pad = y_range == 0 ? 0.05 : 0.05 * y_range
    return (y_min - pad, y_max + pad)
end

function _add_decoding_origin_lines!(ax::Axis, chance_level::Float64, plot_kwargs::Dict)
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

function _plot_error_band!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64},Nothing},
    plot_kwargs::Dict,
)
    band!(ax, times, accuracy .- stderror, accuracy .+ stderror, color = (plot_kwargs[:error_color], plot_kwargs[:error_alpha]))
end

function _plot_accuracy_curve!(ax::Axis, times::Vector{Float64}, accuracy::Vector{Float64}, plot_kwargs::Dict; label::String = "Accuracy")
    line_color = _resolve_theme_linecolor(ax, plot_kwargs[:color])
    line_style = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]
    lines!(ax, times, accuracy, color = line_color, linewidth = plot_kwargs[:linewidth], linestyle = line_style, label = label)
end

function _setup_axis_limits!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64},Nothing},
    plot_kwargs::Dict,
)
    xlim_val = plot_kwargs[:xlim]
    xlims!(ax, isnothing(xlim_val) ? (times[1], times[end]) : xlim_val)

    ylim_val = plot_kwargs[:ylim]
    ylims!(ax, isnothing(ylim_val) ? _decoding_auto_ylim(accuracy, stderror, Bool(plot_kwargs[:show_error])) : ylim_val)
end

function _plot_decoding_to_axis!(
    ax::Axis,
    times::Vector{Float64},
    accuracy::Vector{Float64},
    stderror::Union{Vector{Float64},Nothing},
    chance_level::Float64,
    plot_kwargs::Dict;
    set_limits::Bool = true,
)
    set_limits && _setup_axis_limits!(ax, times, accuracy, stderror, plot_kwargs)
    plot_kwargs[:add_xy_origin] && _add_decoding_origin_lines!(ax, chance_level, plot_kwargs)
    (plot_kwargs[:show_error] && !isnothing(stderror)) && _plot_error_band!(ax, times, accuracy, stderror, plot_kwargs)
    _plot_accuracy_curve!(ax, times, accuracy, plot_kwargs)
    get(plot_kwargs, :legend, false) && axislegend(ax, position = :rt)
end

function _finalize_decoding_figure!(fig::Figure, plot_kwargs::Dict, window_title::String; display_plot::Bool = plot_kwargs[:display_plot])
    if !isempty(plot_kwargs[:figure_title])
        Label(fig[0, :], plot_kwargs[:figure_title], fontsize = plot_kwargs[:figure_title_fontsize], font = :bold, tellwidth = false)
    end
    _set_window_title(window_title)
    display_plot && _display_figure(fig)
    _set_window_title("Makie")
    return nothing
end

"""
    plot_decoding(decoded::DecodedData; kwargs...)

Plot decoding accuracy over time.

Creates a line plot showing classification accuracy at each time point,
with optional error bars and chance level reference line.

# Arguments
- `decoded::DecodedData`: DecodedData object containing decoding results
- `kwargs`: Additional keyword arguments (see PLOT_DECODING_KWARGS)

$(_generate_kwargs_doc(PLOT_DECODING_KWARGS))

# Examples
```julia
# Basic plot
plot_decoding(decoded)

# Custom styling
plot_decoding(decoded, color=:red, linewidth=3, show_error=false)

# With custom title
plot_decoding(decoded, plot_title="Face vs. Object Decoding")
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
function plot_decoding(filepath::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_decoding(data; kwargs...)
    end

    files = _find_batch_files(filepath, input_dir, participant_selection)
    isempty(files) && @minimal_error "No files matching pattern '$filepath' in $input_dir"

    results = NamedTuple[]
    for file in sort(files, by = _natural_sort_key)
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(results, plot_decoding(data; kwargs...))
    end
    return results
end

function plot_decoding(decoded::DecodedData; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)
    title_text = plot_kwargs[:plot_title]

    fig = _eegfun_figure(size = (800, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
        xgridvisible = plot_kwargs[:xgrid],
        ygridvisible = plot_kwargs[:ygrid],
    )

    if plot_kwargs[:show_title]
        if isnothing(title_text)
            ax.title = "Decoding: $(join(decoded.condition_names, " vs "))"
        elseif !isempty(title_text)
            ax.title = title_text
        end
    end

    _plot_decoding_to_axis!(ax, decoded.times, decoded.average_score, decoded.stderror, decoded.parameters.chance_level, plot_kwargs)
    _finalize_decoding_figure!(fig, plot_kwargs, "Decoding")
    return (fig = fig, axes = [ax])
end

function plot_decoding(decoded_list::AbstractVector{<:DecodedData}; kwargs...)
    isempty(decoded_list) && @minimal_error("Cannot plot empty decoded data list")

    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)
    title_text = plot_kwargs[:plot_title]
    show_title = plot_kwargs[:show_title]
    show_err = Bool(plot_kwargs[:show_error])

    n_subjects = length(decoded_list)
    rows, cols = _best_rect(n_subjects)
    fig = _eegfun_figure(size = (400 * cols, 300 * rows))

    time_min, time_max = Inf, -Inf
    y_min, y_max = Inf, -Inf
    for d in decoded_list
        time_min = min(time_min, d.times[1])
        time_max = max(time_max, d.times[end])
        d_min, d_max = _decoding_ylim_bounds(d.average_score, d.stderror, show_err)
        y_min = min(y_min, d_min)
        y_max = max(y_max, d_max)
    end

    y_range = y_max - y_min
    pad = y_range == 0 ? 0.05 : 0.05 * y_range
    shared_xlim = isnothing(plot_kwargs[:xlim]) ? (time_min, time_max) : plot_kwargs[:xlim]
    shared_ylim = isnothing(plot_kwargs[:ylim]) ? (y_min - pad, y_max + pad) : plot_kwargs[:ylim]

    axes = Vector{Axis}(undef, n_subjects)
    for (idx, decoded) in enumerate(decoded_list)
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1

        ax = Axis(
            fig[row, col],
            xlabel = row == rows ? plot_kwargs[:xlabel] : "",
            ylabel = col == 1 ? plot_kwargs[:ylabel] : "",
            xgridvisible = plot_kwargs[:xgrid],
            ygridvisible = plot_kwargs[:ygrid],
        )

        if show_title
            if isnothing(title_text)
                ax.title = isnothing(decoded.file) ? "Subject $idx" : basename(decoded.file)
            elseif !isempty(title_text)
                ax.title = title_text
            end
        end

        xlims!(ax, shared_xlim)
        ylims!(ax, shared_ylim)
        _plot_decoding_to_axis!(
            ax,
            decoded.times,
            decoded.average_score,
            decoded.stderror,
            decoded.parameters.chance_level,
            plot_kwargs;
            set_limits = false,
        )
        axes[idx] = ax
    end

    _finalize_decoding_figure!(fig, plot_kwargs, "Decoding")
    return (fig = fig, axes = axes)
end

function plot_decoding(decoded::DecodedData, stats::DecodingStatisticsResult; kwargs...)
    decoded.times == stats.times || @minimal_error("DecodedData and DecodingStatisticsResult must have matching time vectors")

    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)

    result = plot_decoding(decoded; kwargs..., display_plot = false)
    fig = result.fig
    ax = result.axes[1]

    y_min, y_max =
        !isnothing(plot_kwargs[:ylim]) ? Float64.(plot_kwargs[:ylim]) :
        _decoding_auto_ylim(decoded.average_score, decoded.stderror, Bool(plot_kwargs[:show_error]))
    y_range = y_max - y_min
    y_range == 0 && (y_range = 0.1)

    if plot_kwargs[:show_significance] && any(stats.significant_mask)
        sig_color = plot_kwargs[:sig_color]
        sig_alpha = plot_kwargs[:sig_alpha]
        sig_pos = plot_kwargs[:significance_position]

        bar_y = if sig_pos == :top
            y_max - 0.04 * y_range
        elseif sig_pos == :zero
            decoded.parameters.chance_level
        elseif sig_pos isa Number
            Float64(sig_pos)
        else
            y_min + 0.02 * y_range
        end
        bar_height = y_range * 0.02

        for (t_start, t_end) in _find_continuous_regions(stats.significant_mask, stats.times)
            poly!(ax, Rect(t_start, bar_y, t_end - t_start, bar_height), color = (sig_color, sig_alpha), strokewidth = 0)
        end

        if !isnothing(stats.clusters) && !isempty(stats.clusters)
            for cluster in stats.clusters
                cluster.is_significant || continue
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

    _set_window_title("Decoding")
    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end

"""
    plot_confusion_matrix(decoded::DecodedData; time_point::Union{Float64, Int, Nothing} = nothing, kwargs...)

Plot confusion matrix for decoding results, where rows (`i`) correspond to True class
and columns (`j`) correspond to Predicted class.

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
function plot_confusion_matrix(decoded::DecodedData; time_point::Union{Float64,Int,Nothing} = nothing, kwargs...)
    isnothing(decoded.confusion_matrix) && @minimal_error("No confusion matrix data available in DecodedData")

    plot_kwargs = _merge_plot_kwargs(PLOT_DECODING_KWARGS, kwargs)
    cm = decoded.confusion_matrix

    if isnothing(time_point)
        confusion = dropdims(mean(cm, dims = 1), dims = 1)
        default_title = "Average Confusion Matrix"
    elseif time_point isa Float64
        time_idx = find_closest_time_index(decoded.times, time_point)
        confusion = @view cm[time_idx, :, :]
        default_title = "Confusion Matrix at $(round(time_point, digits=3)) s"
    else
        confusion = @view cm[time_point, :, :]
        default_title = "Confusion Matrix at $(round(decoded.times[time_point], digits=3)) s"
    end

    title_text = plot_kwargs[:plot_title]
    title_str = !plot_kwargs[:show_title] ? "" : (isnothing(title_text) ? default_title : title_text)
    n_classes = size(confusion, 1)

    # yreversed = true puts row 1 (True class 1) at top so diagonal runs top-left to bottom-right
    fig = _eegfun_figure(size = (650, 600))
    ax = Axis(fig[1, 1], title = title_str, xlabel = "Predicted", ylabel = "True", yreversed = true, aspect = DataAspect())

    # Makie heatmap!(ax, xs, ys, Z) maps Z[x, y] to (x = xs[x], y = ys[y]); permutedims maps confusion[i, j] to (x = j, y = i)
    cmap_val = isnothing(plot_kwargs[:colormap]) ? :Blues : plot_kwargs[:colormap]
    cr_val = isnothing(plot_kwargs[:colorrange]) ? (0.0, 1.0) : plot_kwargs[:colorrange]
    hm = heatmap!(ax, 1:n_classes, 1:n_classes, permutedims(confusion), colormap = cmap_val, colorrange = cr_val)

    plot_kwargs[:colorbar_plot] && Colorbar(fig[1, 2], hm, label = "Proportion")

    cr_mid = (cr_val[1] + cr_val[2]) / 2
    for i = 1:n_classes, j = 1:n_classes
        val = confusion[i, j]
        text!(ax, j, i, text = string(round(val, digits = 2)), color = val > cr_mid ? :white : :black, align = (:center, :center))
    end

    ax.xticks = (1:n_classes, decoded.condition_names)
    ax.yticks = (1:n_classes, decoded.condition_names)

    _finalize_decoding_figure!(fig, plot_kwargs, "Confusion Matrix")
    return (fig = fig, axes = [ax])
end
