"""
Time-frequency plotting functions for visualizing TimeFreqData.
"""

# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_TF_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display
    :display_plot => (true, "Whether to display the plot"),
    :theme_fontsize => (24, "Font size for theme"),

    # Layout
    :layout => (:single, "Layout type: :single, :grid, or :topo"),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),

    # Colormap and color range
    :colormap => (:viridis, "Colormap for the heatmap"),
    :colorrange => (nothing, "Color range as (min, max) tuple. If nothing, automatically determined from data"),
    # Colorbar parameters - get all Colorbar attributes with their actual defaults
    [
        Symbol("colorbar_$(attr)") => (get(COLORBAR_DEFAULTS, attr, nothing), "Colorbar $(attr) parameter") for
        attr in propertynames(Colorbar)
    ]...,

    # Specific colorbar overrides for tf
    :colorbar_plot => (true, "Whether to show the colorbar"),
    :colorbar_position => (:right, "Position of the colorbar (:right, :left, :top, :bottom, or tuple)"),
    :colorbar_label => (nothing, "Custom colorbar label. If nothing, automatically determined"),

    # Axis
    :ylogscale => (false, "Whether to use logarithmic scale for the frequency (y) axis"),
    :xlim => (nothing, "X-axis time limits as (min, max) tuple in seconds. If nothing, shows all time points"),
    :interpolate => (false, "Whether to interpolate the heatmap for a smoother appearance"),
    :plot_title => (nothing, "Plot title. If nothing, automatically determined from condition name and channel"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :xticks => (nothing, "Custom x-axis ticks (e.g., -0.2:0.2:1.0)"),
    :yticks => (nothing, "Custom y-axis ticks (e.g., [2, 10, 20, 40, 80])"),
    :time_unit => (
        :s,
        "Time unit for x-axis display (:s or :ms). Only affects axis labels and tick formatting — all intervals remain in seconds.",
    ),

    # Baseline
    :baseline_method => (:db, "Baseline correction method: :db, :absolute, :relative, :relchange, :percent, :zscore"),

    # Layout parameters - dynamically pull all layout options
    [Symbol("layout_$(attr)") => val for (attr, val) in LAYOUT_KWARGS]...,
)


"""
    plot_tf(filepath::String; kwargs...)
    plot_tf(tf::TimeFreqData; channel_selection, baseline_interval, kwargs...)
    plot_tf(tfs::Vector{TimeFreqData}; channel_selection, baseline_interval, kwargs...)

Plot time-frequency data as a heatmap (time × frequency). Supports single condition,
multi-condition grid, and topographic layouts.

# Arguments
- `filepath::String`: Path to a `.jld2` file, or a `TimeFreqData` / `Vector{TimeFreqData}` object
- `channel_selection::Function`: Channel filter (default: all channels)
- `baseline_interval`: Baseline window as `(start, stop)` in seconds (default: `nothing`)

$(_generate_kwargs_doc(PLOT_TF_KWARGS))

# Examples
```julia
plot_tf("tf_morlet_result.jld2")
plot_tf("tf_morlet_result.jld2"; layout=:grid, channel_selection=channels([:Cz, :Pz]))
plot_tf(tf; colormap=:RdBu, ylogscale=true, baseline_interval=(-0.3, 0.0))
```
"""
function plot_tf(filepath::String; kwargs...)
    data = read_data(filepath)
    if isnothing(data)
        @minimal_error "No data found in file: $filepath"
    end
    return plot_tf(data; kwargs...)
end
function plot_tf(
    tf_data::Vector{TimeFreqData};
    channel_selection::Function = channels(),
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    kwargs...,
)
    isempty(tf_data) && error("Empty TimeFreqData vector")

    plot_kwargs = _merge_plot_kwargs(PLOT_TF_KWARGS, kwargs)

    # Extract colorbar/layout kwargs
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    colormap = plot_kwargs[:colormap]
    colorrange = plot_kwargs[:colorrange]
    colorbar = pop!(plot_kwargs, :colorbar_plot, true)
    ylogscale = plot_kwargs[:ylogscale]
    xlim = plot_kwargs[:xlim]
    xticks = plot_kwargs[:xticks]
    yticks = plot_kwargs[:yticks]
    interpolate = plot_kwargs[:interpolate]
    baseline_method = plot_kwargs[:baseline_method]
    time_unit = plot_kwargs[:time_unit]
    grid_dims = get(layout_kwargs, :grid_dims, nothing)

    n = length(tf_data)
    rows, cols = isnothing(grid_dims) ? _best_rect(n) : grid_dims
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])
    fig = Figure(size = (max(600, cols * 400), max(400, rows * 350)), figure_padding = plot_kwargs[:figure_padding])
    _set_window_title("$(basename(first(tf_data).file)) — Time-Frequency ($(n) conditions)")

    # Apply baseline to all conditions
    tf_plots = map(tf_data) do tf
        if !isnothing(baseline_interval) && !isnothing(tf.baseline)
            tf
        elseif !isnothing(baseline_interval)
            tf_baseline(tf, baseline_interval; method = baseline_method)
        else
            tf
        end
    end

    # Resolve channel (same for all conditions — use first condition's channels)
    all_channels = channel_labels(first(tf_plots))
    selected_mask = channel_selection(all_channels)
    selected_channels = all_channels[selected_mask]
    isempty(selected_channels) && error("No channels matched. Available: $(all_channels)")
    channel = first(selected_channels)

    # Shared times/freqs (assume consistent across conditions)
    all_times = Float64.(sort(unique(first(tf_plots).data_power.time)))
    freqs_vec = Float64.(sort(unique(first(tf_plots).data_power.freq)))

    # Pre-compute shared color range across all conditions
    if isnothing(colorrange)
        global_min = Inf
        global_max = -Inf
        for tf in tf_plots
            power_mat = _tf_df_to_matrix(tf.data_power, channel, freqs_vec, all_times)
            valid = filter(!isnan, power_mat)
            if !isempty(valid)
                cmin, cmax = extrema(valid)
                global_min = min(global_min, cmin)
                global_max = max(global_max, cmax)
            end
        end
        shared_colorrange = (global_min, global_max)
    else
        shared_colorrange = colorrange
    end

    # Plot each condition in its grid cell
    axes = Axis[]
    last_hm = nothing
    for (idx, tf) in enumerate(tf_plots)
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)

        hm = _plot_tf_heatmap!(
            ax,
            tf,
            channel,
            freqs_vec,
            all_times;
            colormap = colormap,
            colorrange = shared_colorrange,
            ylogscale = ylogscale,
            xlim = xlim,
            interpolate = interpolate,
            xticks = xticks,
            yticks = yticks,
            time_unit = time_unit,
        )
        last_hm = hm
        ax.title = "$(tf.condition_name) — $channel"

        # Only show y-labels on leftmost, x-labels on bottom
        col != 1 && (ax.ylabel = ""; ax.yticklabelsvisible = false)
        row != rows && (ax.xlabel = ""; ax.xticklabelsvisible = false)
    end

    length(axes) > 1 && linkaxes!(axes...)

    if colorbar && !isnothing(last_hm)
        cb_label =
            isnothing(plot_kwargs[:colorbar_label]) ? _tf_colorbar_label(first(tf_plots), baseline_interval, baseline_method) :
            plot_kwargs[:colorbar_label]

        cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:rows, 1:cols)
        Colorbar(fig[cb_pos...], last_hm; label = cb_label, colorbar_kwargs...)
    end

    _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

function plot_tf(
    tf_data::TimeFreqData;
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    kwargs...,
)

    plot_kwargs = _merge_plot_kwargs(PLOT_TF_KWARGS, kwargs)

    # Extract colorbar/layout kwargs
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    layout_kwargs   = _extract_layout_kwargs(plot_kwargs)

    layout                  = plot_kwargs[:layout]
    colormap                = plot_kwargs[:colormap]
    colorrange              = plot_kwargs[:colorrange]
    title                   = plot_kwargs[:plot_title]
    colorbar                = pop!(plot_kwargs, :colorbar_plot, true)
    ylogscale               = plot_kwargs[:ylogscale]
    xlim                    = plot_kwargs[:xlim]
    xticks                  = plot_kwargs[:xticks]
    yticks                  = plot_kwargs[:yticks]
    interpolate             = plot_kwargs[:interpolate]
    baseline_method         = plot_kwargs[:baseline_method]
    time_unit               = plot_kwargs[:time_unit]
    layout_grid_dims        = get(layout_kwargs, :grid_dims, nothing)
    layout_grid_rowgap      = get(layout_kwargs, :grid_rowgap, 10)
    layout_grid_colgap      = get(layout_kwargs, :grid_colgap, 10)
    layout_topo_plot_width  = get(layout_kwargs, :topo_plot_width, 0.1)
    layout_topo_plot_height = get(layout_kwargs, :topo_plot_height, 0.1)

    # Apply baseline if requested, but only if data hasn't already been baselined
    if !isnothing(baseline_interval) && !isnothing(tf_data.baseline)
        @warn "Data has already been baselined (method: $(tf_data.baseline.method), window: $(tf_data.baseline.window)). " *
              "Ignoring baseline_interval parameter. Use the data as-is or create a new TimeFreqData without baseline."
        tf_plot = tf_data
    elseif !isnothing(baseline_interval)
        tf_plot = tf_baseline(tf_data, baseline_interval; method = baseline_method)
    else
        tf_plot = tf_data
    end

    # Resolve channels via channel_selection predicate
    all_channels = channel_labels(tf_plot)
    selected_mask = channel_selection(all_channels)
    selected_channels = all_channels[selected_mask]
    isempty(selected_channels) && error("No channels matched. Available: $(all_channels)")

    # Apply user-specified plotting order if provided
    if !isnothing(channel_plot_order)
        ordered = Symbol[ch for ch in channel_plot_order if ch in selected_channels]
        if !isempty(ordered)
            selected_channels = ordered
        end
    end

    # For :single layout, take only the first channel
    plot_channels = layout === :single ? [first(selected_channels)] : selected_channels

    # Get unique times and frequencies (shared across all channels)
    all_times = Float64.(sort(unique(tf_plot.data_power.time)))
    freqs_vec = Float64.(sort(unique(tf_plot.data_power.freq)))

    # Create layout and figure
    layout_kwargs = Dict{Symbol,Any}(
        :grid_dims => layout_grid_dims,
        :grid_rowgap => layout_grid_rowgap,
        :grid_colgap => layout_grid_colgap,
        :topo_plot_width => layout_topo_plot_width,
        :topo_plot_height => layout_topo_plot_height,
    )

    # Determine figure size based on layout
    fig_size = if layout === :single
        (800, 500)
    elseif layout === :topo
        (1200, 1000)
    else # :grid
        n = length(plot_channels)
        rows, cols = !isnothing(layout_grid_dims) ? layout_grid_dims : _best_rect(n)
        (max(600, cols * 350), max(400, rows * 300))
    end

    _set_window_title("$(basename(tf_data.file)) — $(tf_data.condition_name) — Time-Frequency")
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])
    fig = Figure(size = fig_size, figure_padding = plot_kwargs[:figure_padding])

    eeg_layout = hasproperty(tf_data, :layout) ? tf_data.layout : nothing
    plot_layout = create_layout(layout, plot_channels, eeg_layout; layout_kwargs...)

    axes, layout_channels = _apply_layout!(fig, plot_layout; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)

    # Pre-compute shared color range across all channels if not specified
    if isnothing(colorrange)
        global_min = Inf
        global_max = -Inf
        for ch in plot_channels
            power_mat = _tf_df_to_matrix(tf_plot.data_power, ch, freqs_vec, all_times)
            valid = filter(!isnan, power_mat)
            if !isempty(valid)
                cmin, cmax = extrema(valid)
                global_min = min(global_min, cmin)
                global_max = max(global_max, cmax)
            end
        end
        shared_colorrange = (global_min, global_max)
    else
        shared_colorrange = colorrange
    end

    # Plot heatmap on each axis
    last_hm = nothing
    for (ax_idx, (ax, channel)) in enumerate(zip(axes, layout_channels))
        hm = _plot_tf_heatmap!(
            ax,
            tf_plot,
            channel,
            freqs_vec,
            all_times;
            colormap = colormap,
            colorrange = shared_colorrange,
            ylogscale = ylogscale,
            xlim = xlim,
            interpolate = interpolate,
            xticks = xticks,
            yticks = yticks,
            time_unit = time_unit,
        )
        last_hm = hm

        # Set title per axis
        if layout === :single
            ax.title = isnothing(title) ? "$(tf_data.condition_name) - $channel" : title
        else
            ax.title = string(channel)
        end
    end

    # Apply layout-specific axis properties (hide redundant labels for grid/topo)
    if plot_layout.type == :grid
        rows, cols = plot_layout.dims
        for (idx, ax) in enumerate(axes)
            row = fld(idx - 1, cols) + 1
            col = mod(idx - 1, cols) + 1
            # Only show y-axis labels on leftmost column
            if col != 1
                ax.ylabel = ""
                ax.yticklabelsvisible = false
            end
            # Only show x-axis labels on bottom row
            if row != rows
                ax.xlabel = ""
                ax.xticklabelsvisible = false
            end
        end
    elseif plot_layout.type == :topo
        for ax in axes
            hidedecorations!(ax, grid = false, ticks = true, ticklabels = true)
            hidespines!(ax)
        end
    end

    # Link axes for consistent zoom
    length(axes) > 1 && linkaxes!(axes...)

    # Add shared colorbar
    if colorbar && !isnothing(last_hm)
        cb_label =
            isnothing(plot_kwargs[:colorbar_label]) ? _tf_colorbar_label(tf_plot, baseline_interval, baseline_method) :
            plot_kwargs[:colorbar_label]

        if layout === :single
            cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
            Colorbar(fig[cb_pos...], last_hm; label = cb_label, colorbar_kwargs...)
        elseif layout === :grid
            rows, cols = plot_layout.dims
            cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:rows, 1:cols)
            Colorbar(fig[cb_pos...], last_hm; label = cb_label, colorbar_kwargs...)
        else # :topo
            # Topo is a bit weird because of the scale axis, let's just stick it on the right
            cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
            Colorbar(fig[cb_pos[1], cb_pos[2]+1], last_hm; label = cb_label, height = Relative(0.5), colorbar_kwargs...)
        end
    end

    _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end


# === INTERNAL HELPERS ===

"""
    _plot_tf_heatmap!(ax, tf_plot, channel, freqs_vec, times; kwargs...)

Plot a single TF heatmap on the given axis. Returns the heatmap object.
"""
function _plot_tf_heatmap!(
    ax::Axis,
    tf_plot::TimeFreqData,
    channel::Symbol,
    freqs_vec::Vector{Float64},
    times::Vector{Float64};
    colormap = :viridis,
    colorrange::Tuple{Real,Real} = (-1.0, 1.0),
    ylogscale::Bool = false,
    xlim::Union{Nothing,Tuple{Real,Real}} = nothing,
    interpolate::Bool = false,
    xticks = nothing,
    yticks = nothing,
    time_unit::Symbol = :s,
)

    # Extract power matrix: [n_freqs × n_times]
    power_mat = _tf_df_to_matrix(tf_plot.data_power, channel, freqs_vec, times)

    # Configure axis
    if time_unit == :ms
        ax.xlabel = "Time (ms)"
        ax.xtickformat = values -> [string(round(Int, v * 1000)) for v in values]
    else
        ax.xlabel = "Time (s)"
    end
    ax.ylabel = "Frequency (Hz)"

    if ylogscale
        ylims!(ax, (minimum(freqs_vec), maximum(freqs_vec)))
        ax.yscale = log10
        ax.ytickformat = values -> [string(round(Int, v)) for v in values]
    end

    !isnothing(xlim) && xlims!(ax, xlim)
    !isnothing(xticks) && (ax.xticks = xticks)
    !isnothing(yticks) && (ax.yticks = yticks)

    # Plot heatmap (expects n_times × n_freqs via zero-allocation lazy transpose)
    hm = heatmap!(
        ax,
        times,
        freqs_vec,
        transpose(power_mat),
        colormap = colormap,
        colorrange = colorrange,
        nan_color = :transparent,
        interpolate = interpolate,
    )

    return hm
end


"""
    _tf_colorbar_label(tf_plot, baseline_interval, baseline_method)

Determine the colorbar label based on baseline information.
"""
function _tf_colorbar_label(tf_plot::TimeFreqData, baseline_interval, baseline_method::Symbol)
    if !isnothing(tf_plot.baseline)
        method = tf_plot.baseline.method
        return method == :db ? "Power (dB)" : method == :percent ? "Power (% change)" : method == :relchange ? "Power (relative)" : "Power"
    elseif !isnothing(baseline_interval)
        return baseline_method == :db ? "Power (dB)" :
               baseline_method == :percent ? "Power (% change)" : baseline_method == :relchange ? "Power (relative)" : "Power"
    else
        return "Power"
    end
end
