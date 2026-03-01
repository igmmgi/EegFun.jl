"""
Time-frequency plotting functions for visualizing TimeFreqData.
"""

"""
    plot_tf(filepath::String; kwargs...)

Load time-frequency data from a JLD2 file and plot it.

# Arguments
- `filepath::String`: Path to JLD2 file containing TimeFreqData

# Examples
```julia
plot_tf("tf_morlet_result.jld2")
plot_tf("tf_morlet_result.jld2"; layout=:grid, channel_selection=channels([:Cz, :Pz]))
```
"""
function plot_tf(filepath::String; kwargs...)
    data = read_data(filepath)
    if isnothing(data)
        @minimal_error "No data found in file: $filepath"
    end
    return plot_tf(data; kwargs...)
end

"""
    plot_tf(tf_data::Vector{TimeFreqData}; kwargs...)

Plot multiple TimeFreqData conditions, each in a separate subplot.
The first selected channel is shown for each condition.

# Examples
```julia
plot_tf(tf_data_vector; channel_selection=channels(:Cz), baseline_interval=(-0.3, 0.0))
```
"""
function plot_tf(
    tf_data::Vector{TimeFreqData};
    channel_selection::Function = channels(),
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    baseline_method::Symbol = :db,
    colormap = :viridis,
    colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
    colorbar::Bool = true,
    ylogscale::Bool = false,
    xlim::Union{Nothing,Tuple{Real,Real}} = nothing,
    interpolate::Bool = false,
    grid_dims = nothing,
)
    isempty(tf_data) && error("Empty TimeFreqData vector")

    n = length(tf_data)
    rows, cols = isnothing(grid_dims) ? _best_rect(n) : grid_dims
    fig = Figure(size = (max(600, cols * 400), max(400, rows * 350)))
    _set_window_title("Time-Frequency — $(n) conditions")

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
        )
        last_hm = hm
        ax.title = "$(tf.condition_name) — $channel"

        # Only show y-labels on leftmost, x-labels on bottom
        col != 1 && (ax.ylabel = ""; ax.yticklabelsvisible = false)
        row != rows && (ax.xlabel = ""; ax.xticklabelsvisible = false)
    end

    length(axes) > 1 && linkaxes!(axes...)

    if colorbar && !isnothing(last_hm)
        cb_label = _tf_colorbar_label(first(tf_plots), baseline_interval, baseline_method)
        Colorbar(fig[1:rows, cols+1], last_hm, label = cb_label)
    end

    display(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end


"""
    plot_tf(tf_data::TimeFreqData;
            layout=:single, channel_selection=channels(),
            baseline_interval=nothing, baseline_method=:db,
            colormap=:viridis, colorrange=nothing,
            title=nothing, colorbar=true, ylogscale=false,
            xlim=nothing)

Plot time-frequency data with configurable layout for single or multi-channel display.

# Arguments
- `tf_data::TimeFreqData`: Time-frequency data

# Keyword Arguments
- `layout::Union{Symbol,PlotLayout}`: Layout type - `:single` (default), `:grid`, `:topo`, or a `PlotLayout` object
  - `:single` — plot first selected channel in a single pane
  - `:grid` — one heatmap per selected channel in a grid arrangement
  - `:topo` — channels positioned by their topographic coordinates
- `channel_selection::Function`: Channel selection predicate (default: all channels)
  - Example: `channel_selection=channels(:Cz)` for a specific channel
  - Example: `channel_selection=channels([:Cz, :Pz, :Oz])` for multiple channels
- `baseline_interval`: Optional baseline interval (start, stop) in seconds
- `baseline_method`: Baseline method if baseline_interval provided (`:db`, `:percent`, `:relchange`)
- `colormap`: Colormap (default: `:viridis`)
- `colorrange`: Color range (default: auto from data)
- `title`: Plot title (default: auto-generated)
- `colorbar`: Show colorbar (default: `true`)
- `ylogscale`: Use logarithmic scale for y-axis (default: `false`)
- `xlim`: X-axis limits as `(min, max)` tuple
- `layout_grid_dims`: Grid dimensions as `(rows, cols)` tuple (auto if nothing)
- `layout_grid_rowgap`: Gap between rows in grid layout (pixels)
- `layout_grid_colgap`: Gap between columns in grid layout (pixels)
- `layout_topo_plot_width`: Width of individual topo plots (fraction)
- `layout_topo_plot_height`: Height of individual topo plots (fraction)

# Returns
- `NamedTuple (fig, axes)`: Makie figure and vector of axes

# Examples
```julia
# Single channel (default)
plot_tf(tf_data; channel_selection=channels(:Cz), baseline_interval=(-0.3, 0.0))

# Grid layout with multiple channels
plot_tf(tf_data; layout=:grid, channel_selection=channels([:Fz, :Cz, :Pz, :Oz]))

# Topographic layout
plot_tf(tf_data; layout=:topo, channel_selection=channels([:Fz, :Cz, :Pz, :Oz]),
        baseline_interval=(-0.3, 0.0), ylogscale=true)

# Load from file
plot_tf("tf_morlet_result.jld2"; channel_selection=channels(:Cz))
```
"""
function plot_tf(
    tf_data::TimeFreqData;
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    baseline_method::Symbol = :db,
    colormap = :viridis,
    colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
    title::Union{Nothing,String} = nothing,
    colorbar::Bool = true,
    ylogscale::Bool = false,
    xlim::Union{Nothing,Tuple{Real,Real}} = nothing,
    interpolate::Bool = false,
    # Layout kwargs (passed through with layout_ prefix stripped)
    layout_grid_dims = nothing,
    layout_grid_rowgap = 10,
    layout_grid_colgap = 10,
    layout_topo_plot_width = 0.10,
    layout_topo_plot_height = 0.10,
)

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

    _set_window_title("$(tf_data.condition_name) — Time-Frequency")
    fig = Figure(size = fig_size)

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
        cb_label = _tf_colorbar_label(tf_plot, baseline_interval, baseline_method)
        if layout === :single
            Colorbar(fig[1, 2], last_hm, label = cb_label)
        elseif layout === :grid
            rows, cols = plot_layout.dims
            Colorbar(fig[1:rows, cols+1], last_hm, label = cb_label)
        else # :topo
            Colorbar(fig[1, 2], last_hm, label = cb_label, height = Relative(0.5))
        end
    end

    display(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

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
)

    # Extract power matrix: [n_freqs × n_times]
    power_mat = _tf_df_to_matrix(tf_plot.data_power, channel, freqs_vec, times)
    # Transpose for Makie heatmap (expects n_times × n_freqs)
    power_mat = power_mat'

    # Configure axis
    ax.xlabel = "Time (s)"
    ax.ylabel = "Frequency (Hz)"

    if ylogscale
        ylims!(ax, (minimum(freqs_vec), maximum(freqs_vec)))
        ax.yscale = log10
        ax.ytickformat = values -> [string(round(Int, v)) for v in values]
    end

    !isnothing(xlim) && xlims!(ax, xlim)

    # Plot heatmap
    hm = heatmap!(
        ax,
        times,
        freqs_vec,
        power_mat,
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
