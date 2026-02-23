"""
Time-frequency plotting functions for visualizing TimeFreqData.
"""

"""
    plot_tf(tf_data::TimeFreqData;
            channel_selection::Function = channels(),
            baseline_interval=nothing, baseline_method=:db,
            colormap=:viridis, colorrange=nothing,
            title=nothing, colorbar=true, ylogscale=false)

Plot time-frequency data for a specific channel.

# Arguments
- `tf_data::TimeFreqData`: Time-frequency data

# Keyword Arguments
- `channel_selection::Function`: Channel selection predicate (default: all channels).
  The first selected channel is plotted.
  - Example: `channel_selection=channels(:Cz)` for a specific channel
  - Example: `channel_selection=channels([:Cz, :Pz])` plots the first match
- `baseline_interval`: Optional baseline interval (start, stop) in seconds
- `baseline_method`: Baseline method if baseline_interval provided
- `colormap`: Colormap (default: :viridis)
- `colorrange`: Color range (default: auto)
- `title`: Plot title (default: auto)
- `colorbar`: Show colorbar (default: true)
- `ylogscale`: Use logarithmic scale for y-axis (frequencies) (default: false)

# Returns
- `(Figure, Axis)`: Makie figure and axis

# Example
```julia
fig, ax = plot_tf(tf_data; channel_selection=channels(:Cz), baseline_interval=(-0.3, 0.0), ylogscale=true)
```
"""
function plot_tf(
    tf_data::TimeFreqData;
    channel_selection::Function = channels(),
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    baseline_method::Symbol = :db,
    colormap = :viridis,
    colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
    title::Union{Nothing,String} = nothing,
    colorbar::Bool = true,
    ylogscale::Bool = false,
)

    # Apply baseline if requested, but only if data hasn't already been baselined
    if !isnothing(baseline_interval) && tf_data.baseline !== nothing
        @warn "Data has already been baselined (method: $(tf_data.baseline.method), window: $(tf_data.baseline.window)). " *
              "Ignoring baseline_interval parameter. Use the data as-is or create a new TimeFreqData without baseline."
        tf_plot = tf_data
    elseif !isnothing(baseline_interval)
        # Apply baseline on-the-fly
        tf_plot = tf_baseline(tf_data, baseline_interval; method = baseline_method)
    else
        tf_plot = tf_data
    end

    # Resolve channel via channel_selection predicate
    all_channels = channel_labels(tf_plot)
    selected_mask = channel_selection(all_channels)
    selected_channels = all_channels[selected_mask]
    isempty(selected_channels) && error("No channels matched. Available: $(all_channels)")
    channel = first(selected_channels)

    # Get unique times and frequencies (sorted)
    times = sort(unique(tf_plot.data_power.time))
    freqs_vec = sort(unique(tf_plot.data_power.freq))
    n_times = length(times)
    n_freqs = length(freqs_vec)

    # Extract power matrix using shared helper: [n_freqs × n_time]
    power_mat = _tf_df_to_matrix(tf_plot.data_power, channel, freqs_vec, times)

    # Transpose for Makie heatmap (expects n_times × n_freqs)
    # After transpose: rows = times, columns = frequencies
    power_mat = power_mat'

    # Create figure
    fig = Figure(size = (800, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
        title = isnothing(title) ? "$(tf_data.condition_name) - $channel" : title,
        yscale = ylogscale ? log10 : identity,
    )

    # Disable scientific notation for y-axis ticks when using log scale
    if ylogscale
        ax.ytickformat = values -> [string(round(Int, v)) for v in values]
    end

    # Determine color range
    if isnothing(colorrange)
        cmin, cmax = extrema(filter(!isnan, power_mat))
        colorrange = (cmin, cmax)
    end

    # Plot heatmap
    # Set NaN color to transparent so edge-filtered regions are not displayed
    hm = heatmap!(ax, times, freqs_vec, power_mat, colormap = colormap, colorrange = colorrange, nan_color = :transparent)

    if colorbar
        # Determine colorbar label based on baseline information
        if tf_plot.baseline !== nothing
            method = tf_plot.baseline.method
            if method == :db
                label = "Power (dB)"
            elseif method == :percent
                label = "Power (% change)"
            elseif method == :relchange
                label = "Power (relative)"
            else
                label = "Power"
            end
        elseif !isnothing(baseline_interval)
            # Baseline was just applied via parameter
            label =
                baseline_method == :db ? "Power (dB)" :
                baseline_method == :percent ? "Power (% change)" : baseline_method == :relchange ? "Power (relative)" : "Power"
        else
            label = "Power"
        end
        Colorbar(fig[1, 2], hm, label = label)
    end

    display(fig)
    return fig, ax
end
