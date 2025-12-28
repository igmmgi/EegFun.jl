"""
Time-frequency plotting functions for visualizing TimeFreqData.
"""

"""
    plot_time_frequency(tf_data::TimeFreqData, channel::Symbol;
                        baseline_window=nothing, baseline_method=:db,
                        colormap=:viridis, colorrange=nothing,
                        title=nothing, colorbar=true)

Plot time-frequency data for a specific channel.

# Arguments
- `tf_data::TimeFreqData`: Time-frequency data
- `channel::Symbol`: Channel to plot

# Keyword Arguments
- `baseline_window`: Optional baseline window (start, stop) in seconds
- `baseline_method`: Baseline method if baseline_window provided
- `colormap`: Colormap (default: :viridis)
- `colorrange`: Color range (default: auto)
- `title`: Plot title (default: auto)
- `colorbar`: Show colorbar (default: true)

# Returns
- `(Figure, Axis)`: Makie figure and axis

# Example
```julia
fig, ax = plot_time_frequency(tf_data, :Cz; baseline_window=(-0.3, 0.0))
```
"""
function plot_time_frequency(tf_data::TimeFreqData, channel::Symbol;
                             baseline_window::Union{Nothing,Tuple{Real,Real}}=nothing,
                             baseline_method::Symbol=:db,
                             colormap=:viridis,
                             colorrange::Union{Nothing,Tuple{Real,Real}}=nothing,
                             title::Union{Nothing,String}=nothing,
                             colorbar::Bool=true)
    
    # Apply baseline if requested
    tf_plot = isnothing(baseline_window) ? tf_data : tf_baseline(tf_data, baseline_window; method=baseline_method)
    
    # Get unique times and freqs
    times = sort(unique(tf_plot.data.time))
    freqs_vec = sort(unique(tf_plot.data.freq))
    n_times = length(times)
    n_freqs = length(freqs_vec)
    
    # Check channel exists
    if !(channel in channel_labels(tf_plot))
        error("Channel $channel not found. Available: $(channel_labels(tf_plot))")
    end
    
    # Reshape to freq × time matrix
    power_mat = zeros(n_freqs, n_times)
    for ti in 1:n_times
        for fi in 1:n_freqs
            row_idx = (ti - 1) * n_freqs + fi
            power_mat[fi, ti] = tf_plot.data[row_idx, channel]
        end
    end
    
    # Create figure
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1],
              xlabel="Time (s)",
              ylabel="Frequency (Hz)",
              title=isnothing(title) ? "$(tf_data.condition_name) - $channel" : title)
    
    # Determine color range
    if isnothing(colorrange)
        cmin, cmax = extrema(filter(!isnan, power_mat))
        colorrange = (cmin, cmax)
    end
    
    # Plot heatmap
    hm = heatmap!(ax, collect(times), collect(freqs_vec), power_mat, 
                  colormap=colormap, colorrange=colorrange)
    
    if colorbar
        label = isnothing(baseline_window) ? "Power" : 
                (baseline_method == :db ? "Power (dB)" : "Power")
        Colorbar(fig[1, 2], hm, label=label)
    end
    
    display(fig)
    return fig, ax
end

# Convenience: plot first channel if not specified
function plot_time_frequency(tf_data::TimeFreqData; kwargs...)
    ch = channel_labels(tf_data)[1]
    return plot_time_frequency(tf_data, ch; kwargs...)
end
