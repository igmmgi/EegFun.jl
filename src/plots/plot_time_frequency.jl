"""
Time-frequency plotting functions for visualizing TimeFreqData.
"""

"""
    plot_time_frequency(tf_data::TimeFreqData, channel::Symbol;
                        baseline_window=nothing, baseline_method=:db,
                        colormap=:viridis, colorrange=nothing,
                        title=nothing, colorbar=true, ylogscale=false)

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
- `ylogscale`: Use logarithmic scale for y-axis (frequencies) (default: false)

# Returns
- `(Figure, Axis)`: Makie figure and axis

# Example
```julia
fig, ax = plot_time_frequency(tf_data, :Cz; baseline_window=(-0.3, 0.0), ylogscale=true)
```
"""
function plot_time_frequency(tf_data::TimeFreqData, channel::Symbol;
                             baseline_window::Union{Nothing,Tuple{Real,Real}}=nothing,
                             baseline_method::Symbol=:db,
                             colormap=:viridis,
                             colorrange::Union{Nothing,Tuple{Real,Real}}=nothing,
                             title::Union{Nothing,String}=nothing,
                             colorbar::Bool=true,
                             ylogscale::Bool=false)
    
    # Apply baseline if requested, but only if data hasn't already been baselined
    if !isnothing(baseline_window) && tf_data.baseline !== nothing
        @warn "Data has already been baselined (method: $(tf_data.baseline.method), window: $(tf_data.baseline.window)). " *
              "Ignoring baseline_window parameter. Use the data as-is or create a new TimeFreqData without baseline."
        tf_plot = tf_data
    elseif !isnothing(baseline_window)
        # Apply baseline on-the-fly
        tf_plot = tf_baseline(tf_data, baseline_window; method=baseline_method)
    else
        tf_plot = tf_data
    end
    
    # Check channel exists
    if !(channel in channel_labels(tf_plot))
        error("Channel $channel not found. Available: $(channel_labels(tf_plot))")
    end
    
    # Get unique times and frequencies (sorted)
    times = sort(unique(tf_plot.data_power.time))
    freqs_vec = sort(unique(tf_plot.data_power.freq))
    n_times = length(times)
    n_freqs = length(freqs_vec)
    
    # Always use lookup approach to be robust to DataFrame structure
    # Create a lookup dictionary: (time, freq) -> power value
    power_dict = Dict{Tuple{Float64,Float64},Float64}()
    for row in eachrow(tf_plot.data_power)
        key = (row.time, row.freq)
        if !haskey(power_dict, key)  # Take first match if duplicates exist
            power_dict[key] = row[channel]
        end
    end
    
    # Build power matrix: rows = frequencies, columns = times
    # power_mat[fi, ti] means: row fi (frequency), column ti (time)
    power_mat = zeros(n_freqs, n_times)
    for (ti, t) in enumerate(times)
        for (fi, f) in enumerate(freqs_vec)
            key = (t, f)
            if haskey(power_dict, key)
                power_mat[fi, ti] = power_dict[key]
            else
                power_mat[fi, ti] = NaN
            end
        end
    end
    
    # Transpose for Makie heatmap (expects n_times × n_freqs)
    # After transpose: rows = times, columns = frequencies
    power_mat = power_mat'
    
    # Create figure
    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1],
              xlabel="Time (s)",
              ylabel="Frequency (Hz)",
              title=isnothing(title) ? "$(tf_data.condition_name) - $channel" : title,
              yscale=ylogscale ? log10 : identity)
    
    # Disable scientific notation for y-axis ticks when using log scale
    if ylogscale
        ax.ytickformat = values -> [string(round(Int, v)) for v in values]
    end
    
    # Determine color range
    if isnothing(colorrange)
        cmin, cmax = extrema(Base.filter(!isnan, power_mat))
        colorrange = (cmin, cmax)
    end
    
    # Plot heatmap
    # Makie heatmap!(ax, x, y, data) expects data to be (length(y), length(x))
    # x = times (n_times), y = freqs_vec (n_freqs)
    # So data should be (n_freqs, n_times)
    # After transpose, power_mat is (n_times, n_freqs), which is what Makie expects
    # Set NaN color to transparent so edge-filtered regions are not displayed
    hm = heatmap!(ax, times, freqs_vec, power_mat, 
                  colormap=colormap, colorrange=colorrange,
                  nan_color=:transparent)
    
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
        elseif !isnothing(baseline_window)
            # Baseline was just applied via parameter
            label = baseline_method == :db ? "Power (dB)" : 
                   baseline_method == :percent ? "Power (% change)" :
                   baseline_method == :relchange ? "Power (relative)" : "Power"
        else
            label = "Power"
        end
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
