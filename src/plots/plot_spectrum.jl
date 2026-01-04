"""
Power spectrum plotting functions for visualizing SpectrumData.
"""

"""
    plot_freq_spectrum(spectrum_data::SpectrumData;
                       channel_selection::Function=channels(),
                       x_scale::Symbol=:linear,
                       y_scale::Symbol=:linear,
                       unit::Symbol=:linear,
                       max_freq::Union{Nothing,Real}=nothing,
                       colormap=:default,
                       title::Union{Nothing,String}=nothing,
                       show_legend::Bool=true,
                       linewidth::Real=2,
                       line_alpha::Real=0.8)

Plot power spectrum data for selected channels.

# Arguments
- `spectrum_data::SpectrumData`: Power spectrum data

# Keyword Arguments
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Example: `channel_selection=channels([:Cz, :Pz])` for multiple channels
- `x_scale::Symbol=:linear`: X-axis scale (`:linear` or `:log10`)
- `y_scale::Symbol=:linear`: Y-axis scale (`:linear` or `:log10`)
- `unit::Symbol=:linear`: Power unit (`:linear` for μV²/Hz or `:dB` for decibels)
- `max_freq::Union{Nothing,Real}=nothing`: Maximum frequency to display in Hz. If `nothing`, uses all frequencies.
- `colormap`: Colormap for lines (default: automatic colors)
- `title::Union{Nothing,String}=nothing`: Plot title (default: auto-generated)
- `show_legend::Bool=true`: Show legend with channel names
- `linewidth::Real=2`: Line width for spectrum lines
- `line_alpha::Real=0.8`: Transparency for spectrum lines

# Returns
- `(Figure, Axis)`: Makie figure and axis

# Example
```julia
# Plot all channels
fig, ax = plot_freq_spectrum(spectrum_data)

# Single channel with log scales
fig, ax = plot_freq_spectrum(spectrum_data; 
    channel_selection=channels(:Cz), 
    x_scale=:log10, 
    y_scale=:log10)

# With dB units
fig, ax = plot_freq_spectrum(spectrum_data; unit=:dB, max_freq=100.0)
```
"""
function plot_freq_spectrum(
    spectrum_data::SpectrumData;
    channel_selection::Function = channels(),
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    unit::Symbol = :linear,
    max_freq::Union{Nothing,Real} = nothing,
    colormap = :default,
    title::Union{Nothing,String} = nothing,
    show_legend::Bool = true,
    linewidth::Real = 2,
    line_alpha::Real = 0.8,
)
    # Get selected channels
    selected_channels = get_selected_channels(spectrum_data, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && error("No channels selected. Available channels: $(channel_labels(spectrum_data))")

    # Filter to only channels that actually exist in the spectrum data
    # Use propertynames to get column names as Symbols
    all_names = propertynames(spectrum_data.data)
    available_channels = [ch for ch in all_names if ch != :freq]
    valid_channels = [ch for ch in selected_channels if ch in available_channels]
    
    if isempty(valid_channels)
        error("No valid channels found in spectrum data. Selected channels: $(selected_channels), Available channels: $(available_channels)")
    end
    
    # Warn about channels that were selected but don't exist
    missing_channels = setdiff(selected_channels, valid_channels)
    if !isempty(missing_channels)
        @warn "Channels $(missing_channels) not found in spectrum data. Available channels: $(available_channels)"
    end

    # Validate scales
    valid_scales = [:linear, :log10]
    if !(x_scale in valid_scales)
        @warn "Invalid x_scale '$x_scale'. Using :linear instead."
        x_scale = :linear
    end
    if !(y_scale in valid_scales)
        @warn "Invalid y_scale '$y_scale'. Using :linear instead."
        y_scale = :linear
    end

    # Validate unit
    valid_units = [:linear, :dB]
    if !(unit in valid_units)
        @warn "Invalid unit '$unit'. Using :linear instead."
        unit = :linear
    end

    # Get frequency and power data
    freqs = spectrum_data.data.freq
    power_data = Dict{Symbol,Vector{Float64}}()

    # Extract power for each valid channel
    for channel in valid_channels
        power_raw = Float64.(spectrum_data.data[!, channel])

        # Convert to requested unit
        if unit == :dB
            # Convert to dB: 10 * log10(power / reference)
            # Reference is 1 μV²/Hz (standard in EEG)
            power = 10.0 .* log10.(max.(power_raw, 1e-10))  # Avoid log(0) with small threshold
        else
            power = power_raw
        end

        power_data[channel] = power
    end

    # Filter by max_freq if specified
    if !isnothing(max_freq)
        mask = freqs .<= max_freq
        freqs = freqs[mask]
        for channel in keys(power_data)
            power_data[channel] = power_data[channel][mask]
        end
    end

    # Create figure
    fig = Figure(size = (800, 500))
    ax = Axis(
        fig[1, 1],
        xlabel = "Frequency (Hz)",
        ylabel = unit == :dB ? "Power Spectral Density (dB)" : "Power Spectral Density (μV²/Hz)",
        title = isnothing(title) ? "$(spectrum_data.condition_name) - Power Spectrum" : title,
        xscale = x_scale == :log10 ? log10 : identity,
        yscale = y_scale == :log10 ? log10 : identity,
    )

    # Plot each valid channel
    for channel in valid_channels
        lines!(ax, freqs, power_data[channel], label = string(channel), linewidth = linewidth, alpha = line_alpha)
    end

    # Add legend if requested and we have plots
    if show_legend && !isempty(valid_channels)
        n_channels = length(valid_channels)
        n_cols = n_channels > 10 ? cld(n_channels, 20) : 1
        axislegend(ax, nbanks = n_cols)
    end

    # Set appropriate limits
    if x_scale == :log10
        xlims!(ax, (max(0.1, minimum(freqs)), maximum(freqs)))
    else
        xlims!(ax, (0, maximum(freqs)))
    end

    if y_scale == :log10
        all_powers = vcat([power_data[ch] for ch in keys(power_data)]...)
        positive_powers = Base.filter(x -> x > 0, all_powers)
        if !isempty(positive_powers)
            log_min = max(0.001, minimum(positive_powers))
            log_max = maximum(positive_powers)
            ylims!(ax, (log_min, log_max))
        end
    else
        all_powers = vcat([power_data[ch] for ch in keys(power_data)]...)
        ylims!(ax, (0.0, max(0.1, maximum(all_powers))))
    end

    display(fig)
    return fig, ax
end

