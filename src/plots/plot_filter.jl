# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_FILTER_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Axis limits
    :ylimit => ((-100, 5), "Y-axis limits in dB as (min, max) tuple"),
    :xlimit => (nothing, "X-axis limits in Hz as (min, max) tuple. If nothing, automatically determined"),
    
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    
    # Plot styling
    :title => ("Filter Frequency Response", "Plot title"),
    :xlabel => ("Frequency (Hz)", "X-axis label"),
    :ylabel => ("Magnitude (dB)", "Y-axis label"),
    
    # Font sizes
    :title_fontsize => (24, "Font size for title"),
    :label_fontsize => (22, "Font size for axis labels"),
    :tick_fontsize => (20, "Font size for tick labels"),
    :legend_fontsize => (32, "Font size for legend"),
    
    # Line styling
    :actual_linewidth => (4, "Line width for actual response"),
    :ideal_linewidth => (3, "Line width for ideal response"),
    :actual_color => (:black, "Color for actual response"),
    :ideal_color => (:green, "Color for ideal response"),
    :ideal_linestyle => (:dash, "Line style for ideal response"),
    
    # Reference lines
    :reference_lines => ([-3, -6], "Reference lines in dB to display"),
    :reference_color => (:gray, "Color for reference lines"),
    :reference_linestyle => (:dash, "Line style for reference lines"),
    
    # Transition region styling
    :transition_alpha => (0.2, "Transparency for transition region shading"),
    :transition_color => (:gray, "Color for transition region shading"),
    :transition_line_color => (:gray, "Color for transition region lines"),
    :transition_line_style => (:dash, "Line style for transition region lines"),
    
    # Plot parameters
    :n_points => (2000, "Number of frequency points for response calculation"),
    :xscale => (Makie.Symlog10(10.0), "X-axis scale type"),
)

"""
    plot_filter_response(filter, sample_rate::Real, filter_freq::Real, transition_band::Real; kwargs...)

Plot the frequency response of a digital filter with ideal response overlay.

# Arguments
- `filter`: A digital filter object (FIR coefficients or DSP.jl filter)
- `sample_rate`: Sampling rate in Hz
- `filter_freq`: Cutoff frequency in Hz
- `transition_band`: Width of transition band in Hz

$(generate_kwargs_doc(PLOT_FILTER_KWARGS))

# Returns
- `fig`: Makie Figure object
- `ax`: Makie Axis object

# Example
```julia
# Basic usage
fig, ax = plot_filter_response(filter, 1000, 50, 10)

# Custom styling
fig, ax = plot_filter_response(filter, 1000, 50, 10;
    title = "Custom Filter Response",
    actual_color = :blue,
    ideal_color = :red,
    reference_lines = [-3, -6, -12],
    display_plot = false)
```
"""
function plot_filter_response(
    filter,
    sample_rate::Real,
    filter_freq::Real,
    transition_band::Real;
    kwargs...
)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_FILTER_KWARGS, kwargs)
    
    # Extract commonly used values
    ylimit = plot_kwargs[:ylimit]
    xlimit = plot_kwargs[:xlimit]
    n_points = plot_kwargs[:n_points]
    # Determine x-axis limits based on filter type and cutoff
    if isnothing(xlimit)
        if filter_freq < 2
            xlimit = (0, filter_freq * 10)
        else
            xlimit = (0, sample_rate / 2)
        end
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
        title = plot_kwargs[:title],
        titlesize = plot_kwargs[:title_fontsize],
        xlabelsize = plot_kwargs[:label_fontsize],
        ylabelsize = plot_kwargs[:label_fontsize],
        xticklabelsize = plot_kwargs[:tick_fontsize],
        yticklabelsize = plot_kwargs[:tick_fontsize],
        xscale = plot_kwargs[:xscale],
        limits = (xlimit, ylimit),
    )

    # Simple logarithmic frequency spacing
    freqs = exp10.(range(log10(0.01), log10(sample_rate/2), length = n_points))
    freqs = [0.0; freqs]  # Add DC point

    w = 2π * freqs / sample_rate

    # Get frequency response
    if filter isa Vector  # FIR filter coefficients
        n = 0:(length(filter)-1)
        resp = [sum(filter .* exp.(-im * w_k * n)) for w_k in w]
    else  # Other filter types
        resp = freqresp(filter, w)
    end

    mag_db = 20 * log10.(abs.(resp))

    # Get filter type from response
    start_response = mag_db[2]  # Use second point to avoid DC issues
    end_response = mag_db[end]

    # Determine filter type
    filter_type = "hp"
    if start_response > -3 && end_response < -20
        filter_type = "lp"
    end

    # Calculate actual stopband attenuation
    if filter_type == "lp"
        stopband_mask = freqs .>= (filter_freq + transition_band)
    elseif filter_type == "hp"
        stopband_mask = freqs .<= (filter_freq - transition_band)
    end

    stopband_db = mag_db[stopband_mask]
    stopband_db = replace(stopband_db, -Inf => -100.0)
    actual_stopband_atten = mean(stopband_db)
    actual_stopband_linear = 10^(actual_stopband_atten/20)

    # Calculate ideal response
    ideal_response = zeros(length(freqs))
    for (i, f) in enumerate(freqs)
        if filter_type == "lp"
            if f <= filter_freq
                ideal_response[i] = 1.0
            elseif f >= filter_freq + transition_band
                ideal_response[i] = actual_stopband_linear
            else
                ideal_response[i] =
                    1.0 * (filter_freq + transition_band - f) / transition_band +
                    actual_stopband_linear * (f - filter_freq) / transition_band
            end
        elseif filter_type == "hp"
            if f >= filter_freq
                ideal_response[i] = 1.0
            elseif f <= filter_freq - transition_band
                ideal_response[i] = actual_stopband_linear
            else
                ideal_response[i] =
                    1.0 * (f - (filter_freq - transition_band)) / transition_band +
                    actual_stopband_linear * (filter_freq - f) / transition_band
            end
        end
    end

    # Define transition regions first
    if filter_type == "lp"
        f_s = filter_freq + transition_band  # Add for lowpass
        transition_regions = [(filter_freq, f_s)]
    elseif filter_type == "hp"
        f_s = filter_freq - transition_band  # Subtract for highpass
        transition_regions = [(f_s, filter_freq)]
    end

    for (start_f, end_f) in transition_regions
        vspan!(ax, start_f, end_f, color = (plot_kwargs[:transition_color], plot_kwargs[:transition_alpha]))
        vlines!(ax, [start_f, end_f], color = plot_kwargs[:transition_line_color], linestyle = plot_kwargs[:transition_line_style])
    end

    # Plot responses
    lines!(ax, freqs, mag_db, label = "Actual", color = plot_kwargs[:actual_color], linewidth = plot_kwargs[:actual_linewidth])
    ideal_mag_db = 20 * log10.(ideal_response)
    lines!(ax, freqs, ideal_mag_db, color = plot_kwargs[:ideal_color], linestyle = plot_kwargs[:ideal_linestyle], label = "Ideal", linewidth = plot_kwargs[:ideal_linewidth])

    # Add reference lines
    hlines!(ax, plot_kwargs[:reference_lines], color = plot_kwargs[:reference_color], linestyle = plot_kwargs[:reference_linestyle])
    for (i, ref_val) in enumerate(plot_kwargs[:reference_lines])
        text!(ax, sample_rate/2, ref_val, text = "$ref_val dB", align = (:right, :center), fontsize = plot_kwargs[:label_fontsize])
    end

    # X-axis ticks based on xlimits
    if xlimit[2] <= 5  # For low frequency (typically highpass) plots
        xticks = [0, 0.1, 0.2, 0.5, 1, 2, filter_freq, f_s, 5]
    else  # For full range plots (typically lowpass)
        xticks = [0, 1, 2, 5, 10, 20, 50, filter_freq, f_s, 100, sample_rate/2]
    end

    # Filter out ticks beyond xlimit
    xticks = [x for x in xticks if x >= xlimit[1] && x <= xlimit[2]]
    ax.xticks = (xticks, string.(round.(xticks, digits = 1)))

    legend_position = filter_type == "lp" ? :lb : :rb
    axislegend(; position = legend_position, labelsize = plot_kwargs[:legend_fontsize])
    
    if plot_kwargs[:display_plot]
        display(fig)
    end

    return fig, ax
end
