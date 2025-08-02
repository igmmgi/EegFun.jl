"""
    plot_filter_response(filter, sample_rate::Real, filter_freq::Real, transition_band::Real; 
                        ylimit::Tuple=(-100, 5), xlimit::Union{Nothing,Tuple{Real,Real}}=nothing)

Plot the frequency response of a digital filter with ideal response overlay.

Arguments:
- `filter`: A digital filter object (FIR coefficients or DSP.jl filter)
- `sample_rate`: Sampling rate in Hz
- `filter_freq`: Cutoff frequency in Hz
- `transition_band`: Width of transition band in Hz
- `ylimit`: Y-axis limits in dB as (min, max) tuple (default: (-100, 5))
- `xlimit`: Optional X-axis limits in Hz as (min, max) tuple. If nothing, automatically determined.

Returns:
- `fig`: Makie Figure object
- `ax`: Makie Axis object
"""
function plot_filter_response(
    filter,
    sample_rate::Real,
    filter_freq::Real,
    transition_band::Real;
    ylimit::Tuple = (-100, 5),
    xlimit::Union{Nothing,Tuple{Real,Real}} = nothing,
)
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
        xlabel = "Frequency (Hz)",
        ylabel = "Magnitude (dB)",
        title = "Filter Frequency Response",
        titlesize = 24,
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xscale = Makie.Symlog10(10.0),
        limits = (xlimit, ylimit),
    )

    # Simple logarithmic frequency spacing
    n_points = 2000
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
        vspan!(ax, start_f, end_f, color = (:gray, 0.2))
        vlines!(ax, [start_f, end_f], color = :gray, linestyle = :dash)
    end

    # Plot responses
    lines!(ax, freqs, mag_db, label = "Actual", color = :black, linewidth = 4)
    ideal_mag_db = 20 * log10.(ideal_response)
    lines!(ax, freqs, ideal_mag_db, color = :green, linestyle = :dash, label = "Ideal", linewidth = 3)

    # Add reference lines
    hlines!(ax, [-3, -6], color = :gray, linestyle = :dash)
    text!(ax, sample_rate/2, -3, text = "-3 dB", align = (:right, :center), fontsize = 22)
    text!(ax, sample_rate/2, -6, text = "-6 dB", align = (:right, :center), fontsize = 22)

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
    axislegend(; position = legend_position, labelsize = 32)
    display(fig)

    return fig, ax
end
