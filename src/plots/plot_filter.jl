# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_FILTER_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Axis limits
    :ylim => ((-100, 5), "Y-axis limits in dB as (min, max) tuple"),
    :xlim => (nothing, "X-axis limits in Hz as (min, max) tuple. If nothing, automatically determined"),

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

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

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

# =============================================================================
# FILTER VISUALIZATION
# =============================================================================

"""
    plot_filter_response(filter_info::FilterInfo; kwargs...)

Plot the frequency response of a digital filter with ideal response overlay.

# Arguments
- `filter_info::FilterInfo`: Filter information struct
- `kwargs...`: Additional keyword arguments

$(generate_kwargs_doc(PLOT_FILTER_KWARGS))

# Returns
- `fig`: Makie Figure object
- `ax`: Tuple of Makie Axis objects (linear, dB, impulse)
"""
function plot_filter_response(
    filter_info::FilterInfo;
    filter_func::Function = filtfilt,  # Default to zero-phase filtering
    kwargs...,
)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_FILTER_KWARGS, kwargs)

    # Determine x-axis limits based on filter type and cutoff
    xlim = plot_kwargs[:xlim]
    if isnothing(xlim)
        if filter_info.cutoff_freq < 2
            xlim = (0, filter_info.cutoff_freq * 10)
        else
            xlim = (0, filter_info.sample_rate / 2)
        end
    end

    fig = Figure()

    # Base axis properties
    base_props = (
        xlabel = plot_kwargs[:xlabel],
        title = plot_kwargs[:title],
        titlesize = plot_kwargs[:title_fontsize],
        xlabelsize = plot_kwargs[:label_fontsize],
        ylabelsize = plot_kwargs[:label_fontsize],
        xticklabelsize = plot_kwargs[:tick_fontsize],
        yticklabelsize = plot_kwargs[:tick_fontsize],
    )

    # Convert xscale symbol to function if needed
    xscale = plot_kwargs[:xscale]
    if xscale == :linear
        xscale = identity
    elseif xscale == :log
        xscale = log10
        # Adjust x-limits for log scale (can't include 0)
        if xlim[1] == 0
            xlim = (0.01, xlim[2])  # Use small positive value instead of 0
        end
    end

    # Add xscale from kwargs
    xscale_props = (xscale = xscale,)

    # Create three axes in a row
    ax1 = Axis(fig[1, 1]; base_props..., xscale_props..., ylabel = "Magnitude (linear)", limits = (xlim, (0, 1.1)))
    ax2 = Axis(
        fig[1, 2];
        base_props...,
        xscale_props...,
        ylabel = plot_kwargs[:ylabel],
        limits = (xlim, plot_kwargs[:ylim]),
    )
    ax3 = Axis(
        fig[1, 3];
        xlabel = "Time (samples)",
        ylabel = "Amplitude",
        title = "Impulse response",
        titlesize = plot_kwargs[:title_fontsize],
        xlabelsize = plot_kwargs[:label_fontsize],
        ylabelsize = plot_kwargs[:label_fontsize],
        xticklabelsize = plot_kwargs[:tick_fontsize],
        yticklabelsize = plot_kwargs[:tick_fontsize],
    )

    # Set grid using the shared function
    _set_axis_grid!(
        ax1;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )
    _set_axis_grid!(
        ax2;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )
    _set_axis_grid!(
        ax3;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    # Simple logarithmic frequency spacing
    n_points = plot_kwargs[:n_points]
    freqs = exp10.(range(log10(0.01), log10(filter_info.sample_rate/2), length = n_points))
    freqs = [0.0; freqs]  # Add DC point

    w = 2π * freqs / filter_info.sample_rate

    # Get frequency response
    if filter_info.filter_object isa Vector  # FIR filter coefficients
        n = 0:(length(filter_info.filter_object)-1)
        resp = [sum(filter_info.filter_object .* exp.(-im * w_k * n)) for w_k in w]
    else  # Other filter types
        resp = freqresp(filter_info.filter_object, w)
    end

    # Calculate magnitude responses
    mag_linear = abs.(resp)
    mag_db = 20 * log10.(mag_linear)

    # Plot actual responses
    lines!(
        ax1,
        freqs,
        mag_linear,
        label = "Actual",
        color = plot_kwargs[:actual_color],
        linewidth = plot_kwargs[:actual_linewidth],
    )
    lines!(
        ax2,
        freqs,
        mag_db,
        label = "Actual",
        color = plot_kwargs[:actual_color],
        linewidth = plot_kwargs[:actual_linewidth],
    )

    # Add vertical line at cutoff frequency to both subplots
    vlines!(ax1, [filter_info.cutoff_freq], color = :red, linestyle = :dash, linewidth = 2)
    vlines!(ax2, [filter_info.cutoff_freq], color = :red, linestyle = :dash, linewidth = 2)

    # Add reference lines
    hlines!(
        ax1,
        [0.707],
        color = plot_kwargs[:reference_color],
        linestyle = plot_kwargs[:reference_linestyle],
        alpha = 0.5,
    )  # -3 dB point
    hlines!(
        ax2,
        plot_kwargs[:reference_lines],
        color = plot_kwargs[:reference_color],
        linestyle = plot_kwargs[:reference_linestyle],
        alpha = 0.5,
    )

    # Calculate and plot impulse response
    # Create a unit impulse with samples before and after
    n_samples_before = 200  # Samples before impulse
    n_samples_after = 200  # Samples after impulse
    n_total = n_samples_before + n_samples_after + 1
    impulse = zeros(n_total)
    impulse[n_samples_before+1] = 1.0  # Unit impulse at t=0

    # Apply the specified filter function to get impulse response
    if filter_info.filter_method == "fir"
        impulse_response = filter_info.filter_object
        time_samples = (-(length(impulse_response)÷2)):((length(impulse_response)-1)÷2)
    else  # IIR filter
        impulse_response = filter_func(filter_info.filter_object, impulse)
        time_samples = (-n_samples_before):n_samples_after
    end

    # Plot impulse response
    lines!(
        ax3,
        time_samples,
        impulse_response,
        color = plot_kwargs[:actual_color],
        linewidth = plot_kwargs[:actual_linewidth],
    )

    # Add zero line for reference
    hlines!(ax3, [0], color = plot_kwargs[:reference_color], linestyle = plot_kwargs[:reference_linestyle], alpha = 0.5)

    plot_kwargs[:display_plot] && display_figure(fig)

    return fig, (ax1, ax2, ax3)
end
