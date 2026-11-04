# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_FILTER_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :ylim => ((-100, 5), "Y-axis limits in dB as (min, max) tuple"),
    :xlim => (nothing, "X-axis limits in Hz as (min, max) tuple. If nothing, automatically determined"),
    :display_plot => (true, "Whether to display the plot"),
    :plot_title => ("Filter Frequency Response", "Plot title"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :xlabel => ("Frequency (Hz)", "X-axis label"),
    :ylabel => ("Magnitude (dB)", "Y-axis label"),
    :plot_title_fontsize => (24, "Font size for title"),
    :label_fontsize => (22, "Font size for axis labels"),
    :tick_fontsize => (20, "Font size for tick labels"),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    :linewidth => (4, "Line width for response curves"),
    :color => (nothing, "Color for response curves"),
    :reference_lines => ([-3, -6], "Reference lines in dB to display"),
    :reference_color => (:gray, "Color for reference lines"),
    :reference_linestyle => (:dash, "Line style for reference lines"),
    :n_points => (2000, "Number of frequency points for response calculation"),
    :xscale => (Makie.Symlog10(10.0), "X-axis scale type"),
)

function _build_freq_grid(n_points::Int, fs::Float64)
    freqs = Vector{Float64}(undef, n_points + 1)
    w = Vector{Float64}(undef, n_points + 1)
    freqs[1] = 0.0
    w[1] = 0.0
    log_grid = range(log10(0.01), log10(fs / 2); length = n_points)
    w_scale = 2π / fs
    @inbounds for i = 1:n_points
        f = exp10(log_grid[i])
        freqs[i+1] = f
        w[i+1] = f * w_scale
    end
    return freqs, w
end

function _fir_mag_db(b::Vector{Float64}, w::Vector{Float64}, is_filtfilt::Bool)::Vector{Float64}
    n_w = length(w)
    n_taps = length(b)
    mag_db = Vector{Float64}(undef, n_w)
    @inbounds for i = 1:n_w
        w_k = w[i]
        acc_re = 0.0
        acc_im = 0.0
        for k = 1:n_taps
            c = cis(-w_k * (k - 1))
            bk = b[k]
            acc_re = muladd(bk, real(c), acc_re)
            acc_im = muladd(bk, imag(c), acc_im)
        end
        m = hypot(acc_re, acc_im)
        is_filtfilt && (m *= m)
        mag_db[i] = 20.0 * log10(m)
    end
    return mag_db
end

function _complex_to_db(resp::AbstractVector, is_filtfilt::Bool)::Vector{Float64}
    mag_db = Vector{Float64}(undef, length(resp))
    @inbounds for i in eachindex(resp)
        m = abs(resp[i])
        is_filtfilt && (m *= m)
        mag_db[i] = 20.0 * log10(m)
    end
    return mag_db
end

"""
    plot_filter_response(filter_info::FilterInfo; filter_func::Function = filtfilt, kwargs...) -> NamedTuple

Plot the frequency response (dB) and impulse response of a digital `FilterInfo`.

$(_generate_kwargs_doc(PLOT_FILTER_KWARGS))
"""
function plot_filter_response(filter_info::FilterInfo; filter_func::Function = filtfilt, kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_FILTER_KWARGS, kwargs)

    xlim = plot_kwargs[:xlim]
    if isnothing(xlim)
        xlim = filter_info.cutoff_freq < 2 ? (0.0, filter_info.cutoff_freq * 10) : (0.0, filter_info.sample_rate / 2)
    end

    xscale = plot_kwargs[:xscale]
    if xscale === :linear
        xscale = identity
    elseif xscale === :log || xscale === :log10
        xscale = log10
        xlim[1] == 0 && (xlim = (0.01, xlim[2]))
    end

    _set_window_title("Filter Response")
    fig = _eegfun_figure(size = (1200, 600))
    show_title = plot_kwargs[:show_title]

    font_props = (
        titlesize = plot_kwargs[:plot_title_fontsize],
        xlabelsize = plot_kwargs[:label_fontsize],
        ylabelsize = plot_kwargs[:label_fontsize],
        xticklabelsize = plot_kwargs[:tick_fontsize],
        yticklabelsize = plot_kwargs[:tick_fontsize],
    )

    ax_db = Axis(
        fig[1, 1];
        font_props...,
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
        title = show_title && !isnothing(plot_kwargs[:plot_title]) ? plot_kwargs[:plot_title] : "",
        xscale = xscale,
        limits = (xlim, plot_kwargs[:ylim]),
    )
    ax_imp = Axis(fig[1, 2]; font_props..., xlabel = "Time (samples)", ylabel = "Amplitude", title = show_title ? "Impulse Response" : "")

    for a in (ax_db, ax_imp)
        _set_axis_grid!(
            a;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )
    end

    freqs, w = _build_freq_grid(Int(plot_kwargs[:n_points]), Float64(filter_info.sample_rate))

    is_filtfilt = filter_func === filtfilt
    fobj = filter_info.filter_object
    if fobj isa Vector{Float64}
        mag_db = _fir_mag_db(fobj, w, is_filtfilt)
    else
        mag_db = _complex_to_db(freqresp(fobj, w), is_filtfilt)
    end

    lines!(ax_db, freqs, mag_db, color = _resolve_theme_linecolor(ax_db, plot_kwargs[:color]), linewidth = plot_kwargs[:linewidth])
    vlines!(ax_db, [filter_info.cutoff_freq], color = :red, linestyle = :dash, linewidth = 2)
    hlines!(
        ax_db,
        plot_kwargs[:reference_lines],
        color = plot_kwargs[:reference_color],
        linestyle = plot_kwargs[:reference_linestyle],
        alpha = 0.5,
    )

    is_fir = filter_info.filter_method === :fir || fobj isa Vector
    if is_fir && !is_filtfilt
        impulse_response = fobj
        time_samples = (-(length(impulse_response)÷2)):((length(impulse_response)-1)÷2)
    else
        n_half = is_fir ? max(200, length(fobj) * 2) : 200
        impulse = zeros(Float64, 2 * n_half + 1)
        impulse[n_half+1] = 1.0
        impulse_response = is_fir ? filtfilt(fobj, impulse) : filter_func(fobj, impulse)
        time_samples = (-n_half):n_half
    end

    lines!(
        ax_imp,
        time_samples,
        impulse_response,
        color = _resolve_theme_linecolor(ax_imp, plot_kwargs[:color]),
        linewidth = plot_kwargs[:linewidth],
    )
    hlines!(ax_imp, [0], color = plot_kwargs[:reference_color], linestyle = plot_kwargs[:reference_linestyle], alpha = 0.5)

    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax_db, ax_imp])
end
