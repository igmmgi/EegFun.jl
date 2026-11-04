# Power spectrum plotting functions for channels, SpectrumData, and ICA components.

const PLOT_POWER_SPECTRUM_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Spectral analysis parameters
    :window_size => (1024, "Size of the FFT window for spectral estimation"),
    :overlap => (0.5, "Overlap between windows for Welch's method (0.0 to 1.0)"),
    :max_freq => (200.0, "Maximum frequency to display in Hz"),
    :window_function => (DSP.hanning, "Window function for spectral estimation"),

    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    :legend => (true, "Whether to show the legend"),
    :show_freq_bands => (true, "Whether to show frequency band indicators"),

    # Axis styling
    :xlabel => ("Frequency (Hz)", "X-axis label"),
    :ylabel => ("Power Spectral Density (μV²/Hz)", "Y-axis label"),
    :plot_title => ("Power Spectrum", "Plot title"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Scale parameters
    :xscale => (:linear, "X-axis scale: :linear or :log10"),
    :yscale => (:log10, "Y-axis scale: :linear or :log10"),
    :unit => (:linear, "Power unit: :linear (μV²/Hz) or :dB (decibels)"),

    # Font sizes
    :plot_title_fontsize => (16, "Font size for title"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :legend_fontsize => (12, "Font size for legend"),

    # Line styling
    :linewidth => (2, "Line width for spectrum lines"),
    :line_alpha => (0.8, "Transparency for spectrum lines"),

    # Frequency band styling
    :freq_band_alpha => (0.3, "Transparency for frequency band indicators"),
    :freq_band_height => (0.1, "Height of frequency band indicators"),

    # Grid styling
    :xgrid => (true, "Whether to show x-axis grid"),
    :ygrid => (true, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
)

const _EEG_FREQ_BANDS = ("δ" => (0.5, 4.0), "θ" => (4.0, 8.0), "α" => (8.0, 13.0), "β" => (13.0, 30.0), "γ" => (30.0, 100.0))
const _EEG_BAND_COLORS = (:lightblue, :lightgreen, :yellow, :orange, :lightcoral)

function _plot_power_spectrum!(
    fig,
    ax,
    signals::AbstractVector{<:AbstractVector},
    channels_to_plot::Vector{Symbol},
    fs::Real;
    line_labels::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
    kwargs...,
)
    isempty(signals) && return

    kwargs_norm = _normalize_plot_kwargs(kwargs)
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs_norm)

    n_samples = length(first(signals))
    effective_window_size = max(4, min(plot_kwargs[:window_size], n_samples))
    noverlap = round(Int, effective_window_size * plot_kwargs[:overlap])
    n_ch = length(channels_to_plot)

    first_sig = first(signals)
    T = float(eltype(first_sig))
    sig_list = all(s -> eltype(s) === T, signals) ? signals : [Float64.(s) for s in signals]
    target_sig = first(sig_list)

    cfg = DSP.WelchConfig(
        target_sig;
        n = effective_window_size,
        noverlap = noverlap,
        fs = Float64(fs),
        window = plot_kwargs[:window_function],
    )

    pgram_first = DSP.welch_pgram(target_sig, cfg)
    freq_data = fill(DSP.freq(pgram_first), n_ch)

    first_power = DSP.power(pgram_first)
    power_raw_data = Vector{typeof(first_power)}(undef, n_ch)
    power_raw_data[1] = first_power
    for i = 2:n_ch
        power_raw_data[i] = DSP.power(DSP.welch_pgram(sig_list[i], cfg))
    end

    _plot_power_spectrum_data!(fig, ax, freq_data, power_raw_data, channels_to_plot; line_labels = line_labels, kwargs_norm...)
end

function _plot_power_spectrum!(fig, ax, df::DataFrame, channels_to_plot::Vector{Symbol}, fs::Real; kwargs...)
    signals = [df[!, ch] for ch in channels_to_plot]
    _plot_power_spectrum!(fig, ax, signals, channels_to_plot, fs; kwargs...)
end

function _calc_power_extrema(power_data::Vector{<:AbstractVector})
    min_val = Inf
    max_val = -Inf
    pos_min = Inf
    @inbounds for vec in power_data
        for x in vec
            x < min_val && (min_val = x)
            x > max_val && (max_val = x)
            (x > 0 && x < pos_min) && (pos_min = x)
        end
    end
    return (min_val, max_val, pos_min)
end

function _apply_spectrum_yscale!(ax::Axis, y_scale::Symbol, unit::Symbol, min_val::Real, max_val::Real, pos_min::Real)
    if y_scale == :log10 && unit != :dB
        log_min, log_max = isinf(pos_min) ? (0.001, 1.0) : (max(0.001, pos_min), max(max(0.001, pos_min) * 10, max_val))
        ylims!(ax, (log_min, log_max))
        ax.yscale = log10
    else
        ax.yscale = identity
        ylims!(ax, unit == :dB ? (min_val - 5, max_val + 5) : (0.0, max(0.1, max_val)))
    end
    return nothing
end

function _plot_power_spectrum_data!(
    fig,
    ax,
    freq_data::Vector{<:AbstractVector},
    power_raw_data::Vector{<:AbstractVector},
    channels_to_plot::Vector{Symbol};
    line_labels::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, _normalize_plot_kwargs(kwargs))

    max_freq = plot_kwargs[:max_freq]
    x_scale = plot_kwargs[:xscale]
    y_scale = plot_kwargs[:yscale]
    unit = plot_kwargs[:unit]

    if !(x_scale in (:linear, :log10))
        @minimal_warning "Invalid xscale '$x_scale'. Using :linear instead."
        x_scale = :linear
    end
    if !(y_scale in (:linear, :log10))
        @minimal_warning "Invalid yscale '$y_scale'. Using :linear instead."
        y_scale = :linear
    end
    if !(unit in (:linear, :dB))
        @minimal_warning "Invalid unit '$unit'. Using :linear instead."
        unit = :linear
    end

    ax.xlabel = plot_kwargs[:xlabel]
    ax.ylabel = unit == :dB ? "Power Spectral Density (dB)" : plot_kwargs[:ylabel]
    ax.title = plot_kwargs[:show_title] && !isnothing(plot_kwargs[:plot_title]) ? plot_kwargs[:plot_title] : ""
    ax.titlesize = plot_kwargs[:plot_title_fontsize]
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    controls_area = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Relative(0.85))
    colsize!(fig.layout, 2, Relative(0.15))
    rowsize!(fig.layout, 1, Relative(0.95))

    x_options = ["Linear (Hz)", "Log10 (Hz)"]
    x_default_idx = x_scale == :log10 ? 2 : 1
    y_options = ["Linear (μV²/Hz)", "Log10 (μV²/Hz)", "Decibels (dB)"]
    y_default_idx = unit == :dB ? 3 : (y_scale == :log10 ? 2 : 1)

    x_scale_obs = Observable(x_scale)
    y_mode_obs = Observable(y_options[y_default_idx])
    band_ax_ref = Ref{Union{Nothing,Axis}}(nothing)

    Label(controls_area[1, 1], "X-Axis:")
    x_menu = Menu(controls_area[1, 2], options = x_options, default = x_options[x_default_idx])
    Label(controls_area[2, 1], "Y-Axis:")
    y_menu = Menu(controls_area[2, 2], options = y_options, default = y_options[y_default_idx])

    on(x_menu.selection) do sel
        x_scale_obs[] = sel == "Log10 (Hz)" ? :log10 : :linear
    end
    on(y_menu.selection) do sel
        y_mode_obs[] = sel
    end

    if x_scale == :log10
        xlims!(ax, (0.1, max_freq))
        ax.xscale = log10
    else
        xlims!(ax, (0, max_freq))
    end

    get_power_data(u) = u == :dB ? [10.0 .* log10.(max.(psd_raw, 1e-10)) for psd_raw in power_raw_data] : power_raw_data
    initial_power_data = get_power_data(unit)

    n_ch = length(channels_to_plot)
    psd_observables = Vector{Observable}(undef, n_ch)
    lw = plot_kwargs[:linewidth]
    lalpha = plot_kwargs[:line_alpha]

    for (i, ch) in enumerate(channels_to_plot)
        psd_obs = Observable(initial_power_data[i])
        psd_observables[i] = psd_obs
        lbl = !isnothing(line_labels) && i <= length(line_labels) ? String(line_labels[i]) : string(ch)
        lines!(ax, freq_data[i], psd_obs, label = lbl, linewidth = lw, alpha = lalpha)
    end

    init_min, init_max, init_pos_min = _calc_power_extrema(initial_power_data)
    _apply_spectrum_yscale!(ax, y_scale, unit, init_min, init_max, init_pos_min)

    if plot_kwargs[:show_freq_bands]
        band_ax = Axis(fig[2, 1], height = 50)
        band_ax_ref[] = band_ax
        linkxaxes!(ax, band_ax)
        xlims!(band_ax, 0.1, max_freq)
        ylims!(band_ax, 0, 1)
        hidedecorations!(band_ax)
        hidespines!(band_ax)

        band_alpha = plot_kwargs[:freq_band_alpha]
        for (i, (band_name, (fmin, fmax))) in enumerate(_EEG_FREQ_BANDS)
            fmax <= max_freq || continue
            lines!(band_ax, [fmin, fmax], [0.5, 0.5], color = _EEG_BAND_COLORS[i], linewidth = 8, alpha = band_alpha)
            text!(band_ax, (fmin + fmax) / 2, 0.5, text = band_name, align = (:center, :center), fontsize = 32, color = :black)
        end
    end

    on(x_scale_obs) do scale
        xlims!(ax, (0.1, max_freq))
        xscale_fn = scale == :log10 ? log10 : identity
        ax.xscale = xscale_fn
        !isnothing(band_ax_ref[]) && (band_ax_ref[].xscale = xscale_fn)
    end

    on(y_mode_obs) do mode
        new_unit = mode == "Decibels (dB)" ? :dB : :linear
        new_y_scale = mode == "Log10 (μV²/Hz)" ? :log10 : :linear
        converted = get_power_data(new_unit)
        for i in eachindex(psd_observables)
            psd_observables[i][] = converted[i]
        end
        ax.ylabel = new_unit == :dB ? "Power Spectral Density (dB)" : plot_kwargs[:ylabel]
        pmin, pmax, ppos_min = _calc_power_extrema(converted)
        _apply_spectrum_yscale!(ax, new_y_scale, new_unit, pmin, pmax, ppos_min)
    end
end

"""
    plot_channel_spectrum(dat; kwargs...)

Plot power spectrum for specified channels from EEG data with interactive controls.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024)
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5)
- `max_freq::Real`: Maximum frequency to display in Hz (default: 200.0)
- `x_scale::Symbol`: Initial scale for x-axis, one of :linear or :log10 (default: :linear)
- `y_scale::Symbol`: Initial scale for y-axis, one of :linear or :log10 (default: :linear)
- `unit::Symbol`: Power unit, one of :linear (μV²/Hz) or :dB (decibels) (default: :linear)
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning)
- `show_freq_bands::Bool`: Whether to show frequency band indicators (default: true)
- `show_legend::Bool`: Whether to show the legend (default: true)
- `display_plot::Bool`: Whether to display the plot (default: true)

# Features
- Interactive checkboxes to toggle between linear and log scales for both axes
- Automatic multi-column legend for many channels
- Frequency band indicators (δ, θ, α, β, γ) below the main plot
- Responsive grid layout with controls in the top-right corner

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot
- `ax::Axis`: The axis containing the plot

# Examples
```julia
# Basic usage with interactive controls
fig, ax = plot_channel_spectrum(dat)

# With specific channels
fig, ax = plot_channel_spectrum(dat, 
    channel_selection = channels([:Fp1, :Fp2]))

# With log scales
fig, ax = plot_channel_spectrum(dat, 
    x_scale = :log10, 
    y_scale = :log10)

# With dB units (EEGLAB-style)
fig, ax = plot_channel_spectrum(dat, unit = :dB)
```
"""
function plot_channel_spectrum(
    dat::SingleDataFrameEeg;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    kwargs...,
)
    _set_window_title(_generate_window_title(dat))

    kwargs_norm = _normalize_plot_kwargs(kwargs)
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs_norm)

    chan_fn = _to_channel_func(channel_selection)
    channels_to_plot = get_selected_channels(dat.layout, chan_fn)

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    n_rows = nrow(dat.data)
    all_rows = combined_sel isa _AllSamples
    selected_samples = all_rows ? (1:n_rows) : get_selected_samples(dat, combined_sel)

    signals = if all_rows || selected_samples == 1:n_rows
        [dat.data[!, ch] for ch in channels_to_plot]
    else
        [dat.data[!, ch][selected_samples] for ch in channels_to_plot]
    end

    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])

    _plot_power_spectrum!(fig, ax, signals, channels_to_plot, sample_rate(dat); kwargs_norm...)

    if plot_kwargs[:legend] && !isempty(channels_to_plot)
        n_channels = length(channels_to_plot)
        axislegend(ax, nbanks = n_channels > 10 ? cld(n_channels, 20) : 1)
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end

function _fill_normalized_ica_row!(
    dat_matrix::Matrix{Float64},
    row_idx::Int,
    col::AbstractVector,
    selected_samples::AbstractVector{Int},
    mean_c::Float64,
    scale_c::Float64,
)
    @inbounds @simd for j in eachindex(selected_samples)
        dat_matrix[row_idx, j] = (Float64(col[selected_samples[j]]) - mean_c) / scale_c
    end
    return nothing
end

"""
    plot_ica_component_spectrum(dat::ContinuousData, ica_result::InfoIca; component_selection::Function = components(), kwargs...)

Plot power spectrum of ICA component(s) with interactive controls for axis scaling.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `ica_result::InfoIca`: The ICA result object.
- `component_selection::Function`: Function that returns boolean vector for component filtering (default: all components).

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Initial scale for x-axis, one of :linear or :log10 (default: :linear).
- `y_scale::Symbol`: Initial scale for y-axis, one of :linear or :log10 (default: :linear).
- `unit::Symbol`: Power unit, one of :linear (μV²/Hz) or :dB (decibels) (default: :linear).
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning).
- `show_freq_bands::Bool`: Whether to show frequency band indicators (default: true).
- `show_legend::Bool`: Whether to show the legend with component variance percentages (default: true).
- `display_plot::Bool`: Whether to display the plot (default: true).

# Returns
- Named tuple `(fig, axes)` containing the Makie `Figure` and `[ax]`.

# Examples
```julia
plot_ica_component_spectrum(dat, ica_result)
plot_ica_component_spectrum(dat, ica_result, component_selection = components([1, 3, 5]))
plot_ica_component_spectrum(dat, ica_result, component_selection = components(1), unit = :dB)
```
"""
function plot_ica_component_spectrum(
    dat::ContinuousData,
    ica_result::InfoIca;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    component_selection = components(),
    kwargs...,
)
    kwargs_norm = _normalize_plot_kwargs(kwargs)
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs_norm)

    selected_components = get_selected_components(ica_result, component_selection)
    if isempty(selected_components)
        @minimal_warning "No components selected for ICA spectrum plot"
        return (fig = _eegfun_figure(), axes = Axis[])
    end

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    n_rows = nrow(dat.data)
    selected_samples = combined_sel isa _AllSamples ? (1:n_rows) : get_selected_samples(dat, combined_sel)

    ica_labels = ica_result.layout.data[!, :label]
    n_ica_ch = length(ica_labels)
    dat_matrix = Matrix{Float64}(undef, n_ica_ch, length(selected_samples))
    scale_val = ica_result.scale
    for (ci, ch) in enumerate(ica_labels)
        _fill_normalized_ica_row!(dat_matrix, ci, dat.data[!, ch], selected_samples, ica_result.mean[ci], scale_val)
    end

    proj = (@view ica_result.unmixing[selected_components, :]) * dat_matrix

    n_sel = length(selected_components)
    components_to_plot = Vector{Symbol}(undef, n_sel)
    component_signals = Vector{Vector{Float64}}(undef, n_sel)
    component_labels = Vector{String}(undef, n_sel)
    for (i, comp_idx) in enumerate(selected_components)
        components_to_plot[i] = Symbol("IC_", comp_idx)
        component_signals[i] = proj[i, :]
        component_labels[i] = @sprintf("IC %d (%.1f%%)", comp_idx, ica_result.variance[comp_idx] * 100)
    end

    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])

    _plot_power_spectrum!(fig, ax, component_signals, components_to_plot, sample_rate(dat); line_labels = component_labels, kwargs_norm...)

    if plot_kwargs[:legend]
        ncols = n_sel > 10 ? cld(n_sel, 10) : 1
        axislegend(ax, nbanks = ncols)
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax])
end

"""
    plot_channel_spectrum(dat::SpectrumData; kwargs...)

Plot a pre-computed power spectrum with interactive controls.
"""
function plot_channel_spectrum(dat::SpectrumData; channel_selection = channels(), kwargs...)
    _set_window_title(_generate_window_title(dat))

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    all_names = propertynames(dat.data)
    valid_channels = [ch for ch in selected_channels if ch != :freq && ch in all_names]
    isempty(valid_channels) && @minimal_error "No valid channels found in spectrum data."

    shared_freq = Float64.(dat.data[!, :freq])
    freq_data = fill(shared_freq, length(valid_channels))
    power_raw_data = Vector{Float64}[Float64.(dat.data[!, ch]) for ch in valid_channels]

    kwargs_norm = _normalize_plot_kwargs(kwargs)
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs_norm)

    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])

    _plot_power_spectrum_data!(fig, ax, freq_data, power_raw_data, valid_channels; kwargs_norm...)

    if plot_kwargs[:legend]
        n_channels = length(valid_channels)
        axislegend(ax, nbanks = n_channels > 10 ? cld(n_channels, 20) : 1)
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end
