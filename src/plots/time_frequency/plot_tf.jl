# Time-frequency plotting functions for visualizing TimeFreqData.

const _SIGNED_TF_BASELINE_METHODS = (:db, :percent, :relchange, :normchange, :vssum, :zscore, :absolute)

const PLOT_TF_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display
    :display_plot => (true, "Whether to display the plot"),
    :theme_fontsize => (nothing, "Font size for theme"),

    # Layout
    :layout => (:single, "Layout type: :single, :grid, or :topo"),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),

    # Colormap and color range
    :colormap => (nothing, "Colormap for the heatmap"),
    :colorrange => (nothing, "Color range as (min, max) tuple. If nothing, automatically determined from data"),
    # Colorbar parameters
    :colorbar_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Colorbar block (see Makie's [Colorbar documentation](https://docs.makie.org/stable/reference/blocks/colorbar/) for available attributes).",
    ),

    # Specific colorbar overrides for tf
    :colorbar_plot => (true, "Whether to show the colorbar"),
    :colorbar_position => (:right, "Position of the colorbar (:right, :left, :top, :bottom, or tuple)"),
    :colorbar_label => (nothing, "Custom colorbar label. If nothing, automatically determined"),

    # Axis
    :xscale => (:linear, "X-axis scale (:linear or :log10)"),
    :yscale => (:linear, "Frequency (y) axis scale (:linear or :log10; legacy alias `ylogscale` also accepted)"),
    :xlim => (nothing, "X-axis time limits as (min, max) tuple in seconds. If nothing, shows all time points"),
    :interpolate => (false, "Whether to interpolate the heatmap for a smoother appearance"),
    :plot_title => (nothing, "Plot title. If nothing, automatically determined from condition name and channel"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :xticks => (nothing, "Custom x-axis ticks (e.g., -0.2:0.2:1.0)"),
    :yticks => (nothing, "Custom y-axis ticks (e.g., [2, 10, 20, 40, 80])"),
    :time_unit => (
        :s,
        "Time unit for x-axis display (:s or :ms). Only affects axis labels and tick formatting — all intervals remain in seconds.",
    ),

    # Baseline
    :baseline_method => (:db, "Baseline correction method: :db, :absolute, :relative, :relchange, :percent, :zscore"),

    # Layout parameters
    :layout_kwargs => (
        NamedTuple(),
        "Additional parameters to configure grid or topographical layouts (e.g., grid_dims, topo_plot_width, grid_rowgap)",
    ),
)

# === SHARED TF EXTRACTION & RANGE HELPERS ===

function _sorted_unique_f64(col::AbstractVector)::Vector{Float64}
    n = length(col)
    n == 0 && return Float64[]
    if issorted(col)
        cnt = 1
        prev = Float64(col[1])
        @inbounds for i = 2:n
            v = Float64(col[i])
            if v != prev
                cnt += 1
                prev = v
            end
        end
        out = Vector{Float64}(undef, cnt)
        prev = Float64(col[1])
        out[1] = prev
        k = 1
        @inbounds for i = 2:n
            v = Float64(col[i])
            if v != prev
                k += 1
                out[k] = v
                prev = v
            end
        end
        return out
    end
    u = unique(col)
    out = Vector{Float64}(undef, length(u))
    @inbounds for i in eachindex(u)
        out[i] = Float64(u[i])
    end
    sort!(out)
    return out
end

struct _TFFastExtractor
    df::DataFrame
    freqs::Vector{Float64}
    times::Vector{Float64}
    canonical::Bool
    row_index::Union{Nothing,Dict{Tuple{Float64,Float64},Int}}
    rounded_freqs::Vector{Float64}
    rounded_times::Vector{Float64}
end

function _is_canonical_tf_grid(freq_col::AbstractVector, time_col::AbstractVector, freqs::Vector{Float64}, times::Vector{Float64})::Bool
    nf = length(freqs)
    nt = length(times)
    (length(freq_col) == nf * nt && length(time_col) == nf * nt) || return false
    idx = 0
    @inbounds for ti = 1:nt
        t = times[ti]
        for fi = 1:nf
            idx += 1
            (abs(Float64(freq_col[idx]) - freqs[fi]) <= 1e-6 && abs(Float64(time_col[idx]) - t) <= 1e-6) || return false
        end
    end
    return true
end

function _tf_make_extractor(df::DataFrame, freqs::Vector{Float64}, times::Vector{Float64})::_TFFastExtractor
    if _is_canonical_tf_grid(df[!, :freq], df[!, :time], freqs, times)
        return _TFFastExtractor(df, freqs, times, true, nothing, Float64[], Float64[])
    end
    return _TFFastExtractor(df, freqs, times, false, _tf_build_row_index(df), round.(freqs, digits = 6), round.(times, digits = 6))
end

function _tf_extract_matrix(ext::_TFFastExtractor, channel::Symbol)::Matrix{Float64}
    nf = length(ext.freqs)
    nt = length(ext.times)
    if ext.canonical
        col = ext.df[!, channel]
        mat = Matrix{Float64}(undef, nf, nt)
        if col isa Vector{Float64}
            copyto!(mat, col)
        else
            @inbounds for i in eachindex(mat, col)
                mat[i] = Float64(col[i])
            end
        end
        return mat
    end
    return _tf_df_to_matrix(
        ext.df,
        channel,
        ext.freqs,
        ext.times,
        ext.row_index::Dict{Tuple{Float64,Float64},Int};
        rounded_freqs = ext.rounded_freqs,
        rounded_times = ext.rounded_times,
    )
end

function _select_tf_channels(all_channels::Vector{Symbol}, channel_selection, channel_plot_order = nothing)
    chan_func = _to_channel_func(channel_selection)
    selected_channels = all_channels[chan_func(all_channels)]
    isempty(selected_channels) && @minimal_error("No channels matched. Available: $(all_channels)")
    if !isnothing(channel_plot_order)
        ordered = Symbol[Symbol(ch) for ch in channel_plot_order if Symbol(ch) in selected_channels]
        !isempty(ordered) && (selected_channels = ordered)
    end
    return selected_channels
end

function _tf_finite_extrema(mats)
    gmin = Inf
    gmax = -Inf
    for mat in mats
        @inbounds for v in mat
            if isfinite(v)
                gmin = min(gmin, v)
                gmax = max(gmax, v)
            end
        end
    end
    return gmin, gmax
end

function _tf_auto_colorrange(gmin::Float64, gmax::Float64, colormap; force_symmetric::Bool = false, has_signed_baseline::Bool = false)
    isinf(gmin) && return (-1.0, 1.0)
    if force_symmetric || (_is_diverging_colormap(colormap) && (gmin < 0 || has_signed_baseline))
        return _symmetric_colorrange(gmin, gmax)
    elseif gmin ≈ gmax
        return (gmin - 1.0, gmax + 1.0)
    else
        return (gmin, gmax)
    end
end

"""
    plot_tf(filepath::String; kwargs...)
    plot_tf(tf::TimeFreqData; channel_selection, baseline_interval, kwargs...)
    plot_tf(tfs::Vector{TimeFreqData}; channel_selection, baseline_interval, kwargs...)

Plot time-frequency data as a heatmap (time × frequency). Supports single condition,
multi-condition grid, and topographic layouts.

# Arguments
- `filepath::String`: Path to a `.jld2` file, or a `TimeFreqData` / `Vector{TimeFreqData}` object
- `channel_selection::Function`: Channel filter (default: all channels)
- `baseline_interval`: Baseline window as `(start, stop)` in seconds (default: `nothing`)

$(_generate_kwargs_doc(PLOT_TF_KWARGS))

# Examples
```julia
plot_tf("tf_morlet_result.jld2")
plot_tf("tf_morlet_result.jld2"; layout=:grid, channel_selection=channels([:Cz, :Pz]))
plot_tf(tf; colormap=:RdBu, ylogscale=true, baseline_interval=(-0.3, 0.0))
```
"""
function plot_tf(filepath::String; kwargs...)
    data = read_data(filepath)
    isnothing(data) && @minimal_error "No data found in file: $filepath"
    return plot_tf(data; kwargs...)
end

function plot_tf(tf_data::AbstractVector{<:TimeFreqData}; channel_selection = channels(), baseline_interval::Interval = nothing, kwargs...)
    isempty(tf_data) && @minimal_error("Empty TimeFreqData vector")

    plot_kwargs = _merge_plot_kwargs(PLOT_TF_KWARGS, _normalize_plot_kwargs(kwargs))
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))
    layout_kwargs = pop!(plot_kwargs, :layout_kwargs, (;))

    colorrange = plot_kwargs[:colorrange]
    colorbar = pop!(plot_kwargs, :colorbar_plot, true)
    xlogscale = plot_kwargs[:xscale] in (:log10, :log, log10)
    ylogscale = plot_kwargs[:yscale] in (:log10, :log, log10)
    show_title = plot_kwargs[:show_title]
    xlim = plot_kwargs[:xlim]
    xticks = plot_kwargs[:xticks]
    yticks = plot_kwargs[:yticks]
    interpolate = plot_kwargs[:interpolate]
    baseline_method = plot_kwargs[:baseline_method]
    time_unit = plot_kwargs[:time_unit]
    grid_dims = get(layout_kwargs, :grid_dims, nothing)

    n = length(tf_data)
    rows, cols = isnothing(grid_dims) ? _best_rect(n) : grid_dims
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig =
        _eegfun_figure(; size = (max(600, cols * 400), max(400, rows * 350)), figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    colormap = _resolve_theme_colormap(fig, plot_kwargs[:colormap])
    _set_window_title("$(basename(first(tf_data).file)) — Time-Frequency ($(n) conditions)")

    tf_plots = map(tf_data) do tf
        (!isnothing(baseline_interval) && isnothing(tf.baseline)) ? tf_baseline(tf, baseline_interval; method = baseline_method) : tf
    end

    channel = first(_select_tf_channels(channel_labels(first(tf_plots)), channel_selection))
    first_df = first(tf_plots).data_power
    all_times = _sorted_unique_f64(first_df[!, :time])
    freqs_vec = _sorted_unique_f64(first_df[!, :freq])

    power_mats = Vector{Matrix{Float64}}(undef, n)
    for (i, tf) in enumerate(tf_plots)
        power_mats[i] = _tf_extract_matrix(_tf_make_extractor(tf.data_power, freqs_vec, all_times), channel)
    end

    has_signed_baseline =
        any(tf -> (!isnothing(tf.baseline) && tf.baseline.method in _SIGNED_TF_BASELINE_METHODS), tf_plots) ||
        (!isnothing(baseline_interval) && baseline_method in _SIGNED_TF_BASELINE_METHODS)

    shared_colorrange = if isnothing(colorrange)
        gmin, gmax = _tf_finite_extrema(power_mats)
        _tf_auto_colorrange(gmin, gmax, colormap; has_signed_baseline = has_signed_baseline)
    else
        colorrange
    end

    axes = Vector{Axis}(undef, n)
    last_hm = nothing
    for (idx, tf) in enumerate(tf_plots)
        row = fld(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1
        ax = Axis(fig[row, col])
        axes[idx] = ax

        last_hm = _plot_tf_heatmap!(
            ax,
            tf,
            channel,
            freqs_vec,
            all_times;
            power_mat = power_mats[idx],
            colormap = colormap,
            colorrange = shared_colorrange,
            xlogscale = xlogscale,
            ylogscale = ylogscale,
            xlim = xlim,
            interpolate = interpolate,
            xticks = xticks,
            yticks = yticks,
            time_unit = time_unit,
        )
        ax.title = show_title ? (isnothing(plot_kwargs[:plot_title]) ? "$(tf.condition_name) — $channel" : plot_kwargs[:plot_title]) : ""

        col != 1 && (ax.ylabel = ""; ax.yticklabelsvisible = false)
        row != rows && (ax.xlabel = ""; ax.xticklabelsvisible = false)
    end

    length(axes) > 1 && linkaxes!(axes...)

    if colorbar && !isnothing(last_hm)
        cb_label =
            isnothing(plot_kwargs[:colorbar_label]) ? _tf_colorbar_label(first(tf_plots), baseline_interval, baseline_method) :
            plot_kwargs[:colorbar_label]
        cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:rows, 1:cols)
        Colorbar(fig[cb_pos...], last_hm; label = cb_label, colorbar_kwargs...)
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

function plot_tf(
    tf_data::TimeFreqData;
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    baseline_interval::Interval = nothing,
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_TF_KWARGS, _normalize_plot_kwargs(kwargs))
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))
    layout_kw_in = pop!(plot_kwargs, :layout_kwargs, (;))

    layout = plot_kwargs[:layout]
    colorrange = plot_kwargs[:colorrange]
    title = plot_kwargs[:plot_title]
    show_title = plot_kwargs[:show_title]
    colorbar = pop!(plot_kwargs, :colorbar_plot, true)
    xlogscale = plot_kwargs[:xscale] in (:log10, :log, log10)
    ylogscale = plot_kwargs[:yscale] in (:log10, :log, log10)
    xlim = plot_kwargs[:xlim]
    xticks = plot_kwargs[:xticks]
    yticks = plot_kwargs[:yticks]
    interpolate = plot_kwargs[:interpolate]
    baseline_method = plot_kwargs[:baseline_method]
    time_unit = plot_kwargs[:time_unit]
    layout_grid_dims = get(layout_kw_in, :grid_dims, nothing)

    if !isnothing(baseline_interval) && !isnothing(tf_data.baseline)
        @minimal_warning "Data has already been baselined (method: $(tf_data.baseline.method), window: $(tf_data.baseline.window)). " *
                         "Ignoring baseline_interval parameter. Use the data as-is or create a new TimeFreqData without baseline."
        tf_plot = tf_data
    elseif !isnothing(baseline_interval)
        tf_plot = tf_baseline(tf_data, baseline_interval; method = baseline_method)
    else
        tf_plot = tf_data
    end

    eff_baseline_method =
        !isnothing(tf_plot.baseline) ? tf_plot.baseline.method : (!isnothing(baseline_interval) ? baseline_method : nothing)
    has_signed_baseline = eff_baseline_method in _SIGNED_TF_BASELINE_METHODS

    selected_channels = _select_tf_channels(channel_labels(tf_plot), channel_selection, channel_plot_order)
    plot_channels = layout === :single ? [first(selected_channels)] : selected_channels

    df_power = tf_plot.data_power
    all_times = _sorted_unique_f64(df_power[!, :time])
    freqs_vec = _sorted_unique_f64(df_power[!, :freq])

    layout_kwargs = Dict{Symbol,Any}(
        :grid_dims => layout_grid_dims,
        :grid_rowgap => get(layout_kw_in, :grid_rowgap, 10),
        :grid_colgap => get(layout_kw_in, :grid_colgap, 10),
        :topo_plot_width => get(layout_kw_in, :topo_plot_width, 0.1),
        :topo_plot_height => get(layout_kw_in, :topo_plot_height, 0.1),
    )

    fig_size = if layout === :single
        (800, 500)
    elseif layout === :topo
        (1200, 1000)
    else
        rows, cols = !isnothing(layout_grid_dims) ? layout_grid_dims : _best_rect(length(plot_channels))
        (max(600, cols * 350), max(400, rows * 300))
    end

    _set_window_title("$(basename(tf_data.file)) — $(tf_data.condition_name) — Time-Frequency")
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; size = fig_size, figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    colormap = _resolve_theme_colormap(fig, plot_kwargs[:colormap])

    eeg_layout = hasproperty(tf_data, :layout) ? tf_data.layout : nothing
    plot_layout = create_layout(layout, plot_channels, eeg_layout; layout_kwargs...)
    axes, layout_channels = _apply_layout!(fig, plot_layout; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)

    ext = _tf_make_extractor(df_power, freqs_vec, all_times)
    cached_mats = Dict{Symbol,Matrix{Float64}}()
    sizehint!(cached_mats, length(plot_channels))
    for ch in plot_channels
        cached_mats[ch] = _tf_extract_matrix(ext, ch)
    end

    shared_colorrange = if isnothing(colorrange)
        gmin, gmax = _tf_finite_extrema(values(cached_mats))
        _tf_auto_colorrange(gmin, gmax, colormap; has_signed_baseline = has_signed_baseline)
    else
        colorrange
    end

    last_hm = nothing
    for (ax, channel) in zip(axes, layout_channels)
        pmat = get!(() -> _tf_extract_matrix(ext, channel), cached_mats, channel)
        last_hm = _plot_tf_heatmap!(
            ax,
            tf_plot,
            channel,
            freqs_vec,
            all_times;
            power_mat = pmat,
            colormap = colormap,
            colorrange = shared_colorrange,
            xlogscale = xlogscale,
            ylogscale = ylogscale,
            xlim = xlim,
            interpolate = interpolate,
            xticks = xticks,
            yticks = yticks,
            time_unit = time_unit,
        )

        ax.title = if !show_title
            ""
        elseif layout === :single
            isnothing(title) ? "$(tf_data.condition_name) - $channel" : title
        else
            isnothing(title) ? string(channel) : title
        end
    end

    if plot_layout.type == :grid
        rows, cols = plot_layout.dims
        for (idx, ax) in enumerate(axes)
            row = fld(idx - 1, cols) + 1
            col = mod(idx - 1, cols) + 1
            col != 1 && (ax.ylabel = ""; ax.yticklabelsvisible = false)
            row != rows && (ax.xlabel = ""; ax.xticklabelsvisible = false)
        end
    elseif plot_layout.type == :topo
        for ax in axes
            hidedecorations!(ax, grid = false, ticks = true, ticklabels = true)
            hidespines!(ax)
        end
    end

    length(axes) > 1 && linkaxes!(axes...)

    if colorbar && !isnothing(last_hm)
        cb_label =
            isnothing(plot_kwargs[:colorbar_label]) ? _tf_colorbar_label(tf_plot, baseline_interval, baseline_method) :
            plot_kwargs[:colorbar_label]

        if layout === :single
            cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
            Colorbar(fig[cb_pos...], last_hm; label = cb_label, colorbar_kwargs...)
        elseif layout === :grid
            rows, cols = plot_layout.dims
            cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:rows, 1:cols)
            Colorbar(fig[cb_pos...], last_hm; label = cb_label, colorbar_kwargs...)
        else
            cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
            Colorbar(fig[cb_pos[1], cb_pos[2]+1], last_hm; label = cb_label, height = Relative(0.5), colorbar_kwargs...)
        end
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

# === INTERNAL HELPERS ===

function _plot_tf_heatmap!(
    ax::Axis,
    tf_plot::TimeFreqData,
    channel::Symbol,
    freqs_vec::Vector{Float64},
    times::Vector{Float64};
    power_mat::Union{Nothing,Matrix{Float64}} = nothing,
    row_index::Union{Nothing,Dict{Tuple{Float64,Float64},Int}} = nothing,
    colormap = nothing,
    colorrange::Tuple{Real,Real} = (-1.0, 1.0),
    xlogscale::Bool = false,
    ylogscale::Bool = false,
    xlim::Union{Nothing,Tuple{Real,Real}} = nothing,
    interpolate::Bool = false,
    xticks = nothing,
    yticks = nothing,
    time_unit::Symbol = :s,
)
    mat = if !isnothing(power_mat)
        power_mat
    elseif !isnothing(row_index)
        _tf_df_to_matrix(tf_plot.data_power, channel, freqs_vec, times, row_index)
    else
        _tf_extract_matrix(_tf_make_extractor(tf_plot.data_power, freqs_vec, times), channel)
    end

    if time_unit == :ms
        ax.xlabel = "Time (ms)"
        ax.xtickformat = values -> [string(round(Int, v * 1000)) for v in values]
    else
        ax.xlabel = "Time (s)"
    end
    ax.ylabel = "Frequency (Hz)"

    if xlogscale
        pos_times = filter(>(0), times)
        if !isempty(pos_times)
            xlims!(ax, (minimum(pos_times), maximum(pos_times)))
            ax.xscale = log10
        end
    end

    if ylogscale
        ylims!(ax, (minimum(freqs_vec), maximum(freqs_vec)))
        ax.yscale = log10
        ax.ytickformat = values -> [string(round(Int, v)) for v in values]
    end

    !isnothing(xlim) && xlims!(ax, xlim)
    !isnothing(xticks) && (ax.xticks = xticks)
    !isnothing(yticks) && (ax.yticks = yticks)

    return heatmap!(
        ax,
        times,
        freqs_vec,
        transpose(mat),
        colormap = _resolve_theme_colormap(ax, colormap),
        colorrange = colorrange,
        nan_color = :transparent,
        interpolate = interpolate,
    )
end

function _tf_baseline_unit_label(method::Symbol)
    method == :db && return "Power (dB)"
    method == :percent && return "Power (% change)"
    method == :relchange && return "Power (relative)"
    return "Power"
end

function _tf_colorbar_label(tf_plot::TimeFreqData, baseline_interval, baseline_method::Symbol)
    !isnothing(tf_plot.baseline) && return _tf_baseline_unit_label(tf_plot.baseline.method)
    !isnothing(baseline_interval) && return _tf_baseline_unit_label(baseline_method)
    return "Power"
end
