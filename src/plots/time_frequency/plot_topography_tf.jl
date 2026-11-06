# Topographic visualizations for Time-Frequency data and TF statistical results.

function _sorted_interval_bounds(iv)
    raw = _normalize_interval(iv)
    return (min(Float64(raw[1]), Float64(raw[2])), max(Float64(raw[1]), Float64(raw[2])))
end

function _tf_select_window_values(
    all_freqs::Vector{Float64},
    all_times::Vector{Float64},
    freq_range::Interval,
    interval_selection::Interval,
)
    fr = _sorted_interval_bounds(freq_range)
    selected_freqs = [f for f in all_freqs if fr[1] <= f <= fr[2]]
    isempty(selected_freqs) &&
        @minimal_error("No frequencies found in range $(freq_range). Data range: $(first(all_freqs)) to $(last(all_freqs)) Hz")

    if isnothing(interval_selection) || interval_selection isa AllSelection
        selected_times = all_times
    else
        t_start = interval_selection isa TimeSelection ? interval_selection.start : interval_selection[1]
        t_stop = interval_selection isa TimeSelection ? interval_selection.stop : interval_selection[2]
        selected_times = [t for t in all_times if t_start <= t <= t_stop]
    end
    isempty(selected_times) && @minimal_error("No time points found in interval. Data range: $(first(all_times)) to $(last(all_times)) s")

    return fr, selected_freqs, selected_times
end

function _prepare_tf_topo_layout(orig_layout::Layout, common_channels::Vector{Symbol}; need_3d::Bool = true)::Layout
    orig_labels = orig_layout.data[!, :label]
    layout = if length(common_channels) < length(orig_labels) || common_channels != orig_labels
        subset_layout(orig_layout; channel_selection = channels(common_channels))
    elseif !has_2d_coords(orig_layout) || (need_3d && !has_3d_coords(orig_layout))
        copy(orig_layout)
    else
        orig_layout
    end
    has_valid_coordinates(layout) || @minimal_error("Cannot create topographic plot: layout has no spatial coordinates.")
    _ensure_coordinates_2d!(layout)
    need_3d && _ensure_coordinates_3d!(layout)
    return layout
end

function _pop_unused_tf_topo_kwargs!(plot_kwargs::Dict)
    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)
    return nothing
end

"""
    plot_topography(tf::TimeFreqData;
                     freq_range::Tuple{Real, Real},
                     interval_selection::Interval = times(),
                     baseline_interval = nothing,
                     baseline_method::Symbol = :db,
                     kwargs...)

Plot a topographic map of average power within a frequency band and time window.

# Arguments
- `tf::TimeFreqData`: Time-frequency data
- `freq_range::Tuple{Real, Real}`: Frequency range to average over (e.g., `(8.0, 12.0)` for alpha)
- `interval_selection::Interval`: Time window (default: all time points)
- `baseline_interval`: Baseline window as `(start, stop)` in seconds (e.g., `(-0.3, 0.0)`). Default: `nothing`
- `baseline_method::Symbol`: Baseline method (default: `:db`). Options: `:db`, `:absolute`, `:relative`, `:relchange`, `:percent`, `:zscore`

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- Named tuple `(fig, axes, ax)`

# Examples
```julia
# Alpha-band power (8-12 Hz) in 300-500ms window with dB baseline
plot_topography(tf, freq_range=(8.0, 12.0), interval_selection=times(0.3, 0.5),
                baseline_interval=(-0.3, 0.0))

# Theta-band power across full time window
plot_topography(tf, freq_range=(4.0, 7.0))
```
"""
function plot_topography(
    tf::TimeFreqData;
    freq_range::Interval,
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    baseline_method::Symbol = :db,
    kwargs...,
)
    if !isnothing(baseline_interval) && !isnothing(tf.baseline)
        @minimal_warning "Data has already been baselined (method: $(tf.baseline.method)). Ignoring baseline_interval parameter."
        tf_plot = tf
    elseif !isnothing(baseline_interval)
        tf_plot = tf_baseline(tf, baseline_interval; method = baseline_method)
    else
        tf_plot = tf
    end

    eff_baseline_method =
        !isnothing(tf_plot.baseline) ? tf_plot.baseline.method : (!isnothing(baseline_interval) ? baseline_method : nothing)
    has_signed_baseline = eff_baseline_method in _SIGNED_TF_BASELINE_METHODS

    orig_layout = tf_plot.layout
    has_valid_coordinates(orig_layout) || @minimal_error("Cannot create topographic plot: layout has no spatial coordinates.")

    df_power = tf_plot.data_power
    all_times = _sorted_unique_f64(df_power[!, :time])
    all_freqs = _sorted_unique_f64(df_power[!, :freq])
    fr, selected_freqs, selected_times = _tf_select_window_values(all_freqs, all_times, freq_range, interval_selection)

    row_indices = _tf_filter_row_indices(df_power, selected_freqs, selected_times)

    channels_list = channel_labels(tf_plot)
    common_channels = intersect(orig_layout.data[!, :label], channels_list)
    isempty(common_channels) && @minimal_error("None of the TF channels $(channels_list) were found in the layout.")
    layout = _prepare_tf_topo_layout(orig_layout, common_channels; need_3d = true)

    layout_labels = layout.data[!, :label]
    layout_values = _tf_extract_channel_means(df_power, layout_labels, row_indices)

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    colorrange = pop!(plot_kwargs, :colorrange)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    plot_title = pop!(plot_kwargs, :plot_title, "")
    show_title = pop!(plot_kwargs, :show_title, true)
    num_levels = pop!(plot_kwargs, :num_levels)
    highlight_channels = pop!(plot_kwargs, :highlight_channels, nothing)

    _pop_unused_tf_topo_kwargs!(plot_kwargs)
    pop!(plot_kwargs, :dims, nothing)

    freq_str = "$(round(Int, fr[1]))-$(round(Int, fr[2])) Hz"
    time_str =
        (isnothing(interval_selection) || interval_selection isa AllSelection) ? "all times" :
        @sprintf("%.0f–%.0f ms", selected_times[1] * 1000, selected_times[end] * 1000)
    default_title = "$(tf_plot.condition_name) — $freq_str, $time_str"
    user_title = !isnothing(plot_title) && !isempty(plot_title) ? plot_title : figure_title
    title = show_title ? (isempty(user_title) ? default_title : user_title) : ""

    fig = _eegfun_figure(size = (400, 400))
    actual_colormap = _resolve_theme_colormap(fig, colormap)

    if isnothing(colorrange)
        data_min, data_max = extrema(layout_values)
        colorrange = _tf_auto_colorrange(data_min, data_max, actual_colormap; has_signed_baseline = has_signed_baseline)
    end

    ax = Axis(fig[1, 1], aspect = DataAspect(), title = title, titlesize = plot_kwargs[:plot_title_fontsize])

    colorbar_kwargs = Dict{Symbol,Any}(pairs(pop!(plot_kwargs, :colorbar_kwargs, (;))))
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    co, _, _ = _render_topo_surface!(
        fig,
        ax,
        layout_values,
        layout;
        method = method,
        gridscale = gridscale,
        colormap = actual_colormap,
        colorrange = colorrange,
        num_levels = num_levels,
        plot_kwargs...,
    )

    !isnothing(highlight_channels) && _draw_highlight_channels!(ax, highlight_channels, layout)

    if colorbar_plot
        cb_label = _tf_colorbar_label(tf_plot, baseline_interval, baseline_method)
        Colorbar(fig[1, 2], co; colorbar_kwargs..., label = cb_label)
    end

    hidedecorations!(ax)
    display_plot && _display_figure(fig)
    return (fig = fig, axes = [ax], ax = ax)
end

function plot_topography(
    tfs::AbstractVector{<:TimeFreqData};
    freq_range::Interval,
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    baseline_method::Symbol = :db,
    kwargs...,
)
    isempty(tfs) && @minimal_error("Cannot plot empty vector of TimeFreqData")

    tf_plots = if !isnothing(baseline_interval)
        map(tfs) do tf
            if !isnothing(tf.baseline)
                @minimal_warning "$(tf.condition_name) already baselined, skipping."
                tf
            else
                tf_baseline(tf, baseline_interval; method = baseline_method)
            end
        end
    else
        tfs
    end

    n = length(tf_plots)
    n == 1 && return plot_topography(tf_plots[1]; freq_range = freq_range, interval_selection = interval_selection, kwargs...)

    orig_layout = tf_plots[1].layout
    has_valid_coordinates(orig_layout) || @minimal_error("Cannot create topographic plot: layout has no spatial coordinates.")

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    colorrange = pop!(plot_kwargs, :colorrange)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    plot_title = pop!(plot_kwargs, :plot_title, "")
    show_title = pop!(plot_kwargs, :show_title, true)
    num_levels = pop!(plot_kwargs, :num_levels)
    dims = pop!(plot_kwargs, :dims, nothing)
    highlight_channels = pop!(plot_kwargs, :highlight_channels, nothing)

    _pop_unused_tf_topo_kwargs!(plot_kwargs)
    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)

    first_df = tf_plots[1].data_power
    all_freqs = _sorted_unique_f64(first_df[!, :freq])
    all_times = _sorted_unique_f64(first_df[!, :time])
    fr, selected_freqs, selected_times = _tf_select_window_values(all_freqs, all_times, freq_range, interval_selection)

    common_channels = copy(orig_layout.data[!, :label])
    for tf in tf_plots
        common_channels = intersect(common_channels, channel_labels(tf))
    end
    isempty(common_channels) && @minimal_error("None of the TF channels were found in the layout.")
    layout = _prepare_tf_topo_layout(orig_layout, common_channels; need_3d = false)

    layout_labels = layout.data[!, :label]
    all_layout_values = Vector{Vector{Float64}}(undef, n)

    global_min, global_max = Inf, -Inf
    first_rows = _tf_filter_row_indices(first_df, selected_freqs, selected_times)
    for (idx, tf) in enumerate(tf_plots)
        df = tf.data_power
        row_indices =
            (idx == 1 || (df[!, :time] === first_df[!, :time] && df[!, :freq] === first_df[!, :freq])) ? first_rows :
            _tf_filter_row_indices(df, selected_freqs, selected_times)
        layout_values = _tf_extract_channel_means(df, layout_labels, row_indices)
        all_layout_values[idx] = layout_values
        lmin, lmax = extrema(layout_values)
        global_min = min(global_min, lmin)
        global_max = max(global_max, lmax)
    end

    n_rows, n_cols = isnothing(dims) ? _best_rect(n) : dims

    can_batch = n > 1 && _topo_method_supported(method)
    batched_interp = nothing
    xb_shared = nothing
    yb_shared = nothing
    if can_batch
        val_matrix = Matrix{Float64}(undef, length(layout_labels), n)
        for j = 1:n
            val_matrix[:, j] = all_layout_values[j]
        end
        batched_interp, xb_shared, yb_shared = _topo_interpolate(val_matrix, layout, gridscale; method = method)
    end

    freq_str = "$(round(Int, fr[1]))-$(round(Int, fr[2])) Hz"
    fig = _eegfun_figure(size = (250 * n_cols + 100, 250 * n_rows + 50))
    actual_colormap = _resolve_theme_colormap(fig, colormap)

    has_signed_baseline =
        any(tf -> (!isnothing(tf.baseline) && tf.baseline.method in _SIGNED_TF_BASELINE_METHODS), tf_plots) ||
        (!isnothing(baseline_interval) && baseline_method in _SIGNED_TF_BASELINE_METHODS)

    if isnothing(colorrange)
        colorrange = _tf_auto_colorrange(global_min, global_max, actual_colormap; has_signed_baseline = has_signed_baseline)
    end

    header_title = !isnothing(plot_title) && !isempty(plot_title) ? plot_title : figure_title
    (show_title && !isempty(header_title)) && Label(fig[0, 1:n_cols], header_title, fontsize = 18, font = :bold)

    colorbar_kwargs = Dict{Symbol,Any}(pairs(pop!(plot_kwargs, :colorbar_kwargs, (;))))
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot, true)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    axes = Vector{Axis}(undef, n)
    for (idx, tf) in enumerate(tf_plots)
        row = div(idx - 1, n_cols) + 1
        col = mod1(idx, n_cols)
        title = show_title ? "$(tf.condition_name) — $freq_str" : ""
        ax = Axis(fig[row, col], aspect = DataAspect(), title = title, titlesize = plot_kwargs[:plot_title_fontsize])
        axes[idx] = ax

        _render_topo_surface!(
            fig,
            ax,
            all_layout_values[idx],
            layout;
            data_interp = can_batch ? batched_interp[idx] : nothing,
            x_bounds = can_batch ? xb_shared : nothing,
            y_bounds = can_batch ? yb_shared : nothing,
            method = method,
            gridscale = gridscale,
            colormap = actual_colormap,
            colorrange = colorrange,
            num_levels = num_levels,
            plot_kwargs...,
        )
        !isnothing(highlight_channels) && _draw_highlight_channels!(ax, highlight_channels, layout)
        hidedecorations!(ax)
    end

    cb = if colorbar_plot
        cb_label = _tf_colorbar_label(tf_plots[1], baseline_interval, baseline_method)
        Colorbar(fig[1:n_rows, n_cols+1]; colorbar_kwargs..., colormap = actual_colormap, colorrange = colorrange, label = cb_label)
    else
        nothing
    end

    display_plot && _display_figure(fig)
    return (fig = fig, axes = axes, colorbar = cb)
end

# === TOPO STATS FOR TF STATISTICAL RESULTS ===

function _tf_bin_tvalue_means(
    t_mat::Array{Float64,3},
    elec_indices::Vector{Int},
    freq_indices::Vector{Int},
    bin_indices::AbstractVector{Int},
)::Vector{Float64}
    n_elec = length(elec_indices)
    out = Vector{Float64}(undef, n_elec)
    inv_n = 1.0 / (length(freq_indices) * length(bin_indices))
    @inbounds for k = 1:n_elec
        ei = elec_indices[k]
        s = 0.0
        for ti in bin_indices, fi in freq_indices
            s += t_mat[ei, fi, ti]
        end
        m = s * inv_n
        out[k] = isfinite(m) ? m : 0.0
    end
    return out
end

function _tf_bin_sig_mask(
    pos::AbstractArray{Bool,3},
    neg::AbstractArray{Bool,3},
    elec_indices::Vector{Int},
    freq_indices::Vector{Int},
    bin_indices::AbstractVector{Int},
    threshold::Real,
)::BitVector
    n_elec = length(elec_indices)
    mask = BitVector(undef, n_elec)
    inv_n = 1.0 / (length(freq_indices) * length(bin_indices))
    @inbounds for k = 1:n_elec
        ei = elec_indices[k]
        cnt = 0
        for ti in bin_indices, fi in freq_indices
            cnt += (pos[ei, fi, ti] | neg[ei, fi, ti])
        end
        mask[k] = (cnt * inv_n) >= threshold
    end
    return mask
end

"""
    plot_topography_stats(result::TFStatsResult;
                     freq_range::Tuple{Real, Real},
                     n_topos::Int = 10,
                     interval_selection::Interval = times(),
                     topo_data::Symbol = :tvalues,
                     highlight_significant::Bool = true,
                     highlight_color = :white,
                     highlight_marker::Symbol = :circle,
                     highlight_size::Real = 8,
                     highlight_threshold::Real = 0.5,
                     kwargs...)

Plot a grid of topographic maps from TF statistical results across time windows
for a selected frequency band.

Each panel shows the spatial distribution of t-statistics (or power difference)
averaged within a time window and frequency band, with significant channels highlighted.

A channel is marked significant when the proportion of significant time×frequency bins
within the panel's window meets or exceeds `highlight_threshold`.

# Arguments
- `result::TFStatsResult`: TF statistical result (from `permutation_test` or `analytic_test`)
- `freq_range::Tuple{Real, Real}`: Frequency range to average over (e.g., `(8.0, 12.0)`)
- `n_topos::Int`: Number of topographic panels (default: 10)
- `interval_selection::Interval`: Time window to display (default: full range). Use `times(start, stop)` to specify
- `topo_data::Symbol`: `:tvalues` (default) or `:difference`
- `highlight_significant::Bool`: Overlay markers on significant channels (default: true)
- `highlight_threshold::Real`: Proportion of time×freq bins that must be significant
  to mark a channel (default: 0.5). `0.0` = any bin, `1.0` = all bins
- `highlight_color`, `highlight_marker`, `highlight_size`: Marker styling

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- Named tuple `(fig, axes, colorbar)`

# Examples
```julia
# Alpha-band statistics
plot_topography_stats(result, freq_range=(8.0, 12.0))

# Theta-band with more panels and lower threshold
plot_topography_stats(result, freq_range=(4.0, 7.0), n_topos=15, highlight_threshold=0.3)

# Specific time window
plot_topography_stats(result, freq_range=(8.0, 12.0), interval_selection=times(0.1, 0.5))
```
"""
function plot_topography_stats(
    result::TFStatsResult;
    freq_range::Interval,
    n_topos::Int = 10,
    interval_selection::Interval = times(),
    topo_data::Symbol = :tvalues,
    highlight_significant::Bool = true,
    highlight_color = :white,
    highlight_marker::Symbol = :circle,
    highlight_size::Real = 8,
    highlight_threshold::Real = 0.5,
    kwargs...,
)
    topo_data in (:tvalues, :difference) || @minimal_error("topo_data must be :tvalues or :difference, got :$topo_data")

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    colorrange = pop!(plot_kwargs, :colorrange)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    plot_title = pop!(plot_kwargs, :plot_title, "")
    show_title = pop!(plot_kwargs, :show_title, true)
    num_levels = pop!(plot_kwargs, :num_levels)
    highlight_channels = pop!(plot_kwargs, :highlight_channels, nothing)

    _pop_unused_tf_topo_kwargs!(plot_kwargs)
    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)

    electrodes = result.electrodes
    all_time_points = result.time_points
    all_frequencies = result.frequencies

    fr = _sorted_interval_bounds(freq_range)
    freq_indices = [i for (i, f) in enumerate(all_frequencies) if fr[1] <= f <= fr[2]]
    isempty(freq_indices) &&
        @minimal_error("No frequencies found in range $(freq_range). Data range: $(first(all_frequencies)) to $(last(all_frequencies)) Hz")

    if isnothing(interval_selection) || interval_selection isa AllSelection
        t_start, t_end = first(all_time_points), last(all_time_points)
    else
        t_start = interval_selection isa TimeSelection ? interval_selection.start : interval_selection[1]
        t_end = interval_selection isa TimeSelection ? interval_selection.stop : interval_selection[2]
    end
    time_indices = [i for (i, t) in enumerate(all_time_points) if t_start <= t <= t_end]
    isempty(time_indices) && @minimal_error(
        "No time points found in range ($t_start, $t_end). Data range: $(first(all_time_points)) to $(last(all_time_points))"
    )

    n_topos = min(n_topos, length(time_indices))
    bins = _partition_indices(time_indices, n_topos)

    orig_layout = result.data[1].layout
    common_channels = intersect(orig_layout.data[!, :label], electrodes)
    isempty(common_channels) && @minimal_error("None of the analyzed electrodes $(electrodes) were found in the layout.")
    layout = _prepare_tf_topo_layout(orig_layout, common_channels; need_3d = true)

    layout_labels = layout.data[!, :label]
    layout_to_elec_idx = Int[findfirst(==(lbl), electrodes)::Int for lbl in layout_labels]

    topo_values = Vector{Vector{Float64}}(undef, n_topos)
    sig_masks = Vector{BitVector}(undef, n_topos)
    bin_time_labels = Vector{String}(undef, n_topos)

    bin_freqs = topo_data == :difference ? all_frequencies[freq_indices] : Float64[]
    for (i, bin_indices) in enumerate(bins)
        bin_time_start = all_time_points[first(bin_indices)]
        bin_time_end = all_time_points[last(bin_indices)]
        bin_time_labels[i] = @sprintf("%.0f – %.0f ms", bin_time_start * 1000, bin_time_end * 1000)

        if topo_data == :tvalues
            topo_values[i] = _tf_bin_tvalue_means(result.stat_matrix.t, layout_to_elec_idx, freq_indices, bin_indices)
        else
            row_indices = _tf_filter_row_indices(result.data[1].data_power, bin_freqs, @view(all_time_points[bin_indices]))
            topo_values[i] =
                _tf_extract_channel_difference_means(result.data[1].data_power, result.data[2].data_power, layout_labels, row_indices)
        end

        sig_masks[i] = _tf_bin_sig_mask(
            result.masks.positive,
            result.masks.negative,
            layout_to_elec_idx,
            freq_indices,
            bin_indices,
            highlight_threshold,
        )
    end

    if isnothing(colorrange)
        max_abs = 0.0
        for vals in topo_values
            @inbounds for v in vals
                isfinite(v) && (max_abs = max(max_abs, abs(v)))
            end
        end
        max_abs ≈ 0.0 && (max_abs = 1.0)
        colorrange = (-max_abs, max_abs)
    end

    dims = pop!(plot_kwargs, :dims, nothing)
    if isnothing(dims)
        n_rows, n_cols = _best_rect(n_topos)
    else
        n_rows, n_cols = dims
        n_rows * n_cols < n_topos &&
            @minimal_error("Grid dimensions ($n_rows × $n_cols) provide $(n_rows * n_cols) cells but need $n_topos")
    end

    test_type = result isa TFAnalyticResult ? "Analytic" : "Permutation"
    data_label = topo_data == :tvalues ? "t-statistic" : "Power Difference"
    freq_str = "$(round(Int, fr[1]))–$(round(Int, fr[2])) Hz"
    user_title = !isnothing(plot_title) && !isempty(plot_title) ? plot_title : figure_title
    fig_title = isempty(user_title) ? "$test_type Test — $data_label ($freq_str)" : user_title

    fig = _eegfun_figure(size = (200 * n_cols + 100, 200 * n_rows + 50))
    (show_title && !isempty(fig_title)) && Label(fig[0, 1:n_cols], fig_title, fontsize = 18, font = :bold)

    colorbar_kwargs = Dict{Symbol,Any}(pairs(pop!(plot_kwargs, :colorbar_kwargs, (;))))
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot, true)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    actual_colormap = _resolve_theme_colormap(fig, colormap)

    can_batch = n_topos > 1 && _topo_method_supported(method)
    batched_interp = nothing
    xb_shared = nothing
    yb_shared = nothing
    if can_batch
        all_layout_mat = Matrix{Float64}(undef, length(layout_labels), n_topos)
        for i = 1:n_topos
            channel_data = topo_values[i]
            @inbounds for k in eachindex(layout_labels)
                v = channel_data[k]
                all_layout_mat[k, i] = isfinite(v) ? v : 0.0
            end
        end
        batched_interp, xb_shared, yb_shared = _topo_interpolate(all_layout_mat, layout, gridscale; method = method)
    end

    axes = Vector{Axis}(undef, n_topos)
    x2_col = layout.data[!, :x2]
    y2_col = layout.data[!, :y2]

    for i = 1:n_topos
        row = div(i - 1, n_cols) + 1
        col = mod1(i, n_cols)

        ax = Axis(
            fig[row, col],
            aspect = DataAspect(),
            title = show_title ? bin_time_labels[i] : "",
            titlesize = plot_kwargs[:plot_title_fontsize],
        )
        axes[i] = ax

        _render_topo_surface!(
            fig,
            ax,
            topo_values[i],
            layout;
            data_interp = can_batch ? batched_interp[i] : nothing,
            x_bounds = can_batch ? xb_shared : nothing,
            y_bounds = can_batch ? yb_shared : nothing,
            method = method,
            gridscale = gridscale,
            colormap = actual_colormap,
            colorrange = colorrange,
            num_levels = num_levels,
            plot_kwargs...,
        )

        !isnothing(highlight_channels) && _draw_highlight_channels!(ax, highlight_channels, layout)

        if highlight_significant && any(sig_masks[i])
            sig_idx = findall(sig_masks[i])
            if !isempty(sig_idx)
                scatter!(
                    ax,
                    x2_col[sig_idx],
                    y2_col[sig_idx],
                    color = highlight_color,
                    marker = highlight_marker,
                    markersize = highlight_size,
                    strokewidth = 1,
                    strokecolor = :black,
                )
            end
        end

        hidedecorations!(ax)
    end

    cb = if colorbar_plot
        cb_label = topo_data == :tvalues ? "t-statistic" : "Power Difference"
        Colorbar(fig[1:n_rows, n_cols+1]; colorbar_kwargs..., colormap = actual_colormap, colorrange = colorrange, label = cb_label)
    else
        nothing
    end

    display_plot && _display_figure(fig)
    return (fig = fig, axes = axes, colorbar = cb)
end

# === TF TOPOGRAPHY EXTRACTION HELPERS ===

function _tf_filter_rows_typed(
    times::AbstractVector,
    freqs::AbstractVector,
    t_min::Float64,
    t_max::Float64,
    f_min::Float64,
    f_max::Float64,
)::Vector{Int}
    row_indices = Int[]
    @inbounds for i in eachindex(times, freqs)
        t = Float64(times[i])
        f = Float64(freqs[i])
        if (t_min <= t <= t_max) && (f_min <= f <= f_max)
            push!(row_indices, i)
        end
    end
    return row_indices
end

function _tf_filter_row_indices(
    df::DataFrame,
    selected_freqs::AbstractVector{Float64},
    selected_times::AbstractVector{Float64},
)::Vector{Int}
    f_min, f_max = extrema(selected_freqs)
    t_min, t_max = extrema(selected_times)
    return _tf_filter_rows_typed(df[!, :time], df[!, :freq], t_min, t_max, f_min, f_max)
end

function _tf_col_mean_finite(col::AbstractVector, row_indices::Vector{Int})::Float64
    sum_val = 0.0
    n_finite = 0
    @inbounds @simd for idx in eachindex(row_indices)
        v = Float64(col[row_indices[idx]])
        if isfinite(v)
            sum_val += v
            n_finite += 1
        end
    end
    return n_finite == 0 ? 0.0 : sum_val / n_finite
end

function _tf_extract_channel_means(df::DataFrame, channels::Vector{Symbol}, row_indices::Vector{Int})::Vector{Float64}
    channel_data = Vector{Float64}(undef, length(channels))
    for (ci, ch) in enumerate(channels)
        channel_data[ci] = _tf_col_mean_finite(df[!, ch], row_indices)
    end
    return channel_data
end

function _tf_col_diff_mean_finite(col_a::AbstractVector, col_b::AbstractVector, row_indices::Vector{Int})::Float64
    sum_val = 0.0
    n_finite = 0
    @inbounds @simd for idx in eachindex(row_indices)
        r = row_indices[idx]
        v = Float64(col_a[r]) - Float64(col_b[r])
        if isfinite(v)
            sum_val += v
            n_finite += 1
        end
    end
    return n_finite == 0 ? 0.0 : sum_val / n_finite
end

function _tf_extract_channel_difference_means(
    df_a::DataFrame,
    df_b::DataFrame,
    channels::Vector{Symbol},
    row_indices::Vector{Int},
)::Vector{Float64}
    diff_data = Vector{Float64}(undef, length(channels))
    for (ci, ch) in enumerate(channels)
        diff_data[ci] = _tf_col_diff_mean_finite(df_a[!, ch], df_b[!, ch], row_indices)
    end
    return diff_data
end
