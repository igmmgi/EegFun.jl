"""
Topographic visualizations for Time-Frequency data.

`plot_topography` for `TimeFreqData` — single topo of average power in a time-freq window.
`plot_topography_stats` for `TFStatsResult` — grid of topos across time for a selected frequency band.
"""


# === TOPOGRAPHY FOR TIME-FREQUENCY DATA ===
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
- Named tuple `(fig, ax)`

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
    freq_range::Tuple{Real,Real},
    interval_selection::Interval = times(),
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    baseline_method::Symbol = :db,
    kwargs...,
)
    # Apply baseline correction if requested
    if !isnothing(baseline_interval) && !isnothing(tf.baseline)
        @warn "Data has already been baselined (method: $(tf.baseline.method)). " * "Ignoring baseline_interval parameter."
        tf_plot = tf
    elseif !isnothing(baseline_interval)
        tf_plot = tf_baseline(tf, baseline_interval; method = baseline_method)
    else
        tf_plot = tf
    end

    layout = tf_plot.layout

    # Validate layout
    if !has_valid_coordinates(layout)
        @minimal_error("Cannot create topographic plot: layout has no spatial coordinates.")
    end

    # Extract unique time and frequency values
    all_times = Float64.(sort(unique(tf_plot.data_power.time)))
    all_freqs = Float64.(sort(unique(tf_plot.data_power.freq)))

    # Apply frequency range filter
    freq_mask = freq_range[1] .<= all_freqs .<= freq_range[2]
    isempty(findall(freq_mask)) &&
        @minimal_error("No frequencies found in range $(freq_range). Data range: $(first(all_freqs)) to $(last(all_freqs)) Hz")

    # Apply time range filter
    if isnothing(interval_selection) || interval_selection isa AllSelection
        time_mask = fill(true, length(all_times))
    else
        t_start = interval_selection isa TimeSelection ? interval_selection.start : interval_selection[1]
        t_stop = interval_selection isa TimeSelection ? interval_selection.stop : interval_selection[2]
        time_mask = t_start .<= all_times .<= t_stop
    end
    isempty(findall(time_mask)) &&
        @minimal_error("No time points found in interval. Data range: $(first(all_times)) to $(last(all_times)) s")

    selected_times = all_times[time_mask]

    # Get channel labels
    channels_list = channel_labels(tf_plot)

    # Compute average power per channel across the selected time-freq window
    channel_data = Float64[]
    for ch in channels_list
        power_mat = _tf_df_to_matrix(tf_plot.data_power, ch, all_freqs, all_times)  # [freqs × time]
        # Subset to selected freq and time ranges; filter non-finite values
        # (dB baseline can produce -Inf for near-zero power)
        vals = filter(isfinite, vec(power_mat[freq_mask, time_mask]))
        avg_power = isempty(vals) ? 0.0 : mean(vals)
        push!(channel_data, avg_power)
    end

    # Map channel data onto layout order
    layout_labels = layout.data.label
    layout_values = Float64[let ch_idx = findfirst(==(lbl), channels_list)
        !isnothing(ch_idx) ? channel_data[ch_idx] : 0.0
    end for lbl in layout_labels]

    # Merge user kwargs with shared topography defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # Extract topography-specific parameters
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    ylim = pop!(plot_kwargs, :ylim)
    # For TF power, use actual data range (not symmetric around 0)
    # since baselined power can be uniformly negative/positive
    if isnothing(ylim)
        data_min, data_max = extrema(layout_values)
        if data_min ≈ data_max
            ylim = (data_min - 1.0, data_max + 1.0)
        else
            ylim = (data_min, data_max)
        end
    end
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    num_levels = pop!(plot_kwargs, :num_levels)

    # Pop unused keys
    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)
    pop!(plot_kwargs, :dims, nothing)

    # Build title
    freq_str = "$(round(Int, freq_range[1]))-$(round(Int, freq_range[2])) Hz"
    time_str = isnothing(interval_selection) ? "all times" : @sprintf("%.0f–%.0f ms", selected_times[1] * 1000, selected_times[end] * 1000)
    default_title = "$(tf_plot.condition_name) — $freq_str, $time_str"
    title = isempty(figure_title) ? default_title : figure_title

    # Create figure
    fig = Figure(size = (400, 400))
    ax = Axis(fig[1, 1], aspect = DataAspect(), title = title, titlesize = plot_kwargs[:plot_title_fontsize])

    # Ensure coordinates
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    # Extract colorbar kwargs
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    # Render
    co, _, _ = _render_topo_surface!(
        fig,
        ax,
        layout_values,
        layout;
        method = method,
        gridscale = gridscale,
        colormap = colormap,
        ylim = ylim,
        num_levels = num_levels,
        plot_kwargs...,
    )

    if colorbar_plot
        cb_label = _tf_colorbar_label(tf_plot, baseline_interval, baseline_method)
        Colorbar(fig[1, 2], co; colorbar_kwargs..., label = cb_label)
    end

    hidedecorations!(ax)

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, ax = ax)
end
function plot_topography(
    tfs::Vector{TimeFreqData};
    freq_range::Tuple{Real,Real},
    interval_selection::Interval = times(),
    baseline_interval::Union{Nothing,Tuple{Real,Real}} = nothing,
    baseline_method::Symbol = :db,
    kwargs...,
)
    isempty(tfs) && @minimal_error("Cannot plot empty vector of TimeFreqData")

    # Apply baseline correction
    if !isnothing(baseline_interval)
        tf_plots = map(tfs) do tf
            if !isnothing(tf.baseline)
                @warn "$(tf.condition_name) already baselined, skipping."
                tf
            else
                tf_baseline(tf, baseline_interval; method = baseline_method)
            end
        end
    else
        tf_plots = tfs
    end

    n = length(tf_plots)
    if n == 1
        return plot_topography(tf_plots[1]; freq_range = freq_range, interval_selection = interval_selection, kwargs...)
    end

    layout = tf_plots[1].layout
    if !has_valid_coordinates(layout)
        @minimal_error("Cannot create topographic plot: layout has no spatial coordinates.")
    end

    # Merge kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    ylim = pop!(plot_kwargs, :ylim)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    num_levels = pop!(plot_kwargs, :num_levels)
    dims = pop!(plot_kwargs, :dims, nothing)

    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)

    # Labels off in grid by default
    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)

    # Extract common TF data
    all_freqs = Float64.(sort(unique(tf_plots[1].data_power.freq)))
    all_times = Float64.(sort(unique(tf_plots[1].data_power.time)))
    freq_mask = freq_range[1] .<= all_freqs .<= freq_range[2]

    if isnothing(interval_selection) || interval_selection isa AllSelection
        time_mask = fill(true, length(all_times))
    else
        t_start = interval_selection isa TimeSelection ? interval_selection.start : interval_selection[1]
        t_stop = interval_selection isa TimeSelection ? interval_selection.stop : interval_selection[2]
        time_mask = t_start .<= all_times .<= t_stop
    end

    # Compute per-condition channel averages
    channels_list = channel_labels(tf_plots[1])
    layout_labels = layout.data.label
    all_layout_values = Vector{Vector{Float64}}(undef, n)

    global_min, global_max = Inf, -Inf
    for (idx, tf) in enumerate(tf_plots)
        channel_data = Float64[]
        for ch in channels_list
            power_mat = _tf_df_to_matrix(tf.data_power, ch, all_freqs, all_times)
            vals = filter(isfinite, vec(power_mat[freq_mask, time_mask]))
            push!(channel_data, isempty(vals) ? 0.0 : mean(vals))
        end
        layout_values = Float64[let ch_idx = findfirst(==(lbl), channels_list)
            !isnothing(ch_idx) ? channel_data[ch_idx] : 0.0
        end for lbl in layout_labels]
        all_layout_values[idx] = layout_values
        global_min = min(global_min, minimum(layout_values))
        global_max = max(global_max, maximum(layout_values))
    end

    # Shared color range (non-symmetric for TF power)
    if isnothing(ylim)
        if global_min ≈ global_max
            ylim = (global_min - 1.0, global_max + 1.0)
        else
            ylim = (global_min, global_max)
        end
    end

    # Grid layout
    if isnothing(dims)
        n_rows, n_cols = _best_rect(n)
    else
        n_rows, n_cols = dims
    end

    # Ensure coordinates
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    freq_str = "$(round(Int, freq_range[1]))-$(round(Int, freq_range[2])) Hz"
    fig = Figure(size = (250 * n_cols + 100, 250 * n_rows + 50))

    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    pop!(plot_kwargs, :colorbar_plot, nothing)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    axes = Axis[]
    for (idx, tf) in enumerate(tf_plots)
        row = div(idx - 1, n_cols) + 1
        col = mod1(idx, n_cols)
        title = "$(tf.condition_name) — $freq_str"
        ax = Axis(fig[row, col], aspect = DataAspect(), title = title, titlesize = plot_kwargs[:plot_title_fontsize])
        push!(axes, ax)

        _render_topo_surface!(
            fig,
            ax,
            all_layout_values[idx],
            layout;
            method = method,
            gridscale = gridscale,
            colormap = colormap,
            ylim = ylim,
            num_levels = num_levels,
            plot_kwargs...,
        )
        hidedecorations!(ax)
    end

    cb_label = _tf_colorbar_label(tf_plots[1], baseline_interval, baseline_method)
    cb = Colorbar(fig[1:n_rows, n_cols+1]; colorbar_kwargs..., colormap = colormap, colorrange = ylim, label = cb_label)

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = axes, colorbar = cb)
end


# === TOPO STATS FOR TF STATISTICAL RESULTS ===
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
    freq_range::Tuple{Real,Real},
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
    # Validate
    topo_data in (:tvalues, :difference) || error("topo_data must be :tvalues or :difference, got :$topo_data")

    # Merge user kwargs with shared topography defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # Extract topography-specific parameters
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    ylim = pop!(plot_kwargs, :ylim)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    num_levels = pop!(plot_kwargs, :num_levels)

    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)

    # Labels and points off by default for grid layout
    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)

    # Get data dimensions
    electrodes = result.electrodes
    all_time_points = result.time_points
    all_frequencies = result.frequencies

    # Find frequency indices in the selected range
    freq_mask = freq_range[1] .<= all_frequencies .<= freq_range[2]
    freq_indices = findall(freq_mask)
    isempty(freq_indices) &&
        error("No frequencies found in range $(freq_range). Data range: $(first(all_frequencies)) to $(last(all_frequencies)) Hz")

    # Find time indices in the selected range
    if isnothing(interval_selection) || interval_selection isa AllSelection
        t_start, t_end = first(all_time_points), last(all_time_points)
    else
        t_start = interval_selection isa TimeSelection ? interval_selection.start : interval_selection[1]
        t_end = interval_selection isa TimeSelection ? interval_selection.stop : interval_selection[2]
    end
    time_mask = t_start .<= all_time_points .<= t_end
    time_indices = findall(time_mask)
    isempty(time_indices) &&
        error("No time points found in range ($t_start, $t_end). Data range: $(first(all_time_points)) to $(last(all_time_points))")

    # Divide time into bins
    n_topos = min(n_topos, length(time_indices))
    bins = _partition_indices(time_indices, n_topos)

    # Get layout and validate
    layout = result.data[1].layout
    if !has_valid_coordinates(layout)
        error("Cannot create topographic plot: layout has no spatial coordinates.")
    end
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    layout_labels = layout.data.label

    # Compute per-bin topo values and significance masks
    topo_values = Vector{Vector{Float64}}(undef, n_topos)
    sig_masks = Vector{BitVector}(undef, n_topos)
    bin_time_labels = Vector{String}(undef, n_topos)

    for (i, bin_indices) in enumerate(bins)
        bin_time_start = all_time_points[first(bin_indices)]
        bin_time_end = all_time_points[last(bin_indices)]
        bin_time_labels[i] = @sprintf("%.0f – %.0f ms", bin_time_start * 1000, bin_time_end * 1000)

        if topo_data == :tvalues
            # Average t-values across selected freqs and time bin: stat_matrix.t is [electrodes × freqs × time]
            raw = vec(mean(result.stat_matrix.t[:, freq_indices, bin_indices], dims = (2, 3)))
            topo_values[i] = replace(raw, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
        else  # :difference
            # Compute power difference from grand average TF data per electrode
            diff_per_electrode = Vector{Float64}(undef, length(electrodes))
            for (j, ch_sym) in enumerate(electrodes)
                power_a = _tf_power_matrix(result.data[1], ch_sym, all_frequencies, all_time_points)
                power_b = _tf_power_matrix(result.data[2], ch_sym, all_frequencies, all_time_points)
                diff_mat = power_a .- power_b  # [freqs × time]
                vals = filter(isfinite, vec(diff_mat[freq_indices, bin_indices]))
                diff_per_electrode[j] = isempty(vals) ? 0.0 : mean(vals)
            end
            topo_values[i] = diff_per_electrode
        end

        # Significance: proportion of significant freq×time bins per channel in this window
        # masks are [electrodes × freqs × time]
        sig_combined = result.masks.positive[:, freq_indices, bin_indices] .| result.masks.negative[:, freq_indices, bin_indices]
        n_bins = length(freq_indices) * length(bin_indices)
        sig_proportion = vec(sum(sig_combined, dims = (2, 3))) ./ n_bins
        sig_masks[i] = sig_proportion .>= highlight_threshold
    end

    # Compute symmetric color limits
    all_vals = vcat(topo_values...)
    max_abs = maximum(abs, all_vals)
    max_abs ≈ 0.0 && (max_abs = 1.0)
    if isnothing(ylim)
        ylim = (-max_abs, max_abs)
    end

    # Grid layout
    dims = pop!(plot_kwargs, :dims, nothing)
    if isnothing(dims)
        n_rows, n_cols = _best_rect(n_topos)
    else
        n_rows, n_cols = dims
        n_rows * n_cols < n_topos && error("Grid dimensions ($n_rows × $n_cols) provide $(n_rows * n_cols) cells but need $n_topos")
    end

    # Create figure
    test_type = isa(result, TFAnalyticResult) ? "Analytic" : "Permutation"
    data_label = topo_data == :tvalues ? "t-statistic" : "Power Difference"
    freq_str = "$(round(Int, freq_range[1]))–$(round(Int, freq_range[2])) Hz"
    fig_title = isempty(figure_title) ? "$test_type Test — $data_label ($freq_str)" : figure_title

    fig = Figure(size = (200 * n_cols + 100, 200 * n_rows + 50))
    Label(fig[0, 1:n_cols], fig_title, fontsize = 18, font = :bold)

    # Extract colorbar kwargs
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    pop!(plot_kwargs, :colorbar_plot, nothing)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    axes = Axis[]

    for i = 1:n_topos
        row = div(i - 1, n_cols) + 1
        col = mod1(i, n_cols)

        ax = Axis(fig[row, col], aspect = DataAspect(), title = bin_time_labels[i], titlesize = plot_kwargs[:plot_title_fontsize])
        push!(axes, ax)

        # Map electrode values to layout order
        channel_data = topo_values[i]
        layout_values = Float64[let ch_idx = findfirst(==(lbl), electrodes)
            v = !isnothing(ch_idx) ? channel_data[ch_idx] : 0.0
            isfinite(v) ? v : 0.0
        end for lbl in layout_labels]

        _render_topo_surface!(
            fig,
            ax,
            layout_values,
            layout;
            method = method,
            gridscale = gridscale,
            colormap = colormap,
            ylim = ylim,
            num_levels = num_levels,
            plot_kwargs...,
        )

        # Highlight significant channels
        if highlight_significant && any(sig_masks[i])
            sig_channel_indices = findall(sig_masks[i])
            sig_x = Float64[]
            sig_y = Float64[]
            for ch_idx in sig_channel_indices
                ch_sym = electrodes[ch_idx]
                layout_idx = findfirst(==(ch_sym), layout_labels)
                if !isnothing(layout_idx)
                    push!(sig_x, layout.data.x2[layout_idx])
                    push!(sig_y, layout.data.y2[layout_idx])
                end
            end
            if !isempty(sig_x)
                scatter!(
                    ax,
                    sig_x,
                    sig_y,
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

    # Shared colorbar
    cb_label = topo_data == :tvalues ? "t-statistic" : "Power Difference"
    cb = Colorbar(fig[1:n_rows, n_cols+1]; colorbar_kwargs..., colormap = colormap, colorrange = ylim, label = cb_label)

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = axes, colorbar = cb)
end
