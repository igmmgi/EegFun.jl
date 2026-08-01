"""
Topographic visualization of statistical test results across time windows.
Creates a grid of topographic maps showing t-statistics or difference amplitudes
with significant channels highlighted — inspired by FieldTrip's ft_clusterplot.
"""


"""
    plot_topography_stats(result::StatsResult;
                     n_topos::Int = 10,
                     interval_selection::Interval = times(),
                     topo_data::Symbol = :tvalues,
                     highlight_significant::Bool = true,
                     highlight_color = :white,
                     highlight_marker::Symbol = :circle,
                     highlight_size::Real = 8,
                     highlight_threshold::Real = 0.5,
                     kwargs...)

Plot a grid of topographic maps showing statistical results across time windows.

Each panel shows the spatial distribution of t-statistics (or difference amplitudes)
averaged within a time window, with significant channels highlighted as markers.

# Arguments
- `result::StatsResult`: Results from `analytic_test` or `permutation_test`
- `n_topos::Int`: Number of topographic panels (default: 10)
- `interval_selection::Interval`: Time window to display (default: full range). Use `times(start, stop)` to specify
- `topo_data::Symbol`: What to display on the topographic maps:
  - `:tvalues` (default) — t-statistics from the statistical test
  - `:difference` — difference wave amplitude (condition A − condition B)
- `highlight_significant::Bool`: Whether to overlay markers on significant channels (default: true)
- `highlight_color`: Color for significance markers (default: `:white`)
- `highlight_marker::Symbol`: Marker symbol for significant channels (default: `:circle`).
  Options: `:circle`, `:cross`, `:diamond`, `:star5`, `:xcross`, `:utriangle`, `:dtriangle`, etc.
- `highlight_size::Real`: Size of significance markers (default: 8)
- `highlight_threshold::Real`: Proportion of time points within a window that must be significant
  to mark a channel (default: 0.5). `0.0` = any single time point (union), `0.5` = majority,
  `1.0` = all time points (intersection)

Additional keyword arguments from PLOT_TOPOGRAPHY_KWARGS are supported:
$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- Named tuple `(fig, axes, colorbar)`

# Examples
```julia
# Basic usage with analytic test
result = analytic_test(prepared, correction_method=:no)
plot_topography_stats(result)

# Show difference wave amplitudes instead of t-values
plot_topography_stats(result, topo_data=:difference)

# Focus on a specific time window with more panels
plot_topography_stats(result, interval_selection=times(0.1, 0.4), n_topos=15)

# Custom marker style for significant channels
plot_topography_stats(result, highlight_marker=:star5, highlight_color=:yellow, highlight_size=12)
```
"""
function plot_topography_stats(
    result::StatsResult;
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
    # Validate topo_data
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

    # Pop ICA-specific keys not used here
    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)

    # Override defaults: labels and points off by default for grid layout
    # (respects explicit user override via kwargs)
    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)

    # Get time points and determine range
    all_time_points = result.time_points
    if isnothing(interval_selection) || interval_selection isa AllSelection
        t_start, t_end = first(all_time_points), last(all_time_points)
    else
        t_start = interval_selection isa TimeSelection ? interval_selection.start : interval_selection[1]
        t_end = interval_selection isa TimeSelection ? interval_selection.stop : interval_selection[2]
    end

    # Find indices within the time range
    time_mask = t_start .<= all_time_points .<= t_end
    time_indices = findall(time_mask)
    if isempty(time_indices)
        error("No time points found in range ($t_start, $t_end). Data range: $(first(all_time_points)) to $(last(all_time_points))")
    end

    # Divide time indices into n_topos equal bins
    n_topos = min(n_topos, length(time_indices))
    bins = _partition_indices(time_indices, n_topos)

    # Get layout and electrodes
    layout = result.data[1].layout
    electrodes = result.electrodes

    # Validate layout
    if !has_valid_coordinates(layout)
        error("Cannot create topographic plot: layout has no spatial coordinates.")
    end

    # Ensure coordinates
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    # Cache layout labels and valid channel set (avoids per-bin propertynames lookups)
    layout_labels = layout.data.label
    if topo_data == :difference
        cond_A = result.data[1]
        cond_B = result.data[2]
        cond_A_props = Set(propertynames(cond_A.data))
        cond_B_props = Set(propertynames(cond_B.data))
    end

    # Compute data for each bin
    topo_values = Vector{Vector{Float64}}(undef, n_topos)
    sig_masks = Vector{BitVector}(undef, n_topos)
    bin_time_labels = Vector{String}(undef, n_topos)

    for (i, bin_indices) in enumerate(bins)
        bin_time_start = all_time_points[first(bin_indices)]
        bin_time_end = all_time_points[last(bin_indices)]
        bin_time_labels[i] = @sprintf("%.0f – %.0f ms", bin_time_start * 1000, bin_time_end * 1000)

        if topo_data == :tvalues
            topo_values[i] = vec(mean(result.stat_matrix.t[:, bin_indices], dims = 2))
        else  # :difference — pre-allocated vector, uses cached channel sets
            diff_per_electrode = Vector{Float64}(undef, length(electrodes))
            for (j, ch_sym) in enumerate(electrodes)
                if ch_sym in cond_A_props && ch_sym in cond_B_props
                    a_vals = cond_A.data[bin_indices, ch_sym]
                    b_vals = cond_B.data[bin_indices, ch_sym]
                    diff_per_electrode[j] = mean(a_vals .- b_vals)
                else
                    diff_per_electrode[j] = 0.0
                end
            end
            topo_values[i] = diff_per_electrode
        end

        # Significance: channel is significant if proportion of significant time points >= threshold
        n_bin_points = length(bin_indices)
        sig_combined = result.masks.positive[:, bin_indices] .| result.masks.negative[:, bin_indices]
        sig_proportion = vec(sum(sig_combined, dims = 2)) ./ n_bin_points
        sig_masks[i] = sig_proportion .>= highlight_threshold
    end

    # Compute symmetric color limits across all panels for consistent scaling
    all_vals = vcat(topo_values...)
    max_abs = maximum(abs, all_vals)
    if max_abs ≈ 0.0
        max_abs = 1.0
    end
    # Use user-provided ylim or compute symmetric limits
    if isnothing(ylim)
        ylim = (-max_abs, max_abs)
    end

    # Grid layout (dims comes from PLOT_TOPOGRAPHY_KWARGS)
    dims = pop!(plot_kwargs, :dims, nothing)
    if isnothing(dims)
        n_rows, n_cols = _best_rect(n_topos)
    else
        n_rows, n_cols = dims
        if n_rows * n_cols < n_topos
            error("Grid dimensions ($n_rows × $n_cols) provide $(n_rows * n_cols) cells but need $n_topos")
        end
    end

    # Create figure
    test_type = isa(result, AnalyticResult) ? "Analytic" : "Permutation"
    data_label = topo_data == :tvalues ? "t-statistic" : "Difference (μV)"
    fig_title = figure_title == "Topography Plot" ? "$test_type Test — $data_label" : figure_title

    fig = Figure(size = (200 * n_cols + 100, 200 * n_rows + 50))
    if plot_kwargs[:show_title]
        Label(fig[0, 1:n_cols], fig_title, fontsize = 18, font = :bold)
    end

    # Extract colorbar kwargs before render loop (removes colorbar_* keys from plot_kwargs)
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

        ax = Axis(fig[row, col], aspect = DataAspect(), title = bin_time_labels[i], titlesize = plot_kwargs[:title_fontsize])
        push!(axes, ax)

        # Map electrode values to layout order
        channel_data = topo_values[i]
        layout_values = Float64[
            let ch_idx = findfirst(==(lbl), electrodes)
                !isnothing(ch_idx) ? channel_data[ch_idx] : 0.0
            end for lbl in layout_labels
        ]

        # Render using shared helper (uses num_levels, not gridscale, for contour levels)
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

    # Add shared colorbar and include in return value
    cb_label = topo_data == :tvalues ? "t-statistic" : "Difference (μV)"
    cb = Colorbar(fig[1:n_rows, n_cols+1]; colorbar_kwargs..., colormap = colormap, colorrange = ylim, label = cb_label)

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = axes, colorbar = cb)
end



"""
    _partition_indices(indices::Vector{Int}, n::Int) -> Vector{UnitRange{Int}}

Divide a vector of sorted, contiguous indices into `n` approximately equal bins.
Returns vector of UnitRanges spanning from `indices[bin_start]` to `indices[bin_end]`.

Note: assumes the input `indices` are sorted and contiguous (sequential integers).
For non-contiguous indices, the returned ranges may span unintended values.
"""
function _partition_indices(indices::Vector{Int}, n::Int)
    len = length(indices)
    bins = Vector{UnitRange{Int}}(undef, n)
    base_size = div(len, n)
    remainder = mod(len, n)

    start = 1
    for i = 1:n
        bin_size = base_size + (i <= remainder ? 1 : 0)
        stop = start + bin_size - 1
        bins[i] = indices[start]:indices[stop]
        start = stop + 1
    end

    return bins
end
