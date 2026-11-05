function _mean_col_diff(col_a::AbstractVector, col_b::AbstractVector, r::AbstractVector{Int})
    s = 0.0
    @inbounds @simd for k in r
        s += Float64(col_a[k]) - Float64(col_b[k])
    end
    return s / length(r)
end

function _partition_indices(indices::Vector{Int}, n::Int)
    len = length(indices)
    bins = Vector{UnitRange{Int}}(undef, n)
    base_size = div(len, n)
    remainder = mod(len, n)
    start = 1
    for i = 1:n
        stop = start + base_size + (i <= remainder ? 1 : 0) - 1
        bins[i] = indices[start]:indices[stop]
        start = stop + 1
    end
    return bins
end

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

Each panel displays t-statistics (`:tvalues`) or difference amplitudes (`:difference`)
averaged within a time bin, with significant channels highlighted as markers.

# Arguments
- `result::StatsResult`: Result from `analytic_test` or `permutation_test`.
- `n_topos::Int`: Number of topographic panels (default: `10`).
- `interval_selection::Interval`: Time window to display (default: `times()`).
- `topo_data::Symbol`: `:tvalues` (default) or `:difference`.
- `highlight_significant::Bool`: Overlay markers on significant channels (default: `true`).
- `highlight_color`: Marker color for significant channels (default: `:white`).
- `highlight_marker::Symbol`: Marker shape for significant channels (default: `:circle`).
- `highlight_size::Real`: Marker size for significant channels (default: `8`).
- `highlight_threshold::Real`: Minimum significant time-point fraction in a bin to highlight a channel (default: `0.5`).

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- Named tuple `(fig, axes, colorbar)`
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
    topo_data in (:tvalues, :difference) || @minimal_error("topo_data must be :tvalues or :difference, got :$topo_data")

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    colorrange = pop!(plot_kwargs, :colorrange)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    num_levels = pop!(plot_kwargs, :num_levels)
    show_title = pop!(plot_kwargs, :show_title, true)

    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)

    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)

    all_time_points = result.time_points
    t_start, t_end = if isnothing(interval_selection) || interval_selection isa AllSelection
        first(all_time_points), last(all_time_points)
    elseif interval_selection isa TimeSelection
        interval_selection.start, interval_selection.stop
    else
        interval_selection[1], interval_selection[2]
    end

    time_indices = findall(t -> t_start <= t <= t_end, all_time_points)
    isempty(time_indices) && @minimal_error(
        "No time points found in range ($t_start, $t_end). Data range: $(first(all_time_points)) to $(last(all_time_points))"
    )

    n_topos = min(n_topos, length(time_indices))
    bins = _partition_indices(time_indices, n_topos)

    orig_layout = result.data[1].layout
    electrodes = result.electrodes
    n_elec = length(electrodes)
    common_channels = intersect(orig_layout.data.label, electrodes)
    isempty(common_channels) && @minimal_error("None of the analyzed electrodes $(electrodes) were found in the layout.")
    layout =
        orig_layout.data.label == common_channels ? copy(orig_layout) :
        subset_layout(orig_layout; channel_selection = channels(common_channels))

    has_valid_coordinates(layout) || @minimal_error("Cannot create topographic plot: layout has no spatial coordinates.")
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    layout_labels = layout.data.label
    layout_to_elec_idx = Int[findfirst(==(lbl), electrodes)::Int for lbl in layout_labels]
    elec_to_layout_idx = Int[something(findfirst(==(ch), layout_labels), 0) for ch in electrodes]

    t_mat = result.stat_matrix.t
    pos_mask = result.masks.positive
    neg_mask = result.masks.negative

    cond_A_df = topo_data == :difference ? result.data[1].data : nothing
    cond_B_df = topo_data == :difference ? result.data[2].data : nothing
    df_row_map = if topo_data == :difference
        df_times = cond_A_df[!, :time]
        Int[argmin(abs.(df_times .- t)) for t in all_time_points]
    else
        Int[]
    end
    has_diff_cols = if topo_data == :difference
        Bool[hasproperty(cond_A_df, ch) && hasproperty(cond_B_df, ch) for ch in electrodes]
    else
        Bool[]
    end

    topo_values = Vector{Vector{Float64}}(undef, n_topos)
    sig_masks = Vector{BitVector}(undef, n_topos)
    bin_time_labels = Vector{String}(undef, n_topos)
    max_abs = 0.0

    for (i, bin_indices) in enumerate(bins)
        bin_time_labels[i] =
            @sprintf("%.0f – %.0f ms", all_time_points[first(bin_indices)] * 1000, all_time_points[last(bin_indices)] * 1000)
        n_bin_points = length(bin_indices)
        vals = Vector{Float64}(undef, n_elec)
        if topo_data == :tvalues
            inv_n = 1.0 / n_bin_points
            @inbounds for j = 1:n_elec
                s = 0.0
                for k in bin_indices
                    s += Float64(t_mat[j, k])
                end
                v = s * inv_n
                vals[j] = v
                isfinite(v) && abs(v) > max_abs && (max_abs = abs(v))
            end
        else
            df_rows = view(df_row_map, bin_indices)
            for (j, ch_sym) in enumerate(electrodes)
                v = has_diff_cols[j] ? _mean_col_diff(cond_A_df[!, ch_sym], cond_B_df[!, ch_sym], df_rows) : NaN
                vals[j] = v
                isfinite(v) && abs(v) > max_abs && (max_abs = abs(v))
            end
        end
        topo_values[i] = vals

        sm = BitVector(undef, n_elec)
        @inbounds for j = 1:n_elec
            cnt = 0
            for k in bin_indices
                (pos_mask[j, k] || neg_mask[j, k]) && (cnt += 1)
            end
            sm[j] = (cnt / n_bin_points) >= highlight_threshold
        end
        sig_masks[i] = sm
    end

    max_abs ≈ 0.0 && (max_abs = 1.0)
    isnothing(colorrange) && (colorrange = (-max_abs, max_abs))

    dims = pop!(plot_kwargs, :dims, nothing)
    n_rows, n_cols = if isnothing(dims)
        _best_rect(n_topos)
    else
        r, c = dims
        r * c < n_topos && @minimal_error("Grid dimensions ($r × $c) provide $(r * c) cells but need $n_topos")
        r, c
    end

    test_type = result isa AnalyticResult ? "Analytic" : "Permutation"
    data_label = topo_data == :tvalues ? "t-statistic" : "Difference (μV)"
    fig_title = isempty(figure_title) ? "$test_type Test — $data_label" : figure_title

    fig = _eegfun_figure(size = (max(220 * n_cols + 140, 600), 220 * n_rows + 70))
    gl_topos = fig[1, 1] = GridLayout()
    show_title && !isempty(fig_title) && Label(gl_topos[0, 1:n_cols], fig_title, fontsize = 16, font = :bold)
    for c = 1:n_cols
        colsize!(gl_topos, c, Aspect(1, 1.0))
    end

    colorbar_kwargs = Dict{Symbol,Any}(pairs(pop!(plot_kwargs, :colorbar_kwargs, (;))))
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot, true)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    axes = Axis[]
    can_batch = n_topos > 1 && _topo_method_supported(method) && length(layout_labels) >= 3
    batched_interp = nothing
    xb_shared = nothing
    yb_shared = nothing

    if can_batch
        all_layout_mat = Matrix{Float64}(undef, length(layout_labels), n_topos)
        @inbounds for i = 1:n_topos
            channel_data = topo_values[i]
            for (k, ch_idx) in enumerate(layout_to_elec_idx)
                all_layout_mat[k, i] = channel_data[ch_idx]
            end
        end
        batched_interp, xb_shared, yb_shared = _topo_interpolate(all_layout_mat, layout, gridscale; method = method)
    end

    x2_col = layout.data.x2
    y2_col = layout.data.y2

    for i = 1:n_topos
        row = div(i - 1, n_cols) + 1
        col = mod1(i, n_cols)
        ax = Axis(
            gl_topos[row, col],
            aspect = DataAspect(),
            title = show_title ? bin_time_labels[i] : "",
            titlesize = plot_kwargs[:plot_title_fontsize],
        )
        push!(axes, ax)

        channel_data = topo_values[i]
        layout_values = Float64[channel_data[ch_idx] for ch_idx in layout_to_elec_idx]

        _render_topo_surface!(
            fig,
            ax,
            layout_values,
            layout;
            data_interp = can_batch ? batched_interp[i] : nothing,
            x_bounds = can_batch ? xb_shared : nothing,
            y_bounds = can_batch ? yb_shared : nothing,
            method = method,
            gridscale = gridscale,
            colormap = colormap,
            colorrange = colorrange,
            num_levels = num_levels,
            plot_kwargs...,
        )

        if highlight_significant && any(sig_masks[i])
            sig_x = Float64[]
            sig_y = Float64[]
            for ch_idx in findall(sig_masks[i])
                l_idx = elec_to_layout_idx[ch_idx]
                if l_idx > 0
                    push!(sig_x, x2_col[l_idx])
                    push!(sig_y, y2_col[l_idx])
                end
            end
            !isempty(sig_x) && scatter!(
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

        hidedecorations!(ax)
    end

    cb = if colorbar_plot
        Colorbar(fig[1, 2]; colorbar_kwargs..., colormap = _resolve_theme_colormap(fig, colormap), colorrange = colorrange, label = data_label)
    else
        nothing
    end

    display_plot && _display_figure(fig)
    return (fig = fig, axes = axes, colorbar = cb)
end
