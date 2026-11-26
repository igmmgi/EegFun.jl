# Statistical t-map topographies and cluster-highlighted scalp maps.

######################################
########### User Functions ###########
######################################

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

    all_time_points = result.time_points
    bins, bin_time_labels = _select_topo_time_bins(all_time_points, interval_selection, n_topos)
    layout, layout_labels, layout_to_elec_idx = _prepare_stats_topo_layout(result.data[1].layout, result.electrodes)

    sig_2d = result.masks.positive .| result.masks.negative
    topo_values, sig_masks = if topo_data == :tvalues
        _bin_2d_topo_stats(result.stat_matrix.t, sig_2d, layout_to_elec_idx, bins, highlight_threshold)
    else
        cond_A_df = result.data[1].data
        cond_B_df = result.data[2].data
        df_times = cond_A_df[!, :time]
        df_row_map = Int[argmin(abs.(df_times .- t)) for t in all_time_points]
        has_cols = Bool[hasproperty(cond_A_df, ch) && hasproperty(cond_B_df, ch) for ch in layout_labels]
        n_bins = length(bins)
        n_ch = length(layout_labels)
        vals_vec = Vector{Vector{Float64}}(undef, n_bins)
        for (i, bin_indices) in enumerate(bins)
            df_rows = view(df_row_map, bin_indices)
            vals = Vector{Float64}(undef, n_ch)
            for (k, ch_sym) in enumerate(layout_labels)
                vals[k] = has_cols[k] ? _mean_col_diff(cond_A_df[!, ch_sym], cond_B_df[!, ch_sym], df_rows) : NaN
            end
            vals_vec[i] = vals
        end
        _, masks_vec = _bin_2d_topo_stats(result.stat_matrix.t, sig_2d, layout_to_elec_idx, bins, highlight_threshold)
        vals_vec, masks_vec
    end

    test_type = result isa AnalyticResult ? "Analytic" : "Permutation"
    data_label = topo_data == :tvalues ? "t-statistic" : "Difference (μV)"
    return _render_topo_grid!(
        layout,
        topo_values,
        bin_time_labels,
        kwargs;
        default_title = "$test_type Test — $data_label",
        colorbar_label = data_label,
        sig_masks = sig_masks,
        highlight_significant = highlight_significant,
        highlight_color = highlight_color,
        highlight_marker = highlight_marker,
        highlight_size = highlight_size,
    )
end

"""
    plot_topography_stats(result::LmmStatsResult;
                          coef = 2,
                          topo_data::Symbol = :tvalues,
                          cluster_thresh::Real = 2.0,
                          threshold_p = nothing,
                          n_topos::Int = 10,
                          interval_selection::Interval = times(),
                          highlight_significant::Bool = true,
                          highlight_color = :white,
                          highlight_marker::Symbol = :circle,
                          highlight_size::Real = 8,
                          highlight_threshold::Real = 0.5,
                          kwargs...)

Plot a grid of topographic maps showing LMM statistical results (`:tvalues` or `:beta`/`:difference`)
across time windows with significant channels highlighted.
"""
function plot_topography_stats(
    result::LmmStatsResult;
    coef = 2,
    coef_idx::Union{Nothing,Int} = nothing,
    topo_data::Symbol = :tvalues,
    cluster_thresh::Real = 2.0,
    threshold_p = nothing,
    n_topos::Int = 10,
    interval_selection::Interval = times(),
    highlight_significant::Bool = true,
    highlight_color = :white,
    highlight_marker::Symbol = :circle,
    highlight_size::Real = 8,
    highlight_threshold::Real = 0.5,
    kwargs...,
)
    topo_data in (:tvalues, :beta, :difference) || @minimal_error("topo_data must be :tvalues, :beta, or :difference, got :$topo_data")
    c_idx = _resolve_lmm_coef_idx(result, something(coef_idx, coef))
    isnothing(result.epochs) &&
        @minimal_error("plot_topography_stats needs a layout: use a result from the EpochData method of fit_mass_lmm.")

    bins, bin_time_labels = _select_topo_time_bins(result.time_points, interval_selection, n_topos)
    layout, _, layout_to_elec_idx = _prepare_stats_topo_layout(result.epochs.layout, result.channels)

    z_matrix = _lmm_stat_slice(result, c_idx)
    val_matrix = topo_data === :tvalues ? z_matrix : @view(result.beta[:, :, c_idx])
    sig_combined_full = if !isnothing(threshold_p)
        p = (!isnothing(result.p_corrected) && c_idx in result.tested) ? result.p_corrected : result.p_parametric
        @view(p[:, :, c_idx]) .<= threshold_p
    elseif !isnothing(result.p_corrected) && c_idx in result.tested
        @view(result.p_corrected[:, :, c_idx]) .<= 0.05
    else
        abs.(z_matrix) .>= cluster_thresh
    end

    topo_values, sig_masks = _bin_2d_topo_stats(val_matrix, sig_combined_full, layout_to_elec_idx, bins, highlight_threshold)

    coef_name = result.coefficients[c_idx]
    data_label = topo_data === :tvalues ? "t-statistic" : "β (μV)"
    return _render_topo_grid!(
        layout,
        topo_values,
        bin_time_labels,
        kwargs;
        default_title = "LMM: $(coef_name) — $data_label",
        colorbar_label = data_label,
        sig_masks = sig_masks,
        highlight_significant = highlight_significant,
        highlight_color = highlight_color,
        highlight_marker = highlight_marker,
        highlight_size = highlight_size,
    )
end

######################################
######### Internal Functions #########
######################################

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

function _select_topo_time_bins(all_time_points::AbstractVector{<:Real}, interval_selection::Interval, n_topos::Int)
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

    n_bins = min(n_topos, length(time_indices))
    bins = _partition_indices(time_indices, n_bins)
    bin_labels = String[@sprintf("%.0f – %.0f ms", all_time_points[first(b)] * 1000, all_time_points[last(b)] * 1000) for b in bins]
    return bins, bin_labels
end

function _prepare_topo_subset_layout(orig_layout::Layout, common_channels::Vector{Symbol}; need_3d::Bool = true)::Layout
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

function _prepare_stats_topo_layout(orig_layout::Layout, electrodes::AbstractVector{Symbol}; need_3d::Bool = true)
    common_channels = intersect(orig_layout.data[!, :label], electrodes)
    isempty(common_channels) && @minimal_error("None of the analyzed electrodes $(electrodes) were found in the layout.")
    layout = _prepare_topo_subset_layout(orig_layout, common_channels; need_3d = need_3d)
    layout_labels = layout.data[!, :label]
    layout_to_elec_idx = Int[findfirst(==(lbl), electrodes)::Int for lbl in layout_labels]
    return layout, layout_labels, layout_to_elec_idx
end

function _bin_2d_topo_stats(
    val_matrix::AbstractMatrix{<:Real},
    sig_matrix::AbstractMatrix{Bool},
    elec_indices::Vector{Int},
    bins::Vector{UnitRange{Int}},
    highlight_threshold::Real,
)
    n_bins = length(bins)
    n_ch = length(elec_indices)
    topo_values = Vector{Vector{Float64}}(undef, n_bins)
    sig_masks = Vector{BitVector}(undef, n_bins)

    for (i, bin_indices) in enumerate(bins)
        n_bin_points = length(bin_indices)
        vals = Vector{Float64}(undef, n_ch)
        sm = BitVector(undef, n_ch)
        @inbounds for k = 1:n_ch
            ei = elec_indices[k]
            s = 0.0
            cnt_f = 0
            cnt_sig = 0
            for t_idx in bin_indices
                v = Float64(val_matrix[ei, t_idx])
                if isfinite(v)
                    s += v
                    cnt_f += 1
                end
                sig_matrix[ei, t_idx] && (cnt_sig += 1)
            end
            vals[k] = cnt_f > 0 ? s / cnt_f : 0.0
            sm[k] = (cnt_sig / n_bin_points) >= highlight_threshold
        end
        topo_values[i] = vals
        sig_masks[i] = sm
    end
    return topo_values, sig_masks
end

function _render_topo_grid!(
    layout::Layout,
    panel_values::Vector{Vector{Float64}},
    panel_titles::Vector{String},
    raw_kwargs;
    default_title::String = "",
    colorbar_label::String = "",
    sig_masks::Union{Nothing,Vector{BitVector}} = nothing,
    highlight_significant::Bool = false,
    highlight_color = :white,
    highlight_marker::Symbol = :circle,
    highlight_size::Real = 8,
    use_subgrid::Bool = true,
    panel_px::NTuple{4,Int} = (220, 140, 220, 70),
    min_width::Int = 600,
    header_fontsize::Int = 16,
    auto_colorrange = nothing,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, raw_kwargs)

    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    colorrange = pop!(plot_kwargs, :colorrange)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title, "")
    plot_title = pop!(plot_kwargs, :plot_title, "")
    show_title = pop!(plot_kwargs, :show_title, true)
    num_levels = pop!(plot_kwargs, :num_levels)
    dims = pop!(plot_kwargs, :dims, nothing)
    highlight_channels = pop!(plot_kwargs, :highlight_channels, nothing)

    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)

    plot_kwargs[:label_plot] = get(raw_kwargs, :label_plot, false)
    plot_kwargs[:point_plot] = get(raw_kwargs, :point_plot, false)

    n_panels = length(panel_values)
    n_rows, n_cols = if isnothing(dims)
        _best_rect(n_panels)
    else
        r, c = dims
        r * c < n_panels && @minimal_error("Grid dimensions ($r × $c) provide $(r * c) cells but need $n_panels")
        r, c
    end

    user_title =
        (!isnothing(plot_title) && !isempty(plot_title)) ? String(plot_title) : (isnothing(figure_title) ? "" : String(figure_title))
    fig_title = isempty(user_title) ? default_title : user_title

    cw, wpad, rh, hpad = panel_px
    fig = _eegfun_figure(size = (max(cw * n_cols + wpad, min_width), rh * n_rows + hpad))
    actual_colormap = _resolve_theme_colormap(fig, colormap)

    if isnothing(colorrange)
        if !isnothing(auto_colorrange)
            colorrange = auto_colorrange(actual_colormap)
        else
            max_abs = 0.0
            for vals in panel_values
                @inbounds for v in vals
                    isfinite(v) && abs(v) > max_abs && (max_abs = abs(v))
                end
            end
            max_abs ≈ 0.0 && (max_abs = 1.0)
            colorrange = (-max_abs, max_abs)
        end
    end

    target_gl = use_subgrid ? (fig[1, 1] = GridLayout()) : fig.layout
    (show_title && !isempty(fig_title)) && Label(target_gl[0, 1:n_cols], fig_title, fontsize = header_fontsize, font = :bold)

    colorbar_kwargs = Dict{Symbol,Any}(pairs(pop!(plot_kwargs, :colorbar_kwargs, (;))))
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot, true)
    pop!(plot_kwargs, :colorbar_position, nothing)
    pop!(plot_kwargs, :colorbar_plot_numbers, nothing)

    n_ch = length(layout.data[!, :label])
    can_batch = n_panels > 1 && _topo_method_supported(method) && n_ch >= 3
    batched_interp = xb_shared = yb_shared = nothing
    if can_batch
        all_layout_mat = Matrix{Float64}(undef, n_ch, n_panels)
        @inbounds for i = 1:n_panels
            vals = panel_values[i]
            for k = 1:n_ch
                all_layout_mat[k, i] = vals[k]
            end
        end
        batched_interp, xb_shared, yb_shared = _topo_interpolate(all_layout_mat, layout, gridscale; method = method)
    end

    axes = Vector{Axis}(undef, n_panels)
    x2_col = layout.data[!, :x2]
    y2_col = layout.data[!, :y2]

    for i = 1:n_panels
        row = div(i - 1, n_cols) + 1
        col = mod1(i, n_cols)
        ax = Axis(
            target_gl[row, col],
            aspect = DataAspect(),
            title = show_title ? panel_titles[i] : "",
            titlesize = plot_kwargs[:plot_title_fontsize],
        )
        axes[i] = ax

        _render_topo_surface!(
            fig,
            ax,
            panel_values[i],
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

        if highlight_significant && !isnothing(sig_masks) && any(sig_masks[i])
            sig_idx = findall(sig_masks[i])
            scatter!(
                ax,
                x2_col[sig_idx],
                y2_col[sig_idx];
                color = highlight_color,
                marker = highlight_marker,
                markersize = highlight_size,
                strokewidth = 1,
                strokecolor = :black,
            )
        end

        hidedecorations!(ax)
    end

    if use_subgrid
        for c = 1:n_cols
            colsize!(target_gl, c, Aspect(1, 1.0))
        end
    end

    cb = if colorbar_plot
        cb_slot = use_subgrid ? fig[1, 2] : fig[1:n_rows, n_cols+1]
        Colorbar(cb_slot; colorbar_kwargs..., colormap = actual_colormap, colorrange = colorrange, label = colorbar_label)
    else
        nothing
    end

    display_plot && _display_figure(fig)
    return (fig = fig, axes = axes, colorbar = cb)
end
