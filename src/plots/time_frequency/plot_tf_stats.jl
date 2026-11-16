# Plotting functions for Time-Frequency statistical test results.

"""
    plot_tf_stats(result::TFStatsResult;
                  layout::Union{Symbol, PlotLayout} = :grid,
                  layout_kwargs = (;),
                  channel_selection = channels(),
                  channel_plot_order::Union{Nothing, AbstractVector} = nothing,
                  content::Symbol = :tvalues,
                  mask_style::Symbol = :alpha,
                  colormap = nothing,
                  colorrange = nothing,
                  xscale::Symbol = :linear,
                  yscale::Symbol = :linear,
                  colorbar_plot::Bool = true,
                  plot_title::Union{Nothing, AbstractString} = nothing,
                  show_title::Bool = true,
                  significance_color = :black,
                  significance_linewidth::Real = 2.0,
                  figure_size = nothing,
                  display_plot::Bool = true)

Plot TF statistical results as heatmaps with significance overlays across `:single`, `:grid`, `:topo`, or custom `PlotLayout` layouts.

Works with both `TFClusterPermutationResult` (from `permutation_test`) and
`TFAnalyticResult` (from `analytic_test`).

# Arguments
- `result::TFStatsResult`: Statistical result to plot
- `layout::Union{Symbol, PlotLayout}`: `:single`, `:grid` (default), `:topo`, or a `PlotLayout`
- `layout_kwargs`: Additional layout configuration parameters (e.g. `grid_dims`, `topo_plot_width`, `topo_plot_height`)
- `channel_selection`: Channel selection predicate (default: all channels)
- `content::Symbol`: What to plot as the heatmap:
  - `:tvalues` (default) - t-statistic values
  - `:difference` - power difference (condition A - B)
  - `:power_a` - grand average power for condition A
  - `:power_b` - grand average power for condition B
- `mask_style::Symbol`: How to visualize significant regions:
  - `:alpha` (default) - dim non-significant regions
  - `:contour` - black contour lines around significant regions
  - `:hide` - hide non-significant regions
  - `:none` - no significance overlay
- `colormap`: Colormap (default: `:RdBu`)
- `colorrange`: Color range tuple or `nothing` for auto (auto-symmetric for t-values)
- `xscale::Symbol`: Time (x) axis scale (`:linear` or `:log10`, default: `:linear`)
- `yscale::Symbol`: Frequency (y) axis scale (`:linear` or `:log10`, default: `:linear`)
- `colorbar_plot::Bool`: Show colorbar (default: `true`)
- `plot_title`: Optional header title override (default: `nothing`)
- `show_title::Bool`: Show header and channel titles (default: `true`)
- `significance_color`: Color for contour lines (default: `:black`)
- `significance_linewidth::Real`: Width of contour lines (default: 2.0)
- `figure_size`: Figure size tuple or `nothing` for auto
- `display_plot::Bool`: Display the plot (default: true)

# Returns
Named tuple `(fig, axes)` with Makie Figure and vector of Axes.

# Examples
```julia
# Basic t-value heatmap with contour significance
result = permutation_test(prepared; n_permutations=1000, cluster_type=:temporal)
plot_tf_stats(result, channel_selection=channels(:Cz))

# Power difference with hide mask
plot_tf_stats(result, content=:difference, mask_style=:hide)

# Multiple channels in grid
plot_tf_stats(result, channel_selection=channels([:Cz, :Fz, :Pz, :Oz]))

# With analytic test results
result_analytic = analytic_test(prepared)
plot_tf_stats(result_analytic, mask_style=:alpha, colormap=:viridis)
```
"""
function plot_tf_stats(
    result::TFStatsResult;
    layout::Union{Symbol,PlotLayout} = :grid,
    layout_kwargs = (;),
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    content::Symbol = :tvalues,
    mask_style::Symbol = :alpha,
    colormap = nothing,
    colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
    xscale::Symbol = :linear,
    yscale::Symbol = :linear,
    colorbar_plot::Bool = true,
    plot_title::Union{Nothing,AbstractString} = nothing,
    show_title::Bool = true,
    significance_color = :black,
    significance_linewidth::Real = 2.0,
    figure_size::Union{Nothing,Tuple{Int,Int}} = nothing,
    display_plot::Bool = true,
    kwargs...,
)
    content in (:tvalues, :difference, :power_a, :power_b) ||
        @minimal_error("content must be :tvalues, :difference, :power_a, or :power_b, got :$content")
    mask_style in (:alpha, :contour, :hide, :none) ||
        @minimal_error("mask_style must be :alpha, :contour, :hide, or :none, got :$mask_style")

    use_xlog = xscale in (:log10, :log, log10)
    use_ylog = yscale in (:log10, :log, log10)
    show_cb = colorbar_plot
    eff_plot_title = plot_title
    eff_show_title = show_title

    all_electrodes = result.electrodes
    selected_channels = _select_tf_channels(all_electrodes, channel_selection, channel_plot_order)
    plot_channels = layout === :single ? [first(selected_channels)] : selected_channels
    frequencies = result.frequencies
    time_points = result.time_points

    ext_a = content in (:difference, :power_a) ? _tf_make_extractor(result.data[1].data_power, frequencies, time_points) : nothing
    ext_b = content in (:difference, :power_b) ? _tf_make_extractor(result.data[2].data_power, frequencies, time_points) : nothing

    channel_entries = Dict{Symbol,Tuple{Int,Matrix{Float64}}}()
    sizehint!(channel_entries, length(plot_channels))
    for ch_sym in plot_channels
        ch_idx = findfirst(==(ch_sym), all_electrodes)
        isnothing(ch_idx) && continue
        channel_entries[ch_sym] = (ch_idx, _extract_tf_stats_matrix(result, ch_idx, content; ext_a = ext_a, ext_b = ext_b))
    end
    isempty(channel_entries) && @minimal_error("No channels matched. Available: $(all_electrodes)")

    valid_plot_channels = [ch for ch in plot_channels if haskey(channel_entries, ch)]
    n_channels = length(valid_plot_channels)
    layout_grid_dims = get(layout_kwargs, :grid_dims, nothing)
    default_grid_dims = if !isnothing(layout_grid_dims)
        layout_grid_dims
    else
        n_c = min(n_channels, 3)
        (cld(n_channels, n_c), n_c)
    end

    eff_layout_kwargs = Dict{Symbol,Any}(
        :grid_dims => default_grid_dims,
        :grid_rowgap => get(layout_kwargs, :grid_rowgap, 10),
        :grid_colgap => get(layout_kwargs, :grid_colgap, 10),
        :topo_plot_width => get(layout_kwargs, :topo_plot_width, 0.1),
        :topo_plot_height => get(layout_kwargs, :topo_plot_height, 0.1),
    )

    eeg_layout = !isempty(result.data) && hasproperty(result.data[1], :layout) ? result.data[1].layout : nothing
    plot_layout = create_layout(layout, valid_plot_channels, eeg_layout; eff_layout_kwargs...)

    if isnothing(figure_size)
        figure_size = if plot_layout.type === :single
            (800, 500)
        elseif plot_layout.type === :topo
            (1200, 1000)
        else
            n_rows, n_cols = plot_layout.dims
            (n_cols * (show_cb ? 500 : 400), n_rows * 350)
        end
    end

    fig = _eegfun_figure(size = figure_size)
    cmap = _resolve_theme_colormap(fig, colormap)

    cr = if isnothing(colorrange)
        gmin, gmax = _tf_finite_extrema((entry[2] for entry in values(channel_entries)))
        _tf_auto_colorrange(gmin, gmax, cmap; force_symmetric = content in (:tvalues, :difference))
    else
        colorrange
    end

    axes, layout_channels = _apply_layout!(fig, plot_layout; xgrid = false, ygrid = false, xminorgrid = false, yminorgrid = false)

    if eff_show_title
        test_type = result isa TFAnalyticResult ? "Analytic" : "Permutation"
        content_label =
            content == :tvalues ? "t-values" :
            content == :difference ? "Power Difference" :
            content == :power_a ? result.data[1].condition_name : result.data[2].condition_name
        header_text = isnothing(eff_plot_title) ? "$test_type Test — $content_label" : String(eff_plot_title)
        if !isempty(header_text)
            if plot_layout.type === :grid
                _, n_cols = plot_layout.dims
                Label(fig[0, 1:n_cols], header_text, fontsize = 16, font = :bold)
            else
                Label(fig[0, 1], header_text, fontsize = 16, font = :bold)
            end
        end
    end

    last_hm = nothing

    for (ax, ch_sym) in zip(axes, layout_channels)
        ch_idx, data_mat = channel_entries[ch_sym]

        ax.xlabel = "Time (s)"
        ax.ylabel = "Frequency (Hz)"
        ax.title = eff_show_title ? string(ch_sym) : ""
        ax.xscale = use_xlog ? log10 : identity
        ax.yscale = use_ylog ? log10 : identity
        use_ylog && (ax.ytickformat = values -> [string(round(Int, v)) for v in values])

        sig_mask = mask_style != :none ? _extract_tf_significance_mask(result, ch_idx) : nothing

        if isnothing(sig_mask) || mask_style == :none
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(data_mat), colormap = cmap, colorrange = cr, nan_color = :transparent)
        elseif mask_style == :alpha
            last_hm = heatmap!(
                ax,
                time_points,
                frequencies,
                transpose(data_mat),
                colormap = cmap,
                colorrange = cr,
                alpha = 0.15,
                nan_color = :transparent,
            )
            if any(sig_mask)
                highlight = _tf_masked_copy(data_mat, sig_mask)
                last_hm =
                    heatmap!(ax, time_points, frequencies, transpose(highlight), colormap = cmap, colorrange = cr, nan_color = :transparent)
            end
        elseif mask_style == :hide
            highlight = _tf_masked_copy(data_mat, sig_mask)
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(highlight), colormap = cmap, colorrange = cr, nan_color = :transparent)
        elseif mask_style == :contour
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(data_mat), colormap = cmap, colorrange = cr, nan_color = :transparent)
            any(sig_mask) && _render_significance_contour!(
                ax,
                time_points,
                frequencies,
                sig_mask;
                color = significance_color,
                linewidth = significance_linewidth,
            )
        end
    end

    if plot_layout.type === :topo
        for ax in axes
            hidedecorations!(ax, grid = false, ticks = true, ticklabels = true)
            hidespines!(ax)
        end
    end

    length(axes) > 1 && linkaxes!(axes...)

    if show_cb && !isnothing(last_hm)
        cb_label = content == :tvalues ? "t-statistic" : content == :difference ? "Power Difference" : "Power"
        if plot_layout.type === :grid
            n_rows, n_cols = plot_layout.dims
            Colorbar(fig[1:n_rows, n_cols+1], last_hm, label = cb_label)
        elseif plot_layout.type === :topo
            Colorbar(fig[1, 2], last_hm, label = cb_label, height = Relative(0.5))
        else
            Colorbar(fig[1, 2], last_hm, label = cb_label)
        end
    end

    display_plot && _display_figure(fig)
    return (fig = fig, axes = axes)
end

# === Internal helpers ===

function _tf_masked_copy(data_mat::Matrix{Float64}, sig_mask::AbstractMatrix{Bool})
    out = similar(data_mat)
    @inbounds for i in eachindex(data_mat, sig_mask)
        out[i] = sig_mask[i] ? data_mat[i] : NaN
    end
    return out
end

function _extract_tf_stats_matrix(
    result::TFStatsResult,
    ch_idx::Int,
    content::Symbol;
    ext_a::Union{Nothing,_TFFastExtractor} = nothing,
    ext_b::Union{Nothing,_TFFastExtractor} = nothing,
)
    content == :tvalues && return result.stat_matrix.t[ch_idx, :, :]
    ch_sym = result.electrodes[ch_idx]
    if content == :difference
        ea = isnothing(ext_a) ? _tf_make_extractor(result.data[1].data_power, result.frequencies, result.time_points) : ext_a
        eb = isnothing(ext_b) ? _tf_make_extractor(result.data[2].data_power, result.frequencies, result.time_points) : ext_b
        mat = _tf_extract_matrix(ea, ch_sym)
        mat .-= _tf_extract_matrix(eb, ch_sym)
        return mat
    elseif content == :power_b
        eb = isnothing(ext_b) ? _tf_make_extractor(result.data[2].data_power, result.frequencies, result.time_points) : ext_b
        return _tf_extract_matrix(eb, ch_sym)
    else
        ea = isnothing(ext_a) ? _tf_make_extractor(result.data[1].data_power, result.frequencies, result.time_points) : ext_a
        return _tf_extract_matrix(ea, ch_sym)
    end
end

@views _extract_tf_significance_mask(result::TFStatsResult, ch_idx::Int) =
    result.masks.positive[ch_idx, :, :] .| result.masks.negative[ch_idx, :, :]

function _render_significance_contour!(
    ax::Axis,
    time_points::Vector{Float64},
    frequencies::Vector{Float64},
    sig_mask::AbstractMatrix;
    color = :black,
    linewidth::Real = 2.0,
)
    any(sig_mask) || return nothing
    sig_float = Matrix{Float64}(undef, size(sig_mask, 2), size(sig_mask, 1))
    @inbounds for j in axes(sig_mask, 2), i in axes(sig_mask, 1)
        sig_float[j, i] = Float64(sig_mask[i, j])
    end
    contour!(ax, time_points, frequencies, sig_float, levels = [0.5], color = color, linewidth = linewidth)
    return nothing
end
