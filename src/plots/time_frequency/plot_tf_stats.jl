# Plotting functions for Time-Frequency statistical test results.

"""
    plot_tf_stats(result::TFStatsResult;
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

Plot TF statistical results as heatmaps with significance overlays.

Works with both `TFClusterPermutationResult` (from `permutation_test`) and
`TFAnalyticResult` (from `analytic_test`).

# Arguments
- `result::TFStatsResult`: Statistical result to plot
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
    frequencies = result.frequencies
    time_points = result.time_points

    ext_a = content in (:difference, :power_a) ? _tf_make_extractor(result.data[1].data_power, frequencies, time_points) : nothing
    ext_b = content in (:difference, :power_b) ? _tf_make_extractor(result.data[2].data_power, frequencies, time_points) : nothing

    channel_entries = Tuple{Symbol,Int,Matrix{Float64}}[]
    sizehint!(channel_entries, length(selected_channels))
    for ch_sym in selected_channels
        ch_idx = findfirst(==(ch_sym), all_electrodes)
        isnothing(ch_idx) && continue
        push!(channel_entries, (ch_sym, ch_idx, _extract_tf_stats_matrix(result, ch_idx, content; ext_a = ext_a, ext_b = ext_b)))
    end
    isempty(channel_entries) && @minimal_error("No channels matched. Available: $(all_electrodes)")

    n_channels = length(channel_entries)
    n_cols = min(n_channels, 3)
    n_rows = cld(n_channels, n_cols)

    if isnothing(figure_size)
        figure_size = (n_cols * (show_cb ? 500 : 400), n_rows * 350)
    end

    fig = _eegfun_figure(size = figure_size)
    cmap = _resolve_theme_colormap(fig, colormap)

    cr = if isnothing(colorrange)
        gmin, gmax = _tf_finite_extrema((entry[3] for entry in channel_entries))
        _tf_auto_colorrange(gmin, gmax, cmap; force_symmetric = content in (:tvalues, :difference))
    else
        colorrange
    end

    if eff_show_title
        test_type = result isa TFAnalyticResult ? "Analytic" : "Permutation"
        content_label =
            content == :tvalues ? "t-values" :
            content == :difference ? "Power Difference" :
            content == :power_a ? result.data[1].condition_name : result.data[2].condition_name
        header_text = isnothing(eff_plot_title) ? "$test_type Test — $content_label" : String(eff_plot_title)
        !isempty(header_text) && Label(fig[0, 1:n_cols], header_text, fontsize = 16, font = :bold)
    end

    axes = Vector{Axis}(undef, n_channels)
    last_hm = nothing

    for (ch_i, (ch_sym, ch_idx, data_mat)) in enumerate(channel_entries)
        row = fld(ch_i - 1, n_cols) + 1
        col = mod(ch_i - 1, n_cols) + 1

        ax = Axis(
            fig[row, col],
            xlabel = "Time (s)",
            ylabel = "Frequency (Hz)",
            title = eff_show_title ? string(ch_sym) : "",
            xscale = use_xlog ? log10 : identity,
            yscale = use_ylog ? log10 : identity,
        )
        axes[ch_i] = ax
        use_ylog && (ax.ytickformat = values -> [string(round(Int, v)) for v in values])

        sig_mask = mask_style != :none ? _extract_tf_significance_mask(result, ch_idx) : nothing
        has_mask = !isnothing(sig_mask) && any(sig_mask)

        if !has_mask || mask_style == :none
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(data_mat), colormap = cmap, colorrange = cr, nan_color = :transparent)
        elseif mask_style == :alpha
            heatmap!(
                ax,
                time_points,
                frequencies,
                transpose(data_mat),
                colormap = cmap,
                colorrange = cr,
                alpha = 0.15,
                nan_color = :transparent,
            )
            highlight = _tf_masked_copy(data_mat, sig_mask)
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(highlight), colormap = cmap, colorrange = cr, nan_color = :transparent)
        elseif mask_style == :hide
            highlight = _tf_masked_copy(data_mat, sig_mask)
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(highlight), colormap = cmap, colorrange = cr, nan_color = :transparent)
        elseif mask_style == :contour
            last_hm =
                heatmap!(ax, time_points, frequencies, transpose(data_mat), colormap = cmap, colorrange = cr, nan_color = :transparent)
            _render_significance_contour!(
                ax,
                time_points,
                frequencies,
                sig_mask;
                color = significance_color,
                linewidth = significance_linewidth,
            )
        end
    end

    length(axes) > 1 && linkaxes!(axes...)

    if show_cb && !isnothing(last_hm)
        cb_label = content == :tvalues ? "t-statistic" : content == :difference ? "Power Difference" : "Power"
        Colorbar(fig[1:n_rows, n_cols+1], last_hm, label = cb_label)
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
