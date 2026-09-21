"""
Plotting functions for Time-Frequency statistical test results.

Visualizes `TFClusterPermutationResult` and `TFAnalyticResult` as heatmaps
with configurable significance overlays.
"""

"""
    plot_tf_stats(result::TFStatsResult;
                  channel_selection::Function = channels(),
                  channel_plot_order::Union{Nothing, Vector{Symbol}} = nothing,
                  content::Symbol = :tvalues,
                  mask_style::Symbol = :alpha,
                  colormap = nothing,
                  colorrange = nothing,
                  ylogscale::Bool = false,
                  colorbar::Bool = true,
                  significance_color = :black,
                  significance_linewidth::Real = 2.0,
                  figure_size = nothing,
                  display_plot::Bool = true)

Plot TF statistical results as heatmaps with significance overlays.

Works with both `TFClusterPermutationResult` (from `permutation_test`) and
`TFAnalyticResult` (from `analytic_test`).

# Arguments
- `result::TFStatsResult`: Statistical result to plot
- `channel_selection::Function`: Channel selection predicate (default: all channels)
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
- `ylogscale::Bool`: Log scale for frequencies (default: false)
- `colorbar::Bool`: Show colorbar (default: true)
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
    channel_selection::Function = channels(),
    channel_plot_order::Union{Nothing,Vector{Symbol}} = nothing,
    content::Symbol = :tvalues,
    mask_style::Symbol = :alpha,
    colormap = nothing,
    colorrange::Union{Nothing,Tuple{Real,Real}} = nothing,
    ylogscale::Bool = false,
    colorbar::Bool = true,
    significance_color = :black,
    significance_linewidth::Real = 2.0,
    figure_size::Union{Nothing,Tuple{Int,Int}} = nothing,
    display_plot::Bool = true,
)
    # Validate arguments
    content in (:tvalues, :difference, :power_a, :power_b) ||
        error("content must be :tvalues, :difference, :power_a, or :power_b, got :$content")
    mask_style in (:alpha, :contour, :hide, :none) ||
        error("mask_style must be :alpha, :contour, :hide, or :none, got :$mask_style")

    # Select channels
    all_electrodes = result.electrodes
    selected_mask = channel_selection(all_electrodes)
    selected_channels = all_electrodes[selected_mask]
    isempty(selected_channels) && error("No channels matched. Available: $(all_electrodes)")

    # Apply user-specified plotting order if provided
    if !isnothing(channel_plot_order)
        ordered = Symbol[ch for ch in channel_plot_order if ch in selected_channels]
        if !isempty(ordered)
            selected_channels = ordered
        end
    end

    n_channels = length(selected_channels)
    frequencies = result.frequencies
    time_points = result.time_points

    # Grid layout: determine rows × cols
    n_cols = min(n_channels, 3)
    n_rows = ceil(Int, n_channels / n_cols)

    # Figure size: auto-scale based on grid
    if isnothing(figure_size)
        panel_w = colorbar ? 500 : 400
        panel_h = 350
        figure_size = (n_cols * panel_w, n_rows * panel_h)
    end

    fig = Figure(size = figure_size)

    # Build title
    test_type = isa(result, TFAnalyticResult) ? "Analytic" : "Permutation"
    content_label =
        content == :tvalues ? "t-values" :
        content == :difference ? "Power Difference" : content == :power_a ? result.data[1].condition_name : result.data[2].condition_name

    Label(fig[0, 1:n_cols], "$test_type Test — $content_label", fontsize = 16, font = :bold)

    axes = Axis[]
    local last_hm  # for colorbar

    for (ch_i, ch_sym) in enumerate(selected_channels)
        ch_idx = findfirst(==(ch_sym), all_electrodes)
        isnothing(ch_idx) && continue

        row = (ch_i - 1) ÷ n_cols + 1
        col = (ch_i - 1) % n_cols + 1

        ax = Axis(
            fig[row, col],
            xlabel = "Time (s)",
            ylabel = "Frequency (Hz)",
            title = string(ch_sym),
            yscale = ylogscale ? log10 : identity,
        )
        push!(axes, ax)

        if ylogscale
            ax.ytickformat = values -> [string(round(Int, v)) for v in values]
        end

        # Extract 2D data matrix [frequencies × time] for this channel
        data_mat = _extract_tf_stats_matrix(result, ch_idx, content)

        # Determine color range
        cr = colorrange
        if isnothing(cr)
            valid = filter(!isnan, vec(data_mat))
            if isempty(valid)
                cr = (-1.0, 1.0)
            elseif content == :tvalues
                # Symmetric range for t-values
                max_abs = maximum(abs, valid)
                cr = (-max_abs, max_abs)
            else
                cr = (minimum(valid), maximum(valid))
            end
        end

        # Plot heatmap based on mask_style
        has_mask = false
        local sig_mask = nothing
        if mask_style != :none
            sig_mask = _extract_tf_significance_mask(result, ch_idx)
            has_mask = !isnothing(sig_mask) && any(sig_mask)
        end
        
        cmap = _resolve_theme_colormap(ax, colormap)
        
        if !has_mask || mask_style == :none
            hm = heatmap!(ax, time_points, frequencies, data_mat', 
                          colormap = cmap, colorrange = cr, nan_color = :transparent)
            last_hm = hm
        elseif mask_style == :alpha
            heatmap!(ax, time_points, frequencies, data_mat', 
                     colormap = cmap, colorrange = cr, alpha = 0.15, nan_color = :transparent)
            
            highlight = copy(data_mat)
            highlight[.!sig_mask] .= NaN
            hm = heatmap!(ax, time_points, frequencies, highlight', 
                          colormap = cmap, colorrange = cr, nan_color = :transparent)
            last_hm = hm
        elseif mask_style == :hide
            highlight = copy(data_mat)
            highlight[.!sig_mask] .= NaN
            hm = heatmap!(ax, time_points, frequencies, highlight', 
                          colormap = cmap, colorrange = cr, nan_color = :transparent)
            last_hm = hm
        elseif mask_style == :contour
            hm = heatmap!(ax, time_points, frequencies, data_mat', 
                          colormap = cmap, colorrange = cr, nan_color = :transparent)
            last_hm = hm
            _render_significance_contour!(ax, time_points, frequencies, sig_mask; color = significance_color, linewidth = significance_linewidth)
        end
    end

    # Link axes for consistent zoom/pan
    length(axes) > 1 && linkaxes!(axes...)

    # Colorbar (on the right of the last column)
    if colorbar && @isdefined(last_hm)
        cb_label = content == :tvalues ? "t-statistic" : content == :difference ? "Power Difference" : "Power"
        Colorbar(fig[1:n_rows, n_cols+1], last_hm, label = cb_label)
    end

    if display_plot
        display(fig)
    end

    return (fig = fig, axes = axes)
end


# === Internal helpers ===

"""
    _extract_tf_stats_matrix(result, ch_idx, content) -> Matrix{Float64}

Extract a 2D matrix [frequencies × time] for a single channel from TF stats result.
"""
function _extract_tf_stats_matrix(result::TFStatsResult, ch_idx::Int, content::Symbol)
    if content == :tvalues
        return result.stat_matrix.t[ch_idx, :, :]
    else
        # Extract power from grand average TimeFreqData
        tf_data = content == :power_b ? result.data[2] : result.data[1]
        ch_sym = result.electrodes[ch_idx]

        # For :difference, compute A - B
        if content == :difference
            tf_a = result.data[1]
            tf_b = result.data[2]
            return _tf_power_matrix(tf_a, ch_sym, result.frequencies, result.time_points) .-
                   _tf_power_matrix(tf_b, ch_sym, result.frequencies, result.time_points)
        else
            return _tf_power_matrix(tf_data, ch_sym, result.frequencies, result.time_points)
        end
    end
end

"""
    _tf_power_matrix(tf_data, channel, frequencies, time_points) -> Matrix{Float64}

Extract a [frequencies × time] power matrix from a TimeFreqData's DataFrame.
Delegates to the shared `_tf_df_to_matrix` helper.
"""
function _tf_power_matrix(tf_data::TimeFreqData, channel::Symbol, frequencies::Vector{Float64}, time_points::Vector{Float64})
    return _tf_df_to_matrix(tf_data.data_power, channel, frequencies, time_points)
end

"""
    _extract_tf_significance_mask(result, ch_idx) -> BitMatrix

Extract a combined 2D significance mask [frequencies × time] for a single channel.
Combines positive and negative masks.
"""
function _extract_tf_significance_mask(result::TFStatsResult, ch_idx::Int)
    return result.masks.positive[ch_idx, :, :] .| result.masks.negative[ch_idx, :, :]
end


# === Significance overlay renderers ===

"""
    _render_significance_contour!(ax, time_points, frequencies, sig_mask; color, linewidth)

Render contour lines around significant regions on a TF heatmap.
"""
function _render_significance_contour!(
    ax::Axis,
    time_points::Vector{Float64},
    frequencies::Vector{Float64},
    sig_mask::AbstractMatrix;
    color = :black,
    linewidth::Real = 2.0,
)
    if !any(sig_mask)
        return
    end

    # Convert BitMatrix to Float64 for contour (Makie expects numeric)
    # sig_mask is [n_freqs × n_time], transpose for contour!(ax, x, y, z)
    sig_float = Float64.(sig_mask)'  # [n_time × n_freqs]
    contour!(ax, time_points, frequencies, sig_float, levels = [0.5], color = color, linewidth = linewidth)
end
