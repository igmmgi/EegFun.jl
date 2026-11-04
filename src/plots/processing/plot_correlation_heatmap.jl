# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_CORRELATION_HEATMAP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :mask_range => (nothing, "Optional tuple (min_val, max_val) to mask correlations within this range"),
    :display_plot => (true, "Whether to display the plot"),
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :plot_title_fontsize => (16, "Font size for the title"),
    :colormap => (nothing, "Colormap for the heatmap"),
    :colorrange => ((-1, 1), "Color range for the heatmap and colorbar (should be -1 to 1 for correlations)"),
    :nan_color => (:transparent, "Color for NaN values"),
    :xlabel => ("", "Label for x-axis"),
    :ylabel => ("", "Label for y-axis"),
    :label_fontsize => (14, "Font size for axis labels"),
    :xtick_rotation => (π / 4, "Rotation angle for x-axis tick labels"),
    :ytick_rotation => (0, "Rotation angle for y-axis tick labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    :colorbar_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Colorbar block (see Makie's [Colorbar documentation](https://docs.makie.org/stable/reference/blocks/colorbar/) for available attributes).",
    ),
    :colorbar_plot => (true, "Whether to display the colorbar"),
    :colorbar_position => (:right, "Position of the colorbar (:right, :left, :top, :bottom, or tuple)"),
    :colorbar_label => ("Correlation", "Label for the colorbar"),
)

function _fill_corr_row!(mat::Matrix{Float32}, j::Int, col::AbstractVector{<:Real})
    @inbounds for i in eachindex(col)
        mat[j, i] = Float32(col[i])
    end
    return nothing
end

function _fill_corr_row_masked!(mat::Matrix{Float32}, j::Int, col::AbstractVector{<:Real}, lo::Float64, hi::Float64)
    @inbounds for i in eachindex(col)
        v = Float64(col[i])
        mat[j, i] = (lo <= v <= hi) ? NaN32 : Float32(v)
    end
    return nothing
end

function _build_transposed_corr_matrix(corr_df::DataFrame, mask_range)::Matrix{Float32}
    n_rows = nrow(corr_df)
    n_cols = ncol(corr_df) - 1
    mat = Matrix{Float32}(undef, n_cols, n_rows)
    if isnothing(mask_range)
        for j = 1:n_cols
            _fill_corr_row!(mat, j, corr_df[!, j+1])
        end
    else
        lo, hi = Float64(mask_range[1]), Float64(mask_range[2])
        for j = 1:n_cols
            _fill_corr_row_masked!(mat, j, corr_df[!, j+1], lo, hi)
        end
    end
    return mat
end

function _plot_correlation_heatmap_impl!(ax::Axis, corr_df::DataFrame, plot_kwargs::Dict{Symbol,Any})
    corr_t = _build_transposed_corr_matrix(corr_df, plot_kwargs[:mask_range])

    colorrange = plot_kwargs[:colorrange]
    if colorrange[1] < -1 || colorrange[2] > 1
        @minimal_warning("Color range ($(colorrange[1]), $(colorrange[2])) extends beyond typical correlation bounds (-1, 1)")
    end

    row_names = String.(corr_df[!, 1])
    col_names = names(corr_df)[2:end]

    ax.xlabel = plot_kwargs[:xlabel]
    ax.ylabel = plot_kwargs[:ylabel]
    ax.xticks = (1:length(col_names), col_names)
    ax.yticks = (1:length(row_names), row_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.yticklabelrotation = plot_kwargs[:ytick_rotation]
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    plot_title = plot_kwargs[:plot_title]
    if plot_kwargs[:show_title] && !isnothing(plot_title) && !isempty(plot_title)
        ax.title = plot_title
        ax.titlesize = plot_kwargs[:plot_title_fontsize]
    end

    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    actual_cmap = _resolve_theme_colormap(ax, plot_kwargs[:colormap])
    heatmap!(ax, corr_t, colormap = actual_cmap, colorrange = colorrange, nan_color = plot_kwargs[:nan_color])

    if plot_kwargs[:colorbar_plot]
        cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
        sg = ax.layoutobservables.gridcontent[].parent
        Colorbar(sg[cb_pos...]; colormap = actual_cmap, limits = colorrange, plot_kwargs[:colorbar_kwargs]...)
    end
    return ax
end

"""
    plot_correlation_heatmap!(fig::Figure, ax::Axis, corr_df::DataFrame; kwargs...) -> Axis
    plot_correlation_heatmap(corr_df::DataFrame; kwargs...) -> NamedTuple

Plot a correlation matrix `DataFrame` (first column row labels, remaining columns correlation values) as a heatmap.

$(_generate_kwargs_doc(PLOT_CORRELATION_HEATMAP_KWARGS))
"""
function plot_correlation_heatmap!(fig::Figure, ax::Axis, corr_df::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_CORRELATION_HEATMAP_KWARGS, _normalize_plot_kwargs(kwargs))
    return _plot_correlation_heatmap_impl!(ax, corr_df, plot_kwargs)
end

"""
    plot_correlation_heatmap(corr_df::DataFrame; kwargs...) -> NamedTuple

Create a figure and plot a correlation matrix `DataFrame` as a heatmap. See [`plot_correlation_heatmap!`](@ref) for keyword options.
"""
function plot_correlation_heatmap(corr_df::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_CORRELATION_HEATMAP_KWARGS, _normalize_plot_kwargs(kwargs))
    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])
    _plot_correlation_heatmap_impl!(ax, corr_df, plot_kwargs)
    if plot_kwargs[:display_plot]
        _set_window_title("Correlation Heatmap")
        _display_figure(fig)
    end
    return (fig = fig, axes = [ax])
end
