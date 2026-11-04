# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_CHANNEL_SUMMARY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Display plot."),
    :sort_values => (false, "Sort the bars by the values in the `col` column in descending order."),
    :average_over => (nothing, "Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI."),
    :dims => (nothing, "Tuple (rows, cols) for grid layout. default = _best_rect(n_columns)."),
    :bar_color => (:steelblue, "Bar color."),
    :bar_width => (0.8, "Bar width."),
    :bar_alpha => (0.7, "Bar alpha."),
    :error_color => (:black, "Error bar color."),
    :error_linewidth => (2, "Error bar linewidth."),
    :xlabel => ("Electrode", "x-axis label."),
    :ylabel => ("", "y-axis label"),
    :plot_title => (nothing, "Plot title."),
    :show_title => (true, "Whether to show the plot title."),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :plot_title_fontsize => (16, "Font size for title."),
    :label_fontsize => (14, "Font size for axis labels."),
    :tick_fontsize => (12, "Font size for tick labels."),
    :xtick_rotation => (π / 4, "Rotation angle for x-axis tick labels."),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
)

function _configure_channel_summary_axis!(
    ax::Axis,
    channel_names::Vector{String},
    col::Symbol,
    plot_kwargs::Dict{Symbol,Any},
    n_epochs::Union{Int,Nothing},
)
    ax.ylabel =
        plot_kwargs[:ylabel] != "" ? plot_kwargs[:ylabel] :
        (!isnothing(plot_kwargs[:average_over]) ? "$(String(col)) (± 95% CI n=$n_epochs)" : String(col))
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.xlabel = plot_kwargs[:xlabel]
    if plot_kwargs[:show_title] && !isnothing(plot_kwargs[:plot_title])
        ax.title = plot_kwargs[:plot_title]
        ax.titlesize = plot_kwargs[:plot_title_fontsize]
    end
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]
    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )
end

function _plot_channel_summary_impl!(ax::Axis, dat::DataFrame, col::Symbol, channel_selection, plot_kwargs::Dict{Symbol,Any})
    (:channel in propertynames(dat) && col in propertynames(dat)) || @minimal_error("DataFrame must contain :channel and :$col columns.")

    raw_channels = dat[!, :channel]
    chan_syms = eltype(raw_channels) === Symbol ? raw_channels : Symbol.(raw_channels)
    mask = _to_channel_func(channel_selection)(chan_syms)
    all_selected = all(mask)

    avg_col = plot_kwargs[:average_over]
    n_epochs = nothing
    margin_of_error = nothing

    if !isnothing(avg_col)
        avg_col in propertynames(dat) || @minimal_error("Column :$avg_col not found in DataFrame for averaging.")
        sub_dat = all_selected ? dat : dat[mask, :]
        n_epochs = length(unique(sub_dat[!, avg_col]))
        grouped = groupby(sub_dat, :channel)
        comb = combine(grouped, col => mean => :mean, col => std => :std, col => length => :n_samples)
        channel_names = String.(comb[!, :channel])
        values_to_plot = comb[!, :mean]
        margin_of_error = 1.96 .* (comb[!, :std] ./ sqrt.(comb[!, :n_samples]))
    else
        ch_sel = all_selected ? raw_channels : raw_channels[mask]
        channel_names = String.(ch_sel)
        raw_vals = dat[!, col]
        values_to_plot = all_selected ? raw_vals : raw_vals[mask]
    end

    if plot_kwargs[:sort_values]
        perm = sortperm(values_to_plot; rev = true)
        channel_names = channel_names[perm]
        values_to_plot = values_to_plot[perm]
        !isnothing(margin_of_error) && (margin_of_error = margin_of_error[perm])
    end

    _configure_channel_summary_axis!(ax, channel_names, col, plot_kwargs, n_epochs)

    xs = 1:length(values_to_plot)
    barplot!(ax, xs, values_to_plot, color = plot_kwargs[:bar_color], width = plot_kwargs[:bar_width], alpha = plot_kwargs[:bar_alpha])
    if !isnothing(margin_of_error)
        errorbars!(ax, xs, values_to_plot, margin_of_error, color = plot_kwargs[:error_color], linewidth = plot_kwargs[:error_linewidth])
    end
    return ax
end

function _plot_multiple_columns!(
    fig::Figure,
    dat::DataFrame,
    cols::AbstractVector{Symbol},
    plot_kwargs::Dict{Symbol,Any};
    channel_selection = channels(),
)
    n_cols = length(cols)
    axes = Axis[]
    n_cols == 0 && return axes

    if isnothing(plot_kwargs[:dims])
        rs, cs = _best_rect(n_cols)
    else
        dims = plot_kwargs[:dims]
        (dims isa Tuple{Int,Int} && dims[1] > 0 && dims[2] > 0) ||
            @minimal_error("dims must be a tuple of two positive integers (rows, cols), got: $dims")
        rs, cs = dims
    end

    idx = 1
    for r = 1:rs, c = 1:cs
        ax = Axis(fig[r, c])
        _plot_channel_summary_impl!(ax, dat, cols[idx], channel_selection, plot_kwargs)
        push!(axes, ax)
        idx += 1
        idx > n_cols && return axes
    end
    return axes
end

"""
    plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col; channel_selection=channels(), kwargs...) -> Axis
    plot_channel_summary(dat::DataFrame, col; channel_selection=channels(), kwargs...) -> NamedTuple

Plot per-channel summary bar charts (with optional `average_over` 95% CI error bars and multi-column grids) from a summary `DataFrame`.

$(_generate_kwargs_doc(PLOT_CHANNEL_SUMMARY_KWARGS))
"""
function plot_channel_summary!(
    fig::Figure,
    ax::Axis,
    dat::DataFrame,
    col::Union{Symbol,AbstractString};
    channel_selection = channels(),
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)
    return _plot_channel_summary_impl!(ax, dat, Symbol(col), channel_selection, plot_kwargs)
end

"""
    plot_channel_summary(dat::DataFrame, col; channel_selection=channels(), kwargs...) -> NamedTuple

Create a figure and plot per-channel summary bar charts (single or multi-column). See [`plot_channel_summary!`](@ref) for keyword options.
"""
function plot_channel_summary(dat::DataFrame, col::Union{Symbol,AbstractString}; channel_selection = channels(), kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)
    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])
    _plot_channel_summary_impl!(ax, dat, Symbol(col), channel_selection, plot_kwargs)
    if plot_kwargs[:display_plot]
        _set_window_title("Channel Summary")
        _display_figure(fig)
    end
    return (fig = fig, axes = [ax])
end

plot_channel_summary(dat::AbstractVector{<:DataFrame}, col::Union{Symbol,AbstractString}; kwargs...) =
    plot_channel_summary.(dat, Symbol(col); kwargs...)

function plot_channel_summary(
    dat::DataFrame,
    col::AbstractVector{<:Union{Symbol,AbstractString}};
    channel_selection = channels(),
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)
    fig = _eegfun_figure()
    axes = _plot_multiple_columns!(fig, dat, Symbol.(col), plot_kwargs; channel_selection = channel_selection)
    if plot_kwargs[:display_plot]
        _set_window_title("Channel Summary")
        _display_figure(fig)
    end
    return (fig = fig, axes = axes)
end

plot_channel_summary(dat::AbstractVector{<:DataFrame}, col::AbstractVector{<:Union{Symbol,AbstractString}}; kwargs...) =
    plot_channel_summary.(dat, Ref(Symbol.(col)); kwargs...)
