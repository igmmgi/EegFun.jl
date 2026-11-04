# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_JOINT_PROBABILITY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Whether to display the plot"),
    :bar_color => (:steelblue, "Color of the bars"),
    :bar_width => (0.8, "Width of bars"),
    :bar_alpha => (0.7, "Transparency of bars"),
    :xlabel => ("Electrode", "Label for x-axis"),
    :ylabel => ("Joint Probability", "Label for y-axis"),
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :xtick_rotation => (π / 4, "Rotation angle for x-axis tick labels"),
    :plot_title_fontsize => (16, "Font size for title"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :xgrid => (true, "Whether to show x-axis grid"),
    :ygrid => (true, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    :sort_values => (false, "If true, sort the bars by joint probability values in descending order"),
)

function _plot_joint_probability_impl!(ax::Axis, dat::DataFrame, plot_kwargs::Dict{Symbol,Any})
    (:channel in propertynames(dat) && :jp in propertynames(dat)) || @minimal_error("DataFrame must contain :channel and :jp columns.")

    raw_ch = dat[!, :channel]
    jp_vals = dat[!, :jp]
    if plot_kwargs[:sort_values]
        perm = sortperm(jp_vals; rev = true)
        channel_names = String.(raw_ch[perm])
        jp_vals = jp_vals[perm]
    else
        channel_names = String.(raw_ch)
    end

    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.xlabel = plot_kwargs[:xlabel]
    ax.ylabel = plot_kwargs[:ylabel]
    if plot_kwargs[:show_title] && !isnothing(plot_kwargs[:plot_title])
        ax.title = plot_kwargs[:plot_title]
    end
    ax.titlesize = plot_kwargs[:plot_title_fontsize]
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

    barplot!(
        ax,
        1:length(jp_vals),
        jp_vals,
        color = plot_kwargs[:bar_color],
        width = plot_kwargs[:bar_width],
        alpha = plot_kwargs[:bar_alpha],
    )
    return ax
end

"""
    plot_joint_probability!(fig::Figure, ax::Axis, dat::DataFrame; kwargs...) -> Axis
    plot_joint_probability(dat::DataFrame; kwargs...) -> NamedTuple

Plot per-channel joint probability (`:jp`) bar charts from `dat`.

$(_generate_kwargs_doc(PLOT_JOINT_PROBABILITY_KWARGS))
"""
function plot_joint_probability!(fig::Figure, ax::Axis, dat::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_JOINT_PROBABILITY_KWARGS, kwargs)
    return _plot_joint_probability_impl!(ax, dat, plot_kwargs)
end

"""
    plot_joint_probability(dat::DataFrame; kwargs...) -> NamedTuple

Create a figure and plot per-channel joint probability (`:jp`) bar charts. See [`plot_joint_probability!`](@ref) for keyword options.
"""
function plot_joint_probability(dat::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_JOINT_PROBABILITY_KWARGS, kwargs)
    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])
    _plot_joint_probability_impl!(ax, dat, plot_kwargs)
    _set_window_title("Joint Probability")
    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end
