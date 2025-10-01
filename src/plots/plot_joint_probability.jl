# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_JOINT_PROBABILITY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Bar styling
    :bar_color => (:steelblue, "Color of the bars"),
    :bar_width => (0.8, "Width of bars"),
    :bar_alpha => (0.7, "Transparency of bars"),

    # Axis styling
    :xlabel => ("Electrode", "Label for x-axis"),
    :ylabel => ("Joint Probability", "Label for y-axis"),
    :title => ("", "Plot title"),
    :xtick_rotation => (π/4, "Rotation angle for x-axis tick labels"),

    # Font sizes
    :title_fontsize => (16, "Font size for title"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),

    # Grid styling
    :grid_visible => (true, "Whether to show grid"),
    :grid_alpha => (0.3, "Transparency of grid"),

    # Sorting
    :sort_values => (false, "If true, sort the bars by joint probability values in descending order"),
)

"""
    plot_joint_probability!(fig::Figure, ax::Axis, dat::DataFrame; kwargs...)

Plot a bar chart of joint probability values per channel on the provided figure and axis.

This is the mutating version that plots directly on the provided `fig` and `ax` objects.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `:jp`: Containing the joint probability values to plot.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `ax::Axis`: The Makie Axis object to plot on
- `dat::DataFrame`: DataFrame with channel and joint probability data.

$(generate_kwargs_doc(PLOT_JOINT_PROBABILITY_KWARGS))

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    # Basic usage
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_joint_probability!(fig, ax, jp_df)
    
    # Custom styling
    plot_joint_probability!(fig, ax, jp_df;
        title = "Custom Joint Probability",
        bar_color = :orange,
        sort_values = true)
"""
function plot_joint_probability!(fig::Figure, ax::Axis, dat::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_JOINT_PROBABILITY_KWARGS, kwargs)

    # Basic validation
    if :channel ∉ propertynames(dat) || :jp ∉ propertynames(dat)
        @minimal_error("DataFrame must contain :channel and :jp columns.")
    end

    # Optionally sort data
    if plot_kwargs[:sort_values]
        dat = sort(dat, :jp, rev = true)
    end

    channel_names = String.(dat.channel)

    # Configure the axis
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.xlabel = plot_kwargs[:xlabel]
    ax.ylabel = plot_kwargs[:ylabel]
    ax.title = plot_kwargs[:title]
    ax.titlesize = plot_kwargs[:title_fontsize]
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    # Configure grid
    ax.xgridvisible = plot_kwargs[:grid_visible]
    ax.ygridvisible = plot_kwargs[:grid_visible]
    ax.xgridwidth = 1
    ax.ygridwidth = 1
    ax.xgridcolor = (:gray, plot_kwargs[:grid_alpha])
    ax.ygridcolor = (:gray, plot_kwargs[:grid_alpha])

    # Create the bar plot
    barplot!(
        ax,
        1:nrow(dat),
        dat[!, :jp],
        color = plot_kwargs[:bar_color],
        width = plot_kwargs[:bar_width],
        alpha = plot_kwargs[:bar_alpha],
    )

    return nothing
end

"""
    plot_joint_probability(dat::DataFrame; kwargs...)

Plot a bar chart of joint probability values per channel.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `:jp`: Containing the joint probability values to plot.

# Arguments
- `dat::DataFrame`: DataFrame with channel and joint probability data.

$(generate_kwargs_doc(PLOT_JOINT_PROBABILITY_KWARGS))

# Returns
- `Figure`: The Makie Figure object.
- `Axis`: The Makie Axis object for the bar plot.

# Example
    # Basic usage
    fig, ax = plot_joint_probability(jp_df)
    
    # Custom styling
    fig, ax = plot_joint_probability(jp_df;
        title = "Joint Probability Analysis",
        bar_color = :red,
        sort_values = true,
        grid_visible = false)
"""
function plot_joint_probability(dat::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_JOINT_PROBABILITY_KWARGS, kwargs)

    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Use the mutating version to plot
    plot_joint_probability!(fig, ax, dat; kwargs...)

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    return fig, ax
end
