
"""
    plot_joint_probability!(fig::Figure, ax::Axis, dat::DataFrame)

Plot a bar chart of joint probability values per channel on the provided figure and axis.

This is the mutating version that plots directly on the provided `fig` and `ax` objects.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `:jp`: Containing the joint probability values to plot.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `ax::Axis`: The Makie Axis object to plot on
- `dat::DataFrame`: DataFrame with channel and joint probability data.

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    # Assuming jp_df has columns :channel and :jp
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_joint_probability!(fig, ax, jp_df)
"""
function plot_joint_probability!(fig::Figure, ax::Axis, dat::DataFrame)
    # Basic validation
    if :channel ∉ propertynames(dat) || :jp ∉ propertynames(dat)
        error("DataFrame must contain :channel and :jp columns.")
    end
    channel_names = String.(dat.channel)

    # Configure the axis
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = pi / 4
    ax.xlabel = "Electrode"
    ax.ylabel = "Joint Probability"

    # Create the bar plot
    barplot!(ax, 1:nrow(dat), dat[!, :jp])

    return nothing
end

"""
    plot_joint_probability(dat::DataFrame)

Plot a bar chart of joint probability values per channel.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `:jp`: Containing the joint probability values to plot.

# Arguments
- `dat::DataFrame`: DataFrame with channel and joint probability data.

# Returns
- `Figure`: The Makie Figure object.
- `Axis`: The Makie Axis object for the bar plot.

# Example
    # Assuming jp_df has columns :channel and :jp
    fig, ax = plot_joint_probability(jp_df)
    # display(fig) # Display if needed
"""
function plot_joint_probability(dat::DataFrame, display_plot::Bool = true)
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Use the mutating version to plot
    plot_joint_probability!(fig, ax, dat)

    if display_plot
        display_figure(fig)
    end

    return fig, ax

end
