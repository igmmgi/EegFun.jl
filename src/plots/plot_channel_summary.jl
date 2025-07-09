"""
    plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; sort_values=false)

Plot a bar chart summarizing a specific metric per channel from a DataFrame on the provided figure and axis.

This is the mutating version that plots directly on the provided `fig` and `ax` objects.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `col`: The column specified by the `col` argument, containing the values to plot.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `ax::Axis`: The Makie Axis object to plot on
- `dat::DataFrame`: DataFrame containing channel summary data.
- `col::Symbol`: The symbol representing the column in `dat` to plot on the y-axis.
- `sort_values::Bool=false`: If true, sort the bars by the values in the `col` column in descending order. Channel names on the x-axis will be reordered accordingly.

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    # Assuming summary_df has columns :channel, :kurtosis, :variance
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, summary_df, :kurtosis)
"""
function plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; sort_values=false)
    # Check if required columns exist
    if :channel ∉ propertynames(dat) || col ∉ propertynames(dat)
        error("DataFrame must contain :channel and :$col columns.")
    end

    # Optionally sort data
    if sort_values 
        dat = sort(dat, col, rev=true)
    end

    channel_names = String.(plot_dat.channel)
    values_to_plot = plot_dat[!, col]

    # Configure the axis
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = pi / 4
    ax.xlabel = "Electrode"
    ax.ylabel = String(col) # Use column name directly

    # Create the bar plot
    barplot!(ax, 1:nrow(plot_dat), values_to_plot)

    return nothing
end

"""
    plot_channel_summary(dat::DataFrame, col::Symbol; sort_values=false)

Plot a bar chart summarizing a specific metric per channel from a DataFrame.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `col`: The column specified by the `col` argument, containing the values to plot.

# Arguments
- `dat::DataFrame`: DataFrame containing channel summary data.
- `col::Symbol`: The symbol representing the column in `dat` to plot on the y-axis.
- `sort_values::Bool=false`: If true, sort the bars by the values in the `col` column in descending order. Channel names on the x-axis will be reordered accordingly.

# Returns
- `Figure`: The Makie Figure object.
- `Axis`: The Makie Axis object for the bar plot.

# Example
    # Assuming summary_df has columns :channel, :kurtosis, :variance
    fig_kurt, ax_kurt = plot_channel_summary(summary_df, :kurtosis)
    # display(fig_kurt)

    fig_var_sorted, ax_var_sorted = plot_channel_summary(summary_df, :variance, sort_values=true)
    # display(fig_var_sorted)
"""
function plot_channel_summary(dat::DataFrame, col::Symbol; sort_values=false, display_plot::Bool = true)
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    
    # Use the mutating version to plot
    plot_channel_summary!(fig, ax, dat, col; sort_values=sort_values)

    if display_plot
        display_figure(fig)
    end

    return fig, ax

end
