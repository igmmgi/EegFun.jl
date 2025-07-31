"""
    plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; sort_values=false, average_over=nothing)

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
- `average_over::Union{Symbol,Nothing}=nothing`: Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI across this grouping variable.

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    # Assuming summary_df has columns :channel, :kurtosis, :variance, :epoch
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, summary_df, :kurtosis, average_over=:epoch)
"""
function plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; sort_values=false, average_over=nothing)

    # Check if required columns exist
    if :channel ∉ propertynames(dat) || col ∉ propertynames(dat)
        @minimal_error("DataFrame must contain :channel and :$col columns.")
    end

    # If averaging is requested, compute mean and std
    if average_over !== nothing
        if average_over ∉ propertynames(dat)
            @minimal_error("Column :$average_over not found in DataFrame for averaging.")
        end
        
        # Group by channel and compute statistics
        grouped = groupby(dat, :channel)
        
        # Compute 95% confidence intervals around the mean
        summary_stats = combine(grouped, 
            col => mean => :mean, 
            col => std => :std,
            col => length => :n_samples
        )
        
        # Calculate margin of error for 95% confidence intervals
        summary_stats.margin_of_error = 1.96 .* (summary_stats.std ./ sqrt.(summary_stats.n_samples))
        
        # Use the averaged data for plotting
        plot_dat = summary_stats
        
        # Update y-axis label to indicate averaging and number of epochs
        # Count unique values in the averaging column to get number of epochs
        n_epochs = length(unique(dat[!, average_over]))
        ax.ylabel = "$(String(col)) (mean ± 95% CI over $(String(average_over)), n=$n_epochs epochs)"
    else
        # Use original data
        plot_dat = dat
        ax.ylabel = String(col)
    end

    # Optionally sort data
    if sort_values 
        sort_column = average_over !== nothing ? :mean : col
        plot_dat = sort(plot_dat, sort_column, rev=true)
    end
    
    # Extract plotting variables after sorting
    channel_names = String.(plot_dat.channel)
    if average_over !== nothing
        values_to_plot = plot_dat.mean
        margin_of_error = plot_dat.margin_of_error
    else
        values_to_plot = plot_dat[!, col]
        margin_of_error = nothing
    end

    # Configure the axis
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = pi / 4
    ax.xlabel = "Electrode"

    # Create the bar plot
    if average_over !== nothing
        # Plot with error bars and whiskers
        bars = barplot!(ax, 1:length(values_to_plot), values_to_plot)
        
        # Add error bars with whiskers
        for i in 1:length(values_to_plot)
            x_pos = i
            y_pos = values_to_plot[i]
            error = margin_of_error[i]
            
            # Draw vertical error bar line
            lines!(ax, [x_pos, x_pos], [y_pos - error, y_pos + error], 
                   color=:black, linewidth=2)
            
            # Draw whisker caps (horizontal lines at top and bottom)
            whisker_length = 0.1
            lines!(ax, [x_pos - whisker_length, x_pos + whisker_length], [y_pos - error, y_pos - error], 
                   color=:black, linewidth=2)
            lines!(ax, [x_pos - whisker_length, x_pos + whisker_length], [y_pos + error, y_pos + error], 
                   color=:black, linewidth=2)
        end
    else
        # Plot without error bars
        barplot!(ax, 1:length(values_to_plot), values_to_plot)
    end

    return nothing
end

"""
    plot_channel_summary(dat::DataFrame, col::Symbol; sort_values=false, average_over=nothing, display_plot::Bool = true)

Plot a bar chart summarizing a specific metric per channel from a DataFrame.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `col`: The column specified by the `col` argument, containing the values to plot.

# Arguments
- `dat::DataFrame`: DataFrame containing channel summary data.
- `col::Symbol`: The symbol representing the column in `dat` to plot on the y-axis.
- `sort_values::Bool=false`: If true, sort the bars by the values in the `col` column in descending order. Channel names on the x-axis will be reordered accordingly.
- `average_over::Union{Symbol,Nothing}=nothing`: Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI across this grouping variable.
- `display_plot::Bool=true`: Whether to display the plot.

# Returns
- `Figure`: The Makie Figure object.
- `Axis`: The Makie Axis object for the bar plot.

# Example
    # Assuming summary_df has columns :channel, :kurtosis, :variance, :epoch
    fig_kurt, ax_kurt = plot_channel_summary(summary_df, :kurtosis, average_over=:epoch)
    # display(fig_kurt)

    fig_var_sorted, ax_var_sorted = plot_channel_summary(summary_df, :variance, sort_values=true, average_over=:epoch)
    # display(fig_var_sorted)
"""
function plot_channel_summary(dat::DataFrame, col::Symbol; sort_values=false, average_over=nothing, display_plot::Bool = true)
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    
    # Use the mutating version to plot
    plot_channel_summary!(fig, ax, dat, col; sort_values=sort_values, average_over=average_over)

    if display_plot
        display_figure(fig)
    end

    return fig, ax

end
