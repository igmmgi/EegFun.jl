# Default parameters for channel summary plots
const DEFAULT_CHANNEL_SUMMARY_KWARGS = Dict(
    :sort_values => false,
    :average_over => nothing,
    :display_plot => true,
    :bar_color => :steelblue,
    :error_color => :black,
    :error_linewidth => 2,
    :xtick_rotation => π/4,
    :xlabel => "Electrode",
    :bar_width => 0.8,
    :bar_alpha => 0.7,
    :grid_visible => true,
    :grid_alpha => 0.3,
    :title => "",
    :title_fontsize => 16,
    :label_fontsize => 14,
    :tick_fontsize => 12,
)

"""
    plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; kwargs...)

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
- `kwargs...`: Keyword arguments for customization (see below for available options).

# Available Keyword Arguments
- `sort_values::Bool=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:sort_values])`: If true, sort the bars by the values in the `col` column in descending order.
- `average_over::Union{Symbol,Nothing}=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:average_over])`: Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI.
- `bar_color::Symbol=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:bar_color])`: Color of the bars.
- `error_color::Symbol=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:error_color])`: Color of error bars.
- `error_linewidth::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:error_linewidth])`: Line width of error bars.
- `xtick_rotation::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:xtick_rotation])`: Rotation angle for x-axis tick labels.
- `xlabel::String="$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:xlabel])"`: Label for x-axis.
- `title::String="$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:title])"`: Plot title.
- `title_fontsize::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:title_fontsize])`: Font size for title.
- `label_fontsize::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:label_fontsize])`: Font size for axis labels.
- `tick_fontsize::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:tick_fontsize])`: Font size for tick labels.
- `bar_width::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:bar_width])`: Width of bars.
- `bar_alpha::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:bar_alpha])`: Transparency of bars.
- `grid_visible::Bool=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:grid_visible])`: Whether to show grid.
- `grid_alpha::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:grid_alpha])`: Transparency of grid.

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    # Assuming summary_df has columns :channel, :kurtosis, :variance, :epoch
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, summary_df, :kurtosis, average_over=:epoch, title="Kurtosis by Channel")
"""
function plot_channel_summary!(
    fig::Figure,
    ax::Axis,
    dat::DataFrame,
    col::Symbol;
    kwargs...
)
    # Merge user kwargs with defaults
    plot_kwargs = merge(DEFAULT_CHANNEL_SUMMARY_KWARGS, Dict(kwargs))

    # Check if required columns exist
    if :channel ∉ propertynames(dat) || col ∉ propertynames(dat)
        @minimal_error("DataFrame must contain :channel and :$col columns.")
    end

    # If averaging is requested, compute mean and std
    n_epochs = nothing
    if plot_kwargs[:average_over] !== nothing

        if plot_kwargs[:average_over] ∉ propertynames(dat)
            @minimal_error("Column :$(plot_kwargs[:average_over]) not found in DataFrame for averaging.")
        end

        # Count unique values in the averaging column to get number of epochs
        n_epochs = length(unique(dat[!, plot_kwargs[:average_over]]))

        # Group by channel and compute statistics
        grouped = groupby(dat, :channel)
        dat = combine(grouped, col => mean => :mean, col => std => :std, col => length => :n_samples)
        dat.margin_of_error = 1.96 .* (dat.std ./ sqrt.(dat.n_samples))
    end

    # Optionally sort data
    if plot_kwargs[:sort_values]
        sort_column = plot_kwargs[:average_over] !== nothing ? :mean : col
        dat = sort(dat, sort_column, rev = true)
    end

    # Extract plotting variables after sorting
    channel_names = String.(dat.channel)
    if plot_kwargs[:average_over] !== nothing
        values_to_plot = dat.mean
        margin_of_error = dat.margin_of_error
    else
        values_to_plot = dat[!, col]
        margin_of_error = nothing
    end

    # Configure the axis
    ax.ylabel = plot_kwargs[:average_over] !== nothing ? "$(String(col)) (± 95% CI n=$n_epochs)" : "$(String(col))"
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.xlabel = plot_kwargs[:xlabel]
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
    barplot!(ax, 1:length(values_to_plot), values_to_plot, 
             color = plot_kwargs[:bar_color], 
             width = plot_kwargs[:bar_width],
             alpha = plot_kwargs[:bar_alpha])
    
    if plot_kwargs[:average_over] !== nothing # add error bars
        errorbars!(ax, 1:length(values_to_plot), values_to_plot, margin_of_error, 
                  color = plot_kwargs[:error_color], 
                  linewidth = plot_kwargs[:error_linewidth])
    end

    return nothing
end

"""
    plot_channel_summary(dat::DataFrame, col::Symbol; kwargs...)

Plot a bar chart summarizing a specific metric per channel from a DataFrame.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `col`: The column specified by the `col` argument, containing the values to plot.

# Arguments
- `dat::DataFrame`: DataFrame containing channel summary data.
- `col::Symbol`: The symbol representing the column in `dat` to plot on the y-axis.
- `kwargs...`: Keyword arguments for customization (see below for available options).

# Available Keyword Arguments
- `sort_values::Bool=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:sort_values])`: If true, sort the bars by the values in the `col` column in descending order.
- `average_over::Union{Symbol,Nothing}=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:average_over])`: Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI.
- `display_plot::Bool=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:display_plot])`: Whether to display the plot.
- `bar_color::Symbol=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:bar_color])`: Color of the bars.
- `error_color::Symbol=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:error_color])`: Color of error bars.
- `error_linewidth::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:error_linewidth])`: Line width of error bars.
- `xtick_rotation::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:xtick_rotation])`: Rotation angle for x-axis tick labels.
- `xlabel::String="$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:xlabel])"`: Label for x-axis.
- `title::String="$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:title])"`: Plot title.
- `title_fontsize::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:title_fontsize])`: Font size for title.
- `label_fontsize::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:label_fontsize])`: Font size for axis labels.
- `tick_fontsize::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:tick_fontsize])`: Font size for tick labels.
- `bar_width::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:bar_width])`: Width of bars.
- `bar_alpha::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:bar_alpha])`: Transparency of bars.
- `grid_visible::Bool=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:grid_visible])`: Whether to show grid.
- `grid_alpha::Real=$(DEFAULT_CHANNEL_SUMMARY_KWARGS[:grid_alpha])`: Transparency of grid.

# Returns
- `Figure`: The Makie Figure object.
- `Axis`: The Makie Axis object for the bar plot.

# Example
    # Assuming summary_df has columns :channel, :kurtosis, :variance, :epoch
    fig_kurt, ax_kurt = plot_channel_summary(summary_df, :kurtosis, average_over=:epoch, title="Kurtosis by Channel")
    # display(fig_kurt)

    fig_var_sorted, ax_var_sorted = plot_channel_summary(summary_df, :variance, sort_values=true, average_over=:epoch, xlabel="Electrode", bar_color=:red)
    # display(fig_var_sorted)
"""
function plot_channel_summary(
    dat::DataFrame,
    col::Symbol;
    kwargs...
)
    # Merge user kwargs with defaults
    plot_kwargs = merge(DEFAULT_CHANNEL_SUMMARY_KWARGS, Dict(kwargs))
    
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, dat, col; kwargs...)

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    return fig, ax
end
