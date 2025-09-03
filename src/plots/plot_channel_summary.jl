# Default parameters for channel summary plots with descriptions
# Dict is used for documentation and for defaults
const PLOT_CHANNEL_SUMMARY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :sort_values => (false, "If true, sort the bars by the values in the `col` column in descending order."),
    :average_over => (nothing, "Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI."),
    :display_plot => (true, "Whether to display the plot."),
    :bar_color => (:steelblue, "Color of the bars."),
    :bar_width => (0.8, "Width of bars."),
    :bar_alpha => (0.7, "Transparency of bars."),
    :error_color => (:black, "Color of error bars."),
    :error_linewidth => (2, "Line width of error bars."),
    :xlabel => ("Electrode", "Label for x-axis."),
    :title => ("", "Plot title."),
    :title_fontsize => (16, "Font size for title."),
    :label_fontsize => (14, "Font size for axis labels."),
    :tick_fontsize => (12, "Font size for tick labels."),
    :xtick_rotation => (π/4, "Rotation angle for x-axis tick labels."),
    :grid_visible => (true, "Whether to show grid."),
    :grid_alpha => (0.3, "Transparency of grid."),
)


"""
    plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; kwargs...)
    plot_channel_summary(dat::DataFrame, col::Symbol; kwargs...)

Plot a bar chart summarizing a specific metric per channel from a DataFrame.

Assumes the DataFrame `dat` contains at least two columns:
- `:channel`: Containing channel names or identifiers (will be used for x-axis labels).
- `col`: The column specified by the `col` argument, containing the values to plot.

$(SIGNATURES)

# Arguments
- `fig::Figure`: The Makie Figure object to plot on (mutating version only)
- `ax::Axis`: The Makie Axis object to plot on (mutating version only)
- `dat::DataFrame`: DataFrame containing channel summary data.
- `col::Symbol`: The symbol representing the column in `dat` to plot on the y-axis.

$(generate_kwargs_doc(PLOT_CHANNEL_SUMMARY_KWARGS))

# Returns
- **Mutating version**: `nothing` (modifies the provided figure and axis in-place)
- **Non-mutating version**: `(fig::Figure, ax::Axis)` - The created figure and axis objects

# Examples
```julia
# Non-mutating version (creates new figure)
fig, ax = plot_channel_summary(summary_df, :kurtosis)

# Mutating version (plots on existing figure)
fig = Figure()
ax = Axis(fig[1, 1])
plot_channel_summary!(fig, ax, summary_df, :kurtosis)

# Customize appearance
fig, ax = plot_channel_summary(summary_df, :kurtosis, 
    bar_color = :red, 
    title = "Custom Title",
    sort_values = true)

# With averaging and error bars
fig, ax = plot_channel_summary(summary_df, :kurtosis,
    average_over = :epoch,
    error_color = :blue,
    error_linewidth = 3)
```
"""
function plot_channel_summary!(
    fig::Figure,
    ax::Axis,
    dat::DataFrame,
    col::Symbol;
    kwargs...
)
    # Merge user kwargs with defaults and validate
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)

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

# Share documentation with the non-mutating version
@doc (@doc plot_channel_summary!) plot_channel_summary
function plot_channel_summary(
    dat::DataFrame,
    col::Symbol;
    kwargs...
)
    # Merge user kwargs with defaults and validate
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)
    
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, dat, col; kwargs...)

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    return fig, ax
end
