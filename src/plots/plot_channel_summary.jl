# Default parameters for channel summary plots with descriptions
const DEFAULT_CHANNEL_SUMMARY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
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

# Helper function to extract just the defaults
function _get_defaults(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::Dict{Symbol,Any}
    return Dict(key => value[1] for (key, value) in kwargs_dict)
end

# Helper function to generate documentation
function generate_kwargs_doc(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::String
    doc_lines = ["# Keyword Arguments"]
    push!(doc_lines, "All keyword arguments below have sensible defaults defined in `DEFAULT_CHANNEL_SUMMARY_KWARGS`.")
    push!(doc_lines, "You can override any of these defaults by passing the corresponding keyword argument.")
    push!(doc_lines, "")
    for (param_name, (default_val, desc)) in kwargs_dict
        type_info = typeof(default_val)
        push!(doc_lines, "- `$(param_name)::$(type_info)=$(default_val)`: $(desc)")
    end
    return join(doc_lines, "\n")
end

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

$(generate_kwargs_doc(DEFAULT_CHANNEL_SUMMARY_KWARGS))

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Examples
```julia
# Basic usage with defaults
fig = Figure()
ax = Axis(fig[1, 1])
plot_channel_summary!(fig, ax, summary_df, :kurtosis)

# Customize appearance
plot_channel_summary!(fig, ax, summary_df, :kurtosis, 
    bar_color = :red, 
    title = "Custom Title",
    sort_values = true)

# With averaging and error bars
plot_channel_summary!(fig, ax, summary_df, :kurtosis,
    average_over = :epoch,
    error_color = :blue,
    error_linewidth = 3)
```

# See Also
- `plot_channel_summary` for the non-mutating version
"""
function plot_channel_summary!(
    fig::Figure,
    ax::Axis,
    dat::DataFrame,
    col::Symbol;
    kwargs...
)
    # Merge user kwargs with defaults
    plot_kwargs = merge(_get_defaults(DEFAULT_CHANNEL_SUMMARY_KWARGS), Dict(kwargs))

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

$(generate_kwargs_doc(DEFAULT_CHANNEL_SUMMARY_KWARGS))

# Returns
- `Figure`: The Makie Figure object.
- `Axis`: The Makie Axis object for the bar plot.

# Examples
```julia
# Basic usage with defaults
fig, ax = plot_channel_summary(summary_df, :kurtosis)

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

# See Also
- `DEFAULT_CHANNEL_SUMMARY_KWARGS` for the complete list of default values
- `plot_channel_summary!` for the mutating version
"""
function plot_channel_summary(
    dat::DataFrame,
    col::Symbol;
    kwargs...
)
    # Merge user kwargs with defaults
    plot_kwargs = merge(_get_defaults(DEFAULT_CHANNEL_SUMMARY_KWARGS), Dict(kwargs))
    
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, dat, col; kwargs...)

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    return fig, ax
end
