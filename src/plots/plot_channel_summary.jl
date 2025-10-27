# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_CHANNEL_SUMMARY_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Display plot."),
    
    :sort_values => (false, "Sort the bars by the values in the `col` column in descending order."),
    :average_over => (nothing, "Column to average over (e.g., :epoch). If specified, will compute mean ± 95% CI."),
    :dims => (nothing, "Tuple (rows, cols) for grid layout. default = best_rect(n_columns)."),
    :bar_color => (:steelblue, "Bar colour."),
    :bar_width => (0.8, "Bar width."),
    :bar_alpha => (0.7, "Bar alpha."),
    :error_color => (:black, "Error bar color."),
    :error_linewidth => (2, "Error bar linewidth."),
    :xlabel => ("Electrode", "x-axis label."),
    :ylabel => ("", "y-axis label"),
    :title => ("", "Plot title."),
    :title_fontsize => (16, "Font size for title."),
    :label_fontsize => (14, "Font size for axis labels."),
    :tick_fontsize => (12, "Font size for tick labels."),
    :xtick_rotation => (π/4, "Rotation angle for x-axis tick labels."),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

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
function plot_channel_summary!(fig::Figure, ax::Axis, dat::DataFrame, col::Symbol; kwargs...)
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
    _configure_axis!(ax, channel_names, col, plot_kwargs, n_epochs)

    # Create the bar plot
    barplot!(
        ax,
        1:length(values_to_plot),
        values_to_plot,
        color = plot_kwargs[:bar_color],
        width = plot_kwargs[:bar_width],
        alpha = plot_kwargs[:bar_alpha],
    )

    if plot_kwargs[:average_over] !== nothing # add error bars
        errorbars!(
            ax,
            1:length(values_to_plot),
            values_to_plot,
            margin_of_error,
            color = plot_kwargs[:error_color],
            linewidth = plot_kwargs[:error_linewidth],
        )
    end

    return nothing
end

# Share documentation with the non-mutating version
@doc (@doc plot_channel_summary!) plot_channel_summary
function plot_channel_summary(dat::DataFrame, col::Symbol; kwargs...)
    # Merge user kwargs with defaults and validate
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)

    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_channel_summary!(fig, ax, dat, col; plot_kwargs...)

    plot_kwargs[:display_plot] && display_figure(fig)
    return fig, ax
end


"""
    plot_channel_summary!(fig, ax, dat::DataFrame, col::Vector{Symbol}; kwargs...)
    plot_channel_summary(dat::DataFrame, col::Vector{Symbol}; kwargs...)

Plot multiple bar charts summarizing different metrics per channel from a DataFrame.

Creates a grid layout with one subplot per column specified in `col`. Each subplot shows
a bar chart for that specific metric across all channels.

# Arguments
- `fig`: The Makie Figure object to plot on (mutating version only)
- `ax`: The Makie Axis object to plot on (mutating version only) - **Note**: This parameter is ignored in the multiple column version as each subplot gets its own axis
- `dat::DataFrame`: DataFrame containing channel summary data
- `col::Vector{Symbol}`: Vector of column symbols to plot, each will get its own subplot

# Keyword Arguments
$(generate_kwargs_doc(PLOT_CHANNEL_SUMMARY_KWARGS))

# Returns
- **Mutating version**: `nothing` (modifies the provided figure in-place)
- **Non-mutating version**: `fig::Figure` - The created figure object

# Examples
```julia
# Plot multiple metrics in a grid
fig = plot_channel_summary(summary_df, [:kurtosis, :variance, :range])

# Custom grid layout
fig = plot_channel_summary(summary_df, [:min, :max, :std], dims = (2, 2))

# With averaging and error bars
fig = plot_channel_summary(summary_df, [:kurtosis, :variance], 
    average_over = :epoch,
    sort_values = true)
```
"""
function plot_channel_summary!(fig, ax, dat::DataFrame, col::Vector{Symbol}; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)
    _plot_multiple_columns!(fig, dat, col, plot_kwargs)
end

# Share documentation with the non-mutating version
@doc (@doc plot_channel_summary!) plot_channel_summary
function plot_channel_summary(dat::DataFrame, col::Vector{Symbol}; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_CHANNEL_SUMMARY_KWARGS, kwargs)
    fig = Figure()
    _plot_multiple_columns!(fig, dat, col, plot_kwargs)
    plot_kwargs[:display_plot] && display_figure(fig)
    return fig
end


"""
    _plot_multiple_columns!(fig::Figure, dat::DataFrame, col::Vector{Symbol}, plot_kwargs::Dict)

Internal helper function to plot multiple channel summary columns in a grid layout.

Creates a grid of subplots where each column in `col` gets its own subplot showing
a bar chart of that metric across all channels.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `dat::DataFrame`: DataFrame containing channel summary data
- `col::Vector{Symbol}`: Vector of column symbols to plot
- `plot_kwargs::Dict`: Dictionary containing plotting parameters

# Grid Layout
- If `plot_kwargs[:dims]` is `nothing`, uses `best_rect()` to determine optimal grid dimensions
- If `plot_kwargs[:dims]` is provided, uses the specified (rows, cols) tuple
- Plots are arranged in row-major order (left to right, top to bottom)
- Stops plotting when all columns are processed, even if grid has empty spaces
"""
function _plot_multiple_columns!(fig::Figure, dat::DataFrame, col::Vector{Symbol}, plot_kwargs::Dict)
    n_cols = length(col)

    if plot_kwargs[:dims] === nothing
        rs, cs = best_rect(n_cols)
    else
        dims = plot_kwargs[:dims]
        if !isa(dims, Tuple) || length(dims) != 2 || !all(isa.(dims, Int)) || any(dims .<= 0)
            @minimal_error("dims must be a tuple of two positive integers (rows, cols), got: $dims")
        end
        rs, cs = dims
    end

    count = 1
    for r = 1:rs
        for c = 1:cs
            ax = Axis(fig[r, c])
            plot_channel_summary!(fig, ax, dat, col[count]; plot_kwargs...)
            count += 1
            if count > n_cols
                break
            end
        end
    end
end


"""
    _configure_axis!(ax::Axis, channel_names::Vector{String}, col::Symbol, plot_kwargs::Dict, n_epochs::Union{Int,Nothing})

Internal helper function to configure axis properties for channel summary plots.

Sets up all axis properties including labels, ticks, titles, fonts, and grid styling
based on the provided plotting parameters.

# Arguments
- `ax::Axis`: The Makie Axis object to configure
- `channel_names::Vector{String}`: Names of channels for x-axis tick labels
- `col::Symbol`: Column name being plotted (used for y-axis label)
- `plot_kwargs::Dict`: Dictionary containing plotting parameters
- `n_epochs::Union{Int,Nothing}`: Number of epochs (used for error bar labels when averaging)

# Configuration Details
- **Y-axis label**: Uses custom label if provided, otherwise generates dynamic label based on column name and averaging settings
- **X-axis ticks**: Sets channel names as tick labels
- **Fonts**: Configures title, axis label, and tick label font sizes
- **Grid**: Sets grid visibility, width, and transparency
- **Rotation**: Applies x-axis tick label rotation
"""
function _configure_axis!(
    ax::Axis,
    channel_names::Vector{String},
    col::Symbol,
    plot_kwargs::Dict,
    n_epochs::Union{Int,Nothing},
)
    # Set ylabel - use custom if provided, otherwise use dynamic default
    if plot_kwargs[:ylabel] != ""
        ax.ylabel = plot_kwargs[:ylabel]
    else
        ax.ylabel = plot_kwargs[:average_over] !== nothing ? "$(String(col)) (± 95% CI n=$n_epochs)" : "$(String(col))"
    end

    # Configure ticks and labels
    ax.xticks = (1:length(channel_names), channel_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.xlabel = plot_kwargs[:xlabel]
    ax.title = plot_kwargs[:title]
    ax.titlesize = plot_kwargs[:title_fontsize]
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    # Configure grid using the new axis styling function
    _set_axis_grid!(ax; 
                     xgrid = plot_kwargs[:xgrid], 
                     ygrid = plot_kwargs[:ygrid],
                     xminorgrid = plot_kwargs[:xminorgrid], 
                     yminorgrid = plot_kwargs[:yminorgrid])
end
