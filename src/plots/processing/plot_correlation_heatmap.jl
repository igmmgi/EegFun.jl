# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_CORRELATION_HEATMAP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Data parameters
    :mask_range => (nothing, "Optional tuple (min_val, max_val) to mask correlations within this range"),

    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Title and labels
    :plot_title => (nothing, "Plot title"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),
    :plot_title_fontsize => (16, "Font size for the title"),

    # Heatmap styling
    :colormap => (:jet, "Colormap for the heatmap"),
    :colorrange => ((-1, 1), "Color range for the heatmap and colorbar (should be -1 to 1 for correlations)"),
    :nan_color => (:transparent, "Color for NaN values"),

    # Axis styling
    :xlabel => ("", "Label for x-axis"),
    :ylabel => ("", "Label for y-axis"),
    :label_fontsize => (14, "Font size for axis labels"),
    :xtick_rotation => (π / 4, "Rotation angle for x-axis tick labels"),
    :ytick_rotation => (0, "Rotation angle for y-axis tick labels"),
    :tick_fontsize => (12, "Font size for tick labels"),

    # Grid and layout
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Colorbar parameters - get all Colorbar attributes with their actual defaults
    [
        Symbol("colorbar_$(attr)") => (get(COLORBAR_DEFAULTS, attr, nothing), "Colorbar $(attr) parameter") for
        attr in propertynames(Colorbar)
    ]...,

    # Specific colorbar overrides for correlation heatmap
    :colorbar_plot => (true, "Whether to display the colorbar"),
    :colorbar_position => (:right, "Position of the colorbar (:right, :left, :top, :bottom, or tuple)"),
    :colorbar_label => ("Correlation", "Label for the colorbar"),
)

"""
    plot_correlation_heatmap!(fig::Figure, ax::Axis, corr_df::DataFrame; kwargs...)

Plot a heatmap of a correlation matrix stored in a DataFrame on the provided figure and axis.

This is the mutating version that plots directly on the provided `fig` and `ax` objects.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `ax::Axis`: The Makie Axis object to plot on
- `corr_df::DataFrame`: DataFrame containing the correlation matrix. First column for row labels, rest for correlation values.

$(_generate_kwargs_doc(PLOT_CORRELATION_HEATMAP_KWARGS))

# Returns
- **Mutating version**: `nothing` (modifies the provided figure and axis in-place)
- **Non-mutating version**: `(fig::Figure, ax::Axis)` - The created figure and axis objects

# Examples
```julia
# Generate sample correlation data
labels = ["Ch1", "Ch2", "Ch3", "Ch4"]
matrix = cor(rand(100, 4))
corr_df = DataFrame(matrix, labels)
insertcols!(corr_df, 1, :row => labels) # Add row label column

# Non-mutating version (creates new figure)
fig, ax = plot_correlation_heatmap(corr_df)

# Mutating version (plots on existing figure)
fig = Figure()
ax = Axis(fig[1, 1])
plot_correlation_heatmap!(fig, ax, corr_df)

# Customize appearance
fig, ax = plot_correlation_heatmap(corr_df; 
    mask_range = (-0.3, 0.3),
    colormap = :plasma,
    colorbar_label = "Custom Label")
```
"""
function plot_correlation_heatmap!(fig::Figure, ax::Axis, corr_df::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_CORRELATION_HEATMAP_KWARGS, kwargs)

    # Extract the correlation matrix (excluding the row names column)
    # Convert to Float64 to handle potential NaNs from masking
    corr_matrix = Matrix{Float64}(corr_df[:, 2:end])

    # Mask values within the specified range
    if !isnothing(plot_kwargs[:mask_range])
        min_val, max_val = plot_kwargs[:mask_range]
        corr_matrix[(corr_matrix .>= min_val) .& (corr_matrix .<= max_val)] .= NaN
    end

    # Use the specified colorrange
    colorrange = plot_kwargs[:colorrange]
    # Validate colorrange for correlations
    min_val, max_val = colorrange
    if min_val < -1 || max_val > 1
        @minimal_warning("Color range ($min_val, $max_val) extends beyond typical correlation bounds (-1, 1)")
    end

    # Extract row and column names
    row_names = String.(corr_df[!, 1]) # Use first column by index
    col_names = String.(propertynames(corr_df)[2:end])

    # Configure the axis
    ax.xlabel = plot_kwargs[:xlabel]
    ax.ylabel = plot_kwargs[:ylabel]
    ax.xticks = (1:length(col_names), col_names)
    ax.yticks = (1:length(row_names), row_names)
    ax.xticklabelrotation = plot_kwargs[:xtick_rotation]
    ax.yticklabelrotation = plot_kwargs[:ytick_rotation]

    # Set font sizes
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    # Set title
    if !isnothing(plot_kwargs[:plot_title]) && !isempty(plot_kwargs[:plot_title])
        ax.title = plot_kwargs[:plot_title]
        ax.titlesize = plot_kwargs[:plot_title_fontsize]
    end

    # Configure grid using the new axis styling function
    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    # Create the heatmap
    # Makie's heatmap! displays matrices transposed (rows become columns), so we transpose
    # to ensure DataFrame rows appear as plot rows and DataFrame columns as plot columns
    heatmap!(ax, corr_matrix', colormap = plot_kwargs[:colormap], colorrange = colorrange, nan_color = plot_kwargs[:nan_color])

    # Add a colorbar if requested
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot, true)

    if colorbar_plot
        cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
        sg = ax.layoutobservables.gridcontent[].parent
        Colorbar(sg[cb_pos...]; colormap = plot_kwargs[:colormap], limits = colorrange, colorbar_kwargs...)
    end

    return nothing
end

"""
    plot_correlation_heatmap(corr_df::DataFrame; kwargs...)

Plot a heatmap of a correlation matrix. Creates a new figure.

See `plot_correlation_heatmap!` for full documentation of arguments and keyword arguments.

# Returns
- `(fig::Figure, ax::Axis)` - The created figure and axis objects
"""
function plot_correlation_heatmap(corr_df::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_CORRELATION_HEATMAP_KWARGS, kwargs)

    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Use the mutating version to plot
    plot_correlation_heatmap!(fig, ax, corr_df; kwargs...)

    if plot_kwargs[:display_plot]
        _set_window_title("Correlation Heatmap")
        _display_figure(fig)
    end

    return (fig = fig, axes = [ax])
end
