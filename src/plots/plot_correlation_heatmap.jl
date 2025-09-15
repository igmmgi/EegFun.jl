# Default parameters for correlation heatmap plots with descriptions
const PLOT_CORRELATION_HEATMAP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Data parameters
    :mask_range => (nothing, "Optional tuple (min_val, max_val) to mask correlations within this range"),
    
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    
    # Heatmap styling
    :colormap => (:viridis, "Colormap for the heatmap"),
    :colorrange => ((-1, 1), "Color range for the heatmap"),
    :nan_color => (:transparent, "Color for NaN values"),
    
    # Axis styling
    :xlabel => ("", "Label for x-axis"),
    :ylabel => ("", "Label for y-axis"),
    :xtick_rotation => (π/4, "Rotation angle for x-axis tick labels"),
    :ytick_rotation => (0, "Rotation angle for y-axis tick labels"),
    
    # Colorbar parameters
    :plot_colorbar => (true, "Whether to display the colorbar"),
    :colorbar_width => (30, "Width of the colorbar"),
    :colorbar_label => ("Correlation", "Label for the colorbar"),
    :colorbar_limits => ((-1, 1), "Limits for the colorbar"),
)

"""
    plot_correlation_heatmap!(fig::Figure, ax::Axis, corr_df::DataFrame; kwargs...)

Plot a heatmap of a correlation matrix stored in a DataFrame on the provided figure and axis.

This is the mutating version that plots directly on the provided `fig` and `ax` objects.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `ax::Axis`: The Makie Axis object to plot on
- `corr_df::DataFrame`: DataFrame containing the correlation matrix. First column for row labels, rest for correlation values.

$(generate_kwargs_doc(PLOT_CORRELATION_HEATMAP_KWARGS))

# Returns
- `nothing` (modifies the provided figure and axis in-place)

# Example
    # Generate sample correlation data
    labels = ["Ch1", "Ch2", "Ch3", "Ch4"]
    matrix = cor(rand(100, 4))
    corr_df = DataFrame(matrix, labels)
    insertcols!(corr_df, 1, :row => labels) # Add row label column

    # Create figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    
    # Plot the heatmap on the existing figure/axis
    plot_correlation_heatmap!(fig, ax, corr_df)
    
    # Plot with custom styling
    plot_correlation_heatmap!(fig, ax, corr_df; 
                             colormap = :plasma, 
                             mask_range = (-0.3, 0.3),
                             xtick_rotation = π/6)
"""
function plot_correlation_heatmap!(
    fig::Figure,
    ax::Axis,
    corr_df::DataFrame;
    kwargs...
)
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

    # Create the heatmap
    heatmap!(ax, corr_matrix, 
             colormap = plot_kwargs[:colormap], 
             colorrange = plot_kwargs[:colorrange], 
             nan_color = plot_kwargs[:nan_color])

    return nothing
end

"""
    plot_correlation_heatmap(corr_df::DataFrame; kwargs...)

Plot a heatmap of a correlation matrix stored in a DataFrame.

Assumes the first column of `corr_df` contains the row labels (e.g., channel names
as Symbols or Strings) and subsequent columns represent the correlations with other
variables (column names).

# Arguments
- `corr_df::DataFrame`: DataFrame containing the correlation matrix. First column for row labels, rest for correlation values.

$(generate_kwargs_doc(PLOT_CORRELATION_HEATMAP_KWARGS))

# Returns
- `Figure`: The Makie Figure object containing the heatmap and colorbar.
- `Axis`: The Makie Axis object for the heatmap itself.

# Example
    # Generate sample correlation data
    labels = ["Ch1", "Ch2", "Ch3", "Ch4"]
    matrix = cor(rand(100, 4))
    corr_df = DataFrame(matrix, labels)
    insertcols!(corr_df, 1, :row => labels) # Add row label column

    # Plot the heatmap with default settings
    fig, ax = plot_correlation_heatmap(corr_df)

    # Plot with custom styling and masking
    fig_masked, ax_masked = plot_correlation_heatmap(corr_df; 
                                                   mask_range = (-0.3, 0.3),
                                                   colormap = :plasma,
                                                   colorbar_label = "Custom Label")
"""
function plot_correlation_heatmap(
    corr_df::DataFrame;
    kwargs...
)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_CORRELATION_HEATMAP_KWARGS, kwargs)
    
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Use the mutating version to plot
    plot_correlation_heatmap!(fig, ax, corr_df; kwargs...)

    # Add a colorbar if requested
    if plot_kwargs[:plot_colorbar]
        Colorbar(fig[1, 2], 
                limits = plot_kwargs[:colorbar_limits], 
                label = plot_kwargs[:colorbar_label],
                width = plot_kwargs[:colorbar_width])
    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    return fig, ax
end
