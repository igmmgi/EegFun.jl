"""
    plot_correlation_heatmap!(fig::Figure, ax::Axis, corr_df::DataFrame, mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing)

Plot a heatmap of a correlation matrix stored in a DataFrame on the provided figure and axis.

This is the mutating version that plots directly on the provided `fig` and `ax` objects.

# Arguments
- `fig::Figure`: The Makie Figure object to plot on
- `ax::Axis`: The Makie Axis object to plot on
- `corr_df::DataFrame`: DataFrame containing the correlation matrix. First column for row labels, rest for correlation values.
- `mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing`: Optional tuple `(min_val, max_val)`. Correlations within this range (inclusive) will be masked (set to NaN) in the heatmap, making them transparent. Defaults to `nothing` (no masking).

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
    
    # Add a colorbar
    Colorbar(fig[1, 2], limits = (-1, 1), label = "Correlation")
"""
function plot_correlation_heatmap!(
    fig::Figure,
    ax::Axis,
    corr_df::DataFrame,
    mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing,
)
    # Extract the correlation matrix (excluding the row names column)
    # Convert to Float64 to handle potential NaNs from masking
    corr_matrix = Matrix{Float64}(corr_df[:, 2:end])

    # Mask values within the specified range
    if !isnothing(mask_range)
        min_val, max_val = mask_range
        corr_matrix[(corr_matrix .>= min_val) .& (corr_matrix .<= max_val)] .= NaN
    end

    # Extract row and column names
    row_names = String.(corr_df[!, 1]) # Use first column by index
    col_names = String.(propertynames(corr_df)[2:end])

    # Configure the axis
    ax.xlabel = ""
    ax.ylabel = ""
    ax.xticks = (1:length(col_names), col_names)
    ax.yticks = (1:length(row_names), row_names)
    ax.xticklabelrotation = pi/4

    # Create the heatmap
    heatmap!(ax, corr_matrix, colormap = :viridis, colorrange = (-1, 1), nan_color = :transparent)

    return nothing
end

"""
    plot_correlation_heatmap(corr_df::DataFrame, mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing)

Plot a heatmap of a correlation matrix stored in a DataFrame.

Assumes the first column of `corr_df` contains the row labels (e.g., channel names
as Symbols or Strings) and subsequent columns represent the correlations with other
variables (column names).

# Arguments
- `corr_df::DataFrame`: DataFrame containing the correlation matrix. First column for row labels, rest for correlation values.
- `mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing`: Optional tuple `(min_val, max_val)`. Correlations within this range (inclusive) will be masked (set to NaN) in the heatmap, making them transparent. Defaults to `nothing` (no masking).

# Returns
- `Figure`: The Makie Figure object containing the heatmap and colorbar.
- `Axis`: The Makie Axis object for the heatmap itself.

# Example
    # Generate sample correlation data
    labels = ["Ch1", "Ch2", "Ch3", "Ch4"]
    matrix = cor(rand(100, 4))
    corr_df = DataFrame(matrix, labels)
    insertcols!(corr_df, 1, :row => labels) # Add row label column

    # Plot the heatmap
    fig, ax = plot_correlation_heatmap(corr_df)
    # display(fig) # Display if needed

    # Plot with masking for weak correlations
    fig_masked, ax_masked = plot_correlation_heatmap(corr_df, (-0.3, 0.3))
    # display(fig_masked)
"""
function plot_correlation_heatmap(
    corr_df::DataFrame,
    mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    display_plot::Bool = true,
)
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Use the mutating version to plot
    plot_correlation_heatmap!(fig, ax, corr_df, mask_range)

    # Add a colorbar
    Colorbar(fig[1, 2], limits = (-1, 1), label = "Correlation")

    if display_plot
        display_figure(fig)
    end

    return fig, ax

end
