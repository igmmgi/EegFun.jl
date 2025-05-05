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
function plot_correlation_heatmap(corr_df::DataFrame, mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing)
    # Extract the correlation matrix (excluding the row names column)
    # Convert to Float64 to handle potential NaNs from masking
    corr_matrix = Matrix{Float64}(corr_df[:, 2:end])

    # Mask values within the specified range
    if !isnothing(mask_range)
        min_val, max_val = mask_range
        corr_matrix[(corr_matrix.>=min_val).&(corr_matrix.<=max_val)] .= NaN
    end

    # Extract row and column names
    row_names = String.(corr_df[!, 1]) # Use first column by index
    col_names = String.(propertynames(corr_df)[2:end])

    # Create the heatmap
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "",
        ylabel = "",
        xticks = (1:length(col_names), col_names),
        yticks = (1:length(row_names), row_names),
        xticklabelrotation = pi/4
    )
    heatmap!(ax, corr_matrix, colormap = :viridis, colorrange = (-1, 1), nan_color = :transparent)

    # Add a colorbar
    Colorbar(fig[1, 2], limits = (-1, 1), label = "Correlation")

    return fig, ax
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
function plot_joint_probability(dat::DataFrame)
    # Basic validation
    if :channel ∉ propertynames(dat) || :jp ∉ propertynames(dat)
        error("DataFrame must contain :channel and :jp columns.")
    end
    channel_names = String.(dat.channel)

    fig = Figure()
    ax = Axis(fig[1, 1], xticks = (1:length(channel_names), channel_names), xticklabelrotation = pi / 4)
    barplot!(ax, 1:nrow(dat), dat[!, :jp])
    ax.xlabel = "Electrode"
    ax.ylabel = "Joint Probability"

    return fig, ax

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
function plot_channel_summary(dat::DataFrame, col::Symbol; sort_values=false)
    # Check if required columns exist
    if :channel ∉ propertynames(dat) || col ∉ propertynames(dat)
        error("DataFrame must contain :channel and :$col columns.")
    end

    # Optionally sort data
    plot_dat = if sort_values
        sort(dat, col, rev=true)
    else
        dat
    end

    channel_names = String.(plot_dat.channel)
    values_to_plot = plot_dat[!, col]

    fig = Figure()
    ax = Axis(fig[1, 1], xticks = (1:length(channel_names), channel_names), xticklabelrotation = pi / 4)
    barplot!(ax, 1:nrow(plot_dat), values_to_plot)
    ax.xlabel = "Electrode"
    ax.ylabel = String(col) # Use column name directly

    return fig, ax

end
