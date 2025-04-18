function plot_correlation_heatmap(corr_df::DataFrame, mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing)
    """
    Plot a heatmap of the correlation matrix using Makie.

    Parameters:
    - corr_df: A DataFrame containing the correlation matrix with row and column names.
    """
    # Extract the correlation matrix (excluding the row names column)
    corr_matrix = Matrix(corr_df[:, 2:end])

    # Mask values within the specified range
    if !isnothing(mask_range)
        min_val, max_val = mask_range
        corr_matrix[(corr_matrix.>=min_val).&(corr_matrix.<=max_val)] .= NaN
    end

    # Extract row and column names
    row_names = String.(corr_df[!, :row])
    col_names = String.(propertynames(corr_df)[2:end])

    # Create the heatmap
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "",
        ylabel = "",
        xticks = (1:length(col_names), col_names),
        yticks = (1:length(row_names), row_names),
    )
    heatmap!(ax, corr_matrix, colormap = :viridis, colorrange = (-1, 1))

    # Add a colorbar
    Colorbar(fig[1, 2], limits = (-1, 1), label = "Correlation")

    # Display the figure
    display(fig)
    return fig, ax
end


function plot_joint_probability(dat::DataFrame)
    fig = Figure()
    ax = Axis(fig[1, 1], xticks = (1:length(dat.channel), String.(dat.channel)), xticklabelrotation = pi / 4)
    barplot!(ax, 1:nrow(dat), dat[!, :jp])
    ax.xlabel = "Electrode"
    ax.ylabel = "Joint Probability"
    display(fig)
    return fig, ax
end


function plot_channel_summary(dat::DataFrame, col::Symbol)
    fig = Figure()
    ax = Axis(fig[1, 1], xticks = (1:length(dat.channel), String.(dat.channel)), xticklabelrotation = pi / 4)
    barplot!(ax, 1:nrow(dat), dat[!, col])
    ax.xlabel = "Electrode"
    ax.ylabel = String(col)
    display(fig)
    return fig, ax
end


