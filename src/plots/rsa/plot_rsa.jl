"""
Plotting functions for Representational Similarity Analysis (RSA) results.

This module provides functions for visualizing RDMs, model correlations,
and other RSA analysis results.
"""
# ==============================================================================
#   DEFAULT KEYWORD ARGUMENTS
# ==============================================================================

const PLOT_RSA_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("RSA Results", "Title for the plot window"),
    :interactive => (true, "Enable interactive features (true/false)"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("Correlation", "Label for y-axis"),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Show title (true/false)"),

    # Line styling
    :linewidth => (2, "Line width"),
    :color => (:blue, "Color"),
    :linestyle => (:solid, "Line style"),

    # Grid
    :xgrid => (true, "Show x-axis grid (true/false)"),
    :ygrid => (true, "Show y-axis grid (true/false)"),

    # RDM heatmap
    :colormap => (:viridis, "Colormap for RDM heatmap"),
    :show_colorbar => (true, "Show colorbar for RDM (true/false)"),
)

"""
    plot_rdm(rsa_data::RsaData; time_point::Union{Float64, Int, Nothing} = nothing, kwargs...)

Plot Representational Dissimilarity Matrix (RDM) as a heatmap.

# Arguments
- `rsa_data::RsaData`: RSA results
- `time_point::Union{Float64, Int, Nothing}`: Time point to plot (Float64 = time in seconds, Int = index, nothing = average across time)
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot RDM at a specific time point
plot_rdm(rsa_result, time_point=0.3)

# Plot average RDM across all time points
plot_rdm(rsa_result)

# Plot RDM at time index 50
plot_rdm(rsa_result, time_point=50)
```
"""
function plot_rdm(rsa_data::RsaData; time_point::Union{Float64,Int,Nothing} = nothing, kwargs...)
    # Merge defaults with user kwargs
    plot_kwargs = merge(PLOT_RSA_KWARGS, Dict(kwargs))

    # Extract parameters (unwrap tuples if needed)
    _get_val(k) = isa(plot_kwargs[k], Tuple) ? plot_kwargs[k][1] : plot_kwargs[k]
    display_plot = _get_val(:display_plot)
    figure_title = _get_val(:figure_title)
    title_text = _get_val(:title)
    show_title = _get_val(:show_title)

    # Determine which RDM to plot
    n_conditions = length(rsa_data.condition_names)
    n_times = length(rsa_data.times)

    if isnothing(time_point)
        # Average across time
        rdm_to_plot = mean(rsa_data.rdm, dims = 1)[1, :, :]
        time_label = "Average"
    elseif isa(time_point, Int)
        # Use time index
        if time_point < 1 || time_point > n_times
            @minimal_error_throw("Time point index $time_point out of range [1, $n_times]")
        end
        rdm_to_plot = rsa_data.rdm[time_point, :, :]
        time_label = "$(round(rsa_data.times[time_point], digits=3)) s"
    else
        # Find closest time point
        time_idx = argmin(abs.(rsa_data.times .- time_point))
        rdm_to_plot = rsa_data.rdm[time_idx, :, :]
        time_label = "$(round(rsa_data.times[time_idx], digits=3)) s"
    end

    # Create figure
    fig = Figure(title = figure_title)

    if isempty(title_text)
        title_str = "RDM at $time_label"
    else
        title_str = title_text
    end

    ax = Axis(
        fig[1, 1],
        title = show_title ? title_str : "",
        xlabel = "Condition",
        ylabel = "Condition",
        xgridvisible = _get_val(:xgrid),
        ygridvisible = _get_val(:ygrid),
        aspect = DataAspect(),
    )

    # Set condition labels
    ax.xticks = (1:n_conditions, rsa_data.condition_names)
    ax.yticks = (1:n_conditions, rsa_data.condition_names)

    # Plot heatmap
    colormap_val = _get_val(:colormap)
    hm = heatmap!(ax, 1:n_conditions, 1:n_conditions, rdm_to_plot, colormap = colormap_val)

    # Add colorbar if requested
    if _get_val(:show_colorbar)
        Colorbar(fig[1, 2], hm, label = "Dissimilarity")
    end

    # Display if requested
    if display_plot
        display(fig)
    end

    return fig
end

"""
    plot_model_correlations(rsa_data::RsaData; kwargs...)

Plot model correlations over time.

# Arguments
- `rsa_data::RsaData`: RSA results with model comparisons
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot correlations with models
plot_model_correlations(rsa_result)

# Custom styling
plot_model_correlations(rsa_result, colors=[:red, :blue], title="Model Comparison")
```
"""
function plot_model_correlations(rsa_data::RsaData; kwargs...)
    if isnothing(rsa_data.model_correlations) || isnothing(rsa_data.model_names)
        @minimal_error_throw("RSA data does not contain model correlations. Run compare_models() first.")
    end

    # Merge defaults with user kwargs
    plot_kwargs = merge(PLOT_RSA_KWARGS, Dict(kwargs))

    # Extract parameters (unwrap tuples if needed)
    _get_val(k) = isa(plot_kwargs[k], Tuple) ? plot_kwargs[k][1] : plot_kwargs[k]
    display_plot = _get_val(:display_plot)
    figure_title = _get_val(:figure_title)
    title_text = _get_val(:title)
    show_title = _get_val(:show_title)

    times = rsa_data.times
    correlations = rsa_data.model_correlations
    model_names = rsa_data.model_names
    n_models = length(model_names)

    # Auto-generate colors if not provided
    colors = get(kwargs, :colors, nothing)
    if isnothing(colors)
        colors = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]
        colors = [colors[mod1(i, length(colors))] for i = 1:n_models]
    elseif colors isa Symbol
        colors = fill(colors, n_models)
    end

    # Create figure
    fig = Figure(title = figure_title)

    if isempty(title_text)
        title_str = "Model Correlations"
    else
        title_str = title_text
    end

    ax = Axis(
        fig[1, 1],
        title = show_title ? title_str : "",
        xlabel = _get_val(:xlabel),
        ylabel = _get_val(:ylabel),
        xgridvisible = _get_val(:xgrid),
        ygridvisible = _get_val(:ygrid),
    )

    # Set axis limits
    xlim_val = _get_val(:xlim)
    if !isnothing(xlim_val)
        xlims!(ax, xlim_val)
    else
        xlims!(ax, (times[1], times[end]))
    end

    ylim_val = _get_val(:ylim)
    if !isnothing(ylim_val)
        ylims!(ax, ylim_val)
    else
        correlations_clean = filter(!isnan, vec(correlations))
        if isempty(correlations_clean)
            ylims!(ax, (-1.0, 1.0))  # Default range if all NaN
        else
            y_min = minimum(correlations_clean)
            y_max = maximum(correlations_clean)
            y_range = y_max - y_min
            if y_range == 0
                ylims!(ax, (y_min - 0.1, y_max + 0.1))
            else
                ylims!(ax, (y_min - 0.1 * y_range, y_max + 0.1 * y_range))
            end
        end
    end

    # Plot correlations for each model
    for (model_idx, model_name) in enumerate(model_names)
        lines!(
            ax,
            times,
            correlations[:, model_idx],
            color = colors[model_idx],
            linewidth = _get_val(:linewidth),
            linestyle = _get_val(:linestyle),
            label = model_name,
        )

        # Add significance markers if p-values available
        if !isnothing(rsa_data.p_values)
            sig_threshold = 0.05
            sig_indices = findall(rsa_data.p_values[:, model_idx] .< sig_threshold)
            if !isempty(sig_indices)
                scatter!(
                    ax,
                    times[sig_indices],
                    correlations[sig_indices, model_idx],
                    color = colors[model_idx],
                    marker = :star,
                    markersize = 10,
                )
            end
        end
    end

    # Add zero line
    hlines!(ax, 0, color = :black, linewidth = 1, linestyle = :dash)

    # Add legend
    axislegend(ax, position = :rt)

    # Display if requested
    if display_plot
        display(fig)
    end

    return fig
end

"""
    plot_rsa(rsa_data::RsaData; kwargs...)

Main plotting function for RSA results.

If model correlations are available, plots them. Otherwise, plots RDM at average time.

# Arguments
- `rsa_data::RsaData`: RSA results
- `kwargs`: Additional keyword arguments (passed to plot_model_correlations or plot_rdm)

# Examples
```julia
# Plot RSA results (automatically chooses best visualization)
plot_rsa(rsa_result)

# Force RDM plot
plot_rsa(rsa_result, plot_type=:rdm, time_point=0.3)
```
"""
function plot_rsa(rsa_data::RsaData; kwargs...)
    plot_type = get(kwargs, :plot_type, nothing)

    if isnothing(plot_type)
        # Auto-detect: prefer model correlations if available
        if !isnothing(rsa_data.model_correlations)
            return plot_model_correlations(rsa_data; kwargs...)
        else
            return plot_rdm(rsa_data; kwargs...)
        end
    elseif plot_type == :rdm
        return plot_rdm(rsa_data; kwargs...)
    elseif plot_type == :correlations
        return plot_model_correlations(rsa_data; kwargs...)
    else
        @minimal_error_throw("Unknown plot_type: $plot_type. Use :rdm or :correlations")
    end
end

