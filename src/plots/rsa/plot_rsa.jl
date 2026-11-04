# Plotting functions for Representational Similarity Analysis (RSA) results.

const PLOT_RSA_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("", "Title drawn at the top of the entire figure canvas"),
    :figure_title_fontsize => (24, "Font size for figure title"),
    :interactive => (true, "Enable interactive features (true/false)"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("Correlation", "Label for y-axis"),

    # Title
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Show plot title (true/false)"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Line styling
    :linewidth => (2, "Line width"),
    :color => (nothing, "Color"),
    :colors => (nothing, "Custom colors for condition pairs or models"),
    :linestyle => (nothing, "Line style"),
    :legend => (true, "Show legend (true/false)"),

    # Grid
    :xgrid => (true, "Show x-axis grid (true/false)"),
    :ygrid => (true, "Show y-axis grid (true/false)"),

    # RDM heatmap
    :colormap => (nothing, "Colormap for RDM heatmap"),
    :colorrange => (nothing, "Color scale limits as (min, max) tuple for RDM heatmap"),
    :colorbar_plot => (true, "Show colorbar for RDM (true/false)"),
    :plot_type => (nothing, "Plot type for plot_rsa (:rdm or :correlations)"),
)

const _RSA_PAIR_PALETTE = (:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :cyan, :magenta)
const _RSA_MODEL_PALETTE = (:blue, :red, :green, :orange, :purple, :brown, :pink, :gray)

function _rsa_title(plot_kwargs::Dict, default_title::String)
    plot_kwargs[:show_title] || return ""
    title_text = plot_kwargs[:plot_title]
    return isnothing(title_text) ? default_title : title_text
end

@inline function _rsa_series_color(colors, palette::Tuple, idx::Int)
    isnothing(colors) && return palette[mod1(idx, length(palette))]
    colors isa Symbol && return colors
    return colors[mod1(idx, length(colors))]
end

function _finalize_rsa_figure!(fig::Figure, plot_kwargs::Dict)
    figure_title = plot_kwargs[:figure_title]
    if !isempty(figure_title)
        Label(fig[0, :], figure_title, fontsize = plot_kwargs[:figure_title_fontsize], font = :bold, tellwidth = false)
    end
    plot_kwargs[:display_plot] && _display_figure(fig)
    return nothing
end

"""
    plot_rdm_heatmap(rsa_data::RsaData; time_point::Union{Float64, Int, Nothing} = nothing, kwargs...)

Plot Representational Dissimilarity Matrix (RDM) as a heatmap.

Shows a matrix visualization of pairwise dissimilarities between conditions
at a specific time point or averaged across time.

# Arguments
- `rsa_data::RsaData`: RSA results
- `time_point::Union{Float64, Int, Nothing}`: Time point to plot (Float64 = time in seconds, Int = index, nothing = average across time)
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot RDM at a specific time point
plot_rdm_heatmap(rsa_result, time_point=0.3)

# Plot average RDM across all time points
plot_rdm_heatmap(rsa_result)

# Plot RDM at time index 50
plot_rdm_heatmap(rsa_result, time_point=50)
```
"""
function plot_rdm_heatmap(rsa_data::RsaData; time_point::Union{Float64,Int,Nothing} = nothing, kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_RSA_KWARGS, kwargs)
    n_conditions = length(rsa_data.condition_names)
    n_times = length(rsa_data.times)

    if isnothing(time_point)
        rdm_to_plot = dropdims(mean(rsa_data.rdm, dims = 1), dims = 1)
        time_label = "Average"
    elseif time_point isa Int
        (time_point < 1 || time_point > n_times) && @minimal_error("Time point index $time_point out of range [1, $n_times]")
        rdm_to_plot = @view rsa_data.rdm[time_point, :, :]
        time_label = "$(round(rsa_data.times[time_point], digits=3)) s"
    else
        time_idx = find_closest_time_index(rsa_data.times, time_point)
        rdm_to_plot = @view rsa_data.rdm[time_idx, :, :]
        time_label = "$(round(rsa_data.times[time_idx], digits=3)) s"
    end

    fig = _eegfun_figure()
    ax = Axis(
        fig[1, 1],
        title = _rsa_title(plot_kwargs, "RDM at $time_label"),
        xlabel = "Condition",
        ylabel = "Condition",
        xgridvisible = plot_kwargs[:xgrid],
        ygridvisible = plot_kwargs[:ygrid],
        aspect = DataAspect(),
    )

    ax.xticks = (1:n_conditions, rsa_data.condition_names)
    ax.yticks = (1:n_conditions, rsa_data.condition_names)

    colormap_val = _resolve_theme_colormap(ax, plot_kwargs[:colormap])
    cr_val = plot_kwargs[:colorrange]
    rdm_perm = permutedims(rdm_to_plot)
    hm =
        isnothing(cr_val) ? heatmap!(ax, 1:n_conditions, 1:n_conditions, rdm_perm, colormap = colormap_val) :
        heatmap!(ax, 1:n_conditions, 1:n_conditions, rdm_perm, colormap = colormap_val, colorrange = cr_val)

    plot_kwargs[:colorbar_plot] && Colorbar(fig[1, 2], hm, label = "Dissimilarity")
    _finalize_rsa_figure!(fig, plot_kwargs)
    return (fig = fig, axes = [ax])
end

"""
    plot_rdm_timecourse(
        rsa_data::RsaData;
        condition_pairs::Union{Vector{Tuple{Int,Int}}, Symbol, Nothing} = :all,
        kwargs...
    )

Plot dissimilarity timecourse for condition pairs.

Shows how dissimilarity between conditions evolves over time as lines.
This helps visualize when conditions become distinguishable.

# Arguments
- `rsa_data::RsaData`: RSA results
- `condition_pairs`: Which condition pairs to plot
  - `:all` - All pairwise comparisons (default)
  - `:upper` - Upper triangle only (same as :all, no duplicates)
  - `Vector{Tuple{Int,Int}}` - Specific pairs, e.g., `[(1,2), (1,3)]`
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot all condition pairs
plot_rdm_timecourse(rsa_result)

# Plot only specific pairs
plot_rdm_timecourse(rsa_result, condition_pairs=[(1,2), (1,3)])

# Custom styling
plot_rdm_timecourse(rsa_result, plot_title="Dissimilarity Over Time")
```
"""
function plot_rdm_timecourse(rsa_data::RsaData; condition_pairs::Union{Vector{Tuple{Int,Int}},Symbol,Nothing} = :all, kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_RSA_KWARGS, kwargs)
    n_conditions = length(rsa_data.condition_names)
    times = rsa_data.times

    if isnothing(condition_pairs) || condition_pairs == :all || condition_pairs == :upper
        pairs = [(i, j) for i = 1:n_conditions for j = (i+1):n_conditions]
    else
        pairs = condition_pairs
        for (i, j) in pairs
            (i < 1 || i > n_conditions || j < 1 || j > n_conditions) &&
                @minimal_error("Invalid condition pair ($i, $j). Conditions range from 1 to $n_conditions")
        end
    end

    fig = _eegfun_figure()
    ax = Axis(
        fig[1, 1],
        title = _rsa_title(plot_kwargs, "RDM Timecourse"),
        xlabel = plot_kwargs[:xlabel],
        ylabel = "Dissimilarity",
        xgridvisible = plot_kwargs[:xgrid],
        ygridvisible = plot_kwargs[:ygrid],
    )

    xlim_val = plot_kwargs[:xlim]
    xlims!(ax, isnothing(xlim_val) ? (times[1], times[end]) : xlim_val)
    !isnothing(plot_kwargs[:ylim]) && ylims!(ax, plot_kwargs[:ylim])

    colors = plot_kwargs[:colors]
    line_style = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]
    lw = plot_kwargs[:linewidth]

    for (pair_idx, (i, j)) in enumerate(pairs)
        lines!(
            ax,
            times,
            @view(rsa_data.rdm[:, i, j]),
            color = _rsa_series_color(colors, _RSA_PAIR_PALETTE, pair_idx),
            linewidth = lw,
            linestyle = line_style,
            label = "$(rsa_data.condition_names[i]) vs $(rsa_data.condition_names[j])",
        )
    end

    hlines!(ax, 0, color = :black, linewidth = 1, linestyle = :dash)
    (plot_kwargs[:legend] && !isempty(pairs)) && axislegend(ax, position = :rt)
    _finalize_rsa_figure!(fig, plot_kwargs)
    return (fig = fig, axes = [ax])
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
plot_model_correlations(rsa_result, colors=[:red, :blue], plot_title="Model Comparison")
```
"""
function plot_model_correlations(rsa_data::RsaData; kwargs...)
    (isnothing(rsa_data.model_correlations) || isnothing(rsa_data.model_names)) &&
        @minimal_error("RSA data does not contain model correlations. Run compare_models() first.")

    plot_kwargs = _merge_plot_kwargs(PLOT_RSA_KWARGS, kwargs)
    times = rsa_data.times
    correlations = rsa_data.model_correlations
    model_names = rsa_data.model_names
    n_models = length(model_names)

    fig = _eegfun_figure()
    ax = Axis(
        fig[1, 1],
        title = _rsa_title(plot_kwargs, "Model Correlations"),
        xlabel = plot_kwargs[:xlabel],
        ylabel = plot_kwargs[:ylabel],
        xgridvisible = plot_kwargs[:xgrid],
        ygridvisible = plot_kwargs[:ygrid],
    )

    xlim_val = plot_kwargs[:xlim]
    xlims!(ax, isnothing(xlim_val) ? (times[1], times[end]) : xlim_val)

    ylim_val = plot_kwargs[:ylim]
    if !isnothing(ylim_val)
        ylims!(ax, ylim_val)
    else
        y_min, y_max = Inf, -Inf
        @inbounds for v in correlations
            if !isnan(v)
                y_min = min(y_min, v)
                y_max = max(y_max, v)
            end
        end
        if isinf(y_min)
            ylims!(ax, (-1.0, 1.0))
        else
            y_range = y_max - y_min
            pad = y_range == 0 ? 0.1 : 0.1 * y_range
            ylims!(ax, (y_min - pad, y_max + pad))
        end
    end

    colors = plot_kwargs[:colors]
    line_style = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]
    lw = plot_kwargs[:linewidth]
    p_values = rsa_data.p_values

    for (model_idx, model_name) in enumerate(model_names)
        model_color = _rsa_series_color(colors, _RSA_MODEL_PALETTE, model_idx)
        corr_col = @view correlations[:, model_idx]
        lines!(ax, times, corr_col, color = model_color, linewidth = lw, linestyle = line_style, label = model_name)

        if !isnothing(p_values)
            sig_indices = findall(<(0.05), @view(p_values[:, model_idx]))
            if !isempty(sig_indices)
                scatter!(ax, times[sig_indices], corr_col[sig_indices], color = model_color, marker = :star5, markersize = 10)
            end
        end
    end

    hlines!(ax, 0, color = :black, linewidth = 1, linestyle = :dash)
    (plot_kwargs[:legend] && n_models > 0) && axislegend(ax, position = :rt)
    _finalize_rsa_figure!(fig, plot_kwargs)
    return (fig = fig, axes = [ax])
end

"""
    plot_rsa(rsa_data::RsaData; kwargs...)

Main plotting function for RSA results.

If model correlations are available, plots them. Otherwise, plots RDM heatmap at average time.

# Arguments
- `rsa_data::RsaData`: RSA results
- `kwargs`: Additional keyword arguments (passed to plot_model_correlations or plot_rdm_heatmap)

# Examples
```julia
# Plot RSA results (automatically chooses best visualization)
plot_rsa(rsa_result)

# Force RDM heatmap plot
plot_rsa(rsa_result, plot_type=:rdm, time_point=0.3)
```

    plot_rsa(filepath::String; input_dir=pwd(), participant_selection=participants(), kwargs...)

Load RSA data and plot. Accepts either a `.jld2` filepath or a pattern
to discover and plot all matching files (one plot per file).

# Examples
```julia
plot_rsa("rsa_results.jld2")
plot_rsa("rsa_results")
```
"""
function plot_rsa(filepath::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_rsa(data; kwargs...)
    end

    files = _find_batch_files(filepath, input_dir, participant_selection)
    isempty(files) && @minimal_error "No files matching pattern '$filepath' in $input_dir"

    results = NamedTuple[]
    for file in sort(files, by = _natural_sort_key)
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(results, plot_rsa(data; kwargs...))
    end
    return results
end

function plot_rsa(rsa_data::RsaData; kwargs...)
    plot_type = get(kwargs, :plot_type, nothing)
    if isnothing(plot_type)
        return !isnothing(rsa_data.model_correlations) ? plot_model_correlations(rsa_data; kwargs...) :
               plot_rdm_heatmap(rsa_data; kwargs...)
    elseif plot_type == :rdm
        return plot_rdm_heatmap(rsa_data; kwargs...)
    elseif plot_type == :correlations
        return plot_model_correlations(rsa_data; kwargs...)
    else
        @minimal_error("Unknown plot_type: $plot_type. Use :rdm or :correlations")
    end
end
