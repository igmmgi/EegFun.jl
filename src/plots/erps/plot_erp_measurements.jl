"""
    plot_erp_measurements(erp_data, analysis_type;
                          analysis_interval, baseline_interval,
                          layout, channel_selection, condition_selection, kwargs...)

Plot ERP data with measurement overlays computed inline.

Computes a single measurement type across all selected channels and conditions,
then overlays the results on the ERP waveforms. Combines measurement computation
and visualization in one call.

# Arguments
- `erp_data::Union{String, ErpData, Vector{ErpData}}`: ERP data (filepath or data object)
- `analysis_type::String`: Measurement type (e.g. "mean_amplitude", "max_peak_latency")

# Keyword Arguments
- `analysis_interval::Tuple{Real,Real}`: Time interval for measurement (default: full range)
- `baseline_interval::Union{Tuple{Real,Real},Nothing}`: Baseline interval for correction (default: nothing)
- `layout::Union{Symbol, PlotLayout}`: Plot layout — :single, :grid, or :topo (default: :single)
- `channel_selection::Function`: Channel selection predicate (default: all channels)
- `condition_selection::Function`: Condition selection predicate (default: all conditions)
- `kwargs...`: Additional arguments passed to `plot_erp`

# Overlay Types
- Peak amplitude → vertical line at peak time, amplitude label
- Peak latency → vertical line at latency, time label
- Peak-to-peak → vertical lines at both peaks
- Mean amplitude → shaded analysis interval + horizontal mean line
- Area/integral → shaded analysis interval
- Fractional latency → vertical line at fractional latency point

# Examples
```julia
dat = EegFun.read_data("erps_good.jld2")

# Mean amplitude in 300-500ms interval
plot_erp_measurements(dat, "mean_amplitude",
    analysis_interval = (0.3, 0.5),
    baseline_interval = (-0.2, 0.0))

# Peak latency
plot_erp_measurements(dat, "max_peak_latency",
    analysis_interval = (0.0, 0.8),
    baseline_interval = (-0.2, 0.0),
    channel_selection = channels([:Cz, :Pz]))

# Load from file path
plot_erp_measurements("erps_good.jld2", "max_peak_amplitude",
    analysis_interval = (0.1, 0.3))
```
"""
function plot_erp_measurements(erp::ErpData, analysis_type::String; kwargs...)
    return plot_erp_measurements([erp], analysis_type; kwargs...)
end

function plot_erp_measurements(
    filepath::String,
    analysis_type::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    kwargs...,
)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_erp_measurements(data, analysis_type; kwargs...)
    else
        files = _find_batch_files(filepath, input_dir, participant_selection)
        isempty(files) && @minimal_error "No files matching pattern '$filepath' in $input_dir"

        results = NamedTuple[]
        for file in sort(files, by = _natural_sort_key)
            file_path = joinpath(input_dir, file)
            @info "Plotting: $file"
            data = read_data(file_path)
            isnothing(data) && continue
            result = plot_erp_measurements(data, analysis_type; kwargs...)
            push!(results, result)
        end
        return results
    end
end

function plot_erp_measurements(
    erp_datasets::Vector{ErpData},
    analysis_type::String;
    analysis_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    condition_selection::Function = conditions(),
    kwargs...,
)

    # Filter by condition
    selected_indices = get_selected_conditions(erp_datasets, condition_selection)
    erp_datasets = erp_datasets[selected_indices]

    isempty(erp_datasets) && @minimal_error "No conditions selected"

    # Default analysis interval to full time range
    if isnothing(analysis_interval)
        time_data = time(erp_datasets[1])
        analysis_interval = (first(time_data), last(time_data))
    end

    # Create the base ERP plot (force interactive for marker linking)
    kwargs_dict = Dict{Symbol,Any}(kwargs)
    kwargs_dict[:interactive] = true
    if !isnothing(baseline_interval)
        kwargs_dict[:baseline_interval] = baseline_interval
    end

    result = plot_erp(
        erp_datasets;
        layout = layout,
        channel_selection = channel_selection,
        condition_selection = conditions(),  # already filtered
        (; kwargs_dict...)...,
    )

    fig = result.fig
    axes = result.axes
    line_refs = result.line_refs

    # Get selected channels
    first_data = erp_datasets[1]
    metadata_cols = meta_labels(first_data)
    all_channels = setdiff(propertynames(first_data.data), metadata_cols)
    selected_channels = all_channels[channel_selection(all_channels)]

    # Build plot_line lookup from line_refs
    plot_lines = _build_plot_line_lookup(line_refs, axes)

    # Measurement kwargs (defaults from canonical kwargs)
    measurement_kwargs = Dict{Symbol,Any}(k => v[1] for (k, v) in ERP_MEASUREMENTS_KWARGS)

    # Compute measurements and overlay for each dataset/channel
    # Pre-compute baseline-corrected data per dataset (avoid redundant copies per channel)
    baseline_cache = Dict{Int,Tuple{ErpData,Vector{Float64}}}()
    for (idx, dataset) in enumerate(erp_datasets)
        erp_data = baseline(dataset, baseline_interval)
        time_data = time(erp_data)
        baseline_cache[idx] = (erp_data, time_data)
    end

    if layout == :single
        for dataset_idx in eachindex(erp_datasets)
            erp_data, time_data = baseline_cache[dataset_idx]
            for channel in selected_channels
                plot_line = get(plot_lines[1], (dataset_idx, channel), nothing)
                _compute_and_overlay!(
                    axes[1],
                    erp_data,
                    time_data,
                    channel,
                    analysis_type,
                    analysis_interval,
                    measurement_kwargs,
                    plot_line,
                )
            end
        end
    else
        # One channel per axis (grid or topo layout)
        for (ax_idx, channel) in enumerate(selected_channels)
            if ax_idx <= length(axes)
                for dataset_idx in eachindex(erp_datasets)
                    erp_data, time_data = baseline_cache[dataset_idx]
                    plot_line = get(plot_lines[ax_idx], (dataset_idx, channel), nothing)
                    _compute_and_overlay!(
                        axes[ax_idx],
                        erp_data,
                        time_data,
                        channel,
                        analysis_type,
                        analysis_interval,
                        measurement_kwargs,
                        plot_line,
                    )
                end
            end
        end
    end

    return (fig = fig)
end


"""Build plot_line lookup dictionary from line_refs structure."""
function _build_plot_line_lookup(line_refs, axes)
    plot_lines = [Dict{Tuple{Int,Symbol},Any}() for _ in axes]
    if isnothing(line_refs)
        return plot_lines
    end
    for (ax_idx, ax_line_refs) in enumerate(line_refs)
        if ax_idx <= length(plot_lines) && !isnothing(ax_line_refs)
            for (dataset_idx, channel_lines) in ax_line_refs
                if !isnothing(channel_lines)
                    for (channel, line_data) in channel_lines
                        if line_data isa Tuple && length(line_data) >= 1
                            plot_lines[ax_idx][(dataset_idx, channel)] = line_data[1]
                        end
                    end
                end
            end
        end
    end
    return plot_lines
end


"""
Compute measurement for a single channel and draw overlay markers.
Uses pre-computed baseline-corrected data to avoid redundant copies.
"""
function _compute_and_overlay!(
    ax::Axis,
    erp_data::ErpData,
    time_data::Vector{Float64},
    channel::Symbol,
    analysis_type::String,
    analysis_interval::Tuple{Real,Real},
    measurement_kwargs::Dict{Symbol,Any},
    plot_line,
)
    channel_data = erp_data.data[!, channel]

    # Get analysis interval mask
    time_mask = (time_data .>= analysis_interval[1]) .& (time_data .<= analysis_interval[2])
    !any(time_mask) && return

    selected_data = channel_data[time_mask]
    selected_times = time_data[time_mask]
    time_min, time_max = extrema(selected_times)

    # Compute measurement value
    value = _compute_measurement(selected_data, selected_times, analysis_type, measurement_kwargs, channel)
    (isnothing(value) || isnan(value)) && return

    # Get marker visibility and color from plot line
    marker_visible, marker_color = _get_marker_style(plot_line)

    # Draw overlay based on analysis type
    _draw_measurement_overlay!(
        ax,
        analysis_type,
        value,
        selected_data,
        selected_times,
        time_data,
        channel_data,
        time_min,
        time_max,
        marker_visible,
        marker_color,
    )
end


"""Extract visibility observable and color from a plot line."""
function _get_marker_style(plot_line)
    if isnothing(plot_line)
        return Observable(true), Observable(:black)
    end
    color = plot_line.color
    marker_color = color isa Observable ? color : Observable(color)
    return plot_line.visible, marker_color
end


"""Draw measurement-specific overlay markers on an axis."""
function _draw_measurement_overlay!(
    ax,
    analysis_type,
    value,
    selected_data,
    selected_times,
    time_data,
    channel_data,
    time_min,
    time_max,
    marker_visible,
    marker_color,
)
    if analysis_type in ["max_peak_amplitude", "min_peak_amplitude"]
        # Value is amplitude; find peak position for vertical line
        peak_idx = analysis_type == "max_peak_amplitude" ? argmax(selected_data) : argmin(selected_data)
        peak_time = selected_times[peak_idx]
        peak_amp = selected_data[peak_idx]

        vlines!(ax, peak_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        text!(
            ax,
            peak_time,
            peak_amp,
            text = Printf.@sprintf("%.2f μV", value),
            align = (:center, :bottom),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type in ["max_peak_latency", "min_peak_latency"]
        # Value is latency; draw vline there
        latency = value
        if time_min <= latency <= time_max
            latency_idx = searchsortedfirst(time_data, latency)
            latency_amp = channel_data[latency_idx]

            vlines!(ax, latency, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)
            text!(
                ax,
                latency,
                latency_amp,
                text = Printf.@sprintf("%.3f s", latency),
                align = (:center, :bottom),
                color = marker_color,
                fontsize = 14,
                visible = marker_visible,
            )
        end

    elseif analysis_type in ["peak_to_peak_amplitude", "peak_to_peak_latency"]
        # Show both peaks
        max_idx = argmax(selected_data)
        min_idx = argmin(selected_data)
        max_time, min_time = selected_times[max_idx], selected_times[min_idx]
        max_amp, min_amp = selected_data[max_idx], selected_data[min_idx]

        vlines!(ax, max_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        vlines!(ax, min_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)

        label =
            analysis_type == "peak_to_peak_amplitude" ? Printf.@sprintf("P2P: %.2f μV", value) : Printf.@sprintf("P2P Lat: %.3f s", value)
        text!(
            ax,
            (max_time + min_time) / 2,
            (max_amp + min_amp) / 2,
            text = label,
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type == "mean_amplitude"
        # Shade analysis interval + horizontal mean line
        y_min, y_max = extrema(selected_data)
        concrete_color = marker_color isa Observable ? marker_color[] : marker_color
        rect = [Point2f(time_min, y_min), Point2f(time_max, y_min), Point2f(time_max, y_max), Point2f(time_min, y_max)]
        poly!(ax, rect, color = (concrete_color, 0.3), strokewidth = 0, visible = marker_visible)
        lines!(ax, [time_min, time_max], [value, value], color = concrete_color, linewidth = 2, linestyle = :dash, visible = marker_visible)
        text!(
            ax,
            (time_min + time_max) / 2,
            value,
            text = Printf.@sprintf("Mean: %.2f μV", value),
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type in ["rectified_area", "integral", "positive_area", "negative_area"]
        # Shade analysis interval
        y_min, y_max = extrema(selected_data)
        concrete_color = marker_color isa Observable ? marker_color[] : marker_color
        rect = [Point2f(time_min, y_min), Point2f(time_max, y_min), Point2f(time_max, y_max), Point2f(time_min, y_max)]
        poly!(ax, rect, color = (concrete_color, 0.15), strokewidth = 0, visible = marker_visible)
        text!(
            ax,
            (time_min + time_max) / 2,
            (y_min + y_max) / 2,
            text = Printf.@sprintf("%s: %.2f μVs", analysis_type, value),
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type == "fractional_area_latency"
        latency = value
        if time_min <= latency <= time_max
            vlines!(ax, latency, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)
            latency_idx = searchsortedfirst(time_data, latency)
            latency_amp = channel_data[latency_idx]
            text!(
                ax,
                latency,
                latency_amp,
                text = Printf.@sprintf("Frac Area: %.3f s", latency),
                align = (:center, :bottom),
                color = marker_color,
                fontsize = 14,
                visible = marker_visible,
            )
        end

    elseif analysis_type == "fractional_peak_latency"
        latency = value
        if time_min <= latency <= time_max
            vlines!(ax, latency, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)
            latency_idx = searchsortedfirst(time_data, latency)
            latency_amp = channel_data[latency_idx]
            text!(
                ax,
                latency,
                latency_amp,
                text = Printf.@sprintf("Frac Peak: %.3f s", latency),
                align = (:center, :bottom),
                color = marker_color,
                fontsize = 14,
                visible = marker_visible,
            )
        end
    end
end
