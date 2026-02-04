"""
    plot_erp_with_measurements(measurements::ErpMeasurementsResult, erp_data::Union{String, ErpData, Vector{ErpData}};
                               analysis_window::Union{Function,Nothing} = nothing,
                               layout::Union{Symbol, PlotLayout} = :single,
                               channel_selection::Function = channels(),
                               condition_selection::Function = conditions(),
                               participant_selection::Function = participants(),
                               kwargs...)

    plot_erp_with_measurements(measurements_df::DataFrame, erp_data::Union{String, ErpData, Vector{ErpData}}, analysis_type::String;
                               analysis_window::Function = samples(),
                               layout::Union{Symbol, PlotLayout} = :single,
                               channel_selection::Function = channels(),
                               condition_selection::Function = conditions(),
                               participant_selection::Function = participants(),
                               kwargs...)

Plot ERP data with measurement overlays from erp_measurements results.

This function takes measurement results (from `erp_measurements`) and overlays
the corresponding measurement markers on ERP plots. The type of overlay depends on `analysis_type`:
- Peak measurements: vertical lines at peak latencies with amplitude labels
- Mean amplitude: shaded region showing analysis window with mean value label
- Area measurements: shaded region showing analysis window
- Fractional latencies: vertical lines at fractional latency points

# Arguments
- `measurements::ErpMeasurementsResult`: Results from `erp_measurements` (preferred - automatically uses correct analysis_window)
- `measurements_df::DataFrame`: DataFrame from `erp_measurements` (legacy - requires analysis_type parameter)
- `erp_data::Union{String, ErpData, Vector{ErpData}}`: ERP data to plot (filepath or data object)
- `analysis_type::String`: Type of measurement to overlay (only needed for DataFrame input)
- `analysis_window::Union{Function,Nothing}`: Analysis window predicate (optional - uses value from ErpMeasurementsResult if not provided)
- `layout::Union{Symbol, PlotLayout}`: Plot layout (same as `plot_erp`)
- `channel_selection::Function`: Channel selection predicate
- `condition_selection::Function`: Condition selection predicate
- `participant_selection::Function`: Participant selection predicate
- `kwargs...`: Additional plotting arguments (passed to `plot_erp`)

# Examples
```julia
# Get measurements (returns ErpMeasurementsResult)
results = erp_measurements("erps_good", "max_peak_amplitude", 
                          analysis_window = samples((0.1, 0.3)))

# Plot with measurement overlays (automatically uses correct analysis_window)
plot_erp_with_measurements(results, "erps_good.jld2")

# Access DataFrame directly
results.data  # or just use results as a DataFrame
```
"""
function plot_erp_with_measurements(
    measurements::ErpMeasurementsResult,
    erp_data::Union{String,ErpData,Vector{ErpData}};
    analysis_window::Union{Function,Nothing} = nothing,
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    condition_selection::Function = conditions(),
    participant_selection::Function = participants(),
    kwargs...,
)
    # Use analysis_window from measurements if not provided
    analysis_window = something(analysis_window, measurements.analysis_window)

    # Extract DataFrame and analysis type from measurements
    measurements_df = measurements.data
    analysis_type = measurements.analysis_type

    return _plot_erp_with_measurements_impl(
        measurements_df,
        erp_data,
        analysis_type,
        analysis_window,
        layout,
        channel_selection,
        condition_selection,
        participant_selection,
        kwargs...,
    )
end

# Backward compatibility: allow DataFrame input (but require analysis_type)
function plot_erp_with_measurements(
    measurements_df::DataFrame,
    erp_data::Union{String,ErpData,Vector{ErpData}},
    analysis_type::String;
    analysis_window::Function = samples(),
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    condition_selection::Function = conditions(),
    participant_selection::Function = participants(),
    kwargs...,
)
    return _plot_erp_with_measurements_impl(
        measurements_df,
        erp_data,
        analysis_type,
        analysis_window,
        layout,
        channel_selection,
        condition_selection,
        participant_selection,
        kwargs...,
    )
end

# Internal implementation function
function _plot_erp_with_measurements_impl(
    measurements_df::DataFrame,
    erp_data::Union{String,ErpData,Vector{ErpData}},
    analysis_type::String,
    analysis_window::Function,
    layout::Union{Symbol,PlotLayout},
    channel_selection::Function,
    condition_selection::Function,
    participant_selection::Function,
    kwargs...,
)
    # Force interactive=true since we need line_refs for marker visibility linking
    kwargs_dict = Dict{Symbol,Any}(kwargs)
    kwargs_dict[:interactive] = true
    kwargs_filtered = (; kwargs_dict...)

    # Load ERP data if filepath provided
    if erp_data isa String
        data = read_data(erp_data)
        if isnothing(data)
            @minimal_error_throw "No data found in file: $erp_data"
        end
        # Dispatch will handle ErpData vs Vector{ErpData} automatically
        if data isa ErpData
            erp_datasets = [data]
        elseif data isa Vector{<:ErpData}
            erp_datasets = data
        else
            @minimal_error_throw "Expected ErpData or Vector{ErpData}, got $(typeof(data))"
        end
    elseif erp_data isa ErpData
        erp_datasets = [erp_data]
    else
        erp_datasets = erp_data
    end

    # Filter datasets by condition and participant
    # Extract participant IDs from filenames or use all if measurements don't have participant column
    if hasproperty(measurements_df, :participant)
        unique_participants = unique(measurements_df.participant)
        participant_mask = participant_selection(unique_participants)
        selected_participants = unique_participants[participant_mask]
    else
        selected_participants = nothing
    end

    # Filter by condition
    selected_indices = get_selected_conditions(erp_datasets, condition_selection)
    erp_datasets = erp_datasets[selected_indices]

    # Create the base plot
    # interactive=true ensures line_refs is populated for marker linking
    # plot_erp returns (fig=fig, axes=axes, line_refs=line_refs)
    result = plot_erp(
        erp_datasets;
        layout = layout,
        channel_selection = channel_selection,
        condition_selection = conditions(),  # Already filtered above
        kwargs_filtered...,  # Contains interactive=true
    )

    fig = result.fig
    axes = result.axes
    line_refs = result.line_refs

    # Check if line_refs is available
    if line_refs === nothing
        @minimal_warning "line_refs is nothing - markers won't link to plot line visibility. Ensure interactive=true."
    end

    # Get selected channels from first dataset
    first_data = erp_datasets[1]
    metadata_cols = meta_labels(first_data)
    all_channels = setdiff(propertynames(first_data.data), metadata_cols)
    channel_mask = channel_selection(all_channels)
    selected_channels = all_channels[channel_mask]

    # Extract plot lines from line_refs structure
    # Structure: line_refs[ax_idx][dataset_idx][channel] = (line, y_obs)
    plot_lines = [Dict{Tuple{Int,Symbol},Any}() for _ in axes]
    if line_refs !== nothing
        for (ax_idx, ax_line_refs) in enumerate(line_refs)
            if ax_idx <= length(plot_lines) && ax_line_refs !== nothing
                for (dataset_idx, channel_lines) in ax_line_refs
                    if channel_lines !== nothing
                        for (channel, line_data) in channel_lines
                            if line_data isa Tuple && length(line_data) >= 1
                                line = line_data[1]  # Extract the line object
                                key = (dataset_idx, channel)
                                plot_lines[ax_idx][key] = line
                            end
                        end
                    end
                end
            end
        end
    else
        @minimal_warning "line_refs is nothing - markers won't be linked to plot lines. Ensure interactive=true."
    end

    # Determine which channels are plotted on which axes
    # For :single layout, all channels on one axis
    # For :grid/:topo, one channel per axis
    if layout == :single
        # All channels on single axis
        for (dataset_idx, dataset) in enumerate(erp_datasets)
            for channel in selected_channels
                plot_line = get(plot_lines[1], (dataset_idx, channel), nothing)
                if plot_line === nothing
                    @debug "No plot line found for dataset $(dataset.condition_name), channel $channel on axis 1"
                end
                _overlay_measurements!(axes[1], measurements_df, dataset, channel, analysis_type, analysis_window, plot_line)
            end
        end
    else
        # One channel per axis (grid or topo layout)
        # Each axis can have multiple conditions plotted on it, so overlay for all conditions
        for (ax_idx, channel) in enumerate(selected_channels)
            if ax_idx <= length(axes)
                # Loop through all datasets (conditions) for this channel
                for (dataset_idx, dataset) in enumerate(erp_datasets)
                    plot_line = get(plot_lines[ax_idx], (dataset_idx, channel), nothing)
                    if plot_line === nothing
                        @debug "No plot line found for dataset $(dataset.condition_name), channel $channel on axis $ax_idx"
                    end
                    _overlay_measurements!(axes[ax_idx], measurements_df, dataset, channel, analysis_type, analysis_window, plot_line)
                end
            end
        end
    end

    return fig
end

"""
Overlay measurement markers on an axis.
"""
function _overlay_measurements!(
    ax::Axis,
    measurements_df::DataFrame,
    dataset::ErpData,
    channel::Symbol,
    analysis_type::String,
    analysis_window,
    plot_line::Union{Any,Nothing} = nothing,
)
    # Find matching row in measurements DataFrame
    # Match by condition and channel
    # Also check for participant if available
    condition_mask = measurements_df.condition .== dataset.condition

    # If measurements_df has participant column and dataset has participant info, filter by that too
    if hasproperty(measurements_df, :participant) && hasproperty(dataset, :participant)
        participant_mask = measurements_df.participant .== dataset.participant
        condition_mask = condition_mask .& participant_mask
    end

    # Check for valid channel data
    channel_mask = hasproperty(measurements_df, channel) .& .!isnan.(measurements_df[!, channel])

    matching_rows = measurements_df[condition_mask .& channel_mask, :]

    if isempty(matching_rows)
        @debug "No matching measurements for condition $(dataset.condition), channel $channel"
        return  # No measurements for this condition/channel
    end

    # Get measurement value for this channel
    if !hasproperty(matching_rows, channel)
        @debug "Channel $channel not found in matching rows for condition $(dataset.condition)"
        return
    end

    measurement_value = matching_rows[1, channel]
    if isnan(measurement_value) || isinf(measurement_value)
        return
    end

    # Get time data
    time_data = dataset.data[!, :time]
    channel_data = dataset.data[!, channel]

    # Apply analysis window to get time range
    time_mask = analysis_window(dataset.data)
    if any(time_mask)
        time_range = extrema(time_data[time_mask])
        time_min, time_max = time_range
    else
        return  # No valid time window
    end

    # Get visibility and color from plot line
    # Pass these directly to marker properties - Makie handles this efficiently
    marker_visible = if plot_line !== nothing
        try
            plot_line.visible  # This should be an Observable that legend clicks modify
        catch
            Observable(true)
        end
    else
        Observable(true)
    end

    # Get color from plot line (may be Observable or direct value)
    marker_color = if plot_line !== nothing
        try
            color_prop = plot_line.color
            # If it's an Observable, use it directly; otherwise wrap it
            color_prop isa Observable ? color_prop : Observable(color_prop)
        catch
            Observable(:black)
        end
    else
        Observable(:black)
    end

    # Overlay based on measurement type
    if analysis_type in ["max_peak_amplitude", "min_peak_amplitude"]
        # Find peak latency in the analysis window
        window_data = channel_data[time_mask]
        window_times = time_data[time_mask]

        if analysis_type == "max_peak_amplitude"
            peak_idx = argmax(window_data)
        else
            peak_idx = argmin(window_data)
        end

        peak_time = window_times[peak_idx]
        peak_amp = window_data[peak_idx]

        # Draw vertical line at peak with linked visibility and color
        vlines!(ax, peak_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)

        # Add text label with linked visibility and color
        text!(
            ax,
            peak_time,
            peak_amp,
            text = Printf.@sprintf("%.2f μV", measurement_value),
            align = (:center, :bottom),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type in ["max_peak_latency", "min_peak_latency"]
        # For latency measurements, the value IS the latency
        latency = measurement_value

        if time_min <= latency <= time_max
            # Draw vertical line at latency with linked visibility and color
            vlines!(ax, latency, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)

            # Get amplitude at this latency
            latency_idx = argmin(abs.(time_data .- latency))
            latency_amp = channel_data[latency_idx]

            # Add text label with linked visibility and color
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

    elseif analysis_type == "peak_to_peak_latency"
        # peak_to_peak_latency is the DURATION between max and min peaks, not a single latency
        # We need to find both peaks and draw lines at both times
        window_data = channel_data[time_mask]
        window_times = time_data[time_mask]

        # Find max and min peaks in the analysis window
        max_idx = argmax(window_data)
        min_idx = argmin(window_data)

        max_time = window_times[max_idx]
        min_time = window_times[min_idx]
        max_amp = window_data[max_idx]
        min_amp = window_data[min_idx]

        # Draw vertical lines at both peak times with linked visibility and color
        vlines!(ax, max_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        vlines!(ax, min_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)

        # Add text label showing the duration (the measurement value)
        mid_time = (max_time + min_time) / 2
        mid_amp = (max_amp + min_amp) / 2
        text!(
            ax,
            mid_time,
            mid_amp,
            text = Printf.@sprintf("P2P Lat: %.3f s", measurement_value),
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type == "mean_amplitude"
        # Shade analysis window and show mean value with linked visibility and color
        # Get the actual data range in the window to make the band visible
        window_data = channel_data[time_mask]
        y_min, y_max = extrema(window_data)

        # Get concrete color value (not Observable) for band and hlines
        # band! and hlines! may not work well with Observable colors
        concrete_color = marker_color isa Observable ? marker_color[] : marker_color

        # Shade the full data range in the analysis window using a rectangle
        # Create rectangle vertices: bottom-left, bottom-right, top-right, top-left
        rect_vertices = [Point2f(time_min, y_min), Point2f(time_max, y_min), Point2f(time_max, y_max), Point2f(time_min, y_max)]

        # Use poly! to create a filled rectangle
        # Pass Observable directly for visibility (same as vlines! and text!)
        poly!(ax, rect_vertices, color = concrete_color, alpha = 0.3, strokewidth = 0, visible = marker_visible)

        # Draw a horizontal line at the mean value
        # Pass Observable directly for visibility
        hlines!(
            ax,
            measurement_value,
            xmin = time_min,
            xmax = time_max,
            color = concrete_color,
            linewidth = 2,
            linestyle = :dash,
            visible = marker_visible,
        )

        # Add text label at center of window with linked visibility and color
        mid_time = (time_min + time_max) / 2
        text!(
            ax,
            mid_time,
            measurement_value,
            text = Printf.@sprintf("Mean: %.2f μV", measurement_value),
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type in ["rectified_area", "integral", "positive_area", "negative_area"]
        # Shade analysis window with linked visibility and color
        y_min, y_max = extrema(channel_data[time_mask])

        # Get concrete color value (not Observable) for poly
        concrete_color = marker_color isa Observable ? marker_color[] : marker_color

        # Create rectangle vertices: bottom-left, bottom-right, top-right, top-left
        rect_vertices = [Point2f(time_min, y_min), Point2f(time_max, y_min), Point2f(time_max, y_max), Point2f(time_min, y_max)]

        # Use poly! to create a filled rectangle (same approach as mean_amplitude)
        poly!(ax, rect_vertices, color = concrete_color, alpha = 0.15, strokewidth = 0, visible = marker_visible)

        # Add text label with linked visibility and color
        mid_time = (time_min + time_max) / 2
        mid_y = (y_min + y_max) / 2
        text!(
            ax,
            mid_time,
            mid_y,
            text = Printf.@sprintf("%s: %.2f μVs", analysis_type, measurement_value),
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )

    elseif analysis_type == "fractional_area_latency"
        # Fractional area latency: time point where a fraction of the total area has been accumulated
        latency = measurement_value

        if time_min <= latency <= time_max
            # Draw vertical line at fractional latency
            # Use solid line with condition-specific color (from marker_color)
            vlines!(ax, latency, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)

            # Get amplitude at this latency for label placement
            latency_idx = argmin(abs.(time_data .- latency))
            latency_amp = channel_data[latency_idx]

            # More descriptive label explaining what this represents
            # Include condition name if available to distinguish between conditions
            label_text = if hasproperty(dataset, :condition_name) && !isempty(dataset.condition_name)
                Printf.@sprintf("%s: %.3f s", dataset.condition_name, latency)
            else
                Printf.@sprintf("Frac Area: %.3f s", latency)
            end

            text!(
                ax,
                latency,
                latency_amp,
                text = label_text,
                align = (:center, :bottom),
                color = marker_color,
                fontsize = 14,
                visible = marker_visible,
            )
        end

    elseif analysis_type == "fractional_peak_latency"
        # Fractional peak latency: time point where signal reaches a fraction of peak amplitude
        latency = measurement_value

        if time_min <= latency <= time_max
            # Find the peak in the window to show context
            # Use robust peak detection to match what was used in the measurement
            window_data = channel_data[time_mask]
            window_times = time_data[time_mask]

            # Find both max and min peaks, use the one with larger absolute value
            max_idx = argmax(window_data)
            min_idx = argmin(window_data)
            max_amp = window_data[max_idx]
            min_amp = window_data[min_idx]

            if abs(max_amp) >= abs(min_amp)
                peak_idx = max_idx
                peak_amp = max_amp
            else
                peak_idx = min_idx
                peak_amp = min_amp
            end
            peak_time = window_times[peak_idx]

            # Get concrete color for peak line (use same color as plot line)
            concrete_color = marker_color isa Observable ? marker_color[] : marker_color

            # Draw vertical line at peak location (for context) - use solid line with same color
            # Make sure it's visible and distinct from fractional latency line
            vlines!(ax, peak_time, color = concrete_color, linewidth = 1.5, linestyle = :solid, alpha = 0.7, visible = marker_visible)

            # Draw vertical line at fractional latency - use dashed line to distinguish
            vlines!(ax, latency, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)

            # Get amplitude at fractional latency for label placement
            latency_idx = argmin(abs.(time_data .- latency))
            latency_amp = channel_data[latency_idx]

            # More descriptive label showing both peak and fractional latency
            if abs(latency - peak_time) > 0.01
                text!(
                    ax,
                    latency,
                    latency_amp,
                    text = Printf.@sprintf("Frac: %.3f s\nPeak: %.3f s", latency, peak_time),
                    align = (:center, :bottom),
                    color = marker_color,
                    fontsize = 11,
                    visible = marker_visible,
                )
            else
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

    elseif analysis_type == "peak_to_peak_amplitude"
        # Show both peaks with linked visibility and color
        window_data = channel_data[time_mask]
        window_times = time_data[time_mask]

        max_idx = argmax(window_data)
        min_idx = argmin(window_data)

        max_time = window_times[max_idx]
        min_time = window_times[min_idx]
        max_amp = window_data[max_idx]
        min_amp = window_data[min_idx]

        # Draw lines at both peaks with linked visibility and color
        vlines!(ax, max_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        vlines!(ax, min_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)

        # Add text label with linked visibility and color
        text!(
            ax,
            (max_time + min_time) / 2,
            (max_amp + min_amp) / 2,
            text = Printf.@sprintf("P2P: %.2f μV", measurement_value),
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )
    end
end
