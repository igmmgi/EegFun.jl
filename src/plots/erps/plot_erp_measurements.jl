"""
    plot_erp_measurements(erp_data, analysis_type;
                          analysis_interval = nothing,
                          baseline_interval = nothing,
                          layout = :single,
                          channel_selection = channels(),
                          condition_selection = conditions(),
                          average_channels = false,
                          kwargs...)

Plot ERP data with measurement overlays computed inline.

Computes a single measurement type across all selected channels and conditions,
then overlays the results on the ERP waveforms.

# Arguments
- `erp_data::Union{String, ErpData, AbstractVector{<:ErpData}}`: ERP data (filepath, pattern, or data object)
- `analysis_type::String`: Measurement type (e.g. `"mean_amplitude"`, `"max_peak_latency"`)

# Keyword Arguments
- `analysis_interval::Interval`: Time interval for measurement (default: full range)
- `baseline_interval::Interval`: Baseline interval for correction (default: `nothing`)
- `layout::Union{Symbol, PlotLayout}`: Plot layout — `:single`, `:grid`, or `:topo` (default: `:single`)
- `channel_selection`: Channel selection predicate (default: `channels()`)
- `condition_selection`: Condition selection predicate (default: `conditions()`)
- `average_channels::Bool`: Evaluate measurement across spatially pooled ROI waveform (default: `false`)
- `kwargs...`: Additional arguments passed to `plot_erp` and `erp_measurements`

$(_generate_kwargs_doc(ERP_MEASUREMENTS_KWARGS))

# Examples
```julia
plot_erp_measurements(dat, "mean_amplitude", analysis_interval = (0.3, 0.5), baseline_interval = (-0.2, 0.0))
plot_erp_measurements("erps_final.jld2", "max_peak_amplitude", analysis_interval = (0.1, 0.3))
```
"""
plot_erp_measurements(erp::ErpData, analysis_type::String; kwargs...) = plot_erp_measurements([erp], analysis_type; kwargs...)

function plot_erp_measurements(
    filepath::String,
    analysis_type::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    kwargs...,
)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_erp_measurements(data, analysis_type; kwargs...)
    end

    files = _require_batch_files(filepath, input_dir, participant_selection)
    results = NamedTuple[]
    for file in files
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(results, plot_erp_measurements(data, analysis_type; kwargs...))
    end
    return results
end

function plot_erp_measurements(
    erp_datasets::AbstractVector{<:ErpData},
    analysis_type::String;
    analysis_interval::Interval = nothing,
    baseline_interval::Interval = nothing,
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection = channels(),
    condition_selection = conditions(),
    average_channels::Bool = false,
    kwargs...,
)
    selected_indices = get_selected_conditions(erp_datasets, condition_selection)
    erp_datasets = erp_datasets[selected_indices]
    isempty(erp_datasets) && @minimal_error "No conditions selected"

    analysis_interval = if isnothing(analysis_interval)
        t0 = time_vector(erp_datasets[1])
        (first(t0), last(t0))
    else
        _normalize_interval(analysis_interval)
    end

    kwargs_dict = _kwargs_dict(kwargs)
    measurement_kwargs = Dict{Symbol,Any}(k => (haskey(kwargs_dict, k) ? pop!(kwargs_dict, k) : v[1]) for (k, v) in ERP_MEASUREMENTS_KWARGS)

    kwargs_dict[:interactive] = true
    !isnothing(baseline_interval) && (kwargs_dict[:baseline_interval] = baseline_interval)
    average_channels && (kwargs_dict[:average_channels] = true)

    result = plot_erp(
        erp_datasets;
        layout = layout,
        channel_selection = channel_selection,
        condition_selection = conditions(),
        (; kwargs_dict...)...,
    )

    fig = result.fig
    axes = result.axes
    plot_lines = _build_plot_line_lookup(result.line_refs, axes)

    first_data = erp_datasets[1]
    all_channels = setdiff(propertynames(first_data.data), meta_labels(first_data))
    chan_func = _to_channel_func(channel_selection)
    selected_channels = all_channels[chan_func(all_channels)]
    plot_channels = average_channels ? [:avg] : selected_channels

    n_ds = length(erp_datasets)
    cached_erps = Vector{ErpData}(undef, n_ds)
    cached_times = Vector{Vector{Float64}}(undef, n_ds)
    ch_sel_avg = average_channels ? [channels(selected_channels)] : nothing
    for (idx, dataset) in enumerate(erp_datasets)
        erp_data = isnothing(baseline_interval) ? (average_channels ? copy(dataset) : dataset) : baseline(dataset, baseline_interval)
        average_channels && channel_average!(erp_data; channel_selections = ch_sel_avg, output_labels = [:avg], reduce = false)
        cached_erps[idx] = erp_data
        cached_times[idx] = time_vector(erp_data)
    end

    if layout == :single
        for dataset_idx = 1:n_ds, channel in plot_channels
            _compute_and_overlay!(
                axes[1],
                cached_erps[dataset_idx],
                cached_times[dataset_idx],
                channel,
                analysis_type,
                analysis_interval,
                measurement_kwargs,
                get(plot_lines[1], (dataset_idx, channel), nothing),
            )
        end
    else
        for (ax_idx, channel) in enumerate(plot_channels)
            ax_idx > length(axes) && break
            for dataset_idx = 1:n_ds
                _compute_and_overlay!(
                    axes[ax_idx],
                    cached_erps[dataset_idx],
                    cached_times[dataset_idx],
                    channel,
                    analysis_type,
                    analysis_interval,
                    measurement_kwargs,
                    get(plot_lines[ax_idx], (dataset_idx, channel), nothing),
                )
            end
        end
    end

    return (fig = fig, axes = axes)
end

function _build_plot_line_lookup(line_refs, axes)
    plot_lines = [Dict{Tuple{Int,Symbol},Any}() for _ in axes]
    isnothing(line_refs) && return plot_lines
    for (ax_idx, ax_line_refs) in enumerate(line_refs)
        (ax_idx > length(plot_lines) || isnothing(ax_line_refs)) && continue
        for (dataset_idx, channel_lines) in ax_line_refs
            isnothing(channel_lines) && continue
            for (channel, line_data) in channel_lines
                line_data isa Tuple && !isempty(line_data) && (plot_lines[ax_idx][(dataset_idx, channel)] = line_data[1])
            end
        end
    end
    return plot_lines
end

function _select_interval_slice(data_vec::AbstractVector, time_data::Vector{Float64}, t_lo::Real, t_hi::Real)
    if issorted(time_data)
        i1 = searchsortedfirst(time_data, t_lo)
        i2 = searchsortedlast(time_data, t_hi)
        i1 > i2 && return nothing
        return @view(data_vec[i1:i2]), @view(time_data[i1:i2])
    end
    idxs = findall(t -> t_lo <= t <= t_hi, time_data)
    isempty(idxs) && return nothing
    return @view(data_vec[idxs]), @view(time_data[idxs])
end

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
    sliced = _select_interval_slice(channel_data, time_data, analysis_interval[1], analysis_interval[2])
    isnothing(sliced) && return
    selected_data, selected_times = sliced
    time_min, time_max = extrema(selected_times)

    value = _compute_measurement(selected_data, selected_times, analysis_type, measurement_kwargs, channel)
    (isnothing(value) || isnan(value)) && return

    marker_visible, marker_color = if isnothing(plot_line)
        Observable(true), Observable(:black)
    else
        col = plot_line.color
        plot_line.visible, (col isa Observable ? col : Observable(col))
    end

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
    if analysis_type in ("max_peak_amplitude", "min_peak_amplitude")
        peak_idx = analysis_type == "max_peak_amplitude" ? argmax(selected_data) : argmin(selected_data)
        peak_time = selected_times[peak_idx]
        vlines!(ax, peak_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        text!(
            ax,
            peak_time,
            selected_data[peak_idx],
            text = Printf.@sprintf("%.2f μV", value),
            align = (:center, :bottom),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )
    elseif analysis_type in ("max_peak_latency", "min_peak_latency", "fractional_area_latency", "fractional_peak_latency")
        if time_min <= value <= time_max
            vlines!(ax, value, color = marker_color, linewidth = 2, linestyle = :solid, visible = marker_visible)
            latency_amp = channel_data[clamp(searchsortedfirst(time_data, value), 1, length(channel_data))]
            lbl = if analysis_type == "fractional_area_latency"
                Printf.@sprintf("Frac Area: %.3f s", value)
            elseif analysis_type == "fractional_peak_latency"
                Printf.@sprintf("Frac Peak: %.3f s", value)
            else
                Printf.@sprintf("%.3f s", value)
            end
            text!(
                ax,
                value,
                latency_amp,
                text = lbl,
                align = (:center, :bottom),
                color = marker_color,
                fontsize = 14,
                visible = marker_visible,
            )
        end
    elseif analysis_type in ("peak_to_peak_amplitude", "peak_to_peak_latency")
        max_idx, min_idx = argmax(selected_data), argmin(selected_data)
        max_time, min_time = selected_times[max_idx], selected_times[min_idx]
        vlines!(ax, max_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        vlines!(ax, min_time, color = marker_color, linewidth = 1, linestyle = :solid, visible = marker_visible)
        lbl = analysis_type == "peak_to_peak_amplitude" ? Printf.@sprintf("P2P: %.2f μV", value) : Printf.@sprintf("P2P Lat: %.3f s", value)
        text!(
            ax,
            (max_time + min_time) / 2,
            (selected_data[max_idx] + selected_data[min_idx]) / 2,
            text = lbl,
            align = (:center, :center),
            color = marker_color,
            fontsize = 14,
            visible = marker_visible,
        )
    elseif analysis_type == "mean_amplitude"
        y_min, y_max = extrema(selected_data)
        concrete_color = marker_color isa Observable ? marker_color[] : marker_color
        poly!(
            ax,
            Point2f[(time_min, y_min), (time_max, y_min), (time_max, y_max), (time_min, y_max)],
            color = (concrete_color, 0.3),
            strokewidth = 0,
            visible = marker_visible,
        )
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
    elseif analysis_type in ("rectified_area", "integral", "positive_area", "negative_area")
        y_min, y_max = extrema(selected_data)
        concrete_color = marker_color isa Observable ? marker_color[] : marker_color
        poly!(
            ax,
            Point2f[(time_min, y_min), (time_max, y_min), (time_max, y_max), (time_min, y_max)],
            color = (concrete_color, 0.15),
            strokewidth = 0,
            visible = marker_visible,
        )
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
    end
end
