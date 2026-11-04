# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_TRIGGERS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :window_size => (10.0, "Size of the time window to display in seconds."),
    :initial_position => (-2.0, "Initial position of the time interval in seconds."),
    :min_window_size => (0.1, "Minimum size of the time interval in seconds."),
    :max_window_size => (100.0, "Maximum size of the time interval in seconds."),
    :label_fontsize => (20, "Font size for labels."),
    :display_plot => (true, "Whether to display the plot."),
    :plot_title => (nothing, "Custom plot title override."),
    :show_title => (true, "Whether to show the plot title."),
    :legend => (true, "Whether to show the trigger count legend."),
    :ignore_triggers => (Int[], "Vector of trigger codes to ignore."),
    :marker_size => (15, "Size of markers in the trigger overview plot."),
    :linewidth_trigger => (1, "Line width for trigger vertical lines."),
    :line_offset => (0.1, "Offset for trigger vertical lines."),
    :timeline_width => (2, "Line width for the main timeline."),
    :event_line_width => (1, "Line width for individual event lines."),
    :font_size => (24, "Font size for text elements."),
    :event_line_height => (0.05, "Height of event lines."),
    :time_label_offset => (-0.001, "Vertical offset for time labels."),
    :y_min_limit => (-0.1, "Minimum y-axis limit."),
    :y_max_limit => (0.15, "Maximum y-axis limit."),
    :window_size_step => (1.0, "Step size for window size slider."),
    :position_step => (0.5, "Step size for position slider."),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
)

function _filter_triggers(trigger_times, trigger_values, ignore_triggers)
    isempty(ignore_triggers) && return trigger_times, trigger_values
    ignore_set = Set(ignore_triggers)
    keep_mask = .!in.(trigger_values, Ref(ignore_set))
    return trigger_times[keep_mask], trigger_values[keep_mask]
end

function _filter_triggers(trigger_times, trigger_codes, trigger_info, ignore_triggers)
    isempty(ignore_triggers) && return trigger_times, trigger_codes, trigger_info
    ignore_set = Set(ignore_triggers)
    keep_mask = .!in.(trigger_codes, Ref(ignore_set))
    return trigger_times[keep_mask], trigger_codes[keep_mask], trigger_info[keep_mask]
end

function _trigger_time_count(time, triggers, ignore_triggers = Int[])
    trigger_times = Float64[]
    trigger_values = Int[]
    has_ignore = !isempty(ignore_triggers)
    ignore_set = has_ignore ? Set(ignore_triggers) : nothing

    @inbounds for i in eachindex(triggers, time)
        v = triggers[i]
        if v != 0 && (!has_ignore || !(v in ignore_set))
            push!(trigger_times, Float64(time[i]))
            push!(trigger_values, Int(v))
        end
    end
    isempty(trigger_values) && return trigger_times, trigger_values, OrderedDict{Int,Int}()

    counts = Dict{Int,Int}()
    for v in trigger_values
        counts[v] = get(counts, v, 0) + 1
    end
    trigger_count = OrderedDict{Int,Int}()
    for k in sort!(collect(keys(counts)))
        trigger_count[k] = counts[k]
    end
    return trigger_times, trigger_values, trigger_count
end

function _trigger_time_count(times::Vector{Float64}, trigger_codes::Vector{Int16}, trigger_info::Vector{String}, ignore_triggers = Int[])
    if !isempty(ignore_triggers)
        times, trigger_codes, trigger_info = _filter_triggers(times, trigger_codes, trigger_info, ignore_triggers)
    end
    trigger_count = OrderedDict{Int,Int}()
    for val in trigger_codes
        k = Int(val)
        trigger_count[k] = get(trigger_count, k, 0) + 1
    end
    return times, Int.(trigger_codes), trigger_count
end

function _plot_trigger_vertical_lines!(
    ax::Axis,
    times::Vector{Float64},
    y_positions::Vector{Int};
    offset::Real = PLOT_TRIGGERS_KWARGS[:line_offset][1],
    linewidth::Real = PLOT_TRIGGERS_KWARGS[:linewidth_trigger][1],
)
    isempty(times) && return
    n = length(times)
    xs = Vector{Float64}(undef, 3 * n)
    ys = Vector{Float64}(undef, 3 * n)
    @inbounds for i = 1:n
        idx = 3 * i - 2
        xpos = times[i]
        ypos = y_positions[i]
        xs[idx] = xpos
        xs[idx+1] = xpos
        xs[idx+2] = NaN
        ys[idx] = ypos - offset
        ys[idx+1] = ypos + offset
        ys[idx+2] = NaN
    end
    lines!(ax, xs, ys, color = _resolve_theme_linecolor(ax, nothing), linewidth = linewidth)
end

function _extract_trigger_cols(trig_col::AbstractVector, time_col::AbstractVector, info_col::Union{Nothing,AbstractVector}, ignore_triggers)
    has_ignore = !isempty(ignore_triggers)
    ignore_set = has_ignore ? Set(ignore_triggers) : nothing
    has_info_col = !isnothing(info_col)

    trigger_codes = Int16[]
    trigger_times = Float64[]
    trigger_info = String[]
    @inbounds for i in eachindex(trig_col, time_col)
        v = trig_col[i]
        if v != 0
            c = Int16(v)
            if !has_ignore || !(c in ignore_set)
                push!(trigger_codes, c)
                push!(trigger_times, Float64(time_col[i]))
                push!(trigger_info, has_info_col ? String(info_col[i]) : "")
            end
        end
    end
    return trigger_codes, trigger_times, trigger_info
end

function _extract_trigger_data(dat::ContinuousData, ignore_triggers = Int[])
    hasproperty(dat.data, :trigger) || @minimal_error("No triggers column found in data. Expected column name: trigger")
    info_col = hasproperty(dat.data, :trigger_info) ? dat.data[!, :trigger_info] : nothing
    return _extract_trigger_cols(dat.data[!, :trigger], dat.data[!, :time], info_col, ignore_triggers)
end

function _count_triggers(trigger_codes::Vector{Int16})
    counts = Dict{Int,Int}()
    for c in trigger_codes
        k = Int(c)
        counts[k] = get(counts, k, 0) + 1
    end
    trigger_count = OrderedDict{Int,Int}()
    for k in sort!(collect(keys(counts)))
        trigger_count[k] = counts[k]
    end
    return trigger_count
end

function _count_triggers(trigger_codes::Vector{Int16}, trigger_info::Vector{String})
    counts = Dict{Int,Int}()
    first_info = Dict{Int,String}()
    n_info = length(trigger_info)
    @inbounds for (idx, c) in enumerate(trigger_codes)
        k = Int(c)
        counts[k] = get(counts, k, 0) + 1
        if !haskey(first_info, k)
            first_info[k] = idx <= n_info ? trigger_info[idx] : ""
        end
    end
    trigger_count = OrderedDict{Int,Int}()
    trigger_labels = OrderedDict{Int,String}()
    for k in sort!(collect(keys(counts)))
        trigger_count[k] = counts[k]
        info = first_info[k]
        trigger_labels[k] = isempty(info) ? string(k) : "$k ($info)"
    end
    return trigger_count, trigger_labels
end

function _add_trigger_legend_entries!(
    ax::Axis,
    trigger_count::OrderedDict{Int,Int},
    trigger_labels::Union{Nothing,OrderedDict{Int,String}} = nothing,
)
    isempty(trigger_count) && return
    for (key, value) in trigger_count
        lbl = (!isnothing(trigger_labels) && haskey(trigger_labels, key)) ? "$(trigger_labels[key]): $value" : "$key: $value"
        scatter!(ax, [-1000.0], [-1000.0], label = lbl, markersize = 0, color = :transparent, alpha = 0)
    end
end

function _render_trigger_overview(
    trigger_times,
    trigger_values,
    trigger_count,
    trigger_labels::Union{Nothing,OrderedDict{Int,String}},
    plot_kwargs::Dict{Symbol,Any},
)
    if isempty(trigger_count)
        @minimal_warning "No triggers found in the data"
        fig = _eegfun_figure()
        ax = Axis(fig[1, 1])
        return (fig = fig, axes = [ax])
    end

    code_keys = collect(keys(trigger_count))
    ytick_labels = isnothing(trigger_labels) ? string.(code_keys) : [trigger_labels[k] for k in code_keys]

    fig = _eegfun_figure()
    ax = Axis(fig[1, 1], yticks = (1:length(code_keys), ytick_labels))
    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    trigger_data = Dict{Int,Vector{Float64}}(Int(k) => Float64[] for k in code_keys)
    @inbounds for (t, v) in zip(trigger_times, trigger_values)
        vec = get(trigger_data, Int(v), nothing)
        !isnothing(vec) && push!(vec, Float64(t))
    end

    ms = plot_kwargs[:marker_size]
    off = plot_kwargs[:line_offset]
    lw = plot_kwargs[:linewidth_trigger]
    for (unique_idx, (key, value)) in enumerate(trigger_count)
        times = trigger_data[Int(key)]
        y_pos = fill(unique_idx, length(times))
        lbl = (!isnothing(trigger_labels) && haskey(trigger_labels, key)) ? "$(trigger_labels[key]): $value" : "$key: $value"
        scatter!(ax, times, y_pos, label = lbl, markersize = ms)
        _plot_trigger_vertical_lines!(ax, times, y_pos; offset = off, linewidth = lw)
    end

    plot_kwargs[:legend] && (fig[1, 2] = Legend(fig, ax))
    ax.ylabel = "Trigger Value"
    ax.xlabel = "Time (s)"
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    if plot_kwargs[:show_title] && !isnothing(plot_kwargs[:plot_title])
        ax.title = plot_kwargs[:plot_title]
    end

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax])
end

"""
    plot_trigger_overview(trigger_times, trigger_values, trigger_count; kwargs...) -> NamedTuple
    plot_trigger_overview(dat::ContinuousData; kwargs...) -> NamedTuple

Plot trigger onset occurrences across time grouped by trigger value with a summary count legend.

$(_generate_kwargs_doc(PLOT_TRIGGERS_KWARGS))
"""
function plot_trigger_overview(trigger_times, trigger_values, trigger_count; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    return _render_trigger_overview(trigger_times, trigger_values, trigger_count, nothing, plot_kwargs)
end

function plot_trigger_overview(dat::ContinuousData; kwargs...)
    @info "Plotting trigger (cleaned) overview for ContinuousData"
    _set_window_title(_generate_window_title(dat))
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    trigger_codes, trigger_times, trigger_info = _extract_trigger_data(dat, plot_kwargs[:ignore_triggers])

    result = if any(!isempty, trigger_info)
        trigger_count, trigger_labels = _count_triggers(trigger_codes, trigger_info)
        _render_trigger_overview(trigger_times, trigger_codes, trigger_count, trigger_labels, plot_kwargs)
    else
        trigger_count = _count_triggers(trigger_codes)
        _render_trigger_overview(trigger_times, trigger_codes, trigger_count, nothing, plot_kwargs)
    end
    _set_window_title("Makie")
    return result
end

function _setup_axis_properties!(ax::Axis, plot_kwargs::Dict = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, (;)))
    ax.xlabel = "Time (s)"
    ax.ylabel = ""
    ax.title = !plot_kwargs[:show_title] ? "" : something(plot_kwargs[:plot_title], "Event Timing and Intervals")
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.titlesize = 22
    hideydecorations!(ax, ticks = true, label = false)
    ylims!(ax, plot_kwargs[:y_min_limit], plot_kwargs[:y_max_limit])
end

function _create_interactive_sliders(fig::Figure, end_time::Float64, plot_kwargs::Dict)
    initial_position = plot_kwargs[:initial_position]
    slider_position =
        Slider(fig[2, 1], range = initial_position:plot_kwargs[:position_step]:end_time, startvalue = initial_position, snap = true)
    slider_size = Slider(
        fig[3, 1],
        range = plot_kwargs[:min_window_size]:plot_kwargs[:window_size_step]:plot_kwargs[:max_window_size],
        startvalue = plot_kwargs[:window_size],
        snap = true,
    )
    Label(fig[2, 2], @lift("Position: $(round($(slider_position.value), digits=1))s"), fontsize = plot_kwargs[:label_fontsize])
    Label(fig[3, 2], @lift("Window Size: $(round($(slider_size.value), digits=1))s"), fontsize = plot_kwargs[:label_fontsize])
    return slider_position, slider_size, initial_position
end

function _plot_trigger_events!(
    ax::Axis,
    trigger_times::AbstractVector{Float64},
    trigger_codes::AbstractVector{Int16},
    trigger_info::Union{Nothing,AbstractVector{String}} = nothing;
    timeline_start::Float64 = trigger_times[1],
    timeline_end::Float64 = trigger_times[end],
    font_size::Real = PLOT_TRIGGERS_KWARGS[:font_size][1],
    event_line_height::Real = PLOT_TRIGGERS_KWARGS[:event_line_height][1],
    event_line_width::Real = PLOT_TRIGGERS_KWARGS[:event_line_width][1],
    time_label_offset::Real = PLOT_TRIGGERS_KWARGS[:time_label_offset][1],
    timeline_width::Real = PLOT_TRIGGERS_KWARGS[:timeline_width][1],
    draw_timeline::Bool = true,
)
    isempty(trigger_times) && (@minimal_error "No triggers to plot"; return)

    line_col = _resolve_theme_linecolor(ax, nothing)
    text_col = _resolve_theme_textcolor(ax, nothing)
    n = length(trigger_times)

    code_strings = Vector{String}(undef, n)
    time_strings = Vector{String}(undef, n)
    xs = Vector{Float64}(undef, 3 * n)
    ys = Vector{Float64}(undef, 3 * n)

    @inbounds for i = 1:n
        c = Int(trigger_codes[i])
        if !isnothing(trigger_info) && !isempty(trigger_info[i])
            code_strings[i] = "$c ($(trigger_info[i]))"
        else
            code_strings[i] = string(c)
        end
        t = trigger_times[i]
        time_strings[i] = string(round(t, digits = 3))

        idx = 3 * i - 2
        xs[idx] = t
        xs[idx+1] = t
        xs[idx+2] = NaN
        ys[idx] = 0.0
        ys[idx+1] = event_line_height
        ys[idx+2] = NaN
    end

    draw_timeline && lines!(ax, [timeline_start, timeline_end], [0.0, 0.0], color = line_col, linewidth = timeline_width)
    lines!(ax, xs, ys, color = line_col, linewidth = event_line_width)
    text!(
        ax,
        trigger_times,
        fill(event_line_height, n);
        text = code_strings,
        align = (:center, :bottom),
        color = text_col,
        fontsize = font_size,
    )
    text!(
        ax,
        trigger_times,
        fill(time_label_offset, n);
        text = time_strings,
        align = (:center, :top),
        color = text_col,
        fontsize = font_size,
    )

    if n > 1
        interval_positions = Vector{Float64}(undef, n - 1)
        interval_strings = Vector{String}(undef, n - 1)
        @inbounds for i = 1:(n-1)
            t1, t2 = trigger_times[i], trigger_times[i+1]
            interval_positions[i] = 0.5 * (t1 + t2)
            interval_strings[i] = string(round(t2 - t1, digits = 3))
        end
        text!(
            ax,
            interval_positions,
            zeros(Float64, n - 1);
            text = interval_strings,
            align = (:center, :top),
            color = text_col,
            fontsize = font_size,
        )
    end
end

function _create_interactive_trigger_plot_impl(
    trigger_codes::Vector{Int16},
    trigger_times::Vector{Float64},
    trigger_info::Union{Nothing,Vector{String}},
    plot_kwargs::Dict{Symbol,Any},
)
    if isempty(trigger_times)
        @minimal_warning "No triggers found in the data"
        fig = _eegfun_figure()
        ax = Axis(fig[1, 1])
        return (fig = fig, axes = [ax])
    end

    end_time = trigger_times[end] + 2.0
    fig = _eegfun_figure()
    ax = Axis(fig[1, 1])
    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    has_info = !isnothing(trigger_info) && any(!isempty, trigger_info)
    if has_info
        trigger_count, trigger_labels = _count_triggers(trigger_codes, trigger_info)
        _add_trigger_legend_entries!(ax, trigger_count, trigger_labels)
    else
        trigger_count = _count_triggers(trigger_codes)
        _add_trigger_legend_entries!(ax, trigger_count)
    end
    plot_kwargs[:legend] && (fig[1, 2] = Legend(fig, ax))

    slider_position, slider_size, _ = _create_interactive_sliders(fig, end_time, plot_kwargs)
    sorted_times = issorted(trigger_times)
    line_col = _resolve_theme_linecolor(ax, nothing)

    function update_plot!()
        empty!(ax)
        current_start = Float64(slider_position.value[])
        current_end = min(current_start + Float64(slider_size.value[]), end_time)

        lines!(ax, [current_start, current_end], [0.0, 0.0], color = line_col, linewidth = plot_kwargs[:timeline_width])

        if sorted_times
            lo = searchsortedfirst(trigger_times, current_start)
            hi = searchsortedlast(trigger_times, current_end)
            if lo <= hi
                rng = lo:hi
                _plot_trigger_events!(
                    ax,
                    @view(trigger_times[rng]),
                    @view(trigger_codes[rng]),
                    has_info ? @view(trigger_info[rng]) : nothing;
                    timeline_start = current_start,
                    timeline_end = current_end,
                    font_size = plot_kwargs[:font_size],
                    event_line_height = plot_kwargs[:event_line_height],
                    event_line_width = plot_kwargs[:event_line_width],
                    time_label_offset = plot_kwargs[:time_label_offset],
                    timeline_width = plot_kwargs[:timeline_width],
                    draw_timeline = false,
                )
            end
        else
            window_mask = (trigger_times .>= current_start) .& (trigger_times .<= current_end)
            if any(window_mask)
                _plot_trigger_events!(
                    ax,
                    trigger_times[window_mask],
                    trigger_codes[window_mask],
                    has_info ? trigger_info[window_mask] : nothing;
                    timeline_start = current_start,
                    timeline_end = current_end,
                    font_size = plot_kwargs[:font_size],
                    event_line_height = plot_kwargs[:event_line_height],
                    event_line_width = plot_kwargs[:event_line_width],
                    time_label_offset = plot_kwargs[:time_label_offset],
                    timeline_width = plot_kwargs[:timeline_width],
                    draw_timeline = false,
                )
            end
        end

        xlims!(ax, current_start, current_end)
    end

    onany(slider_position.value, slider_size.value) do _, _
        update_plot!()
    end

    _setup_axis_properties!(ax, plot_kwargs)
    update_plot!()
    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax])
end

function _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    return _create_interactive_trigger_plot_impl(trigger_codes, trigger_times, nothing, plot_kwargs)
end

function _create_interactive_trigger_plot(
    trigger_codes::Vector{Int16},
    trigger_times::Vector{Float64},
    trigger_info::Vector{String};
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    return _create_interactive_trigger_plot_impl(trigger_codes, trigger_times, trigger_info, plot_kwargs)
end

"""
    plot_trigger_timing(dat::ContinuousData; kwargs...) -> NamedTuple

Plot interactive trigger onset timing and inter-trigger intervals with position and window-size sliders.

$(_generate_kwargs_doc(PLOT_TRIGGERS_KWARGS))
"""
function plot_trigger_timing(dat::ContinuousData; kwargs...)
    _set_window_title(_generate_window_title(dat))
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    trigger_codes, trigger_times, trigger_info = _extract_trigger_data(dat, plot_kwargs[:ignore_triggers])
    result = _create_interactive_trigger_plot_impl(trigger_codes, trigger_times, trigger_info, plot_kwargs)
    _set_window_title("Makie")
    return result
end
