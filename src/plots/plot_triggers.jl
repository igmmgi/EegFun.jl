########################################################
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_TRIGGERS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :window_size => (10.0, "Size of the time window to display in seconds."),
    :initial_position => (-2.0, "Initial position of the time window in seconds."),
    :min_window_size => (0.1, "Minimum size of the time window in seconds."),
    :max_window_size => (100.0, "Maximum size of the time window in seconds."),
    :label_fontsize => (20, "Font size for labels."),
    :display_plot => (true, "Whether to display the plot."),
    :ignore_triggers => (Int[], "Vector of trigger codes to ignore."),
    :marker_size => (15, "Size of markers in the trigger overview plot."),
    :line_width_trigger => (1, "Line width for trigger vertical lines."),
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
)

########################################################
# Shared trigger utilities
########################################################

"""
    _filter_triggers(trigger_times, trigger_values, ignore_triggers)

Filter out specified trigger codes from trigger data.
Note: This function assumes ignore_triggers is not empty.

# Arguments
- `trigger_times`: Vector of trigger times
- `trigger_values`: Vector of trigger values
- `ignore_triggers`: Vector of trigger codes to ignore (assumed non-empty)

# Returns
- `filtered_times`: Filtered trigger times
- `filtered_values`: Filtered trigger values
"""
function _filter_triggers(trigger_times, trigger_values, ignore_triggers)
    keep_mask = .!in.(trigger_values, Ref(ignore_triggers))
    return trigger_times[keep_mask], trigger_values[keep_mask]
end

"""
    _trigger_time_count(time, triggers, ignore_triggers=Int[])

Internal function to process trigger data and count occurrences.

# Arguments
- `time`: Vector of time points
- `triggers`: Vector of trigger values
- `ignore_triggers`: Vector of trigger codes to ignore (optional)

# Returns
- `trigger_times`: Vector of times when triggers occurred
- `trigger_values`: Vector of trigger values at those times
- `trigger_count`: OrderedDict mapping trigger values to their counts
"""
function _trigger_time_count(time, triggers, ignore_triggers = Int[])

    # Since triggers are already cleaned (only onset values), just find non-zero values
    trigger_indices = findall(triggers .!= 0)
    if isempty(trigger_indices)
        return Float64[], Int[], OrderedDict{Int,Int}()
    end

    trigger_values = triggers[trigger_indices]
    trigger_times = time[trigger_indices]

    # Filter out ignored triggers if any are specified
    if !isempty(ignore_triggers)
        trigger_times, trigger_values = _filter_triggers(trigger_times, trigger_values, ignore_triggers)
        if isempty(trigger_values)
            return Float64[], Int[], OrderedDict{Int,Int}()
        end
    end

    trigger_count = OrderedDict(i => 0 for i in sort!(collect(Set(trigger_values))))
    for val in trigger_values
        trigger_count[val] += 1
    end

    return trigger_times, trigger_values, trigger_count
end

"""
    _trigger_time_count(times::Vector{Float64}, trigger_codes::Vector{Int16}, trigger_info::Vector{String}, ignore_triggers=Int[])

Count trigger occurrences with trigger info support.

# Arguments
- `times`: Vector of trigger times
- `trigger_codes`: Vector of trigger codes
- `trigger_info`: Vector of trigger info strings
- `ignore_triggers`: Vector of trigger codes to ignore (optional)

# Returns
- `trigger_times`: Vector of trigger times
- `trigger_values`: Vector of trigger values
- `trigger_count`: OrderedDict mapping trigger codes to counts
"""
function _trigger_time_count(
    times::Vector{Float64},
    trigger_codes::Vector{Int16},
    trigger_info::Vector{String},
    ignore_triggers = Int[],
)
    # Filter out ignored triggers if any are specified
    if !isempty(ignore_triggers)
        times, trigger_codes, trigger_info = _filter_triggers(times, trigger_codes, trigger_info, ignore_triggers)
    end

    # Count occurrences of each trigger code
    trigger_count = OrderedDict{Int,Int}()
    for val in trigger_codes
        trigger_count[Int(val)] = get(trigger_count, Int(val), 0) + 1
    end

    return times, Int.(trigger_codes), trigger_count
end

"""
    _filter_triggers(trigger_times, trigger_codes, trigger_info, ignore_triggers)

Filter out specified trigger codes from trigger data with trigger info.
"""
function _filter_triggers(trigger_times, trigger_codes, trigger_info, ignore_triggers)
    keep_mask = .!in.(trigger_codes, Ref(ignore_triggers))
    return trigger_times[keep_mask], trigger_codes[keep_mask], trigger_info[keep_mask]
end

"""
    _plot_trigger_vertical_lines!(ax::Axis, times::Vector{Float64}, y_pos::Vector{Int})

Helper function to plot vertical lines for trigger events.
"""
function _plot_trigger_vertical_lines!(ax::Axis, times::Vector{Float64}, y_positions::Vector{Int})
    for (xpos, ypos) in zip(times, y_positions)
        lines!(
            ax,
            [xpos, xpos],
            [ypos - PLOT_TRIGGERS_KWARGS[:line_offset][1], ypos + PLOT_TRIGGERS_KWARGS[:line_offset][1]],
            color = :black,
            linewidth = PLOT_TRIGGERS_KWARGS[:line_width_trigger][1],
        )
    end
end

"""
    _plot_single_trigger_line!(ax::Axis, time::Float64)

Helper function to plot a single trigger vertical line.
"""
function _plot_single_trigger_line!(ax::Axis, time::Float64)
    lines!(
        ax,
        [time, time],
        [0, PLOT_TRIGGERS_KWARGS[:event_line_height][1]],
        color = :black,
        linewidth = PLOT_TRIGGERS_KWARGS[:event_line_width][1],
    )
end

"""
    _add_trigger_text!(ax::Axis, x::Float64, y::Float64, text_str::String, align::Tuple)

Helper function to add trigger text with consistent styling.
"""
function _add_trigger_text!(ax::Axis, x::Float64, y::Float64, text_str::String, align::Tuple)
    text!(ax, x, y, text = text_str, align = align, color = :black, fontsize = PLOT_TRIGGERS_KWARGS[:font_size][1])
end


"""
    _extract_trigger_data(dat::ContinuousData, ignore_triggers=Int[])

Extract trigger information from ContinuousData.

# Arguments
- `dat`: ContinuousData object
- `ignore_triggers`: Vector of trigger codes to ignore (optional)

# Returns
- `trigger_codes`: Vector of trigger codes
- `trigger_times`: Vector of trigger times
- `trigger_info`: Vector of trigger info strings (empty strings if not available)
"""
function _extract_trigger_data(dat::ContinuousData, ignore_triggers = Int[])
    trigger_col = :triggers
    if !hasproperty(dat.data, trigger_col)
        @minimal_error("No triggers column found in data. Expected column name: $trigger_col")
    end

    trigger_positions = findall(x -> x != 0, dat.data[!, trigger_col])
    trigger_codes = Int16.(dat.data[trigger_positions, trigger_col])
    trigger_times = dat.data[trigger_positions, :time]

    # Extract trigger info if available
    trigger_info =
        hasproperty(dat.data, :triggers_info) ? dat.data[trigger_positions, :triggers_info] :
        fill("", length(trigger_positions))

    # Filter out ignored triggers if any are specified
    if !isempty(ignore_triggers)
        keep_mask = .!in.(trigger_codes, Ref(ignore_triggers))
        trigger_codes = trigger_codes[keep_mask]
        trigger_times = trigger_times[keep_mask]
        trigger_info = trigger_info[keep_mask]
    end

    return trigger_codes, trigger_times, trigger_info
end



"""
    _count_triggers(trigger_codes::Vector{Int16})

Count occurrences of each trigger code.
"""
function _count_triggers(trigger_codes::Vector{Int16})
    trigger_count = OrderedDict{Int,Int}()
    for code in sort!(collect(Set(trigger_codes)))
        trigger_count[code] = count(x -> x == code, trigger_codes)
    end
    return trigger_count
end

"""
    _count_triggers(trigger_codes::Vector{Int16}, trigger_info::Vector{String})

Count occurrences of each trigger code with info labels.

# Returns
- `trigger_count`: OrderedDict mapping trigger codes to counts
- `trigger_labels`: OrderedDict mapping trigger codes to display labels (code + info if available)
"""
function _count_triggers(trigger_codes::Vector{Int16}, trigger_info::Vector{String})
    trigger_count = OrderedDict{Int,Int}()
    trigger_labels = OrderedDict{Int,String}()

    for code in sort!(collect(Set(trigger_codes)))
        trigger_count[code] = count(x -> x == code, trigger_codes)

        # Find the first occurrence of this trigger code to get its info
        first_idx = findfirst(x -> x == code, trigger_codes)
        if !isnothing(first_idx) && first_idx <= length(trigger_info)
            info = trigger_info[first_idx]
            if !isempty(info)
                trigger_labels[code] = "$code ($info)"
            else
                trigger_labels[code] = string(code)
            end
        else
            trigger_labels[code] = string(code)
        end
    end

    return trigger_count, trigger_labels
end


"""
    _add_trigger_legend_entries!(ax::Axis, trigger_count::OrderedDict{Int,Int})

Add invisible scatter points with labels for legend creation.
"""
function _add_trigger_legend_entries!(ax::Axis, trigger_count::OrderedDict{Int,Int})
    if isempty(trigger_count)
        return
    end
    # Add invisible scatter points with labels for each trigger type
    for (key, value) in trigger_count
        scatter!(ax, [-1000], [-1000], label = "$key: $value", markersize = 0, color = :transparent, alpha = 0)
    end
end

"""
    _add_trigger_legend_entries!(ax::Axis, trigger_count::OrderedDict{Int,Int}, trigger_labels::OrderedDict{Int,String})

Add invisible scatter points with labels for legend creation using trigger info.
"""
function _add_trigger_legend_entries!(
    ax::Axis,
    trigger_count::OrderedDict{Int,Int},
    trigger_labels::OrderedDict{Int,String},
)
    if isempty(trigger_count)
        return
    end
    # Add invisible scatter points with labels for each trigger type using trigger info
    for (key, value) in trigger_count
        label = haskey(trigger_labels, key) ? "$(trigger_labels[key]): $value" : "$key: $value"
        scatter!(ax, [-1000], [-1000], label = label, markersize = 0, color = :transparent, alpha = 0)
    end
end

########################################################
# Trigger overview plotting functions
########################################################

"""
    plot_trigger_overview(trigger_times, trigger_values, trigger_count)

Plot trigger events as a scatter plot with vertical lines.

# Arguments
- `trigger_times`: Vector of times when triggers occurred
- `trigger_values`: Vector of trigger values at those times
- `trigger_count`: OrderedDict mapping trigger values to their counts

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_trigger_overview(trigger_times, trigger_values, trigger_count; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)

    if isempty(trigger_count)
        @minimal_warning "No triggers found in the data"
        fig = Figure()
        ax = Axis(fig[1, 1])
        return fig, ax
    end

    fig = Figure()
    ax = Axis(fig[1, 1], yticks = (1:length(trigger_count.keys), string.(trigger_count.keys)))

    # Pre-compute trigger data for each type to avoid repeated filtering
    trigger_data = Dict{Int,Vector{Float64}}()
    for (key, _) in trigger_count
        trigger_data[key] = trigger_times[trigger_values .== key]
    end

    for (unique, (key, value)) in enumerate(trigger_count)
        times = trigger_data[key]
        y_pos = fill(unique, length(times))
        scatter!(ax, times, y_pos, label = "$key: $(string(value))", markersize = PLOT_TRIGGERS_KWARGS[:marker_size][1])
        _plot_trigger_vertical_lines!(ax, times, y_pos)
    end
    fig[1, 2] = Legend(fig, ax)
    ax.ylabel = "Trigger Value"
    ax.xlabel = "Time (S)"

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    return fig, ax

end



"""
    plot_trigger_overview(dat::ContinuousData; kwargs...)

Plot trigger events from ContinuousData object.

# Arguments
- `dat`: ContinuousData object containing the EEG data
- `ignore_triggers`: Vector of trigger codes to ignore (optional)
- Other plotting parameters (window_size, display_plot, etc.)

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot

# Example
```julia
# Plot all triggers
fig, ax = plot_trigger_overview(dat)

# Ignore specific trigger codes
fig, ax = plot_trigger_overview(dat; ignore_triggers=[1, 255])
```
"""
function plot_trigger_overview(dat::ContinuousData; kwargs...)
    @info "Plotting trigger (cleaned) overview for ContinuousData"

    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    trigger_codes, trigger_times, trigger_info = _extract_trigger_data(dat, plot_kwargs[:ignore_triggers])
    has_info = any(!isempty, trigger_info)

    if has_info
        # Use enhanced plotting with trigger info
        trigger_count, trigger_labels = _count_triggers(trigger_codes, trigger_info)

        # Create the plot manually to avoid duplicate legends
        fig = Figure()
        ax = Axis(fig[1, 1], yticks = (1:length(trigger_count.keys), [trigger_labels[k] for k in trigger_count.keys]))

        # Pre-compute trigger data for each type
        trigger_data = Dict{Int,Vector{Float64}}()
        for (key, _) in trigger_count
            trigger_data[key] = trigger_times[Int.(trigger_codes) .== key]
        end

        for (unique, (key, value)) in enumerate(trigger_count)
            times = trigger_data[key]
            y_pos = fill(unique, length(times))
            scatter!(
                ax,
                times,
                y_pos,
                label = "$(trigger_labels[key]): $(string(value))",
                markersize = PLOT_TRIGGERS_KWARGS[:marker_size][1],
            )
            _plot_trigger_vertical_lines!(ax, times, y_pos)
        end

        fig[1, 2] = Legend(fig, ax)
        ax.ylabel = "Trigger Value"
        ax.xlabel = "Time (S)"

        if plot_kwargs[:display_plot]
            display_figure(fig)
        end

        return fig, ax
    else
        # Fall back to original behavior
        trigger_count = _count_triggers(trigger_codes)
        fig, ax = plot_trigger_overview(trigger_times, Int.(trigger_codes), trigger_count; kwargs...)
        return fig, ax
    end
end

########################################################
# Trigger timing plotting functions
########################################################

"""
    _setup_axis_properties!(ax::Axis)

Set up common axis properties for trigger timing plots.
"""
function _setup_axis_properties!(ax::Axis)
    ax.xlabel = "Time (s)"
    ax.ylabel = ""
    ax.title = "Event Timing and Intervals"
    ax.xlabelsize = 20
    ax.titlesize = 22
    hideydecorations!(ax, ticks = true, label = false)
    ylims!(ax, PLOT_TRIGGERS_KWARGS[:y_min_limit][1], PLOT_TRIGGERS_KWARGS[:y_max_limit][1])
end

"""
    _create_interactive_sliders(fig::Figure, end_time::Float64, plot_kwargs::Dict)

Create interactive sliders for position and window size.
"""
function _create_interactive_sliders(fig::Figure, end_time::Float64, plot_kwargs::Dict)
    initial_position = PLOT_TRIGGERS_KWARGS[:initial_position][1]

    slider_position = Slider(
        fig[2, 1],
        range = initial_position:PLOT_TRIGGERS_KWARGS[:position_step][1]:end_time,
        startvalue = initial_position,
        snap = true,
    )

    slider_size = Slider(
        fig[3, 1],
        range = plot_kwargs[:min_window_size]:PLOT_TRIGGERS_KWARGS[:window_size_step][1]:plot_kwargs[:max_window_size],
        startvalue = plot_kwargs[:window_size],
        snap = true,
    )

    # Create labels for sliders
    position_label = Label(fig[2, 2], @lift("Position: $(round($(slider_position.value), digits=1))s"), fontsize = plot_kwargs[:label_fontsize])
    size_label     = Label(fig[3, 2], @lift("Window Size: $(round($(slider_size.value), digits=1))s"), fontsize = plot_kwargs[:label_fontsize])

    return slider_position, slider_size, initial_position
end

"""
    _plot_trigger_events!(ax::Axis, trigger_times::Vector{Float64}, trigger_codes::Vector{Int16})

Core function to plot trigger events on the given axis.
"""
function _plot_trigger_events!(ax::Axis, trigger_times::Vector{Float64}, trigger_codes::Vector{Int16})
    _plot_trigger_events!(ax, trigger_times, trigger_codes, fill("", length(trigger_codes)))
end

function _plot_trigger_events!(
    ax::Axis,
    trigger_times::Vector{Float64},
    trigger_codes::Vector{Int16},
    trigger_info::Vector{String},
)
    # Early return if no triggers to plot
    if isempty(trigger_times)
        @minimal_error "No triggers to plot"
        return
    end

    # Pre-compute string conversions with optional trigger info
    code_strings = String[]
    for (code, info) in zip(trigger_codes, trigger_info)
        if !isempty(info)
            push!(code_strings, "$(Int(code)) ($info)")
        else
            push!(code_strings, string(Int(code)))
        end
    end

    time_strings = [string(round(time, digits = 2)) for time in trigger_times]

    # Pre-compute intervals and their string representations
    intervals = diff(trigger_times)
    interval_strings = [string(round(interval, digits = 2)) for interval in intervals]
    interval_positions = [(trigger_times[i] + trigger_times[i+1]) / 2 for i = 1:length(intervals)]

    # Plot horizontal timeline
    timeline_start = 0.0
    lines!(
        ax,
        [timeline_start, trigger_times[end]],
        [0, 0],
        color = :black,
        linewidth = PLOT_TRIGGERS_KWARGS[:timeline_width][1],
    )

    # Plot vertical lines for each event
    for (time, code_str, time_str) in zip(trigger_times, code_strings, time_strings)
        _plot_single_trigger_line!(ax, time)

        # Add trigger code (with optional info) at the top of the line
        _add_trigger_text!(ax, time, PLOT_TRIGGERS_KWARGS[:event_line_height][1], code_str, (:center, :bottom))

        # Add time value below the line
        _add_trigger_text!(ax, time, PLOT_TRIGGERS_KWARGS[:time_label_offset][1], time_str, (:center, :top))
    end

    # Plot intervals as text
    for (x, interval_str) in zip(interval_positions, interval_strings)
        text!(
            ax,
            x,
            0,
            text = interval_str,
            align = (:center, :top),
            color = :black,
            fontsize = PLOT_TRIGGERS_KWARGS[:font_size][1],
        )
    end
end


"""
    _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)

Internal function to create interactive trigger timing plot.
"""
function _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)
    _create_interactive_trigger_plot(trigger_codes, trigger_times, fill("", length(trigger_codes)); kwargs...)
end

function _create_interactive_trigger_plot(
    trigger_codes::Vector{Int16},
    trigger_times::Vector{Float64},
    trigger_info::Vector{String};
    kwargs...,
)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, kwargs)
    if isempty(trigger_times)
        @minimal_warning "No triggers found in the data"
        fig = Figure()
        ax = Axis(fig[1, 1])
        return fig, ax
    end

    # Calculate time range
    end_time = trigger_times[end] + 2.0

    # Create figure with grid layout to accommodate legend
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Add trigger count legend entries (with or without trigger info)
    has_info = any(!isempty, trigger_info)
    if has_info
        trigger_count, trigger_labels = _count_triggers(trigger_codes, trigger_info)
        _add_trigger_legend_entries!(ax, trigger_count, trigger_labels)
    else
        trigger_count = _count_triggers(trigger_codes)
        _add_trigger_legend_entries!(ax, trigger_count)
    end
    fig[1, 2] = Legend(fig, ax)

    # Create Observables for reactive plotting
    window_size = Observable(plot_kwargs[:window_size])
    slider_position, slider_size, initial_position = _create_interactive_sliders(fig, end_time, plot_kwargs)
    window_position = Observable(plot_kwargs[:initial_position])

    # Update Observables when sliders change
    onany(slider_position.value, slider_size.value) do pos, size
        window_position[] = pos
        window_size[] = size
    end

    # Reactive plotting function
    function update_plot!()
        empty!(ax)

        # Get current window bounds
        current_start = window_position[]
        current_end = min(current_start + window_size[], end_time)

        # Filter triggers within the current window
        window_mask = (trigger_times .>= current_start) .&& (trigger_times .<= current_end)
        window_times = trigger_times[window_mask]
        window_codes = trigger_codes[window_mask]
        window_info = trigger_info[window_mask]

        if !isempty(window_times)
            if has_info
                _plot_trigger_events!(ax, window_times, window_codes, window_info)
            else
                _plot_trigger_events!(ax, window_times, window_codes)
            end
        end

        # Update only x-axis limits
        xlims!(ax, current_start, current_end)
    end

    # Set up reactive updates
    onany(window_position, window_size) do _, _
        update_plot!()
    end

    # Set axis properties once
    _setup_axis_properties!(ax)

    # Initial plot
    update_plot!()

    # Handle display
    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    return fig, ax

end





"""
    plot_trigger_timing(dat::ContinuousData; kwargs...)

Plot trigger timing with interactive x-axis sliders for scrolling and window size.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing the triggers
- `ignore_triggers`: Vector of trigger codes to ignore (optional)
- Other plotting parameters (window_size, display_plot, etc.)

# Returns
- `fig::Figure`: The Makie figure object
- `ax::Axis`: The Makie axis object

# Example
```julia
# Plot all triggers
fig, ax = plot_trigger_timing(dat)

# Ignore specific trigger codes
fig, ax = plot_trigger_timing(dat; ignore_triggers=[1, 255])
```
"""
function plot_trigger_timing(dat::ContinuousData; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TRIGGERS_KWARGS, Dict(kwargs))
    trigger_codes, trigger_times, trigger_info = _extract_trigger_data(dat, plot_kwargs[:ignore_triggers])
    return _create_interactive_trigger_plot(trigger_codes, trigger_times, trigger_info; kwargs...)
end
