########################################################
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_TRIGGERS_KWARGS = Dict(
    :window_size => 10.0,
    :initial_position => -2.0,
    :min_window_size => 0.1,
    :max_window_size => 100.0,
    :label_fontsize => 20,
    :display_plot => true,
    :ignore_triggers => Int[],  # Vector of trigger codes to ignore
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
function _trigger_time_count(time, triggers, ignore_triggers=Int[])

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
    _plot_trigger_vertical_lines!(ax::Axis, times::Vector{Float64}, y_pos::Vector{Int})

Helper function to plot vertical lines for trigger events.
"""
function _plot_trigger_vertical_lines!(ax::Axis, times::Vector{Float64}, y_positions::Vector{Int})
    for (xpos, ypos) in zip(times, y_positions)
        lines!(
            ax,
            [xpos, xpos],
            [ypos - DEFAULT_LINE_OFFSET, ypos + DEFAULT_LINE_OFFSET],
            color = :black,
            linewidth = DEFAULT_LINE_WIDTH_TRIGGER,
        )
    end
end

"""
    _plot_single_trigger_line!(ax::Axis, time::Float64)

Helper function to plot a single trigger vertical line.
"""
function _plot_single_trigger_line!(ax::Axis, time::Float64)
    lines!(ax, [time, time], [0, EVENT_LINE_HEIGHT], color = :black, linewidth = DEFAULT_EVENT_LINE_WIDTH)
end

"""
    _add_trigger_text!(ax::Axis, x::Float64, y::Float64, text_str::String, align::Tuple)

Helper function to add trigger text with consistent styling.
"""
function _add_trigger_text!(ax::Axis, x::Float64, y::Float64, text_str::String, align::Tuple)
    text!(ax, x, y, text = text_str, align = align, color = :black, fontsize = DEFAULT_FONT_SIZE)
end

"""
    _extract_trigger_data(dat::BiosemiDataFormat.BiosemiData, ignore_triggers=Int[])

Extract trigger information from BioSemi data.

# Arguments
- `dat`: BioSemiData object
- `ignore_triggers`: Vector of trigger codes to ignore (optional)

# Returns
- `trigger_codes`: Vector of trigger codes
- `trigger_times`: Vector of trigger times
"""
function _extract_trigger_data(dat::BiosemiDataFormat.BiosemiData, ignore_triggers=Int[])
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    trigger_positions = findall(x -> x != 0, cleaned_triggers)
    trigger_codes = Int16.(cleaned_triggers[trigger_positions])
    trigger_times = dat.time[trigger_positions]
    
    # Filter out ignored triggers if any are specified
    if !isempty(ignore_triggers)
        keep_mask = .!in.(trigger_codes, Ref(ignore_triggers))
        trigger_codes = trigger_codes[keep_mask]
        trigger_times = trigger_times[keep_mask]
    end
    
    return trigger_codes, trigger_times
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
"""
function _extract_trigger_data(dat::ContinuousData, ignore_triggers=Int[])
    trigger_col = :triggers
    if !hasproperty(dat.data, trigger_col)
        @minimal_error("No triggers column found in data. Expected column name: $trigger_col")
    end

    trigger_positions = findall(x -> x != 0, dat.data[!, trigger_col])
    trigger_codes = Int16.(dat.data[trigger_positions, trigger_col])
    trigger_times = dat.data[trigger_positions, :time]
    
    # Filter out ignored triggers if any are specified
    if !isempty(ignore_triggers)
        keep_mask = .!in.(trigger_codes, Ref(ignore_triggers))
        trigger_codes = trigger_codes[keep_mask]
        trigger_times = trigger_times[keep_mask]
    end
    
    return trigger_codes, trigger_times
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
    plot_kwargs = merge(DEFAULT_TRIGGER_KWARGS, Dict(kwargs))

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
        scatter!(ax, times, y_pos, label = "$key: $(string(value))", markersize = DEFAULT_MARKER_SIZE)
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
    plot_trigger_overview(dat::BiosemiDataFormat.BiosemiData; kwargs...)

Plot trigger events from BioSemi BDF data.

# Arguments
- `dat`: BioSemiData object containing the EEG data
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
function plot_trigger_overview(dat::BiosemiDataFormat.BiosemiData; kwargs...)
    @info "Plotting trigger (raw) overview for BioSemi data"
    # Merge user kwargs with defaults
    plot_kwargs = merge(DEFAULT_TRIGGER_KWARGS, Dict(kwargs))
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.time, dat.triggers.raw, plot_kwargs[:ignore_triggers])
    return plot_trigger_overview(trigger_times, trigger_values, trigger_count; kwargs...)
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
    # Merge user kwargs with defaults
    plot_kwargs = merge(DEFAULT_TRIGGER_KWARGS, Dict(kwargs))
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.data.time, dat.data.triggers, plot_kwargs[:ignore_triggers])
    return plot_trigger_overview(trigger_times, trigger_values, trigger_count; kwargs...)
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
    ylims!(ax, Y_MIN_LIMIT, Y_MAX_LIMIT)
end

"""
    _create_interactive_sliders(fig::Figure, end_time::Float64, plot_kwargs::Dict)

Create interactive sliders for position and window size.
"""
function _create_interactive_sliders(fig::Figure, end_time::Float64, plot_kwargs::Dict)
    initial_position = DEFAULT_TRIGGER_KWARGS[:initial_position] 

    slider_position =
        Slider(fig[2, 1], range = initial_position:POSITION_STEP:end_time, startvalue = initial_position, snap = true)

    slider_size = Slider(
        fig[3, 1],
        range = plot_kwargs[:min_window_size]:WINDOW_SIZE_STEP:plot_kwargs[:max_window_size],
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
function _plot_trigger_events!(
    ax::Axis,
    trigger_times::Vector{Float64},
    trigger_codes::Vector{Int16}
)
    # Early return if no triggers to plot
    if isempty(trigger_times)
        @minimal_error "No triggers to plot"
        return
    end

    # Pre-compute string conversions
    code_strings = [string(Int(code)) for code in trigger_codes]
    time_strings = [string(round(time, digits = 2)) for time in trigger_times]

    # Pre-compute intervals and their string representations
    intervals = diff(trigger_times)
    interval_strings = [string(round(interval, digits = 2)) for interval in intervals]
    interval_positions = [(trigger_times[i] + trigger_times[i+1]) / 2 for i = 1:length(intervals)]

    # Plot horizontal timeline
    timeline_start = 0.0
    lines!(ax, [timeline_start, trigger_times[end]], [0, 0], color = :black, linewidth = DEFAULT_TIMELINE_WIDTH)

    # Plot vertical lines for each event
    for (time, code_str, time_str) in zip(trigger_times, code_strings, time_strings)
        _plot_single_trigger_line!(ax, time)

        # Add trigger code at the top of the line
        _add_trigger_text!(ax, time, EVENT_LINE_HEIGHT, code_str, (:center, :bottom))
        
        # Add time value below the line
        _add_trigger_text!(ax, time, TIME_LABEL_OFFSET, time_str, (:center, :top))
    end

    # Plot intervals as text
    for (x, interval_str) in zip(interval_positions, interval_strings)
        text!(ax, x, 0, text = interval_str, align = (:center, :top), color = :black, fontsize = DEFAULT_FONT_SIZE)
    end
end

"""
    _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)

Internal function to create interactive trigger timing plot.
"""
function _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = merge(DEFAULT_TRIGGER_KWARGS, Dict(kwargs))
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

    # Add trigger count legend entries
    trigger_count = _count_triggers(trigger_codes)
    _add_trigger_legend_entries!(ax, trigger_count)
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

        if !isempty(window_times)
            _plot_trigger_events!(ax, window_times, window_codes)
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
    plot_trigger_timing(dat::BiosemiDataFormat.BiosemiData; kwargs...)

Plot trigger timing with interactive x-axis sliders for scrolling and window size.

# Arguments
- `dat::BiosemiDataFormat.BiosemiData`: The BioSemiData object containing the triggers
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
function plot_trigger_timing(dat::BiosemiDataFormat.BiosemiData; kwargs...)
    plot_kwargs = merge(DEFAULT_TRIGGER_KWARGS, Dict(kwargs))
    trigger_codes, trigger_times = _extract_trigger_data(dat, plot_kwargs[:ignore_triggers])
    return _create_interactive_trigger_plot(trigger_codes, trigger_times; kwargs...)
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
    plot_kwargs = merge(DEFAULT_TRIGGER_KWARGS, Dict(kwargs))
    trigger_codes, trigger_times = _extract_trigger_data(dat, plot_kwargs[:ignore_triggers])
    return _create_interactive_trigger_plot(trigger_codes, trigger_times; kwargs...)
end
