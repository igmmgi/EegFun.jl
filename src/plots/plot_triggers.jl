########################################################
# Trigger plotting functions
########################################################

# Pre-allocated vectors for performance
const VERTICAL_LINE_X = Ref([0.0, 0.0])
const VERTICAL_LINE_Y = Ref([0.0, EVENT_LINE_HEIGHT])

########################################################
# Shared trigger utilities
########################################################

"""
    _trigger_time_count(time, triggers)

Internal function to process trigger data and count occurrences.

# Arguments
- `time`: Vector of time points
- `triggers`: Vector of trigger values

# Returns
- `trigger_times`: Vector of times when triggers occurred
- `trigger_values`: Vector of trigger values at those times
- `trigger_count`: OrderedDict mapping trigger values to their counts
"""
function _trigger_time_count(time, triggers)
    # Since triggers are already cleaned (only onset values), just find non-zero values
    trigger_indices = findall(triggers .!= 0)

    if isempty(trigger_indices)
        return Float64[], Int[], OrderedDict{Int,Int}()
    end

    trigger_values = triggers[trigger_indices]
    trigger_times = time[trigger_indices]
    trigger_count = OrderedDict(i => 0 for i in sort!(collect(Set(trigger_values))))
    for val in trigger_values
        trigger_count[val] += 1
    end

    return trigger_times, trigger_values, trigger_count
end

"""
    _extract_trigger_data(dat::BioSemiBDF.BioSemiData)

Extract trigger information from BioSemi data.
"""
function _extract_trigger_data(dat::BioSemiBDF.BioSemiData)
    raw_triggers = dat.triggers.raw
    cleaned_triggers = _clean_triggers(raw_triggers)
    trigger_positions = findall(x -> x != 0, cleaned_triggers)
    trigger_codes = Int16.(cleaned_triggers[trigger_positions])
    trigger_times = dat.time[trigger_positions]
    return trigger_codes, trigger_times
end

"""
    _extract_trigger_data(dat::ContinuousData)

Extract trigger information from ContinuousData.
"""
function _extract_trigger_data(dat::ContinuousData)
    trigger_col = :triggers
    if !hasproperty(dat.data, trigger_col)
        error("No triggers column found in data. Expected column name: $trigger_col")
    end
    
    trigger_positions = findall(x -> x != 0, dat.data[!, trigger_col])
    trigger_codes = Int16.(dat.data[trigger_positions, trigger_col])
    trigger_times = dat.data[trigger_positions, :time]
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
        # Place invisible points far outside the plot area
        scatter!(ax, [-1000], [-1000], label = "$key: $value", 
                markersize = 0, color = :transparent, alpha = 0)
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
function plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = true)

    if isempty(trigger_count)
        @warn "No triggers found in the data"
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
        # Add vertical lines
        for (t, y) in zip(times, y_pos)
            lines!(
                ax,
                [t, t],
                [y - DEFAULT_LINE_OFFSET, y + DEFAULT_LINE_OFFSET],
                color = :black,
                linewidth = DEFAULT_LINE_WIDTH_TRIGGER,
            )
        end
    end
    fig[1, 2] = Legend(fig, ax)
    ax.ylabel = "Trigger Value"
    ax.xlabel = "Time (S)"

    if display_plot
        display(fig)
    end

    return fig, ax

end

"""
    plot_trigger_overview(dat::BioSemiBDF.BioSemiData)

Plot trigger events from BioSemi BDF data.

# Arguments
- `dat`: BioSemiData object containing the EEG data

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_trigger_overview(dat::BioSemiBDF.BioSemiData; display_plot = true)
    @info "Plotting trigger (raw) overview for BioSemi data"
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.time, dat.triggers.raw)
    return plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = display_plot)
end

"""
    plot_trigger_overview(dat::ContinuousData)

Plot trigger events from ContinuousData object.

# Arguments
- `dat`: ContinuousData object containing the EEG data

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_trigger_overview(dat::ContinuousData; display_plot = true)
    @info "Plotting trigger (cleaned) overview for ContinuousData"
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.data.time, dat.data.triggers)
    return plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = display_plot)
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
    _create_interactive_sliders(fig::Figure, end_time::Float64)

Create interactive sliders for position and window size.
"""
function _create_interactive_sliders(fig::Figure, end_time::Float64)
    initial_position = -2.0
    
    slider_position = Slider(
        fig[2, 1], 
        range = initial_position:POSITION_STEP:end_time, 
        startvalue = initial_position, 
        snap = true
    )
    
    slider_size = Slider(
        fig[3, 1],
        range = MIN_WINDOW_SIZE:WINDOW_SIZE_STEP:MAX_WINDOW_SIZE,
        startvalue = DEFAULT_WINDOW_SIZE,
        snap = true,
    )
    
    # Create labels for sliders
    position_label = Label(fig[2, 2], @lift("Position: $(round($(slider_position.value), digits=1))s"), fontsize = 20)
    size_label = Label(fig[3, 2], @lift("Window Size: $(round($(slider_size.value), digits=1))s"), fontsize = 20)
    
    return slider_position, slider_size, initial_position
end

"""
    _plot_trigger_events!(ax::Axis, trigger_times::Vector{Float64}, trigger_codes::Vector{Int16}; 
                         use_preallocated::Bool=false)

Core function to plot trigger events on the given axis.
"""
function _plot_trigger_events!(ax::Axis, trigger_times::Vector{Float64}, trigger_codes::Vector{Int16}; 
                              use_preallocated::Bool=false)
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
        if use_preallocated
            # Use pre-allocated vectors for performance
            VERTICAL_LINE_X[][1] = VERTICAL_LINE_X[][2] = time
            lines!(ax, VERTICAL_LINE_X[], VERTICAL_LINE_Y[], color = :black, linewidth = DEFAULT_EVENT_LINE_WIDTH)
        else
            # Create new vectors (for non-interactive plotting)
            lines!(ax, [time, time], [0, EVENT_LINE_HEIGHT], color = :black, linewidth = DEFAULT_EVENT_LINE_WIDTH)
        end
        
        # Add trigger code at the top of the line
        text!(
            ax,
            time,
            EVENT_LINE_HEIGHT,
            text = code_str,
            align = (:center, :bottom),
            color = :black,
            fontsize = DEFAULT_FONT_SIZE,
        )
        
        # Add time value below the line
        text!(
            ax,
            time,
            TIME_LABEL_OFFSET,
            text = time_str,
            align = (:center, :top),
            color = :black,
            fontsize = DEFAULT_FONT_SIZE,
        )
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
            fontsize = DEFAULT_FONT_SIZE,
        )
    end
end

"""
    _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)

Internal function to create interactive trigger timing plot.
"""
function _create_interactive_trigger_plot(trigger_codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)
    if isempty(trigger_times)
        @warn "No triggers found in the data"
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
    window_size = Observable(DEFAULT_WINDOW_SIZE)
    slider_position, slider_size, initial_position = _create_interactive_sliders(fig, end_time)
    window_position = Observable(initial_position)
    
    # Update Observables when sliders change
    on(slider_position.value) do pos
        window_position[] = pos
    end
    
    on(slider_size.value) do size
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
            _plot_trigger_events!(ax, window_times, window_codes, use_preallocated=true)
        end
        
        # Update only x-axis limits
        xlims!(ax, current_start, current_end)
    end
    
    # Set up reactive updates
    on(window_position) do _
        update_plot!()
    end
    
    on(window_size) do _
        update_plot!()
    end
    
    # Set axis properties once
    _setup_axis_properties!(ax)
    
    # Initial plot
    update_plot!()
    
    # Handle display
    display_plot = get(kwargs, :display_plot, true)
    if display_plot
        display(fig)
    end
    
    return fig, ax
end

"""
    plot_trigger_timing!(fig::Figure, ax::Axis, times::Vector{Float64}, codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)

Plot event timing with intervals between triggers.

# Arguments
- `fig`: Figure object
- `ax`: Axis object
- `times`: Vector of time intervals between triggers (unused, kept for backward compatibility)
- `codes`: Vector of trigger codes (Int16)
- `trigger_times`: Vector of absolute trigger times
- `kwargs...`: Additional keyword arguments for customization

# Returns
- Figure and Axis objects
"""
function plot_trigger_timing!(
    fig::Figure,
    ax::Axis,
    times::Vector{Float64},
    codes::Vector{Int16},
    trigger_times::Vector{Float64};
    kwargs...,
)
    _plot_trigger_events!(ax, trigger_times, codes, use_preallocated=false)
    _setup_axis_properties!(ax)
    
    # Add trigger count legend if not empty
    if !isempty(trigger_times)
        trigger_count = _count_triggers(codes)
        _add_trigger_legend_entries!(ax, trigger_count)
        fig[1, 2] = Legend(fig, ax)
    end
    
    return fig, ax
end

"""
    plot_trigger_timing(dat::BioSemiBDF.BioSemiData; kwargs...)

Plot trigger timing with interactive x-axis sliders for scrolling and window size.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemiData object containing the triggers

# Returns
- `fig::Figure`: The Makie figure object
- `ax::Axis`: The Makie axis object

# Example
```julia
fig, ax = plot_trigger_timing(dat)
```
"""
function plot_trigger_timing(dat::BioSemiBDF.BioSemiData; kwargs...)
    trigger_codes, trigger_times = _extract_trigger_data(dat)
    return _create_interactive_trigger_plot(trigger_codes, trigger_times; kwargs...)
end

"""
    plot_trigger_timing(dat::ContinuousData; kwargs...)

Plot trigger timing with interactive x-axis sliders for scrolling and window size.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing the triggers
- `kwargs...`: Additional keyword arguments for customization

# Returns
- `fig::Figure`: The Makie figure object
- `ax::Axis`: The Makie axis object

# Example
```julia
fig, ax = plot_trigger_timing(dat)
```
"""
function plot_trigger_timing(dat::ContinuousData; kwargs...)
    trigger_codes, trigger_times = _extract_trigger_data(dat)
    return _create_interactive_trigger_plot(trigger_codes, trigger_times; kwargs...)
end
