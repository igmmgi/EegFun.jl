########################################################
# Module constants
########################################################

# Plotting constants
const DEFAULT_TIMELINE_WIDTH = 2
const DEFAULT_EVENT_LINE_WIDTH = 1
const DEFAULT_FONT_SIZE = 24
const EVENT_LINE_HEIGHT = 0.05
const TIME_LABEL_OFFSET = -0.001
const Y_MIN_LIMIT = -0.1
const Y_MAX_LIMIT = 0.15

# Slider constants
const DEFAULT_WINDOW_SIZE = 10.0  # seconds
const MIN_WINDOW_SIZE = 1.0       # seconds
const MAX_WINDOW_SIZE = 60.0      # seconds
const WINDOW_SIZE_STEP = 1.0      # seconds
const POSITION_STEP = 0.5         # seconds

# Pre-allocated vectors for performance
const VERTICAL_LINE_X = Ref([0.0, 0.0])
const VERTICAL_LINE_Y = Ref([0.0, EVENT_LINE_HEIGHT])

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
function plot_trigger_timing!(fig::Figure, ax::Axis, times::Vector{Float64}, codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)
    # Pre-compute string conversions to avoid repeated operations in loops
    code_strings = [string(Int(code)) for code in codes]
    time_strings = [string(time) for time in trigger_times]
    
    # Pre-compute intervals and their string representations
    intervals = diff(trigger_times)
    interval_strings = [string(interval) for interval in intervals]
    
    # Pre-compute interval positions
    interval_positions = [(trigger_times[i] + trigger_times[i+1]) / 2 for i in 1:length(intervals)]
    
    # Plot horizontal timeline
    lines!(ax, [0, trigger_times[end]], [0, 0], color=:black, linewidth=DEFAULT_TIMELINE_WIDTH)
    
    # Plot vertical lines for each event
    for (i, (time, code_str, time_str)) in enumerate(zip(trigger_times, code_strings, time_strings))
        # Draw vertical line with shorter height
        lines!(ax, [time, time], [0, EVENT_LINE_HEIGHT], 
               color=:black, linewidth=DEFAULT_EVENT_LINE_WIDTH)
        
        # Add trigger code at the top of the line
        text!(ax, time, EVENT_LINE_HEIGHT, 
              text=code_str,
              align=(:center, :bottom),
              color=:black,
              fontsize=DEFAULT_FONT_SIZE)
              
        # Add time value below the line
        text!(ax, time, TIME_LABEL_OFFSET, 
              text=time_str,
              align=(:center, :top),
              color=:black,
              fontsize=DEFAULT_FONT_SIZE)
    end
    
    # Plot intervals as text
    for (x, interval_str) in zip(interval_positions, interval_strings)
        text!(ax, x, 0, 
              text=interval_str,
              align=(:center, :top),
              color=:black,
              fontsize=DEFAULT_FONT_SIZE)
    end
    
    # Customize axis
    ax.xlabel = "Time (s)"
    ax.ylabel = ""
    ax.title = "Event Timing and Intervals"
    ax.xlabelsize = 20  # Larger axis label
    ax.titlesize = 22   # Larger title
    hideydecorations!(ax, ticks=true, label=false)
    ylims!(ax, Y_MIN_LIMIT, Y_MAX_LIMIT)  # Adjusted y limits to accommodate time values
    
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
    # Extract trigger information from BioSemi data
    raw_triggers = dat.triggers.raw
    
    # Clean triggers to detect only onset events
    cleaned_triggers = _clean_triggers(raw_triggers)
    
    # Find trigger positions and codes from cleaned triggers
    trigger_positions = findall(x -> x != 0, cleaned_triggers)
    trigger_codes = Int16.(cleaned_triggers[trigger_positions])
    trigger_times = dat.time[trigger_positions]
    
    if isempty(trigger_times)
        @warn "No triggers found in the data"
        fig = Figure()
        ax = Axis(fig[1, 1])
        return fig, ax
    end
    
    # Calculate time range
    total_time = trigger_times[end] - trigger_times[1]
    end_time = trigger_times[end] + 2.0
    
    # Create figure with grid layout
    fig = Figure()
    
    # Create axis for the main plot
    ax = Axis(fig[1, 1])
    
    # Create Observables for reactive plotting
    window_size = Observable(DEFAULT_WINDOW_SIZE)
    initial_position = -2.0
    window_position = Observable(initial_position)
    
    # Create sliders below the plot
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
        snap = true
    )
    
    # Create labels for sliders
    position_label = Label(fig[2, 2], @lift("Position: $(round($(slider_position.value), digits=1))s"), fontsize=20)
    size_label = Label(fig[3, 2], @lift("Window Size: $(round($(slider_size.value), digits=1))s"), fontsize=20)
    
    # Update Observables when sliders change
    on(slider_position.value) do pos
        window_position[] = pos
    end
    
    on(slider_size.value) do size
        window_size[] = size
    end
    
    # Reactive plotting function
    function update_plot!()
        # Clear the axis
        empty!(ax)
        
        # Get current window bounds
        current_start = window_position[]
        current_end = min(current_start + window_size[], end_time)
        
        # Filter triggers within the current window
        window_mask = (trigger_times .>= current_start) .&& (trigger_times .<= current_end)
        window_times = trigger_times[window_mask]
        window_codes = trigger_codes[window_mask]
        
        if !isempty(window_times)
            # Pre-compute string conversions
            code_strings = [string(Int(code)) for code in window_codes]
            time_strings = [string(round(time, digits=2)) for time in window_times]
            
            # Pre-compute intervals and their string representations
            intervals = diff(window_times)
            interval_strings = [string(round(interval, digits=2)) for interval in intervals]
            interval_positions = [(window_times[i] + window_times[i+1]) / 2 for i in 1:length(intervals)]
            
            # Plot horizontal timeline
            timeline_start = 0.0  # Start from beginning of data
            lines!(ax, [timeline_start, current_end], [0, 0], color=:black, linewidth=DEFAULT_TIMELINE_WIDTH)
            
            # Plot vertical lines for each event
            for (time, code_str, time_str) in zip(window_times, code_strings, time_strings)
                # Draw vertical line using pre-allocated vectors
                VERTICAL_LINE_X[][1] = VERTICAL_LINE_X[][2] = time
                lines!(ax, VERTICAL_LINE_X[], VERTICAL_LINE_Y[], 
                       color=:black, linewidth=DEFAULT_EVENT_LINE_WIDTH)
                
                # Add trigger code at the top of the line
                text!(ax, time, EVENT_LINE_HEIGHT, 
                      text=code_str,
                      align=(:center, :bottom),
                      color=:black,
                      fontsize=DEFAULT_FONT_SIZE)
                      
                # Add time value below the line
                text!(ax, time, TIME_LABEL_OFFSET, 
                      text=time_str,
                      align=(:center, :top),
                      color=:black,
                      fontsize=DEFAULT_FONT_SIZE)
            end
            
            # Plot intervals as text
            for (x, interval_str) in zip(interval_positions, interval_strings)
                text!(ax, x, 0, 
                      text=interval_str,
                      align=(:center, :top),
                      color=:black,
                      fontsize=DEFAULT_FONT_SIZE)
            end
        end
        
        # Update only x-axis limits (other properties set once outside)
        xlims!(ax, current_start, current_end)
    end
    
    # Set up reactive updates
    on(window_position) do _
        update_plot!()
    end
    
    on(window_size) do _
        update_plot!()
    end
    
    # Set axis properties once (outside update function for performance)
    ax.xlabel = "Time (s)"
    ax.ylabel = ""
    ax.title = "Event Timing and Intervals"
    ax.xlabelsize = 20  # Larger axis label
    ax.titlesize = 22   # Larger title
    hideydecorations!(ax, ticks=true, label=false)
    ylims!(ax, Y_MIN_LIMIT, Y_MAX_LIMIT)
    
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
    # Extract trigger information
    trigger_col = :triggers
    if !hasproperty(dat.data, trigger_col)
        error("No triggers column found in data. Expected column name: $trigger_col")
    end
    
    # Since triggers are already cleaned in ContinuousData, just find non-zero values
    trigger_positions = findall(x -> x != 0, dat.data[!, trigger_col])
    trigger_codes = Int16.(dat.data[trigger_positions, trigger_col])
    trigger_times = dat.data[trigger_positions, :time]
    
    if isempty(trigger_times)
        @warn "No triggers found in the data"
        fig = Figure()
        ax = Axis(fig[1, 1])
        return fig, ax
    end
    
    # Calculate time range
    total_time = trigger_times[end] - trigger_times[1]
    end_time = trigger_times[end] + 2.0
    
    # Create figure with grid layout
    fig = Figure()
    
    # Create axis for the main plot
    ax = Axis(fig[1, 1])
    
    # Create Observables for reactive plotting
    window_size = Observable(DEFAULT_WINDOW_SIZE)
    initial_position = -2.0
    window_position = Observable(initial_position)
    
    # Create sliders below the plot
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
        snap = true
    )
    
    # Create labels for sliders
    position_label = Label(fig[2, 2], @lift("Position: $(round($(slider_position.value), digits=1))s"), fontsize=20)
    size_label = Label(fig[3, 2], @lift("Window Size: $(round($(slider_size.value), digits=1))s"), fontsize=20)
    
    # Update Observables when sliders change
    on(slider_position.value) do pos
        window_position[] = pos
    end
    
    on(slider_size.value) do size
        window_size[] = size
    end
    
    # Reactive plotting function
    function update_plot!()
        # Clear the axis
        empty!(ax)
        
        # Get current window bounds
        current_start = window_position[]
        current_end = min(current_start + window_size[], end_time)
        
        # Filter triggers within the current window
        window_mask = (trigger_times .>= current_start) .&& (trigger_times .<= current_end)
        window_times = trigger_times[window_mask]
        window_codes = trigger_codes[window_mask]
        
        if !isempty(window_times)
            # Pre-compute string conversions
            code_strings = [string(Int(code)) for code in window_codes]
            time_strings = [string(round(time, digits=2)) for time in window_times]
            
            # Pre-compute intervals and their string representations
            intervals = diff(window_times)
            interval_strings = [string(round(interval, digits=2)) for interval in intervals]
            interval_positions = [(window_times[i] + window_times[i+1]) / 2 for i in 1:length(intervals)]
            
            # Plot horizontal timeline
            timeline_start = 0.0  # Start from beginning of data or 5s before first trigger
            lines!(ax, [timeline_start, current_end], [0, 0], color=:black, linewidth=DEFAULT_TIMELINE_WIDTH)
            
            # Plot vertical lines for each event
            for (time, code_str, time_str) in zip(window_times, code_strings, time_strings)
                # Draw vertical line using pre-allocated vectors
                VERTICAL_LINE_X[][1] = VERTICAL_LINE_X[][2] = time
                lines!(ax, VERTICAL_LINE_X[], VERTICAL_LINE_Y[], 
                       color=:black, linewidth=DEFAULT_EVENT_LINE_WIDTH)
                
                # Add trigger code at the top of the line
                text!(ax, time, EVENT_LINE_HEIGHT, 
                      text=code_str,
                      align=(:center, :bottom),
                      color=:black,
                      fontsize=DEFAULT_FONT_SIZE)
                      
                # Add time value below the line
                text!(ax, time, TIME_LABEL_OFFSET, 
                      text=time_str,
                      align=(:center, :top),
                      color=:black,
                      fontsize=DEFAULT_FONT_SIZE)
            end
            
            # Plot intervals as text
            for (x, interval_str) in zip(interval_positions, interval_strings)
                text!(ax, x, 0, 
                      text=interval_str,
                      align=(:center, :top),
                      color=:black,
                      fontsize=DEFAULT_FONT_SIZE)
            end
        end
        
        # Update only x-axis limits (other properties set once outside)
        xlims!(ax, current_start, current_end)
    end
    
    # Set up reactive updates
    on(window_position) do _
        update_plot!()
    end
    
    on(window_size) do _
        update_plot!()
    end
    
    # Set axis properties once (outside update function for performance)
    ax.xlabel = "Time (s)"
    ax.ylabel = ""
    ax.title = "Event Timing and Intervals"
    ax.xlabelsize = 20  # Larger axis label
    ax.titlesize = 22   # Larger title
    hideydecorations!(ax, ticks=true, label=false)
    ylims!(ax, Y_MIN_LIMIT, Y_MAX_LIMIT)
    
    # Initial plot
    update_plot!()
    
    # Handle display
    display_plot = get(kwargs, :display_plot, true)
    if display_plot
        display(fig)
    end
    
    return fig, ax
end 