"""
    plot_events_timing!(fig::Figure, ax::Axis, times::Vector{Float64}, codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)

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
function plot_events_timing!(fig::Figure, ax::Axis, times::Vector{Float64}, codes::Vector{Int16}, trigger_times::Vector{Float64}; kwargs...)
    # Plot horizontal timeline
    lines!(ax, [0, trigger_times[end]], [0, 0], color=:black, linewidth=2)
    
    # Plot vertical lines for each event
    for (i, (time, code)) in enumerate(zip(trigger_times, codes))
        # Draw vertical line with shorter height
        lines!(ax, [time, time], [0, 0.05], 
               color=:black, linewidth=1)
        
        # Add trigger code at the top of the line
        text!(ax, time, 0.05, 
              text=string(Int(code)),
              align=(:center, :bottom),
              color=:black,
              fontsize=16)
              
        # Add time value below the line
        text!(ax, time, -0.001, 
              text=string(time),
              align=(:center, :top),
              color=:black,
              fontsize=16)
    end
    
    # Plot intervals as text
    for i in 1:(length(trigger_times)-1)
        x = (trigger_times[i] + trigger_times[i+1]) / 2
        interval = trigger_times[i+1] - trigger_times[i]
        text!(ax, x, 0, 
              text="$interval",
              align=(:center, :top),
              color=:black,
              fontsize=16)
    end
    
    # Customize axis
    ax.xlabel = "Time (s)"
    ax.ylabel = ""
    ax.title = "Event Timing and Intervals"
    hideydecorations!(ax, ticks=true, label=false)
    ylims!(ax, -0.1, 0.15)  # Adjusted y limits to accommodate time values
    
    return fig, ax
end

"""
    plot_events_timing(dat::BioSemiBDF.BioSemiData; kwargs...)

Plot the timing of events from a BioSemiData object.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemiData object containing the triggers

# Returns
- `fig::Figure`: The Makie figure object
- `ax::Axis`: The Makie axis object

# Example
```julia
fig, ax = plot_events_timing(dat)
```
"""
function plot_events_timing(dat::BioSemiBDF.BioSemiData; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_events_timing!(fig, ax, 
                       dat.triggers.time[:, 2],  # times
                       Int16.(dat.triggers.time[:, 1]),  # codes
                       dat.time[dat.triggers.idx];  # trigger_times
                       kwargs...)
    display_plot = get(kwargs, :display_plot, true)
    if display_plot
        display(fig)
    end
    return fig, ax
end

"""
    plot_events_timing(dat::ContinuousData; kwargs...)

Plot the timing of events from a ContinuousData object.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing the triggers

# Returns
- `fig::Figure`: The Makie figure object
- `ax::Axis`: The Makie axis object

# Example
```julia
fig, ax = plot_events_timing(dat)
```
"""
function plot_events_timing(dat::ContinuousData; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    
    # Extract trigger information from the data DataFrame
    trigger_col = :triggers
    if !hasproperty(dat.data, trigger_col)
        error("No triggers column found in data. Expected column name: $trigger_col")
    end
    
    # Find trigger positions and codes
    trigger_positions = findall(x -> x != 0, dat.data[!, trigger_col])
    trigger_codes = Int16.(dat.data[trigger_positions, trigger_col])
    trigger_times = dat.data[trigger_positions, :time]
    
    # Calculate time intervals between triggers
    times = diff(trigger_times)
    
    plot_events_timing!(fig, ax, 
                       times,  # time intervals
                       trigger_codes,  # trigger codes
                       trigger_times;  # absolute times
                       kwargs...)
    display_plot = get(kwargs, :display_plot, true)
    if display_plot
        display(fig)
    end
    return fig, ax
end 