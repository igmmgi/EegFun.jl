"""
    signal_example_1()

Interactive Signal Generator

A Julia/Makie version of the MATLAB signalExample1.m GUI.
Creates an interactive plot with sliders to control:
- Signal duration (0.2 to 5 seconds)
- Frequency (1 to 100 Hz)
- Phase angle (-π to π)
- Sample rate (0 to 500 Hz)

The plot shows both the continuous signal (blue line) and sampled points (red circles).

# Example
```julia
using eegfun
eegfun.signal_example_1()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
"""
function signal_example_1()
    fig = Figure(
        size = (1000, 400),
        title = "Signal Example 1",
        backgroundcolor = :white
    )
    
    ax = Axis(fig[1, 1:6], 
              title = "Interactive Signal Generator",
              xlabel = "Time (s)",
              ylabel = "Amplitude",
              xgridvisible = true,
              ygridvisible = true,
              titlesize = 24,
              xlabelsize = 22,
              ylabelsize = 22,
              xticklabelsize = 18,
              yticklabelsize = 18)
    
    sig_dur = Observable(1.0)
    freq = Observable(10.0)
    phase = Observable(0.0)
    samp_rate = Observable(100.0)
    show_sampled = Observable(false)
    
    base_time = Observable(Float64[])
    samp_time = Observable(Float64[])
    my_sin1 = Observable(Float64[])
    my_sin2 = Observable(Float64[])
    
    # Function to update the signal data
    function update_signal()
        # Generate time vectors
        base_t = 0:0.00001:sig_dur[]  # 100,000 Hz sampling for "continuous" signal
        samp_t = 0:(1/samp_rate[]):sig_dur[]
        
        # Update observables
        base_time[] = collect(base_t)
        samp_time[] = collect(samp_t)
        my_sin1[] = sin.(2 * freq[] * π * base_t .+ phase[])
        my_sin2[] = sin.(2 * freq[] * π * samp_t .+ phase[])
        
        # Update axis limits
        xlims!(ax, 0, sig_dur[])
        ylims!(ax, -1.2, 1.2)
    end
    
    # Connect observables to update function
    for obs in [sig_dur, freq, phase, samp_rate]
        on(obs) do _
            update_signal()
        end
    end
    
    # Generate initial signal data
    update_signal()
    
    # Plot the signals
    continuous_line = lines!(ax, base_time, my_sin1, color = :blue, linewidth = 2, label = "Signal")
    sampled_line = lines!(ax, samp_time, my_sin2, color = :red, linewidth = 2, label = "Sampled Signal", visible = false)
    sample_points = scatter!(ax, samp_time, my_sin2, color = :red, markersize = 16, label = "Sample Points", visible = false)
    
    # Control visibility of sampled signal
    on(show_sampled) do val
        if val
            sampled_line.visible = true
            sample_points.visible = true
            sampled_line.color = :red
            sample_points.color = :red
            sampled_line.alpha = 1.0
            sample_points.alpha = 1.0
        else
            sampled_line.visible = false
            sample_points.visible = false
        end
    end
    
    # Set initial state (sampled signal should be hidden initially)
    show_sampled[] = false
    
    # Add static legend
    legend = axislegend(ax, position = :rt)
    
    # Create sliders and labels in a separate layout
    slider_layout = GridLayout(fig[2, 1:6], 
                              tellheight = false,
                              valign = :top,
                              padding = (1, 5, 1, 5),
                              rowgap = 1,
                              colgap = 5)
    
    # Duration slider
    duration_slider = Slider(slider_layout[1, 1], 
                            range = 0.2:0.1:5.0, 
                            startvalue = 1.0,
                            width = 300,
                            height = 25)
    duration_label = Label(slider_layout[2, 1], 
                          text = "Duration: 1.0 secs",
                          width = 300,
                          fontsize = 30)
    
    # Frequency slider
    freq_slider = Slider(slider_layout[1, 2], 
                        range = 1.0:1.0:100.0, 
                        startvalue = 10.0,
                        width = 300,
                        height = 25)
    freq_label = Label(slider_layout[2, 2], 
                      text = "Frequency: 10 Hz",
                      width = 300,
                      fontsize = 30)
    
    # Phase slider
    phase_slider = Slider(slider_layout[1, 3], 
                         range = -π:π/8:π, 
                         startvalue = 0.0,
                         width = 300,
                         height = 25)
    phase_label = Label(slider_layout[2, 3], 
                       text = "Phase Angle: 0.0",
                       width = 300,
                       fontsize = 30)
    
    # Sample rate slider
    samp_slider = Slider(slider_layout[1, 4], 
                        range = 1.0:1.0:200.0, 
                        startvalue = 100.0,
                        width = 300,
                        height = 25)
    samp_label = Label(slider_layout[2, 4], 
                      text = "Sample Rate: 100 Hz",
                      width = 300,
                      fontsize = 30)
    
    # Checkbox for showing sampled signal
    show_sampled_checkbox = Checkbox(slider_layout[1, 5], 
                                   checked = false,
                                   width = 300,
                                   height = 30)
    show_sampled_label = Label(slider_layout[2, 5], 
                              text = "Show Sampled Signal",
                              width = 300,
                              fontsize = 30)
    
    # Connect sliders to observables and labels
    on(duration_slider.value) do val
        sig_dur[] = val
        duration_label.text = "Duration: $(round(val, digits=1)) secs"
    end
    
    on(freq_slider.value) do val
        freq[] = val
        freq_label.text = "Frequency: $(round(Int, val)) Hz"
    end
    
    on(phase_slider.value) do val
        phase[] = val
        phase_label.text = "Phase Angle: $(round(val, digits=2))"
    end
    
    on(samp_slider.value) do val
        samp_rate[] = val
        samp_label.text = "Sample Rate: $(round(Int, val)) Hz"
    end
    
    # Connect checkbox to observable
    on(show_sampled_checkbox.checked) do val
        show_sampled[] = val
    end
    
    # Set row heights: 90% for plot, 10% for sliders
    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))
    
    # Display the figure
    display(fig)
    
    return fig
end
