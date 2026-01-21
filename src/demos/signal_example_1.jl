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
    fig = Figure(size = (1000, 400), title = "Signal Example 1", backgroundcolor = :white)

    # Set up adaptive font and UI sizing
    function setup_adaptive_sizing(fig)
        # Create observables for different font sizes and UI elements
        title_font = Observable(24)
        label_font = Observable(22)
        tick_font = Observable(18)
        slider_font = Observable(30)
        slider_width = Observable(300)
        slider_height = Observable(25)

        # Update fonts and UI elements when figure is resized
        on(fig.scene.viewport) do area
            scale_factor = area.widths[1] / 4000  # Base on 1000px width
            title_font[] = max(16, round(Int, 24 * scale_factor))
            label_font[] = max(14, round(Int, 22 * scale_factor))
            tick_font[] = max(12, round(Int, 18 * scale_factor))
            slider_font[] = max(20, round(Int, 30 * scale_factor))
            slider_width[] = max(200, round(Int, 300 * scale_factor))
            slider_height[] = max(15, round(Int, 25 * scale_factor))
        end

        return title_font, label_font, tick_font, slider_font, slider_width, slider_height
    end

    title_font, label_font, tick_font, slider_font, slider_width, slider_height = setup_adaptive_sizing(fig)

    ax = Axis(
        fig[1, 1:6],
        title = "Interactive Signal Generator",
        xlabel = "Time (s)",
        ylabel = "Amplitude",
        xgridvisible = true,
        ygridvisible = true,
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
    )

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
    sampled_line =
        lines!(ax, samp_time, my_sin2, color = :red, linewidth = 2, label = "Sampled Signal", visible = false)
    sample_points =
        scatter!(ax, samp_time, my_sin2, color = :red, markersize = 16, label = "Sample Points", visible = false)

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
    slider_layout =
        GridLayout(fig[2, 1:6], tellheight = false, valign = :top, padding = (1, 5, 1, 5), rowgap = 1, colgap = 5)

    # Define slider parameters
    slider_configs = [
        (
            name = "Duration",
            range = 1.0:1.0:10.0,
            startvalue = 1.0,
            label_text = "Duration: 1.0 secs",
            format_func = x -> "Duration: $(round(x, digits=1)) secs",
        ),
        (
            name = "Frequency",
            range = 1.0:1.0:100.0,
            startvalue = 10.0,
            label_text = "Frequency: 10 Hz",
            format_func = x -> "Frequency: $(round(Int, x)) Hz",
        ),
        (
            name = "Phase",
            range = (-π):(π/4):π,
            startvalue = 0.0,
            label_text = "Phase Angle: 0.0",
            format_func = x -> "Phase Angle: $(round(x, digits=2))",
        ),
        (
            name = "Sample Rate",
            range = 1.0:1.0:500.0,
            startvalue = 100.0,
            label_text = "Sample Rate: 100 Hz",
            format_func = x -> "Sample Rate: $(round(Int, x)) Hz",
        ),
    ]

    # Create sliders and labels
    sliders = []
    labels = []

    for (i, config) in enumerate(slider_configs)
        slider = Slider(
            slider_layout[1, i],
            range = config.range,
            startvalue = config.startvalue,
            width = slider_width,
            height = slider_height,
        )
        label = Label(slider_layout[2, i], text = config.label_text, width = slider_width, fontsize = slider_font)
        push!(sliders, slider)
        push!(labels, label)
    end

    # Extract individual sliders and labels for easier reference
    duration_slider, freq_slider, phase_slider, samp_slider = sliders
    duration_label, freq_label, phase_label, samp_label = labels

    # Checkbox for showing sampled signal
    show_sampled_checkbox = Checkbox(slider_layout[1, 5], checked = false, width = slider_width, height = slider_height)
    show_sampled_label =
        Label(slider_layout[2, 5], text = "Show Sampled Signal", width = slider_width, fontsize = slider_font)

    # Connect sliders to observables and labels
    observables = [sig_dur, freq, phase, samp_rate]

    for (i, config) in enumerate(slider_configs)
        on(sliders[i].value) do val
            observables[i][] = val
            labels[i].text = config.format_func(val)
        end
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
