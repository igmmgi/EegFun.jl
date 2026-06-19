"""
    signal_example_sampling()

Interactive Signal Generator — Nyquist Theorem Demo

A Julia/Makie version of the MATLAB signalExample1.m GUI.
Creates an interactive plot with sliders to control signal parameters and two
optional overlays for comparing reconstruction methods.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Duration | 1–10 s | Length of the displayed signal |
| Frequency | 1–100 Hz | Frequency of the sine wave |
| Phase Angle | -π to π | Phase offset of the sine wave |
| Sample Rate | 1–300 Hz | Number of samples per second |
| Linear | checkbox | Show linear (join-the-dots) reconstruction (red) |
| Sinc | checkbox | Show ideal sinc reconstruction (green) |

## Teaching Notes

This demo illustrates the **Nyquist–Shannon sampling theorem**: a signal
sampled at frequency `fs` can be perfectly reconstructed if `fs ≥ 2f`,
where `f` is the signal frequency.

### Linear reconstruction (red)
Joins consecutive sample points with straight lines. This is the simplest
possible reconstruction and introduces high-frequency artifacts at the corners.
To look convincingly smooth, linear reconstruction typically requires ≈5–10×
the signal frequency — well above the Nyquist limit.

### Sinc reconstruction (green)
Uses the Whittaker–Shannon interpolation formula: each sample is represented
by a sinc function (`sin(πt)/(πt)`) centred at its sample time, and all
contributions are summed to recover the original waveform:

```
x(t) = Σ xₙ · sinc((t - nT) / T)
```

This is the *ideal* reconstruction implied by the Nyquist theorem. When the
sample rate is ≥ 2× the signal frequency, the sinc reconstruction overlaps
the original blue signal almost exactly. The green dots show the sample
values the reconstruction is built from.

Try setting `Sample Rate = 2 × Frequency` and toggling each checkbox to see
the difference between linear and ideal reconstruction.

# See Also
- [Sinc Interpolation for Signal Reconstruction (Wolfram Demonstrations)](https://demonstrations.wolfram.com/SincInterpolationForSignalReconstruction/)

# Example
```julia
using EegFun
EegFun.signal_example_sampling()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax::Axis`: The axis object containing the plots
"""
function signal_example_sampling()
    fig = Figure(size = (1000, 800), title = "Signal Example 1")

    # Set up adaptive font and UI sizing
    """Scale font sizes and UI element dimensions when the figure is resized."""
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
            scale_factor = area.widths[1] / 4000
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

    sig_dur      = Observable(1.0)
    freq         = Observable(10.0)
    phase        = Observable(0.0)
    samp_rate    = Observable(100.0)
    show_sampled = Observable(false)
    show_sinc    = Observable(false)

    base_time  = Observable(Float64[])
    samp_time  = Observable(Float64[])
    recon_time = Observable(Float64[])
    my_sin1    = Observable(Float64[])
    my_sin2    = Observable(Float64[])
    my_sin3    = Observable(Float64[])

    # Sinc (Whittaker-Shannon) reconstruction from discrete samples
    """Whittaker–Shannon sinc interpolation: reconstruct continuous signal from discrete samples."""
    function sinc_reconstruct(t_samp, y_samp, t_out, dt)
        out = zeros(length(t_out))
        for i in eachindex(t_out)
            t = t_out[i]
            acc = 0.0
            for j in eachindex(t_samp)
                acc += y_samp[j] * sinc((t - t_samp[j]) / dt)
            end
            out[i] = acc
        end
        return out
    end

    # Function to update the signal data
    """Recompute base signal, sampled signal, and sinc reconstruction from current slider values."""
    function update_signal()
        # Generate time vectors
        base_t = 0:0.0001:sig_dur[]
        samp_t = 0:(1/samp_rate[]):sig_dur[]

        # Update observables
        base_time[] = collect(base_t)
        samp_time[] = collect(samp_t)
        my_sin1[] = sin.(2 * freq[] * π * base_t .+ phase[])
        my_sin2[] = sin.(2 * freq[] * π * samp_t .+ phase[])

        # Sinc reconstruction on a fixed 2000-point grid
        t_recon = range(0, sig_dur[], length = 2000)
        recon_time[] = collect(t_recon)
        dt = 1.0 / samp_rate[]
        my_sin3[] = sinc_reconstruct(samp_time[], my_sin2[], t_recon, dt)

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
    lines!(ax, base_time, my_sin1, color = :blue, linewidth = 2, label = "Signal")
    sampled_line  = lines!(ax, samp_time, my_sin2, color = :red, linewidth = 2, label = "Linear", visible = false)
    sample_points = scatter!(ax, samp_time, my_sin2, color = :red, markersize = 16, visible = false)
    sinc_line     = lines!(ax, recon_time, my_sin3, color = :green, linewidth = 2, label = "Sinc", visible = false)
    sinc_points   = scatter!(ax, samp_time, my_sin2, color = :green, markersize = 16, visible = false)

    # Control visibility of sampled signal (linear interpolation)
    # The checkbox/legend click sets sampled_line.visible; the observer cascades to scatter points
    on(sampled_line.visible) do val
        sample_points.visible = val
    end
    on(show_sampled) do val
        sampled_line.visible = val
    end

    # Control visibility of sinc reconstruction
    on(sinc_line.visible) do val
        sinc_points.visible = val
    end
    on(show_sinc) do val
        sinc_line.visible = val
    end

    # Set initial state
    show_sampled[] = false
    show_sinc[]    = false

    # Add static legend
    axislegend(ax, position = :rt)

    # Create sliders and labels in a separate layout
    slider_layout = GridLayout(fig[2, 1:6], tellheight = false, valign = :top, padding = (1, 5, 1, 5), rowgap = 1, colgap = 5)

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
            range = 1.0:1.0:300.0,
            startvalue = 100.0,
            label_text = "Sample Rate: 100 Hz",
            format_func = x -> "Sample Rate: $(round(Int, x)) Hz",
        ),
    ]

    # Create sliders and labels
    sliders = []
    labels = []

    for (i, config) in enumerate(slider_configs)
        slider =
            Slider(slider_layout[1, i], range = config.range, startvalue = config.startvalue, width = slider_width, height = slider_height)
        label = Label(slider_layout[2, i], text = config.label_text, width = slider_width, fontsize = slider_font)
        push!(sliders, slider)
        push!(labels, label)
    end

    # Checkbox for linear interpolation
    show_sampled_checkbox = Checkbox(slider_layout[1, 5], checked = false, halign = :center)
    Label(slider_layout[2, 5], text = "Linear", fontsize = slider_font, halign = :center)

    # Checkbox for sinc reconstruction
    show_sinc_checkbox = Checkbox(slider_layout[1, 6], checked = false, halign = :center)
    Label(slider_layout[2, 6], text = "Sinc", fontsize = slider_font, halign = :center)

    # Connect sliders to observables and labels
    observables = [sig_dur, freq, phase, samp_rate]

    for (i, config) in enumerate(slider_configs)
        on(sliders[i].value) do val
            observables[i][] = val
            labels[i].text = config.format_func(val)
        end
    end

    # Connect checkboxes to observables
    on(show_sampled_checkbox.checked) do val
        show_sampled[] = val
    end

    on(show_sinc_checkbox.checked) do val
        show_sinc[] = val
    end

    # Set row heights: 90% for plot, 10% for sliders
    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))

    # Display the figure
    display(fig)

    return fig, ax

end
