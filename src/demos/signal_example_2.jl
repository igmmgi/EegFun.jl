"""
    signal_example_2()

Interactive multi-signal composer with filtering and frequency-domain analysis.

Creates a GUI with three independent sine-wave generators, a noise source, optional
low-pass and high-pass filters, and a live power spectrum. All controls update the
plots in real time.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1–3 | Signal 1–3 | Individual sine waves. Y-axes are linked — amplitude changes are visually comparable across all three. |
| 4 | Noise | Additive Gaussian noise. |
| 5 | Combined Signal | Sum of all three signals plus noise. A red overlay shows the filtered version when a filter is active. |
| 6 | Frequency Domain | Power spectrum of the combined (or filtered) signal. |

The x-axes of plots 1–5 are linked: panning or zooming one panel updates all time-domain panels.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Freq (×3) | 0–80 Hz | Frequency of each sine wave |
| Amp (×3) | 0–10 | Amplitude of each sine wave |
| Phase (×3) | −π to π | Phase offset of each sine wave |
| Noise | 0–2 | Standard deviation of additive Gaussian noise |
| ☐ LP Filter | 0–100 Hz | Low-pass filter cutoff; enable with checkbox |
| ☐ HP Filter | 0–2 Hz | High-pass filter cutoff; enable with checkbox |

# Examples
```julia
using EegFun
EegFun.signal_example_2()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `axes::Vector{Axis}`: The six axis objects (rows 1–6)
"""
function signal_example_2()

    fig = Figure(size = (1200, 800), title = "Signal Example 2", backgroundcolor = :white)

    # Set up adaptive font and UI sizing
    function setup_adaptive_sizing(fig)
        # Create observables for different font sizes and UI elements
        title_font = Observable(24)
        label_font = Observable(20)
        tick_font = Observable(18)
        slider_font = Observable(20)
        slider_width = Observable(150)
        slider_height = Observable(20)

        # Update fonts and UI elements when figure is resized
        on(fig.scene.viewport) do area
            scale_factor = area.widths[1] / 4000  # Base on 1200px width
            title_font[] = max(16, round(Int, 24 * scale_factor))
            label_font[] = max(14, round(Int, 20 * scale_factor))
            tick_font[] = max(12, round(Int, 18 * scale_factor))
            slider_font[] = max(12, round(Int, 20 * scale_factor))
            slider_width[] = max(100, round(Int, 150 * scale_factor))
            slider_height[] = max(15, round(Int, 20 * scale_factor))
        end

        return title_font, label_font, tick_font, slider_font, slider_width, slider_height
    end

    title_font, label_font, tick_font, slider_font, slider_width, slider_height = setup_adaptive_sizing(fig)

    # Create 6 subplots
    plot_titles  = ["Signal 1", "Signal 2", "Signal 3", "Noise", "Combined Signal", "Frequency Domain"]
    plot_ylabels = ["Amplitude", "Amplitude", "Amplitude", "Amplitude", "Amplitude", "Power"]
    # Only the last time-domain plot (5) and the frequency plot (6) need x-axis labels
    plot_xlabels = ["", "", "", "", "Time (s)", "Frequency (Hz)"]

    axes = [
        Axis(
            fig[i, 3:6], # Shifted plots to columns 3:6
            title      = plot_titles[i],
            xlabel     = plot_xlabels[i],
            ylabel     = plot_ylabels[i],
            titlesize  = title_font,
            xlabelsize = label_font,
            ylabelsize = label_font,
            # Hide x tick labels on the first 4 plots — they share the same axis
            xticklabelsvisible = i <= 4 ? false : true,
        ) for i = 1:6
    ]

    # Control columns (1:2) take ~20% of figure width total — plots get the rest
    colsize!(fig.layout, 1, Relative(0.10))
    colsize!(fig.layout, 2, Relative(0.10))

    # Extract individual axes for easier reference
    ax1, ax2, ax3, ax4, ax5, ax6 = axes

    # Link x-axes of all 5 time-domain plots so pan/zoom stays in sync
    linkxaxes!(ax1, ax2, ax3, ax4, ax5)

    # Link y-axes of the 3 individual signal plots to keep amplitude scales consistent
    linkyaxes!(ax1, ax2, ax3)

    # Set up observables for parameters
    sig_dur = Observable(5.0)
    samp_rate = Observable(2000.0)

    # Signal 1 parameters
    freq1 = Observable(0.0)
    amp1 = Observable(1.0)
    phase1 = Observable(0.0)

    # Signal 2 parameters
    freq2 = Observable(0.0)
    amp2 = Observable(1.0)
    phase2 = Observable(0.0)

    # Signal 3 parameters
    freq3 = Observable(0.0)
    amp3 = Observable(1.0)
    phase3 = Observable(0.0)

    # Noise parameter
    noise_level = Observable(0.0)

    # Filter parameters
    filter_freq = Observable(0.0)
    low_pass_filter = Observable(false)
    hp_filter_freq = Observable(0.0)
    high_pass_filter = Observable(false)

    # Create the signal data as observables
    time = Observable(Float64[])
    sig1_data = Observable(Float64[])
    sig2_data = Observable(Float64[])
    sig3_data = Observable(Float64[])
    noise_data = Observable(Float64[])
    combined_data = Observable(Float64[])
    filtered_data = Observable(Float64[])
    freq_domain_data = Observable(Float64[])
    freq_domain_data_filtered = Observable(Float64[])
    freq_axis = Observable(Float64[])

    # Function to update all signals
    function update_signals()
        is_filtering = (low_pass_filter[] && filter_freq[] > 0) || (high_pass_filter[] && hp_filter_freq[] > 0)

        # Generate time vector
        t = 0:(1/samp_rate[]):(sig_dur[]-(1/samp_rate[]))
        time[] = collect(t)

        # Generate individual signals
        sig1 = amp1[] * sin.(2 * π * freq1[] * t .+ phase1[])
        sig2 = amp2[] * sin.(2 * π * freq2[] * t .+ phase2[])
        sig3 = amp3[] * sin.(2 * π * freq3[] * t .+ phase3[])

        # Generate noise
        if noise_level[] == 0
            noise = zeros(length(t))
        else
            noise = noise_level[] * randn(length(t))
        end

        # Combined signal
        combined = sig1 + sig2 + sig3 + noise
        filtered = copy(combined)

        # Apply low-pass filter if enabled
        if low_pass_filter[] && filter_freq[] > 0
            filtered = apply_lowpass_filter(filtered, samp_rate[], filter_freq[])
        end

        # Apply high-pass filter if enabled
        if high_pass_filter[] && hp_filter_freq[] > 0
            filtered = apply_highpass_filter(filtered, samp_rate[], hp_filter_freq[])
        end

        # If no filters are actually active, set filtered to zero
        if !is_filtering
            filtered = zeros(length(t))
        end

        # Frequency domain analysis - Power spectrum only
        fft_result = fft(combined)
        n = length(fft_result)
        half_n = n ÷ 2 + 1
        hz = range(0, stop = samp_rate[] / 2, length = half_n)
        freq_axis[] = collect(hz)

        # Calculate power spectrum 
        power_spectrum = (abs.(fft_result ./ n)[1:half_n] .* 2) .^ 2
        freq_domain_data[] = power_spectrum

        # Calculate filtered power spectrum if filter is enabled
        if is_filtering
            fft_filtered = fft(filtered)
            power_spectrum_filtered = (abs.(fft_filtered ./ n)[1:half_n] .* 2) .^ 2
            freq_domain_data_filtered[] = power_spectrum_filtered
        else
            freq_domain_data_filtered[] = zeros(half_n)
        end

        # Update observables
        sig1_data[] = sig1
        sig2_data[] = sig2
        sig3_data[] = sig3
        noise_data[] = noise
        combined_data[] = combined
        filtered_data[] = filtered

        # Update axis limits
        max_amp = maximum([maximum(abs.(sig1)), maximum(abs.(sig2)), maximum(abs.(sig3)), maximum(abs.(noise)), maximum(abs.(combined))])
        if max_amp > 0
            ylims!(ax1, -max_amp, max_amp)
            ylims!(ax2, -max_amp, max_amp)
            ylims!(ax3, -max_amp, max_amp)
            ylims!(ax4, -max_amp, max_amp)
            # Use a minimum range floor for ax5: near-perfect cancellation produces
            # float residuals (~1e-14) that cause "No strict ticks" warnings and
            # misleading zoomed-in noise. Fall back to the shared scale in that case.
            max_combined = maximum(abs.(combined))
            ylims!(ax5, -max(max_combined, max_amp * 1e-10), max(max_combined, max_amp * 1e-10))
        end

        # Update frequency domain limits
        max_f = maximum([freq1[], freq2[], freq3[]])
        if max_f > 0
            xlims!(ax6, 0, max_f * 1.2)
        end

        # Update frequency domain y-axis limits
        max_power = is_filtering ? maximum(freq_domain_data_filtered[]) : maximum(freq_domain_data[])

        if max_power > 0
            ylims!(ax6, 0, max_power * 1.1)
        else
            ylims!(ax6, 0, 1)
        end
    end

    function apply_lowpass_filter(signal, fs, cutoff)
        if cutoff <= 0 || cutoff >= fs / 2
            return signal
        end
        filter_info = create_lowpass_filter(cutoff, fs, order = 4)
        return filtfilt(filter_info.filter_object, signal)
    end

    # Simple high-pass filter implementation
    function apply_highpass_filter(signal, fs, cutoff)
        if cutoff <= 0 || cutoff >= fs / 2
            return signal
        end
        filter_info = create_highpass_filter(cutoff, fs, order = 4)
        return filtfilt(filter_info.filter_object, signal)
    end

    # Connect observables to update function
    for obs in [
        sig_dur,
        samp_rate,
        freq1,
        amp1,
        phase1,
        freq2,
        amp2,
        phase2,
        freq3,
        amp3,
        phase3,
        noise_level,
        filter_freq,
        low_pass_filter,
        hp_filter_freq,
        high_pass_filter,
    ]
        on(obs) do _
            update_signals()
        end
    end

    # Generate initial signal data
    update_signals()

    # Plot the signals
    lines!(ax1, time, sig1_data, color = :blue, linewidth = 2)
    lines!(ax2, time, sig2_data, color = :blue, linewidth = 2)
    lines!(ax3, time, sig3_data, color = :blue, linewidth = 2)
    lines!(ax4, time, noise_data, color = :blue, linewidth = 2)
    lines!(ax5, time, combined_data, color = :blue, linewidth = 2)

    # Add filtered signal if enabled
    filtered_line = lines!(ax5, time, filtered_data, color = :red, linewidth = 3, visible = false)

    # Frequency domain plot - original signal (only visible if no filters)
    original_freq_line =
        stem!(ax6, freq_axis, freq_domain_data, color = :blue, stemwidth = 2, trunkwidth = 2, markersize = 6, label = "Original")

    # Add filtered frequency domain if enabled
    filtered_freq_line = stem!(
        ax6,
        freq_axis,
        freq_domain_data_filtered,
        color = :blue,
        stemwidth = 2,
        trunkwidth = 2,
        markersize = 6,
        visible = false,
        label = "Filtered",
    )

    # Update filtered signal visibility
    on(low_pass_filter) do val
        is_filtering = (val && filter_freq[] > 0) || (high_pass_filter[] && hp_filter_freq[] > 0)
        filtered_line.visible = is_filtering
        filtered_freq_line.visible = is_filtering
        original_freq_line.visible = !is_filtering
    end

    on(high_pass_filter) do val
        is_filtering = (val && hp_filter_freq[] > 0) || (low_pass_filter[] && filter_freq[] > 0)
        filtered_line.visible = is_filtering
        filtered_freq_line.visible = is_filtering
        original_freq_line.visible = !is_filtering
    end

    # Control panels - one for each row to align with plots
    freq_labels = []
    freq_sliders = []
    amp_labels = []
    amp_sliders = []
    phase_labels = []
    phase_sliders = []

    for i = 1:3
        # Controls: label (fixed width) | slider (fills remaining space)
        layout = GridLayout(fig[i, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)

        freq_label  = Label(layout[1, 1], "Freq: 0.0 Hz", fontsize = slider_font, halign = :right, width = 78)
        freq_slider = Slider(layout[1, 2], range = 0.0:0.5:80.0, startvalue = 0.0, height = slider_height)
        push!(freq_labels, freq_label)
        push!(freq_sliders, freq_slider)

        amp_label  = Label(layout[2, 1], "Amp: 1.0", fontsize = slider_font, halign = :right, width = 78)
        amp_slider = Slider(layout[2, 2], range = 0.0:1.0:10.0, startvalue = 1.0, height = slider_height)
        push!(amp_labels, amp_label)
        push!(amp_sliders, amp_slider)

        phase_label  = Label(layout[3, 1], "Phase: 0.0", fontsize = slider_font, halign = :right, width = 78)
        phase_slider = Slider(layout[3, 2], range = (-π):(π/16):π, startvalue = 0.0, height = slider_height)
        push!(phase_labels, phase_label)
        push!(phase_sliders, phase_slider)
    end

    # Noise control (row 4)
    noise_layout = GridLayout(fig[4, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)
    noise_label  = Label(noise_layout[1, 1], "Noise:", fontsize = slider_font, halign = :right, width = 78)
    noise_slider = Slider(noise_layout[1, 2], range = 0.0:0.2:2.0, startvalue = 0.0, height = slider_height)

    # Filter control (row 5)
    filter_outer_layout = GridLayout(fig[5, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)

    # LP Filter: checkbox (col 1) beside bold label (col 2) on the header row
    filter_checkbox = Checkbox(filter_outer_layout[1, 1], checked = false, halign = :right, width = 78)
    Label(filter_outer_layout[1, 2], "LP Filter", fontsize = slider_font, font = :bold, halign = :left)
    filter_label  = Label(filter_outer_layout[2, 1], "0.0 Hz", fontsize = slider_font, halign = :right, width = 78)
    filter_slider = Slider(filter_outer_layout[2, 2], range = 0.0:1.0:100.0, startvalue = 0.0, height = slider_height)

    # HP Filter: same pattern
    hp_filter_checkbox = Checkbox(filter_outer_layout[3, 1], checked = false, halign = :right, width = 78)
    Label(filter_outer_layout[3, 2], "HP Filter", fontsize = slider_font, font = :bold, halign = :left)
    hp_filter_label  = Label(filter_outer_layout[4, 1], "0.0 Hz", fontsize = slider_font, halign = :right, width = 78)
    hp_filter_slider = Slider(filter_outer_layout[4, 2], range = 0.0:0.1:2.0, startvalue = 0.0, height = slider_height)

    # Connect sliders to observables
    # Signal controls (freq, amp, phase for each of 3 signals)
    freq_observables = [freq1, freq2, freq3]
    amp_observables = [amp1, amp2, amp3]
    phase_observables = [phase1, phase2, phase3]

    for i = 1:3
        # Frequency sliders
        on(freq_sliders[i].value) do val
            freq_observables[i][] = val
            freq_labels[i].text = "Freq: $(round(val, digits=1)) Hz"
        end

        # Amplitude sliders
        on(amp_sliders[i].value) do val
            amp_observables[i][] = val
            amp_labels[i].text = "Amp: $(round(val, digits=1))"
        end

        # Phase sliders
        on(phase_sliders[i].value) do val
            phase_observables[i][] = val
            phase_labels[i].text = "Phase: $(round(val, digits=2))"
        end
    end

    # Other controls
    on(noise_slider.value) do val
        noise_level[] = val
        noise_label.text = "Level: $(round(val, digits=1))"
    end

    on(filter_slider.value) do val
        filter_freq[] = val
        filter_label.text = "$(round(val, digits=1)) Hz"
    end
    on(filter_checkbox.checked) do val
        low_pass_filter[] = val
    end

    on(hp_filter_slider.value) do val
        hp_filter_freq[] = val
        hp_filter_label.text = "$(round(val, digits=1)) Hz"
    end
    on(hp_filter_checkbox.checked) do val
        high_pass_filter[] = val
    end

    display(fig)
    return fig, axes
end
