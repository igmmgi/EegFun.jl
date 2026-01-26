"""
    signal_example_2()

Interactive signal generator with multiple sine waves, noise, filtering, and frequency domain analysis.

This function creates a GUI with:
- Three independent sine wave generators (frequency, amplitude, phase controls)
- Noise generator with amplitude control
- Low-pass and High-pass filters
- Combined signal display
- Power spectrum analysis
- Real-time updates when parameters change

# Examples
```julia
using EegFun
EegFun.signal_example_2()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
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
        slider_width = Observable(120)
        slider_height = Observable(20)

        # Update fonts and UI elements when figure is resized
        on(fig.scene.viewport) do area
            scale_factor = area.widths[1] / 4000  # Base on 1200px width
            title_font[] = max(16, round(Int, 24 * scale_factor))
            label_font[] = max(14, round(Int, 20 * scale_factor))
            tick_font[] = max(12, round(Int, 18 * scale_factor))
            slider_font[] = max(12, round(Int, 20 * scale_factor))
            slider_width[] = max(80, round(Int, 120 * scale_factor))
            slider_height[] = max(15, round(Int, 20 * scale_factor))
        end

        return title_font, label_font, tick_font, slider_font, slider_width, slider_height
    end

    title_font, label_font, tick_font, slider_font, slider_width, slider_height = setup_adaptive_sizing(fig)

    # Create 6 subplots
    plot_titles = ["Signal 1", "Signal 2", "Signal 3", "Noise", "Combined Signal", "Frequency Domain"]
    plot_ylabels = ["Amplitude", "Amplitude", "Amplitude", "Amplitude", "Amplitude", "Power"]
    plot_xlabels = ["Time (s)", "Time (s)", "Time (s)", "Time (s)", "Time (s)", "Frequency (Hz)"]

    axes = [
        Axis(
            fig[i, 2:6],
            title = plot_titles[i],
            xlabel = plot_xlabels[i],
            ylabel = plot_ylabels[i],
            titlesize = title_font,
            xlabelsize = label_font,
            ylabelsize = label_font,
        ) for i = 1:6
    ]

    # Extract individual axes for easier reference
    ax1, ax2, ax3, ax4, ax5, ax6 = axes

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
            ylims!(ax5, -maximum(abs.(combined)), maximum(abs.(combined)))
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
        # Create layout for this signal
        layout = GridLayout(fig[i, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)

        # Create controls
        Label(layout[1, 1], "Signal $i", fontsize = slider_font, font = :bold)

        freq_label = Label(layout[2, 1], "Freq: 0.0 Hz", width = slider_width, fontsize = slider_font, height = 10)
        freq_slider = Slider(layout[3, 1], range = 0.0:0.5:80.0, startvalue = 0.0, width = slider_width, height = slider_height)
        push!(freq_labels, freq_label)
        push!(freq_sliders, freq_slider)

        amp_label = Label(layout[4, 1], "Amp: 1.0", width = slider_width, fontsize = slider_font, height = 10)
        amp_slider = Slider(layout[5, 1], range = 0.0:1.0:10.0, startvalue = 1.0, width = slider_width, height = slider_height)
        push!(amp_labels, amp_label)
        push!(amp_sliders, amp_slider)

        phase_label = Label(layout[6, 1], "Phase: 0.0", width = slider_width, fontsize = slider_font, height = 10)
        phase_slider = Slider(layout[7, 1], range = (-π):(π/16):π, startvalue = 0.0, width = slider_width, height = slider_height)
        push!(phase_labels, phase_label)
        push!(phase_sliders, phase_slider)
    end

    # Noise control (row 4)
    noise_layout = GridLayout(fig[4, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(noise_layout[1, 1], "Noise", fontsize = slider_font, font = :bold)
    noise_label = Label(noise_layout[2, 1], "Level: 0.0", width = slider_width, fontsize = slider_font, height = 10)
    noise_slider = Slider(noise_layout[3, 1], range = 0.0:1.0:10.0, startvalue = 0.0, width = slider_width, height = slider_height)

    # Combined Filter control next to Combined Signal (row 5)
    filter_outer_layout = GridLayout(fig[5, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 5)

    # LP Filter section
    filter_label = Label(filter_outer_layout[1, 1], "LP Filter (0.0 Hz)", fontsize = slider_font, font = :bold)
    filter_slider = Slider(filter_outer_layout[2, 1], range = 0.0:1.0:100.0, startvalue = 0.0, width = slider_width, height = slider_height)

    lp_enable_layout = GridLayout(filter_outer_layout[3, 1], tellheight = true)
    Label(lp_enable_layout[1, 1], "Enable LP", fontsize = 14)
    filter_checkbox = Checkbox(lp_enable_layout[1, 2], checked = false, width = 15, height = 15)

    # HP Filter section
    hp_filter_label = Label(filter_outer_layout[4, 1], "HP Filter (0.0 Hz)", fontsize = slider_font, font = :bold)
    hp_filter_slider =
        Slider(filter_outer_layout[5, 1], range = 0.0:0.1:2.0, startvalue = 0.0, width = slider_width, height = slider_height)

    hp_enable_layout = GridLayout(filter_outer_layout[6, 1], tellheight = true)
    Label(hp_enable_layout[1, 1], "Enable HP", fontsize = 14)
    hp_filter_checkbox = Checkbox(hp_enable_layout[1, 2], checked = false, width = 15, height = 15)

    # Output selection (row 6) - next to Frequency Domain plot
    Label(
        GridLayout(fig[6, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)[1, 1],
        "Power\nSpectrum",
        fontsize = slider_font,
        font = :bold,
    )

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
        filter_label.text = "LP Filter ($(round(val, digits=1)) Hz)"
    end
    on(filter_checkbox.checked) do val
        low_pass_filter[] = val
    end

    on(hp_filter_slider.value) do val
        hp_filter_freq[] = val
        hp_filter_label.text = "HP Filter ($(round(val, digits=1)) Hz)"
    end
    on(hp_filter_checkbox.checked) do val
        high_pass_filter[] = val
    end

    display(fig)
    return fig, axes
end
