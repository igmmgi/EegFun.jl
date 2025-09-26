"""
    signal_example_2()

Interactive signal generator with multiple sine waves, noise, filtering, and frequency domain analysis.

This function creates a GUI with:
- Three independent sine wave generators (frequency, amplitude, phase controls)
- Noise generator with amplitude control
- Low-pass filter
- Combined signal display
- Power spectrum or phase spectrum analysis
- Real-time updates when parameters change

# Examples
```julia
using eegfun
eegfun.signal_example_2()
```
"""
function signal_example_2()

    fig = Figure(
        size = (1200, 800),
        title = "Signal Example 2",
        backgroundcolor = :white
    )
    
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
    
    axes = [Axis(fig[i, 2:6], 
                 title = plot_titles[i], 
                 xlabel = plot_xlabels[i], 
                 ylabel = plot_ylabels[i],
                 titlesize = title_font, xlabelsize = label_font, ylabelsize = label_font) 
            for i in 1:6]
    
    # Extract individual axes for easier reference
    ax1, ax2, ax3, ax4, ax5, ax6 = axes
    
    # Set up observables for parameters
    sig_dur = Observable(3.0)
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
    
    # Filter parameter
    filter_freq = Observable(0.0)
    low_pass_filter = Observable(false)
    
    # Output selection (power spectrum only)
    show_power = Observable(true)
    
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
        # Generate time vector
        t = 0:(1/samp_rate[]):sig_dur[]-(1/samp_rate[])
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
        
        # Apply low-pass filter if enabled
        if low_pass_filter[] && filter_freq[] > 0
            # Simple low-pass filter implementation
            # For a more sophisticated filter, you could use DSP.jl
            filtered = apply_lowpass_filter(combined, samp_rate[], filter_freq[])
        else
            filtered = zeros(length(t))  # Empty array of same length as time
        end
        
        # Frequency domain analysis - Power spectrum only
        fft_result = fft(combined)
        n = length(fft_result)
        hz = collect(0:samp_rate[]/n:samp_rate[]*(n-1)/n)
        
        # Calculate power spectrum 
        fft_normalized = fft_result ./ n
        half_n = n ÷ 2 + 1
        power_spectrum = abs.(fft_normalized[1:half_n]) .* 2
        power_spectrum = power_spectrum.^2
        freq_domain_data[] = power_spectrum[1:half_n]
        freq_axis[] = hz[1:half_n]
        
        # Calculate filtered power spectrum if filter is enabled
        if low_pass_filter[] && filter_freq[] > 0
            fft_filtered = fft(filtered)
            fft_filtered_normalized = fft_filtered ./ n
            power_spectrum_filtered = abs.(fft_filtered_normalized[1:half_n]) .* 2
            power_spectrum_filtered = power_spectrum_filtered.^2
            freq_domain_data_filtered[] = power_spectrum_filtered[1:half_n]
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
        max_amp = maximum([maximum(abs.(sig1)), maximum(abs.(sig2)), maximum(abs.(sig3)), 
                          maximum(abs.(noise)), maximum(abs.(combined))])
        if max_amp > 0
            ylims!(ax1, -max_amp, max_amp)
            ylims!(ax2, -max_amp, max_amp)
            ylims!(ax3, -max_amp, max_amp)
            ylims!(ax4, -max_amp, max_amp)
            ylims!(ax5, -maximum(abs.(combined)), maximum(abs.(combined)))
        end
        
        # Update frequency domain limits
        max_freq = maximum([freq1[], freq2[], freq3[]])
        if max_freq > 0
            xlims!(ax6, 0, max_freq + 0.2 * max_freq)
        end
    end
    
    # Simple low-pass filter implementation
    function apply_lowpass_filter(signal, fs, cutoff)
        if cutoff <= 0 || cutoff >= fs/2
            return signal
        end
        filter_info = create_filter("lp", "iir", cutoff, fs, order = 4)

        return filtfilt(filter_info.filter_object, signal)
        
    end
    
    # Connect observables to update function
    for obs in [sig_dur, samp_rate, freq1, amp1, phase1, freq2, amp2, phase2, 
                freq3, amp3, phase3, noise_level, filter_freq, low_pass_filter]
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
    
    # Update filtered line visibility
    on(low_pass_filter) do val
        if val && filter_freq[] > 0
            filtered_line.visible = true
        else
            filtered_line.visible = false
        end
    end
    
    # Frequency domain plot - original signal (always visible)
    stem!(ax6, freq_axis, freq_domain_data, color = :blue, markersize = 4, label = "Original")
    
    # Add filtered frequency domain if enabled
    filtered_freq_line = stem!(ax6, freq_axis, freq_domain_data_filtered, color = :red, markersize = 4, visible = false, label = "Filtered")
    
    # Update y-axis limits when power data changes
    on(freq_domain_data) do data
        if length(data) > 0
            max_power = maximum(data)
            if max_power > 0
                ylims!(ax6, 0, max_power * 1.1)  # Add 10% margin
            else
                ylims!(ax6, 0, 1)  # Default range when all power is zero
            end
        end
    end
    
    # Also update y-axis limits when filtered power data changes
    on(freq_domain_data_filtered) do data
        if length(data) > 0
            max_original = maximum(freq_domain_data[])
            max_filtered = maximum(data)
            max_power = max(max_original, max_filtered)
            if max_power > 0
                ylims!(ax6, 0, max_power * 1.1)  # Add 10% margin
            else
                ylims!(ax6, 0, 1)  # Default range when all power is zero
            end
        end
    end
    
    # Update filtered frequency line visibility
    on(low_pass_filter) do val
        if val && filter_freq[] > 0
            filtered_freq_line.visible = true
        else
            filtered_freq_line.visible = false
        end
    end
    
    # Control panels - one for each row to align with plots
    signal_layouts = []
    freq_labels = []
    freq_sliders = []
    amp_labels = []
    amp_sliders = []
    phase_labels = []
    phase_sliders = []
    
    for i in 1:3
        # Create layout for this signal
        layout = GridLayout(fig[i, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
        push!(signal_layouts, layout)
        
        # Create controls
        Label(layout[1, 1], "Signal $i", fontsize = slider_font, font = :bold)
        
        freq_label = Label(layout[2, 1], "Freq: 0.0 Hz", width = slider_width, fontsize = slider_font, height = 10)
        freq_slider = Slider(layout[3, 1], range = 0.0:1.0:100.0, startvalue = 0.0, width = slider_width, height = slider_height)
        push!(freq_labels, freq_label)
        push!(freq_sliders, freq_slider)
        
        amp_label = Label(layout[4, 1], "Amp: 1.0", width = slider_width, fontsize = slider_font, height = 10)
        amp_slider = Slider(layout[5, 1], range = 0.0:1.0:10.0, startvalue = 1.0, width = slider_width, height = slider_height)
        push!(amp_labels, amp_label)
        push!(amp_sliders, amp_slider)
        
        phase_label = Label(layout[6, 1], "Phase: 0.0", width = slider_width, fontsize = slider_font, height = 10)
        phase_slider = Slider(layout[7, 1], range = -π:π/16:π, startvalue = 0.0, width = slider_width, height = slider_height)
        push!(phase_labels, phase_label)
        push!(phase_sliders, phase_slider)
    end
    
    # Extract individual controls for easier reference
    sig1_layout, sig2_layout, sig3_layout = signal_layouts
    freq1_label, freq2_label, freq3_label = freq_labels
    freq1_slider, freq2_slider, freq3_slider = freq_sliders
    amp1_label, amp2_label, amp3_label = amp_labels
    amp1_slider, amp2_slider, amp3_slider = amp_sliders
    phase1_label, phase2_label, phase3_label = phase_labels
    phase1_slider, phase2_slider, phase3_slider = phase_sliders
    
    # Noise control (row 4)
    noise_layout = GridLayout(fig[4, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(noise_layout[1, 1], "Noise", fontsize = slider_font, font = :bold)
    noise_label = Label(noise_layout[2, 1], "Level: 0.0", width = slider_width, fontsize = slider_font, height = 10)
    noise_slider = Slider(noise_layout[3, 1], range = 0.0:1.0:10.0, startvalue = 0.0, width = slider_width, height = slider_height)
    
    # Filter control (row 5)
    filter_layout = GridLayout(fig[5, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(filter_layout[1, 1], "LP Filter", fontsize = slider_font, font = :bold)
    filter_label = Label(filter_layout[2, 1], "Cutoff: 0.0 Hz", width = slider_width, fontsize = slider_font, height = 10)
    filter_slider = Slider(filter_layout[3, 1], range = 0.0:1.0:100.0, startvalue = 0.0, width = slider_width, height = slider_height)
    filter_checkbox_label = Label(filter_layout[4, 1], "Enable LP", width = slider_width, fontsize = slider_font)
    filter_checkbox = Checkbox(filter_layout[5, 1], checked = false, width = slider_height, height = slider_height)
    
    # Output selection (row 6) - Power spectrum only
    output_layout = GridLayout(fig[6, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(output_layout[1, 1], "Power\nSpectrum", fontsize = slider_font, font = :bold)

    # Connect sliders to observables
    # Signal controls (freq, amp, phase for each of 3 signals)
    freq_observables = [freq1, freq2, freq3]
    amp_observables = [amp1, amp2, amp3]
    phase_observables = [phase1, phase2, phase3]
    
    for i in 1:3
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
        filter_label.text = "Cutoff: $(round(val, digits=1)) Hz"
    end
    on(filter_checkbox.checked) do val
        low_pass_filter[] = val
    end
    
    display(fig)
    return fig
end
