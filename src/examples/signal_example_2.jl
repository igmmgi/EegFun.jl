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
    # Create the main figure
    fig = Figure(
        size = (1200, 800),
        title = "Signal Example 2",
        backgroundcolor = :white
    )
    
    # Create 6 subplots
    ax1 = Axis(fig[1, 2:6], title = "Signal 1", xlabel = "Time (s)", ylabel = "Amplitude",
               titlesize = 22, xlabelsize = 20, ylabelsize = 20)
    ax2 = Axis(fig[2, 2:6], title = "Signal 2", xlabel = "Time (s)", ylabel = "Amplitude",
               titlesize = 22, xlabelsize = 20, ylabelsize = 20)
    ax3 = Axis(fig[3, 2:6], title = "Signal 3", xlabel = "Time (s)", ylabel = "Amplitude",
               titlesize = 22, xlabelsize = 20, ylabelsize = 20)
    ax4 = Axis(fig[4, 2:6], title = "Noise", xlabel = "Time (s)", ylabel = "Amplitude",
               titlesize = 22, xlabelsize = 20, ylabelsize = 20)
    ax5 = Axis(fig[5, 2:6], title = "Combined Signal", xlabel = "Time (s)", ylabel = "Amplitude",
               titlesize = 22, xlabelsize = 20, ylabelsize = 20)
    ax6 = Axis(fig[6, 2:6], title = "Frequency Domain", xlabel = "Frequency (Hz)", ylabel = "Power",
               titlesize = 22, xlabelsize = 20, ylabelsize = 20)
    
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
        
        # Calculate power spectrum - power should scale as amplitude squared
        # Raw power spectrum without normalization
        power_spectrum = abs.(fft_result).^2
        # Take only the first half (up to Nyquist frequency)
        half_n = n ÷ 2 + 1
        # For single-sided spectrum, multiply by 2 (except DC and Nyquist)
        power_spectrum[2:half_n-1] .*= 2
        freq_domain_data[] = power_spectrum[1:half_n]
        freq_axis[] = hz[1:half_n]
        
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
        # Simple moving average filter as a basic low-pass
        # For production use, consider using DSP.jl with proper Butterworth filters
        window_size = max(1, Int(round(fs / cutoff / 4)))
        if window_size >= length(signal)
            return signal
        end
        
        filtered = similar(signal)
        for i in 1:length(signal)
            start_idx = max(1, i - window_size ÷ 2)
            end_idx = min(length(signal), i + window_size ÷ 2)
            filtered[i] = mean(signal[start_idx:end_idx])
        end
        return filtered
    end
    
    # Connect observables to update function
    for obs in [sig_dur, samp_rate, freq1, amp1, phase1, freq2, amp2, phase2, 
                freq3, amp3, phase3, noise_level, filter_freq, low_pass_filter]
        on(obs) do _
            update_signals()
        end
    end
    
    # Frequency domain plot will update automatically via observables
    
    # Generate initial signal data
    update_signals()
    
    # Plot the signals
    lines!(ax1, time, sig1_data, color = :blue, linewidth = 2)
    lines!(ax2, time, sig2_data, color = :blue, linewidth = 2)
    lines!(ax3, time, sig3_data, color = :blue, linewidth = 2)
    lines!(ax4, time, noise_data, color = :blue, linewidth = 2)
    lines!(ax5, time, combined_data, color = :blue, linewidth = 2)
    
    # Add filtered signal if enabled
    filtered_line = lines!(ax5, time, filtered_data, color = :red, linewidth = 2, visible = false)
    
    # Update filtered line visibility
    on(low_pass_filter) do val
        if val && filter_freq[] > 0
            filtered_line.visible = true
        else
            filtered_line.visible = false
        end
    end
    
    # Frequency domain plot - using stem for discrete points
    stem!(ax6, freq_axis, freq_domain_data, color = :blue, markersize = 4)
    
    # Add filtered frequency domain if enabled
    filtered_freq_line = stem!(ax6, freq_axis, freq_domain_data, color = :red, markersize = 4, visible = false)
    
    # Update filtered frequency line visibility
    on(low_pass_filter) do val
        if val && filter_freq[] > 0
            filtered_freq_line.visible = true
        else
            filtered_freq_line.visible = false
        end
    end
    
    # Control panels - one for each row to align with plots
    # Signal 1 controls (row 1)
    sig1_layout = GridLayout(fig[1, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(sig1_layout[1, 1], "Signal 1", fontsize = 14, font = :bold)
    freq1_label = Label(sig1_layout[2, 1], "Freq: 0.0 Hz", width = 120, fontsize = 20)
    freq1_slider = Slider(sig1_layout[3, 1], range = 0.0:0.1:100.0, startvalue = 0.0, width = 120, height = 20)
    amp1_label = Label(sig1_layout[4, 1], "Amp: 1.0", width = 120, fontsize = 20)
    amp1_slider = Slider(sig1_layout[5, 1], range = 0.0:0.1:5.0, startvalue = 1.0, width = 120, height = 20)
    phase1_label = Label(sig1_layout[6, 1], "Phase: 0.0", width = 120, fontsize = 20)
    phase1_slider = Slider(sig1_layout[7, 1], range = -π:π/16:π, startvalue = 0.0, width = 120, height = 20)
    
    # Signal 2 controls (row 2)
    sig2_layout = GridLayout(fig[2, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(sig2_layout[1, 1], "Signal 2", fontsize = 14, font = :bold)
    freq2_label = Label(sig2_layout[2, 1], "Freq: 0.0 Hz", width = 120, fontsize = 20)
    freq2_slider = Slider(sig2_layout[3, 1], range = 0.0:0.1:100.0, startvalue = 0.0, width = 120, height = 20)
    amp2_label = Label(sig2_layout[4, 1], "Amp: 1.0", width = 120, fontsize = 20)
    amp2_slider = Slider(sig2_layout[5, 1], range = 0.0:0.1:5.0, startvalue = 1.0, width = 120, height = 20)
    phase2_label = Label(sig2_layout[6, 1], "Phase: 0.0", width = 120, fontsize = 20)
    phase2_slider = Slider(sig2_layout[7, 1], range = -π:π/16:π, startvalue = 0.0, width = 120, height = 20)
    
    # Signal 3 controls (row 3)
    sig3_layout = GridLayout(fig[3, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(sig3_layout[1, 1], "Signal 3", fontsize = 14, font = :bold)
    freq3_label = Label(sig3_layout[2, 1], "Freq: 0.0 Hz", width = 120, fontsize = 20)
    freq3_slider = Slider(sig3_layout[3, 1], range = 0.0:0.1:100.0, startvalue = 0.0, width = 120, height = 20)
    amp3_label = Label(sig3_layout[4, 1], "Amp: 1.0", width = 120, fontsize = 20)
    amp3_slider = Slider(sig3_layout[5, 1], range = 0.0:0.1:5.0, startvalue = 1.0, width = 120, height = 20)
    phase3_label = Label(sig3_layout[6, 1], "Phase: 0.0", width = 120, fontsize = 20)
    phase3_slider = Slider(sig3_layout[7, 1], range = -π:π/16:π, startvalue = 0.0, width = 120, height = 20)
    
    # Noise control (row 4)
    noise_layout = GridLayout(fig[4, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(noise_layout[1, 1], "Noise", fontsize = 14, font = :bold)
    noise_label = Label(noise_layout[2, 1], "Level: 0.0", width = 120, fontsize = 20)
    noise_slider = Slider(noise_layout[3, 1], range = 0.0:0.1:10.0, startvalue = 0.0, width = 120, height = 20)
    
    # Filter control (row 5)
    filter_layout = GridLayout(fig[5, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(filter_layout[1, 1], "Filter", fontsize = 14, font = :bold)
    filter_label = Label(filter_layout[2, 1], "Cutoff: 0.0 Hz", width = 120, fontsize = 20)
    filter_slider = Slider(filter_layout[3, 1], range = 0.0:1.0:1000.0, startvalue = 0.0, width = 120, height = 20)
    filter_checkbox_label = Label(filter_layout[4, 1], "Enable LP", width = 120, fontsize = 20)
    filter_checkbox = Checkbox(filter_layout[5, 1], checked = false, width = 20, height = 20)
    
    # Output selection (row 6) - Power spectrum only
    output_layout = GridLayout(fig[6, 1], tellheight = false, valign = :center, padding = (5, 5, 5, 5), rowgap = 0)
    Label(output_layout[1, 1], "Power\nSpectrum", fontsize = 14, font = :bold)

    # Connect sliders to observables
    on(freq1_slider.value) do val
        freq1[] = val
        freq1_label.text = "Freq: $(round(val, digits=1)) Hz"
    end
    on(amp1_slider.value) do val
        amp1[] = val
        amp1_label.text = "Amp: $(round(val, digits=1))"
    end
    on(phase1_slider.value) do val
        phase1[] = val
        phase1_label.text = "Phase: $(round(val, digits=2))"
    end
    
    on(freq2_slider.value) do val
        freq2[] = val
        freq2_label.text = "Freq: $(round(val, digits=1)) Hz"
    end
    on(amp2_slider.value) do val
        amp2[] = val
        amp2_label.text = "Amp: $(round(val, digits=1))"
    end
    on(phase2_slider.value) do val
        phase2[] = val
        phase2_label.text = "Phase: $(round(val, digits=2))"
    end
    
    on(freq3_slider.value) do val
        freq3[] = val
        freq3_label.text = "Freq: $(round(val, digits=1)) Hz"
    end
    on(amp3_slider.value) do val
        amp3[] = val
        amp3_label.text = "Amp: $(round(val, digits=1))"
    end
    on(phase3_slider.value) do val
        phase3[] = val
        phase3_label.text = "Phase: $(round(val, digits=2))"
    end
    
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
    
    # Power spectrum is always shown - no callbacks needed
    
    display(fig)
    return fig
end
