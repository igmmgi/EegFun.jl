"""
    signal_example_composition()

Interactive multi-signal composer with filtering and frequency-domain analysis.

Creates a GUI with three sine-wave generators, additive Gaussian noise, optional Butterworth
low-pass and high-pass filters, waveform presets (square, sawtooth, phase scramble), and a
live power spectrum.

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_composition()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `axes::Vector{Axis}`: The six axis objects (rows 1–6)
"""
function signal_example_composition()
    _check_interactive_demo_backend("Signal Composition demo") || return nothing

    fig = _eegfun_figure(size = (1350, 950), title = "Signal Example 2", backgroundcolor = :white)
    title_font, label_font, tick_font, slider_font, slider_height =
        _demo_adaptive_fonts(fig, 1350, (18, 13), (15, 11), (12, 10), (14, 11), (20, 15))

    plot_titles = ("Signal 1", "Signal 2", "Signal 3", "Noise", "Combined Signal", "Frequency Domain")
    plot_ylabels = ("Amplitude", "Amplitude", "Amplitude", "Amplitude", "Amplitude", "Power")
    plot_xlabels = ("", "", "", "", "Time (s)", "Frequency (Hz)")

    axes = [
        Axis(
            fig[i, 3:6],
            title = plot_titles[i],
            xlabel = plot_xlabels[i],
            ylabel = plot_ylabels[i],
            titlesize = title_font,
            xlabelsize = label_font,
            ylabelsize = label_font,
            xticklabelsize = tick_font,
            yticklabelsize = tick_font,
            xticklabelsvisible = i > 4,
        ) for i = 1:6
    ]

    colsize!(fig.layout, 1, Relative(0.10))
    colsize!(fig.layout, 2, Relative(0.10))

    ax1, ax2, ax3, ax4, ax5, ax6 = axes
    linkxaxes!(ax1, ax2, ax3, ax4, ax5)
    linkyaxes!(ax1, ax2, ax3)

    sig_dur = Observable(5.0)
    samp_rate = Observable(2000.0)

    freq1, amp1, phase1 = Observable(10.0), Observable(1.0), Observable(0.0)
    freq2, amp2, phase2 = Observable(0.0), Observable(1.0), Observable(0.0)
    freq3, amp3, phase3 = Observable(0.0), Observable(1.0), Observable(0.0)

    noise_level = Observable(0.4)
    filter_freq = Observable(0.0)
    low_pass_filter = Observable(false)
    hp_filter_freq = Observable(0.0)
    high_pass_filter = Observable(false)

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

    function _power_spectrum_half(sig::Vector{Float64})
        fft_res = fft(sig)
        n = length(fft_res)
        half_n = n ÷ 2 + 1
        scale = 2.0 / n
        out = Vector{Float64}(undef, half_n)
        @inbounds for k = 1:half_n
            out[k] = abs2(fft_res[k] * scale)
        end
        return out
    end

    function update_signals()
        sr = samp_rate[]
        is_filtering = (low_pass_filter[] && filter_freq[] > 0) || (high_pass_filter[] && hp_filter_freq[] > 0)

        t = collect(0.0:(1.0/sr):(sig_dur[]-1.0/sr))
        n = length(t)
        time[] = t

        a1, f1, p1 = amp1[], freq1[], phase1[]
        a2, f2, p2 = amp2[], freq2[], phase2[]
        a3, f3, p3 = amp3[], freq3[], phase3[]
        sig1 = @. a1 * sin(2π * f1 * t + p1)
        sig2 = @. a2 * sin(2π * f2 * t + p2)
        sig3 = @. a3 * sin(2π * f3 * t + p3)
        nl = noise_level[]
        noise = nl == 0.0 ? zeros(Float64, n) : nl .* randn(n)
        combined = @. sig1 + sig2 + sig3 + noise

        if is_filtering
            filtered = copy(combined)
            if low_pass_filter[] && 0.0 < filter_freq[] < sr / 2
                filtered = filtfilt(create_lowpass_filter(filter_freq[], sr; order = 4).filter_object, filtered)
            end
            if high_pass_filter[] && 0.0 < hp_filter_freq[] < sr / 2
                filtered = filtfilt(create_highpass_filter(hp_filter_freq[], sr; order = 4).filter_object, filtered)
            end
        else
            filtered = zeros(Float64, n)
        end

        half_n = n ÷ 2 + 1
        freq_axis[] = collect(range(0.0, stop = sr / 2, length = half_n))
        freq_domain_data[] = _power_spectrum_half(combined)
        freq_domain_data_filtered[] = is_filtering ? _power_spectrum_half(filtered) : zeros(Float64, half_n)

        sig1_data[] = sig1
        sig2_data[] = sig2
        sig3_data[] = sig3
        noise_data[] = noise
        combined_data[] = combined
        filtered_data[] = filtered

        max_amp = max(maximum(abs, sig1), maximum(abs, sig2), maximum(abs, sig3), maximum(abs, noise), maximum(abs, combined))
        if max_amp > 0
            ylims!(ax1, -max_amp, max_amp)
            ylims!(ax2, -max_amp, max_amp)
            ylims!(ax3, -max_amp, max_amp)
            ylims!(ax4, -max_amp, max_amp)
            lim5 = max(maximum(abs, combined), max_amp * 1e-10)
            ylims!(ax5, -lim5, lim5)
        end

        max_f = max(f1, f2, f3)
        max_f > 0 && xlims!(ax6, 0, max_f * 1.2)

        max_power = is_filtering ? maximum(freq_domain_data_filtered[]) : maximum(freq_domain_data[])
        ylims!(ax6, 0, max_power > 0 ? max_power * 1.1 : 1.0)
    end

    is_applying_preset = Ref(false)
    for obs in (
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
    )
        on(obs) do _
            is_applying_preset[] || update_signals()
        end
    end

    update_signals()

    lines!(ax1, time, sig1_data, color = :blue, linewidth = 2)
    lines!(ax2, time, sig2_data, color = :blue, linewidth = 2)
    lines!(ax3, time, sig3_data, color = :blue, linewidth = 2)
    lines!(ax4, time, noise_data, color = :blue, linewidth = 2)
    lines!(ax5, time, combined_data, color = :blue, linewidth = 2)
    filtered_line = lines!(ax5, time, filtered_data, color = :red, linewidth = 3, visible = false)

    original_freq_line =
        stem!(ax6, freq_axis, freq_domain_data, color = :blue, stemwidth = 2, trunkwidth = 2, markersize = 6, label = "Original")
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

    function _sync_filter_vis()
        is_filtering = (low_pass_filter[] && filter_freq[] > 0) || (high_pass_filter[] && hp_filter_freq[] > 0)
        filtered_line.visible = is_filtering
        filtered_freq_line.visible = is_filtering
        original_freq_line.visible = !is_filtering
    end
    on(_ -> _sync_filter_vis(), low_pass_filter)
    on(_ -> _sync_filter_vis(), high_pass_filter)
    on(_ -> _sync_filter_vis(), filter_freq)
    on(_ -> _sync_filter_vis(), hp_filter_freq)

    freq_labels = Label[]
    freq_sliders = Slider[]
    amp_labels = Label[]
    amp_sliders = Slider[]
    phase_labels = Label[]
    phase_sliders = Slider[]

    freq_observables = (freq1, freq2, freq3)
    amp_observables = (amp1, amp2, amp3)
    phase_observables = (phase1, phase2, phase3)

    for i = 1:3
        layout = GridLayout(fig[i, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)
        start_f = i == 1 ? 10.0 : 0.0

        fl = Label(layout[1, 1], "Freq: $(round(start_f, digits=1)) Hz", fontsize = slider_font, halign = :right, width = 78)
        fs = Slider(layout[1, 2], range = 0.0:0.5:80.0, startvalue = start_f, height = slider_height)
        push!(freq_labels, fl)
        push!(freq_sliders, fs)

        al = Label(layout[2, 1], "Amp: 1.0", fontsize = slider_font, halign = :right, width = 78)
        as = Slider(layout[2, 2], range = 0.0:0.1:10.0, startvalue = 1.0, height = slider_height)
        push!(amp_labels, al)
        push!(amp_sliders, as)

        pl = Label(layout[3, 1], "Phase: 0.0", fontsize = slider_font, halign = :right, width = 78)
        ps = Slider(layout[3, 2], range = (-π):(π/16):π, startvalue = 0.0, height = slider_height)
        push!(phase_labels, pl)
        push!(phase_sliders, ps)

        let idx = i
            on(fs.value) do val
                freq_observables[idx][] = val
                fl.text = "Freq: $(round(val, digits=1)) Hz"
            end
            on(as.value) do val
                amp_observables[idx][] = val
                al.text = "Amp: $(round(val, digits=1))"
            end
            on(ps.value) do val
                phase_observables[idx][] = val
                pl.text = "Phase: $(round(val, digits=2))"
            end
        end
    end

    noise_layout = GridLayout(fig[4, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)
    noise_label = Label(noise_layout[1, 1], "Noise:", fontsize = slider_font, halign = :right, width = 78)
    noise_slider = Slider(noise_layout[1, 2], range = 0.0:0.2:2.0, startvalue = 0.4, height = slider_height)
    on(noise_slider.value) do val
        noise_level[] = val
        noise_label.text = "Level: $(round(val, digits=1))"
    end

    filter_outer = GridLayout(fig[5, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)
    filter_checkbox = Checkbox(filter_outer[1, 1], checked = false, halign = :right, width = 78)
    Label(filter_outer[1, 2], "LP Filter", fontsize = slider_font, font = :bold, halign = :left)
    filter_label = Label(filter_outer[2, 1], "0.0 Hz", fontsize = slider_font, halign = :right, width = 78)
    filter_slider = Slider(filter_outer[2, 2], range = 0.0:1.0:100.0, startvalue = 0.0, height = slider_height)

    hp_filter_checkbox = Checkbox(filter_outer[3, 1], checked = false, halign = :right, width = 78)
    Label(filter_outer[3, 2], "HP Filter", fontsize = slider_font, font = :bold, halign = :left)
    hp_filter_label = Label(filter_outer[4, 1], "0.0 Hz", fontsize = slider_font, halign = :right, width = 78)
    hp_filter_slider = Slider(filter_outer[4, 2], range = 0.0:0.1:2.0, startvalue = 0.0, height = slider_height)

    on(filter_slider.value) do val
        filter_freq[] = val
        filter_label.text = "$(round(val, digits=1)) Hz"
    end
    on(val -> (low_pass_filter[] = val), filter_checkbox.checked)
    on(hp_filter_slider.value) do val
        hp_filter_freq[] = val
        hp_filter_label.text = "$(round(val, digits=1)) Hz"
    end
    on(val -> (high_pass_filter[] = val), hp_filter_checkbox.checked)

    preset_layout = GridLayout(fig[6, 1:2], tellheight = false, valign = :center, padding = (4, 4, 4, 4), rowgap = 4, colgap = 6)
    Label(preset_layout[1, 1:2], "Presets", fontsize = slider_font, font = :bold)
    btn_square = Button(preset_layout[2, 1], label = "Square Wave", width = 120)
    btn_saw = Button(preset_layout[2, 2], label = "Sawtooth", width = 120)
    btn_scramble = Button(preset_layout[3, 1:2], label = "Scramble Phase (Keep Power)", width = 240)

    function _apply_preset!(f::F) where {F}
        is_applying_preset[] = true
        try
            f()
        finally
            is_applying_preset[] = false
        end
        update_signals()
    end

    on(btn_square.clicks) do _
        _apply_preset!() do
            for (i, fval, aval) in ((1, 5.0, 5.0), (2, 15.0, 5.0 / 3.0), (3, 25.0, 1.0))
                set_close_to!(freq_sliders[i], fval)
                set_close_to!(amp_sliders[i], aval)
                set_close_to!(phase_sliders[i], 0.0)
            end
            set_close_to!(noise_slider, 0.0)
            filter_checkbox.checked[] = false
        end
    end

    on(btn_saw.clicks) do _
        _apply_preset!() do
            for (i, fval, aval, pval) in ((1, 5.0, 5.0, 0.0), (2, 10.0, 2.5, π), (3, 15.0, 5.0 / 3.0, 0.0))
                set_close_to!(freq_sliders[i], fval)
                set_close_to!(amp_sliders[i], aval)
                set_close_to!(phase_sliders[i], pval)
            end
            set_close_to!(noise_slider, 0.0)
            filter_checkbox.checked[] = false
        end
    end

    on(btn_scramble.clicks) do _
        _apply_preset!() do
            for i = 1:3
                set_close_to!(phase_sliders[i], rand() * 2π - π)
            end
        end
    end

    _display_demo_figure(fig)
    return fig, axes
end
