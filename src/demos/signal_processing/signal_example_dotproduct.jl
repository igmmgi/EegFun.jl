"""
    signal_example_dotproduct()

Interactive dot-product demo illustrating how the Discrete Fourier Transform (DFT) detects
frequency power via sine and complex (sine + cosine) probe sinusoids.

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_dotproduct()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_sig, ax_prod_sin, ax_prod_cos)` — the three axis objects
"""
function signal_example_dotproduct()
    _check_interactive_demo_backend("Signal Dot Product demo") || return nothing

    FS = 1000.0
    EPOCH = 2.0

    fig = _eegfun_figure(size = (1350, 850), figure_padding = (8, 8, 8, 8), title = "Dot Product Demo")
    lbl_sz, tick_sz, title_sz, ctrl_sz = _demo_adaptive_fonts(fig, 1350, (15, 11), (12, 10), (18, 13), (15, 11))

    ax_common = (
        xlabel = "Time (s)",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible = true,
        ygridvisible = true,
    )
    ax_sig = Axis(fig[1, 1]; title = "Signal (blue)  +  Test sinusoid (red)", ylabel = "Amplitude", ax_common...)
    ax_prod_sin = Axis(fig[2, 1]; title = "Sine product  (signal × sin)", ylabel = "", ax_common...)
    ax_prod_cos = Axis(fig[3, 1]; title = "Cosine product  (signal × cos)", ylabel = "", ax_common...)

    sig_freq = Observable(10.0)
    sig_amp = Observable(1.0)
    sig_phase = Observable(0.0)
    noise_lvl = Observable(0.0)
    test_freq = Observable(10.0)
    use_complex = Observable(false)

    t_arr = collect(range(0.0, EPOCH - 1.0 / FS, step = 1.0 / FS))
    n = length(t_arr)

    sig_buf = Vector{Float64}(undef, n)
    sin_buf = Vector{Float64}(undef, n)
    cos_buf = Vector{Float64}(undef, n)
    psin_buf = Vector{Float64}(undef, n)
    pcos_buf = Vector{Float64}(undef, n)
    sig_data = Observable(sig_buf)
    test_sin_data = Observable(sin_buf)
    test_cos_data = Observable(cos_buf)
    prod_sin_data = Observable(psin_buf)
    prod_cos_data = Observable(pcos_buf)
    zero_vec = zeros(Float64, n)

    function update(; do_notify::Bool = true)
        amp = sig_amp[]
        sf = sig_freq[]
        ph = sig_phase[]
        tf = test_freq[]
        nl = noise_lvl[]
        ω_sig = 2π * sf
        ω_test = 2π * tf

        @inbounds for i = 1:n
            ti = t_arr[i]
            sv = amp * sin(muladd(ω_sig, ti, ph))
            if nl > 0.0
                sv = muladd(nl, randn(), sv)
            end
            ts, tc = sincos(ω_test * ti)
            sig_buf[i] = sv
            sin_buf[i] = ts
            cos_buf[i] = tc
            psin_buf[i] = sv * ts
            pcos_buf[i] = sv * tc
        end

        if do_notify
            notify(sig_data)
            notify(test_sin_data)
            notify(test_cos_data)
            notify(prod_sin_data)
            notify(prod_cos_data)
        end

        xlims!(ax_sig, 0.0, EPOCH)
        xlims!(ax_prod_sin, 0.0, EPOCH)
        xlims!(ax_prod_cos, 0.0, EPOCH)
        a = amp + nl
        ylims!(ax_sig, -a * 1.1, a * 1.1)

        for (ax, pdata) in ((ax_prod_sin, psin_buf), (ax_prod_cos, pcos_buf))
            p_lo, p_hi = extrema(pdata)
            p_min = min(p_lo, -0.05)
            p_max = max(p_hi, 0.05)
            p_margin = (p_max - p_min) * 0.12
            ylims!(ax, p_min - p_margin, p_max + p_margin)
        end
    end

    is_applying_preset = Ref(false)
    for obs in (sig_freq, sig_amp, sig_phase, noise_lvl, test_freq)
        on(obs) do _
            is_applying_preset[] || update()
        end
    end
    update(; do_notify = false)

    lines!(ax_sig, t_arr, sig_data, color = :royalblue, linewidth = 2, label = "Signal")
    lines!(ax_sig, t_arr, test_sin_data, color = :crimson, linewidth = 2, label = "Sine test")
    cos_line = lines!(ax_sig, t_arr, test_cos_data, color = :seagreen, linewidth = 2, label = "Cosine test", visible = false)
    axislegend(ax_sig, position = :rt)

    band!(ax_prod_sin, t_arr, zero_vec, prod_sin_data, color = (:gray, 0.4))
    lines!(ax_prod_sin, t_arr, prod_sin_data, color = :black, linewidth = 2)
    hlines!(ax_prod_sin, [0.0], color = :black, linewidth = 1.0, linestyle = :dash)

    band!(ax_prod_cos, t_arr, zero_vec, prod_cos_data, color = (:gray, 0.4))
    lines!(ax_prod_cos, t_arr, prod_cos_data, color = :black, linewidth = 2)
    hlines!(ax_prod_cos, [0.0], color = :black, linewidth = 1.0, linestyle = :dash)

    dp_sin_text = @lift begin
        dp = sum($prod_sin_data) / n * 2.0
        freq_match = abs($sig_freq - $test_freq) < 0.5
        phase_deg = mod(rad2deg($sig_phase), 360.0)
        phase_match = phase_deg < 10.0 || phase_deg > 350.0
        freq_str = freq_match ? "Freq: ✓ Match" : "Freq: ✗ Mismatch"
        phase_str = phase_match ? "Phase: ✓ Match" : "Phase: ✗ Mismatch"
        prefix = $use_complex ? "Sine dp" : "Dot product"
        "$prefix = $(round(dp, digits=3))     $freq_str     $phase_str"
    end
    text!(
        ax_prod_sin,
        0.5,
        0.97,
        text = dp_sin_text,
        space = :relative,
        fontsize = @lift(max(14, round(Int, $lbl_sz * 1.15))),
        color = :black,
        align = (:center, :top),
    )

    dp_cos_text = @lift("Cosine dp = $(round(sum($prod_cos_data) / n * 2.0, digits=3))")
    text!(
        ax_prod_cos,
        0.5,
        0.97,
        text = dp_cos_text,
        space = :relative,
        fontsize = @lift(max(14, round(Int, $lbl_sz * 1.15))),
        color = :black,
        align = (:center, :top),
    )

    mag_text = @lift begin
        $use_complex || return ""
        dp_sin = sum($prod_sin_data) / n * 2.0
        dp_cos = sum($prod_cos_data) / n * 2.0
        "Magnitude = √(sin² + cos²) = $(round(hypot(dp_sin, dp_cos), digits=3))"
    end
    text!(
        ax_prod_sin,
        0.5,
        0.03,
        text = mag_text,
        space = :relative,
        fontsize = @lift(max(14, round(Int, $lbl_sz * 1.25))),
        color = :black,
        font = :bold,
        align = (:center, :bottom),
    )

    on(use_complex) do is_complex
        ax_prod_cos.blockscene.visible[] = is_complex
        cos_line.visible = is_complex
        if is_complex
            rowsize!(fig.layout, 1, Relative(0.26))
            rowsize!(fig.layout, 2, Relative(0.26))
            rowsize!(fig.layout, 3, Relative(0.26))
            rowsize!(fig.layout, 4, Relative(0.15))
            rowsize!(fig.layout, 5, Relative(0.07))
            ax_sig.title = "Signal (blue)  +  Sine test (red)  +  Cosine test (green)"
        else
            rowsize!(fig.layout, 1, Relative(0.35))
            rowsize!(fig.layout, 2, Relative(0.40))
            rowsize!(fig.layout, 3, Relative(0.0))
            rowsize!(fig.layout, 4, Relative(0.18))
            rowsize!(fig.layout, 5, Relative(0.07))
            ax_sig.title = "Signal (blue)  +  Test sinusoid (red)"
        end
        is_applying_preset[] || update()
    end

    ctrl = GridLayout(fig[4, 1], colgap = 20)
    sl_sf, lbl_sf = _demo_labelled_slider(ctrl, 1, "Signal Freq", 1.0:0.1:40.0, 10.0, v -> "$(round(v, digits=1)) Hz", ctrl_sz)
    sl_sa, lbl_sa = _demo_labelled_slider(ctrl, 2, "Amplitude", 0.1:0.1:2.0, 1.0, v -> "$(round(v, digits=1))", ctrl_sz)
    sl_sp, lbl_sp = _demo_labelled_slider(ctrl, 3, "Phase (°)", 0.0:5.0:355.0, 0.0, v -> "$(round(Int, v))°", ctrl_sz)
    sl_tf, lbl_tf = _demo_labelled_slider(ctrl, 4, "Test Freq", 0.5:0.1:60.0, 10.0, v -> "$(round(v, digits=1)) Hz", ctrl_sz)
    sl_nl, lbl_nl = _demo_labelled_slider(ctrl, 5, "Noise", 0.0:0.05:1.0, 0.0, v -> "$(round(v, digits=2))", ctrl_sz)

    Label(ctrl[1, 6], "Complex", fontsize = ctrl_sz, halign = :center)
    tog_complex = Toggle(ctrl[2, 6], active = false, halign = :center)
    Label(ctrl[3, 6], "", fontsize = ctrl_sz, halign = :center)

    on(sl_sf.value) do v
        sig_freq[] = v
        lbl_sf.text = "$(round(v, digits=1)) Hz"
    end
    on(sl_sa.value) do v
        sig_amp[] = v
        lbl_sa.text = "$(round(v, digits=1))"
    end
    on(sl_sp.value) do v
        sig_phase[] = deg2rad(v)
        lbl_sp.text = "$(round(Int, v))°"
    end
    on(sl_tf.value) do v
        test_freq[] = v
        lbl_tf.text = "$(round(v, digits=1)) Hz"
    end
    on(sl_nl.value) do v
        noise_lvl[] = v
        lbl_nl.text = "$(round(v, digits=2))"
    end
    on(v -> (use_complex[] = v), tog_complex.active)

    for col = 1:6
        colsize!(ctrl, col, Relative(1 / 6))
    end

    preset_grid = GridLayout(fig[5, 1], tellwidth = false, colgap = 10, padding = (10, 10, 10, 10))
    Label(preset_grid[1, 1], "Presets:", fontsize = ctrl_sz, font = :bold, halign = :right)
    btn_perfect = Button(preset_grid[1, 2], label = "Perfect Match", width = 150)
    btn_leakage = Button(preset_grid[1, 3], label = "Spectral Leakage", width = 150)
    btn_ortho = Button(preset_grid[1, 4], label = "Phase Orthogonality (90°)", width = 220)

    function _apply_preset!(sf_val, tf_val, sp_val, reset_complex::Bool)
        is_applying_preset[] = true
        try
            set_close_to!(sl_sf, sf_val)
            set_close_to!(sl_tf, tf_val)
            set_close_to!(sl_sp, sp_val)
            set_close_to!(sl_nl, 0.0)
            reset_complex && (tog_complex.active[] = false)
        finally
            is_applying_preset[] = false
        end
        update()
    end

    on(_ -> _apply_preset!(10.0, 10.0, 0.0, true), btn_perfect.clicks)
    on(_ -> _apply_preset!(10.0, 10.25, 0.0, true), btn_leakage.clicks)
    on(_ -> _apply_preset!(10.0, 10.0, 90.0, false), btn_ortho.clicks)

    rowsize!(fig.layout, 1, Relative(0.35))
    rowsize!(fig.layout, 2, Relative(0.40))
    rowsize!(fig.layout, 3, Relative(0.0))
    rowsize!(fig.layout, 4, Relative(0.18))
    rowsize!(fig.layout, 5, Relative(0.07))
    ax_prod_cos.blockscene.visible[] = false

    _display_demo_figure(fig)
    return fig, (ax_sig, ax_prod_sin, ax_prod_cos)
end
