"""
    signal_example_dotproduct()

Interactive Dot Product Demo — How the DFT Detects Frequencies.

Shows the single core idea behind the **Discrete Fourier Transform (DFT)**:
multiplying a signal by a test sinusoid and summing the result. When the test
frequency matches the signal frequency, the product is all-positive and the
sum (dot product) is large. When they do not match, the product alternates
and the sum is near zero.

A **Complex** toggle adds a cosine test sinusoid alongside the sine, showing
how the DFT achieves phase-independence by computing the magnitude from both
projections.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1 | Signal + Test sinusoid(s) | The target signal (blue), sine test (red), and cosine test (green, when Complex is on). |
| 2 | Sine product | Signal × sine test sinusoid, with dot product value. |
| 3 | Cosine product | Signal × cosine test sinusoid (visible only when Complex is on). |

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Signal Freq | 1–40 Hz | Frequency of the signal sine wave |
| Amplitude | 0.1–2 | Signal amplitude — dot product scales proportionally |
| Phase (°) | 0–355° | Signal phase offset |
| Test Freq | 0.5–60 Hz | Frequency of the test sinusoid (the DFT "probe") |
| Noise | 0–1 | Noise standard deviation |
| Complex | toggle | Switch from sine-only to complex (sine + cosine) probing |

## Teaching Notes

The **dot product** measures how similar two signals are — multiply them
point-by-point and sum the results:

    dot product = Σ a[t] × b[t]

The DFT uses exactly this idea to answer "how much of frequency *f* is in
this signal?" by setting one of the signals to a test sinusoid:

    dot product = Σ signal[t] × sin(2πft)

**Freq match, Phase match:** product is all-positive → large sum.

**Freq mismatch:** product alternates → cancels to near zero.

**Phase mismatch (90°):** product alternates even at matching frequency →
dot product ≈ 0.

**The complex solution:** the sine-only dot product is sensitive to phase.
The DFT solves this by computing *two* dot products — one with a sine and
one with a cosine — and combining them:

    magnitude = √(sine_dp² + cosine_dp²)

The magnitude is phase-independent: it captures the signal's power at
frequency *f* regardless of the signal's phase offset. Toggle **Complex**
on and drag Phase to see this in action.

**Amplitude:** the test sinusoid is always ±1; the dot product scales directly
with the signal amplitude.

The FFT runs this operation for every frequency simultaneously. This demo
shows just one bin at a time.

## See Also

- Cohen, M. X. (2014). *Analyzing Neural Time Series Data*. MIT Press. — Chapter 11

# Examples
```julia
using EegFun
EegFun.signal_example_dotproduct()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_sig, ax_prod_sin, ax_prod_cos)` — the three axis objects
"""
function signal_example_dotproduct()

    FS    = 1000.0  # sample rate (Hz) — high enough for smooth rendering up to 40 Hz
    EPOCH = 2.0     # epoch length (s)

    fig = Figure(size = (1300, 700), figure_padding = (8, 8, 8, 8), title = "Dot Product Demo")

    # ── Font sizing ──────────────────────────────────────────────────────────
    lbl_sz   = Observable(17)
    tick_sz  = Observable(14)
    title_sz = Observable(19)
    ctrl_sz  = Observable(17)

    on(fig.scene.viewport) do area
        sf         = area.widths[1] / 3000
        lbl_sz[]   = max(12, round(Int, 17 * sf))
        tick_sz[]  = max(10, round(Int, 14 * sf))
        title_sz[] = max(14, round(Int, 19 * sf))
        ctrl_sz[]  = max(12, round(Int, 17 * sf))
    end

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax_sig = Axis(
        fig[1, 1],
        title          = "Signal (blue)  +  Test sinusoid (red)",
        xlabel         = "Time (s)",
        ylabel         = "Amplitude",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    ax_prod_sin = Axis(
        fig[2, 1],
        title          = "Sine product  (signal × sin)",
        xlabel         = "Time (s)",
        ylabel         = "",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    ax_prod_cos = Axis(
        fig[3, 1],
        title          = "Cosine product  (signal × cos)",
        xlabel         = "Time (s)",
        ylabel         = "",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    # ── Observables ──────────────────────────────────────────────────────────
    sig_freq    = Observable(10.0)
    sig_amp     = Observable(1.0)
    sig_phase   = Observable(0.0)   # radians
    noise_lvl   = Observable(0.0)
    test_freq   = Observable(10.0)
    use_complex = Observable(false)

    t_arr = collect(range(0.0, EPOCH - 1.0 / FS, step = 1.0 / FS))
    n     = length(t_arr)

    sig_data      = Observable(zeros(n))
    test_sin_data = Observable(zeros(n))
    test_cos_data = Observable(zeros(n))
    prod_sin_data = Observable(zeros(n))
    prod_cos_data = Observable(zeros(n))

    # ── Update ───────────────────────────────────────────────────────────────
    """Recompute signal, test sinusoids, and dot products from current parameter values."""
    function update()
        sig = sig_amp[] .* sin.(2π .* sig_freq[] .* t_arr .+ sig_phase[])
        if noise_lvl[] > 0.0
            sig .+= noise_lvl[] .* randn(n)
        end
        test_sin = sin.(2π .* test_freq[] .* t_arr)
        test_cos = cos.(2π .* test_freq[] .* t_arr)
        prod_sin = sig .* test_sin
        prod_cos = sig .* test_cos

        sig_data[]      = sig
        test_sin_data[] = test_sin
        test_cos_data[] = test_cos
        prod_sin_data[] = prod_sin
        prod_cos_data[] = prod_cos

        xlims!(ax_sig, 0.0, EPOCH)
        xlims!(ax_prod_sin, 0.0, EPOCH)
        xlims!(ax_prod_cos, 0.0, EPOCH)
        a = sig_amp[] + noise_lvl[]
        ylims!(ax_sig, -a * 1.1, a * 1.1)

        # Scale product axes to actual data range
        for (ax, pdata) in ((ax_prod_sin, prod_sin), (ax_prod_cos, prod_cos))
            p_min = min(minimum(pdata), -0.05)
            p_max = max(maximum(pdata), 0.05)
            p_margin = (p_max - p_min) * 0.12
            ylims!(ax, p_min - p_margin, p_max + p_margin)
        end
    end

    for obs in (sig_freq, sig_amp, sig_phase, noise_lvl, test_freq)
        on(obs) do _
            update()
        end
    end
    update()

    # ── Visibility management ────────────────────────────────────────────────
    # Hide/show cosine row and update titles based on Complex toggle
    on(use_complex) do is_complex
        ax_prod_cos.blockscene.visible[] = is_complex

        if is_complex
            rowsize!(fig.layout, 1, Relative(0.28))
            rowsize!(fig.layout, 2, Relative(0.28))
            rowsize!(fig.layout, 3, Relative(0.28))
            rowsize!(fig.layout, 4, Relative(0.16))
            ax_sig.title = "Signal (blue)  +  Sine test (red)  +  Cosine test (green)"
        else
            rowsize!(fig.layout, 1, Relative(0.38))
            rowsize!(fig.layout, 2, Relative(0.44))
            rowsize!(fig.layout, 3, Relative(0.0))
            rowsize!(fig.layout, 4, Relative(0.18))
            ax_sig.title = "Signal (blue)  +  Test sinusoid (red)"
        end
    end

    # ── Plots ─────────────────────────────────────────────────────────────────
    lines!(ax_sig, t_arr, sig_data, color = :royalblue, linewidth = 2, label = "Signal")
    lines!(ax_sig, t_arr, test_sin_data, color = :crimson, linewidth = 2, label = "Sine test")
    cos_line = lines!(ax_sig, t_arr, test_cos_data, color = :seagreen, linewidth = 2, label = "Cosine test", visible = false)
    axislegend(ax_sig, position = :rt)

    # Sine product panel
    lines!(ax_prod_sin, t_arr, prod_sin_data, color = :black, linewidth = 2)
    hlines!(ax_prod_sin, [0.0], color = :black, linewidth = 1.0, linestyle = :dash)

    # Cosine product panel
    lines!(ax_prod_cos, t_arr, prod_cos_data, color = :black, linewidth = 2)
    hlines!(ax_prod_cos, [0.0], color = :black, linewidth = 1.0, linestyle = :dash)

    # Dot product annotation — sine panel
    dp_sin_text = @lift begin
        dp          = sum($prod_sin_data) / n * 2.0
        freq_match  = abs($sig_freq - $test_freq) < 0.5
        phase_deg   = mod(rad2deg($sig_phase), 360.0)
        phase_match = phase_deg < 10.0 || phase_deg > 350.0
        freq_str    = freq_match ? "Freq: ✓ Match" : "Freq: ✗ Mismatch"
        phase_str   = phase_match ? "Phase: ✓ Match" : "Phase: ✗ Mismatch"
        if $use_complex
            "Sine dp = $(round(dp, digits=3))     $freq_str     $phase_str"
        else
            "Dot product = $(round(dp, digits=3))     $freq_str     $phase_str"
        end
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

    # Dot product annotation — cosine panel
    dp_cos_text = @lift begin
        dp = sum($prod_cos_data) / n * 2.0
        "Cosine dp = $(round(dp, digits=3))"
    end
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

    # Magnitude annotation — shown on sine panel when complex is active
    mag_text = @lift begin
        if $use_complex
            dp_sin = sum($prod_sin_data) / n * 2.0
            dp_cos = sum($prod_cos_data) / n * 2.0
            mag    = sqrt(dp_sin^2 + dp_cos^2)
            "Magnitude = √(sin² + cos²) = $(round(mag, digits=3))"
        else
            ""
        end
    end
    text!(
        ax_prod_sin,
        0.5,
        0.03,
        text = mag_text,
        space = :relative,
        fontsize = @lift(max(14, round(Int, $lbl_sz * 1.25))),
        color = :purple,
        font = :bold,
        align = (:center, :bottom),
    )

    # Toggle cosine test line visibility
    on(use_complex) do is_complex
        cos_line.visible = is_complex
        update()  # recompute to trigger annotation refresh
    end

    # ── Controls ──────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[4, 1], colgap = 20)

    """Create a slider with a header label above and a value label below."""
    function labelled_slider(parent, col, header, range_vals, startval, fmt)
        Label(parent[1, col], header, fontsize = ctrl_sz, halign = :center)
        sl  = Slider(parent[2, col], range = range_vals, startvalue = startval)
        lbl = Label(parent[3, col], fmt(startval), fontsize = ctrl_sz, halign = :center)
        return sl, lbl
    end

    sl_sf, lbl_sf = labelled_slider(ctrl, 1, "Signal Freq", 1.0:1.0:40.0, 10.0, v -> "$(round(Int, v)) Hz")
    sl_sa, lbl_sa = labelled_slider(ctrl, 2, "Amplitude", 0.1:0.1:2.0, 1.0, v -> "$(round(v, digits=1))")
    sl_sp, lbl_sp = labelled_slider(ctrl, 3, "Phase (°)", 0.0:5.0:355.0, 0.0, v -> "$(round(Int, v))°")
    sl_tf, lbl_tf = labelled_slider(ctrl, 4, "Test Freq", 0.5:0.5:60.0, 10.0, v -> "$(round(v, digits=1)) Hz")
    sl_nl, lbl_nl = labelled_slider(ctrl, 5, "Noise", 0.0:0.05:1.0, 0.0, v -> "$(round(v, digits=2))")

    # Complex toggle
    Label(ctrl[1, 6], "Complex", fontsize = ctrl_sz, halign = :center)
    tog_complex = Toggle(ctrl[2, 6], active = false, halign = :center)
    Label(ctrl[3, 6], "", fontsize = ctrl_sz, halign = :center)

    on(sl_sf.value) do v
        sig_freq[] = v
        lbl_sf.text = "$(round(Int, v)) Hz"
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
    on(tog_complex.active) do v
        use_complex[] = v
    end

    # Force control columns to fill the full figure width
    for col = 1:6
        colsize!(ctrl, col, Relative(1 / 6))
    end

    # ── Row sizing — initial (sine-only, cosine hidden) ───────────────────────
    rowsize!(fig.layout, 1, Relative(0.38))
    rowsize!(fig.layout, 2, Relative(0.44))
    rowsize!(fig.layout, 3, Relative(0.0))
    rowsize!(fig.layout, 4, Relative(0.18))

    # Hide cosine row initially
    ax_prod_cos.blockscene.visible[] = false

    display(fig)
    return fig, (ax_sig, ax_prod_sin, ax_prod_cos)
end
