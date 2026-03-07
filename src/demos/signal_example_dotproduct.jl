"""
    signal_example_dotproduct()

Interactive Dot Product Demo — How the DFT Detects Frequencies.

Shows the single core idea: multiplying a signal by a test sinusoid and
summing the result. When the test frequency matches the signal frequency,
the product is all-positive and the sum (dot product) is large. When
they do not match, the product alternates and the sum is near zero.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1 | Signal + Test sinusoid | The target signal (blue) and the test sinusoid (red) overlaid. |
| 2 | Element-wise product | Signal × test sinusoid, with the dot product value shown. |

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Signal Freq | 1–40 Hz | Frequency of the signal sine wave |
| Amplitude | 0.1–2 | Signal amplitude — dot product scales proportionally |
| Phase (°) | 0–355° | Signal phase offset |
| Test Freq | 0.5–60 Hz | Frequency of the test sinusoid (the DFT "probe") |
| Noise | 0–1 | Noise standard deviation |

## Teaching Notes

The **Discrete Fourier Transform** answers "how much of frequency *f* is in
this signal?" using a single operation:

    dot product = Σ signal[t] × sin(2πft)

**Freq match, Phase match:** product is all-positive → large sum.

**Freq mismatch:** product alternates → cancels to near zero.

**Phase mismatch (90°):** product alternates even at matching frequency →
dot product ≈ 0. This is why the real DFT uses a *complex* sinusoid (sine +
cosine) to be phase-independent.

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
- `axes::Tuple`: `(ax_sig, ax_prod)` — the two axis objects
"""
function signal_example_dotproduct()

    FS    = 1000.0  # sample rate (Hz) — high enough for smooth rendering up to 40 Hz
    EPOCH = 2.0     # epoch length (s)

    fig = Figure(size = (1300, 600), figure_padding = (8, 8, 8, 8), title = "Dot Product Demo")

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

    ax_prod = Axis(
        fig[2, 1],
        title          = "Element-wise product  (signal × test sinusoid)",
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
    sig_freq  = Observable(10.0)
    sig_amp   = Observable(1.0)
    sig_phase = Observable(0.0)   # radians
    noise_lvl = Observable(0.0)
    test_freq = Observable(10.0)

    t_arr = collect(range(0.0, EPOCH - 1.0 / FS, step = 1.0 / FS))
    n     = length(t_arr)

    sig_data  = Observable(zeros(n))
    test_data = Observable(zeros(n))
    prod_data = Observable(zeros(n))

    # ── Update ───────────────────────────────────────────────────────────────
    function update()
        sig = sig_amp[] .* sin.(2π .* sig_freq[] .* t_arr .+ sig_phase[])
        if noise_lvl[] > 0.0
            sig .+= noise_lvl[] .* randn(n)
        end
        test = sin.(2π .* test_freq[] .* t_arr)
        prod = sig .* test

        sig_data[]  = sig
        test_data[] = test
        prod_data[] = prod

        xlims!(ax_sig, 0.0, EPOCH)
        xlims!(ax_prod, 0.0, EPOCH)
        a = sig_amp[] + noise_lvl[]
        ylims!(ax_sig, -a * 1.1, a * 1.1)

        # Scale product axis to actual data range — avoids empty space when
        # the product is all-positive (frequencies match) or all-negative.
        p_min = min(minimum(prod), -0.05)  # always show a sliver of negative axis
        p_max = max(maximum(prod), 0.05)
        p_margin = (p_max - p_min) * 0.12
        ylims!(ax_prod, p_min - p_margin, p_max + p_margin)
    end

    for obs in (sig_freq, sig_amp, sig_phase, noise_lvl, test_freq)
        on(obs) do _
            update()
        end
    end
    update()

    # ── Plots ─────────────────────────────────────────────────────────────────
    lines!(ax_sig, t_arr, sig_data, color = :royalblue, linewidth = 2, label = "Signal")
    lines!(ax_sig, t_arr, test_data, color = :crimson, linewidth = 2, label = "Test")
    axislegend(ax_sig, position = :rt)

    # Product panel — simple line
    lines!(ax_prod, t_arr, prod_data, color = :black, linewidth = 2)
    hlines!(ax_prod, [0.0], color = :black, linewidth = 1.0, linestyle = :dash)

    # Dot product value — large, centred
    dp_text = @lift begin
        dp          = sum($prod_data) / n * 2.0
        freq_match  = abs($sig_freq - $test_freq) < 0.5
        phase_deg   = mod(rad2deg($sig_phase), 360.0)
        phase_match = phase_deg < 10.0 || phase_deg > 350.0
        freq_str    = freq_match ? "Freq: ✓ Match" : "Freq: ✗ Mismatch"
        phase_str   = phase_match ? "Phase: ✓ Match" : "Phase: ✗ Mismatch"
        "Dot product = $(round(dp, digits=3))     $freq_str     $phase_str"
    end
    text!(
        ax_prod,
        0.5,
        0.97,
        text = dp_text,
        space = :relative,
        fontsize = @lift(max(14, round(Int, $lbl_sz * 1.15))),
        color = :black,
        align = (:center, :top),
    )

    # ── Controls ──────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[3, 1], colgap = 20)

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

    # Force control columns to fill the full figure width (avoids the sub-layout
    # shrinking to slider minimum-width and dragging the plot columns with it).
    for col = 1:5
        colsize!(ctrl, col, Relative(1 / 5))
    end

    # ── Row sizing ────────────────────────────────────────────────────────────
    rowsize!(fig.layout, 1, Relative(0.40))
    rowsize!(fig.layout, 2, Relative(0.46))
    rowsize!(fig.layout, 3, Relative(0.14))

    display(fig)
    return fig, (ax_sig, ax_prod)
end
