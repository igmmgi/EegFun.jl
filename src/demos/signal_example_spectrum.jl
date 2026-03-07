"""
    signal_example_spectrum()

Interactive Power Spectrum Demo — From Time Domain to Frequency Domain.

Shows how a time-domain signal maps to a power spectrum via the FFT, and
demonstrates the key spectral concepts relevant to EEG analysis: frequency
resolution, spectral leakage, and windowing.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1 | Time Domain | The composed signal for the full epoch duration. |
| 2 | Power Spectrum | One-sided FFT power spectrum. Dashed vertical lines mark the true signal frequencies; the orange band marks one frequency-resolution bin (Δf = 1/T). |

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Freq 1 | 1–60 Hz | Frequency of the first sine component |
| Amp 1 | 0–2 | Amplitude of the first sine component |
| Freq 2 | 1–60 Hz | Frequency of the second sine component |
| Amp 2 | 0–2 | Amplitude of the second sine component (0 = off) |
| Noise | 0–2 | Standard deviation of additive Gaussian noise |
| Epoch | 1–10 s | Epoch length — controls frequency resolution (Δf = 1/T) |
| Log scale | toggle | Log y-axis on the spectrum — useful for comparing sidelobes |

## Teaching Notes

**Frequency resolution** is determined by the epoch length: Δf = 1/T.
A 1-second epoch gives 1 Hz bins; a 4-second epoch gives 0.25 Hz bins.

**Spectral leakage** occurs when a signal frequency does not fall exactly on
an FFT bin. Energy "leaks" into neighbouring bins, producing sidelobes.
To see this: set Epoch = 1 s, Freq 1 = 10.25 Hz — the peak is smeared.

The orange band on the spectrum shows the current frequency resolution bin.
Note how it shrinks as you increase the epoch length.

## See Also

- Cohen, M. X. (2014). *Analyzing Neural Time Series Data*. MIT Press. — Chapters 10–11

# Examples
```julia
using EegFun
EegFun.signal_example_spectrum()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_time, ax_freq)` — the two axis objects
"""
function signal_example_spectrum()

    FS = 512.0  # fixed sample rate (Hz) — high enough for up to 60 Hz signals

    fig = Figure(size = (1100, 800), title = "Power Spectrum Demo")

    # ── Font / sizing ───────────────────────────────────────────────────────
    lbl_sz   = Observable(18)
    tick_sz  = Observable(15)
    title_sz = Observable(20)
    ctrl_sz  = Observable(18)

    on(fig.scene.viewport) do area
        sf         = area.widths[1] / 3000
        lbl_sz[]   = max(13, round(Int, 18 * sf))
        tick_sz[]  = max(11, round(Int, 15 * sf))
        title_sz[] = max(15, round(Int, 20 * sf))
        ctrl_sz[]  = max(13, round(Int, 18 * sf))
    end

    # ── Axes ────────────────────────────────────────────────────────────────
    ax_t = Axis(
        fig[1, 1],
        title          = "Time Domain",
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

    ax_f = Axis(
        fig[2, 1],
        title          = "Power Spectrum",
        xlabel         = "Frequency (Hz)",
        ylabel         = "Power",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    # ── Observables ─────────────────────────────────────────────────────────
    freq1     = Observable(10.0)
    amp1      = Observable(1.0)
    freq2     = Observable(25.0)
    amp2      = Observable(0.0)    # amplitude 0 → component 2 silent
    noise_lvl = Observable(0.0)
    epoch_len = Observable(2.0)    # seconds
    log_scale = Observable(false)

    # Data observables
    t_data       = Observable(Float64[])
    sig_data     = Observable(Float64[])
    hz_data      = Observable(Float64[])
    pwr_data     = Observable(Float64[])
    rayleigh_end = Observable(0.5)  # = 1/epoch_len, right edge of resolution band

    # ── Update function ──────────────────────────────────────────────────────
    function update()
        T = epoch_len[]
        t = range(0.0, T - 1.0 / FS, step = 1.0 / FS)
        n = length(t)

        # Compose signal
        sig = amp1[] .* sin.(2π .* freq1[] .* t)
        sig += amp2[] .* sin.(2π .* freq2[] .* t)
        if noise_lvl[] > 0.0
            sig += noise_lvl[] .* randn(n)
        end

        # One-sided power spectrum (amplitude-normalised)
        F = fft(sig)
        half = n ÷ 2 + 1
        hz = collect(range(0.0, FS / 2.0, length = half))
        pwr = @. (abs(F[1:half]) / n)^2 * 2.0
        pwr[1] /= 2.0              # DC: don't double
        iseven(n) && (pwr[end] /= 2.0)  # Nyquist: don't double

        # ── Update observables ──
        t_data[]       = collect(t)
        sig_data[]     = sig
        hz_data[]      = hz
        pwr_data[]     = pwr
        rayleigh_end[] = 1.0 / T

        # ── Axis limits ──
        xlims!(ax_t, 0.0, T)
        ylims!(ax_t, -3.5, 3.5)
        xlims!(ax_f, 0.0, 65.0)

        max_p = max(maximum(pwr), 1e-10)
        if log_scale[]
            # Set strictly-positive limits BEFORE switching to log scale —
            # Makie validates the current limits immediately when yscale changes,
            # so a lower bound of 0.0 would error.
            ylims!(ax_f, max(max_p * 1e-6, 1e-14), max_p * 5.0)
            ax_f.yscale = log10
        else
            # Switch to linear first (safe with any limits), then set limits.
            ax_f.yscale = identity
            ylims!(ax_f, 0.0, max_p * 1.25)
        end
    end

    for obs in (freq1, amp1, freq2, amp2, noise_lvl, epoch_len, log_scale)
        on(obs) do _
            update()
        end
    end
    update()

    # ── Plots ────────────────────────────────────────────────────────────────
    lines!(ax_t, t_data, sig_data, color = :royalblue, linewidth = 2)

    # Frequency resolution band (orange) — one frequency bin wide
    vspan!(ax_f, Observable(0.0), rayleigh_end, color = (:orange, 0.15))

    # Spectrum
    lines!(ax_f, hz_data, pwr_data, color = :royalblue, linewidth = 1.5)

    # True frequency markers
    vlines!(ax_f, @lift([$freq1[]]), color = :crimson, linestyle = :dash, linewidth = 1.5, label = "Freq 1")
    vlines!(ax_f, @lift([$freq2[]]), color = :seagreen, linestyle = :dash, linewidth = 1.5, label = "Freq 2")

    axislegend(ax_f, position = :rt)

    # Dynamic annotation showing current resolution
    res_text = @lift("Δf = $(round(1.0/$epoch_len[], digits=2)) Hz  (T = $(round(Int, $epoch_len[])) s)")
    text!(ax_f, 0.02, 0.02, text = res_text, space = :relative, fontsize = lbl_sz, color = :darkorange, align = (:left, :bottom))

    # ── Controls ─────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[3, 1], colgap = 16)

    function labelled_slider(parent, col, header, range_vals, startval, fmt)
        Label(parent[1, col], header, fontsize = ctrl_sz, halign = :center)
        sl  = Slider(parent[2, col], range = range_vals, startvalue = startval)
        lbl = Label(parent[3, col], fmt(startval), fontsize = ctrl_sz, halign = :center)
        return sl, lbl
    end

    sl_f1, lbl_f1 = labelled_slider(ctrl, 1, "Freq 1", 1.0:0.25:60.0, 10.0, v -> "$(round(v, digits=2)) Hz")
    sl_a1, lbl_a1 = labelled_slider(ctrl, 2, "Amp 1", 0.0:0.1:2.0, 1.0, v -> "$(round(v, digits=1))")
    sl_f2, lbl_f2 = labelled_slider(ctrl, 3, "Freq 2", 1.0:0.25:60.0, 25.0, v -> "$(round(v, digits=2)) Hz")
    sl_a2, lbl_a2 = labelled_slider(ctrl, 4, "Amp 2", 0.0:0.1:2.0, 0.0, v -> "$(round(v, digits=1))")
    sl_nl, lbl_nl = labelled_slider(ctrl, 5, "Noise", 0.0:0.1:2.0, 0.0, v -> "$(round(v, digits=1))")
    sl_ep, lbl_ep = labelled_slider(ctrl, 6, "Epoch", 1.0:1.0:10.0, 2.0, v -> "$(round(Int, v)) s  Δf=$(round(1/v, digits=2)) Hz")

    # Log scale toggle
    Label(ctrl[1, 7], "Log y-axis", fontsize = ctrl_sz, halign = :center)
    tog_log = Toggle(ctrl[2, 7], active = false, halign = :center)
    Label(ctrl[3, 7], "", fontsize = ctrl_sz, halign = :center)

    # ── Wire up controls ─────────────────────────────────────────────────────
    on(sl_f1.value) do v
        freq1[] = v
        lbl_f1.text = "$(round(v, digits=2)) Hz"
    end
    on(sl_a1.value) do v
        amp1[] = v
        lbl_a1.text = "$(round(v, digits=1))"
    end
    on(sl_f2.value) do v
        freq2[] = v
        lbl_f2.text = "$(round(v, digits=2)) Hz"
    end
    on(sl_a2.value) do v
        amp2[] = v
        lbl_a2.text = "$(round(v, digits=1))"
    end
    on(sl_nl.value) do v
        noise_lvl[] = v
        lbl_nl.text = "$(round(v, digits=1))"
    end
    on(sl_ep.value) do v
        epoch_len[] = v
        lbl_ep.text = "$(round(Int, v)) s  Δf=$(round(1/v, digits=2)) Hz"
    end
    on(tog_log.active) do v
        log_scale[] = v
    end

    # ── Row sizing ───────────────────────────────────────────────────────────
    rowsize!(fig.layout, 1, Relative(0.28))   # time domain
    rowsize!(fig.layout, 2, Relative(0.55))   # spectrum
    rowsize!(fig.layout, 3, Relative(0.17))   # controls

    display(fig)
    return fig, (ax_t, ax_f)
end
