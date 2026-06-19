"""
    signal_example_convolution()

Interactive Convolution Demo — Sliding a Kernel Across a Signal.

Demonstrates discrete convolution by sliding a kernel (filter) step by step
across a signal, accumulating the output one position at a time. Three kernel
types are provided: Gaussian (low-pass filter), Boxcar (running average), and
Wavelet (Morlet) — which outputs the amplitude envelope at a specific frequency,
bridging filtering and time-frequency analysis.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1 | Signal + Kernel | The input signal with the kernel overlaid at the current position. The shaded region shows where the kernel overlaps the signal. |
| 2 | Output | The convolution result accumulated up to the current position. For the Wavelet kernel this is the amplitude envelope (always ≥ 0, y-axis fixed 0–1.1). |

## Controls

| Control | Applies to | Range | Description |
|---------|------------|-------|-------------|
| Kernel Position | All | 0–2 s | Scrub the kernel across the signal — output builds up as you drag right |
| Kernel Width | Gaussian / Boxcar | 10–2000 ms | Total span of the kernel (label shows "→ eff. X ms" in Wavelet mode) |
| Cycles | Wavelet | 1–10 | Number of oscillations in the wavelet kernel |
| Wavelet Freq | Wavelet | 1–40 Hz | Target frequency to detect |
| Signal Freq | All | 1–20 Hz | Frequency of the input sine wave |
| Noise | All | 0–1 | Additive Gaussian noise standard deviation |
| Kernel type | All | — | Gaussian / Boxcar / Wavelet (Morlet) |

## Teaching Notes

**Convolution** computes, at each time point *t*, the weighted sum of the
signal under the kernel:

    (x * h)[t] = Σ x[τ] · h[t − τ]

This is how FIR filters work. The kernel *h* is the filter's impulse response:
- **Gaussian** → smooth low-pass filter (removes noise, preserves slow waves)
- **Boxcar** → running average (same idea, sharper edges in frequency domain)
- **Wavelet (Morlet)** → Gaussian × cosine at *f* Hz. Frequency-selective:
  output is the amplitude envelope "how much of *f* Hz is present right now?"
  This is one row of a time-frequency spectrogram.

The **Kernel Width** label shows "→ eff. X ms" in Wavelet mode — this is the
equivalent Gaussian width, determined by `cycles / wavelet_freq`. Setting the
Gaussian Kernel Width slider to that value produces the same Gaussian envelope.

## See Also

- Cohen, M. X. (2014). *Analyzing Neural Time Series Data*. MIT Press. — Chapter 11/12

# Examples
```julia
using EegFun
EegFun.signal_example_convolution()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_sig, ax_out)` — the two axis objects
"""
function signal_example_convolution()

    FS    = 512.0   # fixed sample rate (Hz)
    EPOCH = 2.0     # epoch length (s)

    fig = Figure(size = (1300, 650), figure_padding = (8, 8, 8, 8), title = "Convolution Demo")

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
        title          = "Signal + Kernel at current position",
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
    ax_out = Axis(
        fig[2, 1],
        title          = "Convolution Output",
        xlabel         = "Time (s)",
        ylabel         = "Output",
        titlesize      = title_sz,
        xlabelsize     = lbl_sz,
        ylabelsize     = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible   = true,
        ygridvisible   = true,
    )

    # Link x-axes so zoom/pan stays in sync
    linkxaxes!(ax_sig, ax_out)

    # ── Observables ──────────────────────────────────────────────────────────
    kernel_pos     = Observable(0.0)           # fractional position [0,1]
    kernel_width   = Observable(0.1)           # in seconds (Gaussian / Boxcar)
    kernel_type    = Observable(:gaussian)     # :gaussian, :boxcar, or :wavelet
    wavelet_freq   = Observable(8.0)           # Hz
    wavelet_cycles = Observable(5.0)           # number of cycles (Wavelet mode)
    sig_freq       = Observable(4.0)
    noise_lvl      = Observable(0.3)

    # Kernel types for the menu
    kernel_options = ["Gaussian", "Boxcar", "Wavelet (Morlet)"]

    # Time vector (computed once, fixed)
    t_vec = range(0.0, EPOCH - 1.0 / FS, step = 1.0 / FS)
    t_arr = collect(t_vec)
    n     = length(t_arr)

    # Fixed base noise 
    base_noise = randn(n)

    # Data observables
    sig_data  = Observable(Float64[])
    kern_t    = Observable(Float64[])    # kernel x positions (for overlay plot)
    kern_vals = Observable(Float64[])    # kernel y values (rescaled for overlay)
    shd_x     = Observable(Float64[])    # shaded region x
    shd_y     = Observable(Float64[])    # shaded region y (signal under kernel)
    out_t     = Observable(Float64[])    # output time points up to current pos
    out_data  = Observable(Float64[])    # accumulated output

    # Returns (k_real, k_imag, τ). k_imag is nothing for non-wavelet types.
    # width_s = TOTAL span (edge to edge): half = width_s/2, σ = width_s/6 → ±3σ = ±half.
    """Construct a Gaussian, Boxcar, or Morlet wavelet kernel and return `(k_real, k_imag, τ)`."""
    function build_kernel(type, width_s, wfreq, ncycles)
        if type === :wavelet
            width_s = ncycles / wfreq   # total span = ncycles/freq → exactly ncycles visible
        end
        half = round(Int, width_s / 2.0 * FS)
        σ = width_s / 6.0            # ±3σ = ±half so Gaussian tapers to ~1% at edges
        τ = ((-half):half) ./ FS

        if type === :gaussian
            k = exp.(-(τ .^ 2) ./ (2 * σ^2))
            k ./= sum(k)
            return k, nothing, τ
        elseif type === :boxcar
            k = ones(length(τ))
            k ./= sum(k)
            return k, nothing, τ
        else  # :wavelet — Morlet = Gaussian × complex sinusoid
            gauss    = exp.(-(τ .^ 2) ./ (2 * σ^2))
            k_r      = gauss .* cos.(2π .* wfreq .* τ)   # cosine — peaks at centre
            k_i      = gauss .* sin.(2π .* wfreq .* τ)   # sine — zero at centre
            norm_fac = sum(gauss) / 2.0
            return k_r ./ norm_fac, k_i ./ norm_fac, τ
        end
    end

    # ── Full convolution helper ───────────────────────────────────────────────
    """Linear convolution with same-length (centred) output."""
    function full_conv(sig, k)
        # Simple linear convolution, same-length output (centred)
        out = zeros(length(sig))
        h   = (length(k) - 1) ÷ 2
        for i in eachindex(sig)
            acc = 0.0
            for j in eachindex(k)
                idx = i - h + j - 1
                if 1 <= idx <= length(sig)
                    acc += sig[idx] * k[end-j+1]
                end
            end
            out[i] = acc
        end
        return out
    end

    # ── Update function ──────────────────────────────────────────────────────
    """Recompute signal, kernel overlay, and convolution output from current slider values."""
    function update()
        # Build signal
        sig = sin.(2π .* sig_freq[] .* t_vec)
        if noise_lvl[] > 0.0
            sig .+= noise_lvl[] .* base_noise
        end

        # Build kernel
        k_r, k_i, τ = build_kernel(kernel_type[], kernel_width[], wavelet_freq[], wavelet_cycles[])

        # Current kernel centre position in samples
        pos_s   = kernel_pos[] * EPOCH
        pos_idx = clamp(round(Int, pos_s * FS), 1, n)

        # Compute output: amplitude envelope for wavelet, direct for others
        if k_i !== nothing
            out_r    = full_conv(collect(sig), k_r)
            out_i    = full_conv(collect(sig), k_i)
            full_out = sqrt.(out_r .^ 2 .+ out_i .^ 2)
        else
            full_out = full_conv(collect(sig), k_r)
        end

        # Kernel overlay — show real part, rescaled
        kern_centre = t_arr[pos_idx]
        k_t = kern_centre .+ collect(τ)
        k_v = k_r .* (1.0 / max(maximum(abs.(k_r)), 1e-8))

        # Shaded overlap region
        half = (length(k_r) - 1) ÷ 2
        i_lo = max(1, pos_idx - half)
        i_hi = min(n, pos_idx + half)

        # Output up to current position
        out_indices = 1:pos_idx

        # Update observables
        sig_data[]  = collect(sig)
        kern_t[]    = k_t
        kern_vals[] = collect(k_v)
        shd_x[]     = t_arr[i_lo:i_hi]
        shd_y[]     = collect(sig[i_lo:i_hi])
        out_t[]     = t_arr[out_indices]
        out_data[]  = full_out[out_indices]

        # Update output panel title
        ax_out.title =
            kernel_type[] === :wavelet ? "Amplitude Envelope — |Wavelet convolution at $(round(Int, wavelet_freq[])) Hz|" :
            "Convolution Output (built up as kernel slides)"

        xlims!(ax_sig, 0.0, EPOCH)
        xlims!(ax_out, 0.0, EPOCH)
        lim = 1.5 + noise_lvl[] * 1.5
        ylims!(ax_sig, -lim, lim)
        max_o = max(maximum(abs.(full_out)), 0.01)
        if kernel_type[] === :wavelet
            ylims!(ax_out, 0.0, 1.1)   # signal amplitude fixed at 1 → envelope ≤ 1
        else
            ylims!(ax_out, -max_o * 1.2, max_o * 1.2)
        end
    end

    for obs in (kernel_pos, kernel_width, kernel_type, wavelet_freq, wavelet_cycles, sig_freq, noise_lvl)
        on(obs) do _
            update()
        end
    end
    update()

    # ── Plots ────────────────────────────────────────────────────────────────
    # Signal panel
    lines!(ax_sig, t_arr, sig_data, color = :black, linewidth = 2)

    # Shaded kernel footprint — band from axis bottom to signal value, showing signal samples under kernel
    band!(ax_sig, shd_x, @lift(fill(-3.0, length($shd_x))), @lift(fill(3.0, length($shd_x))), color = (:black, 0.07))

    # Kernel shape drawn as a filled + outline at the current position
    # Rescale kernel to amplitude 0.9 so it sits clearly inside ±1 signal
    band!(ax_sig, kern_t, @lift(zeros(length($kern_t))), kern_vals, color = (:black, 0.20))
    lines!(ax_sig, kern_t, kern_vals, color = :black, linewidth = 2, label = "Kernel")

    # Current position marker
    vlines!(ax_sig, @lift([$kernel_pos[] * EPOCH]), color = :black, linestyle = :dash, linewidth = 1.5)

    axislegend(ax_sig, position = :rt)

    # Output panel
    lines!(ax_out, out_t, out_data, color = :black, linewidth = 2)
    hlines!(ax_out, [0.0], color = :black, linewidth = 0.8, linestyle = :dash)
    vlines!(ax_out, @lift([$kernel_pos[] * EPOCH]), color = :black, linestyle = :dash, linewidth = 1.5)

    # Text: current output value at the kernel position
    out_text = @lift begin
        pos_s = $kernel_pos * EPOCH
        "Output at t = $(round(pos_s, digits=2)) s = $(isempty($out_data) ? "–" : round($out_data[end], digits=3))"
    end
    text!(ax_out, 0.5, 0.97, text = out_text, space = :relative, fontsize = lbl_sz, color = :black, align = (:center, :top))

    # ── Controls ─────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[3, 1], colgap = 16)

    # Kernel position slider (full width)
    Label(ctrl[1, 1:4], "Kernel Position", fontsize = ctrl_sz, halign = :center)
    sl_pos = Slider(ctrl[2, 1:4], range = 0.0:0.001:1.0, startvalue = 0.0)
    lbl_pos = Label(ctrl[3, 1:4], "pos: 0.00 s", fontsize = ctrl_sz, halign = :center)

    inner = GridLayout(ctrl[4, 1:4], colgap = 16)

    """Create a slider with a header label above and a value label below."""
    function labelled_slider(parent, col, header, range_vals, startval, fmt)
        Label(parent[1, col], header, fontsize = ctrl_sz, halign = :center)
        sl  = Slider(parent[2, col], range = range_vals, startvalue = startval)
        lbl = Label(parent[3, col], fmt(startval), fontsize = ctrl_sz, halign = :center)
        return sl, lbl
    end

    sl_kw, lbl_kw = labelled_slider(inner, 1, "Kernel Width", 0.01:0.01:2.0, 0.1, v -> "$(round(Int, v * 1000)) ms")
    sl_cy, lbl_cy = labelled_slider(inner, 2, "Cycles", 1.0:1.0:10.0, 5.0, v -> "$(round(Int, v))")
    sl_wf, lbl_wf = labelled_slider(inner, 3, "Wavelet Freq", 1.0:1.0:40.0, 8.0, v -> "$(round(Int, v)) Hz")
    sl_sf, lbl_sf = labelled_slider(inner, 4, "Signal Freq", 1.0:1.0:20.0, 4.0, v -> "$(round(Int, v)) Hz")
    sl_nl, lbl_nl = labelled_slider(inner, 5, "Noise", 0.0:0.05:1.0, 0.3, v -> "$(round(v, digits=2))")

    # Kernel type menu
    Label(inner[1, 6], "Kernel type", fontsize = ctrl_sz, halign = :center)
    menu_kt = Menu(inner[2, 6], options = kernel_options, default = "Gaussian", fontsize = ctrl_sz)
    Label(inner[3, 6], "", fontsize = ctrl_sz, halign = :center)

    # Helper: keep the Kernel Width label linked to effective width in all modes
    """Update the kernel-width label to show effective width in wavelet mode."""
    function update_kw_label()
        if kernel_type[] === :wavelet
            eff_ms = round(Int, wavelet_cycles[] / wavelet_freq[] * 1000)
            lbl_kw.text = "→ eff. $eff_ms ms"
        else
            lbl_kw.text = "$(round(Int, kernel_width[] * 1000)) ms"
        end
    end

    # Wire up
    on(sl_pos.value) do v
        kernel_pos[] = v
        lbl_pos.text = "pos: $(round(v * EPOCH, digits=2)) s"
    end
    on(sl_kw.value) do v
        kernel_width[] = v
        update_kw_label()
    end
    on(sl_cy.value) do v
        wavelet_cycles[] = v
        lbl_cy.text = "$(round(Int, v))"
        update_kw_label()
    end
    on(sl_wf.value) do v
        wavelet_freq[] = v
        lbl_wf.text = "$(round(Int, v)) Hz"
        update_kw_label()
    end
    on(sl_sf.value) do v
        sig_freq[] = v
        lbl_sf.text = "$(round(Int, v)) Hz"
    end
    on(sl_nl.value) do v
        noise_lvl[] = v
        lbl_nl.text = "$(round(v, digits=2))"
    end
    on(menu_kt.selection) do sel
        kernel_type[] = (sel == "Gaussian") ? :gaussian : (sel == "Boxcar") ? :boxcar : :wavelet
        update_kw_label()
    end

    for col = 1:4
        colsize!(ctrl, col, Relative(1 / 4))
    end
    for col = 1:6
        colsize!(inner, col, Relative(1 / 6))
    end

    # ── Presets ──────────────────────────────────────────────────────────────
    preset_grid = GridLayout(fig[4, 1], tellwidth = false, colgap = 10, padding = (10, 10, 10, 10))
    Label(preset_grid[1, 1], "Presets", fontsize = ctrl_sz, font = :bold, halign = :right)
    btn_play = Button(preset_grid[1, 2], label = "▶ Play Animation", width = 180)
    btn_denoise = Button(preset_grid[1, 3], label = "Smooth Denoise (Gaussian)", width = 240)
    btn_wavelet = Button(preset_grid[1, 4], label = "Amplitude Envelope (Wavelet)", width = 260)

    is_playing = Ref(false)
    on(btn_play.clicks) do _
        is_playing[] && return
        is_playing[] = true
        @async begin
            try
                for p = 0.0:0.01:1.0
                    !is_playing[] && break
                    set_close_to!(sl_pos, p)
                    sleep(0.05)   # 5 seconds total
                end
            finally
                is_playing[] = false
            end
        end
    end

    on(btn_denoise.clicks) do _
        is_playing[] = false
        idx = findfirst(x -> x == "Gaussian", kernel_options)
        if !isnothing(idx)
            menu_kt.i_selected[] = idx
        end
        set_close_to!(sl_kw, 0.15)
        set_close_to!(sl_sf, 4.0)
        set_close_to!(sl_nl, 0.5)
        set_close_to!(sl_pos, 0.0)
    end

    on(btn_wavelet.clicks) do _
        is_playing[] = false
        idx = findfirst(x -> x == "Wavelet (Morlet)", kernel_options)
        if !isnothing(idx)
            menu_kt.i_selected[] = idx
        end
        set_close_to!(sl_wf, 8.0)
        set_close_to!(sl_sf, 8.0)
        set_close_to!(sl_nl, 0.0)
        set_close_to!(sl_pos, 0.0)
    end

    # ── Row sizing ────────────────────────────────────────────────────────────
    rowsize!(fig.layout, 1, Relative(0.38))   # signal + kernel
    rowsize!(fig.layout, 2, Relative(0.38))   # output
    rowsize!(fig.layout, 3, Relative(0.18))   # controls
    rowsize!(fig.layout, 4, Relative(0.06))   # presets

    display(fig)
    return fig, (ax_sig, ax_out)
end
