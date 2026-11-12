"""
    signal_example_eda_scr()

Interactive Electrodermal Activity (EDA) simulator illustrating eccrine sweat gland biophysics
(Edelberg's poral valve model), the parallel skin conductance network (\$G_{\\text{total}} = \\sum G_i = 1 / R_{\\text{total}}\$),
conductance (\\mu S) vs. resistance (k\\Omega) non-linearity, and standard event-related
Skin Conductance Response (ER-SCR) morphology and measures (Latency, Rise Time, Amplitude, and 50% Half-Recovery Time).

```julia
using GLMakie, EegFun
EegFun.signal_example_eda_scr()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax_anatomy::Axis`: Eccrine sweat ducts, poral valve & parallel conductance circuit schematic axis
- `ax_scr::Axis`: Single-trial ER-SCR waveform and morphological measurement axis
"""
function signal_example_eda_scr()
    _check_interactive_demo_backend("Interactive EDA SCR demo") || return nothing
    _set_window_title("Interactive EDA Simulator — Sweat Gland Biophysics & SCR Measures")

    fig = _eegfun_figure(
        size = (1420, 880),
        title = "Interactive EDA Simulator — Sweat Gland Biophysics & SCR Measures",
        backgroundcolor = :white,
    )
    title_font, label_font, tick_font, btn_font = _demo_adaptive_fonts(fig, 1420, (16.5, 13.0), (13.0, 10.5), (11.0, 9.5), (12.5, 10.5))

    scl_base = Observable(10.0)
    scr_amp = Observable(0.5)
    scr_lat = Observable(1.5)
    scr_rise = Observable(2.0)
    scr_rec_half = Observable(4.0)
    playback_time = Observable(-1.0)

    t_stim = 2.0
    t_max = 18.0
    t_vec = collect(range(0.0, t_max, length = 720))
    a_half = 1.6783469900166605

    function eval_scr_unit(t::Real, t_on::Real, t_r::Real, t_h::Real)
        u = t - t_on
        if u <= 0.0
            return 0.0
        elseif u <= t_r
            return 0.5 * (1.0 - cos(π * u / t_r))
        else
            x = (u - t_r) / t_h
            return (1.0 + a_half * x) * exp(-a_half * x)
        end
    end

    # =========================================================================
    # LEFT PANEL: Eccrine Sweat Glands, Poral Valve & Parallel Conductance
    # =========================================================================
    ax_anatomy =
        Axis(fig[1, 1], title = "Eccrine Sweat Glands (Poral Valve) & Parallel Conductance", titlesize = title_font, aspect = DataAspect())
    hidedecorations!(ax_anatomy)
    hidespines!(ax_anatomy)
    xlims!(ax_anatomy, -5.2, 5.2)
    ylims!(ax_anatomy, -3.7, 3.3)

    # Skin layers: Epidermis (top) and Dermis (middle/bottom)
    poly!(ax_anatomy, Rect2f(-4.8, 1.15, 9.6, 0.85), color = "#fde68a", strokecolor = "#d97706", strokewidth = 1.5)
    poly!(ax_anatomy, Rect2f(-4.8, -1.95, 9.6, 3.10), color = "#fecdd3", strokecolor = "#e11d48", strokewidth = 1.2)

    text!(ax_anatomy, 4.65, 1.58, text = "Epidermis\n(Corneum)", align = (:right, :center), font = :bold, fontsize = 12, color = "#92400e")
    text!(ax_anatomy, 4.65, -0.35, text = "Dermis", align = (:right, :center), font = :bold, fontsize = 13, color = "#9f1239")

    # Exosomatic Constant-Voltage (0.5 V) Surface Ag/AgCl Electrodes
    poly!(ax_anatomy, Rect2f(-4.45, 2.00, 1.45, 0.34), color = "#334155", strokecolor = "#0f172a", strokewidth = 1.5)
    poly!(ax_anatomy, Rect2f(3.00, 2.00, 1.45, 0.34), color = "#334155", strokecolor = "#0f172a", strokewidth = 1.5)
    text!(ax_anatomy, -3.72, 2.17, text = "Ag/AgCl (+)", align = (:center, :center), font = :bold, fontsize = 11.0, color = :white)
    text!(ax_anatomy, 3.72, 2.17, text = "Ag/AgCl (–)", align = (:center, :center), font = :bold, fontsize = 11.0, color = :white)

    # Top circuit loop connecting electrodes to Constant Voltage source U = 0.5V & Ammeter (G = I/U)
    lines!(ax_anatomy, [Point2f(-3.72, 2.34), Point2f(-3.72, 2.84), Point2f(-1.90, 2.84)], color = "#0f172a", linewidth = 2.0)
    lines!(ax_anatomy, [Point2f(1.90, 2.84), Point2f(3.72, 2.84), Point2f(3.72, 2.34)], color = "#0f172a", linewidth = 2.0)
    poly!(ax_anatomy, Rect2f(-1.90, 2.58, 3.80, 0.52), color = "#f8fafc", strokecolor = "#0284c7", strokewidth = 1.8)
    text!(
        ax_anatomy,
        0.0,
        2.84,
        text = "Exosomatic U = 0.5 V  (G = I / U)",
        align = (:center, :center),
        font = :bold,
        fontsize = 12,
        color = "#0369a1",
    )

    # Sympathetic cholinergic nerve trunk at bottom of dermis
    lines!(ax_anatomy, [Point2f(-3.8, -1.52), Point2f(3.0, -1.52)], color = "#047857", linewidth = 3.0)
    text!(
        ax_anatomy,
        -0.4,
        -1.76,
        text = "Sympathetic ANS Fibers (Acetylcholine)",
        align = (:center, :center),
        font = :bold,
        fontsize = 11.5,
        color = "#065f46",
    )

    duct_xs = [-3.0f0, -1.2f0, 0.6f0, 2.4f0]
    unit_circ = [Point2f(cos(t), sin(t)) for t in range(0, 2π, length = 36)]

    # Draw static gland secretory coils, duct walls, and parallel resistor symbols
    for (idx, dx) in enumerate(duct_xs)
        # Nerve branch to gland
        lines!(ax_anatomy, [Point2f(dx, -1.52), Point2f(dx, -1.15)], color = "#047857", linewidth = 1.8)
        # Secretory coil in deep dermis
        poly!(
            ax_anatomy,
            [Point2f(dx, -0.95) + 0.36f0 * p for p in unit_circ],
            color = "#fda4af",
            strokecolor = "#be123c",
            strokewidth = 1.6,
        )
        # Outer duct walls rising through dermis and epidermis to pore
        poly!(ax_anatomy, Rect2f(dx - 0.16f0, -0.65f0, 0.32f0, 2.65f0), color = "#fff1f2", strokecolor = "#9f1239", strokewidth = 1.4)
        # Parallel resistor label beside each duct
        text!(
            ax_anatomy,
            dx + 0.50f0,
            0.25f0,
            text = "G$idx\n(1/R$idx)",
            align = (:center, :center),
            font = :bold,
            fontsize = 10.5,
            color = "#334155",
        )
    end

    # Dynamic sweat filling level & pore opening (Edelberg's Poral Valve Model)
    duct_state = lift(scl_base, scr_amp, scr_lat, scr_rise, scr_rec_half, playback_time) do base_g, amp_g, lat, rise, rec_h, pt
        t_on = t_stim + lat
        phasic_frac = pt < 0.0 ? 1.0 : eval_scr_unit(pt, t_on, rise, rec_h)
        curr_g = base_g + amp_g * phasic_frac

        # Base filling from tonic SCL (2–20 µS -> 0.22 to 0.68 of duct) + phasic SCR surge (0.1–1.0 µS -> up to +0.35)
        tonic_fill = clamp(0.22 + ((base_g - 2.0) / 18.0) * 0.46, 0.22, 0.68)
        phasic_boost = clamp((amp_g / 1.0) * phasic_frac * 0.35, 0.0, 0.35)

        tonic_rects = Rect2f[]
        phasic_rects = Rect2f[]
        droplets = Point2f[]
        thresholds = (0.0f0, 0.04f0, 0.08f0, 0.12f0)
        for (i, dx) in enumerate(duct_xs)
            f_tonic = clamp(tonic_fill - thresholds[i], 0.16, 0.68)
            f_tot = clamp(f_tonic + phasic_boost, f_tonic, 1.0)
            h_tonic = Float32(2.65 * f_tonic)
            h_phasic = Float32(2.65 * (f_tot - f_tonic))
            push!(tonic_rects, Rect2f(dx - 0.13f0, -0.65f0, 0.26f0, h_tonic))
            push!(phasic_rects, Rect2f(dx - 0.13f0, -0.65f0 + h_tonic, 0.26f0, h_phasic))
            if f_tot >= 0.86
                push!(droplets, Point2f(dx, 2.08f0))
            end
        end

        valve_stage = if phasic_frac < 0.05
            "Stage A: Baseline Tonic Fill (Pores Closed)"
        elseif phasic_frac < 0.65
            "Stage B: Intraductal Pressure Rising"
        elseif phasic_frac >= 0.65 && (pt < 0.0 || pt <= t_on + rise + 0.4)
            "Stage C: Pores Open → Steep SCR Conductance Surge"
        else
            "Stage D: Sweat Reabsorption / Evaporation → Recovery"
        end

        return (; tonic_rects, phasic_rects, droplets, curr_g, valve_stage)
    end

    poly!(ax_anatomy, lift(d -> d.tonic_rects, duct_state), color = "#7dd3fc", strokewidth = 0)
    poly!(ax_anatomy, lift(d -> d.phasic_rects, duct_state), color = "#0284c7", strokewidth = 0)
    scatter!(ax_anatomy, lift(d -> d.droplets, duct_state), color = "#38bdf8", strokecolor = "#0369a1", strokewidth = 1.5, markersize = 13)

    # Poral Valve Stage Banner
    poly!(ax_anatomy, Rect2f(-4.8, -2.62, 9.6, 0.52), color = "#f0f9ff", strokecolor = "#bae6fd", strokewidth = 1.2)
    text!(
        ax_anatomy,
        0.0,
        -2.36,
        text = lift(d -> "Poral Valve (Edelberg, 1993):  " * d.valve_stage, duct_state),
        align = (:center, :center),
        font = :bold,
        fontsize = 12.5,
        color = "#0369a1",
    )

    # Parallel Conductance vs Resistance Readout Box (Slide 34 demonstration)
    readout_text = lift(scl_base, scr_amp) do base_g, amp_g
        peak_g = base_g + amp_g
        r_base = 1000.0 / base_g
        r_peak = 1000.0 / peak_g
        delta_r = r_base - r_peak
        "G_total = ∑ G_i:  $(round(base_g, digits=1)) → $(round(peak_g, digits=1)) µS (ΔG = +$(round(amp_g, digits=2)) µS)   |   " *
        "R_total = 1/G:  $(round(r_base, digits=1)) → $(round(r_peak, digits=1)) kΩ (ΔR = –$(round(delta_r, digits=2)) kΩ)"
    end
    poly!(ax_anatomy, Rect2f(-5.05, -3.55, 10.1, 0.72), color = "#f8fafc", strokecolor = "#cbd5e1", strokewidth = 1.3)
    text!(ax_anatomy, 0.0, -3.19, text = readout_text, align = (:center, :center), font = :bold, fontsize = 12.0, color = "#0f172a")

    # =========================================================================
    # RIGHT PANEL: Single-Trial ER-SCR Morphology & Dawson Measures (Slide 32)
    # =========================================================================
    ax_scr = Axis(
        fig[1:2, 2],
        title = "Single Event-Related Response — Skin Conductance G [µS]",
        xlabel = "Time [s]",
        ylabel = "Skin Conductance G [µS]",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        backgroundcolor = :white,
        xgridcolor = "#e2e8f0",
        ygridcolor = "#e2e8f0",
    )
    xlims!(ax_scr, 0.0, t_max)

    scr_plot_data = lift(scl_base, scr_amp, scr_lat, scr_rise, scr_rec_half) do base_g, amp_g, lat, rise, rec_h
        t_on = t_stim + lat
        t_pk = t_on + rise
        t_50 = t_pk + rec_h

        y_vec = [base_g + amp_g * eval_scr_unit(t, t_on, rise, rec_h) for t in t_vec]
        y_base = base_g
        y_pk = base_g + amp_g
        y_50 = base_g + 0.5 * amp_g

        # Keep a fixed vertical span (1.15 µS) so changing SCR amplitude
        # visually grows/shrinks the peak instead of rescaling the y-axis and moving the bars
        y_lo = y_base
        span = 1.15
        pad_bot = 0.48 * span
        pad_top = 1.22 * span
        ymin = y_lo - pad_bot
        ymax = y_lo + pad_top

        # Stimulus bar near bottom
        stim_rect = Rect2f(Float32(t_stim), Float32(ymin + 0.025 * (ymax - ymin)), 1.0f0, Float32(0.055 * (ymax - ymin)))

        # Reference vertical & horizontal dashed segments (Slide 32 style)
        y_arrow_row = y_lo - 0.15 * span
        ref_segs = Point2f[
            # Baseline SCL horizontal line
            Point2f(0.4, y_base),
            Point2f(t_max - 0.4, y_base),
            # S_onset vertical line
            Point2f(t_stim, y_arrow_row - 0.06 * span),
            Point2f(t_stim, y_base + 0.25 * (y_pk - y_base)),
            # R_onset vertical dashed line
            Point2f(t_on, y_arrow_row - 0.06 * span),
            Point2f(t_on, y_pk),
            # R_peak vertical dashed line
            Point2f(t_pk, y_arrow_row - 0.06 * span),
            Point2f(t_pk, y_pk),
            # 50% recovery horizontal segment & vertical drop line
            Point2f(t_pk, y_50),
            Point2f(t_50 + 0.6, y_50),
            Point2f(t_50, y_arrow_row - 0.06 * span),
            Point2f(t_50, y_50),
        ]

        # Measurement interval lines (Latency, Rise Time, Rec t/2, Amplitude) with end ticks
        tk = Float32(0.022 * span)
        meas_segs = Point2f[
            # Latency horizontal bar: S_onset -> R_onset
            Point2f(t_stim, y_arrow_row),
            Point2f(t_on, y_arrow_row),
            Point2f(t_stim, y_arrow_row - tk),
            Point2f(t_stim, y_arrow_row + tk),
            # Rise time horizontal bar: R_onset -> R_peak
            Point2f(t_on, y_arrow_row),
            Point2f(t_pk, y_arrow_row),
            Point2f(t_on, y_arrow_row - tk),
            Point2f(t_on, y_arrow_row + tk),
            # Rec t/2 horizontal bar: R_peak -> t_50
            Point2f(t_pk, y_arrow_row),
            Point2f(t_50, y_arrow_row),
            Point2f(t_pk, y_arrow_row - tk),
            Point2f(t_pk, y_arrow_row + tk),
            Point2f(t_50, y_arrow_row - tk),
            Point2f(t_50, y_arrow_row + tk),
            # Amplitude vertical bar at R_onset + 0.35 * rise
            Point2f(t_on + 0.32 * rise, y_base),
            Point2f(t_on + 0.32 * rise, y_pk),
        ]

        amp_str = "SCR ampl.\n$(round(amp_g, digits=2)) µS"
        scl_str = "SCL\n($(round(base_g, digits=1)) µS)"

        return (;
            pts = Point2f.(t_vec, y_vec),
            ymin,
            ymax,
            y_base,
            y_pk,
            y_50,
            t_on,
            t_pk,
            t_50,
            y_arrow_row,
            stim_rect,
            ref_segs,
            meas_segs,
            amp_str,
            scl_str,
            span,
        )
    end

    on(scr_plot_data) do d
        ylims!(ax_scr, d.ymin, d.ymax)
    end
    let d0 = scr_plot_data[]
        ylims!(ax_scr, d0.ymin, d0.ymax)
    end

    # Stimulus ON bar
    poly!(ax_scr, lift(d -> d.stim_rect, scr_plot_data), color = "#1e293b")
    text!(
        ax_scr,
        t_stim + 0.5,
        lift(d -> d.ymin + 0.11 * (d.ymax - d.ymin), scr_plot_data),
        text = "Stimulus (S onset)",
        align = (:center, :bottom),
        font = :bold,
        fontsize = 11.5,
        color = "#1e293b",
    )

    # Dashed reference lines and solid measurement bars
    linesegments!(ax_scr, lift(d -> d.ref_segs, scr_plot_data), color = "#64748b", linestyle = :dash, linewidth = 1.4)
    linesegments!(ax_scr, lift(d -> d.meas_segs, scr_plot_data), color = "#0284c7", linewidth = 2.2)

    # Main EDA curve and key landmark markers (R_onset, R-peak, 50% recovery)
    lines!(ax_scr, lift(d -> d.pts, scr_plot_data), color = "#0f172a", linewidth = 3.0)
    scatter!(
        ax_scr,
        lift(d -> [Point2f(d.t_on, d.y_base), Point2f(d.t_pk, d.y_pk), Point2f(d.t_50, d.y_50)], scr_plot_data),
        color = ["#0284c7", "#e11d48", "#059669"],
        strokecolor = :white,
        strokewidth = 1.5,
        markersize = 12,
    )

    # Anatomical labels matching Slide 32
    text!(
        ax_scr,
        0.35,
        lift(d -> d.y_base + 0.03 * d.span, scr_plot_data),
        text = lift(d -> d.scl_str, scr_plot_data),
        align = (:left, :bottom),
        font = :bold,
        fontsize = 12,
        color = "#334155",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(d.t_on, d.y_arrow_row - 0.075 * d.span), scr_plot_data),
        text = "R onset",
        align = (:center, :top),
        font = :bold,
        fontsize = 12,
        color = "#0284c7",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(d.t_pk, d.y_pk + 0.04 * d.span), scr_plot_data),
        text = "R-peak",
        align = (:center, :bottom),
        font = :bold,
        fontsize = 13,
        color = "#be123c",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(0.5 * (t_stim + d.t_on), d.y_arrow_row + 0.025 * d.span), scr_plot_data),
        text = lift(_ -> "Lat. ($(round(scr_lat[], digits=1)) s)", scr_plot_data),
        align = (:center, :bottom),
        font = :bold,
        fontsize = 11.5,
        color = "#0369a1",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(0.5 * (d.t_on + d.t_pk), d.y_arrow_row + 0.025 * d.span), scr_plot_data),
        text = lift(_ -> "Rise ($(round(scr_rise[], digits=1)) s)", scr_plot_data),
        align = (:center, :bottom),
        font = :bold,
        fontsize = 11.5,
        color = "#0369a1",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(0.5 * (d.t_pk + d.t_50), d.y_arrow_row + 0.025 * d.span), scr_plot_data),
        text = lift(_ -> "Rec t/2 ($(round(scr_rec_half[], digits=1)) s)", scr_plot_data),
        align = (:center, :bottom),
        font = :bold,
        fontsize = 11.5,
        color = "#0369a1",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(d.t_on + 0.28 * scr_rise[], 0.5 * (d.y_base + d.y_pk)), scr_plot_data),
        text = lift(d -> d.amp_str, scr_plot_data),
        align = (:right, :center),
        font = :bold,
        fontsize = 11.5,
        color = "#0f172a",
    )
    text!(
        ax_scr,
        lift(d -> Point2f(d.t_50 + 0.25, d.y_50), scr_plot_data),
        text = "50% Recovery",
        align = (:left, :center),
        font = :bold,
        fontsize = 11.5,
        color = "#059669",
    )

    # Playback cursor
    lines!(
        ax_scr,
        lift(pt -> pt < 0.0 ? [Point2f(NaN, NaN)] : [Point2f(pt, scr_plot_data[].ymin), Point2f(pt, scr_plot_data[].ymax)], playback_time),
        color = "#e11d48",
        linewidth = 2.2,
        linestyle = :dash,
    )

    # =========================================================================
    # CONTROLS PANEL (Bottom-Left under Schematic)
    # =========================================================================
    ctrl_layout = fig[2, 1] = GridLayout(tellwidth = false, halign = :center)
    sl_grid = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(sl_grid[1, 1], text = "Baseline SCL:", font = :bold, fontsize = label_font, halign = :right)
    sl_scl = Slider(sl_grid[1, 2], range = 2.0:0.5:20.0, startvalue = 10.0, width = 110)
    Label(sl_grid[1, 3], text = lift(v -> "$(round(v, digits=1)) µS", sl_scl.value), width = 62, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 4], text = "SCR Ampl.:", font = :bold, fontsize = label_font, halign = :right)
    sl_amp = Slider(sl_grid[1, 5], range = 0.1:0.05:1.0, startvalue = 0.5, width = 110)
    Label(sl_grid[1, 6], text = lift(v -> "$(round(v, digits=2)) µS", sl_amp.value), width = 62, fontsize = label_font, halign = :left)

    Label(sl_grid[2, 1], text = "Latency:", font = :bold, fontsize = label_font, halign = :right)
    sl_lat = Slider(sl_grid[2, 2], range = 1.0:0.1:3.0, startvalue = 1.5, width = 110)
    Label(sl_grid[2, 3], text = lift(v -> "$(round(v, digits=1)) s", sl_lat.value), width = 62, fontsize = label_font, halign = :left)

    Label(sl_grid[2, 4], text = "Rise Time:", font = :bold, fontsize = label_font, halign = :right)
    sl_rise = Slider(sl_grid[2, 5], range = 1.0:0.1:3.0, startvalue = 2.0, width = 110)
    Label(sl_grid[2, 6], text = lift(v -> "$(round(v, digits=1)) s", sl_rise.value), width = 62, fontsize = label_font, halign = :left)

    Label(sl_grid[3, 1], text = "Rec t/2 (50%):", font = :bold, fontsize = label_font, halign = :right)
    sl_rec = Slider(sl_grid[3, 2], range = 2.0:0.5:10.0, startvalue = 4.0, width = 110)
    Label(sl_grid[3, 3], text = lift(v -> "$(round(v, digits=1)) s", sl_rec.value), width = 62, fontsize = label_font, halign = :left)

    colgap!(sl_grid, 8)
    rowgap!(sl_grid, 8)
    rowgap!(ctrl_layout, 12)

    btn_grid = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    btn_play = Button(btn_grid[1, 1], label = "▶ Trigger Stimulus", buttoncolor = "#dcfce7", fontsize = btn_font, width = 165, height = 34)
    btn_subj_a = Button(
        btn_grid[1, 2],
        label = "Low Tonic SCL (5 µS, +0.5 µS)",
        buttoncolor = "#e0f2fe",
        fontsize = btn_font,
        width = 225,
        height = 34,
    )
    btn_subj_b = Button(
        btn_grid[1, 3],
        label = "High Tonic SCL (20 µS, +0.5 µS)",
        buttoncolor = "#fef3c7",
        fontsize = btn_font,
        width = 235,
        height = 34,
    )
    colgap!(btn_grid, 10)

    on(val -> (scl_base[] = val), sl_scl.value)
    on(val -> (scr_amp[] = val), sl_amp.value)
    on(val -> (scr_lat[] = val), sl_lat.value)
    on(val -> (scr_rise[] = val), sl_rise.value)
    on(val -> (scr_rec_half[] = val), sl_rec.value)

    on(btn_subj_a.clicks) do _
        set_close_to!(sl_scl, 5.0)
        set_close_to!(sl_amp, 0.5)
    end
    on(btn_subj_b.clicks) do _
        set_close_to!(sl_scl, 20.0)
        set_close_to!(sl_amp, 0.5)
    end

    is_animating = Ref(false)
    on(btn_play.clicks) do _
        is_animating[] && return
        is_animating[] = true
        @async begin
            try
                n_frames = 90
                dt = t_max / n_frames
                for frame = 0:n_frames
                    (!is_animating[] || !events(fig).window_open[]) && break
                    playback_time[] = frame * dt
                    sleep(1.0 / 30.0)
                end
                if events(fig).window_open[]
                    playback_time[] = -1.0
                end
            finally
                is_animating[] = false
            end
        end
    end

    colsize!(fig.layout, 1, Relative(0.48))
    colsize!(fig.layout, 2, Relative(0.52))
    rowsize!(fig.layout, 1, Relative(0.74))
    rowsize!(fig.layout, 2, Relative(0.26))

    _display_demo_figure(fig)
    return fig, ax_anatomy, ax_scr
end
