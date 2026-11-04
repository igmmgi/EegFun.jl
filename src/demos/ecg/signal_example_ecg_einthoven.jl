"""
    signal_example_ecg_einthoven()

Interactive Electrocardiography (ECG) simulator illustrating cardiac depolarization/repolarization
dipoles (P, Q, R, S, T), Einthoven's triangle, Vectorcardiogram loops, and bipolar limb-lead
projections (Lead I, Lead II, Lead III).

```julia
using GLMakie, EegFun
signal_example_ecg_einthoven()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax_torso::Axis`: Anatomical torso, Einthoven's triangle, cardiac dipole & Vectorcardiogram axis
- `ax_ecg::Axis`: Calibrated ECG strip axis
"""
function signal_example_ecg_einthoven()
    _check_interactive_demo_backend("Interactive ECG Einthoven demo") || return nothing
    _set_window_title("Interactive ECG Simulator — Einthoven's Triangle & Vectorcardiogram")

    fig = _eegfun_figure(
        size = (1400, 860),
        title = "Interactive ECG Simulator — Einthoven's Triangle & Vectorcardiogram",
        backgroundcolor = :white,
    )
    title_font, label_font, tick_font, btn_font = _demo_adaptive_fonts(fig, 1400, (16, 13), (13.0, 11), (11.0, 10), (12.5, 11))

    axis_alpha_deg = Observable(60.0)
    heart_rate_bpm = Observable(72.0)
    noise_level = Observable(15.0)
    selected_lead = Observable(:lead2)
    show_vcg_loop = Observable(true)
    show_projection = Observable(true)
    playback_time = Observable(-1.0)
    is_playing = Observable(false)

    pos_ra = Point2f(-2.60, 1.80)
    pos_la = Point2f(2.60, 1.80)
    pos_ll = Point2f(0.00, -2.70)
    pos_rl = Point2f(-1.10, -2.70)
    heart_center = Point2f(0.25, 0.35)
    dipole_scale = 1.35f0

    ax_torso = Axis(fig[1, 1], title = "Einthoven's Triangle & Vectorcardiogram", titlesize = title_font, aspect = DataAspect())
    hidedecorations!(ax_torso)
    hidespines!(ax_torso)
    xlims!(ax_torso, -3.85, 3.85)
    ylims!(ax_torso, -3.65, 3.45)

    torso_half = [
        Point2f(0.55, 3.10),
        Point2f(0.55, 2.70),
        Point2f(1.60, 2.60),
        Point2f(2.90, 2.30),
        Point2f(3.40, 1.80),
        Point2f(3.20, 0.90),
        Point2f(2.55, 0.70),
        Point2f(2.40, -0.40),
        Point2f(2.15, -1.50),
        Point2f(2.45, -2.60),
        Point2f(2.20, -3.30),
        Point2f(0.00, -3.45),
    ]
    torso_pts = [[Point2f(-p[1], p[2]) for p in reverse(torso_half)]; torso_half]
    poly!(ax_torso, torso_pts, color = RGBAf(0.96, 0.97, 0.98, 1.0), strokecolor = "#94a3b8", strokewidth = 1.8)

    lines!(ax_torso, [Point2f(0.2, 2.50), Point2f(1.2, 2.45), Point2f(2.4, 2.20)], color = "#cbd5e1", linewidth = 2.0)
    lines!(ax_torso, [Point2f(-0.2, 2.50), Point2f(-1.2, 2.45), Point2f(-2.4, 2.20)], color = "#cbd5e1", linewidth = 2.0)
    poly!(ax_torso, Rect2f(-0.18, 0.80, 0.36, 1.55), color = "#e2e8f0", strokecolor = "#cbd5e1", strokewidth = 1.0)
    text!(ax_torso, 0.0, 1.55, text = "Sternum", align = (:center, :center), font = :bold, fontsize = 8.5, color = "#94a3b8")

    poly!(
        ax_torso,
        [Point2f(0.12, 1.00), Point2f(0.35, 1.35), Point2f(0.60, 1.15), Point2f(0.40, 0.85)],
        color = "#dc2626",
        strokecolor = "#991b1b",
        strokewidth = 1.0,
    )
    poly!(
        ax_torso,
        [Point2f(-0.05, 0.90), Point2f(0.15, 1.15), Point2f(0.30, 0.80), Point2f(0.10, 0.65)],
        color = "#2563eb",
        strokecolor = "#1d4ed8",
        strokewidth = 1.0,
    )
    poly!(
        ax_torso,
        [Point2f(-0.45, 0.85), Point2f(-0.05, 0.90), Point2f(0.00, 0.55), Point2f(-0.35, 0.50)],
        color = "#fecaca",
        strokecolor = "#f87171",
        strokewidth = 1.2,
    )
    poly!(
        ax_torso,
        [Point2f(0.20, 0.85), Point2f(0.70, 0.80), Point2f(0.65, 0.50), Point2f(0.15, 0.55)],
        color = "#fecaca",
        strokecolor = "#f87171",
        strokewidth = 1.2,
    )

    poly!(
        ax_torso,
        [
            Point2f(-0.35, 0.50),
            Point2f(0.00, 0.55),
            Point2f(0.65, 0.50),
            Point2f(0.75, 0.15),
            Point2f(0.68, -0.30),
            Point2f(0.40, -0.38),
            Point2f(0.05, -0.20),
            Point2f(-0.25, 0.15),
        ],
        color = "#ef4444",
        strokecolor = "#b91c1c",
        strokewidth = 1.8,
    )
    lines!(ax_torso, [Point2f(0.12, 0.52), Point2f(0.45, -0.35)], color = "#991b1b", linewidth = 1.5, linestyle = :dash)

    scatter!(ax_torso, [-0.28], [0.82], color = "#facc15", strokecolor = "#ca8a04", strokewidth = 1.5, markersize = 12)
    text!(ax_torso, -0.32, 0.92, text = "SA Node", align = (:right, :center), font = :bold, fontsize = 8.5, color = "#ca8a04")
    scatter!(ax_torso, [0.08], [0.48], color = "#10b981", strokecolor = "#047857", strokewidth = 1.5, markersize = 10)
    text!(ax_torso, 0.12, 0.38, text = "AV Node", align = (:left, :center), font = :bold, fontsize = 8.5, color = "#047857")

    T0 = 60.0 / 72.0
    cos48, sin48 = cosd(48.0), sind(48.0)

    function _beat_params(bpm_val::Float64, alpha_val::Float64)
        T = 60.0 / bpm_val
        tachy = bpm_val > 72.0
        k = tachy ? (T / T0) : 1.0
        qt_scale = tachy ? k : min(1.35, sqrt(T / T0))
        t_factor = tachy ? k : 1.0
        tT = (0.25 + 0.18 * qt_scale) * t_factor
        wT = 0.045 * qt_scale * t_factor
        sQ, cQ = sincos(deg2rad(-125.0 + (alpha_val - 60.0) * 0.25))
        sR, cR = sincos(deg2rad(alpha_val))
        sS, cS = sincos(deg2rad(-105.0 + (alpha_val - 60.0) * 0.35))
        sT, cT = sincos(deg2rad(45.0 + (alpha_val - 60.0) * 0.70))
        return (;
            T,
            tachy,
            tP = 0.12 * k,
            inv_wP = 1.0 / (0.022 * k),
            tQ = 0.22 * k,
            inv_wQ = 1.0 / (0.010 * k),
            tR = 0.25 * k,
            inv_wR = 1.0 / (0.014 * k),
            tS = 0.28 * k,
            inv_wS = 1.0 / (0.012 * k),
            tT,
            wT,
            inv_wT = 1.0 / wT,
            t_active_end = tT + 3.0 * wT,
            cQ,
            sQ,
            cR,
            sR,
            cS,
            sS,
            cT,
            sT,
        )
    end

    @inline function _dipole_at(t_in_beat::Float64, bp)
        (t_in_beat > bp.t_active_end && !bp.tachy) && return (0.0, 0.0)
        gP = 0.18 * exp(-((t_in_beat - bp.tP) * bp.inv_wP)^2)
        gQ = 0.16 * exp(-((t_in_beat - bp.tQ) * bp.inv_wQ)^2)
        gR = 1.35 * exp(-((t_in_beat - bp.tR) * bp.inv_wR)^2)
        gS = 0.32 * exp(-((t_in_beat - bp.tS) * bp.inv_wS)^2)
        gT = 0.38 * exp(-((t_in_beat - bp.tT) * bp.inv_wT)^2)
        Px = gP * cos48 + gQ * bp.cQ + gR * bp.cR + gS * bp.cS + gT * bp.cT
        Py = -(gP * sin48 + gQ * bp.sQ + gR * bp.sR + gS * bp.sS + gT * bp.sT)
        return Px, Py
    end

    current_dipole = lift(playback_time, heart_rate_bpm, axis_alpha_deg) do pt, bpm, alpha
        bp = _beat_params(bpm, alpha)
        t_in = pt < 0.0 ? bp.tR : mod(pt, bp.T)
        return _dipole_at(t_in, bp)
    end

    vcg_pts = lift(axis_alpha_deg, heart_rate_bpm) do alpha, bpm
        bp = _beat_params(bpm, alpha)
        t_end = min(bp.T, bp.tT + 2.8 * bp.wT)
        n_vcg = 160
        pts = Vector{Point2f}(undef, n_vcg + 2)
        @inbounds for i = 0:n_vcg
            px, py = _dipole_at((i / n_vcg) * t_end, bp)
            pts[i+1] = heart_center + Point2f(dipole_scale * px, dipole_scale * py)
        end
        pts[n_vcg+2] = heart_center
        return pts
    end
    lines!(ax_torso, vcg_pts, color = "#06b6d4", linewidth = 2.0, visible = show_vcg_loop, label = "Vectorcardiogram")

    dipole_tip = lift(((px, py),) -> heart_center + Point2f(dipole_scale * px, dipole_scale * py), current_dipole)
    linesegments!(ax_torso, lift(tip -> [heart_center, tip], dipole_tip), color = "#e11d48", linewidth = 3.2)
    scatter!(
        ax_torso,
        lift(tip -> [tip[1]], dipole_tip),
        lift(tip -> [tip[2]], dipole_tip),
        color = "#e11d48",
        strokecolor = "#9f1239",
        strokewidth = 1.5,
        markersize = 11,
    )

    col_lead1 = lift(s -> (s === :lead1 || s === :all3) ? "#0284c7" : "#cbd5e1", selected_lead)
    col_lead2 = lift(s -> (s === :lead2 || s === :all3) ? "#e11d48" : "#cbd5e1", selected_lead)
    col_lead3 = lift(s -> (s === :lead3 || s === :all3) ? "#10b981" : "#cbd5e1", selected_lead)

    lines!(ax_torso, [pos_ra, pos_la], color = col_lead1, linewidth = lift(s -> (s === :lead1 || s === :all3) ? 3.0 : 1.4, selected_lead))
    lines!(ax_torso, [pos_ra, pos_ll], color = col_lead2, linewidth = lift(s -> (s === :lead2 || s === :all3) ? 3.0 : 1.4, selected_lead))
    lines!(ax_torso, [pos_la, pos_ll], color = col_lead3, linewidth = lift(s -> (s === :lead3 || s === :all3) ? 3.0 : 1.4, selected_lead))

    text!(
        ax_torso,
        0.0,
        pos_la[2] + 0.16,
        text = "Lead I (0°) ➔",
        align = (:center, :bottom),
        font = :bold,
        fontsize = 10,
        color = col_lead1,
    )
    text!(
        ax_torso,
        (pos_ra[1] + pos_ll[1]) / 2 - 0.22,
        (pos_ra[2] + pos_ll[2]) / 2,
        text = "Lead II (+60°) ➔",
        align = (:right, :center),
        font = :bold,
        fontsize = 10,
        color = col_lead2,
    )
    text!(
        ax_torso,
        (pos_la[1] + pos_ll[1]) / 2 + 0.22,
        (pos_la[2] + pos_ll[2]) / 2,
        text = "➔ Lead III (+120°)",
        align = (:left, :center),
        font = :bold,
        fontsize = 10,
        color = col_lead3,
    )

    cos60, sin60 = cosd(60.0), sind(60.0)
    cos120, sin120 = cosd(120.0), sind(120.0)
    u_lead1 = Point2f(1.0, 0.0)
    u_lead2 = Point2f(cos60, -sin60)
    u_lead3 = Point2f(cos120, -sin120)

    lead_unit_vec = lift(s -> s === :lead1 ? u_lead1 : s === :lead3 ? u_lead3 : u_lead2, selected_lead)
    proj_point = lift(current_dipole, lead_unit_vec) do (px, py), u
        v_proj = Float32(px * u[1] + py * u[2])
        return heart_center + (dipole_scale * v_proj) * u
    end

    linesegments!(
        ax_torso,
        lift((tip, pproj) -> [tip, pproj], dipole_tip, proj_point),
        color = "#64748b",
        linewidth = 1.8,
        linestyle = :dash,
        visible = show_projection,
    )
    linesegments!(ax_torso, lift(pproj -> [heart_center, pproj], proj_point), color = "#2563eb", linewidth = 4.0, visible = show_projection)
    scatter!(
        ax_torso,
        lift(p -> [p[1]], proj_point),
        lift(p -> [p[2]], proj_point),
        color = "#2563eb",
        markersize = 9,
        visible = show_projection,
    )

    scatter!(ax_torso, [pos_ra, pos_la, pos_ll], color = "#cbd5e1", markersize = 24)
    scatter!(ax_torso, [pos_rl], color = "#cbd5e1", markersize = 20)
    scatter!(ax_torso, [pos_ra, pos_la, pos_ll], color = ["#ef4444", "#eab308", "#22c55e"], markersize = 14)
    scatter!(ax_torso, [pos_rl], color = "#1e293b", markersize = 12)

    text!(
        ax_torso,
        pos_ra[1] - 0.28,
        pos_ra[2],
        text = lift(s -> s in (:lead1, :lead2, :all3) ? "RA (–)" : "RA", selected_lead),
        align = (:right, :center),
        font = :bold,
        fontsize = 11,
        color = "#991b1b",
    )
    text!(
        ax_torso,
        pos_la[1] + 0.28,
        pos_la[2],
        text = lift(s -> s === :lead1 ? "LA (+)" : s === :lead3 ? "LA (–)" : s === :all3 ? "LA (+/–)" : "LA", selected_lead),
        align = (:left, :center),
        font = :bold,
        fontsize = 11,
        color = "#854d0e",
    )
    text!(
        ax_torso,
        pos_ll[1] + 0.28,
        pos_ll[2],
        text = lift(s -> s in (:lead2, :lead3, :all3) ? "LL (+)" : "LL", selected_lead),
        align = (:left, :center),
        font = :bold,
        fontsize = 11,
        color = "#166534",
    )
    text!(
        ax_torso,
        pos_rl[1] - 0.25,
        pos_rl[2],
        text = "RL (GND)",
        align = (:right, :center),
        font = :bold,
        fontsize = 9.5,
        color = "#475569",
    )

    t_window = 5.0
    n_ecg_samples = 1500
    t_arr = collect(range(0.0, t_window, length = n_ecg_samples))
    noise_template = [0.6 * sin(i * 12.9898) + 0.35 * sin(2π * 50.0 * t_arr[i]) for i = 1:n_ecg_samples]

    ax_ecg = Axis(
        fig[1, 2],
        title = "ECG Recording",
        xlabel = "Time [s]",
        ylabel = "Potential (a.u.)",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        backgroundcolor = :white,
        xgridcolor = "#fecdd3",
        xgridwidth = 1.0,
        ygridcolor = "#fecdd3",
        ygridwidth = 1.0,
        xminorgridcolor = "#fff1f2",
        xminorgridwidth = 0.5,
        yminorgridcolor = "#fff1f2",
        yminorgridwidth = 0.5,
        xminorgridvisible = true,
        yminorgridvisible = true,
    )
    xlims!(ax_ecg, 0.0, t_window)
    ylims!(ax_ecg, -0.75, 1.75)
    hlines!(ax_ecg, [0.0], color = "#f43f5e", linestyle = :dash, linewidth = 1.2)

    function generate_ecg_traces(alpha_val, bpm_val, n_lvl)
        bp = _beat_params(bpm_val, alpha_val)
        pts_I = Vector{Point2f}(undef, n_ecg_samples)
        pts_II = Vector{Point2f}(undef, n_ecg_samples)
        pts_III = Vector{Point2f}(undef, n_ecg_samples)
        n_scale = (n_lvl / 100.0) * 0.30

        @inbounds for i = 1:n_ecg_samples
            t = t_arr[i]
            px, py = _dipole_at(mod(t, bp.T), bp)
            n = n_scale * noise_template[i]
            pts_I[i] = Point2f(t, px + n)
            pts_II[i] = Point2f(t, px * cos60 - py * sin60 + n)
            pts_III[i] = Point2f(t, px * cos120 - py * sin120 + n)
        end
        return (; pts_I, pts_II, pts_III)
    end

    ecg_data = lift(generate_ecg_traces, axis_alpha_deg, heart_rate_bpm, noise_level)
    active_ecg_pts = lift((d, s) -> s === :lead1 ? d.pts_I : s === :lead3 ? d.pts_III : d.pts_II, ecg_data, selected_lead)
    active_ecg_color = lift(s -> s === :lead1 ? "#0284c7" : s === :lead3 ? "#10b981" : "#e11d48", selected_lead)
    is_all3 = lift(s -> s === :all3, selected_lead)

    lines!(ax_ecg, active_ecg_pts, color = active_ecg_color, linewidth = 2.4, visible = lift(!, is_all3))
    lines!(ax_ecg, lift(d -> d.pts_I, ecg_data), color = "#0284c7", linewidth = 2.0, label = "Lead I (0°)", visible = is_all3)
    lines!(ax_ecg, lift(d -> d.pts_II, ecg_data), color = "#e11d48", linewidth = 2.4, label = "Lead II (+60°)", visible = is_all3)
    lines!(ax_ecg, lift(d -> d.pts_III, ecg_data), color = "#10b981", linewidth = 2.0, label = "Lead III (+120°)", visible = is_all3)

    leg_ecg = axislegend(ax_ecg, position = :rt, labelsize = 12, backgroundcolor = RGBAf(1, 1, 1, 0.88))
    leg_ecg.blockscene.visible[] = false
    on(val -> (leg_ecg.blockscene.visible[] = val), is_all3)

    lines!(
        ax_ecg,
        lift(pt -> pt < 0.0 ? [Point2f(NaN, NaN)] : [Point2f(pt, -0.75), Point2f(pt, 1.75)], playback_time),
        color = "#e11d48",
        linewidth = 2.0,
        linestyle = :dash,
    )

    wave_times = lift(bpm -> begin
        bp = _beat_params(bpm, 60.0)
        (bp.tP, bp.tQ, bp.tR, bp.tS, bp.tT)
    end, heart_rate_bpm)
    wave_font = lift(f -> f * 1.20, label_font)
    wave_col = "#1e293b"

    text!(
        ax_ecg,
        lift(w -> w[1], wave_times),
        0.35,
        text = "P",
        align = (:center, :bottom),
        font = :bold,
        fontsize = wave_font,
        color = wave_col,
    )
    text!(
        ax_ecg,
        lift(w -> w[2], wave_times),
        -0.28,
        text = "Q",
        align = (:center, :top),
        font = :bold,
        fontsize = wave_font,
        color = wave_col,
    )
    text!(
        ax_ecg,
        lift(w -> w[3], wave_times),
        1.52,
        text = "R",
        align = (:center, :bottom),
        font = :bold,
        fontsize = lift(f -> f + 1.5, wave_font),
        color = wave_col,
    )
    text!(
        ax_ecg,
        lift(w -> w[4], wave_times),
        -0.45,
        text = "S",
        align = (:center, :top),
        font = :bold,
        fontsize = wave_font,
        color = wave_col,
    )
    text!(
        ax_ecg,
        lift(w -> w[5], wave_times),
        0.55,
        text = "T",
        align = (:center, :bottom),
        font = :bold,
        fontsize = wave_font,
        color = wave_col,
    )

    ctrl_layout = fig[2, 1:2] = GridLayout(tellwidth = false, halign = :center)
    slider_sub = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(slider_sub[1, 1], text = "Heart Axis (α):", font = :bold, fontsize = label_font)
    sl_axis = Slider(slider_sub[1, 2], range = -45.0:1.0:120.0, startvalue = 60.0, width = 110)
    Label(slider_sub[1, 3], text = lift(v -> "$(round(v, digits=1))°", sl_axis.value), width = 60, fontsize = label_font)

    Label(slider_sub[1, 5], text = "HR (BPM):", font = :bold, fontsize = label_font)
    sl_hr = Slider(slider_sub[1, 6], range = 10.0:1.0:220.0, startvalue = 72.0, width = 110)
    Label(slider_sub[1, 7], text = lift(v -> "$(round(Int, v)) BPM", sl_hr.value), width = 75, fontsize = label_font)

    Label(slider_sub[1, 9], text = "Noise:", font = :bold, fontsize = label_font)
    sl_noise = Slider(slider_sub[1, 10], range = 0.0:1.0:200.0, startvalue = 15.0, width = 120)

    colgap!(slider_sub, 8)
    colsize!(slider_sub, 4, Fixed(50))
    colsize!(slider_sub, 8, Fixed(50))
    rowgap!(ctrl_layout, 12)

    btn_sub = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    btn_play = Button(btn_sub[1, 1], label = "Play Heartbeat", buttoncolor = "#dcfce7", fontsize = btn_font, width = 160, height = 36)

    Label(btn_sub[1, 3], text = "Lead:", font = :bold, fontsize = label_font)
    btn_l1 = Button(btn_sub[1, 4], label = "Lead I", fontsize = btn_font, width = 85, height = 36)
    btn_l2 = Button(btn_sub[1, 5], label = "Lead II (Std)", buttoncolor = "#fee2e2", fontsize = btn_font, width = 115, height = 36)
    btn_l3 = Button(btn_sub[1, 6], label = "Lead III", fontsize = btn_font, width = 90, height = 36)
    btn_all = Button(btn_sub[1, 7], label = "All 3 (I, II, III)", fontsize = btn_font, width = 140, height = 36)

    chk_vcg = Checkbox(btn_sub[1, 9], checked = true)
    Label(btn_sub[1, 10], text = "Show Vectorcardiogram", fontsize = label_font)
    chk_proj = Checkbox(btn_sub[1, 12], checked = true)
    Label(btn_sub[1, 13], text = "Show Projection", fontsize = label_font)

    colgap!(btn_sub, 8)
    colsize!(btn_sub, 2, Fixed(24))
    colsize!(btn_sub, 8, Fixed(28))
    colsize!(btn_sub, 11, Fixed(20))

    on(val -> (axis_alpha_deg[] = val), sl_axis.value)
    on(val -> (heart_rate_bpm[] = val), sl_hr.value)
    on(val -> (noise_level[] = val), sl_noise.value)
    on(val -> (show_vcg_loop[] = val), chk_vcg.checked)
    on(val -> (show_projection[] = val), chk_proj.checked)

    function set_lead_button_colors(active_sym)
        btn_l1.buttoncolor[] = active_sym === :lead1 ? "#bae6fd" : "#f1f5f9"
        btn_l2.buttoncolor[] = active_sym === :lead2 ? "#fecdd3" : "#f1f5f9"
        btn_l3.buttoncolor[] = active_sym === :lead3 ? "#a7f3d0" : "#f1f5f9"
        btn_all.buttoncolor[] = active_sym === :all3 ? "#fbcfe8" : "#f1f5f9"
    end
    set_lead_button_colors(:lead2)

    for (btn, sym) in ((btn_l1, :lead1), (btn_l2, :lead2), (btn_l3, :lead3), (btn_all, :all3))
        let s = sym
            on(btn.clicks) do _
                selected_lead[] = s
                set_lead_button_colors(s)
            end
        end
    end

    on(btn_play.clicks) do _
        is_playing[] = !is_playing[]
        if is_playing[]
            btn_play.label[] = "Pause Heartbeat"
            btn_play.buttoncolor[] = "#fee2e2"
            @async begin
                dt = 1.0 / 30.0
                curr_t = playback_time[] < 0.0 ? 0.0 : playback_time[]
                try
                    while is_playing[] && events(fig).window_open[]
                        curr_t += dt
                        curr_t > t_window && (curr_t = 0.0)
                        playback_time[] = curr_t
                        sleep(dt)
                    end
                finally
                    if !events(fig).window_open[]
                        is_playing[] = false
                    elseif !is_playing[]
                        btn_play.label[] = "Play Heartbeat"
                        btn_play.buttoncolor[] = "#dcfce7"
                    end
                end
            end
        else
            btn_play.label[] = "Play Heartbeat"
            btn_play.buttoncolor[] = "#dcfce7"
        end
    end

    rowsize!(fig.layout, 1, Relative(0.80))
    rowsize!(fig.layout, 2, Relative(0.20))
    colsize!(fig.layout, 1, Relative(0.53))
    colsize!(fig.layout, 2, Relative(0.47))

    _display_demo_figure(fig)
    return fig, ax_torso, ax_ecg
end
